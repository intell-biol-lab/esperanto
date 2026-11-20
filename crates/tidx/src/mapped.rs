//! Read-only mapped table ownership and alignment-independent little-endian access.

use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::{Arc, OnceLock};

use memmap2::Mmap;

use super::{le32, le64, parse_index_layout, HEADER_LEN};

/// Borrowed hit view. Encoded records are decoded only as callers iterate them.
#[derive(Debug, Clone, Copy)]
pub struct LookupHits<'a> {
    data: HitData<'a>,
}

#[derive(Debug, Clone, Copy)]
enum HitData<'a> {
    Decoded(&'a [(u32, u32)]),
    Encoded(&'a [u8]),
}

impl<'a> LookupHits<'a> {
    /// Adapt an existing index provider without changing its lookup contract.
    pub fn from_slice(hits: &'a [(u32, u32)]) -> Self {
        Self {
            data: HitData::Decoded(hits),
        }
    }

    /// Number of hits, including every duplicate in original table order.
    pub fn len(&self) -> usize {
        match self.data {
            HitData::Decoded(hits) => hits.len(),
            HitData::Encoded(bytes) => bytes.len() / 8,
        }
    }

    /// Whether the query missed.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterate without allocating or assuming on-disk alignment/host byte order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (u32, u32)> + 'a {
        let data = self.data;
        (0..self.len()).map(move |i| match data {
            HitData::Decoded(hits) => hits[i],
            HitData::Encoded(bytes) => {
                let entry = &bytes[i * 8..i * 8 + 8];
                (le32(&entry[..4]), le32(&entry[4..]))
            }
        })
    }
}

/// Demand-paged transcript index, sharing mapping ownership across clones.
///
/// Index files must remain immutable while open, as for the genome index.
/// Header, section bounds, names and bucket validation match the owning reader.
/// The normal [`Self::lookup_hits`] path never materializes whole entry tables.
#[derive(Debug, Clone)]
pub struct MappedTidx {
    mapping: Arc<Mmap>,
    k: u32,
    tx_count: u32,
    names: Vec<String>,
    bucket_offsets: Vec<u64>,
    key_start: usize,
    payload_start: usize,
    key_width: usize,
    key_mask: u64,
    bucket_shift: u32,
    /// Compatibility only: initialized once if a caller requires tuple slices.
    pub(super) legacy_payloads: Arc<OnceLock<Vec<(u32, u32)>>>,
}

impl MappedTidx {
    /// Open an immutable index using the existing format's validation.
    pub fn open(path: &Path) -> io::Result<Self> {
        let file = File::open(path)?;
        let len = file.metadata()?.len();
        if len < HEADER_LEN as u64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                super::Error::Truncated(format!(
                    "file is {len} bytes, shorter than the {HEADER_LEN}-byte header"
                )),
            ));
        }
        // SAFETY: callers must keep the index inode immutable while open, the
        // same precondition as existing read-only genome-index mappings. Arc
        // owns the map; all returned views borrow self, with no static casts.
        let mapping = Arc::new(unsafe { Mmap::map(&file) }?);
        let parsed = parse_index_layout(&mapping)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        // Layout validation proves these tables are adjacent and end at EOF.
        let payload_start = mapping.len() - parsed.payloads.len();
        let key_start = payload_start - parsed.keys.len();
        Ok(Self {
            mapping: Arc::clone(&mapping),
            k: parsed.k,
            tx_count: parsed.tx_count,
            names: parsed.names,
            bucket_offsets: parsed.bucket_offsets,
            key_start,
            payload_start,
            key_width: parsed.key_width,
            key_mask: parsed.key_mask,
            bucket_shift: parsed.bucket_shift,
            legacy_payloads: Arc::new(OnceLock::new()),
        })
    }

    /// Indexed k-mer length.
    pub fn k(&self) -> u32 {
        self.k
    }

    /// Number of transcripts in the dense id space.
    pub fn tx_count(&self) -> u32 {
        self.tx_count
    }

    /// Transcript display name, or empty string for an unknown id.
    pub fn transcript_name(&self, tx_id: u32) -> &str {
        self.names.get(tx_id as usize).map_or("", String::as_str)
    }

    fn key(&self, i: usize) -> u64 {
        let start = self.key_start + i * self.key_width;
        let bytes = &self.mapping[start..start + self.key_width];
        if self.key_width == 4 {
            u64::from(le32(bytes))
        } else {
            le64(bytes)
        }
    }

    fn range(&self, canonical: u64) -> std::ops::Range<usize> {
        let bucket = (canonical >> self.bucket_shift) as usize;
        if bucket >= self.bucket_offsets.len() - 1 {
            return 0..0;
        }
        let start = self.bucket_offsets[bucket] as usize;
        let end = self.bucket_offsets[bucket + 1] as usize;
        let low = canonical & self.key_mask;
        let partition = |inclusive: bool| {
            let (mut left, mut right) = (start, end);
            while left < right {
                let mid = left + (right - left) / 2;
                let key = self.key(mid);
                if key < low || (inclusive && key == low) {
                    left = mid + 1;
                } else {
                    right = mid;
                }
            }
            left
        };
        partition(false)..partition(true)
    }

    /// All hits, in exactly the owning reader's order, without table allocation.
    pub fn lookup_hits(&self, canonical: u64) -> LookupHits<'_> {
        let range = self.range(canonical);
        LookupHits {
            data: HitData::Encoded(
                &self.mapping
                    [self.payload_start + range.start * 8..self.payload_start + range.end * 8],
            ),
        }
    }

    /// Compatibility slice API. Materializes the payload table once if used.
    ///
    /// Performance-sensitive code should use [`Self::lookup_hits`] instead.
    pub fn lookup(&self, canonical: u64) -> &[(u32, u32)] {
        let range = self.range(canonical);
        let payloads = self.legacy_payloads.get_or_init(|| {
            self.mapping[self.payload_start..]
                .chunks_exact(8)
                .map(|c| (le32(&c[..4]), le32(&c[4..])))
                .collect()
        });
        &payloads[range]
    }
}
