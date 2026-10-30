//! Disposable annotation-only cache. Content verification precedes every hit.
use super::{parse_gtf, Junction, JunctionLib};
use crate::{error::AlignError, fasta::Reference};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

// Bump when parser semantics, reference identity encoding or payload change.
const FORMAT: &str = "esperanto-gtf-junctions-v1";
const MAX_ENTRY: u64 = 128 << 20;
static NEXT: AtomicU64 = AtomicU64::new(0);
type Record = (u32, u32, u32, bool, u32);

/// Cache outcome, reported separately from scientific mapping statistics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CacheStatus {
    /// No cache directory supplied.
    Disabled,
    /// Source, reference, schema and payload were verified.
    Hit,
    /// Parsed source and atomically published a new entry.
    Rebuilt,
    /// Parsed normally; cache could not be read or published.
    Unavailable,
}

#[derive(Serialize, Deserialize)]
struct Entry {
    format: String,
    key: String,
    payload_sha256: String,
    records: Vec<Record>,
}

/// Load annotation with a disposable cache outside the reference directory.
///
/// `None` bypasses cache reads/writes. A hit always hashes the full source;
/// timestamps alone never qualify. Cache errors fall back to normal parsing.
/// Annotation and reference files must remain immutable during an invocation,
/// as for the normal parser and mmap index loaders.
pub fn from_gtf_cached(
    path: &Path,
    reference: &Reference,
    cache_dir: Option<&Path>,
) -> Result<(JunctionLib, CacheStatus), AlignError> {
    let contig_id = |n: &str| reference.contig_index(n.as_bytes()).unwrap_or(u32::MAX);
    let Some(dir) = cache_dir else {
        return super::from_gtf(path, contig_id).map(|lib| (lib, CacheStatus::Disabled));
    };
    let source_sha = hash_file(path)?;
    let mut key_hash = Sha256::new();
    key_hash.update(FORMAT);
    key_hash.update(source_sha);
    key_hash.update(reference.fasta_sha256);
    key_hash.update((reference.contigs.len() as u64).to_le_bytes());
    for c in &reference.contigs {
        key_hash.update((c.name.len() as u64).to_le_bytes());
        key_hash.update(c.name.as_bytes());
        key_hash.update(c.len.to_le_bytes());
    }
    let key = format!("{:x}", key_hash.finalize());
    let entry_path = dir.join(format!("{key}.json"));
    if let Some(lib) = read_entry(&entry_path, &key) {
        return Ok((lib, CacheStatus::Hit));
    }
    let (lib, parsed_sha) = parse_gtf(path, contig_id)?;
    // Do not publish under a digest from a different source read.
    let published = parsed_sha == source_sha && publish(dir, &entry_path, &key, &lib).is_ok();
    Ok((
        lib,
        if published {
            CacheStatus::Rebuilt
        } else {
            CacheStatus::Unavailable
        },
    ))
}

fn hash_file(path: &Path) -> Result<[u8; 32], AlignError> {
    let io_error = |source| AlignError::FastaIo {
        path: path.display().to_string(),
        source,
    };
    let mut f = File::open(path).map_err(io_error)?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = f.read(&mut buf).map_err(io_error)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(h.finalize().into())
}

fn read_entry(path: &Path, key: &str) -> Option<JunctionLib> {
    let f = File::open(path).ok()?;
    if f.metadata().ok()?.len() > MAX_ENTRY {
        return None;
    }
    let mut bytes = Vec::new();
    f.take(MAX_ENTRY + 1).read_to_end(&mut bytes).ok()?;
    if bytes.len() as u64 > MAX_ENTRY {
        return None;
    }
    let entry: Entry = serde_json::from_slice(&bytes).ok()?;
    if entry.format != FORMAT || entry.key != key {
        return None;
    }
    let payload = serde_json::to_vec(&entry.records).ok()?;
    if format!("{:x}", Sha256::digest(&payload)) != entry.payload_sha256 {
        return None;
    }
    let mut items = Vec::with_capacity(entry.records.len());
    for (contig, start, end, minus_strand, count) in entry.records {
        let j = Junction {
            contig,
            start,
            end,
            minus_strand,
        };
        if start >= end || count == 0 || items.last().is_some_and(|(last, _)| *last >= j) {
            return None;
        }
        items.push((j, count));
    }
    Some(JunctionLib::build_with_counts(items))
}

fn publish(dir: &Path, path: &Path, key: &str, lib: &JunctionLib) -> std::io::Result<()> {
    let records: Vec<Record> = lib
        .junctions
        .iter()
        .zip(&lib.counts)
        .map(|(j, &count)| (j.contig, j.start, j.end, j.minus_strand, count))
        .collect();
    let payload = serde_json::to_vec(&records)?;
    let entry = Entry {
        format: FORMAT.to_string(),
        key: key.to_string(),
        payload_sha256: format!("{:x}", Sha256::digest(&payload)),
        records,
    };
    let bytes = serde_json::to_vec(&entry)?;
    if bytes.len() as u64 > MAX_ENTRY {
        return Err(std::io::Error::other("annotation cache entry too large"));
    }
    fs::create_dir_all(dir)?;
    // create_new protects other writers, including PID reuse and concurrent threads.
    let (tmp, mut f) = loop {
        let tmp = dir.join(format!(
            ".{key}.{}.{}.tmp",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        match OpenOptions::new().write(true).create_new(true).open(&tmp) {
            Ok(f) => break (tmp, f),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e),
        }
    };
    let result = (|| {
        f.write_all(&bytes)?;
        f.sync_all()?;
        drop(f);
        fs::rename(&tmp, path)
    })();
    if result.is_err() {
        let _ = fs::remove_file(tmp);
    }
    result
}
