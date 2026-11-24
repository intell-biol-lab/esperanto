use super::*;
use rust_htslib::bam::{self, header::HeaderRecord, record::CigarString};
use std::sync::atomic::{AtomicUsize, Ordering};

struct Fixture {
    dir: std::path::PathBuf,
    bam: std::path::PathBuf,
    fasta: std::path::PathBuf,
}

impl Fixture {
    fn new() -> Self {
        static SEQUENCE: AtomicUsize = AtomicUsize::new(0);
        let dir = std::env::temp_dir().join(format!(
            "esp-reader-reuse-{}-{}",
            std::process::id(),
            SEQUENCE.fetch_add(1, Ordering::SeqCst)
        ));
        std::fs::create_dir(&dir).unwrap();
        let fasta = dir.join("ref.fa");
        let mut sequence = "A".repeat(2000);
        sequence.replace_range(1000..1004, "TCGN");
        std::fs::write(&fasta, format!(">chrTest\n{sequence}\n")).unwrap();
        rust_htslib::faidx::build(&fasta).unwrap();
        let bam = dir.join("test.bam");
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chrTest")
                .push_tag(b"LN", 2000),
        );
        {
            let mut writer = bam::Writer::from_path(&bam, &header, bam::Format::Bam).unwrap();
            for (name, base, bq) in [
                ("conflict/1", b'A', 30),
                ("conflict/2", b'G', 30),
                ("resolved/1", b'A', 30),
                ("resolved/2", b'G', 40),
                ("ref", b'A', 30),
            ] {
                let mut record = bam::Record::new();
                record.set(
                    name.as_bytes(),
                    Some(&CigarString(vec![Cigar::Match(100)])),
                    &[base; 100],
                    &[bq; 100],
                );
                record.set_tid(0);
                record.set_pos(450);
                record.set_flags(0);
                record.set_mapq(60);
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
        Self { dir, bam, fasta }
    }

    fn open(&self) -> Result<(IndexedReader, rust_htslib::faidx::Reader)> {
        Ok((
            IndexedReader::from_path(&self.bam)?,
            rust_htslib::faidx::Reader::from_path(&self.fasta)?,
        ))
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.dir).unwrap();
    }
}

#[test]
fn many_distinct_sites_reuse_readers_instead_of_reopening_each_group() {
    let f = Fixture::new();
    let sites: Vec<_> = (1..=128).map(|i| ("chrTest".to_owned(), i * 4)).collect();
    let opens = AtomicUsize::new(0);
    let output = feature_matrix_with_readers(&AluIndex::default(), &sites, 1, || {
        opens.fetch_add(1, Ordering::SeqCst);
        f.open()
    })
    .unwrap();
    assert_eq!(output.len(), sites.len());
    assert!(output.iter().all(|r| r.is_ok()));
    assert!(
        opens.load(Ordering::SeqCst) <= 8,
        "readers opened {} times for 128 sites",
        opens.load(Ordering::SeqCst)
    );
}

#[test]
fn reused_fetch_preserves_single_site_values_order_errors_and_duplicates() {
    let f = Fixture::new();
    let sites: Vec<_> = [2000, 1, 501, 1001, 1002, 1003, 1004, 501, 201, 450, 550]
        .into_iter()
        .map(|i| ("chrTest".to_owned(), i))
        .collect();
    let alu = AluIndex::default();
    let expected: Vec<_> = sites
        .iter()
        .map(|s| {
            feature_matrix(&f.bam, &f.fasta, &alu, std::slice::from_ref(s), 1)
                .unwrap()
                .remove(0)
        })
        .collect();
    for threads in [1, 2, 4] {
        let actual = feature_matrix(&f.bam, &f.fasta, &alu, &sites, threads).unwrap();
        for (left, right) in expected.iter().zip(actual) {
            match (left, right) {
                (Ok(a), Ok(b)) => {
                    for (a, b) in a.iter().zip(b) {
                        assert!(a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()));
                    }
                }
                (Err(a), Err(b)) => assert_eq!(a, &b),
                other => panic!("row status changed: {other:?}"),
            }
        }
    }
    let opens = AtomicUsize::new(0);
    let empty = feature_matrix_with_readers(&alu, &[], 2, || {
        opens.fetch_add(1, Ordering::SeqCst);
        f.open()
    })
    .unwrap();
    assert!(empty.is_empty());
    assert_eq!(opens.load(Ordering::SeqCst), 0);
}
