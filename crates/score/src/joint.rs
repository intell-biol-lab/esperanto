//! Candidate union and independent molecular recount, without pooled calls.
use anyhow::{ensure, Context, Result};
use rust_htslib::bam::{IndexedReader, Read as _};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
pub struct JointSample {
    pub sample: String,
    pub bam: PathBuf,
    /// Existing scan BED, or None to use the unchanged discovery implementation.
    pub candidates: Option<PathBuf>,
}

pub fn recount(
    samples: &[JointSample],
    fasta: &Path,
    out: &Path,
    min_depth: u64,
    threads: usize,
) -> Result<()> {
    ensure!(
        !samples.is_empty() && min_depth > 0 && threads > 0,
        "samples, min-depth and threads must be positive"
    );
    ensure!(
        !out.exists(),
        "joint output already exists: {}",
        out.display()
    );
    let reference = rust_htslib::faidx::Reader::from_path(fasta)?;
    let fai = std::fs::read_to_string(format!("{}.fai", fasta.display()))?;
    let mut lengths = BTreeMap::new();
    for line in fai.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        ensure!(fields.len() >= 2, "invalid FASTA index");
        lengths.insert(fields[0].to_owned(), fields[1].parse::<u64>()?);
    }
    let mut names = BTreeSet::new();
    for sample in samples {
        ensure!(
            !sample.sample.is_empty()
                && sample.sample.as_bytes()[0].is_ascii_alphanumeric()
                && sample
                    .sample
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b)),
            "unsafe sample identifier"
        );
        ensure!(names.insert(&sample.sample), "duplicate joint sample");
        let bam = IndexedReader::from_path(&sample.bam)?;
        for tid in 0..bam.header().target_count() {
            let chrom = std::str::from_utf8(bam.header().tid2name(tid))?;
            ensure!(
                lengths.get(chrom).copied() == bam.header().target_len(tid),
                "{}: BAM/reference dictionary mismatch at {chrom}",
                sample.sample
            );
        }
    }
    let parent = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    std::fs::create_dir_all(parent)?;
    let stage = parent.join(format!(
        ".{}.joint-{}.tmp",
        out.file_name()
            .context("joint output needs a directory name")?
            .to_string_lossy(),
        std::process::id()
    ));
    std::fs::create_dir(&stage)?;
    // Preserve an incomplete staging directory on error, never publish completion.
    let mut union: BTreeMap<(String, i64), BTreeSet<String>> = BTreeMap::new();
    for sample in samples {
        let candidate = if let Some(path) = &sample.candidates {
            path.clone()
        } else {
            let path = stage.join(format!("{}.discovered.bed", sample.sample));
            crate::unified::discover_sites(&sample.bam, fasta, &path, threads)?;
            path
        };
        for (line_no, line) in std::fs::read_to_string(&candidate)?.lines().enumerate() {
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track ")
                || line.starts_with("browser ")
            {
                continue;
            }
            let fields: Vec<_> = line.split('\t').collect();
            ensure!(
                fields.len() >= 3,
                "{}: invalid BED row {}",
                candidate.display(),
                line_no + 1
            );
            let start: i64 = fields[1].parse()?;
            let end: i64 = fields[2].parse()?;
            ensure!(
                start >= 0 && start.checked_add(1) == Some(end),
                "candidate BED requires one-base intervals"
            );
            let length = lengths
                .get(fields[0])
                .context("candidate contig absent from FASTA")?;
            ensure!(end as u64 <= *length, "candidate outside FASTA");
            union
                .entry((fields[0].to_owned(), end))
                .or_default()
                .insert(sample.sample.clone());
        }
    }
    let sites: Vec<_> = union.keys().cloned().collect();
    let mut alleles = Vec::new();
    for (chrom, pos) in &sites {
        let base = reference.fetch_seq(chrom, (*pos - 1) as usize, (*pos - 1) as usize)?[0]
            .to_ascii_uppercase();
        ensure!(
            matches!(base, b'A' | b'T'),
            "candidate reference must be A or T at {chrom}:{pos}"
        );
        alleles.push((base, if base == b'A' { b'G' } else { b'C' }));
    }
    let mut origins = BufWriter::new(std::fs::File::create(stage.join("candidate_origins.tsv"))?);
    writeln!(origins, "chrom\tpos\tdiscovered_in")?;
    let mut bed = BufWriter::new(std::fs::File::create(stage.join("candidates.bed"))?);
    for ((chrom, pos), sources) in &union {
        writeln!(bed, "{chrom}\t{}\t{pos}", pos - 1)?;
        writeln!(
            origins,
            "{chrom}\t{pos}\t{}",
            sources.iter().cloned().collect::<Vec<_>>().join(",")
        )?;
    }
    origins.flush()?;
    bed.flush()?;
    let mut matrix = BufWriter::new(std::fs::File::create(stage.join("evidence.tsv"))?);
    writeln!(
        matrix,
        "chrom\tpos\tref\talt\tsample\tA\tC\tG\tT\tdepth\talt_count\talt_fraction\tstate"
    )?;
    let mut ordered: Vec<_> = samples.iter().collect();
    ordered.sort_by(|a, b| a.sample.cmp(&b.sample));
    for sample in ordered {
        let counts = crate::unified::count_sites(&sample.bam, &sites)?;
        for (((chrom, pos), &(reference, alt)), counts) in sites.iter().zip(&alleles).zip(counts) {
            let (reference, alt) = (reference as char, alt as char);
            if let Some(counts) = counts {
                let depth: u64 = counts.iter().sum();
                let alt_count = counts[if alt == 'G' { 2 } else { 1 }];
                let fraction = if depth == 0 {
                    "NA".to_owned()
                } else {
                    format!("{:.8}", alt_count as f64 / depth as f64)
                };
                let state = if depth < min_depth {
                    "INSUFFICIENT_COVERAGE"
                } else if alt_count > 0 {
                    "ALT_OBSERVED"
                } else {
                    "NO_ALT_OBSERVED"
                };
                writeln!(matrix,"{chrom}\t{pos}\t{reference}\t{alt}\t{}\t{}\t{}\t{}\t{}\t{depth}\t{alt_count}\t{fraction}\t{state}",
                    sample.sample,counts[0],counts[1],counts[2],counts[3])?;
            } else {
                writeln!(matrix,"{chrom}\t{pos}\t{reference}\t{alt}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tCONTIG_ABSENT",sample.sample)?;
            }
        }
    }
    matrix.flush()?;
    // Close all files before atomic directory publication.
    drop(matrix);
    drop(origins);
    drop(bed);
    std::fs::rename(&stage, out)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{
        self,
        header::HeaderRecord,
        record::{Aux, Cigar, CigarString},
    };
    use std::sync::atomic::{AtomicUsize, Ordering};
    static NEXT: AtomicUsize = AtomicUsize::new(0);

    fn fixture() -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "esp-joint-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("ref.fa"),
            ">chr1\nAAAAAAAAAATTTTTTTTTTAAAAAAAAAA\n>chr2\nAAAAAAAAAA\n",
        )
        .unwrap();
        rust_htslib::faidx::build(dir.join("ref.fa")).unwrap();
        dir
    }

    fn bam(
        dir: &Path,
        name: &str,
        records: &[(&str, i64, &[u8], u8, u16, &str)],
        chr1_len: u64,
    ) -> PathBuf {
        let mut header = bam::Header::new();
        header.push_record(HeaderRecord::new(b"HD").push_tag(b"SO", "coordinate"));
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", chr1_len),
        );
        let path = dir.join(format!("{name}.bam"));
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
        for &(name, pos, seq, mapq, flags, rg) in records {
            let mut r = bam::Record::new();
            let qual = vec![if name == "low-bq" { 0 } else { 30 }; seq.len()];
            r.set(
                name.as_bytes(),
                Some(&CigarString(vec![Cigar::Match(seq.len() as u32)])),
                seq,
                &qual,
            );
            r.set_tid(0);
            r.set_pos(pos);
            r.set_mapq(mapq);
            r.set_flags(flags);
            r.push_aux(b"RG", Aux::String(rg)).unwrap();
            writer.write(&r).unwrap();
        }
        drop(writer);
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
        path
    }

    fn sample(dir: &Path, name: &str, bam: PathBuf, bed: &str) -> JointSample {
        let path = dir.join(format!("{name}.bed"));
        std::fs::write(&path, bed).unwrap();
        JointSample {
            sample: name.to_owned(),
            bam,
            candidates: Some(path),
        }
    }
    fn matrix(dir: &Path) -> Vec<Vec<String>> {
        std::fs::read_to_string(dir.join("evidence.tsv"))
            .unwrap()
            .lines()
            .skip(1)
            .map(|l| l.split('\t').map(str::to_owned).collect())
            .collect()
    }

    #[test]
    fn union_recounts_private_shared_and_uncovered_sites_without_pooling() {
        let dir = fixture();
        let a = bam(
            &dir,
            "A",
            &[
                ("overlap/1", 0, b"G", 60, 0, "r1"),
                ("overlap/2", 0, b"G", 60, 0, "r1"),
                ("overlap", 0, b"G", 60, 0, "r2"),
                ("filtered", 0, b"G", 0, 0, "r1"),
                ("secondary", 0, b"G", 60, 256, "r1"),
                ("supplementary", 0, b"G", 60, 2048, "r1"),
                ("low-bq", 0, b"G", 60, 0, "r1"),
                ("ref", 0, b"A", 60, 0, "r1"),
                ("conflict/1", 0, b"A", 60, 0, "r1"),
                ("conflict/2", 0, b"G", 60, 0, "r1"),
                ("tc", 10, b"C", 60, 0, "r1"),
            ],
            30,
        );
        let b = bam(&dir, "B", &[("ref", 0, b"A", 60, 0, "r1")], 30);
        let samples = vec![
            sample(&dir, "A", a, "chr1\t0\t1\nchr1\t10\t11\n"),
            sample(&dir, "B", b, "chr1\t1\t2\nchr2\t0\t1\n"),
        ];
        let out = dir.join("joint");
        recount(&samples, &dir.join("ref.fa"), &out, 1, 1).unwrap();
        let rows = matrix(&out);
        assert_eq!(rows.len(), 8);
        let row = |chrom: &str, pos: &str, sample: &str| {
            rows.iter()
                .find(|r| r[0] == chrom && r[1] == pos && r[4] == sample)
                .unwrap()
        };
        let ra = row("chr1", "1", "A");
        assert_eq!(&ra[5..12], ["1", "0", "2", "0", "3", "2", "0.66666667"]);
        assert_eq!(ra[12], "ALT_OBSERVED");
        assert_eq!(row("chr1", "1", "B")[12], "NO_ALT_OBSERVED");
        assert_eq!(row("chr1", "2", "A")[11], "NA");
        assert_eq!(row("chr1", "2", "A")[12], "INSUFFICIENT_COVERAGE");
        assert_eq!(row("chr1", "11", "A")[3], "C");
        assert_eq!(row("chr2", "1", "A")[12], "CONTIG_ABSENT");
        let out2 = dir.join("single");
        recount(&samples[..1], &dir.join("ref.fa"), &out2, 1, 1).unwrap();
        for row in matrix(&out2) {
            assert!(rows.contains(&row), "adding B changed A evidence");
        }
    }

    #[test]
    fn low_depth_retains_alt_evidence_and_dictionary_mismatch_fails() {
        let dir = fixture();
        let a = bam(&dir, "A", &[("alt", 0, b"G", 60, 0, "r1")], 30);
        let samples = vec![sample(&dir, "A", a, "chr1\t0\t1\n")];
        let out = dir.join("joint");
        recount(&samples, &dir.join("ref.fa"), &out, 10, 1).unwrap();
        let rows = matrix(&out);
        assert_eq!(rows[0][10], "1");
        assert_eq!(rows[0][12], "INSUFFICIENT_COVERAGE");
        let wrong = bam(&dir, "wrong", &[], 31);
        let bad = vec![sample(&dir, "wrong", wrong, "chr1\t0\t1\n")];
        assert!(recount(&bad, &dir.join("ref.fa"), &dir.join("bad"), 1, 1).is_err());
        assert!(!dir.join("bad/completion.json").exists());
    }

    #[test]
    fn malformed_or_out_of_reference_candidates_are_rejected() {
        let dir = fixture();
        let a = bam(&dir, "A", &[], 30);
        for (i, bed) in [
            "chr1\t-1\t0\n",
            "chr1\t1\t3\n",
            "chr1\t30\t31\n",
            "unknown\t0\t1\n",
        ]
        .iter()
        .enumerate()
        {
            let s = vec![sample(&dir, "A", a.clone(), bed)];
            assert!(recount(&s, &dir.join("ref.fa"), &dir.join(format!("bad-{i}")), 1, 1).is_err());
        }
    }

    #[test]
    fn unchanged_discovery_and_reordered_samples_give_same_matrix() {
        let dir = fixture();
        let a = bam(&dir, "A", &[("private", 0, b"G", 60, 0, "r1")], 30);
        let b = bam(
            &dir,
            "B",
            &[
                ("shared", 0, b"G", 60, 0, "r1"),
                ("private", 10, b"C", 60, 0, "r1"),
            ],
            30,
        );
        let mut samples = vec![
            JointSample {
                sample: "B".into(),
                bam: b,
                candidates: None,
            },
            JointSample {
                sample: "A".into(),
                bam: a,
                candidates: None,
            },
        ];
        recount(&samples, &dir.join("ref.fa"), &dir.join("first"), 1, 2).unwrap();
        samples.reverse();
        recount(&samples, &dir.join("ref.fa"), &dir.join("second"), 1, 1).unwrap();
        assert_eq!(matrix(&dir.join("first")), matrix(&dir.join("second")));
        assert_eq!(matrix(&dir.join("first")).len(), 4);
    }

    #[test]
    fn empty_candidate_union_is_valid_and_duplicate_samples_are_refused() {
        let dir = fixture();
        let a = bam(&dir, "A", &[], 30);
        let samples = vec![sample(&dir, "A", a, "")];
        recount(&samples, &dir.join("ref.fa"), &dir.join("empty"), 1, 1).unwrap();
        assert!(matrix(&dir.join("empty")).is_empty());
        let duplicate = vec![samples[0].clone(), samples[0].clone()];
        assert!(recount(
            &duplicate,
            &dir.join("ref.fa"),
            &dir.join("duplicate"),
            1,
            1
        )
        .is_err());
    }
}
