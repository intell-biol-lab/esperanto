//! Strict sample sheet contract shared by batch execution and joint recount.

use anyhow::{ensure, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashSet};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sample {
    pub sample: String,
    pub r1: Vec<PathBuf>,
    pub r2: Vec<PathBuf>,
    pub bam: Option<PathBuf>,
    pub sites: Option<PathBuf>,
    /// Joint candidate BED, distinct from run's 1-based sites TSV.
    pub candidates: Option<PathBuf>,
    pub metadata: BTreeMap<String, String>,
}

pub fn read_samples(path: &Path) -> Result<Vec<Sample>> {
    let path = path.canonicalize()?;
    let root = path.parent().context("sheet has no parent")?;
    let text = std::fs::read_to_string(&path)?;
    let mut lines = text
        .lines()
        .filter(|s| !s.is_empty() && !s.starts_with('#'));
    let columns: Vec<_> = lines
        .next()
        .context("empty sample sheet")?
        .split('\t')
        .collect();
    let allowed = [
        "sample",
        "r1",
        "r2",
        "bam",
        "sites",
        "candidates",
        "batch",
        "condition",
        "donor",
        "replicate",
        "strandedness",
    ];
    let mut seen = HashSet::new();
    for column in &columns {
        ensure!(allowed.contains(column), "unknown sample column {column}");
        ensure!(seen.insert(*column), "duplicate sample column {column}");
    }
    ensure!(columns.contains(&"sample"), "missing sample column");
    let mut names = HashSet::new();
    let mut inputs = HashSet::new();
    let mut rows = Vec::new();
    for (line_no, line) in lines.enumerate() {
        let values: Vec<_> = line.split('\t').collect();
        ensure!(
            values.len() == columns.len(),
            "wrong field count at data row {}",
            line_no + 1
        );
        let row: BTreeMap<_, _> = columns.iter().copied().zip(values).collect();
        let get = |key| row.get(key).copied().unwrap_or("");
        let name = get("sample");
        ensure!(
            !name.is_empty()
                && name.len() <= 128
                && name.as_bytes()[0].is_ascii_alphanumeric()
                && name
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b)),
            "unsafe sample name {name:?}"
        );
        ensure!(names.insert(name.to_owned()), "duplicate sample {name}");
        let paths = |field| -> Result<Vec<PathBuf>> {
            let value = get(field);
            if value.is_empty() {
                return Ok(Vec::new());
            }
            value
                .split(',')
                .map(|p| {
                    ensure!(
                        !p.is_empty() && p.trim() == p,
                        "invalid path in {field} for {name}"
                    );
                    let resolved = root
                        .join(p)
                        .canonicalize()
                        .with_context(|| format!("{name}: {field}: {p}"))?;
                    ensure!(
                        resolved.is_file(),
                        "input is not a file: {}",
                        resolved.display()
                    );
                    Ok(resolved)
                })
                .collect()
        };
        let r1 = paths("r1")?;
        let r2 = paths("r2")?;
        let bam = paths("bam")?;
        let sites = paths("sites")?;
        let candidates = paths("candidates")?;
        ensure!(
            bam.len() <= 1 && sites.len() <= 1 && candidates.len() <= 1,
            "one BAM/sites file per sample"
        );
        ensure!(!r1.is_empty() ^ !bam.is_empty(), "{name}: supply r1 or bam");
        ensure!(
            r2.is_empty() || (!r1.is_empty() && r1.len() == r2.len()),
            "{name}: unequal lane counts"
        );
        ensure!(
            sites.is_empty() || !bam.is_empty(),
            "{name}: sites require BAM"
        );
        ensure!(
            candidates.is_empty() || !bam.is_empty(),
            "{name}: candidates require BAM"
        );
        for input in r1.iter().chain(&r2).chain(&bam) {
            #[cfg(unix)]
            let identity = {
                use std::os::unix::fs::MetadataExt;
                let m = input.metadata()?;
                format!("{}:{}", m.dev(), m.ino())
            };
            #[cfg(not(unix))]
            let identity = input.to_string_lossy().into_owned();
            ensure!(
                inputs.insert(identity),
                "input reused across mates/lanes/samples: {}",
                input.display()
            );
        }
        let metadata: BTreeMap<String, String> =
            ["batch", "condition", "donor", "replicate", "strandedness"]
                .into_iter()
                .filter_map(|key| row.get(key).map(|v| (key.to_owned(), (*v).to_owned())))
                .collect();
        let strand = get("strandedness");
        ensure!(
            ["", "unstranded", "forward", "reverse"].contains(&strand),
            "{name}: strandedness must be unstranded|forward|reverse"
        );
        rows.push(Sample {
            sample: name.to_owned(),
            r1,
            r2,
            bam: bam.into_iter().next(),
            sites: sites.into_iter().next(),
            candidates: candidates.into_iter().next(),
            metadata,
        });
    }
    ensure!(!rows.is_empty(), "sample sheet has no samples");
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    static NEXT: AtomicUsize = AtomicUsize::new(0);

    fn fixture(text: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "esp-sheet-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        for name in ["a.fq", "b.fq", "c.fq", "d.fq", "a.bam", "a.bed"] {
            std::fs::write(dir.join(name), name).unwrap();
        }
        let sheet = dir.join("samples.tsv");
        std::fs::write(&sheet, text).unwrap();
        sheet
    }

    #[test]
    fn relative_lanes_and_biological_metadata_are_preserved() {
        let path = fixture(
            "sample\tr1\tr2\tbatch\tcondition\tdonor\nS1\ta.fq,b.fq\tc.fq,d.fq\tB1\tcase\tP1\n",
        );
        let rows = read_samples(&path).unwrap();
        assert_eq!(rows[0].r1.len(), 2);
        assert_eq!(rows[0].r1[0], path.parent().unwrap().join("a.fq"));
        assert_eq!(rows[0].metadata["condition"], "case");
        assert_eq!(rows[0].metadata["batch"], "B1");
    }

    #[test]
    fn indexed_bam_entry_can_have_candidate_bed() {
        let path = fixture("sample\tbam\tsites\nS1\ta.bam\ta.bed\n");
        let rows = read_samples(&path).unwrap();
        assert!(rows[0].r1.is_empty());
        assert!(rows[0].sites.is_some());
    }

    #[test]
    fn rejects_invalid_or_ambiguous_sheets() {
        for text in [
            "sample\tr1\n../escape\ta.fq\n",
            "sample\tr1\nS1\ta.fq\nS1\tb.fq\n",
            "sample\tr1\nS1\ta.fq\nS2\t./a.fq\n",
            "sample\tr1\tr2\nS1\ta.fq,b.fq\tc.fq\n",
            "sample\tr1\tbam\nS1\ta.fq\ta.bam\n",
            "sample\tr1\ttyop\nS1\ta.fq\tx\n",
            "sample\tr1\tr1\nS1\ta.fq\tb.fq\n",
            "sample\tr1\nS1\tmissing.fq\n",
            "sample\tr1\nS1\ta.fq\textra\n",
            "sample\tr1\n",
            "sample\tr1\tstrandedness\nS1\ta.fq\tmaybe\n",
            "sample\tr1\tsites\nS1\ta.fq\ta.bed\n",
        ] {
            assert!(read_samples(&fixture(text)).is_err(), "{text}");
        }
    }
}
