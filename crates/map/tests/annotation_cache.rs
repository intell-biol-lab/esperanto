use esperanto_map::fasta;
use esperanto_map::gtf::{self, CacheStatus, JunctionLib};
use std::fs;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture(std::path::PathBuf);
impl Fixture {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "esp-annotation-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn equal(a: &JunctionLib, b: &JunctionLib) {
    assert_eq!(a.junctions, b.junctions);
    assert_eq!(a.counts, b.counts);
    assert_eq!(a.by_end, b.by_end);
}
const GTF: &[u8] = b"# comment\r\nchr1\tx\texon\t1\t10\t.\t+\t.\ttranscript_id \"a\";\textra\r\nchr1\tx\texon\t21\t30\t.\t+\t.\ttranscript_id \"a\";\r\nchr1\tx\texon\t1\t10\t.\t+\t.\ttranscript_id \"b\";\nchr1\tx\texon\t21\t30\t.\t+\t.\ttranscript_id \"b\";\nunknown\tx\texon\t2\t3\t.\t-\t.\ttranscript_id \"u\";\nunknown\tx\texon\t8\t9\t.\t+\t.\ttranscript_id \"u\";\nchr1\tx\texon\t0\t1\t.\t+\t.\tgene_id \"no_tx\";\nchr1\tx\texon\t1\t2\t.\t+\t.\ttranscript_id \"bad\xff\";\nchr2\tx\texon\t5\t6\t.\t-\t.\ttranscript_id \"bad\xff\";";

#[test]
fn verified_cache_hit_content_reference_and_corruption_controls() {
    let f = Fixture::new();
    let path = f.0.join("a.gtf");
    let cache = f.0.join("cache");
    fs::write(&path, GTF).unwrap();
    let reference = fasta::parse_fasta_bytes(b">chr1\nACGT\n>chr2\nTGCA\n").unwrap();
    let direct = gtf::from_gtf(&path, |n| {
        reference.contig_index(n.as_bytes()).unwrap_or(u32::MAX)
    })
    .unwrap();
    assert_eq!(direct.counts.iter().sum::<u32>(), 4);
    let (a, status) = gtf::from_gtf_cached(&path, &reference, Some(&cache)).unwrap();
    assert_eq!(status, CacheStatus::Rebuilt);
    equal(&a, &direct);
    let (b, status) = gtf::from_gtf_cached(&path, &reference, Some(&cache)).unwrap();
    assert_eq!(status, CacheStatus::Hit);
    equal(&a, &b);
    let (_, status) = gtf::from_gtf_cached(&path, &reference, None).unwrap();
    assert_eq!(status, CacheStatus::Disabled);
    let entry = fs::read_dir(&cache)
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    let valid = fs::read(&entry).unwrap();
    fs::write(&entry, b"truncated").unwrap();
    let (c, status) = gtf::from_gtf_cached(&path, &reference, Some(&cache)).unwrap();
    assert_eq!(status, CacheStatus::Rebuilt);
    equal(&a, &c);
    // Valid JSON with a damaged payload must also be refused.
    let mut damaged: serde_json::Value = serde_json::from_slice(&valid).unwrap();
    damaged["records"][0][4] = serde_json::json!(999);
    fs::write(&entry, serde_json::to_vec(&damaged).unwrap()).unwrap();
    let (c, status) = gtf::from_gtf_cached(&path, &reference, Some(&cache)).unwrap();
    assert_eq!(status, CacheStatus::Rebuilt);
    equal(&a, &c);
    // Same byte length/path and original timestamp: metadata is insufficient.
    let original_time = fs::metadata(&path).unwrap().modified().unwrap();
    let mut changed = GTF.to_vec();
    let position = changed.windows(4).position(|w| w == b"\t21\t").unwrap();
    changed[position + 2] = b'2';
    fs::write(&path, &changed).unwrap();
    fs::File::options()
        .write(true)
        .open(&path)
        .unwrap()
        .set_times(fs::FileTimes::new().set_modified(original_time))
        .unwrap();
    let (d, status) = gtf::from_gtf_cached(&path, &reference, Some(&cache)).unwrap();
    assert_eq!(status, CacheStatus::Rebuilt);
    assert_ne!(a.junctions, d.junctions);
    let reordered = fasta::parse_fasta_bytes(b">chr2\nTGCA\n>chr1\nACGT\n").unwrap();
    let (e, status) = gtf::from_gtf_cached(&path, &reordered, Some(&cache)).unwrap();
    assert_eq!(status, CacheStatus::Rebuilt);
    assert_ne!(d.junctions, e.junctions);
    let mut other = fasta::parse_fasta_bytes(b">chr1\nACGT\n>chr2\nTGCA\n").unwrap();
    other.fasta_sha256[0] ^= 1;
    assert_eq!(
        gtf::from_gtf_cached(&path, &other, Some(&cache)).unwrap().1,
        CacheStatus::Rebuilt
    );
}

#[test]
fn cache_unavailable_and_concurrent_publication_keep_exact_library() {
    let f = Fixture::new();
    let path = f.0.join("a.gtf");
    fs::write(&path, GTF).unwrap();
    let reference = fasta::parse_fasta_bytes(b">chr1\nACGT\n>chr2\nTGCA\n").unwrap();
    let (direct, _) = gtf::from_gtf_cached(&path, &reference, None).unwrap();
    let (fallback, status) = gtf::from_gtf_cached(&path, &reference, Some(&path)).unwrap();
    assert_eq!(status, CacheStatus::Unavailable);
    equal(&direct, &fallback);
    let cache = f.0.join("concurrent");
    std::thread::scope(|s| {
        for _ in 0..8 {
            s.spawn(|| {
                let (lib, _) = gtf::from_gtf_cached(&path, &reference, Some(&cache)).unwrap();
                equal(&lib, &direct);
            });
        }
    });
    let (lib, status) = gtf::from_gtf_cached(&path, &reference, Some(&cache)).unwrap();
    assert_eq!(status, CacheStatus::Hit);
    equal(&lib, &direct);
    assert_eq!(fs::read_dir(cache).unwrap().count(), 1);
    fs::write(&path, b"# ignored\nchr1\tx\texon\tbad\t9\t.\t+\t.\ttx\n").unwrap();
    let err = gtf::from_gtf_cached(&path, &reference, Some(&f.0)).unwrap_err();
    assert!(err.to_string().contains("line 2"));
}

use esperanto_map::error::AlignError;
use esperanto_map::gtf::Junction;
use std::collections::HashMap;
use std::path::Path;
fn legacy_from_gtf<F: Fn(&str) -> u32>(
    path: &Path,
    contig_id_fn: F,
) -> Result<JunctionLib, AlignError> {
    let bytes = std::fs::read(path).map_err(|source| AlignError::FastaIo {
        path: path.display().to_string(),
        source,
    })?;
    let text = String::from_utf8_lossy(&bytes);
    let mut junctions: Vec<Junction> = Vec::new();
    // transcript_id -> exons (contig, minus_strand, 0-based half-open list)
    let mut tx: HashMap<String, TranscriptExons> = HashMap::new();

    for (idx, line) in text.lines().enumerate() {
        let lineno = idx + 1;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        if cols[2] != "exon" {
            continue;
        }
        let contig = contig_id_fn(cols[0]);
        let start: u32 = match cols[3].parse() {
            Ok(v) => v,
            Err(_) => {
                return Err(AlignError::FastaFormat {
                    line: lineno,
                    msg: format!("bad exon start '{}'", cols[3]),
                })
            }
        };
        let end: u32 = match cols[4].parse() {
            Ok(v) => v,
            Err(_) => {
                return Err(AlignError::FastaFormat {
                    line: lineno,
                    msg: format!("bad exon end '{}'", cols[4]),
                })
            }
        };
        let minus_strand = cols[6] == "-";
        let transcript = gtf::extract_attr(cols[8], "transcript_id")
            .unwrap_or_else(|| format!("__line_{}", lineno));
        let entry = tx.entry(transcript).or_insert_with(|| TranscriptExons {
            contig,
            minus_strand,
            exons: Vec::new(),
        });
        entry.exons.push((start.saturating_sub(1), end)); // 0-based half-open exon
    }

    for (_tid, t) in tx {
        if t.exons.len() < 2 {
            continue;
        }
        let mut exons = t.exons;
        exons.sort_unstable();
        for w in exons.windows(2) {
            let (.., e1_end) = w[0];
            let (e2_start, ..) = w[1];
            if e2_start > e1_end {
                junctions.push(Junction {
                    contig: t.contig,
                    start: e1_end,
                    end: e2_start,
                    minus_strand: t.minus_strand,
                });
            }
        }
    }

    Ok(JunctionLib::build(junctions))
}

struct TranscriptExons {
    contig: u32,
    minus_strand: bool,
    exons: Vec<(u32, u32)>,
}

#[test]
fn streaming_parser_matches_legacy_bytes_and_errors() {
    let f = Fixture::new();
    let path = f.0.join("parity.gtf");
    let mut cases = vec![
        GTF.to_vec(),
        vec![],
        b"\n\r\n# ignored\r\nshort\tline\n".to_vec(),
    ];
    for ending in [b"\n".as_slice(), b"\r\n", b""] {
        for attribute in [
            b"transcript_id \"t\";".as_slice(),
            b"",
            b"transcript_id \"bad\xff\";",
        ] {
            let mut data = b"chr1\tx\texon\t1\t10\t.\t-\t.\t".to_vec();
            data.extend_from_slice(attribute);
            data.extend_from_slice(ending);
            data.extend_from_slice(b"chr1\tx\texon\t20\t30\t.\t+\t.\t");
            data.extend_from_slice(attribute);
            data.extend_from_slice(ending);
            cases.push(data);
        }
    }
    cases.push(b"#x\nchr1\tx\texon\tbad\t30\t.\t+\t.\tx".to_vec());
    cases.push(b"chr1\tx\texon\t1\tbad\t.\t+\t.\tx".to_vec());
    for bytes in cases {
        fs::write(&path, bytes).unwrap();
        let contig = |name: &str| if name == "chr1" { 0 } else { 1 };
        let old = legacy_from_gtf(&path, contig);
        let new = gtf::from_gtf(&path, contig);
        match (old, new) {
            (Ok(a), Ok(b)) => equal(&a, &b),
            (Err(a), Err(b)) => assert_eq!(a.to_string(), b.to_string()),
            (a, b) => panic!("legacy/streaming mismatch: {a:?} {b:?}"),
        }
    }
}

#[cfg(unix)]
#[test]
fn piped_annotation_is_consumed_once_and_errors_are_preserved() {
    use std::io::Write;
    use std::process::{Command, Stdio};
    for (bytes, expected) in [
        (b"chr1\tx\texon\t1\t10\t.\t+\t.\ttranscript_id \"a\";\nchr1\tx\texon\t21\t30\t.\t+\t.\ttranscript_id \"a\";\n".as_slice(), "valid"),
        (b"# line 1\nchr1\tx\texon\tbad\t10\t.\t+\t.\ttranscript_id \"a\";\n".as_slice(), "error"),
    ] {
        let mut child = Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "piped_annotation_child", "--nocapture"])
            .env("ESP_TEST_PIPE", expected)
            .stdin(Stdio::piped())
            .spawn().unwrap();
        child.stdin.take().unwrap().write_all(bytes).unwrap();
        assert!(child.wait().unwrap().success());
    }
}

#[cfg(unix)]
#[test]
fn piped_annotation_child() {
    let Ok(expected) = std::env::var("ESP_TEST_PIPE") else {
        return;
    };
    let f = Fixture::new();
    let reference = fasta::parse_fasta_bytes(b">chr1\nACGT\n").unwrap();
    let result = gtf::from_gtf_cached(
        Path::new("/dev/stdin"),
        &reference,
        Some(&f.0.join("cache")),
    );
    if expected == "valid" {
        let (lib, status) = result.unwrap();
        assert_eq!(status, CacheStatus::Disabled);
        assert_eq!(lib.junctions.len(), 1);
        assert_eq!(lib.junctions[0].start, 10);
        assert_eq!(lib.junctions[0].end, 20);
    } else {
        assert!(result.unwrap_err().to_string().contains("line 2"));
    }
}

#[test]
fn cache_verification_fields_are_independently_required() {
    use sha2::{Digest, Sha256};
    let f = Fixture::new();
    let path = f.0.join("a.gtf");
    let cache = f.0.join("cache");
    fs::write(&path, GTF).unwrap();
    let reference = fasta::parse_fasta_bytes(b">chr1\nACGT\n>chr2\nTGCA\n").unwrap();
    let (baseline, _) = gtf::from_gtf_cached(&path, &reference, Some(&cache)).unwrap();
    // The invalid-UTF8 transcript begins on chr1/+ but ends on chr2/-.
    assert!(baseline
        .junctions
        .iter()
        .any(|j| j.contig == 0 && j.start == 2 && j.end == 4 && !j.minus_strand));
    let file = fs::read_dir(&cache)
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    let valid: serde_json::Value = serde_json::from_slice(&fs::read(&file).unwrap()).unwrap();
    for field in ["format", "key", "invalid_record"] {
        let mut entry = valid.clone();
        if field == "invalid_record" {
            entry["records"][0][2] = entry["records"][0][1].clone();
            let payload = serde_json::to_vec(&entry["records"]).unwrap();
            entry["payload_sha256"] = format!("{:x}", Sha256::digest(&payload)).into();
        } else {
            entry[field] = "wrong".into();
        }
        fs::write(&file, serde_json::to_vec(&entry).unwrap()).unwrap();
        let (lib, status) = gtf::from_gtf_cached(&path, &reference, Some(&cache)).unwrap();
        assert_eq!(status, CacheStatus::Rebuilt);
        equal(&lib, &baseline);
    }
    for change in ["order", "name", "length"] {
        let mut other = fasta::parse_fasta_bytes(b">chr1\nACGT\n>chr2\nTGCA\n").unwrap();
        match change {
            "order" => other.contigs.swap(0, 1),
            "name" => other.contigs[0].name = "renamed".into(),
            _ => other.contigs[0].len += 1,
        }
        assert_eq!(other.fasta_sha256, reference.fasta_sha256);
        assert_eq!(
            gtf::from_gtf_cached(&path, &other, Some(&cache)).unwrap().1,
            CacheStatus::Rebuilt
        );
    }
    // An ignored comment byte is still part of source identity.
    let mut bytes = GTF.to_vec();
    bytes[2] ^= 1;
    fs::write(&path, bytes).unwrap();
    let (lib, status) = gtf::from_gtf_cached(&path, &reference, Some(&cache)).unwrap();
    assert_eq!(status, CacheStatus::Rebuilt);
    equal(&lib, &baseline);
}
