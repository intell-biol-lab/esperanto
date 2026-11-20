use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Cigar, CigarString},
};
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};
static NEXT: AtomicUsize = AtomicUsize::new(0);

fn fixture() -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "esp-joint-cli-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("ref.fa"), ">chr1\nAAAAAAAAAA\n").unwrap();
    rust_htslib::faidx::build(dir.join("ref.fa")).unwrap();
    let mut header = bam::Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"SO", "coordinate"));
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 10),
    );
    let path = dir.join("A.bam");
    let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
    let mut record = bam::Record::new();
    record.set(
        b"molecule",
        Some(&CigarString(vec![Cigar::Match(2)])),
        b"GA",
        &[30, 30],
    );
    record.set_tid(0);
    record.set_pos(0);
    record.set_flags(0);
    record.set_mapq(60);
    writer.write(&record).unwrap();
    drop(writer);
    bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
    std::fs::write(dir.join("A.bed"), "chr1\t0\t1\nchr1\t2\t3\n").unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "sample\tbam\tcandidates\tcondition\tbatch\nA\tA.bam\tA.bed\tcase\ttechnical1\n",
    )
    .unwrap();
    dir
}
#[test]
fn joint_cli_publishes_independent_evidence_and_completion_without_touching_bam() {
    let dir = fixture();
    let before = std::fs::read(dir.join("A.bam")).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_esperanto"))
        .arg("joint")
        .arg("--samples")
        .arg(dir.join("samples.tsv"))
        .arg("--fasta")
        .arg(dir.join("ref.fa"))
        .arg("--out")
        .arg(dir.join("joint"))
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = std::fs::read_to_string(dir.join("joint/evidence.tsv")).unwrap();
    assert!(text.contains("chr1\t1\tA\tG\tA\t0\t0\t1\t0\t1\t1\t1.00000000\tINSUFFICIENT_COVERAGE"));
    assert!(text.contains("chr1\t3\tA\tG\tA\t0\t0\t0\t0\t0\t0\tNA\tINSUFFICIENT_COVERAGE"));
    assert_eq!(before, std::fs::read(dir.join("A.bam")).unwrap());
    let done: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.join("joint/completion.json")).unwrap()).unwrap();
    assert_eq!(done["status"], "COMPLETE_INDEPENDENT_RECOUNT");
    let metadata = std::fs::read_to_string(dir.join("joint/samples.json")).unwrap();
    assert!(metadata.contains("technical1") && metadata.contains("case"));
}
#[test]
fn joint_rejects_ambiguous_coordinate_format() {
    let dir = fixture();
    std::fs::write(
        dir.join("samples.tsv"),
        "sample\tbam\tsites\nA\tA.bam\tA.bed\n",
    )
    .unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_esperanto"))
        .arg("joint")
        .arg("--samples")
        .arg(dir.join("samples.tsv"))
        .arg("--fasta")
        .arg(dir.join("ref.fa"))
        .arg("--out")
        .arg(dir.join("joint"))
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("candidates BED column"));
    assert!(!dir.join("joint").exists());
}
