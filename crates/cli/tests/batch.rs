#![cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT: AtomicUsize = AtomicUsize::new(0);
fn fixture() -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "esp-batch-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    for file in ["a.fq", "b.fq", "ref.fa", "ref.fa.fai", "ref.paidx"] {
        std::fs::write(dir.join(file), file).unwrap();
    }
    std::fs::write(dir.join("l1"), "index").unwrap();
    std::fs::write(dir.join("l1.tidx"), "transcripts").unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "sample\tr1\tcondition\nA\ta.fq\tcase\nB\tb.fq\tcontrol\n",
    )
    .unwrap();
    let exe = dir.join("mock mapper");
    std::fs::write(
        &exe,
        r#"#!/usr/bin/python3
import sys, pathlib, json, time
a=sys.argv[1:]
out=pathlib.Path(a[a.index('--out')+1]); out.mkdir(parents=True,exist_ok=True)
def event(state):
    with (out.parent.parent/'events.jsonl').open('a') as f:
        f.write(json.dumps(dict(sample=out.name,state=state,t=time.monotonic()))+'\n')
event('start')
time.sleep(.1)
if (out.parent.parent/'fail-B').exists() and out.name=='B':
    print('intentional failure',file=sys.stderr); sys.exit(7)
(out/'raw.bam').write_bytes(b'unchanged test alignment')
(out/'align_qc.json').write_text(json.dumps(a))
event('end')
"#,
    )
    .unwrap();
    std::fs::set_permissions(&exe, std::fs::Permissions::from_mode(0o755)).unwrap();
    dir
}
fn command(dir: &Path) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_esperanto"));
    c.arg("batch")
        .arg("--samples")
        .arg(dir.join("samples.tsv"))
        .args([
            "--mode",
            "map",
            "--total-threads",
            "4",
            "--threads-per-sample",
            "2",
            "--memory-gib",
            "8",
            "--memory-per-sample-gib",
            "4",
        ])
        .arg("--fasta")
        .arg(dir.join("ref.fa"))
        .arg("--index")
        .arg(dir.join("ref.paidx"))
        .arg("--l1-bundle")
        .arg(dir.join("l1"))
        .arg("--executable")
        .arg(dir.join("mock mapper"))
        .arg("--out")
        .arg(dir.join("results"));
    c
}
#[test]
fn bounded_batch_and_verified_resume() {
    let dir = fixture();
    let result = command(&dir).output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let state: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.join("results/batch.json")).unwrap()).unwrap();
    assert_eq!(state["parallel"], 2);
    assert_eq!(state["jobs"][0]["status"], "succeeded");
    assert_eq!(state["jobs"][1]["status"], "succeeded");
    let text = std::fs::read_to_string(dir.join("results/events.jsonl")).unwrap();
    let mut events: Vec<serde_json::Value> = text
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    events.sort_by(|a, b| {
        a["t"]
            .as_f64()
            .unwrap()
            .partial_cmp(&b["t"].as_f64().unwrap())
            .unwrap()
    });
    let (mut active, mut peak) = (0, 0);
    for event in events {
        active += if event["state"] == "start" { 1 } else { -1 };
        peak = peak.max(active);
    }
    assert_eq!(
        peak, 2,
        "workers did not overlap or exceeded the concurrency budget"
    );
    let before = std::fs::read(dir.join("results/samples/A/raw.bam")).unwrap();
    assert!(command(&dir).arg("--resume").status().unwrap().success());
    let state: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.join("results/batch.json")).unwrap()).unwrap();
    assert_eq!(state["jobs"][0]["attempts"], 1);
    assert_eq!(
        before,
        std::fs::read(dir.join("results/samples/A/raw.bam")).unwrap()
    );
    std::fs::write(dir.join("a.fq"), "changed input").unwrap();
    assert!(!command(&dir).arg("--resume").status().unwrap().success());
}
#[test]
fn failed_sample_isolated_and_retried_without_overwriting_attempt() {
    let dir = fixture();
    std::fs::create_dir(dir.join("results")).unwrap();
    std::fs::write(dir.join("results/fail-B"), "").unwrap();
    assert!(!command(&dir).status().unwrap().success());
    assert!(dir.join("results/samples/A/raw.bam").exists());
    std::fs::remove_file(dir.join("results/fail-B")).unwrap();
    assert!(command(&dir).arg("--resume").status().unwrap().success());
    let state: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.join("results/batch.json")).unwrap()).unwrap();
    assert_eq!(state["jobs"][0]["attempts"], 1);
    assert_eq!(state["jobs"][1]["attempts"], 2);
    assert!(dir.join("results/logs/B.attempt-1.stderr.log").exists());
    assert!(dir.join("results/failed-attempts/B.attempt-1").exists());
}
#[test]
fn output_corruption_is_not_silently_skipped() {
    let dir = fixture();
    assert!(command(&dir).status().unwrap().success());
    std::fs::write(dir.join("results/samples/A/raw.bam"), "corrupted").unwrap();
    assert!(!command(&dir).arg("--resume").status().unwrap().success());
}
#[test]
fn impossible_budget_never_launches_children() {
    let dir = fixture();
    assert!(!command(&dir)
        .args(["--jobs", "0"])
        .status()
        .unwrap()
        .success());
    assert!(!dir.join("results/samples").exists());
}

#[test]
fn dry_run_is_reviewable_and_argv_does_not_invoke_a_shell() {
    let dir = fixture();
    std::fs::write(dir.join("read; touch SHOULD_NOT_EXIST.fq"), "read").unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "sample\tr1\nA\tread; touch SHOULD_NOT_EXIST.fq\n",
    )
    .unwrap();
    assert!(command(&dir).arg("--dry-run").status().unwrap().success());
    assert!(!dir.join("results/samples").exists());
    assert!(command(&dir).arg("--resume").status().unwrap().success());
    assert!(!dir.join("SHOULD_NOT_EXIST.fq").exists());
    let args: Vec<String> = serde_json::from_slice(
        &std::fs::read(dir.join("results/samples/A/align_qc.json")).unwrap(),
    )
    .unwrap();
    assert!(args
        .iter()
        .any(|a| a.ends_with("read; touch SHOULD_NOT_EXIST.fq")));
}
