//! Thin joint evidence entry: candidates and independent counts, no pooled call.
use crate::evidence;
use anyhow::{ensure, Result};
use clap::Args;
use esperanto_flow::samples::read_samples;
use esperanto_score::joint::JointSample;
use serde_json::json;
use std::path::PathBuf;

#[derive(Args)]
pub struct JointArgs {
    /// TSV with sample, bam and optional candidates (one-base BED) plus metadata.
    #[arg(long)]
    samples: PathBuf,
    #[arg(long)]
    fasta: PathBuf,
    /// New output directory; existing outputs are refused.
    #[arg(long)]
    out: PathBuf,
    /// Coverage required for an observed/no-ALT evidence state, not an editing call.
    #[arg(long, default_value_t = 10)]
    min_depth: u64,
    /// Threads for independent candidate discovery; recount is sample-serial.
    #[arg(long, default_value_t = 1)]
    threads: usize,
}

pub fn run(a: JointArgs) -> Result<()> {
    ensure!(
        a.min_depth > 0 && a.threads > 0,
        "min-depth and threads must be positive"
    );
    ensure!(!a.out.exists(), "joint output already exists");
    let samples = read_samples(&a.samples)?;
    ensure!(
        samples.iter().all(|s| s.bam.is_some()),
        "joint requires a BAM for every sample"
    );
    ensure!(
        samples.iter().all(|s| s.sites.is_none()),
        "joint uses the candidates BED column; sites is reserved for run's 1-based TSV"
    );
    let fasta = a.fasta.canonicalize()?;
    let mut inputs = vec![
        a.samples.canonicalize()?,
        fasta.clone(),
        PathBuf::from(format!("{}.fai", fasta.display())),
        std::env::current_exe()?.canonicalize()?,
    ];
    let mut joint = Vec::new();
    for sample in &samples {
        let bam = sample.bam.as_ref().unwrap();
        inputs.push(bam.clone());
        if let Some(sites) = &sample.candidates {
            inputs.push(sites.clone());
        }
        for index in [
            PathBuf::from(format!("{}.bai", bam.display())),
            bam.with_extension("bai"),
            PathBuf::from(format!("{}.csi", bam.display())),
            bam.with_extension("csi"),
        ] {
            if index.is_file() {
                inputs.push(index);
            }
        }
        joint.push(JointSample {
            sample: sample.sample.clone(),
            bam: bam.clone(),
            candidates: sample.candidates.clone(),
        });
    }
    let identities = evidence::identities(&inputs)?;
    esperanto_score::joint::recount(&joint, &fasta, &a.out, a.min_depth, a.threads)?;
    ensure!(
        identities == evidence::identities(&inputs)?,
        "joint input changed during recount; output is incomplete"
    );
    evidence::write_json(&a.out.join("samples.json"), &samples)?;
    evidence::write_json(
        &a.out.join("provenance.json"),
        &json!({
            "schema_version":1,"software_version":env!("CARGO_PKG_VERSION"),
            "source_commit":env!("ESP_BUILD_GIT_SHA"),"worktree":env!("ESP_BUILD_WORKTREE"),
            "inputs_sha256":identities,"min_depth":a.min_depth,"discovery_threads":a.threads,
            "mapq":esperanto_score::unified::MIN_MAPQ,"baseq":esperanto_score::unified::MIN_BQ,
            "counts":"informative A/C/G/T observations, RG/template overlap merged by existing unified reader",
            "semantics":"independent sample evidence; no pooled editing calls, differential tests or model changes"
        }),
    )?;
    let outputs = evidence::identities(std::slice::from_ref(&a.out))?;
    evidence::write_json(
        &a.out.join("completion.json"),
        &json!({
            "schema_version":1,"status":"COMPLETE_INDEPENDENT_RECOUNT",
            "outputs_sha256":outputs,"input_identity_verified_after_recount":true
        }),
    )?;
    eprintln!("[joint] independent evidence: {}", a.out.display());
    Ok(())
}
