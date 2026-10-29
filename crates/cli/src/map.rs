//! `esperanto map` — read alignment (→ `esperanto_flow::map_stage`; RNA
//! 2-pass wiring lives in the flow crate, cli is a thin front).

use std::path::PathBuf;

use clap::Args;

#[derive(Args)]
pub struct MapArgs {
    /// R1 FASTQ input.
    #[arg(long)]
    r1: PathBuf,
    /// R2 FASTQ input (omit for single-end).
    #[arg(long)]
    r2: Option<PathBuf>,
    /// paidx index path.
    #[arg(long)]
    index: PathBuf,
    /// Optional GTF (junction library).
    #[arg(long)]
    gtf: Option<PathBuf>,
    /// Optional L1 engine bundle.
    #[arg(long)]
    l1_bundle: Option<PathBuf>,
    /// Output directory (raw.bam / unmapped.fq.gz / align_qc.json).
    #[arg(long)]
    out: PathBuf,
    /// Worker threads (0 = all cores).
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// Write pass2 candidate_diagnostics.jsonl before collapsed rescue (evidence only).
    #[arg(long)]
    candidate_diagnostics: bool,
    /// Disable conservative paired read-through endpoint inference (ablation).
    #[arg(long)]
    no_fragment_endpoints: bool,
    /// Disable additional ordinary high-frequency context proposals (ablation).
    #[arg(long)]
    no_repeat_context: bool,
}

pub fn run(a: MapArgs) -> anyhow::Result<()> {
    let l1 = crate::resolve::l1_bundle(&a.l1_bundle, &a.index)?;
    let (raw, _index, _l1) = esperanto_flow::map_stage_with_context_options(
        &a.index,
        a.gtf.as_deref(),
        Some(l1.as_path()),
        &a.r1,
        a.r2.as_deref(),
        &a.out,
        crate::resolve::threads(a.threads),
        a.candidate_diagnostics,
        !a.no_fragment_endpoints,
        !a.no_repeat_context,
    )?;
    eprintln!("[map] wrote {}", raw.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[derive(Parser)]
    struct Probe {
        #[command(flatten)]
        map: MapArgs,
    }

    #[test]
    fn repeat_context_has_independent_cli_ablation() {
        assert!(Probe::try_parse_from([
            "probe",
            "--r1",
            "r.fq",
            "--index",
            "r.paidx",
            "--out",
            "out",
            "--no-repeat-context",
        ])
        .is_ok());
        let default = Probe::try_parse_from([
            "probe", "--r1", "r.fq", "--index", "r.paidx", "--out", "out",
        ])
        .unwrap();
        assert!(!default.map.no_repeat_context);
        assert!(!default.map.no_fragment_endpoints);
    }
}
