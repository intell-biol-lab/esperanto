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
    /// Disable constant-base background competition for context proposals (ablation).
    #[arg(long)]
    no_context_background: bool,
    /// Disable mate-supported exact short-body qualification (ablation).
    #[arg(long)]
    no_mate_short_body: bool,
    /// Disable endpoint re-arbitration in a pair-supported substitution domain.
    #[arg(long)]
    no_pair_substitution: bool,
    /// Disable projection of invariant columns at uncertain homopolymer boundaries.
    #[arg(long)]
    no_boundary_core: bool,
    /// Experiment: exact MAPQ0 mates may propose a locus without pruning rivals.
    #[arg(long)]
    experimental_weak_mate_proposals: bool,
    /// Enable experimental exact-path arbitration (rejected for general adoption).
    #[arg(long, conflicts_with = "no_exact_path_competition")]
    experimental_exact_path_competition: bool,
    /// Explicitly keep experimental exact-path arbitration disabled (compatibility).
    #[arg(long, conflicts_with = "experimental_exact_path_competition")]
    no_exact_path_competition: bool,
}

pub fn run(a: MapArgs) -> anyhow::Result<()> {
    let l1 = crate::resolve::l1_bundle(&a.l1_bundle, &a.index)?;
    let (raw, _index, _l1) = esperanto_flow::map_stage_with_weak_mate_options(
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
        !a.no_context_background,
        a.experimental_exact_path_competition,
        !a.no_mate_short_body,
        !a.no_pair_substitution,
        !a.no_boundary_core,
        a.experimental_weak_mate_proposals,
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
        assert!(!default.map.no_context_background);
        assert!(!default.map.experimental_weak_mate_proposals);
    }

    #[test]
    fn weak_mate_proposals_require_explicit_opt_in() {
        let parsed = Probe::try_parse_from([
            "probe",
            "--r1",
            "r.fq",
            "--index",
            "r.paidx",
            "--out",
            "out",
            "--experimental-weak-mate-proposals",
        ])
        .unwrap();
        assert!(parsed.map.experimental_weak_mate_proposals);
    }
    #[test]
    fn context_background_has_independent_cli_ablation() {
        assert!(Probe::try_parse_from([
            "probe",
            "--r1",
            "r.fq",
            "--index",
            "r.paidx",
            "--out",
            "out",
            "--no-context-background",
        ])
        .is_ok());
    }

    #[test]
    fn boundary_core_has_independent_cli_ablation() {
        let base = [
            "probe", "--r1", "r.fq", "--index", "r.paidx", "--out", "out",
        ];
        assert!(!Probe::try_parse_from(base).unwrap().map.no_boundary_core);
        let mut args = base.to_vec();
        args.push("--no-boundary-core");
        let off = Probe::try_parse_from(args).unwrap();
        assert!(off.map.no_boundary_core);
        assert!(!off.map.no_mate_short_body);
        assert!(!off.map.no_pair_substitution);
        assert!(!off.map.no_fragment_endpoints);
    }

    #[test]
    fn pair_substitution_has_independent_cli_ablation() {
        let base = [
            "probe", "--r1", "r.fq", "--index", "r.paidx", "--out", "out",
        ];
        assert!(
            !Probe::try_parse_from(base)
                .unwrap()
                .map
                .no_pair_substitution
        );
        let mut args = base.to_vec();
        args.push("--no-pair-substitution");
        let off = Probe::try_parse_from(args).unwrap();
        assert!(off.map.no_pair_substitution);
        assert!(!off.map.no_mate_short_body);
        assert!(!off.map.no_fragment_endpoints);
        assert!(!off.map.experimental_exact_path_competition);
    }

    #[test]
    fn short_body_has_independent_cli_ablation() {
        let base = [
            "probe", "--r1", "r.fq", "--index", "r.paidx", "--out", "out",
        ];
        let default = Probe::try_parse_from(base).unwrap();
        assert!(!default.map.no_mate_short_body);
        let mut off = base.to_vec();
        off.push("--no-mate-short-body");
        let off = Probe::try_parse_from(off).unwrap();
        assert!(off.map.no_mate_short_body);
        assert!(!off.map.no_fragment_endpoints);
        assert!(!off.map.experimental_exact_path_competition);
    }

    #[test]
    fn exact_path_has_independent_cli_ablation() {
        assert!(Probe::try_parse_from([
            "probe",
            "--r1",
            "r.fq",
            "--index",
            "r.paidx",
            "--out",
            "out",
            "--no-exact-path-competition",
        ])
        .is_ok());
    }

    #[test]
    fn default_policy_rejects_conflicting_experiment_flags() {
        assert!(Probe::try_parse_from([
            "probe",
            "--r1",
            "r.fq",
            "--index",
            "r.paidx",
            "--out",
            "out",
            "--experimental-exact-path-competition",
        ])
        .is_ok());
        let conflicting = Probe::try_parse_from([
            "probe",
            "--r1",
            "r.fq",
            "--index",
            "r.paidx",
            "--out",
            "out",
            "--experimental-exact-path-competition",
            "--no-exact-path-competition",
        ])
        .err()
        .expect("conflicting policy requests must fail");
        assert_eq!(conflicting.kind(), clap::error::ErrorKind::ArgumentConflict);
    }

    #[test]
    fn default_policy_cli_and_old_flow_apis_restore_native_mapping() {
        use rust_htslib::bam::Read as _;
        let dir = std::env::temp_dir().join(format!("esp-default-policy-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut state = 0x19eace11u64;
        let mut genome: Vec<u8> = (0..2000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        let query = genome[200..350].to_vec();
        genome[1100..1102].copy_from_slice(&query[148..]);
        let fasta = dir.join("reference.fa");
        std::fs::write(
            &fasta,
            [b">chrT\n".to_vec(), genome, b"\n".to_vec()].concat(),
        )
        .unwrap();
        std::fs::write(dir.join("reference.fa.fai"), b"chrT\t2000\t6\t2000\t2001\n").unwrap();
        let gtf = dir.join("annotation.gtf");
        std::fs::write(
            &gtf,
            concat!(
                "chrT\tx\ttranscript\t151\t450\t.\t+\t.\ttranscript_id \"a\";\n",
                "chrT\tx\texon\t151\t450\t.\t+\t.\ttranscript_id \"a\";\n",
                "chrT\tx\ttranscript\t151\t1200\t.\t+\t.\ttranscript_id \"b\";\n",
                "chrT\tx\texon\t151\t348\t.\t+\t.\ttranscript_id \"b\";\n",
                "chrT\tx\texon\t1101\t1200\t.\t+\t.\ttranscript_id \"b\";\n",
            ),
        )
        .unwrap();
        let index = dir.join("reference.paidx");
        crate::index::build_all(&fasta, Some(&gtf), &[], &index, 15, 5).unwrap();
        let fq = dir.join("reads.fq");
        let reverse: Vec<u8> = query
            .iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!(),
            })
            .collect();
        let mut reads = Vec::new();
        for (name, seq) in [("forward", &query), ("reverse", &reverse)] {
            reads.extend_from_slice(format!("@{name}\n").as_bytes());
            reads.extend_from_slice(seq);
            reads.extend_from_slice(b"\n+\n");
            reads.extend(std::iter::repeat_n(b'I', 150));
            reads.push(b'\n');
        }
        std::fs::write(&fq, reads).unwrap();
        let bundle = dir.join("reference.bndl");
        for mode in 0..8 {
            let out = dir.join(format!("mode-{mode}"));
            if mode <= 2 {
                let mut argv = vec![
                    "probe".to_owned(),
                    "--r1".into(),
                    fq.display().to_string(),
                    "--index".into(),
                    index.display().to_string(),
                    "--gtf".into(),
                    gtf.display().to_string(),
                    "--out".into(),
                    out.display().to_string(),
                    "--threads".into(),
                    "1".into(),
                ];
                match mode {
                    1 => argv.push("--no-exact-path-competition".into()),
                    2 => argv.push("--experimental-exact-path-competition".into()),
                    _ => {}
                }
                run(Probe::try_parse_from(argv).unwrap().map).unwrap();
            } else {
                use esperanto_flow::*;
                match mode {
                    3 => map_stage(&index, Some(&gtf), Some(&bundle), &fq, None, &out, 1),
                    4 => map_stage_with_diagnostics(
                        &index,
                        Some(&gtf),
                        Some(&bundle),
                        &fq,
                        None,
                        &out,
                        1,
                        true,
                    ),
                    5 => map_stage_with_options(
                        &index,
                        Some(&gtf),
                        Some(&bundle),
                        &fq,
                        None,
                        &out,
                        1,
                        false,
                        true,
                    ),
                    6 => map_stage_with_context_options(
                        &index,
                        Some(&gtf),
                        Some(&bundle),
                        &fq,
                        None,
                        &out,
                        1,
                        false,
                        true,
                        true,
                    ),
                    7 => map_stage_with_background_options(
                        &index,
                        Some(&gtf),
                        Some(&bundle),
                        &fq,
                        None,
                        &out,
                        1,
                        false,
                        true,
                        true,
                        true,
                    ),
                    _ => unreachable!(),
                }
                .unwrap();
            }
            let mut reader = rust_htslib::bam::Reader::from_path(out.join("raw.bam")).unwrap();
            let records: Vec<_> = reader.records().map(Result::unwrap).collect();
            assert_eq!(records.len(), 2);
            for record in records {
                assert!(!record.is_unmapped());
                assert_eq!(
                    record.cigar().to_string(),
                    if mode == 2 { "150M" } else { "148M752N2M" },
                    "mode={mode}"
                );
                assert_eq!(record.mapq(), if mode == 2 { 0 } else { 60 }, "mode={mode}");
                assert_eq!(record.seq().as_bytes(), query);
                assert_eq!(record.qual(), &[40; 150]);
                assert_eq!(record.is_reverse(), record.qname() == b"reverse");
            }
        }
    }
}
