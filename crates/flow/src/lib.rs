//! esperanto-flow — pipeline orchestration (entries FASTQ PE/SE / BAM /
//! BAM+sites, species guardrail).
//!
//! Spec: docs/specs/crates/flow.md (single source of truth).
//! Single entry [`run_pipeline`]: wires qc → map → sort → scan → score → vcf
//! by entry type, artifacts landing in `<out_dir>/<stage>/`. Orchestration
//! carries zero scientific semantics: every numeric contract belongs to the
//! stage crates; flow only wires, fails early, and writes deterministically.

pub mod error;
pub mod guard;
pub mod params;
pub mod resume;
pub mod stages;
pub mod vcf;

pub use error::FlowError;
pub use params::{Entry, RunParams};
pub use stages::{
    map_stage, map_stage_with_background_options, map_stage_with_context_options,
    map_stage_with_diagnostics, map_stage_with_exact_path_options, map_stage_with_options,
    map_stage_with_pair_substitution_options, map_stage_with_short_body_options, run_from,
};

/// Run the pipeline selected by `params.entry()`.
pub fn run_pipeline(params: &RunParams) -> Result<(), FlowError> {
    stages::run(params)
}
