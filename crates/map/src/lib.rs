//! esperanto-map — splice-aware genomic aligner.
//!
//! Spec: docs/specs/crates/map.md.

pub mod align;
pub mod bam;
pub mod chain;
pub mod diagnostics;
pub mod error;
pub mod evidence;
pub mod extend;
pub mod fasta;
pub mod fastq;
mod fragment;
pub mod gtf;
pub mod index;
pub mod index_io;
pub mod intron_chain;
pub mod mapq;
pub mod myers_ea;
pub mod pair;
pub mod pipeline;
mod repeat_context;
pub mod seed;
pub mod sinks;
pub mod span;
pub mod splice;
pub mod split;
pub mod stats;

pub use mapq::ReadAlignment;
