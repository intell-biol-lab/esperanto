//! Optional candidate evidence. Never a source of alignments or pileup input.

use std::io::Write;

use serde::Serialize;

use crate::align::{Placement, REFUSE_MIN_ALIGNED_FRAC};
use crate::error::AlignError;
use crate::extend::{CigarOp, EditDir};
use crate::fastq::FastqRecord;
use crate::index::Index;
use crate::seed::Strand;

#[derive(Clone, Debug, Serialize)]
pub struct CandidateEvidence {
    pub rank_in_pool: usize,
    pub reference_name: String,
    pub reference_start0: u32,
    pub strand: &'static str,
    pub cigar: String,
    pub origin: &'static str,
    pub common_score: i32,
    pub editing_score_mode: &'static str,
    pub ordinary_mismatches: u32,
    pub edit_ag: u32,
    pub edit_tc: u32,
    pub informative_query_bases: u32,
    pub eligibility: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct CandidatePool {
    pub stage: &'static str,
    pub scored_candidates: usize,
    pub excluded_nonpositive: usize,
    pub excluded_invalid_geometry: usize,
    pub candidates: Vec<CandidateEvidence>,
}

impl CandidatePool {
    pub(crate) fn capture(
        index: &Index,
        read_len: usize,
        placements: &[Placement],
        stage: &'static str,
    ) -> Self {
        let mut pool = Self {
            stage,
            scored_candidates: placements.len(),
            excluded_nonpositive: 0,
            excluded_invalid_geometry: 0,
            candidates: Vec::new(),
        };
        for (rank, p) in placements.iter().enumerate() {
            if !legal_geometry(index, read_len, p) {
                pool.excluded_invalid_geometry += 1;
                continue;
            }
            if p.strict.score <= 0 {
                pool.excluded_nonpositive += 1;
                continue;
            }
            pool.candidates.push(CandidateEvidence {
                rank_in_pool: rank,
                reference_name: index.reference.contigs[p.contig as usize].name.clone(),
                reference_start0: p.pos,
                strand: if p.strand == Strand::Plus { "+" } else { "-" },
                cigar: cigar_text(&p.cigar),
                origin: p.origin.label(),
                common_score: p.strict.score,
                editing_score_mode: match p.ea_dir {
                    EditDir::Off => "off",
                    EditDir::Ag => "ag",
                    EditDir::Tc => "tc",
                },
                ordinary_mismatches: p.strict.mm,
                edit_ag: p.strict.ea_ag,
                edit_tc: p.strict.ea_tc,
                informative_query_bases: p.informative_query_coverage(),
                eligibility: if p.primary_eligible(read_len) {
                    "eligible"
                } else {
                    "below_coverage"
                },
            });
        }
        pool
    }
}

fn legal_geometry(index: &Index, read_len: usize, p: &Placement) -> bool {
    let Some(contig) = index.reference.contigs.get(p.contig as usize) else {
        return false;
    };
    let mut query = 0u64;
    let mut reference = u64::from(p.pos);
    let mut aligned_query = false;
    for (i, op) in p.cigar.iter().enumerate() {
        let n = match op {
            CigarOp::Match(n)
            | CigarOp::Ins(n)
            | CigarOp::Del(n)
            | CigarOp::RefSkip(n)
            | CigarOp::SoftClip(n) => u64::from(*n),
        };
        if n == 0 {
            return false;
        }
        match op {
            CigarOp::Match(_) => {
                query += n;
                reference += n;
                aligned_query = true;
            }
            CigarOp::Ins(_) => query += n,
            CigarOp::Del(_) => reference += n,
            CigarOp::RefSkip(_) => {
                if !aligned_query
                    || !p.cigar[i + 1..]
                        .iter()
                        .any(|x| matches!(x, CigarOp::Match(_)))
                {
                    return false;
                }
                reference += n;
            }
            CigarOp::SoftClip(_) => {
                if i != 0 && i + 1 != p.cigar.len() {
                    return false;
                }
                query += n;
            }
        }
    }
    aligned_query && query == read_len as u64 && reference <= u64::from(contig.len)
}

fn cigar_text(cigar: &[CigarOp]) -> String {
    let mut text = String::new();
    for op in cigar {
        let (n, code) = match op {
            CigarOp::Match(n) => (n, 'M'),
            CigarOp::Ins(n) => (n, 'I'),
            CigarOp::Del(n) => (n, 'D'),
            CigarOp::RefSkip(n) => (n, 'N'),
            CigarOp::SoftClip(n) => (n, 'S'),
        };
        use std::fmt::Write as _;
        write!(text, "{n}{code}").expect("write to String");
    }
    text
}

pub(crate) fn write_row(
    sink: &mut dyn Write,
    read: &FastqRecord,
    mate: Option<u8>,
    primary_placed: bool,
    independent: CandidatePool,
    mate_attempt: Option<CandidatePool>,
) -> Result<(), AlignError> {
    #[derive(Serialize)]
    struct Row<'a> {
        schema_version: u32,
        phase: &'static str,
        read_name: &'a str,
        mate: Option<u8>,
        read_length: usize,
        minimum_primary_query_bases: usize,
        pass2_primary_placed: bool,
        independent: CandidatePool,
        mate_attempt: Option<CandidatePool>,
    }
    let error = |source| AlignError::FastqIo {
        path: "<candidate_diagnostics.jsonl>".into(),
        source,
    };
    let name = std::str::from_utf8(&read.name)
        .map_err(|e| error(std::io::Error::new(std::io::ErrorKind::InvalidData, e)))?;
    let row = Row {
        schema_version: 1,
        phase: "pass2_before_collapsed_rescue",
        read_name: name,
        mate,
        read_length: read.seq.len(),
        minimum_primary_query_bases: (REFUSE_MIN_ALIGNED_FRAC * read.seq.len() as f64).ceil()
            as usize,
        pass2_primary_placed: primary_placed,
        independent,
        mate_attempt,
    };
    serde_json::to_writer(&mut *sink, &row).map_err(|e| error(std::io::Error::other(e)))?;
    sink.write_all(b"\n").map_err(error)
}

pub(crate) fn flush(sink: &mut dyn Write) -> Result<(), AlignError> {
    sink.flush().map_err(|source| AlignError::FastqIo {
        path: "<candidate_diagnostics.jsonl>".into(),
        source,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct FailingSink {
        fail_flush: bool,
    }
    impl Write for FailingSink {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if self.fail_flush {
                Ok(bytes.len())
            } else {
                Err(std::io::Error::other("write failed"))
            }
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("flush failed"))
        }
    }

    #[test]
    fn diagnostic_io_errors_and_invalid_utf8_are_reported() {
        let empty = || CandidatePool {
            stage: "independent",
            scored_candidates: 0,
            excluded_nonpositive: 0,
            excluded_invalid_geometry: 0,
            candidates: vec![],
        };
        let mut read = FastqRecord {
            name: b"quote\"and\nnewline".to_vec(),
            seq: b"A".to_vec(),
            qual: b"I".to_vec(),
        };
        let mut bytes = vec![];
        write_row(&mut bytes, &read, None, false, empty(), None).unwrap();
        assert_eq!(bytes.iter().filter(|&&x| x == b'\n').count(), 1);
        let row: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(row["read_name"], "quote\"and\nnewline");
        assert!(write_row(
            &mut FailingSink { fail_flush: false },
            &read,
            None,
            false,
            empty(),
            None
        )
        .is_err());
        assert!(flush(&mut FailingSink { fail_flush: true }).is_err());
        read.name = vec![0xff];
        bytes.clear();
        assert!(write_row(&mut bytes, &read, None, false, empty(), None).is_err());
        assert!(bytes.is_empty());
    }
}
