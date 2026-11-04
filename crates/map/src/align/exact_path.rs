//! Full-query exact paths are sequence-equivalent even with different priors.
use super::{footprint_overlap, CigarOp, Placement};

fn is_exact(p: &Placement, read_len: usize) -> bool {
    read_len > 0
        && p.strict.score > 0
        && p.strict.matches as usize == read_len
        && p.strict.aligned as usize == read_len
        && p.strict.mm == 0
        && p.strict.ea_ag == 0
        && p.strict.ea_tc == 0
        && p.strict.ins == 0
        && p.strict.del == 0
        && p.strict.clipped == 0
        && p.cigar.iter().all(|op| match op {
            CigarOp::Match(n) | CigarOp::RefSkip(n) => *n > 0,
            _ => false,
        })
}

/// Ordered M intervals encode the complete query mapping for M/N-only paths.
/// Coalescing adjacent M avoids treating a CIGAR spelling as new evidence.
fn mapping(p: &Placement) -> Vec<(u32, u32)> {
    let mut blocks: Vec<(u32, u32)> = Vec::new();
    let mut reference = p.pos;
    for op in &p.cigar {
        match op {
            CigarOp::Match(n) => {
                if let Some(last) = blocks.last_mut().filter(|last| last.1 == reference) {
                    last.1 += *n;
                } else {
                    blocks.push((reference, reference + *n));
                }
                reference += *n;
            }
            CigarOp::RefSkip(n) => reference += *n,
            _ => unreachable!("mapping is only used for exact M/N paths"),
        }
    }
    blocks
}

pub(super) fn distinct(a: &Placement, b: &Placement, read_len: usize) -> bool {
    is_exact(a, read_len) && is_exact(b, read_len) && mapping(a) != mapping(b)
}

pub(super) fn representative<'a>(
    placements: &'a [Placement],
    leader: &'a Placement,
    read_len: usize,
) -> (&'a Placement, bool) {
    if !is_exact(leader, read_len) {
        return (leader, false);
    }
    let compatible = |p: &&Placement| is_exact(p, read_len) && footprint_overlap(leader, p) >= 0.5;
    let leader_mapping = mapping(leader);
    if !placements
        .iter()
        .filter(compatible)
        .any(|p| mapping(p) != leader_mapping)
    {
        return (leader, false);
    }
    let representative = placements
        .iter()
        .filter(compatible)
        .min_by_key(|p| (p.strict.n_splice, p.origin.rank(), p.pos, mapping(p)))
        .expect("the exact leader belongs to its own local set");
    (representative, true)
}
