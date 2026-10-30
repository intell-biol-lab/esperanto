//! Tail split rescue (STAR-style split discovery).
//!
//! When the primary alignment leaves a soft-clipped end tail of at least
//! `TAIL_MIN` bases, the tail is relocated across a splice junction. Channel A
//! (library-driven, the main channel) enumerates library junctions in the
//! donor/acceptor ±20bp window and retries the tail start ±5bp, scoring both
//! segments by extension. Channels B (direct genomic scan) and C (tail
//! seeding) provide de-novo fallback when no library candidate succeeds.
//!
//! All CIGAR products conserve the read length exactly (`M/I/S` sums).

use crate::extend::{
    extend_hint, push_op, CigarOp, DiagHint, ExtendBuffer, ExtendParams, Extension,
};
use crate::fasta::Reference;
use crate::gtf::{Junction, JunctionLib, RefinedJunction, SpliceSignal};
use crate::index::Index;
use crate::intron_chain::cigar_conserving;
use crate::seed::{minimizers, SeedParams, Strand};
use crate::span::cigar_ref_end;
use crate::splice::{
    junction_signal, lib_window_candidates, refine_junction, splice_signal, splice_support_pass,
    tail_gate, SpliceParams, SPLIT_COST, SPLIT_LIB_WINDOW,
};

/// Library-driven tail lower bound.
pub const TAIL_MIN: usize = 7;

/// Tail-start retry radius (±5bp).
pub const TAIL_RETRY_RADIUS: i64 = 5;

/// Direct-scan channel tail length range (5..19bp; tails ≥ 19 bp reach the
/// seed channel, which needs k+w-1 = 19 bp for one minimizer).
pub const DIRECT_TAIL_MIN: usize = 5;
pub const DIRECT_TAIL_MAX: usize = 19;

/// Direct-scan donor search radius (±12bp).
pub const DIRECT_DONOR_RADIUS: i64 = 12;

/// Direct-scan downstream search span.
pub const DIRECT_SCAN_SPAN: u32 = 50_000;

/// Exact prefix probe length cap.
pub const PROBE_MAX: usize = 10;

/// Seed channel minimum tail (≥ k = 15).
pub const SEED_TAIL_MIN: usize = 15;

/// Seed channel: first N tail minimizers.
pub const SEED_CHANNEL_MINS: usize = 8;

/// Seed channel occupancy cap (`0 < count ≤ 100`).
pub const SEED_CHANNEL_OCC_CAP: u32 = 100;

/// Seed channel positions cap (first 64).
pub const SEED_CHANNEL_POS_CAP: usize = 64;

/// Extra reference margin past the tail length in tail windows.
const TAIL_MARGIN: u32 = 30;

/// Successful split rescue: replacement placement pieces.
#[derive(Clone, Debug)]
pub struct SplitRescue {
    /// New genomic start of the whole stitched alignment.
    pub pos: u32,
    /// Conserving CIGAR over the whole read.
    pub cigar: Vec<CigarOp>,
    /// The junction the tail was placed across.
    pub junction: RefinedJunction,
    /// Rescue score (channel A retry total; tail extension score otherwise).
    pub score: i32,
}

/// Inputs describing the primary alignment of the read being rescued.
pub struct SplitContext<'a> {
    /// Reference genome.
    pub reference: &'a Reference,
    /// Minimizer index (needed by the seed channel; `None` disables it).
    pub index: Option<&'a Index>,
    /// Junction library (empty ⇒ de-novo channels enabled).
    pub lib: &'a JunctionLib,
    /// The read in alignment orientation (anchor `qpos` space).
    pub read: &'a [u8],
    /// Primary alignment contig.
    pub contig: u32,
    /// Primary alignment strand.
    pub strand: Strand,
    /// Primary genomic start.
    pub pos: u32,
    /// Primary genomic end (M/D/N advance), exclusive.
    pub ref_end: u32,
    /// First aligned read base of the primary alignment.
    pub read_start: u32,
    /// One past the last aligned read base of the primary alignment.
    pub read_end: u32,
    /// Primary CIGAR (whole read).
    pub cigar: &'a [CigarOp],
    /// Seed parameters (seed channel).
    pub seed_params: SeedParams,
    /// Splice parameters.
    pub splice_params: SpliceParams,
    /// Extension parameters.
    pub extend_params: ExtendParams,
}

impl SplitContext<'_> {
    fn minus(&self) -> bool {
        self.strand == Strand::Minus
    }

    /// Right-tail library candidates: junction starts near `ref_end` (plus)
    /// or junction ends near `pos` (minus), zipped with support.
    fn right_lib_candidates(&self) -> Vec<(&Junction, u32)> {
        let lo = self.ref_end.saturating_sub(SPLIT_LIB_WINDOW);
        let hi = self.ref_end + SPLIT_LIB_WINDOW;
        lib_window_candidates(self.lib, self.minus(), self.contig, lo, hi)
    }

    /// Left-tail library candidates: junction ends near `pos` (plus) or
    /// junction starts near `ref_end` (minus).
    fn left_lib_candidates(&self) -> Vec<(&Junction, u32)> {
        let lo = self.pos.saturating_sub(SPLIT_LIB_WINDOW);
        let hi = self.pos + SPLIT_LIB_WINDOW;
        lib_window_candidates(self.lib, !self.minus(), self.contig, lo, hi)
    }
}

/// One `retry_tail_start` outcome.
#[derive(Clone, Debug)]
pub struct RetryOutcome {
    /// Winning tail-start offset (split point in read coordinates).
    pub dro: u32,
    /// Left-window genomic start (for `new_pos` reconstruction).
    pub left_win_lo: u32,
    /// Left segment extension (`read[..dro]`).
    pub left: Extension,
    /// Right segment extension (`read[dro..]`).
    pub right: Extension,
    /// `le.score + re.score + signal_score + lib_bonus − SPLIT_COST`.
    pub total: i32,
}

/// Retry the tail start ±5bp around `tail_start` against one junction.
///
/// Per `dro`: the left segment extends in the exon window ending at `j.start`,
/// the right segment in the window starting at `j.end`; each
/// segment must align ≥ `len × 4/5`. The best total wins (ties keep the
/// lowest `dro`).
#[allow(clippy::too_many_arguments)]
pub fn retry_tail_start(
    reference: &Reference,
    read: &[u8],
    tail_start: u32,
    j: &Junction,
    signal: SpliceSignal,
    support: u32,
    params: &SpliceParams,
    extend_params: &ExtendParams,
    buf: &mut ExtendBuffer,
) -> Option<RetryOutcome> {
    let ctg = reference.contigs.get(j.contig as usize)?;
    let read_len = read.len() as u32;
    let lo_dro = tail_start.saturating_sub(TAIL_RETRY_RADIUS as u32).max(1);
    let hi_dro = (tail_start + TAIL_RETRY_RADIUS as u32).min(read_len.saturating_sub(1));
    let signal_score = signal.score(params);
    let lib_bonus = (support.min(10) * 5) as i32;

    let mut best: Option<RetryOutcome> = None;
    for dro in lo_dro..=hi_dro {
        let lseg = &read[..dro as usize];
        let rseg = &read[dro as usize..];
        let llen = lseg.len() as u32;
        let rlen = rseg.len() as u32;

        // Oriented-read geometry is strand-symmetric: the left segment
        // ends at the donor side (j.start), the right segment starts at the
        // acceptor side (j.end).
        let l_hi = j.start.min(ctg.len);
        let l_lo = l_hi.saturating_sub(llen + TAIL_MARGIN);
        let r_lo = j.end;
        let r_hi =
            ((r_lo as u64) + (rlen as u64) + (TAIL_MARGIN as u64)).min(ctg.len as u64) as u32;
        if l_hi <= l_lo || r_hi <= r_lo {
            continue;
        }
        let lwin = ctg.slice_ascii(l_lo, l_hi);
        let rwin = ctg.slice_ascii(r_lo, r_hi);
        let hint = DiagHint {
            offset: 0,
            num: 1,
            den: 1,
        };
        let le = extend_hint(lseg, &lwin, extend_params, buf, hint);
        let re = extend_hint(rseg, &rwin, extend_params, buf, hint);
        let l_aligned = le.read_end.saturating_sub(le.read_start) as usize;
        let r_aligned = re.read_end.saturating_sub(re.read_start) as usize;
        if l_aligned * 5 < lseg.len() * 4 || r_aligned * 5 < rseg.len() * 4 {
            continue;
        }
        // Local extensions may clip either end or stop short of the exon
        // boundary. Only the two OUTER read ends may remain unaligned.
        // Filter each retry before ranking so an invalid high score cannot
        // hide a lower-scoring, correctly anchored split.
        if le.read_end != llen
            || re.read_start != 0
            || re.ref_start != 0
            || cigar_ref_end(l_lo + le.ref_start, &le.cigar) != j.start
        {
            continue;
        }
        let total = le.score + re.score + signal_score + lib_bonus - SPLIT_COST;
        if best.as_ref().is_none_or(|b| total > b.total) {
            best = Some(RetryOutcome {
                dro,
                left_win_lo: l_lo,
                left: le,
                right: re,
                total,
            });
        }
    }
    best
}

/// Channel-A stitch: `[left ext] N(intron) [right ext]`, fully independent of
/// the primary CIGAR; `new_pos` = left segment reference start.
pub fn build_rescue_ext(outcome: &RetryOutcome, refined: &RefinedJunction) -> SplitRescue {
    let intron = refined.junction.end - refined.junction.start;
    let mut cigar = outcome.left.cigar.clone();
    push_op(&mut cigar, CigarOp::RefSkip(intron));
    for op in &outcome.right.cigar {
        push_op(&mut cigar, *op);
    }
    SplitRescue {
        pos: outcome.left_win_lo + outcome.left.ref_start,
        cigar,
        junction: refined.clone(),
        score: outcome.total,
    }
}

/// Right-tail stitch against the primary body: body CIGAR truncated at
/// `body_read_end` read depth (dropping the trailing soft clip) + `N` + the
/// tail extension as a whole. Conserves the read length exactly.
pub fn build_rescue(
    body_cigar: &[CigarOp],
    body_read_end: u32,
    tail_ext: &Extension,
    intron_len: u32,
) -> Vec<CigarOp> {
    let mut cigar = cigar_before_read_depth(body_cigar, body_read_end);
    push_op(&mut cigar, CigarOp::RefSkip(intron_len));
    for op in &tail_ext.cigar {
        push_op(&mut cigar, *op);
    }
    cigar
}

/// Left-tail stitch: the tail extension as a whole + `N` + body CIGAR from
/// `body_read_start` read depth (dropping the leading soft clip). Conserves
/// the read length exactly.
pub fn build_rescue_left(
    tail_ext: &Extension,
    body_cigar: &[CigarOp],
    body_read_start: u32,
    intron_len: u32,
) -> Vec<CigarOp> {
    let mut cigar: Vec<CigarOp> = Vec::new();
    for op in &tail_ext.cigar {
        push_op(&mut cigar, *op);
    }
    push_op(&mut cigar, CigarOp::RefSkip(intron_len));
    for op in cigar_from_read_depth(body_cigar, body_read_start) {
        push_op(&mut cigar, op);
    }
    cigar
}

/// Validate the assembled query span, terminal-only clipping and actual
/// reference coordinates of the newly introduced junction. No CIGAR edits
/// are made here; a failed candidate leaves the caller's primary intact.
pub(crate) fn rescue_cigar_valid(
    pos: u32,
    cigar: &[CigarOp],
    read_len: usize,
    junction: &Junction,
) -> bool {
    if !cigar_conserving(cigar, read_len) {
        return false;
    }
    let mut reference = pos as u64;
    let mut found = false;
    for (i, &op) in cigar.iter().enumerate() {
        let n = match op {
            CigarOp::Match(n)
            | CigarOp::Ins(n)
            | CigarOp::Del(n)
            | CigarOp::RefSkip(n)
            | CigarOp::SoftClip(n) => n,
        };
        if n == 0 {
            return false;
        }
        match op {
            CigarOp::SoftClip(_) if i != 0 && i + 1 != cigar.len() => return false,
            CigarOp::RefSkip(n) => {
                if i == 0 || i + 1 == cigar.len() {
                    return false;
                }
                if reference == junction.start as u64 && reference + n as u64 == junction.end as u64
                {
                    found = true;
                }
                reference += n as u64;
            }
            CigarOp::Match(n) | CigarOp::Del(n) => reference += n as u64,
            _ => {}
        }
    }
    found && reference <= u32::MAX as u64
}

/// Join a locally extended tail only when both sides meet the same query
/// split and the declared genomic boundaries. Indels and outer clips remain
/// valid; reference-only ops at the body split are accounted for explicitly.
fn anchored_tail_rescue(
    ctx: &SplitContext,
    right: bool,
    body_edge: u32,
    win_lo: u32,
    ext: &Extension,
    junction: &Junction,
) -> Option<(Vec<CigarOp>, u32)> {
    let intron_len = junction.end.checked_sub(junction.start)?;
    let (cigar, pos) = if right {
        if ext.read_start != 0 || win_lo + ext.ref_start != junction.end {
            return None;
        }
        (build_rescue(ctx.cigar, body_edge, ext, intron_len), ctx.pos)
    } else {
        let body_prefix = cigar_before_read_depth(ctx.cigar, body_edge);
        if ext.read_end != body_edge || cigar_ref_end(ctx.pos, &body_prefix) != junction.end {
            return None;
        }
        (
            build_rescue_left(ext, ctx.cigar, body_edge, intron_len),
            win_lo + ext.ref_start,
        )
    };
    rescue_cigar_valid(pos, &cigar, ctx.read.len(), junction).then_some((cigar, pos))
}

/// Prefix of `cigar` consuming at most `depth` read bases (mid-op split; ops
/// past the depth, including trailing reference-only ops, are dropped).
pub fn cigar_before_read_depth(cigar: &[CigarOp], depth: u32) -> Vec<CigarOp> {
    let mut out = Vec::new();
    let mut used = 0u32;
    'outer: for &op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => {
                if used >= depth {
                    break 'outer;
                }
                let room = depth - used;
                if n <= room {
                    push_op(&mut out, op);
                    used += n;
                } else {
                    push_op(&mut out, same_kind(op, room));
                    used = depth;
                }
            }
            CigarOp::Del(_) | CigarOp::RefSkip(_) => {
                if used < depth {
                    push_op(&mut out, op);
                }
            }
        }
    }
    out
}

/// Suffix of `cigar` from `depth` read bases on (leading op split if needed).
pub fn cigar_from_read_depth(cigar: &[CigarOp], depth: u32) -> Vec<CigarOp> {
    let mut out = Vec::new();
    let mut used = 0u32;
    for &op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Ins(n) | CigarOp::SoftClip(n) => {
                let end = used + n;
                if end <= depth {
                    used = end;
                    continue;
                }
                if used < depth {
                    let rest = end - depth;
                    push_op(&mut out, same_kind(op, rest));
                } else {
                    push_op(&mut out, op);
                }
                used = end;
            }
            CigarOp::Del(_) | CigarOp::RefSkip(_) => {
                if used >= depth {
                    push_op(&mut out, op);
                }
            }
        }
    }
    out
}

fn same_kind(op: CigarOp, len: u32) -> CigarOp {
    match op {
        CigarOp::Match(_) => CigarOp::Match(len),
        CigarOp::Ins(_) => CigarOp::Ins(len),
        CigarOp::Del(_) => CigarOp::Del(len),
        CigarOp::RefSkip(_) => CigarOp::RefSkip(len),
        CigarOp::SoftClip(_) => CigarOp::SoftClip(len),
    }
}

pub(crate) fn extend_tail(
    ctg: &crate::fasta::Contig,
    tail: &[u8],
    win_lo: u32,
    win_hi: u32,
    params: &ExtendParams,
    buf: &mut ExtendBuffer,
) -> Option<(u32, Extension)> {
    if win_hi <= win_lo {
        return None;
    }
    let win = ctg.slice_ascii(win_lo, win_hi);
    let ext = extend_hint(
        tail,
        &win,
        params,
        buf,
        DiagHint {
            offset: 0,
            num: 1,
            den: 1,
        },
    );
    if !tail_gate(&ext, tail.len()) {
        return None;
    }
    Some((win_lo, ext))
}

/// Channel A (library-driven): best retry over the ±20bp window junctions.
fn channel_a(ctx: &SplitContext, right: bool, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    let cands = if right {
        ctx.right_lib_candidates()
    } else {
        ctx.left_lib_candidates()
    };
    let ctg = ctx.reference.contigs.get(ctx.contig as usize)?;
    let mut best: Option<(RetryOutcome, RefinedJunction)> = None;
    for (j, support) in cands {
        if j.minus_strand != ctx.minus() {
            continue;
        }
        let intron_len = j.end - j.start;
        if intron_len < ctx.splice_params.min_intron || intron_len > ctx.splice_params.max_intron {
            continue;
        }
        let signal = junction_signal(ctg, j);
        if signal == SpliceSignal::NonCanonical && support < 2 {
            continue;
        }
        let tail_start = if right { ctx.read_end } else { ctx.read_start };
        let outcome = match retry_tail_start(
            ctx.reference,
            ctx.read,
            tail_start,
            j,
            signal,
            support,
            &ctx.splice_params,
            &ctx.extend_params,
            buf,
        ) {
            Some(o) => o,
            None => continue,
        };
        // Keep the highest total; ties keep the first candidate (sorted).
        if best.as_ref().is_none_or(|(b, _)| outcome.total > b.total) {
            best = Some((
                outcome,
                RefinedJunction {
                    junction: *j,
                    signal,
                    known_support: support,
                },
            ));
        }
    }
    let (outcome, refined) = best?;
    let rescue = build_rescue_ext(&outcome, &refined);
    if !rescue_cigar_valid(rescue.pos, &rescue.cigar, ctx.read.len(), &refined.junction) {
        return None;
    }
    Some(rescue)
}

/// Dinucleotide ending at `p` (i.e. `[p-2, p)`), read from a decoded window
/// starting at `lo`; `NN` when out of range.
const GT_SET: [[u8; 2]; 3] = [*b"GT", *b"GC", *b"AT"];
const AG_AC_SET: [[u8; 2]; 2] = [*b"AG", *b"AC"];
const CT_GT_SET: [[u8; 2]; 2] = [*b"CT", *b"GT"];
const AC_GC_AT_SET: [[u8; 2]; 3] = [*b"AC", *b"GC", *b"AT"];

fn dinuc_before(win: &[u8], lo: u32, p: u32) -> [u8; 2] {
    if p < lo + 2 {
        return *b"NN";
    }
    let i = (p - lo) as usize;
    if i > win.len() {
        return *b"NN";
    }
    [
        win[i - 2].to_ascii_uppercase(),
        win[i - 1].to_ascii_uppercase(),
    ]
}

fn dinuc_after(win: &[u8], lo: u32, p: u32) -> [u8; 2] {
    if p < lo {
        return *b"NN";
    }
    let i = (p - lo) as usize;
    if i + 2 > win.len() {
        return *b"NN";
    }
    [win[i].to_ascii_uppercase(), win[i + 1].to_ascii_uppercase()]
}

fn matches_at(win: &[u8], lo: u32, p: u32, pat: &[u8]) -> bool {
    if p < lo {
        return false;
    }
    let i = (p - lo) as usize;
    if i + pat.len() > win.len() {
        return false;
    }
    win[i..i + pat.len()]
        .iter()
        .zip(pat.iter())
        .all(|(r, q)| r.eq_ignore_ascii_case(q))
}

/// Channel B (direct scan, empty library only).
///
/// Right tail: donors ±12 around `ref_end`, probe scan ≤50kb downstream for
/// an exact tail-prefix match ending at a plausible acceptor. Left tail:
/// acceptors ±12 around `pos`, probe scan ≤50kb upstream for an exact
/// tail-suffix match at a plausible donor. Candidate score =
/// `sig.score*100 - |body_boundary_deviation| - |intron - 5000|/1000`; the
/// best candidate's whole tail is then extension-validated (≥ 4/5, score > 0).
/// The discovered junction is recorded with the detected splice signal, support 0.
fn channel_b(ctx: &SplitContext, right: bool, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    let ctg = ctx.reference.contigs.get(ctx.contig as usize)?;
    let tl = if right {
        ctx.read.len() - ctx.read_end as usize
    } else {
        ctx.read_start as usize
    };
    if !(DIRECT_TAIL_MIN..DIRECT_TAIL_MAX).contains(&tl) {
        return None;
    }
    let probe_len = tl.min(PROBE_MAX);
    let min_intron = ctx.splice_params.min_intron;
    let span = DIRECT_SCAN_SPAN as u64;

    // Transcript-oriented dinucleotide sets: body-side boundary vs tail-side
    // boundary of the intron (mirrored on minus-strand junctions).
    let (body_set, tail_set): (&[[u8; 2]], &[[u8; 2]]) = match (right, ctx.minus()) {
        (true, false) => (&GT_SET, &AG_AC_SET),
        (true, true) => (&CT_GT_SET, &AC_GC_AT_SET),
        (false, false) => (&AG_AC_SET, &GT_SET),
        (false, true) => (&AC_GC_AT_SET, &CT_GT_SET),
    };

    // One decoded window covering the ±12 body scan and the full search span.
    let (scan_lo, scan_hi) = if right {
        let lo = ctx.ref_end.saturating_sub(DIRECT_DONOR_RADIUS as u32);
        let hi = (ctx.ref_end as u64 + DIRECT_DONOR_RADIUS as u64).min(ctg.len as u64) as u32;
        (lo, ((hi as u64) + span + 2).min(ctg.len as u64) as u32)
    } else {
        let hi = (ctx.pos as u64 + DIRECT_DONOR_RADIUS as u64 + 1).min(ctg.len as u64) as u32;
        let lo = ctx
            .pos
            .saturating_sub(DIRECT_DONOR_RADIUS as u32)
            .saturating_sub(span as u32 + 2);
        (lo, hi)
    };
    if scan_hi <= scan_lo || (scan_hi - scan_lo) as u64 <= min_intron as u64 {
        return None;
    }
    let win = ctg.slice_ascii(scan_lo, scan_hi);
    let body_anchor = if right { ctx.ref_end } else { ctx.pos };

    // Body-overshoot note: a continuous extension can run a few bp INTO the
    // intron (repeat similarity). The true junction boundary then sits `over`
    // bp before the body edge, and the exon fragment on the tail side starts
    // `over` bp earlier in the read. Each boundary candidate computes its own
    // overshoot and uses an overshoot-adjusted probe and tail.
    let mut candidates = Vec::new(); // (score, intron_start, intron_end, over)
    if right {
        let d_lo = ctx.ref_end.saturating_sub(DIRECT_DONOR_RADIUS as u32);
        let d_hi = (ctx.ref_end + DIRECT_DONOR_RADIUS as u32).min(ctg.len.saturating_sub(2));
        for d in d_lo..=d_hi {
            if !body_set.contains(&dinuc_after(&win, scan_lo, d)) {
                continue;
            }
            let dev = (d as i64 - body_anchor as i64).abs() as i32;
            let over = (ctx.ref_end as i64 - d as i64).max(0) as u32;
            if over > ctx.read_end {
                continue;
            }
            let probe_start = (ctx.read_end - over) as usize;
            let probe_d = &ctx.read[probe_start..probe_start + probe_len];
            let q_lo = d as u64 + min_intron as u64;
            let q_hi = (d as u64 + span).min(scan_hi as u64);
            for q in q_lo..=q_hi {
                let qu = q as u32;
                if tail_set.contains(&dinuc_before(&win, scan_lo, qu))
                    && matches_at(&win, scan_lo, qu, probe_d)
                {
                    let intron = qu as i64 - d as i64;
                    let score = 4 * 100 - dev - ((intron - 5000).abs() / 1000) as i32;
                    candidates.push((score, d, qu, over));
                }
            }
        }
    } else {
        let a_lo = ctx.pos.saturating_sub(DIRECT_DONOR_RADIUS as u32);
        let a_hi = (ctx.pos + DIRECT_DONOR_RADIUS as u32).min(ctg.len);
        for a in a_lo..=a_hi {
            if a < min_intron + 2 || !body_set.contains(&dinuc_before(&win, scan_lo, a)) {
                continue;
            }
            let dev = (a as i64 - body_anchor as i64).abs() as i32;
            let over = (a as i64 - ctx.pos as i64).max(0) as u32;
            if ctx.read_start + over > ctx.read.len() as u32 {
                continue;
            }
            let probe_end = (ctx.read_start + over) as usize;
            let probe_d = &ctx.read[probe_end - probe_len..probe_end];
            let s_lo = ((a as i64 - DIRECT_SCAN_SPAN as i64)
                .max(scan_lo as i64)
                .max(0)) as u64;
            let s_hi = a as u64 - min_intron as u64;
            for s in s_lo..=s_hi {
                if s < probe_len as u64 {
                    continue;
                }
                let su = s as u32;
                if tail_set.contains(&dinuc_after(&win, scan_lo, su))
                    && matches_at(&win, scan_lo, su - probe_len as u32, probe_d)
                {
                    let intron = a as i64 - s as i64;
                    let score = 4 * 100 - dev - ((intron - 5000).abs() / 1000) as i32;
                    candidates.push((score, su, a, over));
                }
            }
        }
    }
    // Stable ordering preserves the previous tie policy. Try the next
    // candidate when the highest scan score cannot form an anchored splice.
    candidates.sort_by_key(|c| std::cmp::Reverse(c.0));
    for (scan_score, start, end, over) in candidates {
        // Overshoot-adjusted tail and body edge: the body ran `over` bp into the
        // intron, so the true exon fragment is `over` bp longer on the tail side.
        let (ext_tail, body_edge) = if right {
            let edge = ctx.read_end - over;
            (&ctx.read[edge as usize..], edge)
        } else {
            let edge = ctx.read_start + over;
            (&ctx.read[..edge as usize], edge)
        };
        let etl = ext_tail.len();

        // Whole-tail extension validation.
        let (win_lo, win_hi) = if right {
            let hi =
                ((end as u64) + (etl as u64) + (TAIL_MARGIN as u64)).min(ctg.len as u64) as u32;
            (end, hi)
        } else {
            let hi = start;
            let lo = hi.saturating_sub(etl as u32 + TAIL_MARGIN);
            (lo, hi)
        };
        let Some((_, ext)) = extend_tail(ctg, ext_tail, win_lo, win_hi, &ctx.extend_params, buf)
        else {
            continue;
        };
        let signal = splice_signal(
            dinuc_after(&win, scan_lo, start),
            dinuc_before(&win, scan_lo, end),
            ctx.minus(),
        );
        let junction = RefinedJunction {
            junction: Junction {
                contig: ctx.contig,
                start,
                end,
                minus_strand: ctx.minus(),
            },
            signal,
            known_support: 0,
        };
        let Some((cigar, pos)) =
            anchored_tail_rescue(ctx, right, body_edge, win_lo, &ext, &junction.junction)
        else {
            continue;
        };
        return Some(SplitRescue {
            pos,
            cigar,
            junction,
            score: ext.score + scan_score,
        });
    }
    None
}

/// Channel C (seed-driven fallback): seed the tail, map its first
/// minimizers through the index, keep hits on the same contig/strand falling
/// in the min/max intron window across the body edge, refine with
/// `refine_junction` + splice support gates, extension-validate. First success (deterministic
/// candidate order) wins.
fn channel_c(ctx: &SplitContext, right: bool, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    let index = ctx.index?;
    let ctg = ctx.reference.contigs.get(ctx.contig as usize)?;
    let (tail, body_edge_read) = if right {
        (&ctx.read[ctx.read_end as usize..], ctx.read_end)
    } else {
        (&ctx.read[..ctx.read_start as usize], ctx.read_start)
    };
    let tl = tail.len();
    if tl < SEED_TAIL_MIN || tl < ctx.seed_params.k as usize {
        return None;
    }
    let k = ctx.seed_params.k;
    let mins: Vec<_> = minimizers(tail, ctx.seed_params)
        .into_iter()
        .take(SEED_CHANNEL_MINS)
        .collect();

    // Naive (start, end) intron candidates from seed placements.
    let mut cands: Vec<(u32, u32)> = Vec::new();
    for m in &mins {
        let hit = index.query(m.kmer);
        if hit.count == 0 || hit.count > SEED_CHANNEL_OCC_CAP {
            continue;
        }
        for &p in hit.positions.iter().take(SEED_CHANNEL_POS_CAP) {
            let (contig2, rpos, rstrand) = crate::index::unpack_pos(p);
            if contig2 != ctx.contig {
                continue;
            }
            let astrand = Strand::xor(m.strand, rstrand);
            if astrand != ctx.strand {
                continue;
            }
            let qpos = if astrand == Strand::Plus {
                m.pos
            } else {
                tl as u32 - m.pos - k
            };
            let g = rpos as i64 - qpos as i64; // genomic start of the tail
            let cand = if right {
                (ctx.ref_end as i64, g)
            } else {
                (g + tl as i64, ctx.pos as i64)
            };
            if cand.1 <= cand.0 {
                continue;
            }
            let len = cand.1 - cand.0;
            if len < ctx.splice_params.min_intron as i64
                || len > ctx.splice_params.max_intron as i64
            {
                continue;
            }
            cands.push((cand.0 as u32, cand.1 as u32));
        }
    }
    cands.sort_unstable();
    cands.dedup();

    for (s0, e0) in cands {
        let Some(refined) = refine_junction(
            ctx.reference,
            ctx.contig,
            s0,
            e0,
            ctx.minus(),
            ctx.lib,
            &ctx.splice_params,
        ) else {
            continue;
        };
        let intron_len = refined.junction.end - refined.junction.start;
        if !splice_support_pass(intron_len, refined.signal, refined.known_support) {
            continue;
        }
        let (win_lo, win_hi) = if right {
            let lo = refined.junction.end;
            let hi = ((lo as u64) + (tl as u64) + (TAIL_MARGIN as u64)).min(ctg.len as u64) as u32;
            (lo, hi)
        } else {
            let hi = refined.junction.start;
            let lo = hi.saturating_sub(tl as u32 + TAIL_MARGIN);
            (lo, hi)
        };
        let Some((win_lo, ext)) = extend_tail(ctg, tail, win_lo, win_hi, &ctx.extend_params, buf)
        else {
            continue;
        };
        let Some((cigar, pos)) =
            anchored_tail_rescue(ctx, right, body_edge_read, win_lo, &ext, &refined.junction)
        else {
            continue;
        };
        return Some(SplitRescue {
            pos,
            cigar,
            junction: refined,
            score: ext.score,
        });
    }
    None
}

/// Rescue a 3' (right) soft-clipped tail across a splice junction.
///
/// Try library channel A first, then direct scan (B) and tail seeding (C).
/// Tails shorter than [`TAIL_MIN`] are not rescued.
pub fn rescue_right_tail(ctx: &SplitContext, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    let tail_len = ctx.read.len() - ctx.read_end as usize;
    if tail_len < TAIL_MIN {
        return None;
    }
    if !ctx.lib.is_empty() {
        // Library first (high-confidence); fall back to de-novo scan/seed
        // (B/C) for junctions NOT in the library (novel splices).
        if let Some(a) = channel_a(ctx, true, buf) {
            return Some(a);
        }
    }
    channel_b(ctx, true, buf).or_else(|| channel_c(ctx, true, buf))
}

/// Rescue a 5' (left) soft-clipped tail across a splice junction (mirror of
/// [`rescue_right_tail`]; the alignment start moves left on success).
pub fn rescue_left_tail(ctx: &SplitContext, buf: &mut ExtendBuffer) -> Option<SplitRescue> {
    let tail_len = ctx.read_start as usize;
    if tail_len < TAIL_MIN {
        return None;
    }
    if !ctx.lib.is_empty() {
        if let Some(a) = channel_a(ctx, false, buf) {
            return Some(a);
        }
    }
    channel_b(ctx, false, buf).or_else(|| channel_c(ctx, false, buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn retry(seq: &str, read: &str, split: u32) -> Option<RetryOutcome> {
        let reference =
            crate::fasta::parse_fasta_bytes(format!(">synthetic\n{seq}\n").as_bytes()).unwrap();
        retry_tail_start(
            &reference,
            read.as_bytes(),
            split,
            &Junction {
                contig: 0,
                start: 60,
                end: 160,
                minus_strand: false,
            },
            SpliceSignal::GtAg,
            2,
            &SpliceParams::default(),
            &ExtendParams::default(),
            &mut ExtendBuffer::default(),
        )
    }

    #[test]
    fn retry_rejects_soft_clipping_at_the_internal_seam() {
        let seq = format!(
            "{}{}GT{}AG{}",
            "A".repeat(20),
            "C".repeat(40),
            "T".repeat(96),
            "G".repeat(60)
        );
        let read = format!("{}AAAA{}", "C".repeat(40), "G".repeat(40));
        let outcome = retry(&seq, &read, 44);
        assert!(
            outcome.is_none(),
            "independent local alignments leave unaligned bases at the splice: {outcome:?}"
        );
    }

    #[test]
    fn retry_rejects_reference_gap_before_donor() {
        let seq = format!(
            "{}{}AAAAGT{}AG{}",
            "A".repeat(20),
            "C".repeat(36),
            "T".repeat(96),
            "G".repeat(60)
        );
        let read = format!("{}{}", "C".repeat(40), "G".repeat(40));
        let outcome = retry(&seq, &read, 40);
        assert!(
            outcome.is_none(),
            "left extension stops before donor: {outcome:?}"
        );
    }

    #[test]
    fn retry_rejects_reference_gap_after_acceptor() {
        let seq = format!(
            "{}{}GT{}AGAAAA{}",
            "A".repeat(20),
            "C".repeat(40),
            "T".repeat(96),
            "G".repeat(60)
        );
        let read = format!("{}{}", "C".repeat(40), "G".repeat(40));
        let outcome = retry(&seq, &read, 40);
        assert!(
            outcome.is_none(),
            "right extension starts past acceptor: {outcome:?}"
        );
    }

    #[test]
    fn retry_keeps_exact_splice_and_outer_soft_clips() {
        let seq = format!(
            "{}{}GT{}AG{}",
            "A".repeat(20),
            "C".repeat(40),
            "T".repeat(96),
            "G".repeat(60)
        );
        for (prefix, suffix) in [("", ""), ("TTTT", "AAAA")] {
            let read = format!("{prefix}{}{suffix}", "C".repeat(40) + &"G".repeat(40));
            let o = retry(&seq, &read, prefix.len() as u32 + 40).expect("valid splice");
            assert_eq!(o.left.read_end, o.dro);
            assert_eq!(o.right.read_start, 0);
            assert_eq!(o.right.ref_start, 0);
        }
    }

    #[test]
    fn tail_stitch_preserves_indels_outer_clips_and_checks_both_sides() {
        use CigarOp::{Del as D, Ins as I, Match as M, SoftClip as S};
        let reference = crate::fasta::parse_fasta_bytes(
            format!(">synthetic\n{}\n", "A".repeat(300)).as_bytes(),
        )
        .unwrap();
        let lib = JunctionLib::default();
        let read = vec![b'A'; 84];
        let j = Junction {
            contig: 0,
            start: 60,
            end: 160,
            minus_strand: false,
        };
        let right_body = vec![S(2), M(18), I(2), M(20), S(42)];
        let left_body = vec![S(42), M(18), I(2), M(20), S(2)];
        let mut ctx = SplitContext {
            reference: &reference,
            index: None,
            lib: &lib,
            read: &read,
            contig: 0,
            strand: Strand::Plus,
            pos: 22,
            ref_end: 60,
            read_start: 2,
            read_end: 42,
            cigar: &right_body,
            seed_params: SeedParams { k: 15, w: 5 },
            splice_params: SpliceParams::default(),
            extend_params: ExtendParams::default(),
        };
        // Insertion in the body and deletion in the tail: query and reference
        // lengths differ, but both touch the correct splice boundaries.
        let mut ext = Extension {
            read_start: 0,
            read_end: 40,
            ref_start: 0,
            cigar: vec![M(18), D(2), M(22), S(2)],
            score: 80,
        };
        let (cigar, pos) = anchored_tail_rescue(&ctx, true, 42, 160, &ext, &j).unwrap();
        assert_eq!(pos, 22);
        assert_eq!(cigar_ref_end(pos, &cigar), 202);
        ext.ref_start = 1;
        assert!(anchored_tail_rescue(&ctx, true, 42, 160, &ext, &j).is_none());
        ext.ref_start = 0;
        ctx.pos = 21; // body ends at 59, not donor 60
        assert!(anchored_tail_rescue(&ctx, true, 42, 160, &ext, &j).is_none());

        ctx.cigar = &left_body;
        ctx.pos = 160;
        ctx.ref_end = 198;
        ctx.read_start = 42;
        ctx.read_end = 82;
        ext = Extension {
            read_start: 2,
            read_end: 42,
            ref_start: 0,
            cigar: vec![S(2), M(18), D(2), M(22)],
            score: 80,
        };
        for strand in [Strand::Plus, Strand::Minus] {
            ctx.strand = strand;
            let mut junction = j;
            junction.minus_strand = strand == Strand::Minus;
            let (cigar, pos) = anchored_tail_rescue(&ctx, false, 42, 18, &ext, &junction).unwrap();
            assert_eq!(pos, 18);
            assert_eq!(cigar_ref_end(pos, &cigar), 198);
        }
        ctx.pos = 161; // suffix would start one base past acceptor
        assert!(anchored_tail_rescue(&ctx, false, 42, 18, &ext, &j).is_none());
        ctx.pos = 160;
        ext.read_end = 40;
        ext.cigar = vec![S(2), M(18), D(2), M(20), S(2)];
        assert!(anchored_tail_rescue(&ctx, false, 42, 20, &ext, &j).is_none());
    }
}
