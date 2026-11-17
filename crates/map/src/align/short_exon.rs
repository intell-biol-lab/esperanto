//! Bounded supplementary library-supported short-chain candidates.
use super::*;

/// Keep the ordinary 40-point chain gate unchanged. A 25..39-point exon
/// can propose a splice only when its endpoint is near the junction library.
const MIN_SHORT_CHAIN: i32 = 25;

impl Aligner<'_> {
    pub(super) fn short_exon_placements(
        &mut self,
        seq: &[u8],
        mins: &[Minimizer],
    ) -> (Vec<Placement>, Option<crate::index::AnchorBudget>) {
        let Some(lib) = self.jlib.clone().filter(|lib| !lib.is_empty()) else {
            return (Vec::new(), None);
        };
        let mut anchors = std::mem::take(&mut self.scratch.anchors);
        let budget = crate::index::collect_anchors_budgeted_into(
            self.index,
            mins,
            seq.len() as u32,
            self.config.seed.k,
            self.index.freq_cutoff,
            MID_ROUND_MAX_ANCHORS,
            &mut anchors,
        );
        let mut chains = Vec::new();
        let mut params = self.config.chain;
        params.min_chain_score = MIN_SHORT_CHAIN;
        chain_anchors_into(&mut anchors, &params, &mut chains, &mut self.scratch.chain);
        self.scratch.anchors = anchors;
        chains.retain(|c| {
            if c.score >= self.config.chain.min_chain_score {
                return false;
            }
            let (start, end) = c.ref_span(self.config.seed.k);
            let radius = crate::splice::SPLIT_LIB_WINDOW;
            !lib.range_start(
                c.contig,
                end.saturating_sub(radius),
                end.saturating_add(radius),
            )
            .0
            .is_empty()
                || !lib
                    .range_end(
                        c.contig,
                        start.saturating_sub(radius),
                        start.saturating_add(radius),
                    )
                    .is_empty()
        });
        let rc = self.take_rc(seq);
        let mut remaining = MAX_TIED_CANDIDATES;
        let mut out = Vec::new();
        for chain in chains.iter().take(MAX_TIED_CANDIDATES) {
            let mut bodies = Vec::new();
            self.extend_unspliced_chain(seq, chain, false, &rc, &mut remaining, &mut bodies);
            for body in bodies {
                let before = ReadAlignment {
                    contig: body.contig,
                    pos: body.pos,
                    strand: body.strand,
                    cigar: body.cigar,
                    score: body.ext_score,
                    chain_score: body.chain_score,
                    ..Default::default()
                };
                let proposal =
                    self.propose_split_rescue_with_library_policy(seq, before.clone(), true);
                let rescued = self.select_rescue(seq, before.clone(), proposal);
                if rescued.junctions.is_empty()
                    || rescued.junctions.iter().any(|j| !lib.contains(&j.junction))
                    || !self.short_exon_tail_beats_background(seq, &before, &rescued)
                {
                    continue;
                }
                out.push(Placement::from_genomic(
                    rescued.contig,
                    rescued.pos,
                    rescued.strand,
                    rescued.cigar,
                    rescued.junctions,
                    rescued.chain_score,
                    rescued.score,
                ));
            }
        }
        self.scratch.rc = rc;
        (out, Some(budget))
    }

    /// Library support is not evidence that a constant terminal run belongs
    /// to the exon. Compare sequence gain alone with the existing coordinate-
    /// free background on newly covered query ends. This conservatively also
    /// refuses biological tails with identical sequence evidence.
    fn short_exon_tail_beats_background(
        &mut self,
        seq: &[u8],
        before: &ReadAlignment,
        after: &ReadAlignment,
    ) -> bool {
        if !self.config.context_background {
            return true;
        }
        let clips = |cigar: &[CigarOp]| {
            let left = match cigar.first() {
                Some(CigarOp::SoftClip(n)) => *n as usize,
                _ => 0,
            };
            let right = match cigar.last() {
                Some(CigarOp::SoftClip(n)) => *n as usize,
                _ => 0,
            };
            (left, right)
        };
        let (old_left, old_right) = clips(&before.cigar);
        let (new_left, new_right) = clips(&after.cigar);
        let rc = self.take_rc(seq);
        let query = if before.strand == Strand::Minus {
            &rc
        } else {
            seq
        };
        let mut background = 0;
        if new_left < old_left {
            background += crate::repeat_context::background_score(
                &query[new_left..old_left],
                self.config.extend,
            );
        }
        if new_right < old_right {
            background += crate::repeat_context::background_score(
                &query[query.len() - old_right..query.len() - new_right],
                self.config.extend,
            );
        }
        let mut placement = |a: &ReadAlignment| {
            let mut p = Placement::from_genomic(
                a.contig,
                a.pos,
                a.strand,
                a.cigar.clone(),
                a.junctions.clone(),
                a.chain_score,
                a.score,
            );
            self.rescore_strict(query, &mut p);
            p.strict.score - p.strict.splice_score
        };
        let sequence_gain = placement(after) - placement(before);
        self.scratch.rc = rc;
        sequence_gain > background
    }
}
