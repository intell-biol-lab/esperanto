//! Generated specificity controls; no human reads.
use crate::{
    align::{revcomp, AlignConfig, Aligner},
    fasta::parse_fasta_bytes,
    index::Index,
    seed::SeedParams,
};

fn random(state: &mut u64, n: usize) -> Vec<u8> {
    (0..n)
        .map(|_| {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            b"ACGT"[(*state & 3) as usize]
        })
        .collect()
}

pub(crate) fn fixture(full_target: bool) -> (Index, Vec<u8>) {
    let mut state = 80817;
    let mut body = random(&mut state, 42);
    body[40..].copy_from_slice(b"GT");
    let query = [body, vec![b'A'; 76]].concat();
    let mut homolog = random(&mut state, 42);
    homolog[40..].copy_from_slice(b"GT");
    let mut tail = vec![b'A'; 76];
    // One high-frequency complex/homopolymer boundary kmer per copy.
    // No A15 exists in the homolog; no count cap is manually changed.
    for p in (13..76).step_by(8) {
        tail[p] = b'C';
    }
    homolog.extend(tail);
    let repeats = [vec![b'N'; 32], homolog, vec![b'N'; 32]]
        .concat()
        .repeat(6000);
    let background = random(&mut state, 1_000_000);
    let mut fasta = format!(
        ">repeats\n{}\n>background\n{}\n",
        std::str::from_utf8(&repeats).unwrap(),
        std::str::from_utf8(&background).unwrap()
    );
    if full_target {
        let target = [
            random(&mut state, 120),
            query.clone(),
            random(&mut state, 120),
        ]
        .concat();
        fasta.push_str(&format!(
            ">full_target\n{}\n",
            std::str::from_utf8(&target).unwrap()
        ));
    }
    let index = Index::build(
        parse_fasta_bytes(fasta.as_bytes()).unwrap(),
        SeedParams { k: 15, w: 5 },
    );
    assert!(index.freq_cutoff < 5000);
    (index, query)
}

fn incomplete_refusal(minus: bool) {
    let (index, seq) = fixture(false);
    let seq = if minus { revcomp(&seq) } else { seq };
    let mut off = AlignConfig::rna_default();
    off.repeat_context = false;
    assert!(
        Aligner::new(&index, off).align_read(&seq).is_none(),
        "control must reproduce ordinary candidate omission"
    );
    let actual = Aligner::new(&index, AlignConfig::rna_default()).align_read(&seq);
    assert!(
        actual.is_none(),
        "partial low-complexity-only support must not fabricate primary coverage: {actual:?}"
    );
}

#[test]
fn context_background_rejects_partial_low_complexity_plus() {
    incomplete_refusal(false);
}

#[test]
fn context_background_rejects_partial_low_complexity_minus() {
    incomplete_refusal(true);
}

#[test]
fn context_background_real_full_sequence_with_same_low_complexity_tail_is_retained() {
    let (index, query) = fixture(true);
    for seq in [query.clone(), revcomp(&query)] {
        let actual = Aligner::new(&index, AlignConfig::rna_default())
            .align_read(&seq)
            .expect("do not blacklist an entire low-complexity read");
        assert_eq!((actual.contig, actual.pos, actual.score), (2, 120, 236));
    }
}

#[test]
fn context_background_ablation_and_empty_cache_keep_truthful_limited_evidence() {
    let (index, query) = fixture(false);
    for seq in [query.clone(), revcomp(&query)] {
        let mut off = AlignConfig::rna_default();
        off.context_background = false;
        let (old, old_cache) = Aligner::new(&index, off).align_read_full(&seq);
        let old = old.unwrap();
        assert_eq!(old.score, 108);
        assert_eq!(old.origin.label(), "repeat_context");
        assert_eq!(old.mapq_override, Some(0));
        let b0 = old_cache
            .unwrap()
            .seed_budget
            .unwrap()
            .repeat_context
            .unwrap();
        assert!(b0.background.is_none());
        let (new, cache) = Aligner::new(&index, AlignConfig::rna_default()).align_read_full(&seq);
        assert!(new.is_none());
        let cache = cache.unwrap();
        assert!(cache.placements.is_empty() && cache.seed_search_limited());
        let b1 = cache.seed_budget.unwrap().repeat_context.unwrap();
        let bg = b1.background.unwrap();
        assert!(bg.score >= 152);
        assert!(bg.rejected_candidates > 0);
        assert_eq!(b1.positive_candidates, 0);
        assert_eq!(b0.scanned_postings, b1.scanned_postings);
        assert_eq!(b0.eligible_postings, b1.eligible_postings);
        let row =
            serde_json::to_value(cache.diagnostic_pool(&index, seq.len(), "independent")).unwrap();
        assert_eq!(
            row["seed_budget"]["repeat_context"]["background"]["score"],
            bg.score
        );
        let mut observed = None;
        let fallback = Aligner::new(&index, AlignConfig::rna_default())
            .align_read_with_mate_observed(&seq, 0, 72, Some(&cache), Some(&mut observed));
        if let Some(mut a) = fallback {
            a.mapq_override = Some(60);
            assert_eq!(
                a.mapq_of_checked(),
                0,
                "old mate hypotheses cannot restore confidence"
            );
        }
        let pool = serde_json::to_value(observed.unwrap()).unwrap();
        assert_eq!(
            pool["seed_budget"]["repeat_context"]["background"]["score"],
            bg.score
        );
    }
}

#[test]
fn context_background_full_low_complexity_context_candidates_remain_ambiguous() {
    let mut state = 41081;
    let query = [random(&mut state, 42), vec![b'A'; 76]].concat();
    let repeated = [vec![b'N'; 32], query.clone(), vec![b'N'; 32]]
        .concat()
        .repeat(6000);
    let background = random(&mut state, 1_000_000);
    let index = Index::build(
        parse_fasta_bytes(
            &[
                b">repeats\n".as_slice(),
                &repeated,
                b"\n>background\n",
                &background,
                b"\n",
            ]
            .concat(),
        )
        .unwrap(),
        SeedParams { k: 15, w: 5 },
    );
    assert!(index.freq_cutoff < 5000);
    for seq in [query.clone(), revcomp(&query)] {
        let (a, cache) = Aligner::new(&index, AlignConfig::rna_default()).align_read_full(&seq);
        let a = a.expect("full matching context evidence must beat its background");
        assert_eq!(a.origin.label(), "repeat_context");
        assert_eq!((a.score, a.second_score), (236, 236));
        assert_eq!(a.cigar, vec![crate::extend::CigarOp::Match(118)]);
        assert!(a.n_cands > 1);
        assert_eq!(a.mapq_override, Some(0));
        let b = cache.unwrap().seed_budget.unwrap().repeat_context.unwrap();
        assert!(b.background.unwrap().score < a.score);
        assert!(b.positive_candidates > 8 && b.retained_candidates == 8);
    }
}
