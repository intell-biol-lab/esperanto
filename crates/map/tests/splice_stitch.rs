use esperanto_map::extend::CigarOp;
use esperanto_map::gtf::{Junction, RefinedJunction, SpliceSignal};
use esperanto_map::splice::stitch_cigar;

fn junction(start: u32, end: u32, minus_strand: bool) -> RefinedJunction {
    RefinedJunction {
        junction: Junction {
            contig: 0,
            start,
            end,
            minus_strand,
        },
        signal: SpliceSignal::GtAg,
        known_support: 1,
    }
}

fn spans(cigar: &[CigarOp]) -> (u32, u32) {
    let mut query = 0;
    let mut reference = 0;
    for operation in cigar {
        match operation {
            CigarOp::Match(length) => {
                assert!(*length > 0);
                query += length;
                reference += length;
            }
            CigarOp::Ins(length) | CigarOp::SoftClip(length) => {
                assert!(*length > 0);
                query += length;
            }
            CigarOp::Del(length) | CigarOp::RefSkip(length) => {
                assert!(*length > 0);
                reference += length;
            }
        }
    }
    (query, reference)
}

#[test]
fn exon_boundary_resumption_never_emits_zero_length_operations() {
    use CigarOp::{Del, Ins, Match, RefSkip, SoftClip};
    let bounds = [(0, 100, 110), (10, 210, 230)];
    let cases = [
        (
            vec![Match(10), Ins(1), Match(10)],
            vec![Match(10), Ins(1), RefSkip(100), Match(10)],
        ),
        (
            vec![Match(10), Del(2), Match(8)],
            vec![Match(10), RefSkip(100), Del(2), Match(8)],
        ),
        (
            vec![Match(8), Del(2), Match(10)],
            vec![Match(8), Del(2), RefSkip(100), Match(10)],
        ),
        (
            vec![Match(8), Del(2), Ins(1), Match(10)],
            vec![Match(8), Del(2), Ins(1), RefSkip(100), Match(10)],
        ),
        (
            vec![SoftClip(3), Match(10), Ins(1), Match(10), SoftClip(2)],
            vec![
                SoftClip(3),
                Match(10),
                Ins(1),
                RefSkip(100),
                Match(10),
                SoftClip(2),
            ],
        ),
    ];
    for minus_strand in [false, true] {
        let junctions = [junction(110, 210, minus_strand)];
        for (input, expected) in &cases {
            let observed = stitch_cigar(input, 0, &bounds, &junctions);
            assert_eq!(observed, *expected, "{input:?}, minus={minus_strand}");
            let (query, reference) = spans(input);
            assert_eq!(spans(&observed), (query, reference + 100));
        }
    }
}

#[test]
fn multiple_exon_boundaries_preserve_indels_and_clips() {
    use CigarOp::{Del, Ins, Match, RefSkip, SoftClip};
    let input = [
        SoftClip(3),
        Match(10),
        Ins(1),
        Match(10),
        Del(2),
        Match(8),
        SoftClip(2),
    ];
    let bounds = [(0, 100, 110), (10, 210, 220), (20, 420, 440)];
    let junctions = [junction(110, 210, false), junction(220, 420, false)];
    let observed = stitch_cigar(&input, 0, &bounds, &junctions);
    assert_eq!(
        observed,
        [
            SoftClip(3),
            Match(10),
            Ins(1),
            RefSkip(100),
            Match(10),
            RefSkip(200),
            Del(2),
            Match(8),
            SoftClip(2)
        ]
    );
    let (query, reference) = spans(&input);
    assert_eq!(spans(&observed), (query, reference + 300));
}

#[test]
fn boundary_offsets_and_non_crossing_paths_are_preserved() {
    use CigarOp::{Del, Ins, Match, RefSkip};
    let bounds = [(0, 100, 110), (10, 210, 230)];
    let junctions = [junction(110, 210, false)];
    let cases = [
        (0, vec![Match(20)], vec![Match(10), RefSkip(100), Match(10)]),
        (
            8,
            vec![Del(4), Match(6)],
            vec![Del(2), RefSkip(100), Del(2), Match(6)],
        ),
        (0, vec![Match(10), Ins(1)], vec![Match(10), Ins(1)]),
        (10, vec![Ins(1), Match(10)], vec![Ins(1), Match(10)]),
        (11, vec![Match(5)], vec![Match(5)]),
        (0, vec![Match(8)], vec![Match(8)]),
    ];
    for (start, input, expected) in cases {
        let observed = stitch_cigar(&input, start, &bounds, &junctions);
        assert_eq!(observed, expected);
        spans(&observed);
    }
    let input = [Match(5), Ins(1), Match(5)];
    assert_eq!(stitch_cigar(&input, 0, &[], &[]), input);
}
