use super::*;

fn observation(bq: u8) -> Observation {
    Observation {
        base: b'A',
        bq,
        mapq: 60.0,
        end_fraction: 0.5,
        orientation: 0,
        softclip_fraction: 0.0,
        gap_nearness: 0.0,
    }
}

#[test]
fn center_observations_follow_stable_template_ids_for_every_map_seed() {
    for reversed in [false, true] {
        for _ in 0..8 {
            let mut ids: Vec<u32> = (0..64).collect();
            if reversed {
                ids.reverse();
            }
            let entries: HashMap<_, _> = ids
                .into_iter()
                .map(|i| (i, observation(i as u8 + 13)))
                .collect();
            let window = HashMap::from([(100, entries)]);
            let bqs: Vec<_> = center_observations(&window, 100)
                .iter()
                .map(|o| o.bq)
                .collect();
            assert_eq!(
                bqs,
                (13..77).collect::<Vec<u8>>(),
                "random map order reached floating reductions"
            );
            assert!(center_observations(&window, 101).is_empty());
        }
    }
}

#[test]
fn ordered_center_keeps_template_conflict_merging_and_contents() {
    let mut entries = HashMap::new();
    entries.insert(9, observation(30));
    let mut alternative = observation(30);
    alternative.base = b'G';
    merge_observation(entries.get_mut(&9).unwrap(), alternative);
    entries.insert(2, observation(60));
    let window = HashMap::from([(100, entries)]);
    let refs = center_observations(&window, 100);
    assert_eq!(refs.len(), 2);
    assert_eq!(refs[0].bq, 60);
    assert_eq!(refs[1].bq, 30);
    assert_eq!(refs[1].base, b'N');
}
