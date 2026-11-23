//! Test-only dev44 arithmetic oracle; never compiled into production.
use super::*;

fn parent_log_likelihood(observations: &[&Observation], ref_base: u8, alt: u8, p: f64) -> f64 {
    observations
        .iter()
        .map(|obs| {
            let error = 10_f64.powf(-(obs.bq as f64) / 10.0);
            let (zero, one) = if obs.base == ref_base {
                (1.0 - error, error / 3.0)
            } else if obs.base == alt {
                (error / 3.0, 1.0 - error)
            } else {
                (error / 3.0, error / 3.0)
            };
            ((1.0 - p) * zero + p * one).max(1e-300).ln()
        })
        .sum()
}

fn parent_likelihood_derivative(
    observations: &[&Observation],
    ref_base: u8,
    alt: u8,
    p: f64,
) -> f64 {
    observations
        .iter()
        .map(|obs| {
            let error = 10_f64.powf(-(obs.bq as f64) / 10.0);
            let (zero, one) = if obs.base == ref_base {
                (1.0 - error, error / 3.0)
            } else if obs.base == alt {
                (error / 3.0, 1.0 - error)
            } else {
                (error / 3.0, error / 3.0)
            };
            (one - zero) / ((1.0 - p) * zero + p * one).max(1e-300)
        })
        .sum()
}

fn parent_likelihood_curvature(
    observations: &[&Observation],
    ref_base: u8,
    alt: u8,
    p: f64,
) -> f64 {
    observations
        .iter()
        .map(|obs| {
            let error = 10_f64.powf(-(obs.bq as f64) / 10.0);
            let (zero, one) = if obs.base == ref_base {
                (1.0 - error, error / 3.0)
            } else if obs.base == alt {
                (error / 3.0, 1.0 - error)
            } else {
                (error / 3.0, error / 3.0)
            };
            let denominator = ((1.0 - p) * zero + p * one).max(1e-300);
            ((one - zero) / denominator).powi(2)
        })
        .sum()
}

fn parent_log_uniform_likelihood(observations: &[&Observation], ref_base: u8, alt: u8) -> f64 {
    GL32_X
        .iter()
        .zip(GL32_W.iter())
        .map(|(&x, w)| parent_log_likelihood(observations, ref_base, alt, x) + w.ln())
        .fold(f64::NEG_INFINITY, log_add_exp)
}

fn parent_likelihood_features(observations: &[&Observation], ref_base: u8, alt: u8) -> (f64, f64) {
    if observations.is_empty() {
        return (0.0, 0.0);
    }
    let log_zero = parent_log_likelihood(observations, ref_base, alt, 0.0);
    let log_one = parent_log_likelihood(observations, ref_base, alt, 1.0);
    if observations.len() <= 63 {
        let log_uniform = parent_log_uniform_likelihood(observations, ref_base, alt);
        return (log_uniform - log_zero, log_uniform - log_one);
    }
    let derivative_zero = parent_likelihood_derivative(observations, ref_base, alt, 0.0);
    let derivative_one = parent_likelihood_derivative(observations, ref_base, alt, 1.0);
    let mode = if derivative_zero <= 0.0 {
        0.0
    } else if derivative_one >= 0.0 {
        1.0
    } else {
        let mut lo = 0.0;
        let mut hi = 1.0;
        for _ in 0..100 {
            let mid = (lo + hi) / 2.0;
            if parent_likelihood_derivative(observations, ref_base, alt, mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        (lo + hi) / 2.0
    };
    let peak = parent_log_likelihood(observations, ref_base, alt, mode);
    let scale = if mode == 0.0 || mode == 1.0 {
        let derivative = if mode == 0.0 {
            derivative_zero.abs()
        } else {
            derivative_one.abs()
        };
        // The opposite endpoint can have an arbitrarily steep slope. Near a
        // flat boundary, curvature rather than slope controls the mass width.
        1.0 / derivative
            .max(parent_likelihood_curvature(observations, ref_base, alt, mode).sqrt())
            .max(1.0)
    } else {
        1.0 / parent_likelihood_curvature(observations, ref_base, alt, mode)
            .max(1.0)
            .sqrt()
    };
    let mut points = vec![0.0, 1.0];
    for multiple in [1.0, 4.0, 16.0, 64.0] {
        for point in [mode - multiple * scale, mode + multiple * scale] {
            if point > 0.0 && point < 1.0 {
                points.push(point);
            }
        }
    }
    points.sort_by(f64::total_cmp);
    points.dedup_by(|a, b| a == b);
    let integral = points
        .windows(2)
        .map(|window| {
            integrate(
                |p| (parent_log_likelihood(observations, ref_base, alt, p) - peak).exp(),
                window[0],
                window[1],
            )
        })
        .sum::<f64>();
    (
        peak + integral.ln() - log_zero,
        peak + integral.ln() - log_one,
    )
}

fn owned(groups: &[(u8, u8, usize)]) -> Vec<Observation> {
    groups
        .iter()
        .flat_map(|&(base, bq, n)| {
            std::iter::repeat_n(
                Observation {
                    base,
                    bq,
                    mapq: 60.0,
                    end_fraction: 0.5,
                    orientation: 0,
                    softclip_fraction: 0.0,
                    gap_nearness: 0.0,
                },
                n,
            )
        })
        .collect()
}

#[test]
fn cached_likelihood_preserves_parent_f64_bits() {
    let mut cases = vec![
        vec![],
        vec![(b'A', 30, 9997), (b'G', 30, 3)],
        vec![(b'A', 30, 3), (b'G', 30, 9997)],
        vec![(b'A', 30, 2997), (b'G', 30, 1)],
        vec![(b'C', 93, 60), (b'G', 93, 3)],
        vec![(b'C', 93, 10000)],
        vec![(b'A', 30, 5000), (b'G', 30, 5000)],
    ];
    for bq in [0, 1, 13, 20, 30, 60, 93, 254, 255] {
        for depth in [1, 63, 64, 127] {
            for base in BASES {
                cases.push(vec![(base, bq, depth)]);
            }
            cases.push(vec![(b'A', bq, depth), (b'G', bq, 3), (b'C', bq, 2)]);
        }
    }
    for groups in cases {
        let mut obs = owned(&groups);
        for reversed in [false, true] {
            if reversed {
                obs.reverse();
            }
            let refs: Vec<_> = obs.iter().collect();
            for (reference, alt) in [(b'A', b'G'), (b'T', b'C')] {
                let expected = parent_likelihood_features(&refs, reference, alt);
                let actual = likelihood_features(&refs, reference, alt);
                for (a, b) in [(actual.0, expected.0), (actual.1, expected.1)] {
                    assert!(
                        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
                        "groups={groups:?} reversed={reversed}: {actual:?} != {expected:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn probabilities_are_prepared_once_per_observation() {
    let obs = owned(&[(b'A', 30, 93), (b'G', 30, 7)]);
    let refs: Vec<_> = obs.iter().collect();
    LIKELIHOOD_WORK.with(|work| work.set((0, 0)));
    let _ = likelihood_features(&refs, b'A', b'G');
    let (probabilities, _) = LIKELIHOOD_WORK.with(|work| work.get());
    assert_eq!(
        probabilities,
        obs.len(),
        "BQ probabilities repeated during integration/search"
    );
}

#[test]
fn stagnant_f64_search_does_not_repeat_derivative_evaluations() {
    let obs = owned(&[(b'A', 30, 93), (b'G', 30, 7)]);
    let refs: Vec<_> = obs.iter().collect();
    LIKELIHOOD_WORK.with(|work| work.set((0, 0)));
    let expected = parent_likelihood_features(&refs, b'A', b'G');
    let actual = likelihood_features(&refs, b'A', b'G');
    assert_eq!(
        (actual.0.to_bits(), actual.1.to_bits()),
        (expected.0.to_bits(), expected.1.to_bits())
    );
    let (_, derivatives) = LIKELIHOOD_WORK.with(|work| work.get());
    assert!(
        derivatives < 70,
        "evaluated derivative {derivatives} times after f64 bounds stopped moving"
    );
}

#[test]
fn cached_endpoints_preserve_every_u8_quality_and_mixed_order() {
    let mut obs = Vec::new();
    for bq in 0..=u8::MAX {
        for base in BASES {
            obs.extend(owned(&[(base, bq, 1)]));
        }
    }
    // Unequal class counts and varying BQ exercise nontrivial derivative sums.
    obs.extend(owned(&[(b'A', 13, 57), (b'G', 93, 19)]));
    for reversed in [false, true] {
        if reversed {
            obs.reverse();
        }
        let refs: Vec<_> = obs.iter().collect();
        for (reference, alt) in [(b'A', b'G'), (b'T', b'C')] {
            let cached = Likelihood::new(&refs, reference, alt);
            for p in GL32_X.into_iter().chain([0.0, 0.5, 1.0]) {
                for (a, b) in [
                    (
                        cached.log(p),
                        parent_log_likelihood(&refs, reference, alt, p),
                    ),
                    (
                        cached.derivative(p),
                        parent_likelihood_derivative(&refs, reference, alt, p),
                    ),
                    (
                        cached.curvature(p),
                        parent_likelihood_curvature(&refs, reference, alt, p),
                    ),
                ] {
                    assert!(
                        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
                        "p={p}, reversed={reversed}: {a} != {b}"
                    );
                }
            }
            let actual = likelihood_features(&refs, reference, alt);
            let expected = parent_likelihood_features(&refs, reference, alt);
            assert_eq!(
                (actual.0.to_bits(), actual.1.to_bits()),
                (expected.0.to_bits(), expected.1.to_bits())
            );
        }
    }
}
