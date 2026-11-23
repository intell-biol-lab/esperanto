//! Native reference-aware molecular evidence scoring.
//!
//! This module implements the frozen 27-feature LightGBM contract used by the
//! unified evidence model. It accepts a model directory containing `model.txt`,
//! `freeze.json`, and `alu.bed`, and computes scores directly from
//! a coordinate-sorted indexed BAM.

use anyhow::{bail, ensure, Context, Result};
use rayon::prelude::*;
use rust_htslib::bam::{record::Cigar, IndexedReader, Read as _};
use std::collections::HashMap;
use std::path::Path;

pub const FEATURE_COUNT: usize = 27;
pub const MIN_MAPQ: u8 = 20;
pub const MIN_BQ: u8 = 13;
pub const RADIUS: i64 = 200;

const BASES: [u8; 4] = *b"ACGT";
const EXCLUDE_FLAGS: u16 = 0x4 | 0x100 | 0x800;

#[derive(Clone, Debug)]
struct Observation {
    base: u8,
    bq: u8,
    mapq: f64,
    end_fraction: f64,
    orientation: i8,
    softclip_fraction: f64,
    gap_nearness: f64,
}

type TemplateId = u32;
type Window = HashMap<i64, HashMap<TemplateId, Observation>>;

#[cfg(test)]
#[path = "reader_reuse_tests.rs"]
mod reader_reuse_tests;

#[cfg(test)]
#[path = "likelihood_cache_tests.rs"]
mod likelihood_cache_tests;

#[cfg(test)]
#[path = "deterministic_order_tests.rs"]
mod deterministic_order_tests;

#[cfg(test)]
#[path = "template_merge_tests.rs"]
mod template_merge_tests;

#[cfg(test)]
thread_local! {
    static LIKELIHOOD_WORK: std::cell::Cell<(usize,usize)> = const { std::cell::Cell::new((0,0)) };
}

#[derive(Clone, Debug)]
struct IntervalIndex {
    starts: Vec<i64>,
    ends: Vec<i64>,
}

#[derive(Clone, Debug, Default)]
pub struct AluIndex {
    intervals: HashMap<String, IntervalIndex>,
}

impl AluIndex {
    pub fn from_bed(path: &Path) -> Result<Self> {
        let mut raw: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for (line_no, line) in std::fs::read_to_string(path)?.lines().enumerate() {
            if line.trim().is_empty() || line.starts_with('#') || line.starts_with("track") {
                continue;
            }
            let fields: Vec<_> = line.split_whitespace().collect();
            ensure!(fields.len() >= 3, "invalid Alu BED row {}", line_no + 1);
            raw.entry(fields[0].to_owned()).or_default().push((
                fields[1]
                    .parse()
                    .with_context(|| format!("Alu start at line {}", line_no + 1))?,
                fields[2]
                    .parse()
                    .with_context(|| format!("Alu end at line {}", line_no + 1))?,
            ));
        }
        let intervals = raw
            .into_iter()
            .map(|(chrom, mut values)| {
                values.sort_unstable();
                let starts = values.iter().map(|x| x.0).collect();
                let mut ends = Vec::with_capacity(values.len());
                for (_, end) in values {
                    ends.push(ends.last().copied().unwrap_or(0).max(end));
                }
                (chrom, IntervalIndex { starts, ends })
            })
            .collect();
        Ok(Self { intervals })
    }

    pub fn contains(&self, chrom: &str, pos0: i64) -> bool {
        let Some(index) = self.intervals.get(chrom) else {
            return false;
        };
        let n = index.starts.partition_point(|&start| start <= pos0);
        n > 0 && index.ends[n - 1] > pos0
    }
}

#[derive(Clone, Debug)]
struct Tree {
    feature: Vec<usize>,
    threshold: Vec<f64>,
    decision: Vec<u8>,
    left: Vec<i32>,
    right: Vec<i32>,
    leaf: Vec<f64>,
}

impl Tree {
    fn score(&self, features: &[f32; FEATURE_COUNT]) -> f64 {
        let mut node = 0_i32;
        while node >= 0 {
            let i = node as usize;
            let value = features[self.feature[i]];
            // LightGBM encodes missing-value handling in decision_type. In
            // addition to NaN, type 1 treats an exact zero as missing (the
            // serialized values are 8/10 rather than plain 2).
            let missing_type = self.decision[i] >> 2;
            let missing = value.is_nan() || (missing_type == 1 && value == 0.0);
            let goes_left = if missing {
                self.decision[i] & 2 != 0
            } else {
                f64::from(value) <= self.threshold[i]
            };
            node = if goes_left {
                self.left[i]
            } else {
                self.right[i]
            };
        }
        self.leaf[(-node - 1) as usize]
    }
}

#[derive(Clone, Debug)]
pub struct Model {
    trees: Vec<Tree>,
}

#[derive(Clone, Debug)]
struct ObservationSite {
    chrom: String,
    pos: i64,
    ref_base: u8,
    counts: [u32; 4],
    depth: u32,
}

fn parse_values<T: std::str::FromStr>(block: &str, key: &str) -> Result<Vec<T>>
where
    T::Err: std::fmt::Debug,
{
    let prefix = format!("{key}=");
    block
        .lines()
        .find_map(|line| line.strip_prefix(&prefix))
        .with_context(|| format!("missing LightGBM field {key}"))?
        .split_whitespace()
        .map(|value| {
            value
                .parse()
                .map_err(|e| anyhow::anyhow!("invalid {key}: {e:?}"))
        })
        .collect()
}

impl Model {
    pub fn load(model_dir: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(model_dir.join("model.txt"))?;
        let mut trees = Vec::new();
        for block in text.split("\nTree=").skip(1) {
            if !block.contains("\nsplit_feature=") {
                break;
            }
            let tree = Tree {
                feature: parse_values(block, "split_feature")?,
                threshold: parse_values(block, "threshold")?,
                decision: parse_values(block, "decision_type")?,
                left: parse_values(block, "left_child")?,
                right: parse_values(block, "right_child")?,
                leaf: parse_values(block, "leaf_value")?,
            };
            ensure!(tree.feature.len() == tree.threshold.len());
            ensure!(tree.feature.len() == tree.decision.len());
            ensure!(tree.feature.len() == tree.left.len());
            ensure!(tree.feature.len() == tree.right.len());
            ensure!(tree.leaf.len() == tree.feature.len() + 1);
            ensure!(tree.feature.iter().all(|&x| x < FEATURE_COUNT));
            trees.push(tree);
        }
        ensure!(
            trees.len() == 400,
            "expected 400 LightGBM trees, found {}",
            trees.len()
        );
        Ok(Self { trees })
    }

    pub fn score(&self, features: &[f32; FEATURE_COUNT]) -> f64 {
        let raw: f64 = self.trees.iter().map(|tree| tree.score(features)).sum();
        1.0 / (1.0 + (-raw).exp())
    }
}

fn template_name(qname: &[u8]) -> &[u8] {
    if qname.len() >= 2
        && qname[qname.len() - 2] == b'/'
        && matches!(qname[qname.len() - 1], b'1' | b'2')
    {
        &qname[..qname.len() - 2]
    } else {
        qname
    }
}

fn base_index(base: u8) -> Option<usize> {
    BASES.iter().position(|&value| value == base)
}

fn merge_observation(old: &mut Observation, new: Observation) {
    if new.bq > old.bq {
        *old = new;
    } else if new.bq == old.bq {
        // Reduce all maximum-BQ evidence, not only discordant bases. These
        // fieldwise summaries are a composite template observation, not a
        // selected physical read. N remains absorbing at the same BQ.
        if old.base != new.base {
            old.base = b'N';
        }
        old.mapq = if old.mapq.is_finite() {
            old.mapq.max(new.mapq)
        } else {
            new.mapq
        };
        old.end_fraction = old.end_fraction.max(new.end_fraction);
        old.orientation = if old.orientation == new.orientation {
            old.orientation
        } else {
            -1
        };
        old.softclip_fraction = old.softclip_fraction.min(new.softclip_fraction);
        old.gap_nearness = old.gap_nearness.min(new.gap_nearness);
    }
}

fn read_observations(
    rec: &rust_htslib::bam::Record,
    left: i64,
    right: i64,
) -> Vec<(i64, Observation)> {
    if rec.flags() & EXCLUDE_FLAGS != 0 || rec.mapq() < MIN_MAPQ || rec.pos() < 0 {
        return Vec::new();
    }
    let seq = rec.seq().as_bytes();
    let qual = rec.qual();
    if seq.is_empty() {
        return Vec::new();
    }
    let length = seq.len();
    let soft = rec
        .cigar()
        .iter()
        .filter_map(|c| match c {
            Cigar::SoftClip(n) => Some(*n as f64),
            _ => None,
        })
        .sum::<f64>()
        / length.max(1) as f64;
    let mut gaps = Vec::new();
    let mut q = 0_i64;
    for op in rec.cigar().iter() {
        let n = op.len() as i64;
        match op {
            Cigar::Ins(_) => {
                gaps.push(q as f64 - 0.5);
                gaps.push((q + n) as f64 - 0.5);
                q += n;
            }
            Cigar::Del(_) | Cigar::RefSkip(_) => {
                gaps.push(q as f64 - 0.5);
            }
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) | Cigar::SoftClip(_) => q += n,
            _ => {}
        }
    }
    let orientation = if rec.is_reverse() ^ rec.is_last_in_template() {
        1
    } else {
        0
    };
    let mapq = if rec.mapq() == 255 {
        f64::NAN
    } else {
        f64::from(rec.mapq())
    };
    let mut out = Vec::new();
    let mut qpos = 0_i64;
    let mut rpos = rec.pos();
    for op in rec.cigar().iter() {
        let n = op.len() as i64;
        match op {
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                let lo = rpos.max(left);
                let hi = (rpos + n).min(right);
                for rp in lo..hi {
                    let qp = (qpos + rp - rpos) as usize;
                    let bq = qual.get(qp).copied().unwrap_or(255);
                    if bq == 255 || bq < MIN_BQ {
                        continue;
                    }
                    let Some(&raw_base) = seq.get(qp) else {
                        continue;
                    };
                    let base = raw_base.to_ascii_uppercase();
                    if base_index(base).is_none() {
                        continue;
                    }
                    let end_fraction =
                        (qp.min(length - 1 - qp) as f64) / ((length - 1).max(1) as f64 / 2.0);
                    let gap = gaps
                        .iter()
                        .map(|g| ((qp as f64) - g).abs())
                        .fold(f64::INFINITY, f64::min);
                    out.push((
                        rp,
                        Observation {
                            base,
                            bq,
                            mapq,
                            end_fraction,
                            orientation,
                            softclip_fraction: soft,
                            gap_nearness: if gap.is_finite() {
                                1.0 / (1.0 + gap)
                            } else {
                                0.0
                            },
                        },
                    ));
                }
                qpos += n;
                rpos += n;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) => qpos += n,
            Cigar::Del(_) | Cigar::RefSkip(_) => rpos += n,
            _ => {}
        }
    }
    out
}

fn merge_template<S: std::hash::BuildHasher>(
    observations: &mut HashMap<TemplateId, Observation, S>,
    id: TemplateId,
    observation: Observation,
) {
    match observations.entry(id) {
        std::collections::hash_map::Entry::Occupied(mut entry) => {
            merge_observation(entry.get_mut(), observation);
        }
        std::collections::hash_map::Entry::Vacant(entry) => {
            entry.insert(observation);
        }
    }
}

fn read_window(bam: &mut IndexedReader, chrom: &str, left: i64, right: i64) -> Result<Window> {
    let tid = bam
        .header()
        .tid(chrom.as_bytes())
        .context("contig absent from BAM")?;
    bam.fetch((tid, left, right))?;
    let mut template_ids: HashMap<Vec<u8>, TemplateId> = HashMap::new();
    let mut positions: Window = HashMap::new();
    for record in bam.records() {
        let record = record?;
        let evidence = read_observations(&record, left, right);
        if evidence.is_empty() {
            continue;
        }
        let mut name = record
            .aux(b"RG")
            .ok()
            .and_then(|aux| match aux {
                rust_htslib::bam::record::Aux::String(value) => Some(value.as_bytes().to_vec()),
                _ => None,
            })
            .unwrap_or_default();
        name.push(0);
        name.extend_from_slice(template_name(record.qname()));
        let next =
            TemplateId::try_from(template_ids.len()).context("too many templates in window")?;
        let id = *template_ids.entry(name).or_insert(next);
        for (pos, observation) in evidence {
            let observations = positions.entry(pos).or_default();
            merge_template(observations, id, observation);
        }
    }
    // Keep conflicting intermediate observations until every record has merged:
    // a later higher-quality observation can resolve an earlier conflict.
    for observations in positions.values_mut() {
        observations.retain(|_, observation| base_index(observation.base).is_some());
    }
    positions.retain(|_, observations| !observations.is_empty());
    Ok(positions)
}
/// Independent informative A/C/G/T template counts using the frozen molecular
/// observation reader. None denotes a contig absent from the BAM dictionary.
/// Windows are bounded to 4096 reference bases; output follows input order.
pub fn count_sites(bam_path: &Path, sites: &[(String, i64)]) -> Result<Vec<Option<[u64; 4]>>> {
    let mut bam = IndexedReader::from_path(bam_path)?;
    let mut order: Vec<_> = sites.iter().enumerate().collect();
    for (_, (_, pos)) in &order {
        ensure!(*pos >= 1, "recount positions must be 1-based and positive");
    }
    order.sort_by(|a, b| a.1.cmp(b.1));
    let mut result = vec![None; sites.len()];
    let mut first = 0;
    while first < order.len() {
        let (chrom, pos) = order[first].1;
        let mut end = first + 1;
        while end < order.len() && order[end].1 .0 == *chrom && order[end].1 .1 - pos < 4096 {
            end += 1;
        }
        if let Some(tid) = bam.header().tid(chrom.as_bytes()) {
            let right = order[end - 1].1 .1;
            ensure!(
                right as u64
                    <= bam
                        .header()
                        .target_len(tid)
                        .context("missing BAM contig length")?,
                "recount position outside BAM contig"
            );
            let window = read_window(&mut bam, chrom, pos - 1, right)?;
            for &(i, (_, p)) in &order[first..end] {
                let mut counts = [0u64; 4];
                if let Some(observations) = window.get(&(p - 1)) {
                    for observation in observations.values() {
                        if let Some(b) = base_index(observation.base) {
                            counts[b] += 1;
                        }
                    }
                }
                result[i] = Some(counts);
            }
        }
        first = end;
    }
    Ok(result)
}

fn mean(values: impl Iterator<Item = f64>) -> f64 {
    let values: Vec<_> = values.filter(|value| value.is_finite()).collect();
    if values.is_empty() {
        f64::NAN
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

/// Observation-specific mixture endpoints, prepared once and visited in the
/// original order for every likelihood/derivative/curvature sum. Do not group
/// equal observations: regrouping would change floating-point addition order.
struct Likelihood {
    probabilities: Vec<(f64, f64)>,
}

impl Likelihood {
    fn new(observations: &[&Observation], ref_base: u8, alt: u8) -> Self {
        let probabilities = observations
            .iter()
            .map(|obs| {
                #[cfg(test)]
                LIKELIHOOD_WORK.with(|work| {
                    let (p, d) = work.get();
                    work.set((p + 1, d));
                });
                let error = 10_f64.powf(-(obs.bq as f64) / 10.0);
                if obs.base == ref_base {
                    (1.0 - error, error / 3.0)
                } else if obs.base == alt {
                    (error / 3.0, 1.0 - error)
                } else {
                    (error / 3.0, error / 3.0)
                }
            })
            .collect();
        Self { probabilities }
    }

    fn log(&self, p: f64) -> f64 {
        self.probabilities
            .iter()
            .map(|&(zero, one)| ((1.0 - p) * zero + p * one).max(1e-300).ln())
            .sum()
    }

    fn derivative(&self, p: f64) -> f64 {
        #[cfg(test)]
        LIKELIHOOD_WORK.with(|work| {
            let (p, d) = work.get();
            work.set((p, d + 1));
        });
        self.probabilities
            .iter()
            .map(|&(zero, one)| (one - zero) / ((1.0 - p) * zero + p * one).max(1e-300))
            .sum()
    }

    fn curvature(&self, p: f64) -> f64 {
        self.probabilities
            .iter()
            .map(|&(zero, one)| {
                let denominator = ((1.0 - p) * zero + p * one).max(1e-300);
                ((one - zero) / denominator).powi(2)
            })
            .sum()
    }
}

fn log_add_exp(left: f64, right: f64) -> f64 {
    if !left.is_finite() {
        return right;
    }
    if !right.is_finite() {
        return left;
    }
    let high = left.max(right);
    high + (left.min(right) - high).exp().ln_1p()
}

fn log_uniform_likelihood(likelihood: &Likelihood) -> f64 {
    GL32_X
        .iter()
        .zip(GL32_W.iter())
        .map(|(&x, w)| likelihood.log(x) + w.ln())
        .fold(f64::NEG_INFINITY, log_add_exp)
}

/// 32-point Gauss-Legendre nodes and weights on [0, 1]. Used for every
/// quadrature in this module so that the cost of one site scales linearly
/// with its observation count, regardless of coverage depth.
#[allow(clippy::excessive_precision)]
const GL32_X: [f64; 32] = [
    0.0013680690752592151,
    0.007194244227365809,
    0.017618872206246805,
    0.03254696203113017,
    0.05183942211697395,
    0.07531619313371502,
    0.10275810201602881,
    0.13390894062985514,
    0.16847786653489238,
    0.20614212137961885,
    0.24655004553388532,
    0.28932436193468236,
    0.33406569885893617,
    0.38035631887393145,
    0.42776401920860174,
    0.4758461671561308,
    0.5241538328438692,
    0.5722359807913983,
    0.6196436811260685,
    0.6659343011410639,
    0.7106756380653176,
    0.7534499544661146,
    0.7938578786203812,
    0.8315221334651076,
    0.8660910593701449,
    0.8972418979839711,
    0.924683806866285,
    0.948160577883026,
    0.9674530379688698,
    0.9823811277937532,
    0.9928057557726342,
    0.9986319309247408,
];
#[allow(clippy::excessive_precision)]
const GL32_W: [f64; 32] = [
    0.003509305004735253,
    0.008137197365452872,
    0.01269603265463101,
    0.01713693145651088,
    0.02141794901111342,
    0.02549902963118805,
    0.02934204673926778,
    0.03291111138818084,
    0.03617289705442417,
    0.03909694789353511,
    0.04165596211347335,
    0.04382604650220189,
    0.04558693934788189,
    0.04692219954040226,
    0.04781936003963735,
    0.04827004425736383,
    0.04827004425736383,
    0.04781936003963735,
    0.04692219954040226,
    0.04558693934788189,
    0.04382604650220189,
    0.04165596211347335,
    0.03909694789353511,
    0.03617289705442417,
    0.03291111138818084,
    0.02934204673926778,
    0.02549902986561342,
    0.02141794901111342,
    0.01713693145651088,
    0.01269603265463101,
    0.008137197365452872,
    0.003509305004735253,
];

/// Fixed-node quadrature of `f` over [a, b].
fn integrate<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> f64 {
    if b <= a {
        return 0.0;
    }
    let width = b - a;
    width
        * GL32_X
            .iter()
            .zip(GL32_W.iter())
            .map(|(&x, &w)| w * f(a + width * x))
            .sum::<f64>()
}

fn likelihood_features(observations: &[&Observation], ref_base: u8, alt: u8) -> (f64, f64) {
    if observations.is_empty() {
        return (0.0, 0.0);
    }
    let likelihood = Likelihood::new(observations, ref_base, alt);
    let log_zero = likelihood.log(0.0);
    let log_one = likelihood.log(1.0);
    if observations.len() <= 63 {
        let log_uniform = log_uniform_likelihood(&likelihood);
        return (log_uniform - log_zero, log_uniform - log_one);
    }
    let derivative_zero = likelihood.derivative(0.0);
    let derivative_one = likelihood.derivative(1.0);
    let mode = if derivative_zero <= 0.0 {
        0.0
    } else if derivative_one >= 0.0 {
        1.0
    } else {
        let mut lo = 0.0;
        let mut hi = 1.0;
        for _ in 0..100 {
            let mid = (lo + hi) / 2.0;
            // A midpoint rounded to an endpoint cannot move the final rounded
            // midpoint: later updates either keep these bounds or collapse to
            // that endpoint. Retain the original 100-iteration ceiling.
            if mid == lo || mid == hi {
                break;
            }
            if likelihood.derivative(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        (lo + hi) / 2.0
    };
    let peak = likelihood.log(mode);
    let scale = if mode == 0.0 || mode == 1.0 {
        let derivative = if mode == 0.0 {
            derivative_zero.abs()
        } else {
            derivative_one.abs()
        };
        // The opposite endpoint can have an arbitrarily steep slope. Near a
        // flat boundary, curvature rather than slope controls the mass width.
        1.0 / derivative.max(likelihood.curvature(mode).sqrt()).max(1.0)
    } else {
        1.0 / likelihood.curvature(mode).max(1.0).sqrt()
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
        .map(|window| integrate(|p| (likelihood.log(p) - peak).exp(), window[0], window[1]))
        .sum::<f64>();
    (
        peak + integral.ln() - log_zero,
        peak + integral.ln() - log_one,
    )
}

// Center reductions must have a reproducible observation traversal.
fn center_observations(window: &Window, pos: i64) -> Vec<&Observation> {
    let mut observations: Vec<_> = window
        .get(&pos)
        .into_iter()
        .flat_map(|x| x.iter())
        .collect();
    // IDs are interned in BAM traversal order, independent of HashMap seeds.
    // Sort only the center reduction sequence; retain every merged observation.
    observations.sort_unstable_by_key(|&(id, _)| *id);
    observations
        .into_iter()
        .map(|(_, observation)| observation)
        .collect()
}

fn feature_vector(
    site: &(String, i64),
    ref_seq: &[u8],
    reference_start: i64,
    window: &Window,
    alu: &AluIndex,
) -> Result<[f32; FEATURE_COUNT]> {
    let (chrom, pos1) = site;
    let pos = *pos1 - 1;
    let local = usize::try_from(pos - reference_start)?;
    let ref_base = ref_seq[local].to_ascii_uppercase();
    let alt = match ref_base {
        b'A' => b'G',
        b'T' => b'C',
        _ => bail!("unsupported reference base at {chrom}:{pos1}"),
    };
    let obs = center_observations(window, pos);
    let alt_obs: Vec<_> = obs.iter().copied().filter(|o| o.base == alt).collect();
    let ref_obs: Vec<_> = obs.iter().copied().filter(|o| o.base == ref_base).collect();
    let depth = obs.len() as f64;
    let alt_count = alt_obs.len() as f64;
    let other = depth - alt_count - ref_obs.len() as f64;
    let (bf_error, bf_homo) = likelihood_features(&obs, ref_base, alt);
    let strand = if !alt_obs.is_empty() && !ref_obs.is_empty() {
        let count = |items: &[&Observation], direction: i8| {
            items.iter().filter(|o| o.orientation == direction).count() as f64
        };
        Some(
            ((count(&alt_obs, 0) + 0.5).ln()
                - (count(&alt_obs, 1) + 0.5).ln()
                - (count(&ref_obs, 0) + 0.5).ln()
                + (count(&ref_obs, 1) + 0.5).ln())
            .abs(),
        )
    } else {
        None
    };
    let mut weighted_opportunities = 0.0;
    let mut weighted_depth = 0.0;
    let mut weighted_alt = 0.0;
    let mut weighted_alt_sites = 0.0;
    let mut all_weighted_depth = 0.0;
    let mut all_weighted_other = 0.0;
    let mut template_contributions: HashMap<TemplateId, u32> = HashMap::new();
    for p in (pos - RADIUS).max(reference_start)
        ..=(pos + RADIUS).min(reference_start + ref_seq.len() as i64 - 1)
    {
        if p == pos {
            continue;
        }
        let Some(local_obs) = window.get(&p) else {
            continue;
        };
        let lr = ref_seq[(p - reference_start) as usize].to_ascii_uppercase();
        if base_index(lr).is_none() {
            continue;
        }
        let n = local_obs.len() as f64;
        let w = (-(p - pos).abs() as f64 / RADIUS as f64).exp();
        let local_ref = local_obs.values().filter(|o| o.base == lr).count() as f64;
        let local_alt = if lr == ref_base {
            local_obs.values().filter(|o| o.base == alt).count() as f64
        } else {
            0.0
        };
        let local_other = n - local_ref - local_alt;
        all_weighted_depth += w * n;
        all_weighted_other += w * local_other;
        if lr != ref_base {
            continue;
        }
        weighted_opportunities += w;
        weighted_depth += w * n;
        weighted_alt += w * local_alt;
        if local_alt > 0.0 {
            weighted_alt_sites += w;
            for (name, o) in local_obs {
                if o.base == alt {
                    *template_contributions.entry(*name).or_default() += 1;
                }
            }
        }
    }
    let alt_observations: f64 = template_contributions.values().map(|&n| f64::from(n)).sum();
    let reuse = if alt_observations > 1.0 {
        template_contributions
            .values()
            .map(|&n| {
                let n = f64::from(n);
                n * (n - 1.0)
            })
            .sum::<f64>()
            / (alt_observations * (alt_observations - 1.0))
    } else {
        0.0
    };
    let rate = |n: f64, d: f64| if d > 0.0 { n / d } else { f64::NAN };
    let previous = if pos > reference_start {
        ref_seq[local - 1]
    } else {
        b'N'
    };
    let following = if local + 1 < ref_seq.len() {
        ref_seq[local + 1]
    } else {
        b'N'
    };
    let (previous, following) = if ref_base == b'T' || ref_base == b'G' {
        (complement(following), complement(previous))
    } else {
        (previous, following)
    };
    let mut values = [f64::NAN; FEATURE_COUNT];
    values[..20].copy_from_slice(&[
        depth.max(0.0).ln_1p(),
        alt_count.ln_1p(),
        rate(alt_count, depth),
        rate(other, depth),
        bf_error,
        bf_homo,
        mean(alt_obs.iter().map(|o| f64::from(o.bq))),
        mean(ref_obs.iter().map(|o| f64::from(o.bq))),
        mean(alt_obs.iter().map(|o| o.mapq)),
        mean(alt_obs.iter().map(|o| o.end_fraction)),
        strand.unwrap_or(f64::NAN),
        mean(alt_obs.iter().map(|o| o.softclip_fraction)),
        mean(alt_obs.iter().map(|o| o.gap_nearness)),
        rate(weighted_alt, weighted_depth),
        rate(weighted_alt_sites, weighted_opportunities),
        rate(all_weighted_other, all_weighted_depth),
        (template_contributions.len() as f64).ln_1p(),
        reuse,
        weighted_opportunities.ln_1p(),
        weighted_depth.ln_1p(),
    ]);
    for (i, base) in [previous, following].into_iter().enumerate() {
        for (j, target) in (*b"ACG").into_iter().enumerate() {
            values[20 + i * 3 + j] = if base == target {
                1.0
            } else if base_index(base).is_some() {
                0.0
            } else {
                f64::NAN
            };
        }
    }
    values[26] = if alu.contains(chrom, pos) { 1.0 } else { 0.0 };
    Ok(values.map(|value| value as f32))
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

pub fn score_sites(
    bam_path: &Path,
    fasta_path: &Path,
    model_dir: &Path,
    alu: &AluIndex,
    sites: &[(String, i64)],
    threads: usize,
) -> Result<Vec<f64>> {
    let model = Model::load(model_dir)?;
    let mut order: Vec<usize> = (0..sites.len()).collect();
    order.sort_by(|&a, &b| sites[a].cmp(&sites[b]));
    let mut groups = Vec::new();
    let mut start = 0;
    while start < order.len() {
        let first = &sites[order[start]];
        let mut end = start + 1;
        let mut right = first.1 - 1 + RADIUS + 1;
        while end < order.len() {
            let next = &sites[order[end]];
            let left = next.1 - 1 - RADIUS;
            if next.0 != first.0 || left > right {
                break;
            }
            right = right.max(next.1 - 1 + RADIUS + 1);
            end += 1;
        }
        groups.push((start, end, first.0.clone(), first.1 - 1 - RADIUS));
        start = end;
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()?;
    let scores = pool.install(|| {
        groups
            .par_iter()
            .map(|&(lo, hi, ref chrom, left)| -> Result<Vec<(usize, f64)>> {
                let right = sites[order[hi - 1]].1 - 1 + RADIUS + 1;
                let mut bam = IndexedReader::from_path(bam_path)?;
                let fasta = rust_htslib::faidx::Reader::from_path(fasta_path)?;
                let ref_start = left.max(0);
                let reference = fasta.fetch_seq(chrom, ref_start as usize, (right - 1) as usize)?;
                let window = read_window(&mut bam, chrom, ref_start, right)?;
                let mut result = Vec::with_capacity(hi - lo);
                for &i in &order[lo..hi] {
                    let features = feature_vector(&sites[i], &reference, ref_start, &window, alu)?;
                    result.push((i, model.score(&features)));
                }
                Ok(result)
            })
            .collect::<Vec<_>>()
    });
    let mut output = vec![0.0; sites.len()];
    for result in scores {
        for (i, score) in result? {
            output[i] = score;
        }
    }
    Ok(output)
}

/// Feature vector or per-site failure, as produced by [`feature_matrix`].
pub type FeatureRow = std::result::Result<[f32; FEATURE_COUNT], String>;

/// Export raw 27-feature vectors for arbitrary sites, mirroring the grouping
/// and window logic of `score_sites`. Per-site failures (e.g. unsupported
/// reference base) are captured as `Err(message)` instead of aborting.
///
/// Groups are additionally capped in size/span for rayon load balance; this
/// does not change any feature value (a site's window depends only on its own
/// +-RADIUS neighbourhood). BAM and FASTA handles are reused within each rayon
/// task; every group still fetches its original independent window.
pub fn feature_matrix(
    bam_path: &Path,
    fasta_path: &Path,
    alu: &AluIndex,
    sites: &[(String, i64)],
    threads: usize,
) -> Result<Vec<FeatureRow>> {
    feature_matrix_with_readers(alu, sites, threads, || {
        Ok((
            IndexedReader::from_path(bam_path)?,
            rust_htslib::faidx::Reader::from_path(fasta_path)?,
        ))
    })
}

fn feature_matrix_with_readers(
    alu: &AluIndex,
    sites: &[(String, i64)],
    threads: usize,
    open: impl Fn() -> Result<(IndexedReader, rust_htslib::faidx::Reader)> + Sync,
) -> Result<Vec<FeatureRow>> {
    let mut order: Vec<usize> = (0..sites.len()).collect();
    order.sort_by(|&a, &b| sites[a].cmp(&sites[b]));
    const MAX_GROUP_SITES: usize = 512;
    const MAX_GROUP_SPAN: i64 = 401;
    let mut groups = Vec::new();
    let mut start = 0;
    while start < order.len() {
        let first = &sites[order[start]];
        let mut end = start + 1;
        let mut right = first.1 - 1 + RADIUS + 1;
        while end < order.len() {
            let next = &sites[order[end]];
            let left = next.1 - 1 - RADIUS;
            if next.0 != first.0 || left > right {
                break;
            }
            if end - start >= MAX_GROUP_SITES
                || next.1 - 1 + RADIUS + 1 - (first.1 - 1 - RADIUS) > MAX_GROUP_SPAN
            {
                break;
            }
            right = right.max(next.1 - 1 + RADIUS + 1);
            end += 1;
        }
        groups.push((start, end, first.0.clone(), first.1 - 1 - RADIUS));
        start = end;
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()?;
    let results = pool.install(|| {
        groups
            .par_iter()
            .map_init(
                || None,
                |readers, &(lo, hi, ref chrom, left)| -> Result<Vec<(usize, FeatureRow)>> {
                    let right = sites[order[hi - 1]].1 - 1 + RADIUS + 1;
                    if readers.is_none() {
                        *readers = Some(open()?);
                    }
                    let (bam, fasta) = readers.as_mut().expect("readers initialized above");
                    let ref_start = left.max(0);
                    let reference =
                        fasta.fetch_seq(chrom, ref_start as usize, (right - 1) as usize)?;
                    let window = read_window(bam, chrom, ref_start, right)?;
                    let mut result = Vec::with_capacity(hi - lo);
                    for &i in &order[lo..hi] {
                        result.push((
                            i,
                            feature_vector(&sites[i], &reference, ref_start, &window, alu)
                                .map_err(|e| e.to_string()),
                        ));
                    }
                    Ok(result)
                },
            )
            .collect::<Vec<_>>()
    });
    let mut output = vec![Err("not computed".to_string()); sites.len()];
    for result in results {
        for (i, value) in result? {
            output[i] = value;
        }
    }
    Ok(output)
}

fn observe_block(
    bam_path: &Path,
    fasta_path: &Path,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<Vec<ObservationSite>> {
    let fasta = rust_htslib::faidx::Reader::from_path(fasta_path)?;
    let reference = fasta.fetch_seq(chrom, start as usize, (end - 1) as usize)?;
    let mut bam = IndexedReader::from_path(bam_path)?;
    let tid = bam
        .header()
        .tid(chrom.as_bytes())
        .context("contig absent from BAM")?;
    bam.fetch((tid, start, end))?;
    let mut depth = vec![0u32; reference.len()];
    let mut counts = vec![[0u32; 4]; reference.len()];
    for record in bam.records() {
        let record = record?;
        for (pos, observation) in read_observations(&record, start, end) {
            let local = usize::try_from(pos - start)?;
            let Some(base) = base_index(observation.base) else {
                continue;
            };
            if base_index(reference[local]).is_none() {
                continue;
            }
            depth[local] += 1;
            counts[local][base] += 1;
        }
    }
    let mut sites = Vec::new();
    for (local, &ref_base) in reference.iter().enumerate() {
        let (alt, alt_i) = match ref_base.to_ascii_uppercase() {
            b'A' => (b'G', 2),
            b'T' => (b'C', 1),
            _ => continue,
        };
        if counts[local][alt_i] == 0 {
            continue;
        }
        sites.push(ObservationSite {
            chrom: chrom.to_owned(),
            pos: start + local as i64 + 1,
            ref_base,
            counts: counts[local],
            depth: depth[local],
        });
        let _ = alt;
    }
    Ok(sites)
}

pub fn discover_sites(
    bam_path: &Path,
    fasta_path: &Path,
    out: &Path,
    threads: usize,
) -> Result<usize> {
    let fai = std::fs::read_to_string(format!("{}.fai", fasta_path.display()))?;
    let lengths: HashMap<String, i64> = fai
        .lines()
        .filter_map(|line| {
            let mut f = line.split('\t');
            Some((f.next()?.to_owned(), f.next()?.parse().ok()?))
        })
        .collect();
    let reader = IndexedReader::from_path(bam_path)?;
    let mut contigs = Vec::new();
    for tid in 0..reader.header().target_count() {
        let chrom = String::from_utf8(reader.header().tid2name(tid).to_vec())?;
        if let Some(&len) = lengths.get(&chrom) {
            contigs.push((chrom, len));
        }
    }
    let mut blocks = Vec::new();
    for (chrom, len) in contigs {
        let mut start = 0;
        while start < len {
            blocks.push((chrom.clone(), start, (start + 8_000_000).min(len)));
            start += 8_000_000;
        }
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()?;
    let mut all = pool
        .install(|| {
            blocks
                .par_iter()
                .map(|(chrom, start, end)| observe_block(bam_path, fasta_path, chrom, *start, *end))
                .collect::<Result<Vec<_>>>()
        })?
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    all.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.pos.cmp(&b.pos)));
    let mut text = String::new();
    for site in &all {
        let alt = if site.ref_base == b'A' { b'G' } else { b'C' };
        let freq = site.counts[if alt == b'G' { 2 } else { 1 }] as f64 / site.depth.max(1) as f64;
        use std::fmt::Write as _;
        writeln!(
            text,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{freq:.8}\t",
            site.chrom,
            site.pos - 1,
            site.pos,
            site.ref_base as char,
            site.depth,
            site.counts[0],
            site.counts[1],
            site.counts[2],
            site.counts[3]
        )?;
    }
    std::fs::write(out, text)?;
    Ok(all.len())
}

#[cfg(test)]
mod template_tie_regression {
    use super::*;

    fn obs(base: u8, bq: u8, mapq: f64, orientation: i8) -> Observation {
        Observation {
            base,
            bq,
            mapq,
            end_fraction: mapq / 100.0,
            orientation,
            softclip_fraction: (100.0 - mapq) / 100.0,
            gap_nearness: (100.0 - mapq) / 100.0,
        }
    }

    fn signature(value: &Observation) -> Vec<u64> {
        vec![
            value.base as u64,
            value.bq as u64,
            if value.mapq.is_nan() {
                u64::MAX
            } else {
                value.mapq.to_bits()
            },
            value.end_fraction.to_bits(),
            (value.orientation as i64) as u64,
            value.softclip_fraction.to_bits(),
            value.gap_nearness.to_bits(),
        ]
    }

    #[test]
    fn equal_bq_same_base_metadata_is_commutative_and_associative() {
        let records = [
            obs(b'G', 30, 20.0, 0),
            obs(b'G', 30, 60.0, 1),
            obs(b'G', 30, f64::NAN, 1),
        ];
        let mut expected = None;
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    continue;
                }
                let k = 3 - i - j;
                let mut left = records[i].clone();
                merge_observation(&mut left, records[j].clone());
                merge_observation(&mut left, records[k].clone());
                let mut right = records[j].clone();
                merge_observation(&mut right, records[k].clone());
                let mut grouped = records[i].clone();
                merge_observation(&mut grouped, right);
                assert_eq!(signature(&left), signature(&grouped));
                if let Some(ref e) = expected {
                    assert_eq!(&signature(&left), e, "input-order metadata change");
                } else {
                    expected = Some(signature(&left));
                }
                assert_eq!(left.mapq, 60.0);
                assert_eq!(left.orientation, -1);
            }
        }
    }

    #[test]
    fn conflict_is_absorbing_until_a_strictly_higher_bq_resolves_it() {
        let mut ambiguous = obs(b'A', 30, 20.0, 0);
        merge_observation(&mut ambiguous, obs(b'G', 30, 60.0, 1));
        assert_eq!(ambiguous.base, b'N');
        merge_observation(&mut ambiguous, obs(b'A', 30, 40.0, 0));
        assert_eq!(ambiguous.base, b'N');
        let winner = obs(b'G', 35, 40.0, 0);
        merge_observation(&mut ambiguous, winner.clone());
        assert_eq!(signature(&ambiguous), signature(&winner));
        // Lower-BQ ties/metadata cannot contaminate the maximum-BQ winner.
        merge_observation(&mut ambiguous, obs(b'A', 30, 60.0, 1));
        assert_eq!(signature(&ambiguous), signature(&winner));
    }
}

#[cfg(test)]
mod likelihood_numerical_regression {
    use super::*;

    fn observations(groups: &[(u8, u8, usize)]) -> Vec<Observation> {
        groups
            .iter()
            .flat_map(|&(base, bq, count)| {
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
                    count,
                )
            })
            .collect()
    }

    fn check(groups: &[(u8, u8, usize)], expected: (f64, f64)) {
        let owned = observations(groups);
        let refs: Vec<_> = owned.iter().collect();
        let actual = likelihood_features(&refs, b'A', b'G');
        for (a, e) in [(actual.0, expected.0), (actual.1, expected.1)] {
            assert!(
                a.is_finite() && (a - e).abs() < 2e-6,
                "groups={groups:?}: actual={actual:?}, expected={expected:?}"
            );
        }
    }

    #[test]
    fn low_depth_third_bases_do_not_underflow() {
        // Independent integral: all third-base factors are constant in p.
        check(&[(b'C', 93, 63)], (0.0, 0.0));
        // Three ALT factors integrate to ((1-e)^4-(e/3)^4)/(4*(1-4e/3)).
        check(
            &[(b'C', 93, 60), (b'G', 93, 3)],
            (66.15166659808187, -1.386294360952661),
        );
    }

    #[test]
    fn boundary_mass_is_not_lost_at_high_depth() {
        // Independently checked by scipy quadrature AND incomplete beta form.
        check(
            &[(b'A', 30, 9997), (b'G', 30, 3)],
            (-8.254418713486338, 79997.38405204924),
        );
        check(
            &[(b'A', 30, 3), (b'G', 30, 9997)],
            (79997.38405204924, -8.254418713486338),
        );
        // Almost-flat boundary: derivative alone is not a local width.
        check(
            &[(b'A', 30, 2997), (b'G', 30, 1)],
            (-7.3125532755095435, 23976.767180405215),
        );
    }

    #[test]
    fn mixed_quality_and_third_base_oracle() {
        check(
            &[
                (b'A', 20, 6000),
                (b'A', 41, 3980),
                (b'G', 20, 10),
                (b'G', 41, 5),
                (b'C', 30, 5),
            ],
            (4.272870200591683, 76002.77702414071),
        );
    }

    #[test]
    fn pure_ref_and_alt_have_closed_form_boundary_integrals() {
        for n in [63, 64, 1000, 10000] {
            for q in [13, 20, 30, 60, 93] {
                let e = 10_f64.powf(-(q as f64) / 10.0);
                let high = 1.0 - e;
                let low = e / 3.0;
                // Integral (high - (high-low)*p)^n dp on [0,1].
                let log_z = (n as f64 + 1.0) * high.ln()
                    + (-(low / high).powi(n as i32 + 1)).ln_1p()
                    - (high - low).ln()
                    - (n as f64 + 1.0).ln();
                let expected = (log_z - n as f64 * high.ln(), log_z - n as f64 * low.ln());
                check(&[(b'A', q, n)], expected);
                check(&[(b'G', q, n)], (expected.1, expected.0));
            }
        }
    }
}

#[cfg(test)]
mod memory_regression {
    use super::*;
    #[test]
    fn position_keys_do_not_own_repeated_template_names() {
        let observation = Observation {
            base: b'G',
            bq: 30,
            mapq: 60.0,
            end_fraction: 0.5,
            orientation: 0,
            softclip_fraction: 0.0,
            gap_nearness: 0.0,
        };
        let mut window = Window::new();
        window.entry(10).or_default().insert(0, observation.clone());
        let key = window.get(&10).unwrap().keys().next().unwrap();
        assert!(
            std::mem::size_of_val(key) <= 8,
            "each position stores an owned template name; compact identity required"
        );
    }
}

#[cfg(test)]
mod window_semantic_regression {
    use super::*;
    use rust_htslib::bam::{
        self,
        header::HeaderRecord,
        record::{Cigar, CigarString},
        Header,
    };
    #[test]
    fn interning_retains_rg_mates_conflicts_and_late_resolution() {
        let path =
            std::env::temp_dir().join(format!("esp-window-intern-{}.bam", std::process::id()));
        let mut header = Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 1000),
        );
        {
            let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            for (name, rg, base, quality, flag) in [
                ("same/1", "A", b'A', 30, 0),
                ("same/2", "A", b'G', 30, 0),
                ("same/1", "A", b'G', 40, 0),
                ("same/1", "B", b'A', 30, 0),
                ("conflict/1", "A", b'A', 30, 0),
                ("conflict/2", "A", b'G', 30, 0),
                ("secondary", "A", b'G', 40, 256),
            ] {
                let mut record = bam::Record::new();
                record.set(
                    name.as_bytes(),
                    Some(&CigarString(vec![Cigar::Match(2)])),
                    &[base; 2],
                    &[quality; 2],
                );
                record.set_tid(0);
                record.set_pos(100);
                record.set_mapq(60);
                record.set_flags(flag);
                record
                    .push_aux(b"RG", bam::record::Aux::String(rg))
                    .unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
        let mut bam = IndexedReader::from_path(&path).unwrap();
        let window = read_window(&mut bam, "chr1", 99, 104).unwrap();
        for position in [100, 101] {
            let observations = &window[&position];
            assert_eq!(observations.len(), 2);
            assert_eq!(
                observations
                    .values()
                    .filter(|o| o.base == b'G' && o.bq == 40)
                    .count(),
                1
            );
            assert_eq!(observations.values().filter(|o| o.base == b'A').count(), 1);
        }
        std::fs::remove_file(&path).unwrap();
        std::fs::remove_file(format!("{}.bai", path.display())).unwrap();
    }
}
