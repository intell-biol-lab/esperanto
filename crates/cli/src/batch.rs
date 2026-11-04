//! Bounded external sample workers; scientific entry points stay unchanged.
use crate::evidence;
use anyhow::{ensure, Context, Result};
use clap::{Args, ValueEnum};
use esperanto_flow::samples::{read_samples, Sample};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum Mode {
    Run,
    Map,
}

#[derive(Args)]
pub struct BatchArgs {
    /// TSV sample sheet; relative inputs are resolved against this file.
    #[arg(long)]
    samples: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, value_enum, default_value = "run")]
    mode: Mode,
    /// Total worker-thread budget (0 = available parallelism).
    #[arg(long, default_value_t = 0)]
    total_threads: usize,
    /// Complete-content input hash readers, capped by total threads and file count.
    /// Validation finishes before sample workers; each reader uses a 1 MiB buffer.
    #[arg(long, default_value_t = 4)]
    validation_threads: usize,
    #[arg(long, default_value_t = 4)]
    threads_per_sample: usize,
    /// Total estimated memory budget; reservations are not OS-enforced limits.
    #[arg(long)]
    memory_gib: usize,
    #[arg(long)]
    memory_per_sample_gib: usize,
    /// Optional additional concurrency cap.
    #[arg(long)]
    jobs: Option<usize>,
    #[arg(long)]
    index: Option<PathBuf>,
    #[arg(long)]
    fasta: Option<PathBuf>,
    #[arg(long)]
    gtf: Option<PathBuf>,
    #[arg(long)]
    bundle: Option<PathBuf>,
    #[arg(long)]
    l1_bundle: Option<PathBuf>,
    /// Explicit worker executable (default: this exact executable).
    #[arg(long)]
    executable: Option<PathBuf>,
    /// Full-pipeline stop stage; not supported with --mode map.
    #[arg(long)]
    stop_after: Option<String>,
    /// Verify pinned inputs/configuration and completed outputs before continuing.
    #[arg(long)]
    resume: bool,
    /// Publish the execution plan without starting workers; continue with --resume.
    #[arg(long)]
    dry_run: bool,
}

#[derive(Serialize, Deserialize)]
struct Job {
    sample: String,
    command: Vec<String>,
    status: String,
    attempts: usize,
    started_unix_seconds: Option<f64>,
    wall_seconds: Option<f64>,
    exit_code: Option<i32>,
    outputs: BTreeMap<PathBuf, String>,
    error: Option<String>,
}

#[derive(Serialize, Deserialize)]
struct State {
    schema_version: usize,
    identity: Value,
    parallel: usize,
    jobs: Vec<Job>,
}

struct Running {
    child: Child,
    job: usize,
    start: Instant,
}
impl Drop for Running {
    fn drop(&mut self) {
        // On coordinator errors, do not leave writers running behind the lock.
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

fn parallel(
    total: usize,
    each: usize,
    memory: usize,
    per: usize,
    jobs: Option<usize>,
) -> Result<usize> {
    ensure!(
        total > 0 && each > 0 && memory > 0 && per > 0,
        "resource budgets must be positive"
    );
    ensure!(jobs != Some(0), "--jobs must be positive");
    let count = (total / each)
        .min(memory / per)
        .min(jobs.unwrap_or(usize::MAX));
    ensure!(count > 0, "one sample exceeds total CPU or memory budget");
    Ok(count)
}

fn absolute(path: &mut Option<PathBuf>) -> Result<()> {
    if let Some(p) = path {
        *p = p.canonicalize()?;
    }
    Ok(())
}

fn add(args: &mut Vec<String>, flag: &str, path: &Path) -> Result<()> {
    args.extend([
        flag.to_owned(),
        path.to_str().context("non-UTF8 worker path")?.to_owned(),
    ]);
    Ok(())
}

fn command(a: &BatchArgs, sample: &Sample, executable: &Path) -> Result<Vec<String>> {
    let mut args = vec![
        executable
            .to_str()
            .context("non-UTF8 executable")?
            .to_owned(),
        if a.mode == Mode::Map { "map" } else { "run" }.to_owned(),
    ];
    for (flag, paths) in [("--r1", &sample.r1), ("--r2", &sample.r2)] {
        for path in paths {
            add(&mut args, flag, path)?;
        }
    }
    for (flag, path) in [
        ("--index", a.index.as_deref()),
        ("--gtf", a.gtf.as_deref()),
        ("--l1-bundle", a.l1_bundle.as_deref()),
    ] {
        if let Some(path) = path {
            add(&mut args, flag, path)?;
        }
    }
    let root = a.out.join("samples");
    if a.mode == Mode::Run {
        add(&mut args, "--fasta", a.fasta.as_deref().unwrap())?;
        add(&mut args, "--bundle", a.bundle.as_deref().unwrap())?;
        if let Some(path) = &sample.bam {
            add(&mut args, "--bam", path)?;
        }
        if let Some(path) = &sample.sites {
            add(&mut args, "--sites", path)?;
        }
        args.extend(["--sample".to_owned(), sample.sample.clone()]);
        if let Some(stop) = &a.stop_after {
            args.extend(["--stop-after".to_owned(), stop.clone()]);
        }
        add(&mut args, "--out", &root)?;
    } else {
        add(&mut args, "--out", &root.join(&sample.sample))?;
    }
    args.extend(["--threads".to_owned(), a.threads_per_sample.to_string()]);
    Ok(args)
}

pub fn run(mut a: BatchArgs) -> Result<()> {
    let invocation_start = Instant::now();
    let total = crate::resolve::threads(a.total_threads);
    ensure!(
        a.validation_threads > 0,
        "--validation-threads must be positive"
    );
    let concurrency = parallel(
        total,
        a.threads_per_sample,
        a.memory_gib,
        a.memory_per_sample_gib,
        a.jobs,
    )?;
    let samples = read_samples(&a.samples)?;
    ensure!(
        a.mode == Mode::Run || a.stop_after.is_none(),
        "--stop-after requires --mode run"
    );
    if let Some(stop) = &a.stop_after {
        ensure!(
            ["qc", "map", "sort", "scan", "score", "vcf"].contains(&stop.as_str()),
            "invalid stop stage"
        );
    }
    if a.mode == Mode::Map {
        ensure!(
            samples.iter().all(|s| s.r1.len() == 1 && s.bam.is_none()),
            "map mode requires one R1 per sample; use run for lanes/BAM"
        );
    }
    // Resolve discovery once and pin explicit paths in every command.
    if let Some(refs) = crate::resolve::refs() {
        refs.fill_fasta(&mut a.fasta);
        refs.fill_index(&mut a.index);
        refs.fill_gtf(&mut a.gtf);
    }
    ensure!(a.fasta.is_some(), "--fasta is required");
    if samples.iter().any(|s| !s.r1.is_empty()) {
        let index = a.index.as_ref().context("FASTQ samples require --index")?;
        a.l1_bundle = Some(crate::resolve::l1_bundle(&a.l1_bundle, index)?);
    }
    if a.mode == Mode::Run {
        a.bundle = Some(crate::resolve::bundle(&a.bundle)?);
    }
    for path in [
        &mut a.fasta,
        &mut a.index,
        &mut a.gtf,
        &mut a.bundle,
        &mut a.l1_bundle,
    ] {
        absolute(path)?;
    }
    let executable = a
        .executable
        .clone()
        .unwrap_or(std::env::current_exe()?)
        .canonicalize()?;
    std::fs::create_dir_all(&a.out)?;
    a.out = a.out.canonicalize()?;
    let _lock = esperanto_flow::resume::RunLock::acquire(&a.out)?;
    let mut paths = vec![
        a.samples.canonicalize()?,
        executable.clone(),
        std::env::current_exe()?.canonicalize()?,
    ];
    for sample in &samples {
        paths.extend(
            sample
                .r1
                .iter()
                .chain(&sample.r2)
                .chain(sample.bam.iter())
                .chain(sample.sites.iter())
                .chain(sample.candidates.iter())
                .cloned(),
        );
        if let Some(bam) = &sample.bam {
            for idx in [
                PathBuf::from(format!("{}.bai", bam.display())),
                bam.with_extension("bai"),
                PathBuf::from(format!("{}.csi", bam.display())),
                bam.with_extension("csi"),
            ] {
                if idx.is_file() {
                    paths.push(idx);
                }
            }
        }
    }
    for path in [&a.fasta, &a.index, &a.gtf, &a.bundle, &a.l1_bundle]
        .into_iter()
        .flatten()
    {
        paths.push(path.clone());
    }
    paths.push(PathBuf::from(format!(
        "{}.fai",
        a.fasta.as_ref().unwrap().display()
    )));
    if let Some(l1) = &a.l1_bundle {
        paths.push(l1.with_extension("tidx"));
    }
    if let Some(index) = &a.index {
        for kind in ["ag.paidx", "tc.paidx"] {
            let p = index.with_extension(kind);
            if p.exists() {
                paths.push(p);
            }
        }
    }
    let controls: BTreeMap<String, String> = std::env::vars()
        .filter(|(key, _)| {
            key.starts_with("ESP_")
                || key.starts_with("ESPERANTO_")
                || [
                    "HOME",
                    "PATH",
                    "XDG_DATA_HOME",
                    "XDG_CACHE_HOME",
                    "OMP_NUM_THREADS",
                    "RAYON_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                ]
                .contains(&key.as_str())
        })
        .collect();
    let commands = samples
        .iter()
        .map(|s| command(&a, s, &executable))
        .collect::<Result<Vec<_>>>()?;
    let hash_start = Instant::now();
    let hashes = evidence::identities_parallel(&paths, a.validation_threads.min(total))?;
    let input_hash_seconds = hash_start.elapsed().as_secs_f64();
    let input_bytes: u64 = hashes
        .keys()
        .map(|p| std::fs::metadata(p).map(|m| m.len()))
        .collect::<std::io::Result<Vec<_>>>()?
        .into_iter()
        .sum();
    let validation_threads = a.validation_threads.min(total).min(hashes.len());
    let input_files = hashes.len();
    let identity = json!({
        "samples": samples, "sha256": hashes,
        "commands": commands, "environment": controls,
        "worker_environment_overrides":{
            "OMP_NUM_THREADS":a.threads_per_sample.to_string(),
            "RAYON_NUM_THREADS":a.threads_per_sample.to_string()
        },
        "mode": format!("{:?}",a.mode), "total_threads":total,
        "threads_per_sample":a.threads_per_sample, "memory_gib":a.memory_gib,
        "memory_per_sample_gib":a.memory_per_sample_gib, "parallel":concurrency,
        "memory_policy":"estimated reservation, not enforced peak"
    });
    let state_path = a.out.join("batch.json");
    let mut state = if a.resume {
        let old: State = serde_json::from_slice(&std::fs::read(&state_path)?)?;
        ensure!(
            old.schema_version == 1 && old.identity == identity,
            "batch identity changed; use a new output directory"
        );
        ensure!(
            old.parallel == concurrency && old.jobs.len() == samples.len(),
            "batch state structure changed"
        );
        for ((job, sample), command) in old.jobs.iter().zip(&samples).zip(&commands) {
            ensure!(
                job.sample == sample.sample && job.command == *command,
                "batch job definition changed"
            );
        }
        for job in &old.jobs {
            if job.status == "succeeded" {
                let current = evidence::identities(&[a.out.join("samples").join(&job.sample)])?;
                ensure!(
                    current == job.outputs,
                    "completed output changed for {}; use a new directory",
                    job.sample
                );
            }
        }
        old
    } else {
        ensure!(!state_path.exists(), "batch already exists; use --resume");
        ensure!(
            !a.out.join("samples").exists(),
            "sample output directory already exists"
        );
        State {
            schema_version: 1,
            identity,
            parallel: concurrency,
            jobs: samples
                .iter()
                .zip(commands)
                .map(|(s, c)| Job {
                    sample: s.sample.clone(),
                    command: c,
                    status: "pending".to_owned(),
                    attempts: 0,
                    started_unix_seconds: None,
                    wall_seconds: None,
                    exit_code: None,
                    outputs: BTreeMap::new(),
                    error: None,
                })
                .collect(),
        }
    };
    evidence::write_json(&state_path, &state)?;
    let mut attempt = 1;
    while a.out.join(format!("validation-{attempt}.json")).exists() {
        attempt += 1;
    }
    let receipt_path = a.out.join(format!("validation-{attempt}.json"));
    let mut receipt = json!({
        "schema_version":1, "status":"prepared",
        "validation_threads":validation_threads,
        "requested_validation_threads":a.validation_threads,
        "total_threads":total, "input_files":input_files, "input_bytes":input_bytes,
        "input_hash_seconds":input_hash_seconds,
        "preparation_seconds":invocation_start.elapsed().as_secs_f64(),
        "output_hash_seconds":0.0,
        "policy":"complete-content SHA256; readers joined before sample workers; output hashing serial"
    });
    evidence::write_json(&receipt_path, &receipt)?;
    if a.dry_run {
        receipt["status"] = json!("dry_run");
        receipt["wall_seconds"] = json!(invocation_start.elapsed().as_secs_f64());
        evidence::write_json(&receipt_path, &receipt)?;
        return Ok(());
    }
    let worker_phase = Instant::now();
    let mut output_hash_seconds = 0.0;
    std::fs::create_dir_all(a.out.join("logs"))?;
    std::fs::create_dir_all(a.out.join("samples"))?;
    let pending: Vec<_> = state
        .jobs
        .iter()
        .enumerate()
        .filter(|(_, j)| j.status != "succeeded")
        .map(|(i, _)| i)
        .collect();
    let mut next = 0;
    let mut running: Vec<Running> = Vec::new();
    while next < pending.len() || !running.is_empty() {
        while next < pending.len() && running.len() < concurrency {
            let i = pending[next];
            next += 1;
            let job = &mut state.jobs[i];
            let dest = a.out.join("samples").join(&job.sample);
            let mut args = job.command[1..].to_vec();
            if dest.exists() {
                if a.mode == Mode::Run && dest.join("run.json").exists() && a.stop_after.is_none() {
                    args = vec!["resume".to_owned(), dest.to_str().unwrap().to_owned()];
                } else {
                    let archive = a.out.join("failed-attempts");
                    std::fs::create_dir_all(&archive)?;
                    std::fs::rename(
                        &dest,
                        archive.join(format!("{}.attempt-{}", job.sample, job.attempts)),
                    )?;
                }
            }
            job.attempts += 1;
            let prefix = a
                .out
                .join("logs")
                .join(format!("{}.attempt-{}", job.sample, job.attempts));
            let stdout = std::fs::File::create(
                prefix.with_extension(format!("attempt-{}.stdout.log", job.attempts)),
            )?;
            let stderr = std::fs::File::create(
                prefix.with_extension(format!("attempt-{}.stderr.log", job.attempts)),
            )?;
            evidence::write_json(
                &a.out.join("logs").join(format!(
                    "{}.attempt-{}.command.json",
                    job.sample, job.attempts
                )),
                &json!({"executable":executable,"args":args}),
            )?;
            job.started_unix_seconds =
                Some(SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs_f64());
            job.status = "running".to_owned();
            job.error = None;
            job.exit_code = None;
            job.wall_seconds = None;
            job.outputs.clear();
            let start = Instant::now();
            match Command::new(&executable)
                .args(&args)
                .env("OMP_NUM_THREADS", a.threads_per_sample.to_string())
                .env("RAYON_NUM_THREADS", a.threads_per_sample.to_string())
                .stdout(Stdio::from(stdout))
                .stderr(Stdio::from(stderr))
                .spawn()
            {
                Ok(child) => running.push(Running {
                    child,
                    job: i,
                    start,
                }),
                Err(error) => {
                    job.status = "failed".to_owned();
                    job.error = Some(error.to_string());
                    evidence::write_json(
                        &a.out.join("logs").join(format!(
                            "{}.attempt-{}.execution.json",
                            job.sample, job.attempts
                        )),
                        &job,
                    )?;
                }
            }
            evidence::write_json(&state_path, &state)?;
        }
        let mut n = 0;
        while n < running.len() {
            if let Some(status) = running[n].child.try_wait()? {
                let active = running.swap_remove(n);
                let job = &mut state.jobs[active.job];
                job.wall_seconds = Some(active.start.elapsed().as_secs_f64());
                job.exit_code = status.code();
                job.status = if status.success() {
                    "succeeded"
                } else {
                    "failed"
                }
                .to_owned();
                if status.success() {
                    let output_start = Instant::now();
                    match evidence::identities(&[a.out.join("samples").join(&job.sample)]) {
                        Ok(outputs) if !outputs.is_empty() => job.outputs = outputs,
                        outcome => {
                            job.status = "failed".to_owned();
                            job.error = Some(match outcome {
                                Err(error) => error.to_string(),
                                _ => "worker exited successfully without output".to_owned(),
                            });
                        }
                    }
                    output_hash_seconds += output_start.elapsed().as_secs_f64();
                }
                eprintln!("[batch] {}: {}", job.sample, job.status);
                evidence::write_json(
                    &a.out.join("logs").join(format!(
                        "{}.attempt-{}.execution.json",
                        job.sample, job.attempts
                    )),
                    &job,
                )?;
                evidence::write_json(&state_path, &state)?;
            } else {
                n += 1;
            }
        }
        if !running.is_empty() {
            std::thread::sleep(Duration::from_millis(50));
        }
    }
    let succeeded = state.jobs.iter().all(|j| j.status == "succeeded");
    receipt["status"] = json!(if succeeded { "succeeded" } else { "failed" });
    receipt["worker_phase_seconds"] = json!(worker_phase.elapsed().as_secs_f64());
    receipt["output_hash_seconds"] = json!(output_hash_seconds);
    receipt["wall_seconds"] = json!(invocation_start.elapsed().as_secs_f64());
    evidence::write_json(&receipt_path, &receipt)?;
    ensure!(
        succeeded,
        "one or more samples failed; inspect batch.json and logs, then --resume"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reservations_bound_cpu_and_memory_independently() {
        assert_eq!(parallel(32, 4, 96, 12, None).unwrap(), 8);
        assert_eq!(parallel(32, 4, 48, 12, None).unwrap(), 4);
        assert_eq!(parallel(32, 4, 96, 12, Some(2)).unwrap(), 2);
        for (t, e, m, p, j) in [
            (0, 4, 96, 12, None),
            (32, 0, 96, 12, None),
            (2, 4, 96, 12, None),
            (32, 4, 1, 12, None),
            (32, 4, 96, 0, None),
            (32, 4, 96, 12, Some(0)),
        ] {
            assert!(parallel(t, e, m, p, j).is_err());
        }
    }
}
