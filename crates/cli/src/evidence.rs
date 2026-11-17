//! Read-only content identities and atomic evidence publication.
use anyhow::{ensure, Result};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

pub fn hash(path: &Path) -> Result<String> {
    let mut input = std::fs::File::open(path)?;
    let mut digest = Sha256::new();
    let mut buffer = vec![0u8; 1024 * 1024];
    loop {
        let n = input.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        digest.update(&buffer[..n]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

pub fn files(path: &Path) -> Result<Vec<PathBuf>> {
    if path.is_file() {
        return Ok(vec![path.to_owned()]);
    }
    ensure!(path.is_dir(), "missing resource {}", path.display());
    let mut result = Vec::new();
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        ensure!(
            !entry.file_type()?.is_symlink(),
            "directory resource contains symlink: {}",
            entry.path().display()
        );
        result.extend(files(&entry.path())?);
    }
    result.sort();
    Ok(result)
}

pub fn identities(paths: &[PathBuf]) -> Result<BTreeMap<PathBuf, String>> {
    let mut result = BTreeMap::new();
    for path in paths {
        for file in files(path)? {
            let file = file.canonicalize()?;
            if !result.contains_key(&file) {
                result.insert(file.clone(), hash(&file)?);
            }
        }
    }
    Ok(result)
}

/// Complete-content fingerprints with a bounded number of independent readers.
/// Errors are returned in canonical-path order after all readers have joined.
pub fn identities_parallel(paths: &[PathBuf], threads: usize) -> Result<BTreeMap<PathBuf, String>> {
    identities_using(paths, threads, &hash)
}

fn identities_using(
    paths: &[PathBuf],
    threads: usize,
    hash_file: &(impl Fn(&Path) -> Result<String> + Sync),
) -> Result<BTreeMap<PathBuf, String>> {
    ensure!(threads > 0, "validation threads must be positive");
    let mut unique = BTreeSet::new();
    for path in paths {
        for file in files(path)? {
            unique.insert(file.canonicalize()?);
        }
    }
    let files: Vec<_> = unique.into_iter().collect();
    let next = AtomicUsize::new(0);
    let mut results = std::thread::scope(|scope| -> Result<Vec<_>> {
        let handles: Vec<_> = (0..threads.min(files.len()))
            .map(|_| {
                scope.spawn(|| {
                    let mut result = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(path) = files.get(i) else { break };
                        result.push((i, hash_file(path)));
                    }
                    result
                })
            })
            .collect();
        let mut results = Vec::new();
        let mut panicked = false;
        for handle in handles {
            match handle.join() {
                Ok(part) => results.extend(part),
                Err(_) => panicked = true,
            }
        }
        ensure!(!panicked, "validation worker panicked");
        Ok(results)
    })?;
    results.sort_by_key(|(i, _)| *i);
    let mut output = BTreeMap::new();
    for (i, result) in results {
        output.insert(files[i].clone(), result?);
    }
    Ok(output)
}

pub fn write_json(path: &Path, value: &impl Serialize) -> Result<()> {
    let temp = path.with_extension("json.tmp");
    let mut data = serde_json::to_vec_pretty(value)?;
    data.push(b'\n');
    std::fs::write(&temp, data)?;
    std::fs::rename(temp, path)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn fixture() -> PathBuf {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let dir = std::env::temp_dir().join(format!(
            "esp-hashes-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        for n in 0..17 {
            std::fs::write(dir.join(format!("{n:02}.bin")), vec![n; 1031 + n as usize]).unwrap();
        }
        dir
    }

    #[test]
    fn parallel_digests_preserve_identity_and_detect_content_mutation() {
        let dir = fixture();
        let paths = vec![dir.clone(), dir.join("00.bin"), dir.join("./01.bin")];
        let expected = identities(&paths).unwrap();
        assert_eq!(expected.len(), 17);
        for cap in [1, 2, 4, 64] {
            assert_eq!(identities_parallel(&paths, cap).unwrap(), expected);
        }
        std::fs::write(dir.join("00.bin"), vec![42; 1031]).unwrap();
        assert_ne!(identities_parallel(&paths, 4).unwrap(), expected);
        assert!(identities_parallel(&paths, 0).is_err());
        assert!(identities_parallel(&[dir.join("missing")], 4).is_err());
        assert!(identities_parallel(&[], 4).unwrap().is_empty());
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn hash_workers_are_bounded_and_failures_are_deterministic() {
        let dir = fixture();
        let active = AtomicUsize::new(0);
        let peak = AtomicUsize::new(0);
        let completed = AtomicUsize::new(0);
        let result = identities_using(std::slice::from_ref(&dir), 3, &|path| {
            let n = active.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(n, Ordering::SeqCst);
            std::thread::sleep(std::time::Duration::from_millis(10));
            active.fetch_sub(1, Ordering::SeqCst);
            completed.fetch_add(1, Ordering::SeqCst);
            anyhow::bail!("failed {}", path.file_name().unwrap().to_string_lossy())
        });
        assert!(result.unwrap_err().to_string().contains("00.bin"));
        assert!(peak.load(Ordering::SeqCst) <= 3);
        assert!(peak.load(Ordering::SeqCst) > 1);
        assert_eq!(active.load(Ordering::SeqCst), 0);
        assert_eq!(completed.load(Ordering::SeqCst), 17);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn directory_symlink_is_rejected_for_parallel_validation() {
        let dir = fixture();
        std::os::unix::fs::symlink(dir.join("00.bin"), dir.join("link")).unwrap();
        assert!(identities_parallel(std::slice::from_ref(&dir), 2).is_err());
        std::fs::remove_dir_all(dir).unwrap();
    }
}
