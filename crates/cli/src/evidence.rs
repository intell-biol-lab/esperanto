//! Read-only content identities and atomic evidence publication.
use anyhow::{ensure, Result};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

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

pub fn write_json(path: &Path, value: &impl Serialize) -> Result<()> {
    let temp = path.with_extension("json.tmp");
    let mut data = serde_json::to_vec_pretty(value)?;
    data.push(b'\n');
    std::fs::write(&temp, data)?;
    std::fs::rename(temp, path)?;
    Ok(())
}
