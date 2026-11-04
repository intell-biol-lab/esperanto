use std::path::Path;
use std::process::Command;

fn git(root: &Path, args: &[&str]) -> Option<String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(args)
        .output()
        .ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_owned())
}

fn main() {
    let manifest = std::env::var_os("CARGO_MANIFEST_DIR").unwrap();
    let root = Path::new(&manifest).parent().unwrap().parent().unwrap();
    // Watch source and Git metadata, not build output. This also works with
    // linked worktrees where .git is a file and the index lives elsewhere.
    for path in [
        ".cargo",
        ".github",
        ".gitignore",
        "AGENTS.md",
        "CHANGELOG.md",
        "CONTRIBUTING.md",
        "Cargo.lock",
        "Cargo.toml",
        "LICENSE",
        "README.md",
        "bundle",
        "crates",
        "docs",
        "install.sh",
        "rust-toolchain.toml",
        "scripts",
    ] {
        println!("cargo:rerun-if-changed={}", root.join(path).display());
    }
    for path in ["HEAD", "index"] {
        if let Some(path) = git(root, &["rev-parse", "--git-path", path]) {
            println!("cargo:rerun-if-changed={}", root.join(path).display());
        }
    }
    if let Some(branch) = git(root, &["symbolic-ref", "-q", "HEAD"]) {
        if let Some(path) = git(root, &["rev-parse", "--git-path", &branch]) {
            println!("cargo:rerun-if-changed={}", root.join(path).display());
        }
    }
    if let Some(path) = git(root, &["rev-parse", "--git-path", "packed-refs"]) {
        let path = root.join(path);
        if path.exists() {
            println!("cargo:rerun-if-changed={}", path.display());
        }
    }
    // Do not attribute a source archive to an unrelated enclosing repository.
    let is_checkout = git(root, &["rev-parse", "--show-toplevel"])
        .is_some_and(|top| Path::new(&top).canonicalize().ok() == root.canonicalize().ok());
    let sha = is_checkout
        .then(|| git(root, &["rev-parse", "HEAD"]))
        .flatten()
        .unwrap_or_else(|| "unknown".to_owned());
    let state = if is_checkout {
        match git(root, &["status", "--porcelain", "--untracked-files=normal"]) {
            Some(status) if status.is_empty() => "clean",
            Some(_) => "dirty",
            None => "unknown",
        }
    } else {
        "unknown"
    };
    println!("cargo:rustc-env=ESP_BUILD_GIT_SHA={sha}");
    println!("cargo:rustc-env=ESP_BUILD_WORKTREE={state}");
}
