use std::process::Command;

fn esperanto() -> Command {
    Command::new(env!("CARGO_BIN_EXE_esperanto"))
}

#[test]
fn long_version_identifies_the_repository_and_exact_source_commit() {
    let output = esperanto().arg("--version").output().unwrap();
    assert!(output.status.success());
    let version = String::from_utf8(output.stdout).unwrap();
    assert!(version.lines().next().unwrap().contains("-dev."));
    assert!(version.contains("repository: https://github.com/intell-biol-lab/esperanto"));
    assert!(
        version.contains("worktree: clean")
            || version.contains("worktree: dirty")
            || version.contains("worktree: unknown")
    );
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let git = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(root)
        .output();
    if let Ok(git) = git {
        if git.status.success() {
            let sha = String::from_utf8(git.stdout).unwrap();
            assert!(version.contains(&format!("commit: {}", sha.trim())));
            return;
        }
    }
    assert!(version.contains("commit: unknown"));
}

#[test]
fn development_update_refuses_even_with_yes_and_an_unreachable_network() {
    for args in [vec!["update", "--yes"], vec!["update", "--check"]] {
        let output = esperanto()
            .args(args)
            .env("https_proxy", "http://127.0.0.1:9")
            .env("HTTPS_PROXY", "http://127.0.0.1:9")
            .env("http_proxy", "http://127.0.0.1:9")
            .env("HTTP_PROXY", "http://127.0.0.1:9")
            .env("NO_PROXY", "")
            .env("no_proxy", "")
            .output()
            .unwrap();
        assert!(!output.status.success());
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(
            stderr.contains("self-update is disabled for development builds"),
            "{stderr}"
        );
        assert!(stderr.contains("https://github.com/intell-biol-lab/esperanto"));
    }
}
