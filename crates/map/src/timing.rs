//! Opt-in wall-clock profiling at stage and chunk boundaries.
//!
//! `ESP_PROFILE_STAGES=1` enables stderr records. The setting is cached on the
//! first timer start and applies for the rest of the process; changing the
//! environment afterward has no effect. Disabled timers never read the clock.
//! Keep guards outside per-read loops so both checks and output stay bounded
//! by the number of stages/chunks. Nested timings are inclusive, not additive.

use std::io::Write;
use std::sync::OnceLock;
use std::time::Instant;

static ENABLED: OnceLock<bool> = OnceLock::new();

/// Logs elapsed wall seconds on drop, including early returns and unwinding.
#[must_use = "keep the timer guard alive until the stage ends"]
pub struct StageTimer {
    name: &'static str,
    started: Option<Instant>,
}

impl StageTimer {
    /// Start a named stage if process-scoped profiling is enabled.
    pub fn start(name: &'static str) -> Self {
        let enabled = *ENABLED
            .get_or_init(|| std::env::var_os("ESP_PROFILE_STAGES").is_some_and(|v| v == "1"));
        Self {
            name,
            started: enabled.then(Instant::now),
        }
    }
}

impl Drop for StageTimer {
    fn drop(&mut self) {
        if let Some(started) = self.started {
            let seconds = started.elapsed().as_secs_f64();
            // Profiling must not turn a stderr failure into a mapping failure.
            let _ = writeln!(
                std::io::stderr().lock(),
                "[esp-profile] {} {:.6}",
                self.name,
                seconds
            );
        }
    }
}
