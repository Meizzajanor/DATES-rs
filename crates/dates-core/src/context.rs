//! Machine-readable run context used by `xtask`.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

/// Recorded command execution in the run artifact.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RunCommand {
    pub program: String,
    pub args: Vec<String>,
}

/// Recorded verification check in the run artifact.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RunCheck {
    pub name: String,
    pub status: String,
    pub details: String,
}

/// Persistent last-run metadata that agents are required to keep current.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RunArtifact {
    pub date: String,
    pub status: String,
    pub touched_modules: Vec<String>,
    pub commands: Vec<RunCommand>,
    pub checks: Vec<RunCheck>,
    pub documentation: Vec<String>,
    pub parity_state: String,
    pub open_gaps: Vec<String>,
}

impl RunArtifact {
    /// Load the last-run artifact from disk.
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let raw = fs::read_to_string(path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        serde_json::from_str(&raw).with_context(|| format!("failed to parse {}", path.display()))
    }

    /// Persist the last-run artifact to disk.
    pub fn write(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        let raw = serde_json::to_string_pretty(self)?;
        fs::write(path, raw).with_context(|| format!("failed to write {}", path.display()))?;
        Ok(())
    }
}
