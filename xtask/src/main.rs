use std::path::Path;

use anyhow::{Result, bail};
use dates_core::context::RunArtifact;

fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    match args.next().as_deref() {
        Some("verify") => verify(),
        Some(other) => bail!("unknown xtask command: {other}"),
        None => verify(),
    }
}

fn verify() -> Result<()> {
    let artifact = RunArtifact::load(Path::new("context/last-run.json"))?;
    if artifact.documentation.is_empty() {
        bail!("run artifact must list updated documentation");
    }
    if artifact.parity_state.trim().is_empty() {
        bail!("run artifact must record parity_state");
    }
    if !Path::new("AGENTS.md").exists() {
        bail!("AGENTS.md is required");
    }
    if !Path::new("docs/last-run.md").exists() {
        bail!("docs/last-run.md is required");
    }
    println!("xtask verify: ok");
    Ok(())
}
