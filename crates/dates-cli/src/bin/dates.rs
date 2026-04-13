use std::path::PathBuf;

use anyhow::Result;
use dates_core::dates::run_dates;

fn main() -> Result<()> {
    let mut par = None::<PathBuf>;
    let mut verbose = false;
    let mut show_version = false;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-p" => par = args.next().map(PathBuf::from),
            "-V" => verbose = true,
            "-v" => show_version = true,
            other => anyhow::bail!("unknown argument: {other}"),
        }
    }
    if show_version {
        println!("version: 753");
        if par.is_none() {
            return Ok(());
        }
    }
    run_dates(
        par.as_deref()
            .ok_or_else(|| anyhow::anyhow!("missing -p"))?,
        verbose,
    )
}
