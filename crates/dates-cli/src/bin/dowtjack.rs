use std::path::PathBuf;

use anyhow::Result;
use dates_core::workflow::run_dowtjack;

fn main() -> Result<()> {
    let mut input = None::<PathBuf>;
    let mut output = None::<PathBuf>;
    let mut mean = None::<f64>;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-i" => input = args.next().map(PathBuf::from),
            "-o" => output = args.next().map(PathBuf::from),
            "-m" => mean = args.next().map(|value| value.parse()).transpose()?,
            other => anyhow::bail!("unknown argument: {other}"),
        }
    }
    let summary = run_dowtjack(
        input
            .as_deref()
            .ok_or_else(|| anyhow::anyhow!("missing -i"))?,
        output
            .as_deref()
            .ok_or_else(|| anyhow::anyhow!("missing -o"))?,
        mean.ok_or_else(|| anyhow::anyhow!("missing -m"))?,
    )?;
    eprintln!("jackknife mean={} stderr={}", summary.mean, summary.std_err);
    Ok(())
}
