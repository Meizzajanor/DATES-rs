use std::path::PathBuf;

use anyhow::Result;
use dates_core::workflow::run_simpjack2;

fn main() -> Result<()> {
    let mut input = Some(PathBuf::from("outf3"));
    let mut mean = None::<f64>;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-i" => input = args.next().map(PathBuf::from),
            "-m" => mean = args.next().map(|value| value.parse()).transpose()?,
            other => anyhow::bail!("unknown argument: {other}"),
        }
    }
    let (_, line) = run_simpjack2(
        input
            .as_deref()
            .ok_or_else(|| anyhow::anyhow!("missing -i"))?,
        mean,
    )?;
    println!("{line}");
    Ok(())
}
