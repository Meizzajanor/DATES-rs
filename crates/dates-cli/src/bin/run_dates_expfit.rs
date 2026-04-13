use std::path::PathBuf;

use anyhow::Result;
use dates_core::workflow::run_dates_expfit_from_par;

fn main() -> Result<()> {
    let mut par = None::<PathBuf>;
    let mut data_col = 3usize;
    let mut low_cm = 0.45f64;
    let mut seed = 0u64;
    let mut admix = None::<String>;
    let mut affine = false;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-p" => par = args.next().map(PathBuf::from),
            "-c" => {
                data_col = args
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("missing -c value"))?
                    .parse()?
            }
            "-l" => {
                low_cm = args
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("missing -l value"))?
                    .parse()?
            }
            "-r" => {
                seed = args
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("missing -r value"))?
                    .parse()?
            }
            "-z" => admix = args.next(),
            "-a" => affine = true,
            other => anyhow::bail!("unknown argument: {other}"),
        }
    }
    let output = run_dates_expfit_from_par(
        par.as_deref()
            .ok_or_else(|| anyhow::anyhow!("missing -p"))?,
        data_col,
        low_cm,
        affine,
        seed,
        admix.as_deref(),
    )?;
    println!("{}", output.display());
    Ok(())
}
