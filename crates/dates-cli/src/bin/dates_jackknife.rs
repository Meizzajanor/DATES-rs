use std::path::PathBuf;

use anyhow::Result;
use dates_core::workflow::{DatesJackknifeRequest, run_dates_jackknife};

fn main() -> Result<()> {
    let mut par = None::<PathBuf>;
    let mut data_col = 3usize;
    let mut low_cm = 0.45f64;
    let mut high_cm = 20.0f64;
    let mut seed = 0u64;
    let mut map = None::<PathBuf>;
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
            "-h" => {
                high_cm = args
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("missing -h value"))?
                    .parse()?
            }
            "-r" => {
                seed = args
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("missing -r value"))?
                    .parse()?
            }
            "-m" => map = args.next().map(PathBuf::from),
            "-z" => admix = args.next(),
            "-a" => affine = true,
            other => anyhow::bail!("unknown argument: {other}"),
        }
    }
    let summary = run_dates_jackknife(&DatesJackknifeRequest {
        par_path: par.ok_or_else(|| anyhow::anyhow!("missing -p"))?,
        data_col,
        low_cm,
        high_cm,
        snp_override: map,
        admix_override: admix,
        output_dir: None,
        prefix_override: None,
        affine,
        seed,
    })?;
    println!("{:9.3}{:9.3}", summary.mean, summary.std_err);
    Ok(())
}
