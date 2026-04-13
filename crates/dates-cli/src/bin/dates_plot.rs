use anyhow::Result;
use dates_core::workflow::run_dates_plot;

fn main() -> Result<()> {
    let mut prefix = None::<String>;
    let mut data_col = 3usize;
    let mut low_cm = 0.45f64;
    let mut high_cm = 20.0f64;
    let mut step = 0.001f64;
    let mut seed = 77u64;
    let mut affine = false;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-i" => prefix = args.next(),
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
            "-s" => {
                step = args
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("missing -s value"))?
                    .parse()?
            }
            "-r" => {
                seed = args
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("missing -r value"))?
                    .parse()?
            }
            "-a" => affine = true,
            other => anyhow::bail!("unknown argument: {other}"),
        }
    }
    run_dates_plot(
        prefix
            .as_deref()
            .ok_or_else(|| anyhow::anyhow!("missing -i"))?,
        data_col,
        low_cm,
        high_cm,
        step,
        affine,
        seed,
    )?;
    Ok(())
}
