use std::path::PathBuf;

use anyhow::Result;
use dates_core::dataset::FitRequest;
use dates_core::workflow::run_fit;

fn main() -> Result<()> {
    let mut input = None::<PathBuf>;
    let mut output = None::<PathBuf>;
    let mut num_exp = 2usize;
    let mut data_col = 1usize;
    let mut low_cm = -1.0e20;
    let mut high_cm = 1.0e20;
    let mut step = None::<f64>;
    let mut add_x = 0.0f64;
    let mut affine = false;
    let mut seed = 0u64;
    let mut help = false;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-i" => input = args.next().map(PathBuf::from),
            "-o" => output = args.next().map(PathBuf::from),
            "-n" => {
                num_exp = args
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("missing -n value"))?
                    .parse()?
            }
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
                step = Some(
                    args.next()
                        .ok_or_else(|| anyhow::anyhow!("missing -s value"))?
                        .parse()?,
                )
            }
            "-x" => {
                add_x = args
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("missing -x value"))?
                    .parse()?
            }
            "-a" => affine = true,
            "-r" => {
                seed = args
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("missing -r value"))?
                    .parse()?
            }
            "-m" => help = true,
            "-V" => {}
            other => anyhow::bail!("unknown argument: {other}"),
        }
    }
    if help || input.is_none() {
        println!("expfit:");
        println!(" -i iname   ## input");
        println!(" -o oname   ## output");
        println!(" -n numexp  ## number of exponentials");
        println!(" -c col     ## data column  (0 is xval)");
        println!(" -l loval   ## lowest x value");
        println!(" -h hival   ## highest x value");
        println!(" -s stepsize   ## step size (Morgans)");
        println!(" -x val     ## value to add to x (deprecated)");
        println!(" -r ran     ## seed for random generator");
        println!(" -a         ## affine mode (add constant)");
        println!(" -V         ## verbose mode");
        println!(" -m         ## print help menu and quit");
        return Ok(());
    }
    let request = FitRequest {
        input: input.unwrap(),
        output,
        num_exp,
        data_col,
        low_cm,
        high_cm,
        step_morgans: step,
        add_x,
        affine,
        seed,
    };
    let (_, stdout) = run_fit(&request, "dates_expfit")?;
    print!("{stdout}");
    Ok(())
}
