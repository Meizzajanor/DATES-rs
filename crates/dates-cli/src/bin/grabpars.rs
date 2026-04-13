use std::path::PathBuf;

use anyhow::{Result, bail};
use dates_core::workflow::grab_parameter;

fn main() -> Result<()> {
    let mut parname = None::<PathBuf>;
    let mut key = None::<String>;
    let mut verbose = false;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-p" => parname = args.next().map(PathBuf::from),
            "-x" => key = args.next(),
            "-V" => verbose = true,
            other => bail!("unknown argument: {other}"),
        }
    }
    let parname = parname.ok_or_else(|| anyhow::anyhow!("missing -p"))?;
    let key = key.ok_or_else(|| anyhow::anyhow!("missing -x"))?;
    let value = grab_parameter(&parname, &key)?;
    if verbose {
        eprintln!("{} {}", parname.display(), key);
    }
    println!("{value}");
    Ok(())
}
