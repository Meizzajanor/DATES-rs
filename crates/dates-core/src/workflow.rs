//! Legacy-compatible helper workflows for the installed command surface.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow};

use crate::config::{LegacyParamFile, infer_output_prefix};
use crate::dataset::{FitRequest, FitResult, JackknifeSummary, chromosome_counts};
use crate::fit::{fit_request, render_fit_stdout};
use crate::jackknife::{dowtjack, simpjack2};
use crate::plot::{PlotSpec, write_pdf, write_ps, write_xtxt};

/// Return the resolved value for `key`, or `NOTFOUND` when the key is absent.
pub fn grab_parameter(par_path: &Path, key: &str) -> Result<String> {
    let params = LegacyParamFile::load(par_path)?;
    let normalized = key.trim().trim_end_matches(':');
    Ok(params.get(normalized).unwrap_or("NOTFOUND").to_owned())
}

/// Run the `dates_expfit` workflow on one input file.
pub fn run_fit(request: &FitRequest, program_name: &str) -> Result<(FitResult, String)> {
    let result = fit_request(request)?;
    let stdout = render_fit_stdout(program_name, &result, request.affine);
    Ok((result, stdout))
}

/// Run the `dates_plot` helper workflow.
pub fn run_dates_plot(
    prefix: &str,
    data_col: usize,
    low_cm: f64,
    high_cm: f64,
    step_morgans: f64,
    affine: bool,
    seed: u64,
) -> Result<FitResult> {
    let input = PathBuf::from(format!("{prefix}.out"));
    let output = PathBuf::from(format!("{prefix}.fit"));
    let log_path = PathBuf::from(format!("{prefix}_expfit.log"));
    let request = FitRequest {
        input: input.clone(),
        output: Some(output.clone()),
        num_exp: 1,
        data_col,
        low_cm,
        high_cm,
        step_morgans: Some(step_morgans),
        add_x: 0.0,
        affine,
        seed,
    };
    let (result, stdout) = run_fit(&request, "dates_expfit")?;
    fs::write(&log_path, stdout)
        .with_context(|| format!("failed to write {}", log_path.display()))?;
    let spec = PlotSpec::from_fit(format!("DATES: {prefix}"), &result.rows);
    write_xtxt(Path::new(&format!("{prefix}.xtxt")), &spec, &output)?;
    write_ps(Path::new(&format!("{prefix}.ps")), &spec)?;
    write_pdf(Path::new(&format!("{prefix}.pdf")), &spec)?;
    Ok(result)
}

/// Run `run_dates_expfit`, the non-jackknife helper.
pub fn run_dates_expfit_from_par(
    par_path: &Path,
    data_col: usize,
    low_cm: f64,
    affine: bool,
    seed: u64,
    admix_override: Option<&str>,
) -> Result<PathBuf> {
    let params = LegacyParamFile::load(par_path)?;
    let prefix = infer_output_prefix(&params, admix_override)?;
    let binsize = params
        .get("binsize")
        .unwrap_or(".001")
        .parse::<f64>()
        .context("invalid binsize")?;
    let output = PathBuf::from(format!("{prefix}:expfit.out"));
    let request = FitRequest {
        input: PathBuf::from(format!("{prefix}.out")),
        output: Some(output.clone()),
        num_exp: 1,
        data_col,
        low_cm,
        high_cm: f64::INFINITY,
        step_morgans: Some(binsize),
        add_x: 0.0,
        affine,
        seed,
    };
    let (_, stdout) = run_fit(&request, "dates_expfit")?;
    let log_path = PathBuf::from(format!("{prefix}_expfit.log"));
    fs::write(log_path, stdout)?;
    Ok(output)
}

/// Request object for the `dates_jackknife` helper.
#[derive(Clone, Debug)]
pub struct DatesJackknifeRequest {
    pub par_path: PathBuf,
    pub data_col: usize,
    pub low_cm: f64,
    pub high_cm: f64,
    pub snp_override: Option<PathBuf>,
    pub admix_override: Option<String>,
    pub affine: bool,
    pub seed: u64,
}

/// Run `dates_jackknife`.
pub fn run_dates_jackknife(request: &DatesJackknifeRequest) -> Result<JackknifeSummary> {
    let params = LegacyParamFile::load(&request.par_path)?;
    let prefix = infer_output_prefix(&params, request.admix_override.as_deref())?;
    let snp_name = request
        .snp_override
        .clone()
        .or_else(|| params.get("snpname").map(PathBuf::from))
        .ok_or_else(|| anyhow!("missing snpname"))?;
    let binsize = params
        .get("binsize")
        .unwrap_or(".001")
        .parse::<f64>()
        .context("invalid binsize")?;
    let num_chrom = params
        .get("numchrom")
        .unwrap_or("22")
        .parse::<i32>()
        .context("invalid numchrom")?;
    let counts = chromosome_counts(&snp_name)?;
    let mut jin = String::new();
    let per_log = PathBuf::from(format!("expfit_{prefix}.log"));
    let full_log = PathBuf::from(format!("expfit_{prefix}.flog"));
    let mut per_log_body = String::new();
    for chrom in 1..=num_chrom {
        let input = PathBuf::from(format!("{prefix}.out:{chrom}"));
        let output = PathBuf::from(format!("{prefix}.fit"));
        let fit_request = FitRequest {
            input,
            output: Some(output),
            num_exp: 1,
            data_col: request.data_col,
            low_cm: request.low_cm,
            high_cm: request.high_cm,
            step_morgans: Some(binsize),
            add_x: 0.0,
            affine: request.affine,
            seed: request.seed + chrom as u64,
        };
        let (fit, stdout) = run_fit(&fit_request, "dates_expfit")?;
        per_log_body.push_str(&stdout);
        jin.push_str(&format!(
            "{:3}  {:12.6}  {:12.6}\n",
            chrom,
            counts.get(&chrom).copied().unwrap_or_default() as f64,
            fit.mean_generations[0]
        ));
    }
    fs::write(&per_log, per_log_body)?;
    let overall_request = FitRequest {
        input: PathBuf::from(format!("{prefix}.out")),
        output: Some(PathBuf::from(format!("{prefix}.fit"))),
        num_exp: 1,
        data_col: request.data_col,
        low_cm: request.low_cm,
        high_cm: request.high_cm,
        step_morgans: Some(binsize),
        add_x: 0.0,
        affine: request.affine,
        seed: request.seed,
    };
    let (overall_fit, stdout) = run_fit(&overall_request, "dates_expfit")?;
    fs::write(&full_log, stdout)?;
    let jin_path = PathBuf::from(format!("{prefix}.jin"));
    fs::write(&jin_path, &jin)?;
    let jout_path = PathBuf::from(format!("{prefix}.jout"));
    let summary = dowtjack(&jin_path, &jout_path, overall_fit.mean_generations[0])?;
    run_dates_plot(
        &prefix,
        request.data_col,
        request.low_cm,
        request.high_cm,
        binsize,
        request.affine,
        request.seed,
    )?;
    Ok(summary)
}

/// Run `dowtjack`.
pub fn run_dowtjack(input: &Path, output: &Path, mean: f64) -> Result<JackknifeSummary> {
    dowtjack(input, output, mean)
}

/// Run `simpjack2`.
pub fn run_simpjack2(input: &Path, mean: Option<f64>) -> Result<(JackknifeSummary, String)> {
    simpjack2(input, mean)
}
