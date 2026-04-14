//! Legacy-compatible helper workflows for the installed command surface.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};

use crate::config::{LegacyParamFile, OutputPrefix, resolve_output_prefix, resolve_param_path};
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
    let output_base = std::env::current_dir()?;
    let prefix = OutputPrefix::resolve(prefix, &output_base)?;
    run_dates_plot_with_prefix(
        &prefix,
        data_col,
        low_cm,
        high_cm,
        step_morgans,
        affine,
        seed,
    )
}

fn run_dates_plot_with_prefix(
    prefix: &OutputPrefix,
    data_col: usize,
    low_cm: f64,
    high_cm: f64,
    step_morgans: f64,
    affine: bool,
    seed: u64,
) -> Result<FitResult> {
    let input = prefix.out_path();
    let output = prefix.fit_path();
    let log_path = prefix.expfit_log_path();
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
    let spec = PlotSpec::from_fit(
        format!("DATES: {}", prefix.raw()),
        &result.rows,
        low_cm,
        high_cm,
    );
    write_xtxt(&prefix.xtxt_path(), &spec, &output)?;
    write_ps(&prefix.ps_path(), &spec)?;
    write_pdf(&prefix.pdf_path(), &spec)?;
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
    run_dates_expfit_from_par_with_paths(&DatesExpfitRequest {
        par_path,
        data_col,
        low_cm,
        affine,
        seed,
        admix_override,
        output_dir: None,
        prefix_override: None,
    })
}

/// Request object for the `run_dates_expfit` helper.
#[derive(Clone, Copy, Debug)]
pub(crate) struct DatesExpfitRequest<'a> {
    pub par_path: &'a Path,
    pub data_col: usize,
    pub low_cm: f64,
    pub affine: bool,
    pub seed: u64,
    pub admix_override: Option<&'a str>,
    pub output_dir: Option<&'a Path>,
    pub prefix_override: Option<&'a str>,
}

/// Run `run_dates_expfit` with an explicit derived-output location.
pub(crate) fn run_dates_expfit_from_par_with_paths(
    request: &DatesExpfitRequest<'_>,
) -> Result<PathBuf> {
    let params = LegacyParamFile::load(request.par_path)?;
    reject_runmode_two(&params)?;
    let prefix = resolve_workflow_prefix(
        &params,
        request.admix_override,
        request.output_dir,
        request.prefix_override,
    )?;
    let binsize = params
        .get("binsize")
        .unwrap_or(".001")
        .parse::<f64>()
        .context("invalid binsize")?;
    let output = prefix.expfit_out_path();
    let fit_request = FitRequest {
        input: prefix.out_path(),
        output: Some(output.clone()),
        num_exp: 1,
        data_col: request.data_col,
        low_cm: request.low_cm,
        high_cm: f64::INFINITY,
        step_morgans: Some(binsize),
        add_x: 0.0,
        affine: request.affine,
        seed: request.seed,
    };
    let (_, stdout) = run_fit(&fit_request, "dates_expfit")?;
    let log_path = prefix.expfit_log_path();
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
    pub output_dir: Option<PathBuf>,
    pub prefix_override: Option<String>,
    pub affine: bool,
    pub seed: u64,
}

/// Run `dates_jackknife`.
pub fn run_dates_jackknife(request: &DatesJackknifeRequest) -> Result<JackknifeSummary> {
    let params = LegacyParamFile::load(&request.par_path)?;
    reject_runmode_two(&params)?;
    let prefix = resolve_workflow_prefix(
        &params,
        request.admix_override.as_deref(),
        request.output_dir.as_deref(),
        request.prefix_override.as_deref(),
    )?;
    let snp_name = request
        .snp_override
        .clone()
        .or_else(|| {
            params
                .get("snpname")
                .map(|path| resolve_param_path(&request.par_path, path))
        })
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
    let per_log = prefix.jackknife_log_path();
    let full_log = prefix.jackknife_full_log_path();
    let mut per_log_body = String::new();
    for chrom in 1..=num_chrom {
        let input = prefix.chrom_out_path(chrom);
        let output = prefix.fit_path();
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
        input: prefix.out_path(),
        output: Some(prefix.fit_path()),
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
    let jin_path = prefix.jin_path();
    fs::write(&jin_path, &jin)?;
    let jout_path = prefix.jout_path();
    let summary = dowtjack(&jin_path, &jout_path, overall_fit.mean_generations[0])?;
    run_dates_plot_with_prefix(
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

fn resolve_workflow_prefix(
    params: &LegacyParamFile,
    admix_override: Option<&str>,
    output_dir: Option<&Path>,
    prefix_override: Option<&str>,
) -> Result<OutputPrefix> {
    let output_base = match output_dir {
        Some(path) => path.to_path_buf(),
        None => std::env::current_dir()?,
    };
    Ok(match prefix_override {
        Some(prefix) => OutputPrefix::resolve(prefix, &output_base)?,
        None => resolve_output_prefix(params, admix_override, &output_base)?,
    })
}

fn reject_runmode_two(params: &LegacyParamFile) -> Result<()> {
    if params.get("runmode") == Some("2") {
        bail!("runmode 2 is not yet supported end-to-end in DATES-rs");
    }
    Ok(())
}
