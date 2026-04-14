//! Exponential fitting used by `dates_expfit` and the helper workflows.
//!
//! This module translates the active fitting path from `src/dates_expfit.c`,
//! `src/fitexp.c`, `src/gslfit.c`, and `src/regsubs.c`. The original code uses
//! random initialization followed by Nelder-Mead over decay bases in `(0, 1)`.
//! The Rust port keeps the same model shape and initialization style while
//! using a bounded local search over the same parameterization.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};
use nalgebra::{DMatrix, DVector};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::dataset::{FitRequest, FitResult, FitRow};

/// Run the DATES exponential fitting workflow on a whitespace-delimited table.
pub fn fit_request(request: &FitRequest) -> Result<FitResult> {
    let mut rows = load_fit_rows(&request.input, request.data_col)?;
    if rows.len() < 2 {
        bail!("fit input must contain at least two usable rows");
    }
    rows.retain(|(distance_cm, _)| {
        *distance_cm >= request.low_cm && *distance_cm <= request.high_cm
    });
    if rows.is_empty() {
        bail!("no data after applying fit range");
    }
    if rows.len() < 2 {
        bail!("fit range must retain at least two usable rows");
    }
    let mut step = request
        .step_morgans
        .unwrap_or_else(|| (rows[1].0 - rows[0].0) / 100.0);
    if step < 0.0 {
        bail!("step negative");
    }
    if step == 0.0 {
        step = 1.0;
    }
    let xbase = rows
        .iter()
        .map(|(distance, _)| *distance)
        .collect::<Vec<_>>();
    let observed = rows
        .iter()
        .map(|(_, value)| *value + request.add_x)
        .collect::<Vec<_>>();
    let init_iter = ((100.0 * (2.0f64).powi(request.num_exp as i32)).round() as usize).max(10);
    let basis = optimize_basis(
        &observed,
        request.num_exp.max(1),
        request.affine,
        request.seed,
        init_iter,
    )?;
    let scored = score_with_basis(&observed, &basis, request.affine)?;
    let halflives = basis
        .iter()
        .map(|value| -std::f64::consts::LN_2 / value.ln())
        .collect::<Vec<_>>();
    let mean_generations = halflives
        .iter()
        .map(|half| std::f64::consts::LN_2 / (half * step))
        .collect::<Vec<_>>();
    let mut coefficients = scored.coefficients.clone();
    for (index, mean) in mean_generations.iter().enumerate() {
        coefficients[index] *= (*mean * (xbase[0] * 0.01 - step)).exp();
    }
    let rows = xbase
        .iter()
        .zip(scored.fitted.iter())
        .zip(observed.iter())
        .map(|((distance_cm, fitted), observed)| FitRow {
            distance_cm: *distance_cm,
            observed: *observed,
            fitted: *fitted,
            residual: *observed - *fitted,
        })
        .collect::<Vec<_>>();
    if let Some(output) = &request.output {
        write_fit_rows(output, &rows)?;
    }
    Ok(FitResult {
        rows,
        error_sd: scored.mse.sqrt(),
        halflives,
        mean_generations,
        basis,
        coefficients,
        step_morgans: step,
    })
}

/// Render the legacy stdout summary for `dates_expfit`.
pub fn render_fit_stdout(program_name: &str, result: &FitResult, affine: bool) -> String {
    let mut out = String::new();
    out.push_str(&format!("{program_name} version: 200\n"));
    out.push_str(&format!("step (Morgans) :: {:12.6}\n", result.step_morgans));
    if affine {
        out.push_str(&format!(
            "fitting {} exponentials + affine\n",
            result.basis.len()
        ));
    } else {
        out.push_str(&format!("fitting {} exponentials\n", result.basis.len()));
    }
    out.push_str(&format!("error sd: {:12.6}\n", result.error_sd));
    out.push_str("halflife:");
    for value in &result.halflives {
        out.push_str(&format!(" {:12.6}", value));
    }
    out.push('\n');
    out.push_str("mean (generations):");
    for value in &result.mean_generations {
        out.push_str(&format!(" {:12.6}", value));
    }
    out.push('\n');
    out.push_str("coefficients:");
    for value in &result.coefficients {
        out.push_str(&format!(" {:12.6}", value));
    }
    out.push('\n');
    out.push_str("##end of run\n");
    out
}

fn load_fit_rows(path: &Path, data_col: usize) -> Result<Vec<(f64, f64)>> {
    let raw =
        fs::read_to_string(path).with_context(|| format!("failed to read {}", path.display()))?;
    let mut rows = Vec::new();
    for line in raw.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<_> = trimmed.split_whitespace().collect();
        if parts.len() <= data_col {
            bail!("fit row has too few columns: {trimmed}");
        }
        rows.push((parts[0].parse::<f64>()?, parts[data_col].parse::<f64>()?));
    }
    Ok(rows)
}

fn write_fit_rows(path: &Path, rows: &[FitRow]) -> Result<()> {
    let mut out = String::new();
    let label = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("fit");
    out.push_str(&format!("##fit: {label}\n"));
    for row in rows {
        out.push_str(&format!(
            "{:12.6} {:12.6} {:12.6} {:12.6}\n",
            row.distance_cm, row.observed, row.fitted, row.residual
        ));
    }
    fs::write(path, out).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn optimize_basis(
    observed: &[f64],
    num_exp: usize,
    affine: bool,
    seed: u64,
    init_iter: usize,
) -> Result<Vec<f64>> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed.max(1));
    let mut best = vec![0.5; num_exp];
    let mut best_score = f64::INFINITY;
    let mut mul = 1.0;
    let inner = init_iter / 10 + 10;
    for iter in 1..=init_iter {
        let mut trial = (0..num_exp)
            .map(|_| rng.random_range(0.000_001f64..0.999_999f64))
            .collect::<Vec<_>>();
        if best_score.is_finite() {
            for (value, base) in trial.iter_mut().zip(best.iter()) {
                *value = *value * mul + *base * (1.0 - mul);
            }
        }
        trial.sort_by(|lhs, rhs| lhs.partial_cmp(rhs).unwrap());
        let score = score_with_basis(observed, &trial, affine)?.mse;
        if score < best_score {
            best = trial;
            best_score = score;
        }
        if iter % inner == 0 {
            mul *= 0.5;
        }
    }
    let mut step = 0.05;
    while step > 1.0e-6 {
        let mut improved = false;
        for index in 0..best.len() {
            let current = best[index];
            for direction in [-1.0, 1.0] {
                let mut candidate = best.clone();
                candidate[index] = (current + direction * step).clamp(1.0e-9, 0.999_999_999);
                candidate.sort_by(|lhs, rhs| lhs.partial_cmp(rhs).unwrap());
                let score = score_with_basis(observed, &candidate, affine)?.mse;
                if score < best_score {
                    best = candidate;
                    best_score = score;
                    improved = true;
                }
            }
        }
        if !improved {
            step *= 0.5;
        }
    }
    Ok(best)
}

struct ScoredFit {
    mse: f64,
    fitted: Vec<f64>,
    coefficients: Vec<f64>,
}

fn score_with_basis(observed: &[f64], basis: &[f64], affine: bool) -> Result<ScoredFit> {
    let n = observed.len();
    let cols = basis.len() + usize::from(affine);
    let mut eq = Vec::with_capacity(n * cols);
    let mut powers = vec![1.0; cols];
    for _ in 0..n {
        for index in 0..basis.len() {
            powers[index] *= basis[index];
        }
        eq.extend_from_slice(&powers);
    }
    let x = DMatrix::from_row_slice(n, cols, &eq);
    let y = DVector::from_column_slice(observed);
    let coeffs = x
        .clone()
        .svd(true, true)
        .solve(&y, 1.0e-12)
        .map_err(|_| anyhow!("regression solve failed"))?;
    let fitted = x * coeffs.clone();
    let mse = observed
        .iter()
        .zip(fitted.iter())
        .map(|(left, right)| {
            let residual = left - right;
            residual * residual
        })
        .sum::<f64>()
        / observed.len() as f64;
    Ok(ScoredFit {
        mse,
        fitted: fitted.iter().copied().collect(),
        coefficients: coeffs.iter().copied().collect(),
    })
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;

    #[test]
    fn fit_request_runs_on_small_table() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("toy.out");
        fs::write(
            &input,
            "0.5 0.050\n0.6 0.040\n0.7 0.033\n0.8 0.026\n0.9 0.021\n",
        )
        .unwrap();
        let output = dir.path().join("toy.fit");
        let result = fit_request(&FitRequest {
            input,
            output: Some(output),
            num_exp: 1,
            data_col: 1,
            low_cm: 0.5,
            high_cm: 5.0,
            step_morgans: Some(0.001),
            add_x: 0.0,
            affine: true,
            seed: 77,
        })
        .unwrap();
        assert_eq!(result.basis.len(), 1);
        assert!(!result.rows.is_empty());
        let _ = PathBuf::from("unused");
    }

    #[test]
    fn fit_request_normalizes_zero_step_before_generation_conversion() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("toy.out");
        fs::write(
            &input,
            "0.5 0.050\n0.6 0.040\n0.7 0.033\n0.8 0.026\n0.9 0.021\n",
        )
        .unwrap();
        let result = fit_request(&FitRequest {
            input,
            output: None,
            num_exp: 1,
            data_col: 1,
            low_cm: 0.5,
            high_cm: 5.0,
            step_morgans: Some(0.0),
            add_x: 0.0,
            affine: true,
            seed: 77,
        })
        .unwrap();
        assert_eq!(result.step_morgans, 1.0);
        assert!(
            result
                .mean_generations
                .iter()
                .all(|value| value.is_finite())
        );
        assert!(result.coefficients.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn fit_request_rejects_one_row_after_filtering() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("toy.out");
        fs::write(&input, "0.5 0.050\n0.6 0.040\n0.7 0.033\n").unwrap();
        let err = fit_request(&FitRequest {
            input,
            output: None,
            num_exp: 1,
            data_col: 1,
            low_cm: 0.7,
            high_cm: 5.0,
            step_morgans: Some(0.001),
            add_x: 0.0,
            affine: true,
            seed: 77,
        })
        .unwrap_err()
        .to_string();
        assert!(
            err.contains("fit range must retain at least two usable rows"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn fit_request_rejects_empty_window_after_filtering() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("toy.out");
        fs::write(&input, "0.5 0.050\n0.6 0.040\n0.7 0.033\n").unwrap();
        let err = fit_request(&FitRequest {
            input,
            output: None,
            num_exp: 1,
            data_col: 1,
            low_cm: 1.0,
            high_cm: 5.0,
            step_morgans: Some(0.001),
            add_x: 0.0,
            affine: true,
            seed: 77,
        })
        .unwrap_err()
        .to_string();
        assert!(
            err.contains("no data after applying fit range"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn fit_request_infers_step_from_retained_rows() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("toy.out");
        fs::write(
            &input,
            "0.5 0.050\n0.6 0.040\n1.5 0.033\n1.7 0.026\n1.9 0.021\n",
        )
        .unwrap();
        let result = fit_request(&FitRequest {
            input,
            output: None,
            num_exp: 1,
            data_col: 1,
            low_cm: 1.5,
            high_cm: 5.0,
            step_morgans: None,
            add_x: 0.0,
            affine: true,
            seed: 77,
        })
        .unwrap();
        assert!((result.step_morgans - 0.002).abs() < 1.0e-12);
    }
}
