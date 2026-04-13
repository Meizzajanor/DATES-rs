//! Helper programs for weighted jackknife summaries.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::corr::weighted_jackknife;
use crate::dataset::{JackknifeSample, JackknifeSummary};

/// Reimplementation of `dowtjack`.
pub fn dowtjack(input: &Path, output: &Path, mean: f64) -> Result<JackknifeSummary> {
    let samples = parse_weighted_samples(input)?;
    let weights = samples
        .iter()
        .map(|sample| sample.weight)
        .collect::<Vec<_>>();
    let estimates = samples
        .iter()
        .map(|sample| sample.estimate)
        .collect::<Vec<_>>();
    let (jack_mean, std_err) = weighted_jackknife(mean, &estimates, &weights)?;
    fs::write(output, format!("{jack_mean:9.3}{std_err:9.3}\n"))
        .with_context(|| format!("failed to write {}", output.display()))?;
    Ok(JackknifeSummary {
        mean: jack_mean,
        std_err,
    })
}

/// Reimplementation of `simpjack2`.
pub fn simpjack2(input: &Path, override_mean: Option<f64>) -> Result<(JackknifeSummary, String)> {
    let rows = parse_numeric_rows(input)?;
    let mut global_mean = override_mean;
    let mut samples = Vec::new();
    for row in rows {
        if row.len() < 3 {
            continue;
        }
        let block = row[0] as i32;
        if block == 0 && global_mean.is_none() {
            global_mean = Some(row[2]);
            continue;
        }
        samples.push(JackknifeSample {
            block,
            weight: row[1],
            estimate: row[2],
        });
    }
    let mean = global_mean.ok_or_else(|| anyhow::anyhow!("simpjack2 requires a global mean"))?;
    let weights = samples
        .iter()
        .map(|sample| sample.weight)
        .collect::<Vec<_>>();
    let estimates = samples
        .iter()
        .map(|sample| sample.estimate)
        .collect::<Vec<_>>();
    let (jack_mean, std_err) = weighted_jackknife(mean, &estimates, &weights)?;
    let line = format!(
        "## simpjack2: {} {:12.6} {:12.6} {:12.6} {:9.3}",
        input.display(),
        jack_mean,
        mean,
        std_err,
        jack_mean / std_err.max(1.0e-20)
    );
    Ok((
        JackknifeSummary {
            mean: jack_mean,
            std_err,
        },
        line,
    ))
}

fn parse_weighted_samples(path: &Path) -> Result<Vec<JackknifeSample>> {
    let rows = parse_numeric_rows(path)?;
    let mut samples = Vec::new();
    for row in rows {
        if row.len() < 3 {
            continue;
        }
        samples.push(JackknifeSample {
            block: row[0] as i32,
            weight: row[1],
            estimate: row[2],
        });
    }
    if samples.len() < 2 {
        bail!("not enough weighted jackknife samples");
    }
    Ok(samples)
}

fn parse_numeric_rows(path: &Path) -> Result<Vec<Vec<f64>>> {
    let raw =
        fs::read_to_string(path).with_context(|| format!("failed to read {}", path.display()))?;
    let mut rows = Vec::new();
    for line in raw.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let row = trimmed
            .split_whitespace()
            .map(str::parse::<f64>)
            .collect::<std::result::Result<Vec<_>, _>>()?;
        rows.push(row);
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dowtjack_and_simpjack2_work() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("jin.txt");
        fs::write(&input, "1 10 12.0\n2 9 11.5\n3 11 12.4\n").unwrap();
        let output = dir.path().join("jout.txt");
        let summary = dowtjack(&input, &output, 12.0).unwrap();
        assert!(summary.mean.is_finite());
        let simp_input = dir.path().join("simp.txt");
        fs::write(&simp_input, "0 0 12.0\n1 10 12.0\n2 9 11.5\n3 11 12.4\n").unwrap();
        let (_, line) = simpjack2(&simp_input, None).unwrap();
        assert!(line.contains("simpjack2"));
    }
}
