//! Correlation and weighted jackknife helpers.
//!
//! These routines are direct Rust translations of the sufficient-statistics
//! helpers in `src/ldsubs.c`, `src/qpsubs.c`, and `src/nicksrc/statsubs.c`.

use anyhow::{Result, bail};

/// Sufficient statistics for the covariance and correlation summaries written by DATES.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Corr {
    pub s0: f64,
    pub s1: f64,
    pub s2: f64,
    pub s11: f64,
    pub s12: f64,
    pub s22: f64,
    pub m1: f64,
    pub m2: f64,
    pub v11: f64,
    pub v12: f64,
    pub v22: f64,
    pub corr: f64,
    pub z: f64,
}

impl Corr {
    /// Reset all sufficient statistics.
    pub fn clear(&mut self) {
        *self = Self::default();
    }

    /// Add one observation pair.
    pub fn add(&mut self, x1: f64, x2: f64) {
        self.s0 += 1.0;
        self.s1 += x1;
        self.s2 += x2;
        self.s11 += x1 * x1;
        self.s12 += x1 * x2;
        self.s22 += x2 * x2;
    }

    /// Add pre-aggregated sufficient statistics.
    pub fn add_summary(&mut self, x0: f64, x1: f64, x2: f64, x12: f64, x11: f64, x22: f64) {
        self.s0 += x0;
        self.s1 += x1;
        self.s2 += x2;
        self.s12 += x12;
        self.s11 += x11;
        self.s22 += x22;
    }

    /// Combine two correlation accumulators by addition.
    pub fn plus(self, other: Self) -> Self {
        Self {
            s0: self.s0 + other.s0,
            s1: self.s1 + other.s1,
            s2: self.s2 + other.s2,
            s11: self.s11 + other.s11,
            s12: self.s12 + other.s12,
            s22: self.s22 + other.s22,
            ..Self::default()
        }
    }

    /// Combine two correlation accumulators by subtraction.
    pub fn minus(self, other: Self) -> Result<Self> {
        let out = Self {
            s0: self.s0 - other.s0,
            s1: self.s1 - other.s1,
            s2: self.s2 - other.s2,
            s11: self.s11 - other.s11,
            s12: self.s12 - other.s12,
            s22: self.s22 - other.s22,
            ..Self::default()
        };
        if out.s0 < -1.0e-6 {
            bail!("invalid negative sufficient-statistics count");
        }
        Ok(out)
    }

    /// Recompute mean, variance, correlation, and Z-score.
    pub fn calculate(&mut self, zero_centered: bool, z_transform: bool) -> Result<bool> {
        self.corr = 0.0;
        self.z = 0.0;
        if self.s0 < 0.5 {
            return Ok(false);
        }
        let yn = self.s0;
        self.m1 = self.s1 / self.s0;
        self.m2 = self.s2 / self.s0;
        let (m1, m2) = if zero_centered {
            (0.0, 0.0)
        } else {
            (self.m1, self.m2)
        };
        self.v11 = (self.s11 - yn * m1 * m1) / yn;
        self.v12 = (self.s12 - yn * m1 * m2) / yn;
        self.v22 = (self.s22 - yn * m2 * m2) / yn;
        self.corr = self.v12 / (self.v11 * self.v22 + 1.0e-20).sqrt();
        self.z = yn.sqrt() * self.corr;
        if z_transform {
            if yn < 4.0 {
                return Ok(false);
            }
            let clipped = self.corr.clamp(-0.9, 0.9);
            let r = 0.5 * ((1.0 + clipped) / (1.0 - clipped)).ln();
            self.z = (yn - 3.0).sqrt() * r;
        }
        Ok(true)
    }
}

/// Weighted jackknife estimate and standard error.
pub fn weighted_jackknife(
    mean: f64,
    leave_one_out_means: &[f64],
    weights: &[f64],
) -> Result<(f64, f64)> {
    let filtered: Vec<(f64, f64)> = leave_one_out_means
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .filter(|(_, weight)| *weight >= 1.0e-6)
        .collect();
    if filtered.len() <= 1 {
        bail!("weighted_jackknife requires at least two blocks");
    }
    let yn: f64 = filtered.iter().map(|(_, weight)| weight).sum();
    let jackest = filtered.iter().map(|(value, _)| mean - value).sum::<f64>()
        + filtered
            .iter()
            .map(|(value, weight)| weight * value)
            .sum::<f64>()
            / yn;
    let mut yvar = 0.0;
    for (value, weight) in &filtered {
        let hh = yn / weight;
        let w1 = hh - 1.0;
        let mut xtau = hh * mean - w1 * value;
        xtau -= jackest;
        yvar += (xtau * xtau) / w1;
    }
    yvar /= filtered.len() as f64;
    Ok((jackest, yvar.max(0.0).sqrt()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weighted_jackknife_handles_simple_input() {
        let (est, sig) = weighted_jackknife(10.0, &[9.0, 11.0, 10.5], &[5.0, 4.0, 6.0]).unwrap();
        assert!(est.is_finite());
        assert!(sig.is_finite());
    }
}
