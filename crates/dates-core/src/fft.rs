//! FFT-backed correlation helpers.
//!
//! These mirror the original `src/fftsubs.c` logic closely enough to preserve
//! the positive-lag autocorrelation and cross-correlation behavior used by the
//! q-bin execution path in `dates`.

use num_complex::Complex64;
use rustfft::FftPlanner;

fn fft_len(m: usize) -> usize {
    if m == 0 {
        return 2;
    }
    m.next_power_of_two() * 2
}

/// Autocorrelation for non-negative lags.
pub fn auto_positive(a: &[f64], max_lag: usize) -> Vec<f64> {
    let n = fft_len(a.len());
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(n);
    let ifft = planner.plan_fft_inverse(n);

    let mut buffer = vec![Complex64::new(0.0, 0.0); n];
    for (index, value) in a.iter().enumerate() {
        buffer[index].re = *value;
    }
    fft.process(&mut buffer);
    for value in &mut buffer {
        *value *= value.conj();
    }
    ifft.process(&mut buffer);

    let scale = 1.0 / n as f64;
    let out_max = max_lag.min(n.saturating_sub(1));
    (0..=out_max).map(|lag| buffer[lag].re * scale).collect()
}

/// Cross-correlation for non-negative lags, `sum_i a[i] * b[i + lag]`.
pub fn cross_positive(a: &[f64], b: &[f64], max_lag: usize) -> Vec<f64> {
    let n = fft_len(a.len().max(b.len()));
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(n);
    let ifft = planner.plan_fft_inverse(n);

    let mut left = vec![Complex64::new(0.0, 0.0); n];
    let mut right = vec![Complex64::new(0.0, 0.0); n];
    for (index, value) in a.iter().enumerate() {
        left[index].re = *value;
    }
    for (index, value) in b.iter().enumerate() {
        right[n - index - 1].re = *value;
    }
    fft.process(&mut left);
    fft.process(&mut right);
    let mut product = left
        .iter()
        .zip(right.iter())
        .map(|(lhs, rhs)| lhs * rhs)
        .collect::<Vec<_>>();
    ifft.process(&mut product);

    let scale = 1.0 / n as f64;
    (0..=max_lag)
        .map(|lag| {
            let t = n.saturating_sub(lag + 1);
            product.get(t).map(|value| value.re * scale).unwrap_or(0.0)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn autocorrelation_matches_small_example() {
        let out = auto_positive(&[1.0, 2.0, 3.0], 2);
        assert_eq!(out.len(), 3);
        assert!((out[0] - 14.0).abs() < 1.0e-6);
        assert!((out[1] - 8.0).abs() < 1.0e-6);
        assert!((out[2] - 3.0).abs() < 1.0e-6);
    }

    #[test]
    fn cross_correlation_matches_small_example() {
        let out = cross_positive(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0], 2);
        assert_eq!(out.len(), 3);
        assert!((out[0] - 32.0).abs() < 1.0e-6);
        assert!((out[1] - 17.0).abs() < 1.0e-6);
        assert!((out[2] - 6.0).abs() < 1.0e-6);
    }
}
