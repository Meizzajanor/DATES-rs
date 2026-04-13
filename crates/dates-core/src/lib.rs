//! Pure-Rust runtime support for the DATES toolchain.
//!
//! This crate intentionally mirrors the supported execution surface of the
//! original C/Perl distribution while exposing a reusable Rust API for testing,
//! documentation, and future maintenance.

pub mod config;
pub mod context;
pub mod corr;
pub mod dataset;
pub mod dates;
pub mod fft;
pub mod fit;
pub mod jackknife;
pub mod plot;
pub mod workflow;

pub use config::{DatesParams, LegacyParamFile};
pub use context::{RunArtifact, RunCheck, RunCommand};
pub use corr::{Corr, weighted_jackknife};
pub use dataset::{
    AdmixJob, CovarianceBin, Dataset, FitRequest, FitResult, FitRow, Individual, JackknifeSample,
    JackknifeSummary, Snp,
};
