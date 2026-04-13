# DATES-rs

`DATES-rs` is a pure-Rust port of the supported DATES runtime toolchain. It is structured as a Cargo workspace with a reusable core crate, CLI binaries that preserve the legacy command names, and an `xtask` verifier for parity and documentation hygiene.

## Workspace
- `crates/dates-core`: config parsing, Eigenstrat I/O, covariance kernels, fitting, jackknife utilities, plotting, and workflow orchestration.
- `crates/dates-cli`: binaries named `dates`, `dates_expfit`, `grabpars`, `dowtjack`, `simpjack2`, `dates_jackknife`, `dates_plot`, and `run_dates_expfit`.
- `xtask`: repository maintenance and verification commands.

## Verification
```bash
cargo fmt --check
cargo clippy --workspace --all-targets
cargo test --workspace
cargo doc --workspace --no-deps
cargo run -p xtask -- verify
```

## Compatibility Notes
- The Rust port preserves the supported legacy CLI surface and parameter-file syntax.
- Documentation, parity status, and last-run context live under `docs/` and `context/`.
- Current validation scope is recorded in [docs/parity-status.md](docs/parity-status.md).

## Fixture Corpus
- `fixtures/toy` contains a self-contained text-Eigenstrat corpus and golden outputs for the full Rust CLI surface.
- The current integration suite validates all eight installed binaries against that corpus.
