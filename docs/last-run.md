# Last Run

- Date: 2026-04-13
- Status: implementation verified
- Touched modules: workspace manifests, core library (config, dates, dataset), xtask manifest, docs, agent context
- Commands run:
  `cargo build --workspace`
  `cargo fmt`
  `cargo fmt --check`
  `cargo clippy --workspace`
  `cargo test --workspace`
- Tests run:
  `cargo test --workspace`
  8 CLI integration tests
  9 core unit tests
- Parity state: validated against the self-contained `fixtures/toy` golden corpus for the installed Rust CLI surface
- Documentation updated: `docs/last-run.md`
- Changes in this run:
  Removed unused `thiserror` dependency from `dates-core`
  Removed unused `clap` dependency from `dates-cli`
  Removed unused `serde`/`serde_json` from `xtask`
  Added `numchrom` field to `DatesParams` for configurable chromosome count
  Used `numchrom` for chromosome filtering, accumulator sizing, and jackknife output bounds
  Rejected `timeoffsetname` when `qbin > 0` with explicit error message
- Remaining gaps:
  packed Eigenstrat support
  large-example parity corpus
  production-scale performance validation
