# Last Run

- Date: 2026-04-13
- Status: implementation verified
- Touched modules: `dates-core` (`dates`, `workflow`), docs
- Commands run:
  `cargo fmt`
  `cargo test --workspace`
  `cargo clippy`
  `cargo run -p xtask -- verify`
- Tests run:
  `cargo test --workspace`
  37 total tests passed
  `cargo clippy`
  workspace lint pass with no warnings
- Parity state: unchanged from the prior verified state; helper workflow behavior remains validated against the self-contained `fixtures/toy` golden corpus and targeted regressions for helper path resolution, `numchrom > 22`, filtered fit windows, non-default plot ranges, helper failure without `cwd` leakage, and explicit `runmode == 2` rejection
- Documentation updated: `docs/last-run.md`, `context/last-run.json`
- Changes in this run:
  Replaced the internal `run_dates_expfit_from_par_with_paths` argument list with a dedicated request struct to satisfy `clippy::too_many_arguments`
  Updated the runtime call site to pass the request object without changing the public CLI-facing behavior
- Remaining gaps:
  packed Eigenstrat support
  `runmode == 2` end-to-end compatibility
  large-example parity corpus
  production-scale runtime and numerical validation
