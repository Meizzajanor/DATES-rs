# Last Run

- Date: 2026-04-13
- Status: implementation verified
- Touched modules: `dates-core` (`config`, `dataset`, `dates`, `fit`, `plot`, `workflow`), `dates-cli` (tests/fixtures), docs
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
  Replaced the internal `run_dates_expfit_from_par_with_paths` argument list with a dedicated `DatesExpfitRequest` struct to satisfy `clippy::too_many_arguments`
  Updated the runtime call site to pass the request object without changing the public CLI-facing behavior
  Introduced `OutputPrefix` handling so helper-generated outputs consistently derive related artifact paths from the configured prefix
  Resolved parameter-file-relative input and output paths in the helper workflow so `.par` execution matches the legacy toolchain more closely when invoked outside the parameter file directory
  Added SNP ordering validation to reject incompatible or unsorted marker input earlier in the workflow instead of silently proceeding with invalid assumptions
  Adjusted plot range handling so non-default fit and plotting windows propagate correctly through the runtime and plotting helpers
  Expanded regression coverage for helper path resolution, `numchrom > 22`, filtered fit windows, non-default plot ranges, helper failure without `cwd` leakage, and explicit `runmode == 2` rejection
  Made `OutputPrefix::resolve` and `from_resolved` return `Result` to surface a clear error when the prefix path has no valid file name component
  Renamed the local `FitRequest` binding in `run_dates_expfit_from_par_with_paths` from `request` to `fit_request` to eliminate shadowing of the outer `DatesExpfitRequest` parameter
- Remaining gaps:
  packed Eigenstrat support
  `runmode == 2` end-to-end compatibility
  large-example parity corpus
  production-scale runtime and numerical validation
