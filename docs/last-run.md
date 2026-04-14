# Last Run

- Date: 2026-04-13
- Status: implementation verified
- Touched modules: dates-core (dataset, dates, config)
- Commands run:
  `cargo build --workspace`
  `cargo test --workspace`
- Tests run:
  `cargo test --workspace`
  8 CLI integration tests
  11 core unit tests
- Parity state: validated against the self-contained `fixtures/toy` golden corpus for the installed Rust CLI surface
- Documentation updated: `docs/last-run.md`, `context/last-run.json`
- Changes in this run:
  Replaced wholesale `fs::read` + `String::from_utf8` in `load_text_genotypes` with `BufReader` line-by-line streaming to reduce peak memory on large `.geno` files
  Refactored `build_present_snps` and `build_weighted_values` to accept pre-allocated `&mut Vec` buffers, lifted allocations outside per-individual loops in `run_direct_mode` and `run_qbin_mode`
  Eliminated intermediate `Vec` allocations inside `build_weighted_values` by computing dot products inline
  Removed now-unused `dot()` helper
  Added explanatory comment on the `for _ in 0..8` circuit-breaker loop in `resolve_entries`
- Remaining gaps:
  packed Eigenstrat support
  large-example parity corpus
  CLI argument parsing robustness (optional)
