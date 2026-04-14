# Parity Status

## Supported Surface
- `grabpars`: implemented in Rust.
- `dowtjack`: implemented in Rust.
- `simpjack2`: implemented in Rust.
- `dates_expfit`: implemented in Rust.
- `dates_plot`: implemented in Rust.
- `run_dates_expfit`: implemented in Rust.
- `dates_jackknife`: implemented in Rust.
- `dates`: implemented for the supported text-Eigenstrat workflow used by the Rust fixture corpus.

## Verified Corpus
- `fixtures/toy` exercises the complete installed CLI surface and is covered by integration tests.
- Verification compares text artifacts exactly for the current golden corpus.
- Regression tests now also cover parameter-file-relative helper inputs, `numchrom > 22`, non-default plot ranges, filtered fit windows, helper-failure `cwd` safety, and explicit `runmode == 2` rejection.

## Known Gaps
- Packed-Eigenstrat input compatibility is not yet validated.
- `runmode == 2` is explicitly rejected because the legacy reduced-column output is not yet supported end-to-end by the helper workflow chain.
- Large-example parity against the original `family_packed.geno` workflow is not yet recorded in this repo.
- The current parity corpus is intentionally small and does not yet benchmark runtime or numerical stability on production-scale data.
