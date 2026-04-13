# C-to-Rust Traceability

## Supported Runtime Surface
- `src/dates.c` -> `dates_core::dates`
- `src/dates_expfit.c`, `src/fitexp.c`, `src/gslfit.c`, `src/regsubs.c` -> `dates_core::fit`
- `src/dowtjack.c`, `src/simpjack2.c`, `src/qpsubs.c`, `src/nicksrc/statsubs.c` -> `dates_core::jackknife` and `dates_core::corr`
- `src/grabpars.c`, `src/perlsrc/dates_jackknife`, `src/perlsrc/dates_plot`, `src/perlsrc/run_dates_expfit` -> `dates_core::workflow`
- `src/ldsubs.c` -> `dates_core::corr`
- `src/fftsubs.c` -> `dates_core::fft`
- `src/mcio.c`, `src/admutils.c`, `src/egsubs.c` -> `dates_core::dataset` and `dates_core::config`

