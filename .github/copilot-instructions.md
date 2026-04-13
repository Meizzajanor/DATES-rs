# DATES-rs Agent Context

## Project Goal
- `DATES-rs` is the pure-Rust runtime port of the supported DATES toolchain.
- Preserve the user-facing command surface for `dates`, `dates_expfit`, `grabpars`, `dowtjack`, `simpjack2`, `dates_jackknife`, `dates_plot`, and `run_dates_expfit`.
- Prefer behavioral and file compatibility with the original C/Perl workflow, with documented numeric tolerances where exact reproduction is not practical.

## Required Maintenance Rule
- After every implementation run, update this file if the working agreement changed.
- After every implementation run, update `docs/last-run.md` and `context/last-run.json`.
- The run report must list touched modules, commands executed, tests executed, parity state, documentation touched, and open gaps.
- If code and documentation disagree, update the documentation in the same run.

## Current Architecture
- `crates/dates-core`: reusable core logic, parsers, numerical routines, workflow helpers, plotting, and runtime orchestration.
- `crates/dates-cli`: the installed command surface, one Rust binary per legacy command name.
- `xtask`: repository verification helpers, including doc/context consistency checks.

## Guardrails
- Preserve ASCII text outputs unless the original workflow emits binary artifacts.
- Keep comments focused on scientific intent, compatibility assumptions, or non-obvious translation details.
- Record every unsupported legacy edge case in `docs/parity-status.md` until it is implemented.

