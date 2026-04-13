# Architecture

## Core Modules
- `config`: legacy parameter-file parsing with variable substitution and typed extraction.
- `dataset`: text Eigenstrat loading and the typed data model used by the runtime surface.
- `corr`: sufficient-statistics correlation helpers and weighted jackknife routines translated from the original codebase.
- `fit`: exponential fitting for `dates_expfit` and the helper workflows.
- `plot`: generation of `.xtxt`, `.ps`, and `.pdf` artifacts from one plot model.
- `workflow`: legacy-compatible orchestration for helper commands.
- `dates`: runtime engine for DATES covariance accumulation and output generation.

## Translation Strategy
- Keep public behavior compatibility at the CLI and file level.
- Replace the C global-state style with typed Rust request/response flows.
- Document source provenance for every non-trivial numerical kernel.

