# Development validation — 2026-10-17

The package remains version 0.0.1.9000 pending remote checks and release review.

## Local checks

- R 4.6.1, macOS Tahoe 26.6.1, arm64.
- `devtools::test()`: 509 assertions, 0 failures/warnings/skips, with every
  optional adapter installed. Tests include independent Gaussian penalised
  refits, NCV/QNCV, Poisson and sampled discrete bam.
- `covr::package_coverage()`: 97.11% overall line coverage, 96.34% across
  folds/indices/nei/utils. `tools/coverage-gate.R` passes.
- `R CMD build .`: all five HTML vignettes built.
- `R CMD check --as-cran`: 0 errors, 0 warnings, 1 incoming NOTE (new submission,
  deliberate development version and not-yet-published site/workflow URLs).
  Includes PDF/HTML manuals, examples, tests and vignette rebuilding. No
  `--no-manual`, `--no-build-vignettes` or force-Suggests suppression was used.
- `pkgdown::build_site()`: reference pages and all five articles built. Installed
  rmarkdown/pkgdown emit Pandoc deprecation messages for `--mathml`/`--mathjax`;
  these are upstream command-line compatibility messages, not failed builds.
- `git diff --check`: clean.

Initial sandboxed network checks could not resolve remote hosts; the successful
full check and site build used approved network access. Initial source-check
issues were fixed before the successful run.

## Performance evidence

All scripts use R and are in `benchmarks/`. Heap figures below are R GC
high-water estimates, not operating-system peak resident memory. Timings depend
on hardware and dependencies; they are observations rather than guarantees.

| Fixture | Observations | Memberships | Construction | Packed output |
|---|---:|---:|---:|---:|
| Irregular unequal temporal series, sparse | 1,000,000 | 3,796,000 | ~0.30 s | 27.18 MB estimated |
| Temporal whole-series exclusions | 1,000,000 | 1,078,979,374 | ~0.26 s | 4.33 GB estimated; guard rejects |
| 1,000 fixed stations × 1,000 times | 1,000,000 | 14,972,012 | 14.37 s | 71,888,736 bytes measured |

The space-time fixture exported `nei` in 47.95 s, with export R heap high-water
976.9 MiB. Five boundary/interior focal sets and their packed endpoints matched
independent full-data scans. Moving-location, polygon and dense shared-station
fixtures also passed direct references.

`benchmarks/approximation.R` at n=120 found Gaussian fixed-basis NCV and exact
refits agreed to <=3.3e-13 relative discrepancy across the recorded widths.
Poisson NCV discrepancies were approximately 0.23–0.49% on that fixture. These
examples do not establish accuracy for other families or large omissions.

Use `fit$gcv.ubre` for the fitted criterion in these comparisons: mgcv can update
`fit$NCV` during covariance calculations with its default neighbourhood policy.

## Remaining release gates

Run the supplied GitHub OS/R matrix and minimal-install job, configure Pages and
Codecov, verify public URLs and uploads, review scientific examples, update the
version, and review CRAN submission metadata. See `tools/RELEASE.md`. At the time
of the local validation above, no remote push, release tag, CRAN submission or
deployment had been performed. Subsequent workflow and deployment results are
available in the repository's GitHub Actions history.
