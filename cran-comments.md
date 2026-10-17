# 0.1.0 release candidate — not submitted

## Local test environment

* macOS Tahoe 26.6.1, arm64, R 4.6.1
* All optional integrations installed, including mgcv 1.9-4, sf, units, clock,
  tsibble, cubble, ape, rsample and ggplot2

## R CMD check

`R CMD build .`, then `R CMD check --as-cran naboer_0.1.0.tar.gz`:

**0 errors, 0 warnings, 1 NOTE.**

The only incoming NOTE is **New submission**, expected for a package not yet on
CRAN. There are no development-version or unpublished-URL diagnostics.
PDF/HTML manuals, examples, tests and all five vignettes (including rebuilding)
were checked. No `--no-manual`, `--no-build-vignettes`, forced Suggests suppression
or incoming-check suppression was used.

## Other validation

651 local assertions pass without warnings or skips. Coverage is 97.18% overall
and 96.54% core. The reference site and articles build. Million-row temporal and
space-time benchmarks pass direct references; the latter includes packed export
and operating-system peak-memory measurement. Weighted Gaussian tests check
block assessments and repeated scoring. Large-omission approximation results
and scientific limits are documented in `tools/REVIEW.md`.

## Remote checks

Pre-review implementation `92db241` passed release R checks on Linux, Windows
and macOS, Linux devel/oldrel/R4.4, and the minimal-dependency job. Reviewed
candidate runs will be recorded after pushing.

Standard remote checks use `--no-manual`, omitting the PDF manual; they build and
check vignettes. The full local candidate check includes the manual. The minimal
job deliberately omits optional integrations and tests core/minimal contexts
with `_R_CHECK_FORCE_SUGGESTS_=false`.

No release tag, GitHub release or CRAN submission has been created.
