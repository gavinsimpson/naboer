This is a development build (0.0.1.9000), not a CRAN submission.

## Local test environment

* macOS Tahoe 26.6.1, arm64, R 4.6.1
* All optional adapters installed, including sf, cubble, tsibble, ape and rsample

## R CMD check

`R CMD check --as-cran`: 0 errors, 0 warnings, 1 NOTE.
Manual generation, examples, tests and vignette rebuilding were included.

The incoming NOTE identifies a new submission, the deliberately unreleased
version component 9000, and URLs for the not-yet-published pkgdown site and
GitHub workflow. Resolve the public URLs and select version 0.1.0 before
submission. These diagnostics are release prerequisites, not test failures.

## Other validation

509 test assertions pass without warnings or skips. Line coverage: 97.11%
overall and 96.34% across core storage, construction, export and validation.
The package website and all five articles build locally. Million-observation
temporal and fixed-station spatiotemporal benchmarks pass independent membership
checks, including end-to-end packed export for the latter.

Linux, Windows, R devel, oldrel and R 4.4 results must be recorded from GitHub
Actions before submission. The workflows are provided but have not been run
remotely from this working tree.
