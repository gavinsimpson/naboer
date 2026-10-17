# Release-candidate validation — 2026-10-17

Candidate **0.1.0**. No release tag or CRAN submission has been made. See
`tools/REVIEW.md` for the API/scientific review and documented limits.

## Local candidate checks

- R 4.6.1, macOS Tahoe 26.6.1, arm64; Apple M1 Pro, 32 GiB RAM.
- Optional adapters: mgcv 1.9-4, sf 1.1-3, clock 0.7.4, ape 5.8-1,
  rsample 1.3.2, tsibble 1.2.0 and cubble 1.0.0. Spatial libraries:
  GEOS 3.13.0, GDAL 3.8.5 and PROJ 9.5.1.
- `devtools::test()`: **651 assertions, 0 failures/warnings/skips**, with all
  optional integrations installed. Weighted Gaussian refits check block
  assessments, repeated scoring and large omissions independently.
- `covr::package_coverage()` (covr 3.6.5): **97.20% overall**, **96.61% core**
  across folds/indices/nei/utils. `tools/coverage-gate.R` passes.
- `R CMD build .`: all five vignettes built. Candidate tarball
  `naboer_0.1.0.tar.gz` has MD5 `76ba4b53f8e377e5293bc7e3e4d3d24f`.
- `R CMD check --as-cran naboer_0.1.0.tar.gz`: **0 errors, 0 warnings, 1 NOTE**.
  The only incoming NOTE says **New submission**. No development-version or
  unpublished-URL diagnostic remains. PDF/HTML manuals, examples, tests and
  vignette rebuilding were included. No `--no-manual`, `--no-build-vignettes`,
  force-Suggests suppression or incoming-check suppression was used.
- `pkgdown::build_site()`: reference pages and all five articles built. Upstream
  Pandoc deprecation messages for `--mathml`/`--mathjax` did not prevent building;
  they are not R CMD check warnings. A representative manual page was visually
  inspected.
- All README R examples were separately executed successfully against 0.1.0,
  including the fitted GAM and combined spatial/temporal workflow.
- Added R code/examples use quoted package names in `library()` calls.
- `git diff --check`: clean.

Full check/site commands used approved network access for incoming and URL
checks. Formatter caching was disabled because the sandbox could not write the
default user cache; this did not affect package checks.

## Million-row candidate benchmarks

Both R scripts ran with `/usr/bin/time -l`. Successful timing runs used approved
access because the sandbox blocked the utility's macOS `sysctl` query. The two
processes briefly overlapped; timings are observations, not controlled speed
comparisons or guarantees.

| Fixture | Observations | Omission memberships | Construction | Packed output |
|---|---:|---:|---:|---:|
| Irregular unequal temporal series, sparse | 1,000,000 | 3,796,000 | 0.398 s | 27,184,000 bytes estimated |
| Temporal whole-series exclusions | 1,000,000 | 1,078,979,374 | 0.338 s | 4,327,917,496 bytes estimated; guard rejects |
| 1,000 stations × 1,000 times | 1,000,000 | 14,972,012 | 13.131 s | 71,888,736 bytes measured |

Temporal process maximum RSS was **473,169,920 bytes**, peak footprint
**456,721,320 bytes**. R heap high-water was 200.2 MiB for sparse and 272.1 MiB
for dense construction. Sampled memberships matched direct scans.

Space-time export took **45.612 s**. Process maximum RSS was **1,429,880,832
bytes**, peak footprint **1,301,382,824 bytes**. Export R heap high-water was
**922.8 MiB**. Compact object size was 564,000,952 bytes, including canonical
metadata. Five boundary/interior folds and packed endpoints matched full-data
scans. No observation-level dense distance matrix was used.

RSS, macOS footprint, R heap high-water and final packed bytes are different
measurements. The 1 GiB guard limits estimated exported arrays, not total process
memory. Independent review probes also matched 24 projected/geographic joint
queries and 30 hex-grid subset fixtures against direct references.

## Statistical evidence

`benchmarks/approximation.R` compared eight radii for each of two families at
n=120. `tools/approximation-review.csv` records all results; `tools/REVIEW.md`
interprets them. All refits converged. Gaussian NCV/QNCV matched exact
fixed-basis refits to <=3.3e-13 relative criterion discrepancy. At the widest
Poisson window (maximum omission 80.8%), NCV differed by 15.36% and QNCV by
68.98%, with maximum link prediction error 0.805. These are stress-test results,
not general accuracy guarantees or smoothing/interval coverage validation.

Use `fit$gcv.ubre` for the fitted criterion: mgcv may update `fit$NCV` during
covariance calculations under its default neighbourhood policy. The study keeps
the full-data basis and penalty fixed.

## Remote validation

The pre-review implementation at `92db241` passed the complete push matrix and
minimal-dependency job. Codecov processed coverage, and the root site deployed:

- [R-CMD-check](https://github.com/gavinsimpson/naboer/actions/runs/38979324226)
- [Coverage](https://github.com/gavinsimpson/naboer/actions/runs/38979324131)
- [pkgdown](https://github.com/gavinsimpson/naboer/actions/runs/38979324237)
- [Pages](https://github.com/gavinsimpson/naboer/actions/runs/38979553126)

Reviewed candidate `617f710e95582a4159abf69d55fcda0edefc5ddd` passed:

- [R-CMD-check 39008110486](https://github.com/gavinsimpson/naboer/actions/runs/39008110486):
  release R on Linux/macOS/Windows, Linux devel/oldrel/R4.4 and minimal job.
- [Coverage 39008110491](https://github.com/gavinsimpson/naboer/actions/runs/39008110491):
  covr 97.20% overall, 96.61% core; OIDC upload accepted and Codecov state complete
  for the exact candidate SHA. Codecov's display truncates overall to 97.19%.
- [pkgdown 39008110499](https://github.com/gavinsimpson/naboer/actions/runs/39008110499)
  and [Pages 39008233996](https://github.com/gavinsimpson/naboer/actions/runs/39008233996):
  successful. Live root site returns HTTP 200 and shows 0.1.0; updated temporal,
  forecasting API, approximation and spatial guidance were checked.

Standard CI uses `--no-manual`, skipping the PDF manual, and builds/checks
vignettes. The minimal job uses `_R_CHECK_FORCE_SUGGESTS_=false` and intentionally
tests core/minimal contexts without optional adapters. The full release job is
dispatch-only and was not triggered. The successful full local candidate check
includes the manual and all optional integrations. Subsequent evidence-only
edits affect files excluded from the validated source package.
