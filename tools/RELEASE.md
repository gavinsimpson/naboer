# Release checklist for 0.1.0

Candidate **0.1.0**, reviewed 2026-10-17. **Do not tag or submit to CRAN without
a new explicit instruction.** See `tools/REVIEW.md` for API/scientific findings
and `tools/VALIDATION.md` for evidence and validation scope.

## Candidate review and local validation

- [x] Review API, canonical order, assessment/omission semantics, prediction
      targets, domains and optional adapters; fix identified issues.
- [x] Document mgcv 1.9-4 compatibility limits: at least two exported folds;
      explicit sampled fold vectors for nonsingleton omission sets.
- [x] Review large-omission approximation against exact fixed-basis refits;
      distinguish fixture evidence from general accuracy/uncertainty guarantees.
- [x] Full optional-dependency tests pass: 651 assertions, no warnings or skips.
- [x] Coverage gates pass: 97.20% overall and 96.61% core.
- [x] README examples, five vignettes, reference manual and pkgdown build.
- [x] Million-row temporal and fixed-station space-time fixtures pass direct
      references; record hardware, runtime, output size, R heap and process peak
      memory separately.
- [x] Build and full `R CMD check --as-cran naboer_0.1.0.tar.gz` include manuals
      and vignettes: 0 errors/warnings, only expected New submission NOTE.
- [x] Review NEWS classifications; amend unreleased feature entries once and
      update DESCRIPTION/NEWS to 0.1.0.
- [x] Review public links and incoming checks; update `cran-comments.md` with
      actual local results and the expected NOTE.

## Remote gates

- [x] Baseline `92db241`: platform matrix, minimal job, OIDC Codecov uploads,
      pkgdown and actual GitHub Pages deployment all verified.
- [x] Reviewed candidate `617f710`: Linux/macOS/Windows release, Linux devel,
      oldrel and R4.4 all passed in
      [R-CMD-check 39008110486](https://github.com/gavinsimpson/naboer/actions/runs/39008110486).
- [x] Candidate minimal-dependency job passed in the same run.
- [x] [Coverage 39008110491](https://github.com/gavinsimpson/naboer/actions/runs/39008110491)
      passed; Codecov fully processed the report: 97.20% overall, 96.61% core.
- [x] [pkgdown 39008110499](https://github.com/gavinsimpson/naboer/actions/runs/39008110499)
      and [Pages 39008233996](https://github.com/gavinsimpson/naboer/actions/runs/39008233996)
      passed. Live version 0.1.0 and revised articles/API pages verified.

Subsequent evidence-only edits affect files excluded by `.Rbuildignore`, not the
validated package payload. Their push workflows are still monitored.

Standard push CI uses `--no-manual` and omits the PDF manual; vignettes are
built and checked. The minimal job intentionally excludes optional adapters
with `_R_CHECK_FORCE_SUGGESTS_=false`. The full release job is dispatch-only.
No such exclusions were used for the full local candidate check.

## Deliberately deferred

- [ ] Maintainer approval of the candidate and documented scientific limits.
- [ ] Create a tag/GitHub release — not authorised in this task.
- [ ] Submit to CRAN — not authorised in this task.

Window choice and model-specific approximation checks remain the analyst's
responsibility. Later source changes require proportionate revalidation.
