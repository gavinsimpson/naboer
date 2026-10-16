# Release checklist for 0.1.0

The version remains 0.0.1.9000 until these gates are independently verified.
Local checks alone do not certify GitHub's platform matrix.

- [x] Local full dependency tests pass with no skipped optional integrations.
- [ ] Checks pass on Linux, macOS, Windows, R devel, oldrel and R 4.4.
- [ ] Minimal installation gives clear optional-dependency errors.
- [x] Overall line coverage >=90%; core storage/export/validation >=95%.
- [x] README examples, vignettes and pkgdown build locally.
- [x] Million-row temporal and fixed-station spatiotemporal fixtures pass;
      records include hardware, runtime, output size and peak memory.
- [ ] Approximation study reviewed for large omission fractions.
- [ ] Full `R CMD build .` and `R CMD check --as-cran naboer_0.1.0.tar.gz`
      include manuals/vignettes and have no unexplained diagnostics.
- [ ] GitHub Pages configured for gh-pages; Codecov credentials configured if
      needed; remote workflow results and uploads verified.
- [ ] NEWS classifications reviewed and DESCRIPTION updated to 0.1.0.
- [ ] Links and CRAN incoming checks reviewed; cran-comments.md records actual
      platforms, results and expected notes.

Standard PR CI uses `--no-manual` while building: it omits the PDF reference
manual, so it cannot replace the full release check. The minimal dependency job
uses `_R_CHECK_FORCE_SUGGESTS_=false` and intentionally excludes optional adapters.
Always distinguish these exclusions from package failures in reports.

Local validation does not submit to CRAN, tag a release or publish remotely.
