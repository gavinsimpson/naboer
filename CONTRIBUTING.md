# Contributing

Use R for implementation, automation, benchmarks and supporting analysis. Quote
package names in every `library()` call. Follow tidyverse style, keep functions
pipeable with data first, and document public arguments with roxygen2.

Install development dependencies, then run:

```r
roxygen2::roxygenise()
devtools::test()
devtools::check()
pkgdown::build_site()
```

Tests compare membership with independently computed sets and statistical
criteria with exact Gaussian refits. Keep simulations and million-observation
benchmarks in `benchmarks/`, outside routine CRAN tests. When reporting check
results, name any options that caused skipped validation or expected notes.

An optional package belongs in Suggests. Check availability with a clear error
before using it; never attach it as a side effect. Do not silently reorder model
observations, repair geometry, fill time gaps, change global spatial options,
change the caller's RNG for an explicit seed, or reinterpret keys as independence.

Release notes classify new capabilities under New features, changes to existing
behaviour under User visible changes, and corrections under Bug fixes. Amend an
existing unreleased feature entry when completing it; record each change once.
