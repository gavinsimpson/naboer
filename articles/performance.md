# Storage, diagnostics and performance

The fold object keeps canonical data and compressed indices. Contiguous
temporal memberships use shared ordering and range endpoints; repeated
spatial/taxon sets share membership definitions. Fixed-station
space-time queries search temporal intervals for neighbouring station
geometries. Analysis sets are complements and are not stored as
duplicated vectors.

``` r

x <- tibble::tibble(time = seq_len(1000))
folds <- neighbourhoods(x,
  index = time,
  temporal = temporal_window(time_span(3, "seconds"))
)
print(folds)
#> <naboer_folds: neighbourhood>
#> 1,000 observations; 1,000 folds
tibble::as_tibble(folds, folds = c(1, 500, 1000))
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                1          4        996        3            0.004
#> 2 500              1          7        993        6            0.007
#> 3 1000             1          4        996        3            0.004
diagnostics <- summary(folds)
diagnostics[c("assessment_size", "omission_size", "omitted_fraction", "unassessed")]
#> $assessment_size
#> [1] 1 1
#> 
#> $omission_size
#> [1] 4 7
#> 
#> $omitted_fraction
#> [1] 0.004 0.007
#> 
#> $unassessed
#> [1] 0
```

`as_tibble()` omits membership list-columns unless `membership = TRUE`.
[`assessment_indices()`](https://gavinsimpson.github.io/naboer/reference/assessment_indices.md),
[`omission_indices()`](https://gavinsimpson.github.io/naboer/reference/assessment_indices.md)
and
[`analysis_indices()`](https://gavinsimpson.github.io/naboer/reference/assessment_indices.md)
expand only requested folds; asking for all folds is explicit but can
consume large amounts of memory.
[`summary()`](https://rdrr.io/r/base/summary.html) reports assessment
coverage without retaining expanded memberships. Plot selected folds
with `ggplot2::autoplot(folds, folds = 1)`.

## Materialising nei

`mgcv` consumes packed observation indices. This output is proportional
to the total membership count even when the fold object is compact.
Export checks the integer-index limit and estimates output bytes before
allocation. The default 1 GiB guard can be adjusted explicitly:

``` r

nei <- as_nei(folds, max_bytes = 100 * 1024^2)
sum(vapply(nei, length, integer(1)))
#> [1] 9988
```

The estimate covers packed integer arrays, not the total R process or
temporary allocation peaks. Use a limit below available memory. Sampling
bam neighbourhood criteria does not remove the requirement to define all
neighbourhoods. Forecast designs retain compressed future exclusions
until export, which can itself be quadratic in the number of
observations/origins.

## Reproducible benchmarks

The source repository’s `benchmarks/` directory contains R scripts for
million-row multiple-series data, fixed monitoring stations, moving
locations, polygons and approximation studies. They are excluded from
routine CRAN tests. Records should include input dimensions,
memberships, compressed/output sizes, runtime, R heap high-water or
process peak memory (clearly distinguished), hardware, R/dependency
versions and brute-force agreement on sampled focal observations.

Dense windows are deliberately tested with preallocation guards. A
sparse million-row design and a dense million-row design have very
different output costs; construction runtime does not predict packed
size.
