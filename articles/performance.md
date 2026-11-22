# Storage, diagnostics and performance

``` r

library("naboer")
library("tibble")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library("ggplot2")
}
```

Fold objects store the data, the observations to assess and the
observations to omit from training. Observations used for training are
called the analysis set. `naboer` works out that set when you ask for it
by taking the observations that remain after omission.

To save memory, a time window can be stored as start and end positions
in a shared time ordering. Folds with the same set of neighbours can
share that set, for example when several observations come from the same
site or taxon. For fixed stations, when a window requires observations
to be close in both distance and elapsed time, the design stores time
ranges for each nearby station.

## Inspecting folds

Create 1,000 observations one second apart. The time window leaves the
assessment observation and any observations within three seconds out of
training.

``` r

x <- tibble(time = seq_len(1000))
folds <- neighbourhoods(
  x,
  index = time,
  temporal = temporal_window(time_span(3, "seconds"))
)
```

Inspect folds at the start, middle and end of the series. Each fold
assesses one observation. The middle fold omits seven observations from
training, while the folds at the ends omit four each.

``` r

as_tibble(folds, folds = c(1, 500, 1000))
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                1          4        996        3            0.004
#> 2 500              1          7        993        6            0.007
#> 3 1000             1          4        996        3            0.004
```

Use [`summary()`](https://rdrr.io/r/base/summary.html) to check the
range of assessment and omission sizes, the fraction of observations
omitted, and how many observations are never assessed. This avoids
making a separate list of row numbers for every fold.

``` r

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

If you need the actual row numbers, use `membership = TRUE` with
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
for assessment and omission, or call
[`assessment_indices()`](https://gavinsimpson.github.io/naboer/reference/assessment_indices.md),
[`omission_indices()`](https://gavinsimpson.github.io/naboer/reference/assessment_indices.md)
or
[`analysis_indices()`](https://gavinsimpson.github.io/naboer/reference/assessment_indices.md)
directly. Supply `folds` to limit the result to the folds you need.
Asking for every fold can use much more memory than storing the design.
With `ggplot2` loaded, you can plot the observation roles in a single
fold with `autoplot(folds, folds = 1)`.

## Exporting to mgcv

`mgcv` takes a list called `nei`. It contains all the omission and
assessment row numbers in long integer vectors, along with two vectors
marking where each fold ends. The size depends on the total number of
entries across folds, so a compact fold object can still produce a large
export.
[`as_nei()`](https://gavinsimpson.github.io/naboer/reference/as_nei.md)
checks the integer index limit and estimates the size of these vectors
before creating them.

Export the example with a limit of 100 MiB, then count the entries in
each vector. The default limit is 1 GiB and can be changed with the
`naboer.max_bytes` option or for one call with `max_bytes`.

``` r

nei <- as_nei(folds, max_bytes = 100 * 1024^2)
lengths(nei)
#>    a   ma    d   md 
#> 6988 1000 1000 1000
```

The estimate covers the integer vectors. Leave room for the rest of the
R session and for temporary memory used during export. Sampling
neighbourhoods for the `bam()` criterion still requires the complete
`nei` structure. A forecasting fold can omit every observation after its
forecast origin. The same observations can be omitted at many origins,
so the export can grow with both the number of observations and the
number of origins.

Export and model fitting have separate costs. With
`bam(discrete = TRUE)` and neighbourhood cross-validation, fitting
requires a matrix inversion for each neighbourhood. The matrix gets
larger as more observations are omitted, so large blocks can be slow to
fit even when the exported vectors fit in memory. ([mgcv documentation
for
`bam()`](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/bam.html))

## Exporting to rsample

[`as_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
exports the training and assessment row numbers for each split. In a
leave-one-out design with one million rows, every split has 999,999
training rows. The training indices alone would require about 4 TB of
memory.
[`as_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
checks its own `max_bytes` limit, which also defaults to 1 GiB. The
estimate covers the index vectors. The `rsample` split objects and data
need additional memory.

## Reproducible benchmarks

The benchmark scripts cover million-row time series, fixed monitoring
stations, moving locations, polygon data and comparisons between
approximate and exact refits. These scripts run separately from routine
CRAN tests. ([Benchmark
scripts](https://github.com/gavinsimpson/naboer/tree/main/benchmarks))

When comparing runs, record the data size, how many observations each
fold includes, the sizes of the fold object and export, and the elapsed
time. Record the hardware and the versions of R and its packages.
Specify whether the peak is for memory managed by R or for the whole
process. Where possible, check neighbourhoods for a sample of
observations against a direct calculation.

The benchmarks also test whether memory limits catch designs with large
neighbourhoods before allocation. Two designs can have the same number
of rows and very different export sizes, depending on how many
neighbours each observation has.
