# Assessment, omission and prediction targets

``` r

library("naboer")
library("tibble")
```

naboer builds the `nei` argument that mgcv uses to select smoothness by
neighbourhood cross-validation (NCV) or its quadratic variant (QNCV).
Start with a worked example that makes a neighbourhood for each
observation, fits a smooth and compares it with leave-one-out selection.
[Soil moisture over
time](https://gavinsimpson.github.io/naboer/articles/temporal.html#choose-smoothness-for-a-time-series).
[Daily maximum temperature across weather
stations](https://gavinsimpson.github.io/naboer/articles/spatial.html#choose-smoothness-for-spatial-data).

Suppose you want to predict a measurement at a new location. If
measurements from nearby locations are similar, fitting the model with
those neighbours may make the prediction task easier than the one you
care about. A cross-validation fold can leave out the target location
and its neighbours, then score the prediction only at the target.

naboer distinguishes three sets of rows in each fold.

- **Assessment** rows are the ones whose predictions contribute to the
  score.
- **Omission** rows are left out when fitting. They include every
  assessment row and can also include neighbours, observations from the
  future, or observations outside a fixed-length training window.
- **Analysis** rows are the remaining rows used to fit the model.

## Inspect a fold

Create eight observations in time order. We will score rows 3 and 6 in
separate folds, leaving out one extra row on either side of each target.

``` r

x <- tibble(time = 1:8, y = c(1, 2, 2, 4, 3, 5, 7, 6))
```

Use
[`custom_folds()`](https://gavinsimpson.github.io/naboer/reference/custom_folds.md)
to specify the two sets of row numbers for each fold.

``` r

folds <- custom_folds(x, assessment = list(3L, 6L), omission = list(2:4, 5:7))
```

The fold table shows how many rows are in each set. The `n_buffer`
column counts the extra omissions. They appear as `buffer_size` in
[`summary()`](https://rdrr.io/r/base/summary.html) and as `omitted` in
plots.

``` r

as_tibble(folds)
#> # A tibble: 2 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                1          3          5        2            0.375
#> 2 2                1          3          5        2            0.375
```

Inspect the row numbers for the first fold. Row 3 is scored, rows 2 and
4 are also left out, and rows 1 and 5 to 8 are used for fitting.

``` r

assessment_indices(folds, folds = 1)
#> [[1]]
#> [1] 3
omission_indices(folds, folds = 1)
#> [[1]]
#> [1] 2 3 4
analysis_indices(folds, folds = 1)
#> [[1]]
#> [1] 1 5 6 7 8
```

## Choose a design for the prediction task

Start with where, when or for whom the predictions will be made. The
table shows common starting points.

| Prediction task | Design to consider |
|----|----|
| New observations from already sampled units | Leave-one-out or random V-fold, if neighbours in training reflect how predictions will be made |
| Predictions with no nearby training observations | Leave out a neighbourhood around each target, or add a buffer around each held-out block |
| New locations in the sampled region | Hold out spatial blocks |
| New periods | Hold out temporal blocks |
| Future observations given data available at the time | Use forecast origins, a horizon and, if needed, a gap |
| New groups | Hold out entire groups or combinations of grouping levels |
| New environments | Hold out groups defined by environmental predictors |
| New taxa | Hold out clades, distance-based blocks or neighbourhoods on a phylogenetic tree |

Spatial and temporal blocks can also separate the environmental
conditions of the training and assessment data. Larger errors may then
reflect the difficulty of predicting under unfamiliar conditions. Choose
a design that represents the range of conditions you expect to predict
in. For neighbourhood designs, choose the radius using knowledge of the
system or a separate analysis of dependence. Check that the model can be
fitted and can make the required predictions in each fold, especially if
assessment rows contain factor levels absent from training.

Use a key to identify repeated units or series. Use `.by` when the data
contain independent domains and automatically generated neighbourhoods
should stay within each domain. Custom fold membership and existing
exclusions are kept as supplied. A V-fold design can allocate whole
blocks from several domains to the same fold. When forecasting,
observations from other domains, including later observations, can
remain available for fitting. Define domains according to which
information would actually be available when making the prediction.

## Fit a model with mgcv

Prepare the data you will use in the model before constructing folds.
Row numbers in the design refer to the data returned by
[`fold_data()`](https://gavinsimpson.github.io/naboer/reference/fold_data.md).
Use that data for fitting, or supply your model data to
`as_nei(data = ...)` to check alignment. If you set `id` when
constructing the folds, the check compares the identifiers in order.
Otherwise, it compares all columns, including their order, types and
values. Rebuild the folds if you change the observations.

Pass the exported design as `nei` to mgcv’s `gam()` with
`method = "NCV"` or `method = "QNCV"`. Without `nei`, mgcv uses
leave-one-out cross-validation. With `nei`, mgcv scores predictions for
the supplied assessment rows under the supplied omissions. The final
model uses all observations, with smoothness chosen by that criterion.
QNCV uses a quadratic approximation to the score. For Gaussian models
with the identity link, as in the worked examples, the loss is already
quadratic. [mgcv NCV
documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/NCV.html).

## Read the exported indices

To see how larger omissions are represented, export the original
two-fold design.

``` r

as_nei(folds)
#> $a
#> [1] 2 3 4 5 6 7
#> 
#> $ma
#> [1] 3 6
#> 
#> $d
#> [1] 3 6
#> 
#> $md
#> [1] 1 2
```

The vector `a` lists the omitted row numbers for each fold in turn, and
`ma` marks where each fold ends in that vector. Here `ma` is `c(3, 6)`,
so the first fold uses positions 1 to 3 and the second uses positions 4
to 6. The pair `d` and `md` works the same way for assessment. If a row
is assessed in several folds, each assessment contributes to the
criterion without automatic reweighting.

### Sample folds in bam

`gam()` supports neighbourhood cross-validation (NCV) and its quadratic
variant (QNCV). `bam()` supports NCV with `discrete = TRUE`. To compute
the `bam()` criterion on a sample of folds, keep the full design and
pass `sample` to
[`as_nei()`](https://gavinsimpson.github.io/naboer/reference/as_nei.md).
For example, create one fold per row that omits its immediate neighbours
in time, then choose three folds explicitly.

``` r

buffered <- neighbourhoods(x,
  index = time, temporal = temporal_window(time_span(1))
)
sampled_nei <- as_nei(buffered, sample = c(1L, 3L, 5L))
sampled_nei$sample
#> [1] 1 3 5
```

When any fold omits more than one row, mgcv 1.9-4 requires explicit fold
numbers for sampling. From mgcv 1.9-5, you can also supply a count, such
as `sample = 3L`. To use every fold, leave `sample` unset. Export for
mgcv requires at least two folds. You can still inspect a single fold in
naboer or use it for a separate refit.

### Covariance and standard errors

Wood’s dependence-adjusted covariance assumes one assessment fold per
observation, reciprocal neighbourhoods, and negligible residual
correlation outside each neighbourhood. Reciprocity means that whenever
row `j` is in row `i`’s neighbourhood, row `i` is also in row `j`’s.
Asymmetric windows can still define a cross-validation score.

For mgcv through the tested 1.9-6 development version,
[`as_nei()`](https://gavinsimpson.github.io/naboer/reference/as_nei.md)
puts folds in data row order when each observation has exactly one fold
assessing it alone. It also translates any explicit sample numbers to
that order. On these versions, requesting the raw jackknife covariance
with `jackknife = TRUE` requires that same assessment layout.

In mgcv 1.9-5 and the tested 1.9-6 development version, the default
dependence-adjusted covariance in `gam()` can count a model offset
twice. Use `as_nei(folds, jackknife = FALSE)` if standard Bayesian
covariance is appropriate for your model. For discrete `bam()` on these
versions, check that prediction standard errors are finite and that the
covariance matrix has no substantially negative eigenvalues. A sample of
one fold can also trigger a warning about taking a maximum over an empty
set.

`bam()` applies NCV to a penalised linear approximation of the model. An
exact non-Gaussian refit uses the full loss function. The implementation
in `bam()` is designed for small neighbourhoods and can be expensive
when folds omit large blocks. [mgcv NCV
documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/NCV.html),
[mgcv bam
documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/bam.html).

## Check the approximation for large omissions

For non-quadratic loss functions, NCV approximates what happens to the
fitted model when observations are left out. QNCV also uses a quadratic
approximation to the score, which keeps the criterion finite for
supported models. The two methods can differ in how closely they match
an exact refit. Wood’s accuracy result as sample size increases assumes
the neighbourhood size stays fixed. Check accuracy separately for
designs that omit large blocks or much of the future record at each
forecast origin.

The repository’s benchmark compares both criteria with exact refits. It
holds the smoothing parameter and the basis for the smooth fixed to
isolate the effect of the approximation. [Approximation
script](https://github.com/gavinsimpson/naboer/blob/main/benchmarks/approximation.R).

For 120 observations, seed 20261017 and a smoothing parameter of 0.2,
relative differences in the Gaussian scores were below `1e-12`. In the
Poisson example, the differences varied with the radius of the temporal
window.

| Temporal radius | Maximum omitted fraction | NCV relative criterion error | QNCV relative criterion error |
|---:|---:|---:|---:|
| 7 | 12.5% | 0.23% | 0.47% |
| 15 | 25.8% | 6.04% | 15.10% |
| 30 | 50.8% | 0.37% | 25.37% |
| 48 | 80.8% | 15.36% | 68.98% |

These results describe one simulated data set. The error can rise or
fall as the window gets wider, and an accurate total score can hide
discrepancies in individual predictions. The script also checks whether
refits converge and records the largest NCV prediction difference on the
scale of the model’s linear predictor. For your own model, compare both
scores and predictions across representative folds and smoothing
parameters. If the approximation is unsuitable, use
[`as_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
to export the folds for explicit refitting.

## References

Roberts et al. (2017), *Cross-validation strategies for data with
temporal, spatial, hierarchical, or phylogenetic structure*, Ecography
40, 913–929.
[doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881).

Wood (2026), *On Neighbourhood Cross Validation*, Statistical Science
41, 827–846.
[doi:10.1214/25-STS1015](https://doi.org/10.1214/25-STS1015).
