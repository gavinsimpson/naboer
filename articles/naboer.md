# Assessment, omission and prediction targets

A fold assesses observations whose predictions contribute to a criterion
and omits observations unavailable during fitting. Assessment is a
subset of omission; the difference contains exclusion-only observations.
These may be dependence buffers, unavailable future observations, or
observations outside a sliding training window. Analysis is the
complement of omission. Diagnostics retain the names `n_buffer` and
`buffer_size` for these exclusion-only counts; plots label them
`omitted`.

``` r

x <- tibble::tibble(time = 1:8, y = c(1, 2, 2, 4, 3, 5, 7, 6))
folds <- custom_folds(x, assessment = list(3L, 6L), omission = list(2:4, 5:7))
tibble::as_tibble(folds, membership = TRUE)
#> # A tibble: 2 × 8
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction assessment
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl> <list>    
#> 1 1                1          3          5        2            0.375 <int [1]> 
#> 2 2                1          3          5        2            0.375 <int [1]> 
#> # ℹ 1 more variable: omission <list>
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

`nei$a` concatenates omission indices; `ma` contains cumulative
endpoints. `d` and `md` encode assessment. The first endpoint starts at
position one; there is no leading zero. Repeated assessment across folds
is retained, including its effect on the criterion. No automatic fold
weighting is applied.

| Prediction target | Design |
|----|----|
| Baseline prediction among sampled units | LOO or random V-fold, where appropriate for the dependence and target |
| Beyond local residual dependence | Focal neighbourhoods, optionally buffered blocks |
| New locations in the sampled region | Spatial block holdouts |
| New periods | Temporal blocks |
| Future observations | Forecast origins, horizon and gap |
| New groups | Whole-group or hierarchical holdouts |
| New environments | Predictor-space clusters or strata |
| New taxa | Phylogenetic clades, distance blocks or neighbourhoods |

These designs pose different questions. Blocked validation can require
extrapolation in predictors; a large error does not automatically imply
a poor interpolation model. Choose dependence radii using scientific
knowledge or separate diagnostics; naboer does not infer them
automatically.

Keys identify repeated series or observational units; they do not
declare independence. Use `.by` only for domains across which generated
dependence queries should never pass. Custom membership and existing
exclusions are retained as supplied. Whole-block V-fold allocation can
combine several domains in one fold; `.by` does not request a separate
V-fold schedule per domain. Observations in other domains remain
available for fitting, including future observations in forecasting
designs. If information would not be available at prediction time, do
not declare it to be an independent domain merely to obtain a particular
fold layout.

Blocking in space or time can also separate the training and assessment
ranges of environmental predictors. Match the design to the intended
deployment: interpolation among sampled units, transfer to new units,
and extrapolation to new environments are different estimands. A valid
fold design does not ensure that every training model is identifiable or
can predict previously unseen factor levels. Check those properties for
the model and prediction target.

## Model alignment and mgcv

Prepare complete modelling data first. Constructors preserve canonical
order.
[`fold_data()`](https://gavinsimpson.github.io/naboer/reference/fold_data.md)
returns that order; `as_nei(data = ...)` can check it. With `id`,
alignment checks compare identifiers; otherwise all canonical columns
must match. Rebuild folds after changing observations.

``` r

all_folds <- loo_folds(x)
fit <- mgcv::gam(y ~ s(time, k = 4),
  data = fold_data(all_folds),
  method = "NCV", nei = as_nei(all_folds), na.action = stats::na.fail
)
fit$sp
#>  s(time) 
#> 28199.07
```

`gam()` supports NCV and QNCV; `bam()` supports NCV with
`discrete = TRUE`. Leave `jackknife` unspecified to preserve mgcv’s
covariance default. For sampled bam criteria, keep every neighbourhood
and supply `sample` through
[`as_nei()`](https://gavinsimpson.github.io/naboer/reference/as_nei.md).
Discarding unsampled definitions changes the information available for
covariance. For buffered or blocked designs, use at least two explicit
sampled fold numbers, for example
`as_nei(all_folds, sample = c(1L, 3L, 5L))`. In the tested mgcv 1.9-4
implementation, scalar counts sample packed endpoints rather than fold
numbers, so naboer rejects that form unless all omission sets are
singletons. Export also requires at least two folds: a single native
fold is useful for inspection or explicit rsample refitting, but the
tested `gam()` path rejects it. The initial integration baseline is mgcv
1.9-4. For `bam()`, NCV is applied to the working penalised linear
model; its score should not be interpreted as an exact non-Gaussian
refit criterion. Large blocked omissions can be particularly expensive
in `bam()`; its NCV path is intended for small neighbourhoods.

## Approximation accuracy and large omissions

A valid design does not guarantee approximation accuracy. Except for
quadratic loss, NCV uses a local approximation to the omitted-data fit.
QNCV also approximates the criterion; it provides a finite criterion
where supported, not a guarantee of agreement with exact refits or of
improved accuracy over NCV. The asymptotic argument in the NCV paper
assumes neighbourhood size does not increase with sample size. Large
blocks and forecast omissions need separate checks rather than an appeal
to that small-neighbourhood result.

`benchmarks/approximation.R` compares both criteria with exact refits
using the same full-data basis and fixed smoothing parameter. Rebuilding
a smooth basis on each training subset can change its centring or
penalty scaling, confounding the comparison. With 120 observations, seed
20261017 and smoothing parameter 0.2, Gaussian relative criterion
discrepancies were below `3.3e-13`. The corresponding Poisson stress
test gave:

| Temporal radius | Maximum omitted fraction | NCV relative criterion error | QNCV relative criterion error |
|---:|---:|---:|---:|
| 7 | 12.5% | 0.23% | 0.47% |
| 15 | 25.8% | 6.04% | 15.10% |
| 30 | 50.8% | 0.37% | 25.37% |
| 48 | 80.8% | 15.36% | 68.98% |

These are one-fixture results, not accuracy bounds or recommended radii.
Errors need not grow monotonically with window width, and total
criterion agreement can hide differences between individual predictions.
The script also records the largest link-scale prediction discrepancy
and checks refit convergence. It does not test the accuracy of optimised
smoothing parameters or interval coverage. For a scientific application,
compare approximate and exact refits on representative folds and
candidate smoothing parameters. Use explicit refitting through
[`as_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
when the approximation is unsuitable.

See the [mgcv NCV
documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/NCV.html)
and [bam
documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/bam.html)
for the model-specific interpretation and computational limits.

## References

Roberts et al. (2017), *Cross-validation strategies for data with
temporal, spatial, hierarchical, or phylogenetic structure*, Ecography
40:913–929.
[doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881).

Wood, *On Neighbourhood Cross Validation*,
[arXiv:2404.16490](https://arxiv.org/abs/2404.16490).
