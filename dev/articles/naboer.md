# Assessment, omission and prediction targets

A fold assesses observations whose predictions contribute to a criterion
and omits observations unavailable during fitting. Assessment is a
subset of omission; the difference is an exclusion buffer. Analysis is
the complement of omission.

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

## Model alignment and mgcv

Prepare complete modelling data first. Constructors preserve canonical
order.
[`fold_data()`](https://gavinsimpson.github.io/naboer/dev/reference/fold_data.md)
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
[`as_nei()`](https://gavinsimpson.github.io/naboer/dev/reference/as_nei.md).
Discarding unsampled definitions changes the information available for
covariance.

A valid design does not guarantee approximation accuracy. Compare with
exact refits for representative models, particularly for large omission
fractions. `benchmarks/approximation.R` provides a reproducible
fixed-basis study.

## References

Roberts et al. (2017), *Cross-validation strategies for data with
temporal, spatial, hierarchical, or phylogenetic structure*, Ecography
40:913–929.
[doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881).

Wood, *On Neighbourhood Cross Validation*,
[arXiv:2404.16490](https://arxiv.org/abs/2404.16490).
