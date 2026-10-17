# Construct environmental cross-validation folds

Select numeric predictors for k-means clustering, or supply existing
strata. Clustering is performed separately in each independent domain.
Standardising predictors gives each nonconstant variable equal scale.
Constant predictors are reported and excluded within each domain. For
clustering, `v` is the number of clusters per domain, and each cluster
forms one holdout. Clusters can have unequal sizes. `balance` only
affects allocation of preexisting `strata`, not k-means clustering.

## Usage

``` r
environmental_folds(
  data,
  predictors = NULL,
  v = NULL,
  strata = NULL,
  standardise = TRUE,
  balance = TRUE,
  seed = NULL,
  .by = NULL,
  id = NULL
)
```

## Arguments

- data:

  An observation data frame.

- predictors:

  Numeric predictor columns, using tidy selection.

- v:

  Number of k-means clusters, or folds of whole existing strata.

- strata:

  Existing environmental strata columns, using tidy selection.

- standardise:

  Standardise numeric predictors before clustering.

- balance:

  Balance observation counts when allocating groups to folds.

- seed:

  Optional local random seed.

- .by:

  Columns defining independent domains.

- id:

  Optional unique observation identifier column.
