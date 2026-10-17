# Random V-fold cross-validation

Random V-fold cross-validation

## Usage

``` r
vfold_folds(
  data,
  v = 10L,
  group = NULL,
  balance = TRUE,
  seed = NULL,
  .by = NULL,
  id = NULL
)
```

## Arguments

- data:

  A data frame, tibble, sf object, tsibble or cubble.

- v:

  Number of folds, at least two.

- group:

  Optional tidy selection of columns kept together during allocation.

- balance:

  Balance observation counts rather than group counts.

- seed:

  Optional local random seed; does not change the caller's RNG state.

- .by:

  Tidy selection of independent domain columns. These constrain
  generated dependence queries and qualify group/block identities.
  Custom indices and existing omissions are retained as supplied; V-fold
  allocation can combine whole blocks from several domains in one fold.

- id:

  Optional tidy selection of a unique observation identifier.

## Value

A compact `naboer_folds` object.

## Details

`.by` qualifies group identities; it does not create a separate V-fold
schedule within each domain. Every observation is assessed once across
the requested `v` folds.
