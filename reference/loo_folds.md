# Leave-one-observation-out folds

Leave-one-observation-out folds

## Usage

``` r
loo_folds(data, .by = NULL, id = NULL)
```

## Arguments

- data:

  A data frame, tibble, sf object, tsibble or cubble.

- .by:

  Tidy selection of independent domain columns. These constrain
  generated dependence queries and qualify group/block identities.
  Custom indices and existing omissions are retained as supplied; V-fold
  allocation can combine whole blocks from several domains in one fold.

- id:

  Optional tidy selection of a unique observation identifier.

## Value

A compact `naboer_folds` object.
