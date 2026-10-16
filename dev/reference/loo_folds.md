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

  Tidy selection of independent domain columns.

- id:

  Optional tidy selection of a unique observation identifier.

## Value

A compact `naboer_folds` object.
