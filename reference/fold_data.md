# Prepare observations in the order used by a fold object

Prepare observations in the order used by a fold object

## Usage

``` r
fold_data(x, data)
```

## Arguments

- x:

  A `naboer_folds` object.

- data:

  The observation data used to construct the folds, in their original
  row order.

## Value

A tibble in the exact order used by exported indices. For cubbles,
`attr(result, "naboer_source")` maps observations to source rows.

## Details

Supply the original data. The observation count is checked after
preparation, and the row order is assumed. For cubbles, preparation
unfolds the observations and attaches their spatial data.
