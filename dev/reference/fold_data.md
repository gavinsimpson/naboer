# Retrieve the canonical observation data

Retrieve the canonical observation data

## Usage

``` r
fold_data(x)
```

## Arguments

- x:

  A `naboer_folds` object.

## Value

A tibble in the exact order used by exported indices. For cubbles,
`attr(result, "naboer_source")` maps observations to source rows.
