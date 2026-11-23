# Plot selected cross-validation folds

By default, plots use observation row numbers. Supplying the original
data lets the plot use spatial geometries or the temporal index when
available. Only requested folds are expanded. The `omitted` role denotes
omission-only observations; these can be dependence buffers, unavailable
future observations, or other exclusions.

## Usage

``` r
# S3 method for class 'naboer_folds'
autoplot(object, folds = 1L, ..., data = NULL)
```

## Arguments

- object:

  A `naboer_folds` object.

- folds:

  Integer fold selection; defaults to the first fold.

- ...:

  Unused; must be empty.

- data:

  Optional observation data used to construct the folds, in their
  original row order.

## Value

A ggplot object.
