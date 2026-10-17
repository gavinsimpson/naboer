# Plot selected cross-validation folds

Spatial geometries are used when available, otherwise the temporal index
or observation row is plotted. Only requested folds are expanded. The
`omitted` role denotes omission-only observations; these can be
dependence buffers, unavailable future observations, or other
exclusions.

## Usage

``` r
# S3 method for class 'naboer_folds'
autoplot(object, folds = 1L, ...)
```

## Arguments

- object:

  A `naboer_folds` object.

- folds:

  Integer fold selection; defaults to the first fold.

- ...:

  Unused; must be empty.

## Value

A ggplot object.
