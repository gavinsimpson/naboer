# Plot selected cross-validation folds

Spatial geometries are used when available, otherwise the temporal index
or observation row is plotted. Only requested folds are expanded.

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
