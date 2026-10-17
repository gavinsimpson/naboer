# Summarise a cross-validation design

Summarise a cross-validation design

## Usage

``` r
# S3 method for class 'naboer_folds'
summary(object, ...)
```

## Arguments

- object:

  A `naboer_folds` object.

- ...:

  Unused; must be empty.

## Value

A list with design counts, fold-size ranges and an observation-level
assessment-frequency tibble. Summary never expands membership
list-columns. `buffer_size` counts omission-only observations, including
forecasting exclusions outside the assessment horizon or training
window.
