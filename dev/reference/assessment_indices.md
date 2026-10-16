# Retrieve observation indices for folds

Retrieve observation indices for folds

## Usage

``` r
assessment_indices(x, folds = NULL)

omission_indices(x, folds = NULL)

analysis_indices(x, folds = NULL)
```

## Arguments

- x:

  A `naboer_folds` object.

- folds:

  Optional integer vector selecting folds. Defaults to all folds.

## Value

A list of integer vectors. Expanding all large folds can use substantial
memory.
