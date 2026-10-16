# Inspect folds as a tibble

Inspect folds as a tibble

## Usage

``` r
# S3 method for class 'naboer_folds'
as_tibble(
  x,
  ...,
  folds = NULL,
  membership = FALSE,
  .rows = NULL,
  .name_repair = "check_unique"
)
```

## Arguments

- x:

  A `naboer_folds` object.

- ...:

  Unused; must be empty.

- folds:

  Optional integer fold selection, defaulting to all folds.

- membership:

  Include assessment and omission list-columns. Defaults to FALSE so
  large designs are not expanded accidentally.

- .rows, .name_repair:

  Reserved tibble arguments.

## Value

A tibble with fold identifiers and analysis/assessment/omission counts.
