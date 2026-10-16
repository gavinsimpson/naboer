# Construct custom cross-validation folds

Assessment observations are scored; omission observations are excluded
from fitting. Every assessment set must be a subset of its omission set.
Indices refer to
[`fold_data()`](https://gavinsimpson.github.io/naboer/reference/fold_data.md)
order and must be unique within a set.

## Usage

``` r
custom_folds(data, assessment, omission = assessment, .by = NULL, id = NULL)
```

## Arguments

- data:

  A data frame, tibble, sf object, tsibble or cubble.

- assessment, omission:

  Lists of observation indices, one vector per fold. Omission defaults
  to assessment.

- .by:

  Tidy selection of independent domain columns.

- id:

  Optional tidy selection of a unique observation identifier.

## Value

A compact `naboer_folds` object.

## Examples

``` r
x <- tibble::tibble(time = 1:6, y = c(1, 3, 2, 4, 6, 5))
folds <- custom_folds(x, list(2L, 5L), list(1:3, 4:6))
as_nei(folds)
#> $a
#> [1] 1 2 3 4 5 6
#> 
#> $ma
#> [1] 3 6
#> 
#> $d
#> [1] 2 5
#> 
#> $md
#> [1] 1 2
#> 
```
