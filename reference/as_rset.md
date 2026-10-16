# Convert between naboer folds and rsample resamples

Analysis sets are the complement of omission sets. Assessment membership
is preserved independently, including exclusion-only observations in
buffered designs. Bootstrap resamples with duplicated analysis rows are
not supported.

## Usage

``` r
as_rset(x)

from_rset(x, id = NULL)
```

## Arguments

- x:

  A `naboer_folds` object for `as_rset()`, or an rsample `rset` for
  `from_rset()`.

- id:

  Optional unique observation identifier column, using tidy selection.

## Value

`as_rset()` returns an rsample manual rset; `from_rset()` returns a
`naboer_folds` object.
