# Convert between naboer folds and rsample resamples

Analysis sets are the complement of omission sets. Assessment membership
is preserved independently, including exclusion-only observations in
buffered designs. Bootstrap resamples with duplicated analysis rows are
not supported. Resamples must index individual observations. Convert
nested cubbles to their temporal face before creating resamples;
spatial-row indices cannot safely be reused after expanding nested
temporal observations. Exported manual rsets retain column metadata for
series keys, index, independent domains and observation identifiers.
Import restores it unless explicit tidy selections override it. External
rsets can supply these selections directly.

## Usage

``` r
as_rset(x, max_bytes = getOption("naboer.max_bytes", 1024^3))

from_rset(x, id = NULL, key = NULL, index = NULL, .by = NULL)
```

## Arguments

- x:

  A `naboer_folds` object for `as_rset()`, or an rsample `rset` for
  `from_rset()`.

- max_bytes:

  Maximum estimated storage for analysis and assessment indices exported
  to rsample, in bytes. Defaults to option `naboer.max_bytes`, or 1 GiB.
  The estimate excludes R object overhead. Use `Inf` to disable the
  memory guard after checking the design size.

- id:

  Optional unique observation identifier column, using tidy selection.

- key, index, .by:

  Tidy selections of series keys, temporal index and independent
  domains. Omitted selections restore metadata saved by `as_rset()`, or
  use native tsibble/cubble metadata where available.

## Value

`as_rset()` returns an rsample manual rset; `from_rset()` returns a
`naboer_folds` object.
