# Construct contiguous temporal blocks

Construct contiguous temporal blocks

## Usage

``` r
temporal_folds(
  data,
  size,
  index = NULL,
  key = NULL,
  .by = NULL,
  id = NULL,
  origin = NULL,
  within = NULL,
  v = NULL,
  balance = TRUE,
  seed = NULL
)
```

## Arguments

- data:

  A data frame, tsibble, or cubble.

- size:

  A positive
  [`time_span()`](https://gavinsimpson.github.io/naboer/dev/reference/time_span.md)
  specifying block width.

- index, key, .by, id:

  Columns identifying time, series, independent domains, and unique
  observations. Use unquoted column names or tidy selections.

- origin:

  A scalar index value defining block boundaries. By default use the
  earliest observation within each independent domain or series.

- within:

  Choose `"key"` or `"domain"` for multiple series. Count spans require
  `"key"` and count rows in stable time order.

- v:

  Optional number of folds to which whole blocks are allocated.

- balance:

  Balance observation counts when assigning blocks to folds.

- seed:

  Optional local random seed.

## Value

A `naboer_folds` object.
