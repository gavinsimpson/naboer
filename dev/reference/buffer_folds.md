# Add dependence buffers around fold assessment observations

Existing omissions are retained. Only assessment observations define
buffer centres; buffers do not recursively expand from other omitted
observations.

## Usage

``` r
buffer_folds(
  x,
  spatial = NULL,
  temporal = NULL,
  phylogenetic = NULL,
  combine = NULL,
  key = NULL,
  index = NULL,
  within = NULL
)
```

## Arguments

- x:

  A `naboer_folds` object.

- spatial:

  A
  [`spatial_window()`](https://gavinsimpson.github.io/naboer/dev/reference/spatial_window.md)
  specification.

- temporal:

  A
  [`temporal_window()`](https://gavinsimpson.github.io/naboer/dev/reference/temporal_window.md)
  specification.

- phylogenetic:

  A
  [`phylogenetic_window()`](https://gavinsimpson.github.io/naboer/dev/reference/phylogenetic_window.md)
  specification.

- combine:

  Either `"intersection"` or `"union"`; required for multiple windows.

- key, index:

  Optional tidy selections overriding the stored series key and temporal
  index for buffer construction.

- within:

  For temporal searches, `"key"` or `"domain"`. Required when purely
  temporal data contain multiple series. Defaults to `"domain"` for
  spatiotemporal searches.

## Value

A `naboer_folds` object with expanded omission membership.
