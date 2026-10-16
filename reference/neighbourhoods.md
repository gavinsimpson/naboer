# Construct focal-observation neighbourhoods

Each observation is assessed once, omitting its dependence
neighbourhood. Multiple windows require explicit union or intersection.
Series keys identify observations; `.by` prevents dependence across
independent domains.

## Usage

``` r
neighbourhoods(
  data,
  spatial = NULL,
  temporal = NULL,
  phylogenetic = NULL,
  combine = NULL,
  key = NULL,
  index = NULL,
  .by = NULL,
  within = NULL,
  id = NULL
)
```

## Arguments

- data:

  A data frame, tibble, sf object, tsibble or cubble.

- spatial:

  A
  [`spatial_window()`](https://gavinsimpson.github.io/naboer/reference/spatial_window.md)
  specification.

- temporal:

  A
  [`temporal_window()`](https://gavinsimpson.github.io/naboer/reference/temporal_window.md)
  specification.

- phylogenetic:

  A
  [`phylogenetic_window()`](https://gavinsimpson.github.io/naboer/reference/phylogenetic_window.md)
  specification.

- combine:

  Either `"intersection"` or `"union"`; required for multiple windows.

- key, index, .by, id:

  Tidy selections of series keys, temporal index, independent domains
  and observation identifiers.

- within:

  For temporal searches, `"key"` or `"domain"`. Required when purely
  temporal data contain multiple series. Defaults to `"domain"` for
  spatiotemporal searches.

## Value

A compact `naboer_folds` object.
