# Construct spatiotemporal block folds

Construct spatiotemporal block folds

## Usage

``` r
spatiotemporal_folds(
  data,
  spatial_size = NULL,
  temporal_size = NULL,
  mode = c("crossed", "space", "time"),
  shape = c("square", "hexagonal"),
  regions = NULL,
  assignment = NULL,
  spatial_origin = NULL,
  temporal_origin = NULL,
  grid_crs = NULL,
  within = "domain",
  v = NULL,
  balance = TRUE,
  seed = NULL,
  key = NULL,
  index = NULL,
  .by = NULL,
  id = NULL
)
```

## Arguments

- data:

  A data frame, tibble, sf object, tsibble, or cubble.

- spatial_size:

  Spatial grid width with physical distance units.

- temporal_size:

  A time span specifying temporal block width.

- mode:

  Cross space and time blocks, or hold out space/time blocks alone.

- shape:

  Square or hexagonal grid cells.

- regions:

  Optional sf/sfc polygon regions instead of a generated grid.

- assignment:

  Optional tidy selection of a column containing integer region row
  numbers. Required to resolve overlapping custom regions.

- spatial_origin, temporal_origin:

  Origins for the respective block systems.

- grid_crs:

  Optional projected CRS for generating grid cells.

- within:

  Temporal grouping, `"key"` or `"domain"`; the default is `"domain"` so
  temporal blocks align across sites.

- v:

  Optional number of folds to assign whole blocks to.

- balance:

  Balance observation counts when assigning whole blocks.

- seed:

  Optional local random seed.

- key, index, .by, id:

  Tidy selections of series keys, temporal index, independent domains,
  and observation identifiers.

## Value

A `naboer_folds` object.
