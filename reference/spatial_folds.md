# Construct spatial cross-validation folds

Blocks are defined independently of buffers and fold counts. Shared
boundaries are assigned to the first region in deterministic grid/region
order. Features intersecting several blocks are also assigned to the
first intersected block; their original geometries are preserved.
Geographic grids use a local azimuthal equidistant projection centred on
the data bounding box unless `grid_crs` is supplied. Custom regions
retain their own geometry.

## Usage

``` r
spatial_folds(
  data,
  size = NULL,
  shape = c("square", "hexagonal"),
  regions = NULL,
  assignment = NULL,
  origin = NULL,
  grid_crs = NULL,
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

- size:

  A positive scalar with physical distance units, defining grid width.

- shape:

  Square or hexagonal grid cells.

- regions:

  Optional sf/sfc polygon regions instead of a generated grid.

- assignment:

  Optional tidy selection of a column containing integer region row
  numbers. Required to resolve overlapping custom regions.

- origin:

  Optional numeric x/y grid anchor in the grid CRS. The grid extends in
  both directions from the anchor: square cell boundaries pass through
  it, and a hexagonal cell is centred on it. The default is the
  coordinate-wise floor of the bounding-box minimum divided by cell
  width.

- grid_crs:

  Optional projected CRS for generating grid cells.

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
