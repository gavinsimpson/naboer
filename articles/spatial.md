# Spatial and spatiotemporal data

Use `sf` geometry columns and physical distance units. Projected queries
use indexed candidates with exact geometry-to-geometry distance checks.
Geographic queries respect the current sf spherical/ellipsoidal engine.
Constructors never change `sf_use_s2()`, repair invalid geometry or
assume a missing CRS.

Radius thresholds and block assignments use the stored coordinates. An
`sf` precision attribute is retained in
[`fold_data()`](https://gavinsimpson.github.io/naboer/reference/fold_data.md)
but does not implicitly round coordinates during construction. Round the
coordinates explicitly beforehand if the intended design uses rounded
locations. This also keeps indexed candidates consistent with exact
geometry distances.

``` r

sites <- sf::st_as_sf(
  tibble::tibble(site = letters[1:6], east = c(0, 100, 500, 1000, 2000, 3000), north = 0),
  coords = c("east", "north"), crs = 32632
)
f <- neighbourhoods(sites, spatial = spatial_window(units::set_units(200, "m")))
tibble::as_tibble(f, folds = 1:3, membership = TRUE)
#> # A tibble: 3 × 8
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction assessment
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl> <list>    
#> 1 1                1          2          4        1            0.333 <int [1]> 
#> 2 2                1          2          4        1            0.333 <int [1]> 
#> 3 3                1          1          5        0            0.167 <int [1]> 
#> # ℹ 1 more variable: omission <list>
```

Square and hexagonal block widths are distinct from buffers. Custom
polygon regions are also supported. A shared boundary belongs to the
first region in deterministic region order. Overlapping region interiors
require an explicit integer `assignment` column. Non-point observations
retain their original geometry; block assignment uses the first
intersecting block rather than silently substituting a centroid.

Use `origin = c(x, y)` to keep a grid anchored across datasets or
subsets, including negative coordinates. Square boundaries pass through
the anchor; a hexagonal cell is centred on it. Hexagonal `size` is the
distance between opposite edges, not the side length. Shared-boundary
ties use cells ordered south to north, then west to east. With the same
explicit anchor, changing the data extent does not move the grid. The
default anchor is calculated from each dataset’s extent.

Grid creation checks its estimated allocation against
`options(naboer.max_bytes)` (1 GiB by default) before allocating polygon
cells. Increasing the block size or limiting the spatial extent can
reduce the allocation.

``` r

blocks <- spatial_folds(sites, size = units::set_units(500, "m"))
buffered <- buffer_folds(blocks, spatial = spatial_window(units::set_units(100, "m")))
tibble::as_tibble(buffered)
#> # A tibble: 5 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                1          2          4        1            0.333
#> 2 2                2          3          3        1            0.5  
#> 3 3                1          1          5        0            0.167
#> 4 4                1          1          5        0            0.167
#> 5 5                1          1          5        0            0.167
```

Generated geographic grids use a local azimuthal equidistant projection.
For large extents, select an appropriate explicit projected `grid_crs`.
Physical distance radii remain separate from that grid representation.
Geographic extents spanning more than 180 degrees of longitude require
an explicit projected `grid_crs`, so an automatic local projection
cannot silently split an antimeridian study area. Geographic distance
windows still work across the antimeridian. A projected grid width
describes distances in that projection; it need not be the same surface
distance everywhere over a large study area.

## Repeated station observations

``` r

observations <- sites[rep(seq_len(nrow(sites)), each = 10), ]
observations$date <- rep(as.Date("2024-01-01") + 0:9, nrow(sites))
joint <- neighbourhoods(
  observations,
  key = site, index = date,
  spatial = spatial_window(units::set_units(200, "m")),
  temporal = temporal_window(time_span(1, "days")),
  combine = "intersection"
)
tibble::as_tibble(joint, folds = c(1, 15, 60))
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                1          4         56        3           0.0667
#> 2 15               1          6         54        5           0.1   
#> 3 60               1          2         58        1           0.0333
```

Intersection omits observations close in both space and time. Union
omits observations close in either dimension. These represent different
assumptions about residual dependence. `.by` imposes hard
independent-domain boundaries in either case; keys identify stations and
need not imply independent residuals. For fixed stations, geometry
relationships are reused across times. Moving locations use each
observation’s geometry.

## Cubble

Both nested spatial and long temporal cubbles are accepted.
[`fold_data()`](https://gavinsimpson.github.io/naboer/reference/fold_data.md)
gives one canonical row per temporal observation. Its `naboer_source`
attribute maps canonical rows back to the nested spatial row and
within-series row, or to the source long row. Actual spatial `sf`
geometries are restored for temporal cubbles; polygon geometry is
preserved. The active geometry is preserved even when the spatial
component contains several `sfc` columns. Missing coordinates do not
prevent a purely temporal design: these rows remain in
[`fold_data()`](https://gavinsimpson.github.io/naboer/reference/fold_data.md).
Spatial constructors require complete, valid geometries and a CRS.

Cubble 1.0.0 itself accepts one key column. For data naturally
identified by several columns, ordinary `sf` objects or tibbles support
composite `key` and `.by` selections directly, or a composite identifier
can be created before constructing a cubble.

``` r

values <- tibble::tibble(
  site = rep(letters[1:6], each = 10),
  date = rep(as.Date("2024-01-01") + 0:9, 6), value = seq_len(60)
)
cube <- cubble::make_cubble(spatial = sites, temporal = values, key = site, index = date)
cube_folds <- neighbourhoods(cube,
  spatial = spatial_window(units::set_units(200, "m")),
  temporal = temporal_window(time_span(1, "days")), combine = "intersection"
)
head(fold_data(cube_folds))
#> # A tibble: 6 × 6
#>   site  date       value     x     y    geometry
#>   <chr> <date>     <int> <dbl> <dbl> <POINT [m]>
#> 1 a     2024-01-01     1     0     0       (0 0)
#> 2 a     2024-01-02     2     0     0       (0 0)
#> 3 a     2024-01-03     3     0     0       (0 0)
#> 4 a     2024-01-04     4     0     0       (0 0)
#> 5 a     2024-01-05     5     0     0       (0 0)
#> 6 a     2024-01-06     6     0     0       (0 0)
```

Use `spatiotemporal_folds(mode = "crossed")` for joint block holdouts,
`mode = "space"` for holding locations out across all times, and
`mode = "time"` for holding periods out across stations.
