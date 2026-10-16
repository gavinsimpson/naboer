# Specify a spatial neighbourhood

Distances are measured between the original geometries, including
polygons. The current sf spherical/ellipsoidal setting is respected and
never changed.

## Usage

``` r
spatial_window(radius)
```

## Arguments

- radius:

  A nonnegative scalar with physical distance units, for example
  `units::set_units(5, "km")`.

## Value

A spatial window specification.
