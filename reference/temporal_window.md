# Describe an asymmetric temporal neighbourhood

Window boundaries are inclusive. Elapsed windows compare the time
difference to the requested span, allowing for floating-point roundoff
at the boundary.

## Usage

``` r
temporal_window(before, after = before)
```

## Arguments

- before, after:

  Temporal spans made with
  [`time_span()`](https://gavinsimpson.github.io/naboer/reference/time_span.md).

## Value

A temporal window specification.
