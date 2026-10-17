# Export mgcv neighbourhood structures

The result is a plain list of one-based indices, with cumulative
endpoints. No rescaling of repeated assessments or covariance settings
is performed. Neighbourhood sampling preserves the complete structure
needed by bam. Export requires at least two folds for compatibility with
the tested mgcv integration. Single native folds remain usable for
inspection and rsample. In mgcv 1.9-4, scalar sampling incorrectly
samples packed endpoints for nonsingleton omissions. For these designs,
supply at least two explicit fold numbers, or omit `sample` to use every
fold. This compatibility restriction is enforced without requiring mgcv
to be installed.

## Usage

``` r
as_nei(
  x,
  sample = NULL,
  jackknife = NULL,
  max_bytes = getOption("naboer.max_bytes", 1024^3),
  data = NULL
)
```

## Arguments

- x:

  A `naboer_folds` object.

- sample:

  Optional bam sampling specification: one positive count, or a vector
  of neighbourhood numbers. A length-one value always means a count and
  is supported only when every omission set is a singleton.

- jackknife:

  Optional logical gam covariance setting. `NULL` preserves mgcv's
  default by omitting the field.

- max_bytes:

  Maximum estimated packed output size in bytes; defaults to option
  `naboer.max_bytes`, or 1 GiB. Use `Inf` to disable the memory guard.

- data:

  Optional data used for alignment checking. Supply unchanged canonical
  data, or data containing the original unique identifier column.

## Value

A plain list with `a`, `ma`, `d`, `md`, and optional settings.
