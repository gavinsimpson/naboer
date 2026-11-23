# Export mgcv neighbourhood structures

The result is a plain list of one-based indices, with cumulative
endpoints. Repeated assessments retain their original weight.
Neighbourhood sampling preserves the complete structure needed by bam.
Export requires at least two folds for compatibility with the tested
mgcv integration. Single native folds remain usable for inspection and
rsample.

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
  of input fold numbers. A length-one value always means a count.

- jackknife:

  Optional logical gam covariance setting. `NULL` preserves mgcv's
  default by omitting the field. On versions that require focal folds in
  observation order, `TRUE` requires exactly one singleton assessment
  fold per observation.

- max_bytes:

  Maximum estimated packed output size in bytes; defaults to option
  `naboer.max_bytes`, or 1 GiB. Use `Inf` to disable the memory guard.

- data:

  Optional modelling data for a row-count check. Supply the data in the
  observation order used to build the folds. For cubbles, use the
  unfolded observations returned by
  [`fold_data()`](https://gavinsimpson.github.io/naboer/reference/fold_data.md).

## Value

A plain list with `a`, `ma`, `d`, `md`, and optional settings.

## Details

For compatibility with covariance calculations in affected versions of
mgcv, designs with one singleton assessment fold per observation are
packed in observation order. Explicit sample indices are translated from
input fold numbers to packed fold numbers. On these versions, a sample
containing every fold is exported without a sampling field.

For nonsingleton omissions, sampling by count requires mgcv 1.9-5 or
later because 1.9-4 samples packed endpoints instead of fold numbers.
With 1.9-4, use at least two explicit fold numbers or use all folds.
