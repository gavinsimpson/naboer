# Package index

## Neighbourhoods and windows

- [`neighbourhoods()`](https://gavinsimpson.github.io/naboer/dev/reference/neighbourhoods.md)
  : Construct focal-observation neighbourhoods
- [`spatial_window()`](https://gavinsimpson.github.io/naboer/dev/reference/spatial_window.md)
  : Specify a spatial neighbourhood
- [`temporal_window()`](https://gavinsimpson.github.io/naboer/dev/reference/temporal_window.md)
  : Describe an asymmetric temporal neighbourhood
- [`time_span()`](https://gavinsimpson.github.io/naboer/dev/reference/time_span.md)
  : Describe a temporal span
- [`phylogenetic_window()`](https://gavinsimpson.github.io/naboer/dev/reference/phylogenetic_window.md)
  : Specify a phylogenetic neighbourhood
- [`buffer_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/buffer_folds.md)
  : Add dependence buffers around fold assessment observations

## Cross-validation strategies

- [`custom_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/custom_folds.md)
  : Construct custom cross-validation folds
- [`loo_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/loo_folds.md)
  : Leave-one-observation-out folds
- [`vfold_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/vfold_folds.md)
  : Random V-fold cross-validation
- [`spatial_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/spatial_folds.md)
  : Construct spatial cross-validation folds
- [`temporal_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/temporal_folds.md)
  : Construct contiguous temporal blocks
- [`spatiotemporal_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/spatiotemporal_folds.md)
  : Construct spatiotemporal block folds
- [`forecast_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/forecast_folds.md)
  : Construct rolling-origin forecasting folds
- [`group_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/group_folds.md)
  : Hold out observational groups
- [`environmental_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/environmental_folds.md)
  : Construct environmental cross-validation folds
- [`phylogenetic_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/phylogenetic_folds.md)
  : Construct phylogenetic cross-validation folds

## Inspect and export

- [`fold_data()`](https://gavinsimpson.github.io/naboer/dev/reference/fold_data.md)
  : Retrieve the canonical observation data
- [`assessment_indices()`](https://gavinsimpson.github.io/naboer/dev/reference/assessment_indices.md)
  [`omission_indices()`](https://gavinsimpson.github.io/naboer/dev/reference/assessment_indices.md)
  [`analysis_indices()`](https://gavinsimpson.github.io/naboer/dev/reference/assessment_indices.md)
  : Retrieve observation indices for folds
- [`as_tibble(`*`<naboer_folds>`*`)`](https://gavinsimpson.github.io/naboer/dev/reference/as_tibble.naboer_folds.md)
  : Inspect folds as a tibble
- [`summary(`*`<naboer_folds>`*`)`](https://gavinsimpson.github.io/naboer/dev/reference/summary.naboer_folds.md)
  : Summarise a cross-validation design
- [`autoplot(`*`<naboer_folds>`*`)`](https://gavinsimpson.github.io/naboer/dev/reference/autoplot.naboer_folds.md)
  : Plot selected cross-validation folds
- [`as_nei()`](https://gavinsimpson.github.io/naboer/dev/reference/as_nei.md)
  : Export mgcv neighbourhood structures
- [`validate_folds()`](https://gavinsimpson.github.io/naboer/dev/reference/validate_folds.md)
  : Validate fold membership
- [`validate_nei()`](https://gavinsimpson.github.io/naboer/dev/reference/validate_nei.md)
  : Validate a packed mgcv neighbourhood structure
- [`as_rset()`](https://gavinsimpson.github.io/naboer/dev/reference/as_rset.md)
  [`from_rset()`](https://gavinsimpson.github.io/naboer/dev/reference/as_rset.md)
  : Convert between naboer folds and rsample resamples

## Package

- [`naboer`](https://gavinsimpson.github.io/naboer/dev/reference/naboer-package.md)
  [`naboer-package`](https://gavinsimpson.github.io/naboer/dev/reference/naboer-package.md)
  : Dependence-aware cross-validation neighbourhoods
