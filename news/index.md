# Changelog

## naboer 0.1.0

### New features

- Added cross-validation folds with separate assessment and omission
  membership, validation, diagnostics, tibble views and plots. Designs
  export complete `nei` structures for mgcv NCV/QNCV, including optional
  neighbourhood sampling and covariance settings. Export requires at
  least two folds. Sampling by count for nonsingleton omissions requires
  mgcv 1.9-5 or later.
- Added spatial, temporal and spatiotemporal neighbourhoods and block
  designs, asymmetric elapsed/calendar/count windows, buffers and
  rolling-origin forecasts across multiple series and independent
  domains. Shared forecasting cutoffs exclude future observations across
  dependent series; calendar shifts preserve zero-width instants and use
  the index time zone.
- Added grouped and hierarchical holdouts, environmental clustering and
  strata, phylogenetic distance neighbourhoods and clade/cluster
  holdouts, and faithful rsample import/export with preserved key/domain
  metadata and guarded analysis expansion. Imports require resamples
  indexing individual observations.
- Added optional sf, tsibble and cubble adapters, indexed spatial
  queries, compressed interval storage, allocation guards and
  reproducible R benchmarks. Explicit spatial grid anchors remain stable
  across data subsets, and radius queries consistently use stored
  coordinates. The approximation study includes large-omission stress
  tests and independent weighted Gaussian refits.

### User visible changes

- The minimum supported R version is now 4.4.0.
- A default installation now includes mgcv for NCV integration, clock
  for calendar windows, blocks and forecasts, and units for distance
  specifications.

### Bug fixes

- On affected mgcv releases, export complete focal designs in
  observation order for covariance calculations. Preserve which folds
  are sampled and omit the sampling field when all folds are selected.
  Reject explicit jackknife covariance for unsupported assessment
  layouts on those releases.
- Preserve geometry type and coordinate precision when reusing spatial
  neighbourhoods. Fold plots use the active geometry column.
- Keep cubble source row identifiers aligned when temporal conversion
  sorts observations, so rsample import rejects splits indexed in a
  different order.
- Preserve time of day and fractional seconds under calendar shifts. Use
  the initial span’s daylight-saving rules for the first forecast origin
  and include elapsed-window boundaries consistently in both directions.
