# naboer 0.1.0

## New features

* Added cross-validation folds with separate assessment and omission
  membership, validation, diagnostics, tibble views and plots. Designs export
  complete `nei` structures for mgcv NCV/QNCV, including optional neighbourhood
  sampling and covariance settings. Compatibility checks require at least two
  folds and explicit sampled fold numbers for nonsingleton omissions under the
  tested mgcv 1.9-4 integration.
* Added spatial, temporal and spatiotemporal neighbourhoods and block designs,
  asymmetric elapsed/calendar/count windows, buffers and rolling-origin
  forecasts across multiple series and independent domains. Shared forecasting
  cutoffs exclude future observations across dependent series; calendar shifts
  preserve zero-width instants and use the index time zone.
* Added grouped and hierarchical holdouts, environmental clustering and strata,
  phylogenetic distance neighbourhoods and clade/cluster holdouts, and faithful
  rsample import/export with preserved key/domain metadata and guarded analysis
  expansion. Imports require resamples indexing individual observations.
* Added optional sf, tsibble and cubble adapters, indexed spatial queries,
  compressed interval storage, allocation guards and reproducible R benchmarks.
  Explicit spatial grid anchors remain stable across data subsets, and radius
  queries consistently use stored coordinates. The approximation study includes
  large-omission stress tests and independent weighted Gaussian refits.

## User visible changes

* The minimum supported R version is now 4.4.0.
