# naboer 0.0.1.9000

## New features

* Added cross-validation folds with separate assessment and omission
  membership, validation, diagnostics, tibble views and plots. Designs export
  complete `nei` structures for mgcv NCV/QNCV, including optional neighbourhood
  sampling and covariance settings.
* Added spatial, temporal and spatiotemporal neighbourhoods and block designs,
  asymmetric elapsed/calendar/count windows, buffers and rolling-origin
  forecasts across multiple series and independent domains.
* Added grouped and hierarchical holdouts, environmental clustering and strata,
  phylogenetic distance neighbourhoods and clade/cluster holdouts, and faithful
  rsample import/export.
* Added optional sf, tsibble and cubble adapters, indexed spatial queries,
  compressed interval storage, allocation guards and reproducible R benchmarks.

## User visible changes

* The minimum supported R version is now 4.4.0.
