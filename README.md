# naboer

[![R-CMD-check](https://github.com/gavinsimpson/naboer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gavinsimpson/naboer/actions/workflows/R-CMD-check.yaml)
[![Coverage](https://codecov.io/gh/gavinsimpson/naboer/branch/main/graph/badge.svg)](https://app.codecov.io/gh/gavinsimpson/naboer)

naboer constructs dependence-aware cross-validation designs for ecological and
environmental data. It distinguishes observations **assessed** by a fold from
observations **omitted** during fitting, including dependence buffers. Export
these designs to `mgcv` for NCV/QNCV smoothness selection, or to `rsample` for
explicit resampling workflows.

This is the 0.1.0 release candidate, pending release approval. `tools/RELEASE.md`
records completed validation and remaining release gates. It has not been tagged
or submitted to CRAN.

## Temporal neighbourhoods

```r
library("naboer")

observations <- tibble::tibble(
  site = rep(c("north", "south"), each = 30),
  time = rep(1:30, 2),
  response = sin(time / 5) + rep(c(0, 0.5), each = 30)
)

folds <- observations |>
  neighbourhoods(
    key = site, index = time, within = "key",
    temporal = temporal_window(time_span(2, "seconds"))
  )
tibble::as_tibble(folds, folds = 1:3, membership = TRUE)
nei <- as_nei(folds)

if (requireNamespace("mgcv", quietly = TRUE)) {
  fit <- mgcv::gam(
    response ~ s(time, k = 8) + site,
    data = fold_data(folds), method = "NCV", nei = nei,
    na.action = stats::na.fail
  )
}
```

Numeric time indices are seconds. For dates, use `Date`/`POSIXct` with elapsed
or calendar spans. `time_span(1, "months", "calendar")` follows the calendar;
`time_span(3, "observations")` counts rows within an ordered series. Tsibble
keys/index and cubble metadata are inferred automatically.

## Space and time together

```r
if (requireNamespace("sf", quietly = TRUE) &&
    requireNamespace("units", quietly = TRUE)) {
  stations <- sf::st_as_sf(
    tibble::tibble(site = c("a", "b", "c"), x = c(0, 1000, 5000), y = 0),
    coords = c("x", "y"), crs = 32632
  )
  space <- neighbourhoods(
    stations, spatial = spatial_window(units::set_units(2, "km"))
  )
  observations <- stations[rep(seq_len(nrow(stations)), each = 10), ]
  observations$date <- rep(as.Date("2026-01-01") + 0:9, nrow(stations))
  space_time <- observations |>
    neighbourhoods(
      key = site, index = date,
      spatial = spatial_window(units::set_units(2, "km")),
      temporal = temporal_window(time_span(1, "days")),
      combine = "intersection"
    )
}
```

Supply both window arguments and explicitly choose `combine = "intersection"`
or `"union"` for spatiotemporal neighbourhoods. Keys identify stations or
series; `.by` defines independent domains. Multiple temporal series require
`within = "key"` or `within = "domain"`. Neither keys nor `group_by()` silently
declare independence.
For forecasting, use shared elapsed/calendar cutoffs across dependent series.
Count-based horizons require one series per explicitly independent domain.

## Choose a prediction task

| Task | Constructors |
|---|---|
| Predict observations beyond local dependence | `neighbourhoods()`, `buffer_folds()` |
| Predict new places or periods | `spatial_folds()`, `temporal_folds()`, `spatiotemporal_folds()` |
| Forecast using only information available at an origin | `forecast_folds()` |
| Predict new groups, environments or taxa | `group_folds()`, `environmental_folds()`, `phylogenetic_folds()` |
| Supply a custom or external design | `custom_folds()`, `from_rset()` |

Inspect with `summary()`, `tibble::as_tibble()` and `ggplot2::autoplot()`. Large
memberships remain compressed until requested. Export estimates packed size and
defaults to a 1 GiB limit, controlled by `max_bytes` or `naboer.max_bytes`.
Million-observation sparse designs are benchmarked; wide windows can still
require enormous output. `as_rset()` also guards analysis-complement expansion,
which can be much larger than `nei` export.

Prepare complete modelling data before constructing folds. Rebuild the design
after changing observations or their order. A valid design specifies a prediction
task; it does not establish NCV approximation accuracy or automatically select
an appropriate dependence radius.

R >= 4.4 is required. Spatial, temporal-class, phylogenetic, plotting and rsample
integrations use optional dependencies. See the [package site](https://gavinsimpson.github.io/naboer/)
for worked workflows and the statistical contract.
