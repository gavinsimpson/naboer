# naboer

[![R-CMD-check](https://github.com/gavinsimpson/naboer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gavinsimpson/naboer/actions/workflows/R-CMD-check.yaml)
[![Coverage](https://codecov.io/gh/gavinsimpson/naboer/branch/main/graph/badge.svg)](https://app.codecov.io/gh/gavinsimpson/naboer)

naboer constructs dependence-aware cross-validation designs for
ecological and environmental data. It distinguishes observations
**assessed** by a fold from observations **omitted** during fitting,
including dependence buffers. Export these designs to `mgcv` for
NCV/QNCV smoothness selection, or to `rsample` for explicit resampling
workflows.

This is the 0.1.0 release candidate, pending release approval.
`tools/RELEASE.md` records completed validation and remaining release
gates. It has not been tagged or submitted to CRAN.

## Temporal neighbourhoods

``` r

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

fit <- mgcv::gam(
  response ~ s(time, k = 8) + site,
  data = observations, method = "NCV", nei = nei,
  na.action = stats::na.fail
)
```

Keep the observation table separately from the fold object. For model
fitting, use the same observations in the same row order as the folds.
If you supply `data` to
[`as_nei()`](https://gavinsimpson.github.io/naboer/reference/as_nei.md),
it checks the row count.

Numeric time indices are seconds. For dates, use `Date`/`POSIXct` with
elapsed or calendar spans. `time_span(1, "months", "calendar")` follows
the calendar; `time_span(3, "observations")` counts rows within an
ordered series. Tsibble keys/index and cubble metadata are inferred
automatically.

## Space and time together

Install `sf` for site geometries. Use `units`, which comes with a
standard naboer installation, to specify distances for windows and
buffers. For spatiotemporal data, use a cubble to link repeated
measurements to site geometries. This example also requires `cubble`.

``` r

stations <- sf::st_as_sf(
  tibble::tibble(site = c("a", "b", "c"), x = c(0, 1000, 5000), y = 0),
  coords = c("x", "y"), crs = 32632
)
space <- neighbourhoods(
  stations, spatial = spatial_window(units::set_units(2, "km"))
)
measurements <- tibble::tibble(
  site = rep(stations$site, each = 10),
  date = rep(as.Date("2026-01-01") + 0:9, nrow(stations)),
  value = seq_len(nrow(stations) * 10)
)
observations <- cubble::make_cubble(
  spatial = stations, temporal = measurements, key = site, index = date
)
space_time <- observations |>
  neighbourhoods(
    spatial = spatial_window(units::set_units(2, "km")),
    temporal = temporal_window(time_span(1, "days")),
    combine = "intersection"
  )
```

Use `fold_data(space_time, observations)` to unfold the cubble into the
row order used by the folds for modelling.

Supply both window arguments and explicitly choose
`combine = "intersection"` or `"union"` for spatiotemporal
neighbourhoods. Keys identify stations or series; `.by` defines
independent domains. Multiple temporal series require `within = "key"`
or `within = "domain"`. Neither keys nor `group_by()` silently declare
independence. For forecasting, use shared elapsed/calendar cutoffs
across dependent series. Count-based horizons require one series per
explicitly independent domain.

## Choose a prediction task

| Task | Constructors |
|----|----|
| Predict observations beyond local dependence | [`neighbourhoods()`](https://gavinsimpson.github.io/naboer/reference/neighbourhoods.md), [`buffer_folds()`](https://gavinsimpson.github.io/naboer/reference/buffer_folds.md) |
| Predict new places or periods | [`spatial_folds()`](https://gavinsimpson.github.io/naboer/reference/spatial_folds.md), [`temporal_folds()`](https://gavinsimpson.github.io/naboer/reference/temporal_folds.md), [`spatiotemporal_folds()`](https://gavinsimpson.github.io/naboer/reference/spatiotemporal_folds.md) |
| Forecast using only information available at an origin | [`forecast_folds()`](https://gavinsimpson.github.io/naboer/reference/forecast_folds.md) |
| Predict new groups, environments or taxa | [`group_folds()`](https://gavinsimpson.github.io/naboer/reference/group_folds.md), [`environmental_folds()`](https://gavinsimpson.github.io/naboer/reference/environmental_folds.md), [`phylogenetic_folds()`](https://gavinsimpson.github.io/naboer/reference/phylogenetic_folds.md) |
| Supply a custom or external design | [`custom_folds()`](https://gavinsimpson.github.io/naboer/reference/custom_folds.md), [`from_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md) |

Inspect with [`summary()`](https://rdrr.io/r/base/summary.html),
[`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
and
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html).
Supply `data` to `autoplot()` to show locations or times. By default the
plot shows observation positions. Pass the original observations to
[`buffer_folds()`](https://gavinsimpson.github.io/naboer/reference/buffer_folds.md)
and
[`as_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md).
Large memberships remain compressed until requested. Export estimates
packed size and defaults to a 1 GiB limit, controlled by `max_bytes` or
`naboer.max_bytes`. Million-observation sparse designs are benchmarked;
wide windows can still require enormous output.
[`as_rset()`](https://gavinsimpson.github.io/naboer/reference/as_rset.md)
also guards analysis-complement expansion, which can be much larger than
`nei` export.

Prepare complete modelling data before constructing folds. Rebuild the
design after changing observations or their order. A valid design
specifies a prediction task; it does not establish NCV approximation
accuracy or automatically select an appropriate dependence radius.

R 4.4 or later is required. A standard installation includes `mgcv` for
model fitting, `clock` for calendar spans and `units` for distances.
Install `sf` for spatial neighbourhoods, `cubble` for observations
stored in a cubble, `tsibble` for indexed time series, `ape` for
phylogenetic neighbourhoods or `rsample` to import or export resamples.
Diagnostic plots use `ggplot2`. [Package site and worked
examples](https://gavinsimpson.github.io/naboer/).
