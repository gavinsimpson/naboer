# Temporal neighbourhoods and forecasting

``` r

library("naboer")
library("tibble")
library("mgcv")
library("datasets")

have_gratia <- requireNamespace("gratia", quietly = TRUE)
if (have_gratia) {
  library("gratia")
  library("ggplot2")
}

have_tsibble <- requireNamespace("tsibble", quietly = TRUE)
if (have_tsibble) {
  library("tsibble")
}
```

## Choose smoothness for a time series

We want a smooth that describes how soil moisture changed through 2019
at Balruddery, near Dundee. Soil moisture on successive days is often
similar. A leave-one-out fit can use neighbouring days to predict a
held-out day, which can favour a smooth that follows short-term
fluctuations. We will compare that fit with one that also leaves out
nearby days during smoothness selection.

Load the 365 daily measurements from the COSMOS-UK monitoring network
with [`data()`](https://rdrr.io/r/utils/data.html). The `date` column
defines the neighbourhoods. The number of days since 1 January is the
predictor in the model, and `vwc` is the volumetric water content as a
percentage. Contains data supplied by UK Centre for Ecology & Hydrology.
[COSMOS-UK data and
licence](https://doi.org/10.5285/dde10e2b-ee1a-4b2b-a247-176db2f7bb31).

``` r

data("balruddery_2019", package = "naboer")
soil <- balruddery_2019
soil$day <- as.numeric(soil$date - min(soil$date))
```

Make one fold per day. In each fold, leave out the target day and up to
20 days on either side, and score the prediction for the target day
alone. This window follows the example in David L. Miller’s article. The
appropriate window for another series depends on how long residual
dependence persists. [Miller’s time-series
example](https://calgary.converged.yt/articles/ncv_timeseries.html#what-about-restricting-the-prediction-neighbourhood).

``` r

soil_folds <- neighbourhoods(soil,
  index = date, temporal = temporal_window(time_span(20, "days"))
)
```

Check a fold in the middle of the year. To score 2 July, the fit leaves
out 12 June to 22 July inclusive.

``` r

july_fold <- match(as.Date("2019-07-02"), soil$date)
range(soil$date[omission_indices(soil_folds, july_fold)[[1]]])
#> [1] "2019-06-12" "2019-07-22"
```

Fit the same smooth twice with mgcv. With `nei` omitted,
`method = "NCV"` uses leave-one-out cross-validation. For the second
fit, export the temporal folds with
[`as_nei()`](https://gavinsimpson.github.io/naboer/reference/as_nei.md)
and supply that design as `nei`. Both fits use the same data and basis,
so the change is in how the smoothing parameter is selected.

``` r

soil_loo <- gam(vwc ~ s(day, k = 20),
  data = soil, method = "NCV", na.action = na.fail
)

soil_nei <- as_nei(soil_folds)
soil_ncv <- gam(vwc ~ s(day, k = 20),
  data = soil, method = "NCV", nei = soil_nei, na.action = na.fail
)
```

Use gratia’s
[`edf()`](https://gavinsimpson.github.io/gratia/reference/edf.html) to
compare the effective degrees of freedom (EDF) of the smooth in each
fit. Higher EDF allows the curve to follow more fluctuations in the
data.

``` r

c(
  `Leave-one-out` = edf(soil_loo)$.edf,
  `Temporal neighbourhoods` = edf(soil_ncv)$.edf
)
#>           Leave-one-out Temporal neighbourhoods 
#>               17.397113                3.171146
```

Use gratia’s
[`add_fitted()`](https://gavinsimpson.github.io/gratia/reference/add_fitted.html)
to add each model’s fitted moisture values to the observed dates. Both
final models use all 365 observations.

``` r

soil_curves <- rbind(
  transform(add_fitted(soil, soil_loo, value = "fitted_vwc"),
    scheme = "Leave-one-out"
  ),
  transform(add_fitted(soil, soil_ncv, value = "fitted_vwc"),
    scheme = "Temporal neighbourhoods"
  )
)
```

Overlay the two fitted curves on the observations in grey. The
neighbourhood fit follows the broader changes in moisture over the year.

``` r

ggplot(soil, aes(x = date)) +
  geom_line(aes(y = vwc), colour = "grey60", linewidth = 0.4) +
  geom_line(
    data = soil_curves,
    aes(y = fitted_vwc, colour = scheme, linetype = scheme),
    linewidth = 0.9
  ) +
  scale_colour_manual(values = c(
    "Leave-one-out" = "#D55E00", "Temporal neighbourhoods" = "#0072B2"
  )) +
  scale_linetype_manual(values = c(
    "Leave-one-out" = "dashed", "Temporal neighbourhoods" = "solid"
  )) +
  labs(
    x = "Date in 2019", y = "Volumetric water content (%)",
    colour = NULL, linetype = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
```

![Daily soil moisture at Balruddery in grey and two fitted curves. The
temporal neighbourhood fit is smoother than the leave-one-out
fit.](temporal_files/figure-html/soil-curves-1.png)

The two-sided window selects a smooth for the observed period.
Forecasting uses an information cut-off, described later in this
article.

## Multiple time series

When you have several time series, decide whether observations at the
same time in different series should be held out together. For example,
two monitoring sites might respond to the same weather. You might want
to hold out both sites to assess predictions for a new period.

### Set up the time series

Create daily records for two sites. In these examples, `site` is the key
that identifies each series and `date` is the time index. Supplying a
key keeps track of the series. Use `.by` separately when the data
contain independent domains and the design should keep neighbourhoods
within them.

``` r

x <- tibble(
  site = rep(c("a", "b"), each = 12),
  date = rep(as.Date("2024-01-01") + 0:11, 2)
)
```

### Leave out observations near each target

Use
[`neighbourhoods()`](https://gavinsimpson.github.io/naboer/reference/neighbourhoods.md)
to assess one row per fold. A two-day window leaves out the target and
any observations up to two days before or after it. Here
`within = "key"` restricts that window to the target’s site. Use
`within = "domain"` to include nearby dates from both sites when they
belong to the same domain.

``` r

f <- neighbourhoods(x,
  key = site, index = date, within = "key",
  temporal = temporal_window(time_span(2, "days"))
)
```

Inspect the rows omitted in fold 3. The target is 3 January at site a,
so the model leaves out 1 to 5 January at that site and scores only 3
January.

``` r

x[omission_indices(f, 3)[[1]], ]
#> # A tibble: 5 × 2
#>   site  date      
#>   <chr> <date>    
#> 1 a     2024-01-01
#> 2 a     2024-01-02
#> 3 a     2024-01-03
#> 4 a     2024-01-04
#> 5 a     2024-01-05
```

### Choose how to measure a span

- **Elapsed time** is measured in seconds. An elapsed day is 86,400
  seconds. Numeric time indices are also interpreted in seconds.
- **Calendar time** follows the local calendar. Across a daylight-saving
  change, a calendar day can be 23 or 25 hours. When shifting by months
  or years, a date that would fall past the end of a month moves to its
  last day. Calendar shifts preserve the time of day, including
  fractional seconds, subject to the chosen daylight-saving rules.
  naboer uses `clock` for this calendar arithmetic.
- **Observation counts** measure a span in recorded rows. Blocks start
  at the first row of each series, and neighbourhoods count away from
  the target. Use `time_span(2, "observations")` for a span of two rows
  and set `within = "key"`. Resolve tied times in each series before
  using counts. Time-based windows can include observations with the
  same timestamp.

Temporal windows include observations on both boundaries. Elapsed
windows allow for rounding when calculating time differences at a
boundary. When calendar arithmetic would reach a repeated or skipped
local clock time, the calculation fails by default. Use the `ambiguous`
or `nonexistent` argument of
[`time_span()`](https://gavinsimpson.github.io/naboer/reference/time_span.md)
to choose how to resolve those times. A zero-length shift preserves the
original timestamp, including which occurrence of a repeated time it
represents. For date-times, calendar arithmetic uses the time zone of
the index even if supplied origins are displayed in a different zone.

## Hold out whole periods

To assess a new period at both sites, use
[`temporal_folds()`](https://gavinsimpson.github.io/naboer/reference/temporal_folds.md)
with `within = "domain"`. Start the four-day blocks on 1 January.

``` r

blocks <- temporal_folds(x, time_span(4, "days"),
  key = site, index = date,
  within = "domain", origin = as.Date("2024-01-01")
)
```

The fold summary shows three blocks of eight assessment rows each, four
from each site.

``` r

as_tibble(blocks)
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                8          8         16        0            0.333
#> 2 2                8          8         16        0            0.333
#> 3 3                8          8         16        0            0.333
```

A time-based block includes its start and excludes its end, so 5 January
belongs to the second block. If you omit `origin`, each domain or series
starts at its earliest observation. Count-based blocks start at the
first observation of each series and take no `origin` argument.

## Assess forecasts using information available at the time

A forecasting fold has a cut-off called the *origin*. Here the model can
use records through the origin. The gap is the period after the origin
when no predictions are scored. The horizon is the period immediately
after the gap when predictions are scored. Later records, including
those beyond the horizon, are left out of training.

Generate the first origin three elapsed days after the first
observation. Move the origin forward by two days at a time, the default
step when the horizon is two days. With multiple series in one domain,
time-based forecasts use a shared origin so neither site’s future
records enter the fit.

``` r

forecast <- forecast_folds(x,
  initial = time_span(3, "days"),
  horizon = time_span(2, "days"), gap = time_span(1, "days"),
  key = site, index = date
)
```

Inspect the first fold. Both sites contribute training records from 1 to
4 January. The records on 5 January form the gap, and the predictions
for 6 and 7 January are scored.

``` r

x[analysis_indices(forecast, 1)[[1]], ]
#> # A tibble: 8 × 2
#>   site  date      
#>   <chr> <date>    
#> 1 a     2024-01-01
#> 2 a     2024-01-02
#> 3 a     2024-01-03
#> 4 a     2024-01-04
#> 5 b     2024-01-01
#> 6 b     2024-01-02
#> 7 b     2024-01-03
#> 8 b     2024-01-04
x[assessment_indices(forecast, 1)[[1]], ]
#> # A tibble: 4 × 2
#>   site  date      
#>   <chr> <date>    
#> 1 a     2024-01-06
#> 2 a     2024-01-07
#> 3 b     2024-01-06
#> 4 b     2024-01-07
```

Training starts with the first observation by default. For a fixed
training history, set `window = "sliding"` and supply a `lookback` span.
Training then uses times strictly after `origin - lookback` and through
the origin.

Generated origins require the horizon to end on or before the last
observed time in the domain. Series can end earlier or have gaps, and
generated horizons with no observations to assess are skipped. To choose
cut-off dates yourself, supply `origins` in place of `initial`. Each
explicit origin must leave at least one training and one assessment row
in every domain.

When both `initial` and `step` are calendar spans measured in months or
years, origins are calculated from the first observation’s date to keep
month-end and leap-day adjustments from carrying over to later origins.
The same anchoring applies when both spans use calendar days or weeks.
The daylight-saving rules in `initial` apply to the first origin, and
those in `step` apply to later origins.

### Forecast across separate domains

Observation-count forecasts require a single series per domain. For
related series recorded at different times, use elapsed or calendar
spans so there is a common time cut-off. Forecasting with
`within = "key"` also requires a single series per domain.

Setting `.by = site` gives each site its own forecast schedule and
leaves records from other sites available for fitting, including their
later records. If the same cut-off should apply to both sites, keep them
in one domain as above.

### Refit models at forecast origins

Use the `nottem` series to forecast monthly mean air temperatures at
Nottingham Castle. The record covers 1920 to 1939 and reports
temperatures in degrees Fahrenheit. Date each observation at the end of
its month. [R documentation for
`nottem`](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/nottem.html).

``` r

nottingham <- tibble(
  date = seq(as.Date("1920-02-01"), by = "month", length.out = length(nottem)) - 1,
  month = as.integer(cycle(nottem)),
  temperature = as.numeric(nottem)
)
```

Set forecast origins at the end of 1929, 1931 and 1933. Each fold trains
on observations through December of that year and assesses all 12 months
of the following year. Get the row numbers for fitting and scoring each
model.

``` r

cutoffs <- as.Date(c("1929-12-31", "1931-12-31", "1933-12-31"))
nottingham_folds <- forecast_folds(nottingham,
  index = date, origins = cutoffs,
  horizon = time_span(1, "years", type = "calendar")
)
training_rows <- analysis_indices(nottingham_folds)
assessment_rows <- assessment_indices(nottingham_folds)
```

Fit a seasonal GAM using only the training rows. The month of the year
is known in advance, so the model needs no temperature records beyond
the cut-off to make its forecasts. A cyclic smooth connects December to
January, with knots at 0.5 and 12.5 defining the ends of the annual
cycle. Score each set of forecasts with root mean squared error (RMSE).

``` r

forecast_rmse <- function(training, assessment) {
  fit <- gam(temperature ~ s(month, bs = "cc", k = 8),
    data = training, method = "REML", knots = list(month = c(0.5, 12.5))
  )
  predicted <- predict(fit, newdata = assessment["month"])
  sqrt(mean((assessment$temperature - predicted)^2))
}
```

Refit the GAM separately for each forecast origin. Each row of the
result summarises forecast errors, in degrees Fahrenheit, over the
following year.

``` r

rmse <- vapply(seq_along(training_rows), function(i) {
  forecast_rmse(
    nottingham[training_rows[[i]], ],
    nottingham[assessment_rows[[i]], ]
  )
}, numeric(1))
tibble(origin = cutoffs, rmse_f = rmse)
#> # A tibble: 3 × 2
#>   origin     rmse_f
#>   <date>      <dbl>
#> 1 1929-12-31   1.63
#> 2 1931-12-31   1.87
#> 3 1933-12-31   2.71
```

When forecasting with other predictors, use only information available
at the origin to construct them or to fill in missing values. Forecast
folds often leave out many observations. If you use neighbourhood
cross-validation (NCV) in place of explicit refitting, compare its
approximation with exact refits for your model and forecast folds. See
the [comparison with exact
refits](https://gavinsimpson.github.io/naboer/articles/naboer.html#check-the-approximation-for-large-omissions).

## Use a tsibble

With an ordinary tibble, supply `key` and `index` explicitly. With a
tsibble, naboer reads them from the object. If you have `tsibble`
installed, convert the example data and use `within = "key"` to make
four-day blocks separately for each site.

``` r

series <- as_tsibble(x, key = site, index = date)
temporal_folds(series, time_span(4, "days"), within = "key")
#> <naboer_folds: temporal>
#> 24 observations; 6 folds
```

For tsibble indices that represent weeks, months or quarters, temporal
calculations use the first date of each period. The supplied tsibble
keeps its original index class.
