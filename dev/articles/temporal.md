# Multiple time series and forecasting

Keys identify series. `.by` limits dependence to independent domains.
Explicitly choose whether temporal windows span all series in a domain
or stay inside each series. Ordinary tibbles need explicit metadata;
tsibbles supply it.

``` r

x <- tibble::tibble(
  site = rep(c("a", "b"), each = 12),
  date = rep(as.Date("2024-01-01") + 0:11, 2)
)
f <- neighbourhoods(x,
  key = site, index = date, within = "key",
  temporal = temporal_window(time_span(2, "days"))
)
tibble::as_tibble(f, folds = 1:3)
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                1          3         21        2            0.125
#> 2 2                1          4         20        3            0.167
#> 3 3                1          5         19        4            0.208
```

Elapsed days are exactly 86,400 seconds. Calendar days follow local
civil time and can span 23 or 25 hours across daylight-saving changes.
Calendar months/years clamp invalid month-end dates.
Ambiguous/nonexistent times error unless explicitly resolved with
[`time_span()`](https://gavinsimpson.github.io/naboer/dev/reference/time_span.md)
arguments. Numeric indices are seconds. Tsibble year-week/month/quarter
indices use period-start dates for queries;
[`fold_data()`](https://gavinsimpson.github.io/naboer/dev/reference/fold_data.md)
preserves their original class.

Observation-count spans count actual rows and require `within = "key"`.
They do not infer missing observations. Resolve tied times within a
counted series first. Duration windows can include simultaneous
observations.

``` r

blocks <- temporal_folds(x, time_span(4, "days"),
  key = site, index = date,
  within = "domain", origin = as.Date("2024-01-01")
)
tibble::as_tibble(blocks)
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                8          8         16        0            0.333
#> 2 2                8          8         16        0            0.333
#> 3 3                8          8         16        0            0.333
```

Windows include both boundaries; adjacent blocks are half-open. Default
block origins are the earliest observation in each selected domain or
series.

## Forecasting

An origin defines what was available for fitting. Exclude the gap, then
score the assessment horizon. All later observations remain omitted,
including those beyond the horizon. Origins span all participating
series in a domain by default.

``` r

forecast <- forecast_folds(x,
  initial = time_span(3, "days"),
  horizon = time_span(2, "days"), gap = time_span(1, "days"),
  key = site, index = date
)
tibble::as_tibble(forecast)
#> # A tibble: 3 × 6
#>   fold  n_assessment n_omission n_analysis n_buffer omitted_fraction
#>   <chr>        <dbl>      <dbl>      <dbl>    <dbl>            <dbl>
#> 1 1                4         16          8       12            0.667
#> 2 2                4         12         12        8            0.5  
#> 3 3                4          8         16        4            0.333
analysis_indices(forecast, 1)
#> [[1]]
#> [1]  1  2  3  4 13 14 15 16
assessment_indices(forecast, 1)
#> [[1]]
#> [1]  6  7 18 19
```

Use `window = "sliding"` with `lookback` for finite training history.
Generated origins require complete horizons; empty generated assessment
windows are skipped. Explicit `origins` require nonempty training and
assessment in every participating domain. Construct predictors and
imputations using only information available at each origin in
explicit-refit workflows.

``` r

series <- tsibble::as_tsibble(x, key = site, index = date)
temporal_folds(series, time_span(4, "days"), within = "key")
#> <naboer_folds: temporal>
#> 24 observations; 6 folds
```
