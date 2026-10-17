# Construct rolling-origin forecasting folds

Assessment observations are strictly after the origin and gap and no
later than the horizon endpoint. Every future observation, including
observations beyond the horizon, is omitted. Sliding windows also omit
old observations. Other independent domains remain available for
fitting. Forecasts across several series in the same domain must use
time spans and `within = "domain"`, so that every series shares the
cutoff. Count forecasts require one series per independent domain;
`within = "key"` is rejected when it would leave another series' future
observations available in that domain.

## Usage

``` r
forecast_folds(
  data,
  initial = NULL,
  horizon,
  step = horizon,
  gap = NULL,
  window = c("expanding", "sliding"),
  lookback = NULL,
  origins = NULL,
  index = NULL,
  key = NULL,
  .by = NULL,
  id = NULL,
  within = NULL
)
```

## Arguments

- data:

  A data frame, tsibble, or cubble.

- initial:

  Initial training span, measured from the first observation. Required
  unless `origins` is supplied.

- horizon:

  Positive span assessed after each gap.

- step:

  Positive distance between generated origins; defaults to horizon.

- gap:

  Nonnegative span between training and assessment; defaults to zero in
  the same units and type as `horizon`.

- window:

  `"expanding"` or `"sliding"` training window.

- lookback:

  Positive training span for sliding windows. Training includes
  observations strictly after `origin - lookback` and through the
  origin.

- origins:

  Optional explicit index values. Generated origins require a complete
  horizon within the observed range; explicit origins require at least
  one assessment and one training observation in each domain.

- index, key, .by, id:

  Columns identifying time, series, independent domains, and unique
  observations. Use unquoted column names or tidy selections.

- within:

  Choose `"key"` or `"domain"` for multiple series. Count spans require
  `"key"` and count rows in stable time order.

## Value

A `naboer_folds` object. Generated empty assessment windows are skipped.
