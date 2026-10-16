# Describe a temporal span

Elapsed spans use seconds (a day is 86400 seconds); calendar spans use
local calendar arithmetic. Calendar month and year arithmetic clamps
invalid dates to the final day of the month. Count spans count ordered
observations. Numeric time indices are interpreted in seconds. Tsibble
year-week, year-month and year-quarter indices use the first date of
each period for elapsed and calendar arithmetic; source index columns
are preserved.

## Usage

``` r
time_span(
  n,
  unit = "seconds",
  type = "elapsed",
  ambiguous = "error",
  nonexistent = "error"
)
```

## Arguments

- n:

  Nonnegative span length. Calendar and count spans require integers.

- unit:

  A time unit, or `"observations"` for a count span.

- type:

  One of `"elapsed"`, `"calendar"`, or `"count"` (also accepted as
  `"observations"`). The default for `unit = "observations"` is
  `"count"`.

- ambiguous, nonexistent:

  Resolution passed to
  [`clock::add_days()`](https://clock.r-lib.org/reference/clock-arithmetic.html)
  for ambiguous or nonexistent local times; both default to `"error"`.

## Value

A temporal span specification.
