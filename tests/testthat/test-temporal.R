test_that("temporal spans reject invalid and incompatible specifications", {
  expect_error(time_span(-1, "days"), class = "naboer_error")
  expect_error(time_span(0.5, "months", "calendar"), "whole-number")
  expect_error(time_span(1, "months"), "not supported")
  expect_error(time_span(1, "seconds", "calendar"), "not supported")
  expect_error(time_span(1, "days", ambiguous = "guess"), "ambiguous")
  expect_equal(time_span(2, "observation")$type, "count")
  expect_equal(time_span(2, "observations", "observations")$type, "count")
  expect_error(temporal_window(1), "time_span")
  expect_error(temporal_window(time_span(1, "days"), time_span(1, "observations")), "both")
})

test_that("sorted temporal searches agree with brute-force inclusive windows", {
  x <- tibble::tibble(
    series = rep(c("a", "b"), each = 5),
    domain = rep(c("north", "south"), 5),
    time = c(8, 1, 2, 5, 3, 5, 0, 2, 6, 4)
  )
  ctx <- prepare_data(x, key = rlang::quo(series), index = rlang::quo(time), by = rlang::quo(domain))
  w <- temporal_window(time_span(2, "seconds"), time_span(1, "seconds"))
  expect_error(temporal_membership(ctx, w), "Multiple series")
  for (within in c("key", "domain")) {
    membership <- temporal_membership(ctx, w, within)
    domains <- domain_ids(ctx, within)
    expected <- lapply(seq_len(nrow(x)), function(i) {
      which(domains == domains[i] & x$time >= x$time[i] - 2 & x$time <= x$time[i] + 1)
    })
    actual <- lapply(seq_len(nrow(x)), function(i) sort(index_get(membership, i)))
    expect_equal(actual, expected)
  }
})

test_that("count windows use stable series order and reject tied timestamps", {
  x <- tibble::tibble(
    series = c("a", "b", "a", "b", "a", "b"),
    time = c(100, 6, 0, 1, 1, 30)
  )
  ctx <- prepare_data(x, key = rlang::quo(series), index = rlang::quo(time))
  w <- temporal_window(time_span(1, "observations"), time_span(0, "observations"))
  expect_error(temporal_membership(ctx, w, "domain"), "require")
  z <- temporal_membership(ctx, w, "key")
  expect_equal(sort(index_get(z, 1)), c(1L, 5L))
  expect_equal(index_get(z, 3), 3L)
  expect_equal(sort(index_get(z, 6)), c(2L, 6L))
  x$time[5] <- 0
  ctx <- prepare_data(x, key = rlang::quo(series), index = rlang::quo(time))
  expect_error(temporal_membership(ctx, w, "key"), "unique index")
})

test_that("calendar months and years clamp month ends and leap dates", {
  skip_if_not_installed("clock")
  month <- time_span(1, "months", "calendar")
  year <- time_span(1, "years", "calendar")
  expect_equal(temporal_shift(as.Date("2024-01-31"), month), as.Date("2024-02-29"))
  expect_equal(temporal_shift(as.Date("2024-02-29"), year), as.Date("2025-02-28"))
  x <- tibble::tibble(time = as.Date(c("2024-03-30", "2024-01-31", "2024-03-31", "2024-02-29")))
  ctx <- prepare_data(x, index = rlang::quo(time))
  block <- temporal_blocks(ctx, month, origin = as.Date("2024-01-31"))
  expect_equal(block[1], block[4])
  expect_false(block[3] == block[4])
  z <- temporal_membership(ctx, temporal_window(month), "domain")
  expect_true(4L %in% index_get(z, 2))
})

test_that("calendar days differ from elapsed days across DST", {
  skip_if_not_installed("clock")
  x <- as.POSIXct("2024-03-30 12:00:00", tz = "Europe/Copenhagen")
  calendar <- temporal_shift(x, time_span(1, "days", "calendar"))
  elapsed <- temporal_shift(x, time_span(1, "days"))
  expect_equal(as.numeric(difftime(calendar, x, units = "hours")), 23)
  expect_equal(as.numeric(difftime(elapsed, x, units = "hours")), 24)
  nonexistent <- as.POSIXct("2024-03-30 02:30:00", tz = "Europe/Copenhagen")
  expect_error(temporal_shift(nonexistent, time_span(1, "days", "calendar")), class = "naboer_error_time")
  resolved <- temporal_shift(nonexistent, time_span(1, "days", "calendar", nonexistent = "roll-forward"))
  expect_equal(format(resolved, "%H:%M"), "03:00")
  ambiguous <- as.POSIXct("2024-10-26 02:30:00", tz = "Europe/Copenhagen")
  expect_error(temporal_shift(ambiguous, time_span(1, "days", "calendar")), class = "naboer_error_time")
})

test_that("temporal blocks are half open and retain whole blocks in V folds", {
  x <- tibble::tibble(time = c(4, 0, 1, 2, 3, 5, 6, 7))
  f <- temporal_folds(x, time_span(2, "seconds"), index = time, origin = 0)
  expect_equal(assessment_indices(f), list(c(1L, 6L), c(2L, 3L), c(4L, 5L), c(7L, 8L)))
  f2 <- temporal_folds(x, time_span(2, "seconds"), index = time, origin = 0, v = 2, seed = 99)
  sets <- assessment_indices(f2)
  expect_true(all(vapply(assessment_indices(f), function(b) any(vapply(sets, function(s) all(b %in% s), logical(1))), logical(1))))
  expect_error(temporal_folds(x, time_span(0, "days"), index = time), "positive")
  expect_error(temporal_folds(x, time_span(2, "observations"), index = time, within = "key", origin = 1), "do not supply")
})

test_that("count blocks and irregular elapsed blocks differ", {
  x <- tibble::tibble(time = c(0, 1, 100, 101, 102, 200))
  counts <- temporal_folds(x, time_span(2, "observations"), index = time, within = "key")
  expect_equal(assessment_indices(counts), list(1:2, 3:4, 5:6))
  elapsed <- temporal_folds(x, time_span(50, "seconds"), index = time)
  expect_equal(assessment_indices(elapsed), list(1:2, 3:5, 6L))
})

test_that("forecasts omit gaps and all future observations", {
  x <- tibble::tibble(time = c(6, 0, 4, 1, 5, 2, 3, 7))
  f <- forecast_folds(x,
    initial = time_span(2), horizon = time_span(2),
    step = time_span(2), gap = time_span(1), index = time
  )
  expect_equal(f$spec$origins, c(2, 4))
  expect_equal(assessment_indices(f)[[1]], c(3L, 5L))
  expect_equal(analysis_indices(f)[[1]], c(2L, 4L, 6L))
  expect_equal(omission_indices(f)[[1]], c(1L, 3L, 5L, 7L, 8L))
  sliding <- forecast_folds(x,
    origins = 4, horizon = time_span(2),
    window = "sliding", lookback = time_span(2), index = time
  )
  expect_equal(analysis_indices(sliding)[[1]], c(3L, 7L))
  expect_equal(assessment_indices(sliding)[[1]], c(1L, 5L))
  expect_s3_class(sliding$omission, "naboer_temporal_complement")
})

test_that("forecast origins span unequal series together and respect domains", {
  x <- tibble::tibble(
    network = c(rep("a", 8), rep("b", 4)),
    series = c(rep("x", 4), rep("y", 4), rep("z", 4)),
    time = c(0, 1, 3, 5, 0, 2, 4, 6, 0, 2, 4, 6)
  )
  f <- forecast_folds(x,
    horizon = time_span(3), origins = 2,
    index = time, key = series, .by = network
  )
  expect_equal(assessment_indices(f), list(c(3L, 4L, 7L), 11L))
  expect_equal(omission_indices(f), list(c(3L, 4L, 7L, 8L), c(11L, 12L)))
  expect_true(all(9:12 %in% analysis_indices(f)[[1]]))
  expect_true(all(1:8 %in% analysis_indices(f)[[2]]))
})

test_that("count forecasting and explicit-origin errors are clear", {
  x <- tibble::tibble(time = c(0, 1, 20, 50, 70, 90))
  f <- forecast_folds(x,
    initial = time_span(2, "observations"),
    horizon = time_span(2, "observations"), index = time
  )
  expect_equal(assessment_indices(f), list(3:4, 5:6))
  expect_equal(analysis_indices(f), list(1:2, 1:4))
  expect_error(forecast_folds(x, initial = time_span(1000), horizon = time_span(2), index = time), "No forecast")
  expect_error(forecast_folds(x, origins = -1, horizon = time_span(2), index = time), "nonempty")
  expect_error(forecast_folds(x, origins = 90, horizon = time_span(2), index = time), "nonempty")
  expect_error(forecast_folds(x, origins = 1, initial = time_span(1), horizon = time_span(2), index = time), "either")
})

test_that("tsibble metadata supplies temporal index and keys", {
  skip_if_not_installed("tsibble")
  x <- tsibble::as_tsibble(
    tibble::tibble(
      site = rep(c("a", "b"), each = 4),
      date = rep(as.Date("2024-01-01") + 0:3, 2)
    ),
    key = site, index = date
  )
  expect_error(temporal_folds(x, time_span(2, "days")), "Multiple series")
  f <- temporal_folds(x, time_span(2, "days"), within = "key")
  expect_equal(length(assessment_indices(f)), 4L)
  expect_equal(f$ctx$key, "site")
  expect_equal(f$ctx$index, "date")
})

test_that("period indexes use period-start dates and preserve source metadata", {
  skip_if_not_installed("tsibble")
  skip_if_not_installed("clock")
  starts <- as.Date(c("2023-12-01", "2024-01-01", "2024-02-01", "2024-03-01", "2024-04-01", "2024-05-01"))
  for (period in c("yearweek", "yearmonth", "yearquarter")) {
    constructor <- get(period, asNamespace("tsibble"))
    dates <- if (period == "yearquarter") seq(as.Date("2023-01-01"), by = "3 months", length.out = 6) else starts
    index <- constructor(dates)
    x <- tsibble::as_tsibble(tibble::tibble(time = index), index = time)
    ctx <- prepare_data(x)
    expect_s3_class(temporal_check_index(ctx), "Date")
    span <- if (period == "yearquarter") time_span(3, "months", "calendar") else time_span(1, "months", "calendar")
    folded <- temporal_folds(x, span)
    reference <- temporal_folds(tibble::tibble(time = as.Date(index)), span, index = time)
    expect_equal(assessment_indices(folded), assessment_indices(reference))
    expect_equal(fold_data(folded)$time, index)
    forecast <- forecast_folds(x, horizon = span, origins = index[3])
    expect_equal(assessment_indices(forecast)[[1]], 4L)
  }
})

test_that("numeric indices use seconds and calendar operations reject them", {
  x <- tibble::tibble(time = c(0, 60, 120, 3600))
  ctx <- prepare_data(x, index = rlang::quo(time))
  z <- temporal_membership(ctx, temporal_window(time_span(1, "minutes")), "domain")
  expect_equal(index_get(z, 2), 1:3)
  expect_error(temporal_blocks(ctx, time_span(1, "months", "calendar")), "Calendar")
  expect_error(time_span(1, "fortnights"), "not supported")
  dates <- tibble::tibble(time = as.Date("2024-01-01") + 0:3)
  expect_error(temporal_folds(dates, time_span(1, "days"), index = time, origin = 0), "same temporal index class")
})

test_that("calendar forecast origins avoid month-end drift", {
  skip_if_not_installed("clock")
  x <- tibble::tibble(time = as.Date(c("2024-01-31", "2024-02-29", "2024-03-29", "2024-03-31", "2024-04-29", "2024-04-30", "2024-05-31")))
  span <- time_span(1, "months", "calendar")
  f <- forecast_folds(x, origins = as.Date("2024-02-29"), horizon = span, index = time)
  expect_equal(analysis_indices(f)[[1]], 1:2)
  expect_equal(assessment_indices(f)[[1]], 3L)
  generated <- forecast_folds(x, initial = span, horizon = span, index = time)
  expect_equal(generated$spec$origins, as.Date(c("2024-02-29", "2024-03-31")))
})

test_that("temporal specification errors are actionable", {
  for (bad in list(NA_real_, Inf, numeric(), c(1, 2), "1")) {
    expect_error(time_span(bad), class = "naboer_error")
  }
  expect_error(time_span(1, NA_character_), "unit")
  expect_error(time_span(1, c("days", "weeks")), "unit")
  expect_error(time_span(1, "days", nonexistent = "guess"), "resolution")
  expect_error(temporal_window(time_span(1, "observations"), time_span(1)), "both")
  x <- tibble::tibble(time = c(1, 2, Inf))
  ctx <- prepare_data(x, index = rlang::quo(time))
  expect_error(temporal_membership(ctx, temporal_window(time_span(1)), "domain"), "finite")
  x <- tibble::tibble(time = letters[1:3])
  ctx <- prepare_data(x, index = rlang::quo(time))
  expect_error(temporal_membership(ctx, temporal_window(time_span(1)), "domain"), "finite")
  ctx <- prepare_data(tibble::tibble(time = 1:3))
  expect_error(temporal_membership(ctx, temporal_window(time_span(1))), "exactly one")
  expect_error(temporal_membership(ctx, time_span(1)), "temporal_window")
  expect_equal(temporal_seconds(time_span(1, "weeks")), 604800)
  expect_equal(temporal_seconds(time_span(1, "hours")), 3600)
})

test_that("calendar block origins and boundaries work before the origin", {
  skip_if_not_installed("clock")
  x <- tibble::tibble(time = as.Date(c("2024-01-01", "2024-01-30", "2024-01-31", "2024-02-01", "2024-02-29", "2024-03-30", "2024-03-31")))
  f <- temporal_folds(x, time_span(1, "months", "calendar"), index = time, origin = as.Date("2024-01-31"))
  expect_equal(assessment_indices(f), list(1:2, 3:4, 5:6, 7L))
  expect_error(temporal_folds(x, time_span(1, "days"), index = time, origin = as.Date(NA)), "nonmissing")
  expect_error(temporal_folds(x, time_span(1, "days"), index = time, origin = x$time[1:2]), "one nonmissing")
  x <- tibble::tibble(time = as.POSIXct(c("2024-03-30 12:00:00", "2024-03-31 11:00:00", "2024-03-31 12:00:00", "2024-04-01 12:00:00"), tz = "Europe/Copenhagen"))
  f <- temporal_folds(x, time_span(1, "days", "calendar"), index = time)
  expect_equal(assessment_indices(f), list(1:2, 3L, 4L))
  weeks <- temporal_folds(tibble::tibble(time = as.Date("2024-01-01") + 0:20), time_span(1, "weeks", "calendar"), index = time)
  expect_equal(assessment_indices(weeks), list(1:7, 8:14, 15:21))
  years <- temporal_folds(tibble::tibble(time = as.Date(c("2024-02-29", "2025-02-27", "2025-02-28", "2026-02-28"))), time_span(1, "years", "calendar"), index = time)
  expect_equal(assessment_indices(years), list(1:2, 3L, 4L))
})

test_that("forecast argument validation avoids ambiguous or leaking designs", {
  x <- tibble::tibble(time = 1:10)
  expect_error(forecast_folds(x, horizon = time_span(2), index = time), "initial")
  expect_error(forecast_folds(x, horizon = time_span(2), origins = 3, index = time, lookback = time_span(2)), "only used")
  expect_error(forecast_folds(x, horizon = time_span(2), origins = 3, index = time, window = "sliding"), "lookback")
  expect_error(forecast_folds(x, horizon = time_span(2), origins = 3, index = time, step = time_span(1, "observations")), "all count")
  expect_error(forecast_folds(x, horizon = time_span(2), origins = numeric(), index = time), "finite")
  expect_error(forecast_folds(x, horizon = time_span(2), origins = c(3, 3), index = time), "unique")
  expect_error(forecast_folds(x, horizon = time_span(2), origins = as.Date("2024-01-01"), index = time), "same temporal")
  expect_error(forecast_folds(x, horizon = time_span(20), initial = time_span(2), index = time), "No forecast")
  expect_error(forecast_folds(x, horizon = time_span(2, "observations"), origins = 2.5, index = time), "must occur")
  f <- forecast_folds(x,
    horizon = time_span(2, "observations"), origins = c(4, 6), index = time,
    window = "sliding", lookback = time_span(2, "observations")
  )
  expect_equal(analysis_indices(f), list(3:4, 5:6))
  expect_equal(assessment_indices(f), list(5:6, 7:8))
  expect_error(forecast_folds(x, horizon = time_span(2, "observations"), initial = time_span(100, "observations"), index = time), "No forecast")
})

test_that("POSIXlt indices retain calendar timezone semantics", {
  skip_if_not_installed("clock")
  original <- as.POSIXlt(as.POSIXct(c("2024-03-30 12:00:00", "2024-03-31 12:00:00", "2024-04-01 12:00:00"), tz = "Europe/Copenhagen"))
  ctx <- prepare_data(tibble::tibble(time = original), index = rlang::quo(time))
  expect_s3_class(temporal_check_index(ctx), "POSIXct")
  blocks <- temporal_blocks(ctx, time_span(1, "days", "calendar"))
  expect_equal(blocks, 1:3)
})

test_that("zero calendar shifts preserve already resolved repeated clock times", {
  skip_if_not_installed("clock")
  instants <- as.POSIXct(c("2024-10-27 00:30:00", "2024-10-27 01:30:00", "2024-10-28 01:30:00"), tz = "UTC")
  attr(instants, "tzone") <- "Europe/Copenhagen"
  for (resolution in c("error", "earliest", "latest")) {
    zero <- time_span(0, "days", "calendar", ambiguous = resolution)
    expect_identical(temporal_shift(instants, zero), instants)
    folded <- neighbourhoods(tibble::tibble(time = instants),
      temporal = temporal_window(zero), index = time
    )
    expect_equal(omission_indices(folded), as.list(1:3))
  }
  shifts <- temporal_shift(rep(instants[2], 2), time_span(1, "days", "calendar"), c(0, 1))
  expect_equal(as.numeric(shifts), as.numeric(instants[c(2, 3)]))
  # A known origin in the second repeated hour remains a valid cutoff when
  # the default zero gap accompanies a positive calendar horizon.
  forecast <- forecast_folds(tibble::tibble(time = instants),
    origins = instants[2], horizon = time_span(1, "days", "calendar"), index = time
  )
  expect_equal(analysis_indices(forecast), list(1:2))
  expect_equal(assessment_indices(forecast), list(3L))
})

test_that("forecasting never leaves another dependent series' future available", {
  data <- tibble::tibble(
    site = rep(c("a", "b"), each = 4),
    time = c(0, 1, 3, 5, 0, 2, 4, 6)
  )
  for (origins in list(NULL, 2)) {
    if (is.null(origins)) {
      expect_error(
        forecast_folds(data,
          initial = time_span(1),
          horizon = time_span(2), key = site, index = time, within = "key"
        ),
        "one series per independent domain",
        class = "naboer_error"
      )
    } else {
      expect_error(
        forecast_folds(data,
          origins = origins,
          horizon = time_span(2), key = site, index = time, within = "key"
        ),
        "one series per independent domain",
        class = "naboer_error"
      )
    }
  }
  expect_error(
    forecast_folds(data,
      initial = time_span(2, "observations"),
      horizon = time_span(1, "observations"), key = site, index = time
    ),
    "one series per independent domain"
  )
  expect_error(
    forecast_folds(data,
      origins = 2, horizon = time_span(2),
      key = site, index = time, within = "k"
    ),
    "one series per independent domain"
  )
  shared <- forecast_folds(data,
    origins = 2, horizon = time_span(2),
    key = site, index = time
  )
  expect_equal(analysis_indices(shared), list(c(1L, 2L, 5L, 6L)))
  expect_equal(omission_indices(shared), list(c(3L, 4L, 7L, 8L)))
  expect_equal(assessment_indices(shared), list(c(3L, 7L)))
  # Separating series through .by is an explicit declaration of independence.
  independent <- forecast_folds(data,
    initial = time_span(2, "observations"),
    horizon = time_span(1, "observations"), key = site, index = time, .by = site
  )
  expect_equal(analysis_indices(independent, 1)[[1]], c(1L, 2L, 5L, 6L, 7L, 8L))
  expect_equal(omission_indices(independent, 1)[[1]], 3:4)
})

test_that("temporal origins reject coercions and unrepresentable spans", {
  data <- tibble::tibble(time = 1:5)
  expect_error(temporal_folds(data, time_span(1), index = time, origin = "1"),
    "same temporal index class",
    class = "naboer_error"
  )
  expect_error(forecast_folds(data, origins = "2", horizon = time_span(1), index = time),
    "same temporal index class",
    class = "naboer_error"
  )
  expect_error(forecast_folds(data, origins = factor("2"), horizon = time_span(1), index = time),
    "same temporal index class",
    class = "naboer_error"
  )
  expect_error(neighbourhoods(data, temporal = temporal_window(time_span(1e308, "weeks")), index = time),
    "too large",
    class = "naboer_error"
  )
  matrix_data <- tibble::tibble(time = matrix(1:10, ncol = 2))
  expect_error(neighbourhoods(matrix_data, temporal = temporal_window(time_span(1)), index = time),
    "finite numeric",
    class = "naboer_error"
  )
})

test_that("forecast origins must advance at the index's numerical resolution", {
  data <- tibble::tibble(time = c(1e9, 1e9 + 1, 1e9 + 2))
  expect_error(forecast_folds(data,
    initial = time_span(1), horizon = time_span(0.1),
    step = time_span(1e-20), index = time
  ), "does not advance", class = "naboer_error_time")
})

test_that("calendar origins use the data index timezone rather than their display zone", {
  skip_if_not_installed("clock")
  dates <- as.POSIXct(c(
    "2024-03-30 12:00:00", "2024-03-31 12:00:00",
    "2024-03-31 12:30:00", "2024-04-01 12:00:00"
  ), tz = "Europe/Copenhagen")
  utc_origin <- dates[1]
  attr(utc_origin, "tzone") <- "UTC"
  data <- tibble::tibble(time = dates)
  blocks <- temporal_folds(data, time_span(1, "days", "calendar"),
    index = time, origin = utc_origin
  )
  expect_equal(assessment_indices(blocks), list(1L, 2:3, 4L))
  forecast <- forecast_folds(data,
    origins = utc_origin,
    horizon = time_span(1, "days", "calendar"), index = time
  )
  expect_equal(assessment_indices(forecast), list(2L))
  expect_equal(omission_indices(forecast), list(2:4))
  expect_equal(attr(forecast$spec$origins, "tzone"), "Europe/Copenhagen")
})

test_that("equivalent calendar origin steps retain the same leap-day anchor", {
  skip_if_not_installed("clock")
  data <- tibble::tibble(time = as.Date(c(
    "2024-02-29", "2025-02-28",
    "2026-02-28", "2027-02-28", "2028-02-29", "2029-02-28"
  )))
  years <- forecast_folds(data,
    initial = time_span(1, "years", "calendar"),
    step = time_span(1, "years", "calendar"), horizon = time_span(1, "years", "calendar"), index = time
  )
  months <- forecast_folds(data,
    initial = time_span(12, "months", "calendar"),
    step = time_span(1, "years", "calendar"), horizon = time_span(1, "years", "calendar"), index = time
  )
  expect_equal(months$spec$origins, years$spec$origins)
  # The 2027 origin has no assessment before its 2028-02-28 endpoint.
  expect_equal(months$spec$origins, as.Date(c("2025-02-28", "2026-02-28", "2028-02-29")))
})
