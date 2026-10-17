#' Describe a temporal span
#'
#' Elapsed spans use seconds (a day is 86400 seconds); calendar spans use local
#' calendar arithmetic. Calendar month and year arithmetic clamps invalid dates
#' to the final day of the month. Count spans count ordered observations.
#' Numeric time indices are interpreted in seconds. Tsibble year-week,
#' year-month and year-quarter indices use the first date of each period for
#' elapsed and calendar arithmetic; source index columns are preserved.
#'
#' @param n Nonnegative span length. Calendar and count spans require integers.
#' @param unit A time unit, or `"observations"` for a count span.
#' @param type One of `"elapsed"`, `"calendar"`, or `"count"` (also accepted as
#'   `"observations"`). The default for `unit = "observations"` is `"count"`.
#' @param ambiguous,nonexistent Resolution passed to [clock::add_days()] for
#'   ambiguous or nonexistent local times; both default to `"error"`.
#' @return A temporal span specification.
#' @export
time_span <- function(n, unit = "seconds", type = "elapsed",
                      ambiguous = "error", nonexistent = "error") {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) || n < 0) {
    abort_naboer("`n` must be one finite nonnegative number.")
  }
  if (!is.character(unit) || length(unit) != 1L || is.na(unit)) {
    abort_naboer("`unit` must be one time-unit string.")
  }
  aliases <- c(
    second = "seconds", minute = "minutes", hour = "hours",
    day = "days", week = "weeks", month = "months", year = "years",
    observation = "observations"
  )
  if (unit %in% names(aliases)) unit <- unname(aliases[[unit]])
  if (missing(type) && unit == "observations") type <- "count"
  if (identical(type, "observations")) type <- "count"
  type <- match.arg(type, c("elapsed", "calendar", "count"))
  allowed <- switch(type,
    elapsed = c("seconds", "minutes", "hours", "days", "weeks"),
    calendar = c("days", "weeks", "months", "years"),
    count = "observations"
  )
  if (!unit %in% allowed) {
    abort_naboer(paste0("Unit `", unit, "` is not supported for `", type, "` spans."))
  }
  if (type != "elapsed" && n != floor(n)) {
    abort_naboer("Calendar and count spans require a whole-number `n`.")
  }
  if (!is.character(ambiguous) || length(ambiguous) != 1L ||
    !ambiguous %in% c("error", "earliest", "latest")) {
    abort_naboer("`ambiguous` must be 'error', 'earliest', or 'latest'.")
  }
  if (!is.character(nonexistent) || length(nonexistent) != 1L ||
    !nonexistent %in% c("error", "roll-forward", "roll-backward", "shift-forward", "shift-backward")) {
    abort_naboer("Unsupported `nonexistent` time resolution.")
  }
  structure(list(
    n = n, unit = unit, type = type, ambiguous = ambiguous,
    nonexistent = nonexistent
  ), class = "naboer_time_span")
}

#' Describe an asymmetric temporal neighbourhood
#' @param before,after Temporal spans made with [time_span()].
#' @return A temporal window specification.
#' @export
temporal_window <- function(before, after = before) {
  temporal_check_span(before, "before")
  temporal_check_span(after, "after")
  if (xor(before$type == "count", after$type == "count")) {
    abort_naboer("`before` and `after` must both count observations, or both measure time.")
  }
  structure(list(before = before, after = after), class = "naboer_temporal_window")
}

temporal_check_span <- function(x, arg, positive = FALSE) {
  if (!inherits(x, "naboer_time_span")) {
    abort_naboer(paste0("`", arg, "` must be constructed with `time_span()`."))
  }
  if (positive && x$n <= 0) abort_naboer(paste0("`", arg, "` must be positive."))
  invisible(x)
}

temporal_check_index <- function(ctx) {
  if (length(ctx$index) != 1L) abort_naboer("Supply exactly one temporal `index` column.")
  x <- temporal_normalize_index(ctx$data[[ctx$index]])
  if (!(inherits(x, "Date") || inherits(x, "POSIXt") || is.numeric(x)) ||
    !is.null(dim(x)) || length(x) != ctx$n ||
    anyNA(x) || any(!is.finite(as.numeric(x)))) {
    abort_naboer("The temporal index must contain finite numeric, Date, or POSIXct values without missing values.")
  }
  x
}

temporal_normalize_index <- function(x) {
  if (inherits(x, c("yearweek", "yearmonth", "yearquarter"))) {
    return(as.Date(x))
  }
  if (inherits(x, "POSIXlt")) {
    return(as.POSIXct(x))
  }
  x
}

temporal_domain_ids <- function(ctx, within = NULL, count = FALSE) {
  nkeys <- if (length(ctx$key)) length(unique(group_ids(ctx$data, c(ctx$by, ctx$key)))) else 1L
  if (is.null(within)) {
    if (nkeys > length(unique(ctx$domains))) {
      abort_naboer("Multiple series require explicit `within = 'key'` or `within = 'domain'`.")
    }
    within <- "domain"
  }
  within <- match.arg(within, c("key", "domain"))
  if (count && within != "key") {
    abort_naboer("Observation-count windows require `within = 'key'`.")
  }
  domain_ids(ctx, within)
}

temporal_seconds <- function(span) {
  seconds <- span$n * switch(span$unit,
    seconds = 1,
    minutes = 60,
    hours = 3600,
    days = 86400,
    weeks = 604800
  )
  if (!is.finite(seconds)) {
    abort_naboer("The elapsed span is too large to represent in seconds.")
  }
  seconds
}

temporal_numeric <- function(x) {
  if (inherits(x, "Date")) as.numeric(x) * 86400 else as.numeric(x)
}

temporal_shift <- function(x, span, multiplier = 1) {
  if (span$type == "elapsed") {
    seconds <- temporal_seconds(span) * multiplier
    if (inherits(x, "Date")) {
      return(x + seconds / 86400)
    }
    return(x + seconds)
  }
  if (span$type == "count") abort_naboer("Count spans require an ordered observation index.")
  if (!inherits(x, c("Date", "POSIXt"))) {
    abort_naboer("Calendar spans require a Date or POSIXct index.")
  }
  n <- span$n * multiplier
  recycled <- vctrs::vec_recycle_common(x = x, n = n)
  x <- recycled$x
  n <- recycled$n
  changed <- n != 0
  # The input already identifies a real instant, including which occurrence
  # of a repeated clock time it represents. A zero shift must preserve it.
  if (!any(changed)) {
    return(x)
  }
  need_package("clock")
  args <- list(x = x[changed], n = n[changed])
  if (span$unit %in% c("months", "years")) args$invalid <- "previous"
  if (inherits(x, "POSIXt")) {
    args$ambiguous <- span$ambiguous
    args$nonexistent <- span$nonexistent
  }
  fun <- switch(span$unit,
    days = clock::add_days,
    weeks = clock::add_weeks,
    months = clock::add_months,
    years = clock::add_years
  )
  result <- tryCatch(do.call(fun, args), error = function(e) {
    abort_naboer(paste0("Calendar arithmetic failed: ", conditionMessage(e)), "time")
  })
  x[changed] <- result
  x
}

temporal_compatible_origin <- function(x, origin) {
  (inherits(origin, "Date") || inherits(origin, "POSIXt") || is.numeric(origin)) &&
    is.null(dim(origin)) &&
    inherits(x, "Date") == inherits(origin, "Date") &&
    inherits(x, "POSIXt") == inherits(origin, "POSIXt")
}

temporal_origin_zone <- function(x, origin) {
  if (inherits(x, "POSIXt") && inherits(origin, "POSIXt")) {
    # Origins identify instants, while calendar arithmetic follows the data's
    # civil-time zone even if an origin was supplied in another display zone.
    attr(origin, "tzone") <- attr(x, "tzone")
  }
  origin
}

# Compact focal membership in the permutation sorted by domain and time.
temporal_membership <- function(ctx, window, within = NULL) {
  if (!inherits(window, "naboer_temporal_window")) {
    abort_naboer("`temporal` must be constructed with `temporal_window()`.")
  }
  x <- temporal_check_index(ctx)
  count <- window$before$type == "count"
  domains <- temporal_domain_ids(ctx, within, count)
  if (count) temporal_unique_times(ctx, domains)
  ord <- order(domains, x, seq_len(ctx$n))
  runs <- split(seq_along(ord), domains[ord])
  start <- end <- integer(ctx$n)
  for (positions in runs) {
    rows <- ord[positions]
    n <- length(rows)
    if (count) {
      lo <- pmax(1L, seq_len(n) - window$before$n)
      hi <- pmin(n, seq_len(n) + window$after$n)
    } else {
      values <- temporal_numeric(x[rows])
      lower <- temporal_numeric(temporal_shift(x[rows], window$before, -1))
      upper <- temporal_numeric(temporal_shift(x[rows], window$after, 1))
      lo <- findInterval(lower, values, left.open = TRUE) + 1L
      hi <- findInterval(upper, values)
    }
    start[rows] <- positions[1L] + lo - 1L
    end[rows] <- positions[1L] + hi - 1L
  }
  range_index(start, end, order = ord)
}

#' Construct contiguous temporal blocks
#' @param data A data frame, tsibble, or cubble.
#' @param size A positive [time_span()] specifying block width.
#' @param index,key,.by,id Columns identifying time, series, independent domains,
#'   and unique observations. Use unquoted column names or tidy selections.
#' @param origin A scalar index value defining block boundaries. By default use
#'   the earliest observation within each independent domain or series.
#' @param within Choose `"key"` or `"domain"` for multiple series. Count spans
#'   require `"key"` and count rows in stable time order.
#' @param v Optional number of folds to which whole blocks are allocated.
#' @param balance Balance observation counts when assigning blocks to folds.
#' @param seed Optional local random seed.
#' @return A `naboer_folds` object.
#' @export
temporal_folds <- function(data, size, index = NULL, key = NULL, .by = NULL,
                           id = NULL, origin = NULL, within = NULL, v = NULL,
                           balance = TRUE, seed = NULL) {
  ctx <- prepare_data(data,
    key = rlang::enquo(key), index = rlang::enquo(index),
    by = rlang::enquo(.by), id = rlang::enquo(id)
  )
  blocks <- temporal_blocks(ctx, size, origin, within)
  group_folds_from_blocks(ctx, blocks,
    v = v, balance = balance, seed = seed,
    spec = list(
      strategy = "temporal", size = size,
      origin = origin, within = within
    )
  )
}

temporal_blocks <- function(ctx, size, origin = NULL, within = NULL) {
  temporal_check_span(size, "size", positive = TRUE)
  x <- temporal_check_index(ctx)
  domains <- temporal_domain_ids(ctx, within, size$type == "count")
  if (size$type == "count") temporal_unique_times(ctx, domains)
  origin <- temporal_normalize_index(origin)
  if (!is.null(origin) && (length(origin) != 1L || anyNA(origin))) {
    abort_naboer("`origin` must be one nonmissing temporal index value.")
  }
  if (!is.null(origin) &&
    (!temporal_compatible_origin(x, origin) ||
      !is.finite(as.numeric(origin)))) {
    abort_naboer("`origin` must be finite and have the same temporal index class as the data.")
  }
  if (size$type == "count" && !is.null(origin)) {
    abort_naboer("Count blocks start at the first observation; do not supply `origin`.")
  }
  origin <- temporal_origin_zone(x, origin)
  block <- numeric(ctx$n)
  for (rows in split(seq_len(ctx$n), domains)) {
    if (size$type == "count") {
      ordered <- rows[order(x[rows], rows)]
      block[ordered] <- (seq_along(ordered) - 1L) %/% size$n
      next
    }
    anchor <- if (is.null(origin)) min(x[rows]) else origin
    if (size$type == "elapsed") {
      block[rows] <- floor((temporal_numeric(x[rows]) - temporal_numeric(anchor)) /
        temporal_seconds(size))
    } else {
      if (!inherits(anchor, class(x)[1L])) {
        abort_naboer("Calendar `origin` must have the same index class as the data.")
      }
      local_date <- function(y) {
        if (inherits(y, "POSIXt")) {
          as.Date(y, tz = if (is.null(attr(y, "tzone"))) "" else attr(y, "tzone"))
        } else {
          as.Date(y)
        }
      }
      if (size$unit %in% c("months", "years")) {
        xx <- as.POSIXlt(x[rows])
        aa <- as.POSIXlt(anchor)
        delta <- (xx$year - aa$year) * 12 + xx$mon - aa$mon
        divisor <- size$n * if (size$unit == "years") 12 else 1
        k <- floor(delta / divisor)
      } else {
        delta <- as.numeric(local_date(x[rows]) - local_date(anchor))
        divisor <- size$n * if (size$unit == "weeks") 7 else 1
        k <- floor(delta / divisor)
      }
      values <- temporal_numeric(x[rows])
      # Correct estimated positions against direct shifts from the anchor,
      # avoiding cumulative month-end drift and retaining half-open bins.
      repeat {
        lower <- temporal_numeric(temporal_shift(rep(anchor, length(k)), size, k))
        bad <- values < lower
        if (!any(bad)) break
        k[bad] <- k[bad] - 1
      }
      repeat {
        upper <- temporal_numeric(temporal_shift(rep(anchor, length(k)), size, k + 1))
        bad <- values >= upper
        if (!any(bad)) break
        k[bad] <- k[bad] + 1
      }
      block[rows] <- k
    }
  }
  group_ids(data.frame(domain = domains, block = block), c("domain", "block"))
}

#' Construct rolling-origin forecasting folds
#'
#' Assessment observations are strictly after the origin and gap and no later
#' than the horizon endpoint. Every future observation, including observations
#' beyond the horizon, is omitted. Sliding windows also omit old observations.
#' Other independent domains remain available for fitting.
#' Forecasts across several series in the same domain must use time spans and
#' `within = "domain"`, so that every series shares the cutoff. Count forecasts
#' require one series per independent domain; `within = "key"` is rejected when
#' it would leave another series' future observations available in that domain.
#'
#' @inheritParams temporal_folds
#' @param initial Initial training span, measured from the first observation.
#'   Required unless `origins` is supplied.
#' @param horizon Positive span assessed after each gap.
#' @param step Positive distance between generated origins; defaults to horizon.
#' @param gap Nonnegative span between training and assessment; defaults to zero
#'   in the same units and type as `horizon`.
#' @param window `"expanding"` or `"sliding"` training window.
#' @param lookback Positive training span for sliding windows. Training includes
#'   observations strictly after `origin - lookback` and through the origin.
#' @param origins Optional explicit index values. Generated origins require a
#'   complete horizon within the observed range; explicit origins require at
#'   least one assessment and one training observation in each domain.
#' @return A `naboer_folds` object. Generated empty assessment windows are skipped.
#' @export
forecast_folds <- function(data, initial = NULL, horizon, step = horizon,
                           gap = NULL, window = c("expanding", "sliding"),
                           lookback = NULL, origins = NULL, index = NULL,
                           key = NULL, .by = NULL, id = NULL, within = NULL) {
  window <- match.arg(window)
  temporal_check_span(horizon, "horizon", positive = TRUE)
  temporal_check_span(step, "step", positive = TRUE)
  if (is.null(gap)) {
    gap <- horizon
    gap$n <- 0
  }
  temporal_check_span(gap, "gap")
  if (is.null(origins)) {
    temporal_check_span(initial, "initial", positive = TRUE)
  } else if (!is.null(initial)) abort_naboer("Supply either `initial` or `origins`, not both.")
  if (window == "sliding") {
    temporal_check_span(lookback, "lookback", positive = TRUE)
  } else if (!is.null(lookback)) abort_naboer("`lookback` is only used for sliding windows.")
  spans <- Filter(Negate(is.null), list(initial, horizon, step, gap, lookback))
  count <- horizon$type == "count"
  if (any(vapply(spans, function(s) (s$type == "count") != count, logical(1)))) {
    abort_naboer("Forecast spans must all count observations or all measure time.")
  }
  ctx <- prepare_data(data,
    key = rlang::enquo(key), index = rlang::enquo(index),
    by = rlang::enquo(.by), id = rlang::enquo(id)
  )
  x <- temporal_check_index(ctx)
  # Forecasting defaults to shared origins across the independent domain.
  if (is.null(within)) within <- if (count) "key" else "domain"
  within <- match.arg(within, c("key", "domain"))
  if (identical(within, "key") && length(ctx$key) &&
    length(unique(group_ids(ctx$data, c(ctx$by, ctx$key)))) >
      length(unique(ctx$domains))) {
    abort_naboer(paste0(
      "`within = 'key'` forecasting requires one series per independent domain. ",
      "Use elapsed/calendar spans with `within = 'domain'` for shared cutoffs; ",
      "use `.by` to separate series only when they are independent."
    ))
  }
  domains <- temporal_domain_ids(ctx, within, count)
  if (count) temporal_unique_times(ctx, domains)
  if (!is.null(origins)) {
    origins <- temporal_normalize_index(origins)
    if (!temporal_compatible_origin(x, origins)) {
      abort_naboer("`origins` must have the same temporal index class as the data.")
    }
    if (!length(origins) || anyNA(origins) || any(!is.finite(as.numeric(origins)))) {
      abort_naboer("`origins` must contain finite nonmissing index values.")
    }
    if (anyDuplicated(origins)) abort_naboer("`origins` must be unique.")
    origins <- temporal_origin_zone(x, origins)
    origins <- sort(origins)
  }
  ord <- order(domains, x, seq_len(ctx$n))
  astart <- aend <- tstart <- tend <- dstart <- dend <- integer()
  origin_values <- list()
  fold_domains <- integer()
  for (positions in split(seq_along(ord), domains[ord])) {
    rows <- ord[positions]
    xx <- x[rows]
    values <- temporal_numeric(xx)
    n <- length(rows)
    explicit <- !is.null(origins)
    if (count) {
      if (explicit) {
        offsets <- match(origins, xx)
        if (anyNA(offsets)) abort_naboer("Count-based `origins` must occur in every selected series.")
      } else {
        last <- n - gap$n - horizon$n
        offsets <- if (initial$n <= last) seq(initial$n, last, by = step$n) else numeric()
      }
      origin <- xx[offsets]
      lower <- pmin(n, offsets + gap$n)
      upper <- pmin(n, offsets + gap$n + horizon$n)
      train_lo <- if (window == "expanding") rep.int(1, length(offsets)) else pmax(1, offsets - lookback$n + 1)
      train_hi <- offsets
    } else {
      if (explicit) {
        origin <- origins
      } else {
        first <- temporal_shift(xx[1L], initial)
        if (temporal_numeric(first) > values[n]) next
        # Shift directly from an anchor. Equivalent calendar units share an
        # anchor (years/months; weeks/days), avoiding month-end drift.
        out <- list()
        k <- 0L
        previous <- -Inf
        repeat {
          month_units <- c("months", "years")
          day_units <- c("days", "weeks")
          if (initial$type == "calendar" && step$type == "calendar" &&
            ((initial$unit %in% month_units && step$unit %in% month_units) ||
              (initial$unit %in% day_units && step$unit %in% day_units))) {
            distance <- step
            distance$unit <- if (step$unit %in% month_units) "months" else "days"
            multiplier_for <- function(unit) {
              switch(unit,
                years = 12,
                weeks = 7,
                1
              )
            }
            distance$n <- initial$n * multiplier_for(initial$unit) +
              k * step$n * multiplier_for(step$unit)
            candidate <- temporal_shift(xx[1L], distance)
          } else {
            candidate <- temporal_shift(first, step, k)
          }
          candidate_numeric <- temporal_numeric(candidate)
          if (candidate_numeric <= previous) {
            abort_naboer(paste0(
              "The forecast `step` does not advance to a distinct instant. ",
              "Increase its width or supply explicit distinct `origins`."
            ), "time")
          }
          previous <- candidate_numeric
          endpoint <- temporal_shift(temporal_shift(candidate, gap), horizon)
          if (temporal_numeric(endpoint) > values[n]) break
          out[[length(out) + 1L]] <- candidate
          k <- k + 1L
          if (k > .Machine$integer.max - 1L) abort_naboer("Too many forecast origins.")
        }
        origin <- if (length(out)) do.call(c, out) else xx[FALSE]
      }
      lower <- findInterval(temporal_numeric(temporal_shift(origin, gap)), values)
      upper <- findInterval(temporal_numeric(temporal_shift(temporal_shift(origin, gap), horizon)), values)
      train_hi <- findInterval(temporal_numeric(origin), values)
      train_lo <- if (window == "expanding") {
        rep.int(1, length(origin))
      } else {
        findInterval(temporal_numeric(temporal_shift(origin, lookback, -1)), values) + 1L
      }
    }
    valid <- upper > lower & train_hi >= train_lo
    if (explicit && any(!valid)) {
      abort_naboer("Each explicit origin must have nonempty training and assessment sets in every domain.")
    }
    if (!any(valid)) next
    offset <- positions[1L] - 1L
    astart <- c(astart, offset + lower[valid] + 1L)
    aend <- c(aend, offset + upper[valid])
    tstart <- c(tstart, offset + train_lo[valid])
    tend <- c(tend, offset + train_hi[valid])
    dstart <- c(dstart, rep.int(positions[1L], sum(valid)))
    dend <- c(dend, rep.int(positions[length(positions)], sum(valid)))
    origin_values[[length(origin_values) + 1L]] <- origin[valid]
    fold_domains <- c(fold_domains, rep.int(domains[rows[1L]], sum(valid)))
  }
  if (!length(astart)) abort_naboer("No forecast folds have nonempty training and assessment sets; reduce `initial`, `gap`, or `horizon`.")
  new_folds(ctx, range_index(astart, aend, order = ord),
    temporal_complement_index(tstart, tend, dstart, dend, ord),
    spec = list(
      strategy = "forecast", window = window, initial = initial,
      horizon = horizon, gap = gap, step = step, lookback = lookback,
      origins = do.call(c, origin_values), domains = fold_domains,
      within = within
    )
  )
}

temporal_complement_index <- function(start, end, domain_start, domain_end, order) {
  structure(
    list(
      start = as.integer(start), end = as.integer(end),
      domain_start = as.integer(domain_start), domain_end = as.integer(domain_end),
      order = as.integer(order)
    ),
    class = c("naboer_temporal_complement", "naboer_index")
  )
}

temporal_unique_times <- function(ctx, domains) {
  key_time <- data.frame(domain = domains, time = ctx$data[[ctx$index]])
  if (anyDuplicated(key_time)) {
    abort_naboer("Observation-count spans require unique index values within each series; resolve tied times first.")
  }
}
