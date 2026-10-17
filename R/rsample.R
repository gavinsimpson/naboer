#' Convert between naboer folds and rsample resamples
#'
#' Analysis sets are the complement of omission sets. Assessment membership is
#' preserved independently, including exclusion-only observations in buffered
#' designs. Bootstrap resamples with duplicated analysis rows are not supported.
#' Resamples must index individual observations. Convert nested cubbles to
#' their temporal face before creating resamples; spatial-row indices cannot
#' safely be reused after expanding nested temporal observations.
#' Exported manual rsets retain column metadata for series keys, index,
#' independent domains and observation identifiers. Import restores it unless
#' explicit tidy selections override it. External rsets can supply these
#' selections directly.
#'
#' @param x A `naboer_folds` object for `as_rset()`, or an rsample `rset` for
#'   `from_rset()`.
#' @param id Optional unique observation identifier column, using tidy selection.
#' @param key,index,.by Tidy selections of series keys, temporal index and
#'   independent domains. Omitted selections restore metadata saved by
#'   `as_rset()`, or use native tsibble/cubble metadata where available.
#' @param max_bytes Maximum estimated storage for analysis and assessment indices
#'   exported to rsample, in bytes. Defaults to option `naboer.max_bytes`, or
#'   1 GiB. The estimate excludes R object overhead. Use `Inf` to disable the
#'   memory guard after checking the design size.
#' @return `as_rset()` returns an rsample manual rset; `from_rset()` returns
#'   a `naboer_folds` object.
#' @export
as_rset <- function(x, max_bytes = getOption("naboer.max_bytes", 1024^3)) {
  need_package("rsample")
  validate_folds(x)
  if (!is.numeric(max_bytes) || length(max_bytes) != 1L || is.na(max_bytes) ||
    max_bytes < 0) {
    abort_naboer("`max_bytes` must be nonnegative or Inf.")
  }
  if (x$ctx$n > .Machine$integer.max) {
    abort_naboer("Observation indices exceed rsample's integer index limit.", "overflow")
  }
  # A compact LOO design has tiny omissions but quadratically many analysis
  # indices. Check before either accessor materialises its membership lists.
  count <- sum(x$ctx$n - index_sizes(x$omission)) + sum(index_sizes(x$assessment))
  bytes <- 4 * count
  if (bytes > max_bytes) {
    abort_naboer(sprintf(paste0(
      "Exporting rsample indices requires approximately %.1f MiB before R ",
      "object overhead, exceeding `max_bytes`. Reduce the design or explicitly ",
      "raise the limit."
    ), bytes / 1024^2), "memory")
  }
  data <- fold_data(x)
  analysis <- analysis_indices(x)
  assessment <- assessment_indices(x)
  splits <- lapply(seq_along(assessment), function(i) {
    rsample::make_splits(list(
      analysis = as.integer(analysis[[i]]),
      assessment = as.integer(assessment[[i]])
    ), data)
  })
  out <- rsample::manual_rset(splits, x$labels)
  attr(out, "naboer_columns") <- x$ctx[c("key", "index", "by", "id")]
  out
}

#' @rdname as_rset
#' @export
from_rset <- function(x, id = NULL, key = NULL, index = NULL, .by = NULL) {
  need_package("rsample")
  if (!inherits(x, "rset") || !nrow(x) || !is.list(x$splits)) {
    abort_naboer("`x` must be a nonempty rsample rset.")
  }
  if (!all(vapply(x$splits, inherits, logical(1), "rsplit"))) {
    abort_naboer("Every resample must be an rsample rsplit.")
  }
  data <- x$splits[[1L]]$data
  if (!all(vapply(x$splits, function(split) identical(split$data, data), logical(1)))) {
    abort_naboer("All resamples must refer to identical data in identical row order.")
  }
  selections <- list(
    key = rlang::enquo(key), index = rlang::enquo(index),
    by = rlang::enquo(.by), id = rlang::enquo(id)
  )
  columns <- attr(x, "naboer_columns", exact = TRUE)
  if (!is.null(columns)) {
    fields <- names(selections)
    if (!is.list(columns) || !setequal(names(columns), fields) ||
      anyDuplicated(names(columns)) || !all(vapply(columns, function(z) {
      is.character(z) && !anyNA(z) && !anyDuplicated(z) && all(z %in% names(data))
    }, logical(1)))) {
      abort_naboer("Stored resample column metadata is invalid or refers to missing columns.", "alignment")
    }
    for (field in fields) {
      if (rlang::quo_is_null(selections[[field]])) {
        selections[[field]] <- rlang::new_quosure(
          rlang::expr(tidyselect::all_of(!!columns[[field]])), rlang::base_env()
        )
      }
    }
  }
  ctx <- prepare_data(data,
    key = selections$key, index = selections$index,
    by = selections$by, id = selections$id
  )
  if (ctx$n != nrow(data) || (!is.null(ctx$source) &&
    !identical(ctx$source$.source_row, seq_len(nrow(data))))) {
    abort_naboer(paste0(
      "Resamples must index individual observations, but preparing the data ",
      "changes its rows. Use `cubble::face_temporal()` or `fold_data()` before ",
      "creating rsample splits."
    ), "alignment")
  }
  assessment <- omission <- vector("list", nrow(x))
  all_rows <- seq_len(ctx$n)
  for (i in seq_len(nrow(x))) {
    analysis <- as.integer(x$splits[[i]], data = "analysis")
    assess <- as.integer(x$splits[[i]], data = "assessment")
    if (anyDuplicated(analysis)) {
      abort_naboer("Duplicated analysis rows, including bootstrap samples, cannot be represented by `nei`.")
    }
    for (indices in list(analysis, assess)) {
      if (!is.numeric(indices) || !length(indices) || anyNA(indices) ||
        any(!is.finite(indices)) || any(indices != floor(indices)) ||
        any(indices < 1L | indices > ctx$n) ||
        anyDuplicated(indices)) {
        abort_naboer("Resample indices must be nonempty, unique observation indices.")
      }
    }
    if (length(intersect(analysis, assess))) {
      abort_naboer("Resample analysis and assessment rows must not overlap.")
    }
    assessment[[i]] <- sort.int(assess)
    omission[[i]] <- setdiff(all_rows, analysis)
  }
  out <- new_folds(ctx, list_index(assessment), list_index(omission),
    spec = list(strategy = "rsample")
  )
  id_cols <- setdiff(names(x), "splits")
  if (length(id_cols)) {
    out$labels <- make.unique(do.call(paste, c(as.list(x[id_cols]), sep = "/")))
  }
  out
}
