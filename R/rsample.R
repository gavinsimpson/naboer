#' Convert between naboer folds and rsample resamples
#'
#' Analysis sets are the complement of omission sets. Assessment membership is
#' preserved independently, including exclusion-only observations in buffered
#' designs. Bootstrap resamples with duplicated analysis rows are not supported.
#'
#' @param x A `naboer_folds` object for `as_rset()`, or an rsample `rset` for
#'   `from_rset()`.
#' @param id Optional unique observation identifier column, using tidy selection.
#' @return `as_rset()` returns an rsample manual rset; `from_rset()` returns
#'   a `naboer_folds` object.
#' @export
as_rset <- function(x) {
  need_package("rsample")
  validate_folds(x)
  data <- fold_data(x)
  analysis <- analysis_indices(x)
  assessment <- assessment_indices(x)
  splits <- lapply(seq_along(assessment), function(i) {
    rsample::make_splits(list(
      analysis = as.integer(analysis[[i]]),
      assessment = as.integer(assessment[[i]])
    ), data)
  })
  rsample::manual_rset(splits, x$labels)
}

#' @rdname as_rset
#' @export
from_rset <- function(x, id = NULL) {
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
  ctx <- prepare_data(data, id = rlang::enquo(id))
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
