#' Export mgcv neighbourhood structures
#'
#' The result is a plain list of one-based indices, with cumulative endpoints.
#' No rescaling of repeated assessments or covariance settings is performed.
#' Neighbourhood sampling preserves the complete structure needed by bam.
#' @param x A `naboer_folds` object.
#' @param sample Optional bam sampling specification: one positive count, or a
#'   vector of neighbourhood numbers. A length-one value always means a count.
#' @param jackknife Optional logical gam covariance setting. `NULL` preserves
#'   mgcv's default by omitting the field.
#' @param max_bytes Maximum estimated packed output size in bytes; defaults to
#'   option `naboer.max_bytes`, or 1 GiB. Use `Inf` to disable the memory guard.
#' @param data Optional data used for alignment checking. Supply unchanged
#'   canonical data, or data containing the original unique identifier column.
#' @return A plain list with `a`, `ma`, `d`, `md`, and optional settings.
#' @export
as_nei <- function(x, sample = NULL, jackknife = NULL,
                   max_bytes = getOption("naboer.max_bytes", 1024^3), data = NULL) {
  validate_folds(x)
  if (!is.null(data)) {
    if (!is.data.frame(data) || nrow(data) != x$ctx$n) abort_naboer("Data rows do not match fold observations.", "alignment")
    aligned <- if (length(x$ctx$id)) {
      identical(data[[x$ctx$id]], x$ctx$data[[x$ctx$id]])
    } else {
      identical(tibble::as_tibble(data), x$ctx$data)
    }
    if (!aligned) abort_naboer("Data order or identifiers do not match the fold data. Rebuild folds after changing observations.", "alignment")
  }
  m <- index_nfold(x$assessment)
  if (!is.null(jackknife) && (!is.logical(jackknife) || length(jackknife) != 1L || is.na(jackknife))) abort_naboer("`jackknife` must be NULL, TRUE or FALSE.")
  if (!is.null(sample)) {
    if (!is.numeric(sample) || !length(sample) || anyNA(sample) || any(!is.finite(sample)) ||
      any(sample != floor(sample)) || any(sample < 1 | sample > m) ||
      (length(sample) > 1L && anyDuplicated(sample))) {
      abort_naboer("`sample` must be a positive fold count or unique valid fold indices.")
    }
    sample <- as.integer(sample)
  }
  if (!is.numeric(max_bytes) || length(max_bytes) != 1L || is.na(max_bytes) || max_bytes < 0) abort_naboer("`max_bytes` must be nonnegative or Inf.")
  lengths <- c(sum(index_sizes(x$assessment)), sum(index_sizes(x$omission)))
  if (any(lengths > .Machine$integer.max) || m > .Machine$integer.max) abort_naboer("The packed structure exceeds mgcv's integer index limit. Reduce the number or width of folds.", "overflow")
  bytes <- 4 * (sum(lengths) + 2 * m)
  if (bytes > max_bytes) abort_naboer(sprintf("Packed nei requires approximately %.1f MiB, exceeding `max_bytes`. Reduce the design or explicitly raise the limit.", bytes / 1024^2), "memory")
  a <- pack_index(x$omission)
  d <- pack_index(x$assessment)
  out <- list(a = a$values, ma = a$ends, d = d$values, md = d$ends)
  if (!is.null(sample)) out$sample <- sample
  if (!is.null(jackknife)) out$jackknife <- jackknife
  validate_nei(out, x$ctx$n)
  out
}

#' Validate a packed mgcv neighbourhood structure
#' @param x A list with `a`, `ma`, `d` and `md` components.
#' @param n Number of observations; required for checking upper bounds.
#' @return `x`, invisibly, or a structured error.
#' @export
validate_nei <- function(x, n) {
  n <- check_count(n, "n", 2L)
  required <- c("a", "ma", "d", "md")
  if (!is.list(x) || !all(required %in% names(x))) abort_naboer("Supply all of a, ma, d and md.", "structure")
  for (nm in required) {
    z <- x[[nm]]
    if (!is.numeric(z) || !length(z) || anyNA(z) || any(!is.finite(z)) ||
      any(z < 1 | z > .Machine$integer.max | z != floor(z))) {
      abort_naboer(paste0("Invalid `", nm, "` indices."), "indices")
    }
  }
  if (length(x$ma) != length(x$md) || any(diff(c(0, x$ma)) < 1) ||
    any(diff(c(0, x$md)) < 1) || tail(x$ma, 1) != length(x$a) ||
    tail(x$md, 1) != length(x$d) || any(x$a > n) || any(x$d > n)) {
    abort_naboer("Packed endpoints or observation bounds are invalid.", "indices")
  }
  for (i in seq_along(x$ma)) {
    a <- x$a[seq.int(if (i == 1L) 1L else x$ma[[i - 1L]] + 1L, x$ma[[i]])]
    d <- x$d[seq.int(if (i == 1L) 1L else x$md[[i - 1L]] + 1L, x$md[[i]])]
    if (anyDuplicated(a) || anyDuplicated(d) || !all(d %in% a) || length(a) >= n) abort_naboer("Invalid assessment/omission membership or empty analysis set.", "indices")
  }
  if (!is.null(x$jackknife) && (!is.logical(x$jackknife) || length(x$jackknife) != 1L || is.na(x$jackknife))) abort_naboer("Invalid jackknife setting.")
  if (!is.null(x$sample)) {
    s <- x$sample
    if (!is.numeric(s) || !length(s) || anyNA(s) || any(!is.finite(s)) || any(s != floor(s)) ||
      any(s < 1 | s > length(x$ma)) || (length(s) > 1L && anyDuplicated(s))) {
      abort_naboer("Invalid sample setting.")
    }
  }
  invisible(x)
}
