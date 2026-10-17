index_frequency <- function(x, n) {
  if (inherits(x, "naboer_ranges")) {
    active <- x$end >= x$start
    delta <- tabulate(x$start[active], n + 1L) - tabulate(x$end[active] + 1L, n + 1L)
    freq <- cumsum(delta)[seq_len(n)]
    if (!is.null(x$order)) freq <- freq[x$inverse]
    return(freq)
  }
  if (inherits(x, "naboer_shared")) {
    weights <- tabulate(x$lookup, index_nfold(x$members))
    out <- numeric(n)
    for (i in which(weights > 0L)) {
      rows <- index_get(x$members, i)
      out[rows] <- out[rows] + weights[[i]]
    }
    return(out)
  }
  out <- numeric(n)
  for (i in seq_len(index_nfold(x))) {
    rows <- index_get(x, i)
    out[rows] <- out[rows] + 1
  }
  out
}

#' Inspect folds as a tibble
#' @param x A `naboer_folds` object.
#' @param ... Unused; must be empty.
#' @param folds Optional integer fold selection, defaulting to all folds.
#' @param membership Include assessment and omission list-columns. Defaults to
#'   FALSE so large designs are not expanded accidentally.
#' @param .rows,.name_repair Reserved tibble arguments.
#' @return A tibble with fold identifiers and analysis/assessment/omission counts.
#'   `n_buffer` counts all omission-only observations, including unavailable
#'   future or old observations in forecasting designs, not just spatial or
#'   temporal dependence buffers.
#' @exportS3Method tibble::as_tibble
as_tibble.naboer_folds <- function(x, ..., folds = NULL, membership = FALSE,
                                   .rows = NULL, .name_repair = "check_unique") {
  rlang::check_dots_empty()
  selected <- fold_selection(x, folds)
  aa <- index_sizes(x$assessment)[selected]
  oo <- index_sizes(x$omission)[selected]
  out <- tibble::tibble(
    fold = x$labels[selected], n_assessment = aa,
    n_omission = oo, n_analysis = x$ctx$n - oo,
    n_buffer = oo - aa, omitted_fraction = oo / x$ctx$n
  )
  if (isTRUE(membership)) {
    out$assessment <- assessment_indices(x, selected)
    out$omission <- omission_indices(x, selected)
  }
  out
}

#' Summarise a cross-validation design
#' @param object A `naboer_folds` object.
#' @param ... Unused; must be empty.
#' @return A list with design counts, fold-size ranges and an observation-level
#'   assessment-frequency tibble. Summary never expands membership list-columns.
#'   `buffer_size` counts omission-only observations, including forecasting
#'   exclusions outside the assessment horizon or training window.
#' @export
summary.naboer_folds <- function(object, ...) {
  rlang::check_dots_empty()
  aa <- index_sizes(object$assessment)
  oo <- index_sizes(object$omission)
  freq <- index_frequency(object$assessment, object$ctx$n)
  list(
    strategy = object$spec$strategy, observations = object$ctx$n,
    folds = length(aa), assessment_size = range(aa), omission_size = range(oo),
    omitted_fraction = range(oo / object$ctx$n), buffer_size = range(oo - aa),
    unassessed = sum(freq == 0L), repeated_assessment = sum(freq > 1L),
    assessment_frequency = tibble::tibble(observation = seq_len(object$ctx$n), frequency = freq)
  )
}

#' @export
print.naboer_folds <- function(x, ...) {
  cat("<naboer_folds: ", x$spec$strategy, ">\n", sep = "")
  cat(
    format(x$ctx$n, big.mark = ","), "observations;",
    format(index_nfold(x$assessment), big.mark = ","), "folds\n"
  )
  invisible(x)
}

#' Plot selected cross-validation folds
#'
#' Spatial geometries are used when available, otherwise the temporal index or
#' observation row is plotted. Only requested folds are expanded.
#' The `omitted` role denotes omission-only observations; these can be dependence
#' buffers, unavailable future observations, or other exclusions.
#' @param object A `naboer_folds` object.
#' @param folds Integer fold selection; defaults to the first fold.
#' @param ... Unused; must be empty.
#' @return A ggplot object.
#' @exportS3Method ggplot2::autoplot
autoplot.naboer_folds <- function(object, folds = 1L, ...) {
  need_package("ggplot2")
  rlang::check_dots_empty()
  selected <- fold_selection(object, folds)
  if (!length(selected)) abort_naboer("Select at least one fold to plot.")
  geom <- names(object$ctx$data)[vapply(object$ctx$data, inherits, logical(1), "sfc")]
  pieces <- lapply(selected, function(i) {
    data <- object$ctx$data
    role <- rep("analysis", object$ctx$n)
    role[index_get(object$omission, i)] <- "omitted"
    role[index_get(object$assessment, i)] <- "assessment"
    data$.naboer_role <- factor(role, c("analysis", "omitted", "assessment"))
    data$.naboer_fold <- object$labels[[i]]
    data$.naboer_observation <- seq_len(object$ctx$n)
    if (length(object$ctx$index)) data$.naboer_time <- data[[object$ctx$index]]
    data
  })
  data <- dplyr::bind_rows(pieces)
  if (length(geom)) {
    need_package("sf")
    data <- sf::st_as_sf(data, sf_column_name = geom[[1L]])
    p <- ggplot2::ggplot(data) +
      ggplot2::geom_sf(ggplot2::aes(colour = .data$.naboer_role))
  } else {
    if (!length(object$ctx$index)) data$.naboer_time <- data$.naboer_observation
    p <- ggplot2::ggplot(data, ggplot2::aes(
      x = .data$.naboer_time,
      y = .data$.naboer_role, colour = .data$.naboer_role
    )) +
      ggplot2::geom_point()
  }
  p + ggplot2::facet_wrap(ggplot2::vars(.data$.naboer_fold)) +
    ggplot2::scale_colour_manual(values = c(analysis = "#737373", omitted = "#E69F00", assessment = "#0072B2"), drop = FALSE) +
    ggplot2::labs(colour = "Role", x = if (length(object$ctx$index)) object$ctx$index else "Observation", y = NULL) +
    ggplot2::theme_minimal()
}
