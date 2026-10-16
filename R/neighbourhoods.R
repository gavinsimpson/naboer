#' Construct focal-observation neighbourhoods
#'
#' Each observation is assessed once, omitting its dependence neighbourhood.
#' Multiple windows require explicit union or intersection. Series keys identify
#' observations; `.by` prevents dependence across independent domains.
#' @param data A data frame, tibble, sf object, tsibble or cubble.
#' @param spatial A [spatial_window()] specification.
#' @param temporal A [temporal_window()] specification.
#' @param phylogenetic A [phylogenetic_window()] specification.
#' @param combine Either `"intersection"` or `"union"`; required for multiple windows.
#' @param key,index,.by,id Tidy selections of series keys, temporal index,
#'   independent domains and observation identifiers.
#' @param within For temporal searches, `"key"` or `"domain"`. Required when
#'   purely temporal data contain multiple series. Defaults to `"domain"` for
#'   spatiotemporal searches.
#' @return A compact `naboer_folds` object.
#' @export
neighbourhoods <- function(data, spatial = NULL, temporal = NULL, phylogenetic = NULL,
                           combine = NULL, key = NULL, index = NULL, .by = NULL,
                           within = NULL, id = NULL) {
  ctx <- prepare_data(data, rlang::enquo(key), rlang::enquo(index), rlang::enquo(.by), rlang::enquo(id))
  count <- sum(!vapply(list(spatial, temporal, phylogenetic), is.null, logical(1)))
  if (!count) abort_naboer("Supply at least one neighbourhood window; use loo_folds() for leave-one-out.")
  if (count > 1L && is.null(combine)) abort_naboer("Choose `combine = \"intersection\"` or `\"union\"` for multiple windows.")
  if (!is.null(combine)) combine <- match.arg(combine, c("intersection", "union"))
  parts <- list()
  if (!is.null(spatial)) parts$spatial <- spatial_membership(ctx, spatial)
  if (!is.null(temporal)) {
    if (is.null(within) && !is.null(spatial)) within <- "domain"
    parts$temporal <- temporal_membership(ctx, temporal, within)
  }
  if (!is.null(phylogenetic)) parts$phylogenetic <- phylogenetic_membership(ctx, phylogenetic)
  omit <- if (length(parts) == 1L) {
    parts[[1L]]
  } else if (
    identical(names(parts), c("spatial", "temporal")) && identical(combine, "intersection") &&
      temporal$before$type != "count"
  ) {
    spatiotemporal_membership(ctx, parts$spatial, temporal, within)
  } else {
    combined_index(parts, combine)
  }
  new_folds(
    ctx, point_index(ctx$n), omit,
    list(
      strategy = "neighbourhood", spatial = spatial, temporal = temporal,
      phylogenetic = phylogenetic, combine = combine, within = within,
      spatial_engine = if (!is.null(spatial)) attr(parts$spatial, "spatial_engine")
    )
  )
}

#' Add dependence buffers around fold assessment observations
#'
#' Existing omissions are retained. Only assessment observations define buffer
#' centres; buffers do not recursively expand from other omitted observations.
#' @param x A `naboer_folds` object.
#' @inheritParams neighbourhoods
#' @param key,index Optional tidy selections overriding the stored series key
#'   and temporal index for buffer construction.
#' @return A `naboer_folds` object with expanded omission membership.
#' @export
buffer_folds <- function(x, spatial = NULL, temporal = NULL, phylogenetic = NULL,
                         combine = NULL, key = NULL, index = NULL, within = NULL) {
  validate_folds(x)
  ctx <- x$ctx
  if (!rlang::quo_is_null(rlang::enquo(key))) ctx$key <- select_cols(ctx$data, rlang::enquo(key), "key")
  if (!rlang::quo_is_null(rlang::enquo(index))) ctx$index <- select_cols(ctx$data, rlang::enquo(index), "index")
  count <- sum(!vapply(list(spatial, temporal, phylogenetic), is.null, logical(1)))
  if (!count) abort_naboer("Supply at least one buffer window.")
  if (count > 1L && is.null(combine)) abort_naboer("Specify how multiple buffer windows combine.")
  if (!is.null(combine)) combine <- match.arg(combine, c("intersection", "union"))
  parts <- list()
  if (!is.null(spatial)) parts$spatial <- spatial_membership(ctx, spatial)
  if (!is.null(temporal)) {
    if (is.null(within) && !is.null(spatial)) within <- "domain"
    parts$temporal <- temporal_membership(ctx, temporal, within)
  }
  if (!is.null(phylogenetic)) parts$phylogenetic <- phylogenetic_membership(ctx, phylogenetic)
  neighbourhood <- if (length(parts) == 1L) parts[[1L]] else combined_index(parts, combine)
  sets <- lapply(seq_len(index_nfold(x$assessment)), function(i) {
    centres <- index_get(x$assessment, i)
    additions <- unlist(lapply(centres, function(j) index_get(neighbourhood, j)), use.names = FALSE)
    sort(unique(c(index_get(x$omission, i), additions)))
  })
  spec <- x$spec
  spec$buffer <- list(spatial = spatial, temporal = temporal, phylogenetic = phylogenetic, combine = combine)
  new_folds(ctx, x$assessment, sets, spec)
}
