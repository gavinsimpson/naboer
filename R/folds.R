#' Construct custom cross-validation folds
#'
#' Assessment observations are scored; omission observations are excluded from
#' fitting. Every assessment set must be a subset of its omission set. Indices
#' refer to [fold_data()] order and must be unique within a set.
#' @param data A data frame, tibble, sf object, tsibble or cubble.
#' @param assessment,omission Lists of observation indices, one vector per fold.
#'   Omission defaults to assessment.
#' @param .by Tidy selection of independent domain columns. These constrain
#'   generated dependence queries and qualify group/block identities. Custom
#'   indices and existing omissions are retained as supplied; V-fold allocation
#'   can combine whole blocks from several domains in one fold.
#' @param id Optional tidy selection of a unique observation identifier.
#' @return A compact `naboer_folds` object.
#' @export
#' @examples
#' x <- tibble::tibble(time = 1:6, y = c(1, 3, 2, 4, 6, 5))
#' folds <- custom_folds(x, list(2L, 5L), list(1:3, 4:6))
#' as_nei(folds)
custom_folds <- function(data, assessment, omission = assessment, .by = NULL, id = NULL) {
  ctx <- prepare_data(data, by = rlang::enquo(.by), id = rlang::enquo(id))
  check_sets <- function(sets) {
    if (!is.list(sets)) abort_naboer("Supply membership as a list of numeric index vectors.", "indices")
    lapply(sets, function(z) {
      if (!is.numeric(z) || anyNA(z) || any(!is.finite(z)) ||
        any(z != floor(z)) || any(z < 1 | z > ctx$n) || anyDuplicated(z)) {
        abort_naboer("Indices must be unique whole numbers in the observation range.", "indices")
      }
      sort(as.integer(z))
    })
  }
  new_folds(ctx, check_sets(assessment), check_sets(omission), list(strategy = "custom"))
}

new_folds <- function(ctx, assessment, omission, spec = list()) {
  assessment <- list_index(assessment)
  x <- structure(
    list(
      ctx = ctx, assessment = assessment,
      omission = list_index(omission), spec = spec,
      labels = as.character(seq_len(index_nfold(assessment)))
    ),
    class = "naboer_folds"
  )
  validate_folds(x)
  x
}

#' Validate fold membership
#' @param x A `naboer_folds` object.
#' @return `x`, invisibly, or a structured error.
#' @export
validate_folds <- function(x) {
  if (!inherits(x, "naboer_folds") || !is.list(x$ctx) ||
    !inherits(x$assessment, "naboer_index") || !inherits(x$omission, "naboer_index")) {
    abort_naboer("Expected a `naboer_folds` object.", "structure")
  }
  if (!is.numeric(x$ctx$n) || length(x$ctx$n) != 1L || is.na(x$ctx$n) ||
    !is.finite(x$ctx$n) || x$ctx$n < 2 || x$ctx$n != floor(x$ctx$n) ||
    !is.data.frame(x$ctx$data) || nrow(x$ctx$data) != x$ctx$n) {
    abort_naboer("Fold observations and their stored row count do not match.", "structure")
  }
  validate_index(x$assessment, x$ctx$n)
  validate_index(x$omission, x$ctx$n)
  m <- index_nfold(x$assessment)
  if (!m || m != index_nfold(x$omission) || length(x$labels) != m) abort_naboer("Fold counts do not match.", "structure")
  assessment_sizes <- index_sizes(x$assessment)
  omission_sizes <- index_sizes(x$omission)
  if (anyNA(c(assessment_sizes, omission_sizes)) ||
    any(!is.finite(c(assessment_sizes, omission_sizes))) ||
    any(c(assessment_sizes, omission_sizes) < 0)) {
    abort_naboer("Membership sizes must be finite and nonnegative.", "indices")
  }
  if (any(assessment_sizes < 1)) abort_naboer("Assessment sets must be nonempty.", "indices")
  if (any(omission_sizes >= x$ctx$n)) abort_naboer("Every fold must leave a nonempty analysis set.", "indices")
  # Compact stores have vectorised validation; general sets are checked fold
  # by fold. No expanded list is retained during validation.
  a <- x$assessment
  o <- x$omission
  if (!validate_compact_membership(x)) {
    for (i in seq_len(m)) {
      aa <- index_get(a, i)
      oo <- index_get(o, i)
      if (anyNA(c(aa, oo)) || any(c(aa, oo) < 1L | c(aa, oo) > x$ctx$n) ||
        anyDuplicated(aa) || anyDuplicated(oo)) {
        abort_naboer("Invalid or duplicated membership indices.", "indices")
      }
      if (!all(aa %in% oo)) abort_naboer("Assessment must be included in omission.", "indices")
    }
  }
  invisible(x)
}

#' Leave-one-observation-out folds
#' @inheritParams custom_folds
#' @return A compact `naboer_folds` object.
#' @export
loo_folds <- function(data, .by = NULL, id = NULL) {
  ctx <- prepare_data(data, by = rlang::enquo(.by), id = rlang::enquo(id))
  new_folds(ctx, point_index(ctx$n), point_index(ctx$n), list(strategy = "leave-one-out"))
}

#' Random V-fold cross-validation
#' @inheritParams custom_folds
#' @param v Number of folds, at least two.
#' @param group Optional tidy selection of columns kept together during allocation.
#' @param balance Balance observation counts rather than group counts.
#' @param seed Optional local random seed; does not change the caller's RNG state.
#' @details `.by` qualifies group identities; it does not create a separate
#'   V-fold schedule within each domain. Every observation is assessed once
#'   across the requested `v` folds.
#' @return A compact `naboer_folds` object.
#' @export
vfold_folds <- function(data, v = 10L, group = NULL, balance = TRUE, seed = NULL,
                        .by = NULL, id = NULL) {
  ctx <- prepare_data(data, by = rlang::enquo(.by), id = rlang::enquo(id))
  cols <- select_cols(ctx$data, rlang::enquo(group), "group")
  if (length(cols) && anyNA(ctx$data[cols])) abort_naboer("Group columns contain missing values.", "missing")
  block <- if (length(cols)) group_ids(ctx$data, cols) else seq_len(ctx$n)
  group_folds_from_blocks(
    ctx, block, v, balance, seed,
    list(strategy = "vfold", group = cols)
  )
}

group_folds_from_blocks <- function(ctx, block, v = NULL, balance = TRUE, seed = NULL, spec = list()) {
  if (!is.logical(balance) || length(balance) != 1L || is.na(balance)) {
    abort_naboer("`balance` must be TRUE or FALSE.")
  }
  if (length(block) != ctx$n || anyNA(block)) abort_naboer("Blocks must identify every observation.", "indices")
  labels <- group_ids(tibble::tibble(domain = ctx$domains, block = block), c("domain", "block"))
  units <- split(seq_len(ctx$n), labels)
  if (!is.null(v)) {
    v <- check_count(v, "v", 2L)
    if (v > length(units)) abort_naboer("`v` exceeds the number of blocks or groups.")
    assignment <- with_seed(seed, {
      sizes <- lengths(units)
      ord <- sample.int(length(units))
      if (isTRUE(balance)) ord <- ord[order(-sizes[ord], method = "radix")]
      load <- numeric(v)
      out <- integer(length(units))
      for (j in ord) {
        k <- which.min(load)
        out[[j]] <- k
        load[[k]] <- load[[k]] + if (isTRUE(balance)) sizes[[j]] else 1
      }
      out
    })
    sets <- lapply(split(seq_along(units), assignment), function(j) sort(unlist(units[j], use.names = FALSE)))
  } else {
    sets <- unname(units)
  }
  spec$blocks <- block
  spec$v <- v
  spec$balance <- balance
  if (!is.null(seed) || is.null(spec$seed)) spec$seed <- seed
  new_folds(ctx, sets, sets, spec)
}

#' Retrieve the canonical observation data
#' @param x A `naboer_folds` object.
#' @return A tibble in the exact order used by exported indices. For cubbles,
#'   `attr(result, "naboer_source")` maps observations to source rows.
#' @export
fold_data <- function(x) {
  if (!inherits(x, "naboer_folds")) abort_naboer("Expected a `naboer_folds` object.", "structure")
  out <- x$ctx$data
  if (!is.null(x$ctx$source)) attr(out, "naboer_source") <- x$ctx$source
  out
}

fold_selection <- function(x, folds) {
  if (is.null(folds)) {
    return(seq_len(index_nfold(x$assessment)))
  }
  if (!is.numeric(folds) || anyNA(folds) || any(folds != floor(folds)) ||
    any(folds < 1 | folds > index_nfold(x$assessment))) {
    abort_naboer("Invalid fold selection.", "indices")
  }
  as.integer(folds)
}

#' Retrieve observation indices for folds
#' @param x A `naboer_folds` object.
#' @param folds Optional integer vector selecting folds. Defaults to all folds.
#' @return A list of integer vectors. Expanding all large folds can use substantial memory.
#' @export
assessment_indices <- function(x, folds = NULL) {
  lapply(fold_selection(x, folds), function(i) sort(index_get(x$assessment, i)))
}

#' @rdname assessment_indices
#' @export
omission_indices <- function(x, folds = NULL) {
  lapply(fold_selection(x, folds), function(i) sort(index_get(x$omission, i)))
}

#' @rdname assessment_indices
#' @export
analysis_indices <- function(x, folds = NULL) {
  lapply(fold_selection(x, folds), function(i) setdiff(seq_len(x$ctx$n), index_get(x$omission, i)))
}
