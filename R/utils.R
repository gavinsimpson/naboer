abort_naboer <- function(message, class = "argument") {
  rlang::abort(message, class = c(paste0("naboer_error_", class), "naboer_error"))
}

need_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    abort_naboer(paste0("Install the optional package `", package, "` to use this feature."), "dependency")
  }
}

check_count <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
    x != floor(x) || x < minimum || x > .Machine$integer.max) {
    abort_naboer(paste0("`", name, "` must be a whole number >= ", minimum, "."))
  }
  as.integer(x)
}

select_cols <- function(data, quo, arg, allow_empty = TRUE) {
  out <- names(tidyselect::eval_select(quo, data, allow_rename = FALSE))
  if (!allow_empty && !length(out)) abort_naboer(paste0("Select at least one column for `", arg, "`."))
  out
}

group_ids <- function(data, columns) {
  if (!length(columns)) {
    return(rep.int(1L, nrow(data)))
  }
  as.integer(vctrs::vec_group_id(data[columns]))
}

domain_ids <- function(ctx, within = NULL) {
  group_ids(ctx$data, unique(c(ctx$by, if (identical(within, "key")) ctx$key)))
}

with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  seed <- check_count(seed, "seed", 0L)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)
  force(code)
}

prepare_data <- function(data, key = rlang::quo(NULL), index = rlang::quo(NULL),
                         by = rlang::quo(NULL), id = rlang::quo(NULL)) {
  if (!is.data.frame(data)) abort_naboer("`data` must be a data frame.")
  if (inherits(data, "grouped_df") && rlang::quo_is_null(by)) {
    abort_naboer("Specify independent domains with `.by`, or `ungroup()` the input.")
  }
  inferred_key <- inferred_index <- character()
  source <- NULL
  if (inherits(data, "cubble_df")) {
    need_package("cubble")
    cube <- cubble_data(data)
    data <- cube$data
    inferred_key <- cube$key
    inferred_index <- cube$index
    source <- cube$source
  } else if (inherits(data, "tbl_ts")) {
    need_package("tsibble")
    inferred_key <- tsibble::key_vars(data)
    inferred_index <- tsibble::index_var(data)
  }
  data <- tibble::as_tibble(data)
  if (nrow(data) < 2L) abort_naboer("`data` must contain at least two observations.")
  if (anyDuplicated(names(data))) abort_naboer("Column names must be unique.")
  keys <- if (rlang::quo_is_null(key)) inferred_key else select_cols(data, key, "key")
  times <- if (rlang::quo_is_null(index)) inferred_index else select_cols(data, index, "index")
  bys <- select_cols(data, by, ".by")
  ids <- select_cols(data, id, "id")
  ctx <- list(
    data = data, n = nrow(data), key = keys, index = times, by = bys,
    id = ids, domains = group_ids(data, bys), source = source
  )
  validate_context_columns(ctx)
  ctx
}

validate_context_columns <- function(ctx) {
  if (length(ctx$index) > 1L || length(ctx$id) > 1L) {
    abort_naboer("`index` and `id` must each select at most one column.")
  }
  for (nm in unique(c(ctx$key, ctx$index, ctx$by, ctx$id))) {
    if (anyNA(ctx$data[[nm]])) {
      abort_naboer(paste0("Membership column `", nm, "` contains missing values."), "missing")
    }
  }
  if (length(ctx$id) && anyDuplicated(ctx$data[[ctx$id]])) {
    abort_naboer("`id` must uniquely identify observations.", "alignment")
  }
  invisible(ctx)
}
