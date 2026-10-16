test_that("corrupted observation permutations are rejected before export", {
  x <- loo_folds(tibble::tibble(y = 1:6))
  x$omission <- range_index(1:6, 1:6, order = 1:6)
  expect_invisible(validate_folds(x))
  x$omission$order[[2]] <- 1L
  expect_error(validate_folds(x), class = "naboer_error_indices")
  x$omission <- range_index(1:6, 1:6, order = 6:1)
  expect_error(validate_folds(x), "included", class = "naboer_error_indices")
  x$omission <- range_index(1:6, 1:6)
  x$omission$start[[1]] <- 0L
  expect_error(validate_folds(x), class = "naboer_error_indices")
  x$omission$start[[1]] <- NA_integer_
  expect_error(validate_folds(x), class = "naboer_error_indices")
})

test_that("shared memberships retain validation for every focal row", {
  x <- loo_folds(tibble::tibble(y = 1:6))
  x$omission <- shared_index(list_index(list(1:3, 4:6)), c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_invisible(validate_folds(x))
  x$omission$lookup[[1]] <- 2L
  expect_error(validate_folds(x), "included", class = "naboer_error_indices")
  x$omission$lookup[[1]] <- 3L
  expect_error(validate_folds(x), class = "naboer_error_indices")
  x$omission <- shared_index(list(c(1L, 1L, 2L, 3L), 4:6), c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_error(validate_folds(x), class = "naboer_error_indices")
})

test_that("multirange validation detects overlap, missing focal rows and stale sizes", {
  folds <- loo_folds(tibble::tibble(y = 1:6))
  positions <- c(1L, 4L, 2L, 5L, 3L, 6L, 1L, 4L, 2L, 5L, 3L, 6L)
  folds$omission <- structure(
    list(
      start = positions, end = positions, ends = seq.int(2L, 12L, by = 2L),
      sizes = rep(2, 6), order = 1:6
    ),
    class = c("naboer_multiranges", "naboer_index")
  )
  expect_invisible(validate_folds(folds))
  expect_equal(omission_indices(folds, 1), list(c(1L, 4L)))
  bad <- folds
  bad$omission$ends[[2]] <- 2L
  expect_error(validate_folds(bad), class = "naboer_error_indices")
  bad <- folds
  bad$omission$start[[2]] <- 1L
  expect_error(validate_folds(bad), class = "naboer_error_indices")
  bad <- folds
  bad$omission$start[[1]] <- 2L
  bad$omission$end[[1]] <- 2L
  expect_error(validate_folds(bad), "included", class = "naboer_error_indices")
  bad <- folds
  bad$omission$sizes[[1]] <- 1
  expect_error(validate_folds(bad), class = "naboer_error_indices")
  bad <- folds
  bad$omission$order <- c(1L, 1L, 3:6)
  expect_error(validate_folds(bad), class = "naboer_error_indices")
})

test_that("forecast compact validation prevents training-assessment overlap", {
  x <- tibble::tibble(time = seq_len(12))
  folds <- forecast_folds(x,
    index = time, initial = time_span(4, "observations", "count"),
    horizon = time_span(2, "observations", "count"),
    step = time_span(2, "observations", "count"), within = "key"
  )
  expect_invisible(validate_folds(folds))
  bad <- folds
  bad$omission$order <- rev(bad$omission$order)
  expect_error(validate_folds(bad), class = "naboer_error_indices")
  bad <- folds
  bad$omission$domain_start[[1]] <- 2L
  expect_error(validate_folds(bad), class = "naboer_error_indices")
  bad <- folds
  bad$assessment$start[[1]] <- bad$omission$start[[1]]
  expect_error(validate_folds(bad), "within its domain", class = "naboer_error_indices")
})

test_that("dependencies and imported corrupted indices fail clearly", {
  expect_error(need_package("naboer.optional.package.not.installed"), "Install the optional",
    class = "naboer_error_dependency"
  )
  skip_if_not_installed("rsample")
  split <- rsample::make_splits(
    list(analysis = 1:3, assessment = 4:6),
    tibble::tibble(y = 1:6)
  )
  split$in_id[[1]] <- 1.5
  resamples <- rsample::manual_rset(list(split), "corrupted")
  expect_error(from_rset(resamples), "indices", class = "naboer_error_argument")
})

test_that("a local random seed does not create global random state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) saved <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", saved, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
  x <- tibble::tibble(group = rep(letters[1:3], each = 2), time = 1:6)
  expect_s3_class(vfold_folds(x, v = 3, seed = 123), "naboer_folds")
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_error(group_folds(x, tidyselect::all_of(character())), "Select at least one")
  x$group[[1]] <- NA_character_
  expect_error(vfold_folds(x, v = 2, group = group), "missing")
})

test_that("temporal-phylogenetic combinations and buffers follow set semantics", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):2,c:1);")
  x <- tibble::tibble(taxon = rep(c("a", "b", "c"), each = 2), time = 1:6)
  phylo <- phylogenetic_window(tree, taxon, 0)
  time <- temporal_window(time_span(1, "seconds"))
  union <- neighbourhoods(x,
    phylogenetic = phylo, temporal = time,
    index = time, combine = "union"
  )
  intersection <- neighbourhoods(x,
    phylogenetic = phylo, temporal = time,
    index = time, combine = "intersection"
  )
  expect_equal(omission_indices(union, 2), list(1:3))
  expect_equal(omission_indices(intersection, 2), list(1:2))
  buffered <- buffer_folds(loo_folds(x),
    temporal = time, phylogenetic = phylo,
    index = time, combine = "union"
  )
  expect_equal(omission_indices(buffered), omission_indices(union))
  expect_error(buffer_folds(loo_folds(x),
    temporal = time, phylogenetic = phylo,
    index = time
  ), "multiple buffer")
  expect_error(neighbourhoods(x), "at least one")
  expect_error(neighbourhoods(x,
    temporal = time, phylogenetic = phylo,
    index = time
  ), "combine")
  expect_error(
    neighbourhoods(x, phylogenetic = phylogenetic_window(tree, c(taxon, time), 0)),
    "exactly one.*taxon"
  )
  clades <- phylogenetic_folds(x, tree, taxon, clades = as.list(1:3))
  expect_equal(assessment_indices(clades), list(1:2, 3:4, 5:6))
})

test_that("packed and combined stores cannot export corrupt membership", {
  folds <- custom_folds(tibble::tibble(y = 1:6), list(c(1L, 3L), c(4L, 6L)))
  expect_s3_class(folds$assessment, "naboer_packed")
  for (bad in list(NA_integer_, 0L, 7L, 1.5)) {
    broken <- folds
    broken$assessment$values[[1]] <- bad
    expect_error(as_nei(broken), class = "naboer_error_indices")
  }
  broken <- folds
  broken$assessment$ends[[2]] <- 100L
  expect_error(as_nei(broken), class = "naboer_error_indices")
  broken <- folds
  broken$assessment$ends <- rev(broken$assessment$ends)
  expect_error(as_nei(broken), class = "naboer_error_indices")
  broken <- loo_folds(tibble::tibble(y = 1:6))
  broken$omission <- combined_index(list(point_index(6), point_index(6)), "union")
  expect_equal(as_nei(broken), as_nei(loo_folds(tibble::tibble(y = 1:6))))
  broken$omission$combine <- "unknown"
  expect_error(as_nei(broken), class = "naboer_error_indices")
  broken$omission$combine <- "union"
  broken$omission$indices[[2]] <- point_index(5)
  expect_error(as_nei(broken), class = "naboer_error_indices")
  broken$omission <- structure(list(), class = "naboer_index")
  expect_error(as_nei(broken), class = "naboer_error_indices")
})
