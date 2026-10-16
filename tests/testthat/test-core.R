test_that("custom assessment and omission are distinct and export correctly", {
  x <- tibble::tibble(id = letters[1:8], y = 1:8)
  folds <- custom_folds(x, list(c(3L, 1L), 6L), list(c(4L, 1L, 3L), 5:7), id = id)
  expect_s3_class(folds, "naboer_folds")
  expect_invisible(validate_folds(folds))
  expect_equal(fold_data(folds), x)
  expect_equal(assessment_indices(folds), list(c(1L, 3L), 6L))
  expect_equal(omission_indices(folds, c(2L, 1L)), list(5:7, c(1L, 3L, 4L)))
  expect_equal(analysis_indices(folds), list(c(2L, 5:8), c(1:4, 8L)))
  expect_equal(assessment_indices(folds, integer()), list())
  nei <- as_nei(folds)
  expect_identical(nei, list(
    a = c(1L, 3L, 4L, 5:7), ma = c(3L, 6L),
    d = c(1L, 3L, 6L), md = c(2L, 3L)
  ))
  expect_invisible(validate_nei(nei, nrow(x)))
})

test_that("assessment repetition across folds preserves its criterion weight", {
  x <- tibble::tibble(y = 1:5)
  folds <- custom_folds(x, list(2L, 2L), list(1:2, 2:3))
  nei <- as_nei(folds)
  expect_identical(nei$d, c(2L, 2L))
  expect_identical(nei$md, 1:2)
})

test_that("LOO and random V-fold designs cover observations exactly once", {
  x <- tibble::tibble(group = rep(letters[1:6], each = 2), y = 1:12)
  loo <- loo_folds(x)
  expect_equal(assessment_indices(loo), as.list(1:12))
  expect_equal(as_nei(loo), list(a = 1:12, ma = 1:12, d = 1:12, md = 1:12))
  set.seed(100)
  old <- .Random.seed
  a <- vfold_folds(x, v = 3, group = group, seed = 12)
  expect_identical(.Random.seed, old)
  b <- vfold_folds(x, v = 3, group = group, seed = 12)
  expect_equal(assessment_indices(a), assessment_indices(b))
  expect_equal(sort(unlist(assessment_indices(a), use.names = FALSE)), 1:12)
  expect_equal(lengths(assessment_indices(a)), rep(4L, 3))
  expect_true(all(vapply(assessment_indices(a), function(i) {
    all(table(x$group[i]) == 2L)
  }, logical(1))))
  expect_length(assessment_indices(vfold_folds(x, v = 3, balance = FALSE, seed = 4)), 3)
  expect_error(vfold_folds(x, v = 0), "at least|>=")
  expect_error(vfold_folds(x, v = 7, group = group), "exceeds")
  expect_error(vfold_folds(x, v = 3, balance = NA), "TRUE or FALSE")
  expect_error(vfold_folds(x, v = 3, seed = -1), "seed")
})

test_that("construction and selection reject invalid memberships", {
  x <- tibble::tibble(y = 1:5)
  expect_error(custom_folds(x, 1:2), "list", class = "naboer_error_indices")
  for (bad in list(c(1L, 1L), 0L, 6L, 1.5, NA_integer_, Inf, "1")) {
    expect_error(custom_folds(x, list(bad)), class = "naboer_error_indices")
  }
  expect_error(custom_folds(x, list(integer())), "nonempty")
  expect_error(custom_folds(x, list(1L), list(2L)), "included")
  expect_error(custom_folds(x, list(1L), list(1:5)), "analysis")
  expect_error(custom_folds(x, list(1L), list(1L, 2L)), "counts")
  expect_error(custom_folds(x, list()), "counts")
  folds <- loo_folds(x)
  for (bad in list(0L, 6L, 1.5, NA_integer_, Inf, "a")) {
    expect_error(assessment_indices(folds, bad), class = "naboer_error_indices")
  }
  expect_error(fold_data(list()), class = "naboer_error_structure")
  expect_error(validate_folds(list()), class = "naboer_error_structure")
  invalid <- folds
  invalid$ctx$n <- NULL
  expect_error(validate_folds(invalid), class = "naboer_error_structure")
  invalid <- folds
  invalid$labels <- "one"
  expect_error(validate_folds(invalid), "counts")
})

test_that("data domains, identifiers and grouped input are explicit", {
  x <- tibble::tibble(id = letters[1:5], domain = c(1L, 1L, 2L, 2L, 2L), y = 1:5)
  expect_error(loo_folds(x[1, , drop = FALSE]), "two observations")
  grouped <- dplyr::group_by(x, domain)
  expect_error(loo_folds(grouped), "ungroup")
  expect_s3_class(loo_folds(grouped, .by = domain), "naboer_folds")
  bad <- x
  bad$id[[2]] <- bad$id[[1]]
  expect_error(loo_folds(bad, id = id), "uniquely", class = "naboer_error_alignment")
  bad <- x
  bad$domain[[2]] <- NA_integer_
  expect_error(loo_folds(bad, .by = domain), "missing", class = "naboer_error_missing")
  expect_error(loo_folds(x, id = c(id, y)), "at most one")
  folds <- loo_folds(x, id = id)
  expect_identical(as_nei(folds, data = x[c("id", "y")]), as_nei(folds))
  expect_error(as_nei(folds, data = x[5:1, ]), "order", class = "naboer_error_alignment")
  expect_error(as_nei(folds, data = x[-1, ]), "rows", class = "naboer_error_alignment")
  expect_error(as_nei(loo_folds(x), data = transform(x, y = y + 1)), "order")
})

test_that("packed export guards allocation and preserves optional settings", {
  folds <- loo_folds(tibble::tibble(y = 1:8))
  plain <- as_nei(folds)
  expect_null(plain$sample)
  expect_null(plain$jackknife)
  sampled <- as_nei(folds, sample = c(2L, 5L), jackknife = FALSE)
  expect_identical(sampled[c("a", "ma", "d", "md")], plain)
  expect_identical(sampled$sample, c(2L, 5L))
  expect_false(sampled$jackknife)
  expect_identical(as_nei(folds, sample = 3)$sample, 3L)
  expect_true(as_nei(folds, jackknife = TRUE)$jackknife)
  expect_error(as_nei(folds, max_bytes = 1), "exceed", class = "naboer_error_memory")
  expect_equal(as_nei(folds, max_bytes = Inf), plain)
  for (bad in list(-1, NA_real_, c(1, 2), "a")) {
    expect_error(as_nei(folds, max_bytes = bad), "max_bytes")
  }
  for (bad in list(0, NA_real_, c(1, 1), 9, 1.5, Inf, character())) {
    expect_error(as_nei(folds, sample = bad), "sample")
  }
  for (bad in list(1, NA, c(TRUE, FALSE))) {
    expect_error(as_nei(folds, jackknife = bad), "jackknife")
  }
  n <- 60000L
  starts <- rep(1L, n)
  ends <- rep(n - 1L, n)
  starts[[n]] <- 2L
  ends[[n]] <- n
  large <- new_folds(
    prepare_data(tibble::tibble(y = seq_len(n))), point_index(n),
    range_index(starts, ends)
  )
  expect_error(as_nei(large, max_bytes = Inf), "integer index", class = "naboer_error_overflow")
})

test_that("packed structure validation rejects broken endpoints and settings", {
  good <- as_nei(loo_folds(tibble::tibble(y = 1:5)))
  expect_error(validate_nei(list(), 5), "Supply all")
  expect_error(validate_nei(good, 1), "n")
  for (bad in list(NA_integer_, 0L, Inf, 1.5, "a")) {
    broken <- good
    broken$a <- bad
    expect_error(validate_nei(broken, 5), "Invalid.*indices")
  }
  broken <- good
  broken$ma[[2]] <- 1L
  expect_error(validate_nei(broken, 5), "endpoints")
  broken <- good
  broken$a[[1]] <- 6L
  expect_error(validate_nei(broken, 5), "bounds")
  broken <- good
  broken$d[[1]] <- 2L
  expect_error(validate_nei(broken, 5), "membership")
  broken <- good
  broken$jackknife <- NA
  expect_error(validate_nei(broken, 5), "jackknife")
  broken <- good
  broken$sample <- c(1L, 1L)
  expect_error(validate_nei(broken, 5), "sample")
})
