test_that("rsample round trips preserve exclusion-only observations", {
  skip_if_not_installed("rsample")
  x <- tibble::tibble(y = 1:8)
  original <- custom_folds(x, list(2L, 6L), list(1:3, 5:7))
  resamples <- as_rset(original)
  expect_s3_class(resamples, "manual_rset")
  expect_equal(as.integer(resamples$splits[[1]], data = "analysis"), 4:8)
  expect_equal(as.integer(resamples$splits[[1]], data = "assessment"), 2L)
  imported <- from_rset(resamples)
  expect_equal(assessment_indices(imported), assessment_indices(original))
  expect_equal(omission_indices(imported), omission_indices(original))
  expect_equal(analysis_indices(imported), analysis_indices(original))
  expect_equal(fold_data(imported), fold_data(original))
})

test_that("ordinary V-fold and Monte Carlo resamples import faithfully", {
  skip_if_not_installed("rsample")
  x <- tibble::tibble(id = letters[1:12], y = 1:12)
  set.seed(5)
  original <- rsample::vfold_cv(x, v = 3)
  imported <- from_rset(original, id = id)
  expect_equal(imported$labels, original$id)
  for (i in seq_len(nrow(original))) {
    expect_equal(
      analysis_indices(imported, i)[[1]],
      sort(as.integer(original$splits[[i]], data = "analysis"))
    )
    expect_equal(
      assessment_indices(imported, i)[[1]],
      sort(as.integer(original$splits[[i]], data = "assessment"))
    )
  }
  repeated <- rsample::vfold_cv(x, v = 3, repeats = 2)
  imported <- from_rset(repeated)
  expect_length(unique(imported$labels), 6L)
})

test_that("unrepresentable and invalid resamples are rejected", {
  skip_if_not_installed("rsample")
  x <- tibble::tibble(y = 1:8)
  split <- function(analysis, assessment, data = x) {
    rsample::make_splits(list(analysis = analysis, assessment = assessment), data)
  }
  rset <- function(...) rsample::manual_rset(list(...), as.character(seq_along(list(...))))
  expect_error(from_rset(rset(split(c(1L, 1L, 2L), 3:8))), "Duplicated.*bootstrap")
  expect_error(from_rset(rset(split(1:4, c(5L, 5L)))), "unique")
  expect_error(from_rset(rset(split(1:4, 4:8))), "overlap")
  expect_error(from_rset(rset(split(1:4, 9L))), "indices")
  expect_error(from_rset(rset(split(1:4, integer()))), "nonempty")
  expect_error(from_rset(rset(
    split(1:4, 5:8),
    split(1:4, 5:8, tibble::tibble(y = 8:1))
  )), "identical")
  expect_error(from_rset(x), "nonempty rsample")
  malformed <- rset(split(1:4, 5:8))
  malformed$splits[[1]] <- list()
  expect_error(from_rset(malformed), "rsplit")
})
