# These tests require only package Imports and testthat. Dependency checks are
# evaluated with the real helper body and a local namespace-availability shim;
# this avoids changing installed packages or global base namespace bindings.
test_that("core folds and temporal elapsed workflows need no specialist packages", {
  data <- tibble::tibble(time = 1:8, y = c(1, 3, 2, 4, 5, 7, 6, 8))
  expect_equal(as_nei(loo_folds(data))$a, 1:8)
  f <- neighbourhoods(data, temporal = temporal_window(time_span(1)), index = time)
  expect_equal(omission_indices(f, 4)[[1]], 3:5)
  expect_equal(fold_data(f), data)
  expect_s3_class(tibble::as_tibble(f), "tbl_df")
})

test_that("absent optional dependencies produce install guidance", {
  dependency_check <- need_package
  dependency_environment <- new.env(parent = environment(dependency_check))
  dependency_environment$requireNamespace <- function(package, quietly = TRUE) FALSE
  environment(dependency_check) <- dependency_environment
  testthat::local_mocked_bindings(need_package = dependency_check)
  data <- tibble::tibble(time = 1:4)
  expect_error(spatial_window(1), "Install.*units", class = "naboer_error_dependency")
  window <- structure(list(radius = 1), class = "naboer_spatial_window")
  expect_error(neighbourhoods(data, spatial = window), "Install.*sf", class = "naboer_error_dependency")
  expect_error(phylogenetic_window(NULL, taxon, 1), "Install.*ape", class = "naboer_error_dependency")
  expect_error(
    temporal_folds(tibble::tibble(date = as.Date("2024-01-01") + 0:3),
      time_span(1, "months", "calendar"),
      index = date
    ),
    "Install.*clock",
    class = "naboer_error_dependency"
  )
  expect_error(as_rset(loo_folds(data)), "Install.*rsample", class = "naboer_error_dependency")
  expect_error(autoplot.naboer_folds(loo_folds(data)), "Install.*ggplot2", class = "naboer_error_dependency")
  fake_tsibble <- data
  class(fake_tsibble) <- c("tbl_ts", class(data))
  expect_error(loo_folds(fake_tsibble), "Install.*tsibble", class = "naboer_error_dependency")
  fake_cubble <- data
  class(fake_cubble) <- c("cubble_df", class(data))
  expect_error(loo_folds(fake_cubble), "Install.*cubble", class = "naboer_error_dependency")
})

test_that("a genuinely absent namespace raises a structured dependency error", {
  expect_error(need_package("naboerDeliberatelyAbsentTestDependency"),
    "Install.*naboerDeliberatelyAbsentTestDependency",
    class = "naboer_error_dependency"
  )
})
