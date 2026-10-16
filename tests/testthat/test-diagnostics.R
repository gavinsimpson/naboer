test_that("diagnostics count assessment coverage and buffer exclusions", {
  x <- tibble::tibble(time = 1:7)
  f <- custom_folds(x, list(c(2L, 3L), 3L), list(1:4, 2:4))
  s <- summary(f)
  expect_equal(s$unassessed, 5L)
  expect_equal(s$repeated_assessment, 1L)
  expect_equal(s$assessment_frequency$frequency, c(0, 1, 2, 0, 0, 0, 0))
  view <- tibble::as_tibble(f, membership = TRUE, folds = 2L)
  expect_equal(view$n_buffer, 2)
  expect_equal(view$assessment, list(3L))
  expect_false("assessment" %in% names(tibble::as_tibble(f)))
  expect_output(print(f), "7 observations; 2 folds")
  g <- neighbourhoods(x, index = time, temporal = temporal_window(time_span(1)))
  expect_equal(summary(g)$assessment_frequency$frequency, rep(1, 7))
  expect_error(tibble::as_tibble(g, unused = 1), "must be empty")
})

test_that("buffers use assessment centres and preserve existing exclusions", {
  x <- tibble::tibble(time = 1:10)
  f <- custom_folds(x, list(3L), list(c(3L, 8L)))
  b <- buffer_folds(f, index = time, temporal = temporal_window(time_span(1)))
  expect_equal(omission_indices(b), list(c(2:4, 8L)))
  expect_equal(assessment_indices(b), list(3L))
  expect_error(buffer_folds(f), "at least one")
})

test_that("plotting returns valid temporal and spatial plots", {
  skip_if_not_installed("ggplot2")
  f <- loo_folds(tibble::tibble(time = 1:5))
  expect_s3_class(ggplot2::autoplot(f, folds = c(1, 2)), "ggplot")
  expect_error(ggplot2::autoplot(f, folds = integer()), "at least one")
  g <- temporal_folds(tibble::tibble(time = 1:6), time_span(2), index = time)
  expect_silent(ggplot2::ggplot_build(ggplot2::autoplot(g)))
  skip_if_not_installed("sf")
  spatial <- sf::st_as_sf(data.frame(x = 1:4, y = 0), coords = c("x", "y"), crs = 32632)
  p <- ggplot2::autoplot(loo_folds(spatial))
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("single-window and combined neighbourhood contracts are explicit", {
  x <- tibble::tibble(time = c(3, 1, 5, 2, 4))
  expect_error(neighbourhoods(x), "at least one")
  expect_error(neighbourhoods(x,
    temporal = temporal_window(time_span(1)),
    spatial = list()
  ), "combine")
  f <- neighbourhoods(x, index = time, temporal = temporal_window(time_span(1)))
  expect_equal(omission_indices(f), lapply(x$time, function(t) which(abs(x$time - t) <= 1)))
  expect_equal(as_nei(f)$md, 1:5)
})

test_that("spatial-only public designs validate and export shared station ranges", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  x <- sf::st_as_sf(tibble::tibble(
    station = rep(1:4, each = 3), time = rep(1:3, 4),
    x = rep(c(0, 100, 200, 1000), each = 3), y = 0
  ), coords = c("x", "y"), crs = 32632)
  folds <- neighbourhoods(x, spatial = spatial_window(units::set_units(100, "m")))
  expected <- lapply(rep(c(0, 100, 200, 1000), each = 3), function(value) {
    which(abs(rep(c(0, 100, 200, 1000), each = 3) - value) <= 100)
  })
  expect_equal(omission_indices(folds), expected)
  expect_invisible(validate_folds(folds))
  packed <- as_nei(folds)
  expect_equal(packed$a, as.integer(unlist(expected, use.names = FALSE)))
  expect_equal(summary(folds)$assessment_frequency$frequency, rep(1, 12))
  broken <- folds
  broken$omission$lookup[1:3] <- broken$omission$lookup[10]
  expect_error(validate_folds(broken), "included", class = "naboer_error_indices")
})

test_that("forecast summaries respect unsorted multiple assessment ranges", {
  x <- tibble::tibble(time = c(9, 1, 7, 2, 8, 3, 6, 4, 5, 10))
  f <- forecast_folds(x, index = time, initial = time_span(3), horizon = time_span(2))
  expected <- tabulate(unlist(assessment_indices(f), use.names = FALSE), nrow(x))
  expect_equal(summary(f)$assessment_frequency$frequency, expected)
})
