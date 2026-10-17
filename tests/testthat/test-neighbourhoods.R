test_that("public temporal neighbourhoods require explicit series semantics", {
  data <- tibble::tibble(
    network = rep(1:2, each = 6),
    site = rep(rep(c("a", "b"), each = 3), 2),
    time = rep(c(0, 1, 3), 4)
  )
  window <- temporal_window(time_span(1))
  expect_error(neighbourhoods(data, temporal = window, key = site, index = time), "Multiple series")
  by_series <- neighbourhoods(data,
    temporal = window, key = site, index = time,
    .by = network, within = "key"
  )
  by_domain <- neighbourhoods(data,
    temporal = window, key = site, index = time,
    .by = network, within = "domain"
  )
  expect_equal(omission_indices(by_series, 1)[[1]], 1:2)
  expect_equal(omission_indices(by_domain, 1)[[1]], c(1L, 2L, 4L, 5L))
  expect_equal(assessment_indices(by_domain), as.list(1:12))
  expect_equal(fold_data(by_domain), data)
})

test_that("public phylogenetic neighbourhoods and buffers respect taxon domains", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):2);")
  data <- tibble::tibble(
    taxon = rep(c("a", "b", "c", "d"), 2),
    network = rep(1:2, each = 4)
  )
  window <- phylogenetic_window(tree, taxon, 2)
  folds <- neighbourhoods(data, phylogenetic = window, .by = network)
  expect_equal(omission_indices(folds), list(1:2, 1:2, 3:4, 3:4, 5:6, 5:6, 7:8, 7:8))
  initial <- custom_folds(data, list(1L), list(c(1L, 7L)), .by = network)
  buffered <- buffer_folds(initial, phylogenetic = window)
  expect_equal(assessment_indices(buffered), list(1L))
  expect_equal(omission_indices(buffered), list(c(1L, 2L, 7L)))
  expect_equal(omission_indices(initial), list(c(1L, 7L)))
})

test_that("buffer centres are assessment observations and existing exclusions persist", {
  data <- tibble::tibble(time = 0:8, site = rep("a", 9))
  initial <- custom_folds(data, list(3L), list(c(3L, 7L)))
  buffered <- buffer_folds(initial,
    temporal = temporal_window(time_span(1)),
    key = site, index = time, within = "key"
  )
  expect_equal(omission_indices(buffered), list(c(2L, 3L, 4L, 7L)))
  expect_equal(assessment_indices(buffered), list(3L))
  expect_equal(buffered$ctx$key, "site")
  expect_equal(buffered$ctx$index, "time")
  expect_equal(analysis_indices(buffered), list(c(1L, 5L, 6L, 8L, 9L)))
})

test_that("buffer metadata overrides receive constructor validation", {
  data <- tibble::tibble(time = 1:6, site = c("a", NA, "a", "b", "b", "b"))
  initial <- loo_folds(data)
  expect_error(
    buffer_folds(initial,
      temporal = temporal_window(time_span(1)),
      key = site, index = time, within = "key"
    ),
    "Membership column `site` contains missing values",
    class = "naboer_error_missing"
  )
  expect_error(
    buffer_folds(initial,
      temporal = temporal_window(time_span(1)),
      index = c(time, site)
    ),
    "at most one",
    class = "naboer_error_argument"
  )
  expect_length(initial$ctx$index, 0)
  expect_length(initial$ctx$key, 0)
})

test_that("combined spatial and temporal neighbourhoods implement set operations", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  data <- sf::st_as_sf(
    tibble::tibble(
      network = rep(1:2, each = 4), site = rep(letters[1:4], 2),
      time = rep(c(0, 2, 3, 8), 2), x = rep(c(0, 1, 3, 6), 2), y = 0
    ),
    coords = c("x", "y"), crs = 3857
  )
  spatial <- spatial_window(units::set_units(1.1, "m"))
  temporal <- temporal_window(time_span(2))
  intersection <- neighbourhoods(data,
    spatial = spatial, temporal = temporal,
    combine = "intersection", key = site, index = time,
    .by = network
  )
  union <- neighbourhoods(data,
    spatial = spatial, temporal = temporal,
    combine = "union", key = site, index = time, .by = network
  )
  expect_equal(omission_indices(intersection, 2)[[1]], 1:2)
  expect_equal(omission_indices(union, 2)[[1]], 1:3)
  expect_equal(omission_indices(union, 6)[[1]], 5:7)
  expect_equal(intersection$spec$within, "domain")
  expect_equal(intersection$spec$spatial_engine, "GEOS")
  initial <- custom_folds(data, list(2L), list(c(2L, 4L)), .by = network)
  both <- buffer_folds(initial,
    spatial = spatial, temporal = temporal,
    combine = "intersection", key = site, index = time
  )
  either <- buffer_folds(initial,
    spatial = spatial, temporal = temporal,
    combine = "union", key = site, index = time
  )
  expect_equal(omission_indices(both), list(c(1L, 2L, 4L)))
  expect_equal(omission_indices(either), list(1:4))
  # Only the first independent domain is omitted; the second remains analysis.
  expect_equal(analysis_indices(either), list(5:8))
  space_only <- buffer_folds(initial, spatial = spatial)
  expect_equal(omission_indices(space_only), list(c(1L, 2L, 4L)))
})

test_that("all three dependence dimensions combine without losing exclusions", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):2);")
  data <- sf::st_as_sf(
    tibble::tibble(
      network = rep(1:2, each = 4), taxon = rep(letters[1:4], 2),
      time = rep(c(0, 2, 3, 8), 2), x = rep(c(0, 1, 3, 6), 2), y = 0
    ),
    coords = c("x", "y"), crs = 3857
  )
  spatial <- spatial_window(units::set_units(1.1, "m"))
  temporal <- temporal_window(time_span(2))
  phylogenetic <- phylogenetic_window(tree, taxon, 2)
  intersection <- neighbourhoods(data,
    spatial = spatial, temporal = temporal,
    phylogenetic = phylogenetic, combine = "intersection",
    index = time, .by = network
  )
  union <- neighbourhoods(data,
    spatial = spatial, temporal = temporal,
    phylogenetic = phylogenetic, combine = "union",
    index = time, .by = network
  )
  expect_equal(omission_indices(intersection, 3)[[1]], 3L)
  expect_equal(omission_indices(union, 3)[[1]], 2:4)
  initial <- custom_folds(data, list(3L), list(c(3L, 8L)), .by = network)
  buffered <- buffer_folds(initial,
    spatial = spatial, temporal = temporal,
    phylogenetic = phylogenetic, combine = "union", index = time
  )
  expect_equal(omission_indices(buffered), list(c(2L, 3L, 4L, 8L)))
  # Temporal+phylogenetic composition also works without a spatial adapter.
  temporal_tree <- neighbourhoods(sf::st_drop_geometry(data),
    temporal = temporal,
    phylogenetic = phylogenetic, combine = "intersection",
    index = time, .by = network
  )
  expect_equal(omission_indices(temporal_tree, 3)[[1]], 3L)
})

test_that("count windows support public spatial composition", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  data <- sf::st_as_sf(
    tibble::tibble(
      site = rep(c("a", "b"), each = 3),
      time = rep(c(0, 4, 20), 2),
      x = rep(c(0, 10), each = 3), y = 0
    ),
    coords = c("x", "y"), crs = 3857
  )
  folded <- neighbourhoods(data,
    spatial = spatial_window(units::set_units(1, "m")),
    temporal = temporal_window(time_span(1, "observations")),
    combine = "intersection", index = time, key = site,
    within = "key"
  )
  expect_equal(omission_indices(folded, 1)[[1]], 1:2)
  expect_equal(omission_indices(folded, 5)[[1]], 4:6)
})

test_that("public window and combination errors give useful messages", {
  data <- tibble::tibble(time = 1:6)
  window <- temporal_window(time_span(1))
  expect_error(neighbourhoods(data), "at least one")
  expect_error(buffer_folds(loo_folds(data)), "at least one")
  expect_error(neighbourhoods(data, temporal = window, spatial = 1, index = time), "combine")
  expect_error(buffer_folds(loo_folds(data), temporal = window, spatial = 1), "combine")
  expect_error(neighbourhoods(data, temporal = window, combine = "invalid", index = time), "arg")
  expect_error(buffer_folds(loo_folds(data), temporal = window, combine = "invalid", index = time), "arg")
  expect_error(neighbourhoods(data, temporal = 1, index = time), "temporal_window")
  expect_error(neighbourhoods(data, spatial = 1), "spatial_window")
  expect_error(neighbourhoods(data, phylogenetic = 1), "phylogenetic_window")
  expect_error(buffer_folds(loo_folds(data), temporal = 1, index = time), "temporal_window")
  expect_error(buffer_folds(loo_folds(data), spatial = 1), "spatial_window")
  expect_error(buffer_folds(loo_folds(data), phylogenetic = 1), "phylogenetic_window")
})

test_that("forecast validation keeps compressed future exclusions", {
  x <- tibble::tibble(time = seq_len(10000L))
  f <- forecast_folds(x,
    horizon = time_span(100), initial = time_span(1000),
    step = time_span(100), index = time
  )
  expect_s3_class(f$omission, "naboer_temporal_complement")
  # Validation must not expand the large future exclusions to check subsets.
  testthat::local_mocked_bindings(index_get = function(...) stop("Unexpected expansion"))
  expect_silent(validate_folds(f))
})
