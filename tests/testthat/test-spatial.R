spatial_points <- function(x, y = rep(0, length(x)), crs = 3857, ...) {
  sf::st_as_sf(data.frame(x = x, y = y, ...), coords = c("x", "y"), crs = crs)
}

spatial_sets <- function(ctx, radius) {
  x <- spatial_membership(ctx, spatial_window(radius))
  lapply(seq_len(ctx$n), function(i) index_get(x, i))
}

test_that("spatial windows validate physical units", {
  skip_if_not_installed("units")
  expect_error(spatial_window(1), "units")
  expect_error(spatial_window(units::set_units(1, "s")), "distance")
  expect_error(spatial_window(units::set_units(-1, "m")), "nonnegative")
  expect_error(spatial_window(units::set_units(c(1, 2), "m")), "one finite")
})

test_that("indexed projected queries agree with brute-force distances", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  set.seed(421)
  points <- spatial_points(runif(30, -3, 3), runif(30, -3, 3),
    domain = rep(c("a", "b"), 15)
  )
  ctx <- prepare_data(points, by = rlang::quo(domain))
  radius <- units::set_units(1, "m")
  actual <- spatial_sets(ctx, radius)
  distances <- sf::st_distance(points)
  expected <- lapply(seq_len(nrow(points)), function(i) {
    unname(which(distances[i, ] <= radius & points$domain == points$domain[[i]]))
  })
  expect_equal(actual, expected)
  old <- options(naboer.spatial_batch_size = 3L)
  on.exit(options(old))
  expect_equal(spatial_sets(ctx, radius), expected)
  options(naboer.spatial_batch_size = 0)
  expect_error(spatial_sets(ctx, radius), "positive")
})

test_that("station caching preserves domains, duplicates, and inclusive boundaries", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(1, 0, 0, 2), domain = c(1, 1, 2, 1))
  ctx <- prepare_data(points, by = rlang::quo(domain))
  expect_equal(
    spatial_sets(ctx, units::set_units(1, "m")),
    list(c(1L, 2L, 4L), c(1L, 2L), 3L, c(1L, 4L))
  )
  expect_equal(
    spatial_sets(ctx, units::set_units(0, "m")),
    as.list(seq_len(4))
  )
})

test_that("polygons use boundary distance rather than centroids", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  square <- function(x) {
    sf::st_polygon(list(rbind(
      c(x, 0), c(x + 10, 0),
      c(x + 10, 1), c(x, 1), c(x, 0)
    )))
  }
  geom <- sf::st_sfc(square(0), square(11), square(100), crs = 3857)
  data <- sf::st_sf(site = 1:3, geometry = geom)
  expect_equal(
    spatial_sets(prepare_data(data), units::set_units(1, "m")),
    list(1:2, 1:2, 3L)
  )
  expect_equal(sf::st_geometry(data), geom)
})

test_that("geographic queries preserve the user's distance engine", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(0, 0.005, 1), crs = 4326)
  state <- sf::sf_use_s2()
  expect_equal(
    spatial_sets(prepare_data(points), units::set_units(1, "km")),
    list(1:2, 1:2, 3L)
  )
  expect_identical(sf::sf_use_s2(), state)
})

test_that("spatial inputs are never silently repaired", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(0:2, crs = NA)
  expect_error(spatial_sets(prepare_data(points), units::set_units(1, "m")), "CRS")
  points <- spatial_points(0:2)
  sf::st_geometry(points)[[2]] <- sf::st_point()
  expect_error(spatial_sets(prepare_data(points), units::set_units(1, "m")), "nonempty")
  expect_error(spatial_sets(
    prepare_data(tibble::tibble(x = 1:3)),
    units::set_units(1, "m")
  ), "geometry column")
})

test_that("square and hexagonal blocks are deterministic and respect domains", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(0, 1, 2, 4), domain = c(1, 1, 2, 2))
  for (shape in c("square", "hexagonal")) {
    folds <- spatial_folds(points, units::set_units(2, "m"),
      shape = shape,
      .by = domain
    )
    expect_s3_class(folds, "naboer_folds")
    expect_equal(
      folds$spec$blocks,
      spatial_folds(points, units::set_units(2, "m"),
        shape = shape,
        .by = domain
      )$spec$blocks
    )
    for (i in seq_len(index_nfold(folds$assessment))) {
      expect_length(unique(points$domain[index_get(folds$assessment, i)]), 1L)
    }
  }
  expect_error(spatial_folds(points, units::set_units(0, "m")), "positive")
  expect_error(spatial_folds(points, units::set_units(1, "m"), origin = 1), "x and y")
  expect_error(spatial_folds(points, units::set_units(1, "m"), origin = c(Inf, 1)), "x and y")
})

test_that("custom regions resolve shared boundaries but reject interior overlaps", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  rect <- function(left, right) {
    sf::st_polygon(list(rbind(
      c(left, -1), c(right, -1),
      c(right, 1), c(left, 1), c(left, -1)
    )))
  }
  regions <- sf::st_sfc(rect(-1, 1), rect(1, 3), crs = 3857)
  points <- spatial_points(c(0, 1, 2))
  folds <- spatial_folds(points, regions = regions)
  expect_equal(index_get(folds$assessment, 1), 1:2)
  overlap <- sf::st_sfc(rect(-1, 2), rect(1, 3), crs = 3857)
  expect_error(spatial_folds(points, regions = overlap), "overlap")
  points$region <- c(1L, 2L, 2L)
  assigned <- spatial_folds(points, regions = overlap, assignment = region)
  expect_equal(index_get(assigned$assessment, 2), 2:3)
  points$region <- c(3L, 2L, 2L)
  expect_error(spatial_folds(points, regions = overlap, assignment = region), "row numbers")
  expect_error(spatial_folds(spatial_points(c(0, 20)), regions = regions), "do not intersect")
})

test_that("spatiotemporal block modes have distinct prediction targets", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(0, 0, 3, 3), date = as.Date("2020-01-01") + c(0, 3, 0, 3))
  args <- list(
    data = points, spatial_size = units::set_units(2, "m"),
    temporal_size = time_span(2, "days"), index = "date"
  )
  crossed <- do.call(spatiotemporal_folds, args)
  space <- do.call(spatiotemporal_folds, c(args, list(mode = "space")))
  time <- do.call(spatiotemporal_folds, c(args, list(mode = "time")))
  expect_equal(index_nfold(crossed$assessment), 4L)
  expect_equal(index_nfold(space$assessment), 2L)
  expect_equal(index_nfold(time$assessment), 2L)
})

test_that("joint indexed searches match brute-force windows across multiple series", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(0, 1, 0, 0, 1, 2, 0, 2),
    date = as.Date("2020-01-01") + c(0, 0, 3, 1, 3, 0, 2, 3),
    station = c("a", "b", "a", "a", "b", "c", "a", "c"),
    domain = c(1, 1, 1, 2, 1, 1, 1, 1)
  )
  ctx <- prepare_data(points,
    key = rlang::quo(station), index = rlang::quo(date),
    by = rlang::quo(domain)
  )
  radius <- units::set_units(1, "m")
  spatial <- spatial_membership(ctx, spatial_window(radius))
  window <- temporal_window(time_span(1, "days"), time_span(2, "days"))
  for (within in c("domain", "key")) {
    joint <- spatiotemporal_membership(ctx, spatial, window, within)
    distances <- sf::st_distance(points)
    expected <- lapply(seq_len(nrow(points)), function(i) {
      same <- points$domain == points$domain[[i]]
      if (within == "key") same <- same & points$station == points$station[[i]]
      unname(which(distances[i, ] <= radius & same &
        points$date >= points$date[[i]] - 1 & points$date <= points$date[[i]] + 2))
    })
    actual <- lapply(seq_len(ctx$n), function(i) index_get(joint, i))
    expect_equal(actual, expected)
    expect_equal(index_sizes(joint), as.numeric(lengths(expected)))
    expect_equal(index_nfold(joint), nrow(points))
  }
})

test_that("geographic grids, explicit projected grids, and feet CRSs retain units", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  geographic <- spatial_points(c(12, 12.01, 12.02), c(55, 55.01, 55.02), crs = 4326)
  expect_s3_class(spatial_folds(geographic, units::set_units(1, "km")), "naboer_folds")
  expect_s3_class(spatial_folds(geographic, units::set_units(1, "km"),
    grid_crs = 32633
  ), "naboer_folds")
  expect_error(spatial_folds(geographic, units::set_units(1, "km"),
    grid_crs = 4326
  ), "projected")
  feet <- spatial_points(0:2, crs = 2263)
  radius <- units::set_units(1.01, "US_survey_foot")
  expect_equal(spatial_sets(prepare_data(feet), radius), list(1:2, 1:3, 2:3))
})

test_that("ellipsoidal distances preserve sf state and return exact geometry distances", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  skip_if_not_installed("lwgeom")
  old <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old)))
  suppressMessages(sf::sf_use_s2(FALSE))
  points <- spatial_points(c(0, 0.005, 1), crs = 4326)
  expect_equal(
    spatial_sets(prepare_data(points), units::set_units(1, "km")),
    list(1:2, 1:2, 3L)
  )
  expect_false(sf::sf_use_s2())
})

test_that("joint calendar and count windows respect temporal semantics", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  skip_if_not_installed("clock")
  points <- spatial_points(c(0, 1, 0, 1, 0, 1),
    date = as.Date(c(
      "2020-01-31", "2020-01-31", "2020-02-29",
      "2020-02-29", "2020-03-31", "2020-03-31"
    )),
    station = rep(c("a", "b"), 3)
  )
  ctx <- prepare_data(points, key = rlang::quo(station), index = rlang::quo(date))
  spatial <- spatial_membership(ctx, spatial_window(units::set_units(1, "m")))
  calendar <- temporal_window(time_span(1, "months", "calendar"), time_span(0, "days"))
  joint <- spatiotemporal_membership(ctx, spatial, calendar)
  expect_equal(index_get(joint, 5), 3:6)
  count <- temporal_window(time_span(1, "observations"))
  joint <- spatiotemporal_membership(ctx, spatial, count, "key")
  expect_equal(index_get(joint, 3), c(1L, 3L, 5L))
  expect_error(spatiotemporal_membership(ctx, spatial, count, "domain"), "within")
})

test_that("repeated-site spatial memberships retain compact ranges for dense outputs", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  data <- spatial_points(rep(0:9, 100), domain = rep(1:2, each = 500))
  ctx <- prepare_data(data, by = rlang::quo(domain))
  membership <- spatial_membership(ctx, spatial_window(units::set_units(100, "m")))
  expect_s3_class(membership$members, "naboer_multiranges")
  expect_true(as.numeric(object.size(membership)) < 100000)
  expect_equal(index_sizes(membership), rep(500, 1000))
  expect_equal(index_get(membership, 1), 1:500)
  expect_equal(index_get(membership, 1000), 501:1000)
})

test_that("geographic custom regions distinguish boundary contact from interior overlap", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  skip_if_not(sf::sf_use_s2(), "Spherical geometry disabled in the caller")
  rect <- function(left, right) {
    sf::st_polygon(list(rbind(
      c(left, 54), c(right, 54),
      c(right, 56), c(left, 56), c(left, 54)
    )))
  }
  points <- spatial_points(c(11.5, 12, 12.5), rep(55, 3), crs = 4326)
  regions <- sf::st_sfc(rect(11, 12), rect(12, 13), crs = 4326)
  expect_s3_class(spatial_folds(points, regions = regions), "naboer_folds")
  overlap <- sf::st_sfc(rect(11, 12.5), rect(12, 13), crs = 4326)
  expect_error(spatial_folds(points, regions = overlap), "overlap")
})

test_that("sf's active geometry is respected when several sfc columns exist", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(0, 1, 10))
  points$alternative <- sf::st_geometry(spatial_points(c(0, 100, 200)))
  expect_equal(
    spatial_sets(prepare_data(points), units::set_units(1, "m")),
    list(1:2, 1:2, 3L)
  )
})

test_that("precision metadata cannot discard focal observations or change radius distances", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(0.51, 0.51, 1.49, 10))
  sf::st_geometry(points) <- sf::st_set_precision(sf::st_geometry(points), 1)
  original <- sf::st_geometry(points)
  for (width in c(0, 0.1, 1)) {
    radius <- units::set_units(width, "m", mode = "standard")
    folds <- neighbourhoods(points, spatial = spatial_window(radius))
    distances <- sf::st_distance(points)
    expected <- lapply(seq_len(nrow(points)), function(i) unname(which(distances[i, ] <= radius)))
    expect_equal(omission_indices(folds), expected)
    expect_identical(fold_data(folds)$geometry, original)
  }
  blocks <- spatial_folds(points, units::set_units(1, "m"), origin = c(0, 0))
  expect_equal(index_get(blocks$assessment, 1), 1:2)
  expect_equal(index_get(blocks$assessment, 2), 3L)
})

test_that("fixed grid anchors extend into negative coordinates without changing cells", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(-3.5, -2.5, -1.5, 0.5, 1.5, 4.5))
  folds <- spatial_folds(points, units::set_units(2, "m"), origin = c(0, 0))
  expect_equal(assessment_indices(folds), list(1:2, 3L, 4:5, 6L))
  translated <- spatial_folds(points, units::set_units(2, "m"), origin = c(20, 20))
  expect_equal(assessment_indices(translated), assessment_indices(folds))
  expect_equal(folds$spec$origin, c(0, 0))
  # With an explicit anchor, subsetting cannot move the grid underneath the data.
  for (shape in c("square", "hexagonal")) {
    all <- spatial_folds(points, units::set_units(2, "m"), origin = c(0, 0), shape = shape)
    subset <- spatial_folds(points[-1, ], units::set_units(2, "m"), origin = c(0, 0), shape = shape)
    same <- outer(all$spec$blocks[-1], all$spec$blocks[-1], `==`)
    expect_equal(outer(subset$spec$blocks, subset$spec$blocks, `==`), same)
  }
})

test_that("grid-edge observations keep their cell after subsetting", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(-2, 0, 0.1, 2, 2.1, 4), c(0, 0, 0, 0, 0, 0))
  for (shape in c("square", "hexagonal")) {
    full <- spatial_folds(points, units::set_units(2, "m"), origin = c(0, 0), shape = shape)
    reduced <- spatial_folds(points[-1, ], units::set_units(2, "m"), origin = c(0, 0), shape = shape)
    expected <- outer(full$spec$blocks[-1], full$spec$blocks[-1], `==`)
    expect_equal(outer(reduced$spec$blocks, reduced$spec$blocks, `==`), expected)
  }
})

test_that("hexagonal lattice phase stays fixed across changing extents", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  set.seed(314)
  for (iteration in seq_len(10)) {
    xy <- matrix(runif(60, -10, 10), ncol = 2)
    points <- spatial_points(xy[, 1], xy[, 2])
    full <- spatial_folds(points, units::set_units(3, "m"),
      origin = c(0.5, -0.7), shape = "hexagonal"
    )
    keep <- sort(sample.int(30, 20))
    subset <- spatial_folds(points[keep, ], units::set_units(3, "m"),
      origin = c(0.5, -0.7), shape = "hexagonal"
    )
    expect_equal(
      outer(subset$spec$blocks, subset$spec$blocks, `==`),
      outer(full$spec$blocks[keep], full$spec$blocks[keep], `==`)
    )
  }
  # Distance across opposite edges is the specified width, not hex side length.
  cells <- spatial_hex_grid(sf::st_geometry(points), 3, c(0, 0))
  expected_area <- sqrt(3) / 2 * 3^2
  expect_equal(as.numeric(sf::st_area(cells)), rep(expected_area, length(cells)), tolerance = 1e-12)
})

test_that("grid allocations are checked before generating excessive cell counts", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(0, 2), c(0, 2))
  for (shape in c("square", "hexagonal")) {
    expect_error(spatial_folds(points, units::set_units(1e-12, "m"), shape = shape),
      class = "naboer_error_overflow"
    )
  }
  old <- options(naboer.max_bytes = 1)
  on.exit(options(old))
  for (shape in c("square", "hexagonal")) {
    expect_error(spatial_folds(points, units::set_units(1, "m"), shape = shape),
      class = "naboer_error_memory"
    )
  }
})

test_that("antimeridian distances work and geographic grids need an explicit projection", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(179.9, -179.9, 178), crs = 4326)
  expect_equal(
    spatial_sets(prepare_data(points), units::set_units(30, "km")),
    list(1:2, 1:2, 3L)
  )
  expect_error(spatial_folds(points, units::set_units(50, "km")), "grid_crs.*antimeridian")
  projection <- "+proj=aeqd +lon_0=180 +lat_0=0 +datum=WGS84 +units=m"
  expect_s3_class(
    spatial_folds(points, units::set_units(50, "km"), grid_crs = projection),
    "naboer_folds"
  )
})

test_that("spatiotemporal designs retain the parameters needed to reproduce blocks", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(-1, -1, 3, 3),
    date = as.Date("2020-01-01") + c(0, 4, 0, 4)
  )
  time_origin <- as.Date("2020-01-01")
  folds <- spatiotemporal_folds(points,
    spatial_size = units::set_units(2, "m"), temporal_size = time_span(2, "days"),
    index = date, spatial_origin = c(0, 0), temporal_origin = time_origin,
    shape = "hexagonal", within = "domain"
  )
  expect_equal(folds$spec$shape, "hexagonal")
  expect_equal(folds$spec$spatial_origin, c(0, 0))
  expect_equal(folds$spec$temporal_origin, time_origin)
  expect_equal(folds$spec$within, "domain")
  expect_equal(folds$spec$grid_crs, sf::st_crs(points))
})

test_that("block metadata report the engine and CRS actually used after projection", {
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  points <- spatial_points(c(12, 12.01, 12.02), c(55, 55.01, 55.02), crs = 4326)
  points$date <- as.Date("2020-01-01") + 0:2
  for (folds in list(
    spatial_folds(points, units::set_units(500, "m")),
    spatiotemporal_folds(points, units::set_units(500, "m"), time_span(1, "days"), index = date)
  )) {
    expect_identical(folds$spec$spatial_engine, "GEOS")
    expect_false(sf::st_is_longlat(folds$spec$grid_crs))
  }
})
