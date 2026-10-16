test_that("cubble adapters infer keys/index and preserve unequal series mappings", {
  skip_if_not_installed("cubble")
  skip_if_not_installed("tsibble")
  skip_if_not_installed("sf")
  data <- cubble::climate_mel
  data$ts <- lapply(seq_len(nrow(data)), function(i) data$ts[[i]][seq_len(i + 1L), ])
  nested <- cubble_data(data)
  long_data <- cubble::face_temporal(data)
  long <- cubble_data(long_data)
  expect_equal(nested$data, long$data)
  expect_identical(nested$key, "id")
  expect_identical(nested$index, "date")
  expect_equal(nested$source$.source_row, rep(1:3, 2:4))
  expect_equal(nested$source$.source_subrow, sequence(2:4))
  expect_equal(long$source$.source_row, seq_len(9))
  ctx <- prepare_data(data)
  expect_identical(ctx$key, "id")
  expect_identical(ctx$index, "date")
  expect_equal(ctx$source, nested$source)
  expect_equal(ctx$data$date, unlist(lapply(data$ts, `[[`, "date")) |> as.Date(origin = "1970-01-01"))
})

test_that("long cubble row order survives canonicalisation", {
  skip_if_not_installed("cubble")
  skip_if_not_installed("sf")
  long <- cubble::face_temporal(cubble::climate_mel)
  order <- c(30:20, 1:19)
  long <- long[order, ]
  converted <- cubble_data(long)
  expect_equal(converted$data$id, long$id)
  expect_equal(converted$data$date, long$date)
})

test_that("cubble's original polygons survive nested and temporal forms", {
  skip_if_not_installed("cubble")
  skip_if_not_installed("sf")
  skip_if_not_installed("units")
  polygon <- function(x) {
    sf::st_polygon(list(rbind(
      c(x, 0), c(x + 10, 0),
      c(x + 10, 1), c(x, 1), c(x, 0)
    )))
  }
  stations <- sf::st_sf(
    station = c("a", "b"), x = c(5, 16), y = 0.5,
    geometry = sf::st_sfc(polygon(0), polygon(11), crs = 3857)
  )
  times <- tibble::tibble(
    station = c("a", "a", "b"),
    date = as.Date("2020-01-01") + c(0, 2, 0), value = 1:3
  )
  cube <- cubble::make_cubble(stations, times,
    key = station, index = date,
    coords = c(x, y)
  )
  nested <- cubble_data(cube)
  temporal <- cubble_data(cubble::face_temporal(cube))
  expect_equal(nested$data, temporal$data)
  geometry <- nested$data[[which(vapply(nested$data, inherits, logical(1), "sfc"))]]
  expect_true(all(sf::st_geometry_type(geometry) == "POLYGON"))
  expect_equal(geometry, sf::st_geometry(stations)[c(1, 1, 2)])
  context <- prepare_data(cube)
  membership <- spatial_membership(context, spatial_window(units::set_units(1, "m")))
  expect_equal(index_get(membership, 1), 1:3)
})

test_that("one-station nested cubbles are valid multiple-observation data", {
  skip_if_not_installed("cubble")
  skip_if_not_installed("sf")
  cube <- cubble::climate_mel[1, ]
  folds <- temporal_folds(cube, size = time_span(2, "days"))
  expect_equal(nrow(fold_data(folds)), 10L)
  expect_equal(index_nfold(folds$assessment), 5L)
})
