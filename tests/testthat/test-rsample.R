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

test_that("nested cubble source-row resamples cannot silently change observation units", {
  skip_if_not_installed("rsample")
  skip_if_not_installed("cubble")
  skip_if_not_installed("sf")
  nested <- cubble::climate_mel
  split <- rsample::make_splits(list(analysis = 2:3, assessment = 1L), nested)
  resamples <- rsample::manual_rset(list(split), "station1")
  expect_error(from_rset(resamples), "individual observations.*changes its rows",
    class = "naboer_error_alignment"
  )
  # Equal source/canonical row counts can still hide a nonidentity mapping.
  nested$ts <- list(nested$ts[[1]][1:2, ], nested$ts[[2]][integer(), ], nested$ts[[3]][1, ])
  split <- rsample::make_splits(list(analysis = 2:3, assessment = 1L), nested)
  expect_error(from_rset(rsample::manual_rset(list(split), "station1")),
    "face_temporal",
    class = "naboer_error_alignment"
  )
})

test_that("observation-level long cubble and canonical exports import faithfully", {
  skip_if_not_installed("rsample")
  skip_if_not_installed("cubble")
  skip_if_not_installed("sf")
  long <- cubble::face_temporal(cubble::climate_mel)
  long <- long[rev(seq_len(nrow(long))), ]
  split <- rsample::make_splits(list(analysis = 4:30, assessment = 1:2), long)
  imported <- from_rset(rsample::manual_rset(list(split), "buffered"))
  expect_equal(assessment_indices(imported), list(1:2))
  expect_equal(omission_indices(imported), list(1:3))
  expect_equal(analysis_indices(imported), list(4:30))
  expect_equal(fold_data(imported)$date, long$date)
  original <- temporal_folds(cubble::climate_mel, time_span(2, "days"), within = "key")
  restored <- from_rset(as_rset(original))
  expect_equal(assessment_indices(restored), assessment_indices(original))
  expect_equal(analysis_indices(restored), analysis_indices(original))
  expect_equal(fold_data(restored), fold_data(original))
})

test_that("rsample conversion guards large analysis complements before expansion", {
  skip_if_not_installed("rsample")
  # Native LOO storage is linear, but these analysis complements would require
  # about 3.6 GB of integer indices. The error must occur without allocating it.
  large <- loo_folds(tibble::tibble(y = seq_len(30000L)))
  expect_error(as_rset(large), "before R object overhead.*exceeding",
    class = "naboer_error_memory"
  )
  small <- loo_folds(tibble::tibble(y = 1:8))
  expect_error(as_rset(small, max_bytes = 1), "Reduce the design",
    class = "naboer_error_memory"
  )
  old <- options(naboer.max_bytes = 1)
  on.exit(options(old))
  expect_error(as_rset(small), class = "naboer_error_memory")
  expect_equal(
    analysis_indices(from_rset(as_rset(small, max_bytes = Inf))),
    analysis_indices(small)
  )
  for (invalid in list(NA_real_, -1, "unlimited", c(1, 2))) {
    expect_error(as_rset(small, max_bytes = invalid), "max_bytes")
  }
})

test_that("round trips preserve independent domains and temporal column metadata", {
  skip_if_not_installed("rsample")
  x <- tibble::tibble(
    id = 1:8, network = rep(c("north", "south"), each = 4),
    station = "a", time = rep(1:4, 2), other_time = rep(2:5, 2)
  )
  original <- neighbourhoods(x,
    temporal = temporal_window(time_span(0)),
    key = station, index = time, .by = network, id = id, within = "key"
  )
  resamples <- as_rset(original)
  restored <- from_rset(resamples)
  expect_equal(
    restored$ctx[c("key", "index", "by", "id", "domains")],
    original$ctx[c("key", "index", "by", "id", "domains")]
  )
  buffered <- buffer_folds(restored, temporal = temporal_window(time_span(1)), within = "key")
  expect_equal(omission_indices(buffered, 2), list(1:3))
  explicit <- from_rset(resamples, index = other_time, .by = station)
  expect_identical(explicit$ctx$index, "other_time")
  expect_identical(explicit$ctx$by, "station")
  expect_equal(explicit$ctx$domains, rep(1L, 8))
  attr(resamples, "naboer_columns") <- NULL
  external <- from_rset(resamples, id = id, key = station, index = time, .by = network)
  expect_equal(
    external$ctx[c("key", "index", "by", "id", "domains")],
    original$ctx[c("key", "index", "by", "id", "domains")]
  )
})

test_that("stale resample metadata and invalid recovered identifiers are rejected", {
  skip_if_not_installed("rsample")
  x <- tibble::tibble(id = 1:8, time = 1:8)
  resamples <- as_rset(loo_folds(x, id = id))
  missing <- resamples
  attr(missing, "naboer_columns")$by <- "removed_network"
  expect_error(from_rset(missing), "metadata.*missing columns", class = "naboer_error_alignment")
  invalid <- resamples
  attr(invalid, "naboer_columns") <- list(key = 1L)
  expect_error(from_rset(invalid), "metadata", class = "naboer_error_alignment")
  # Same modified source data in every split still must satisfy saved IDs.
  for (i in seq_along(resamples$splits)) resamples$splits[[i]]$data$id[[2]] <- 1L
  expect_error(from_rset(resamples), "uniquely identify", class = "naboer_error_alignment")
})
