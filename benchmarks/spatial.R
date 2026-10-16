# Smaller moving-location, polygon, and dense-output scenarios.
# Rscript benchmarks/spatial.R
if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(quiet = TRUE)
} else {
  library("naboer")
}
library("sf")
library("units")
library("tibble")

run_spatial_case <- function(label, data, radius, temporal = NULL) {
  invisible(gc(reset = TRUE))
  elapsed <- system.time({
    folds <- neighbourhoods(data,
      spatial = spatial_window(radius),
      temporal = temporal, index = "time",
      .by = "domain",
      combine = if (!is.null(temporal)) "intersection" else NULL
    )
  })
  peak <- gc()
  focal <- unique(as.integer(c(1, nrow(data) %/% 2, nrow(data))))
  for (i in focal) {
    distances <- sf::st_distance(sf::st_geometry(data)[rep(i, nrow(data))],
      sf::st_geometry(data),
      by_element = TRUE
    )
    same <- data$domain == data$domain[[i]]
    if (!is.null(temporal)) same <- same & abs(data$time - data$time[[i]]) <= 2
    expected <- unname(which(distances <= radius & same))
    stopifnot(identical(omission_indices(folds, i)[[1L]], as.integer(expected)))
  }
  tibble::tibble(
    scenario = label, observations = nrow(data), elapsed_seconds = elapsed[["elapsed"]],
    R_heap_high_water_MiB = sum(peak[, ncol(peak)]),
    memberships = sum(naboer:::index_sizes(folds$omission)),
    object_bytes = as.numeric(object.size(folds))
  )
}

n <- 1000L
station <- rep(1:20, length.out = n)
time <- (seq_len(n) - 1L) %/% 20L
moving <- sf::st_as_sf(tibble::tibble(
  time = time, domain = 1L, x = station * 100 + sin(time), y = time * 0.1
), coords = c("x", "y"), crs = 3857)

polygon <- function(x, y) {
  sf::st_polygon(list(rbind(
    c(x, y), c(x + 1, y), c(x + 1, y + 1),
    c(x, y + 1), c(x, y)
  )))
}
polygons <- sf::st_sf(
  time = seq_len(500L), domain = 1L,
  geometry = sf::st_sfc(lapply(seq_len(500L), function(i) {
    polygon((i %% 25L) * 3, (i %/% 25L) * 3)
  }), crs = 3857)
)
dense <- sf::st_as_sf(tibble::tibble(
  time = seq_len(n), domain = rep(1:2, each = n / 2),
  x = rep(1:10, length.out = n), y = 0
), coords = c("x", "y"), crs = 3857)

results <- rbind(
  run_spatial_case(
    "moving locations", moving, units::set_units(1, "m"),
    temporal_window(time_span(2, "seconds"))
  ),
  run_spatial_case("polygon distances", polygons, units::set_units(2, "m")),
  run_spatial_case("dense spatial domains", dense, units::set_units(100, "m"))
)
print(Sys.info()[c("sysname", "release", "machine")])
print(R.version.string)
print(results, width = Inf)
