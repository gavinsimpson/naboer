# Run after installing naboer, or from the package root with pkgload available:
# Rscript benchmarks/spatiotemporal.R 1000000 --export
# Use /usr/bin/time -l on macOS (or /usr/bin/time -v on Linux) for process peak RSS.
if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(quiet = TRUE)
} else {
  library("naboer")
}
library("sf")
library("units")
library("tibble")

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args)) as.integer(args[[1L]]) else 1000000L
stopifnot(!is.na(n), n >= 100L)
stations <- min(1000L, floor(sqrt(n)))
station <- rep(seq_len(stations), length.out = n)
time <- (seq_len(n) - 1L) %/% stations
observations <- sf::st_as_sf(
  tibble::tibble(station = station, time = time, x = station * 100, y = 0),
  coords = c("x", "y"), crs = 3857
)
# One hundred metres spacing means at most three neighbouring stations. A
# two-second radius gives at most five observation times at each neighbour.
radius <- units::set_units(100, "m")
window <- temporal_window(time_span(2, "seconds"))
invisible(gc(reset = TRUE))
runtime <- system.time({
  folds <- neighbourhoods(observations,
    spatial = spatial_window(radius),
    temporal = window, key = station, index = time,
    combine = "intersection"
  )
})
peak <- gc()
counts <- naboer:::index_sizes(folds$omission)
# Independent direct reference over all rows for a few boundary/interior folds.
focal <- unique(as.integer(c(1L, stations, stations + 1L, n %/% 2L, n)))
for (i in focal) {
  expected <- which(abs(station - station[[i]]) <= 1L &
    abs(time - time[[i]]) <= 2)
  stopifnot(identical(omission_indices(folds, i)[[1L]], as.integer(expected)))
}
report <- list(
  hardware = Sys.info()[c("sysname", "release", "version", "machine")],
  R = R.version.string,
  package_versions = vapply(
    c("naboer", "sf", "units"),
    function(p) as.character(utils::packageVersion(p)), character(1)
  ),
  observations = n, stations = stations, runtime_seconds = unclass(runtime),
  R_heap_high_water_MiB = sum(peak[, ncol(peak)]),
  # gc() column names differ across R releases; print the complete matrix too.
  gc = peak,
  memberships = sum(counts),
  compact_object_bytes = as.numeric(object.size(folds)),
  estimated_nei_integer_bytes = 4 * (sum(counts) + n + 2 * n),
  checked_focal_rows = focal
)
if ("--export" %in% args) {
  invisible(gc(reset = TRUE))
  export_time <- system.time(nei <- as_nei(folds))
  export_peak <- gc()
  for (i in focal) {
    expected <- which(abs(station - station[[i]]) <= 1L &
      abs(time - time[[i]]) <= 2)
    begin <- if (i == 1L) 1L else nei$ma[[i - 1L]] + 1L
    actual <- nei$a[seq.int(begin, nei$ma[[i]])]
    stopifnot(identical(actual, as.integer(expected)))
    stopifnot(identical(nei$d[[nei$md[[i]]]], i))
  }
  report$export_runtime_seconds <- unclass(export_time)
  report$nei_object_bytes <- as.numeric(object.size(nei))
  report$export_R_heap_high_water_MiB <- sum(export_peak[, ncol(export_peak)])
  report$export_gc <- export_peak
}
print(report)
