# Rscript benchmarks/temporal.R 1000000
# Run with /usr/bin/time -l (macOS) or /usr/bin/time -v (Linux) for peak RSS.
# The gc report measures R heap high-water marks, not operating-system RSS.
if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(quiet = TRUE)
} else {
  library("naboer")
}
library("tibble")

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args)) as.integer(args[[1L]]) else 1000000L
stopifnot(!is.na(n), n >= 100L)
series <- min(1000L, n %/% 20L)
# Deliberately unequal lengths, deterministic irregular spacing and mixed rows.
lengths <- seq_len(series) %% 37L + 20L
lengths <- as.integer(floor(lengths / sum(lengths) * n))
lengths[1L] <- lengths[1L] + n - sum(lengths)
site <- rep(seq_len(series), lengths)
position <- sequence(lengths)
time <- position * 60 + (position %% 5L) * 3
row_order <- order(position, site)
observations <- tibble::tibble(site = site[row_order], time = time[row_order])

measure <- function(label, width) {
  invisible(gc(reset = TRUE))
  timing <- system.time({
    folds <- neighbourhoods(observations,
      temporal = temporal_window(time_span(width, "seconds")),
      key = site, index = time, within = "key"
    )
  })
  heap <- gc()
  focal <- unique(as.integer(c(1L, n %/% 3L, n %/% 2L, n)))
  for (i in focal) {
    expected <- which(observations$site == observations$site[i] &
      abs(observations$time - observations$time[i]) <= width)
    stopifnot(identical(omission_indices(folds, i)[[1L]], as.integer(expected)))
  }
  counts <- naboer:::index_sizes(folds$omission)
  estimated <- 4 * (sum(counts) + 3 * n)
  guard <- if (estimated > 1024^3) {
    result <- tryCatch(as_nei(folds, max_bytes = 1024^3),
      naboer_error_memory = function(e) class(e)[1L],
      naboer_error_overflow = function(e) class(e)[1L]
    )
    stopifnot(is.character(result), length(result) == 1L)
    result
  } else {
    "below default materialisation limit"
  }
  list(
    design = label, radius_seconds = width,
    runtime_seconds = unclass(timing), gc = heap,
    R_heap_high_water_MiB = sum(heap[, ncol(heap)]),
    memberships = sum(counts), compact_bytes = as.numeric(object.size(folds)),
    packed_nei_integer_bytes = estimated, materialisation_guard = guard,
    checked_focal_rows = focal
  )
}

report <- list(
  hardware = Sys.info()[c("sysname", "release", "version", "machine")],
  R = R.version.string,
  naboer = as.character(getNamespaceVersion("naboer")),
  observations = n, series = series,
  sparse = measure("sparse irregular multiple series", 120),
  dense = measure("entire-series omissions retained as ranges", max(time))
)
print(report)
