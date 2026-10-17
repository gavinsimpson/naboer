# Reproducible fixed-smoothing NCV approximation study. This is deliberately
# outside routine CRAN tests because it repeatedly refits complete models.
# Rscript benchmarks/approximation.R [n=120] [output.csv]
# Gaussian NCV should match full refits at fixed smoothing; Poisson differences
# measure the local approximation and may grow with omission width. The largest
# radius omits up to 80% of the data. This is a stress test, not a recommended
# neighbourhood size or a guarantee for other data, bases, or families.
if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(quiet = TRUE)
} else {
  library("naboer")
}
library("mgcv")
library("tibble")

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args)) as.integer(args[[1L]]) else 120L
stopifnot(!is.na(n), n >= 40L)
set.seed(20261017)
latent <- as.numeric(stats::arima.sim(model = list(ar = 0.65), n = n, sd = 0.25))
x <- seq(0, 1, length.out = n)
data <- tibble::tibble(
  time = seq_len(n), x = x,
  gaussian = sin(2 * pi * x) + latent,
  poisson = stats::rpois(n, exp(0.7 + sin(2 * pi * x) + latent))
)
widths <- sort(unique(c(0, 1, 3, 7, floor(n * c(0.05, 0.125, 0.25, 0.4)))))
smoothing <- 0.2
rows <- list()
for (response in c("gaussian", "poisson")) {
  family <- if (response == "gaussian") stats::gaussian() else stats::poisson()
  formula <- stats::as.formula(paste(response, "~ s(x, k = 10)"))
  setup <- mgcv::gam(formula,
    family = family, data = data, sp = smoothing,
    fit = FALSE, na.action = stats::na.fail
  )
  for (width in widths) {
    folds <- neighbourhoods(data, temporal = temporal_window(time_span(width)), index = time)
    timing <- system.time({
      approximate <- mgcv::gam(formula,
        family = family, data = data,
        method = "NCV", sp = smoothing,
        nei = as_nei(folds), na.action = stats::na.fail
      )
      approximate_q <- mgcv::gam(formula,
        family = family, data = data,
        method = "QNCV", sp = smoothing,
        nei = as_nei(folds), na.action = stats::na.fail
      )
      # Reuse the full-data setup and set omitted observation weights to zero.
      # Rebuilding smooths on training subsets changes basis/penalty scaling
      # and would compare a different model rather than the NCV approximation.
      exact <- 0
      exact_links <- numeric(n)
      converged <- TRUE
      for (i in seq_len(n)) {
        refit_setup <- setup
        refit_setup$w[omission_indices(folds, i)[[1L]]] <- 0
        refit <- mgcv::gam(G = refit_setup, method = "REML")
        converged <- converged && isTRUE(refit$converged)
        assess <- assessment_indices(folds, i)[[1L]]
        exact_links[assess] <- stats::predict(refit, newdata = data[assess, ], type = "link")
        prediction <- family$linkinv(exact_links[assess])
        exact <- exact + sum(family$dev.resids(data[[response]][assess], prediction, rep(1, length(assess))))
      }
    })
    # $NCV may be overwritten by mgcv's covariance calculations under its
    # default jackknife policy; $gcv.ubre retains the fitted criterion.
    score <- as.numeric(approximate$gcv.ubre)
    score_q <- as.numeric(approximate_q$gcv.ubre)
    approximate_links <- as.numeric(attr(approximate$gcv.ubre, "eta.cv"))
    stopifnot(
      converged, isTRUE(approximate$converged), isTRUE(approximate_q$converged),
      length(approximate_links) == n, all(is.finite(c(exact, score, score_q, approximate_links)))
    )
    relative <- abs(score - exact) / max(abs(exact), .Machine$double.eps)
    if (response == "gaussian") {
      stopifnot(relative < 1e-6, isTRUE(all.equal(score, score_q, tolerance = 1e-6)))
    }
    omission_sizes <- naboer:::index_sizes(folds$omission)
    rows[[length(rows) + 1L]] <- tibble::tibble(
      family = response, observations = n, radius = width, smoothing = smoothing,
      memberships = sum(omission_sizes),
      minimum_omitted_fraction = min(omission_sizes) / n,
      maximum_omitted_fraction = max(omission_sizes) / n,
      all_refits_converged = converged,
      exact = exact, NCV = score, QNCV = score_q,
      maximum_absolute_link_error = max(abs(approximate_links - exact_links)),
      NCV_signed_relative_error = (score - exact) / max(abs(exact), .Machine$double.eps),
      NCV_relative_error = relative,
      QNCV_relative_error = abs(score_q - exact) / max(abs(exact), .Machine$double.eps),
      elapsed_seconds = unname(timing[["elapsed"]])
    )
  }
}
result <- do.call(rbind, rows)
print(result, width = Inf)
print(list(seed = 20261017L, session = utils::sessionInfo(), hardware = Sys.info()))
if (length(args) > 1L) utils::write.csv(result, args[[2L]], row.names = FALSE)
