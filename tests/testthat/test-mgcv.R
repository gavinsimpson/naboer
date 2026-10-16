test_that("Gaussian NCV equals independent fixed-penalty leave-out refits", {
  skip_if_not_installed("mgcv", minimum_version = "1.9.4")
  n <- 40L
  x <- data.frame(x = seq(0, 1, length.out = n))
  x$y <- sin(2 * pi * x$x) + 0.1 * cos(seq_len(n) * 2)
  omitted <- lapply(seq_len(n), function(i) seq.int(max(1L, i - 1L), min(n, i + 1L)))
  folds <- custom_folds(x, as.list(seq_len(n)), omitted)
  nei <- as_nei(folds)
  # Reusing the full basis keeps penalty normalisation and centring fixed.
  # Rebuilding a separate basis on each training sample changes the problem.
  setup <- mgcv::gam(y ~ s(x, k = 6, bs = "cr"), data = x, fit = FALSE)
  sp <- 0.7
  penalty <- matrix(0, ncol(setup$X), ncol(setup$X))
  j <- seq.int(setup$off[[1]], length.out = nrow(setup$S[[1]]))
  penalty[j, j] <- sp * setup$S[[1]]
  predictions <- vapply(seq_len(n), function(i) {
    training <- setdiff(seq_len(n), omitted[[i]])
    design <- setup$X[training, , drop = FALSE]
    beta <- solve(crossprod(design) + penalty, crossprod(design, setup$y[training]))
    drop(setup$X[i, , drop = FALSE] %*% beta)
  }, numeric(1))
  fit <- mgcv::gam(
    G = setup, sp = sp, method = "NCV", nei = nei,
    na.action = stats::na.fail
  )
  # mgcv's optimised criterion is a sum of squared prediction errors for a
  # Gaussian response; gcv.ubre holds it, including the fold prediction values.
  expect_equal(as.numeric(fit$gcv.ubre), sum((setup$y - predictions)^2), tolerance = 1e-6)
  expect_equal(as.numeric(attr(fit$gcv.ubre, "eta.cv")), predictions, tolerance = 1e-6)
  quadratic <- mgcv::gam(
    G = setup, sp = sp, method = "QNCV", nei = nei,
    na.action = stats::na.fail
  )
  expect_equal(as.numeric(quadratic$gcv.ubre), as.numeric(fit$gcv.ubre), tolerance = 1e-6)
  expect_equal(stats::coef(quadratic), stats::coef(fit), tolerance = 1e-6)
})

test_that("gam optimises NCV and QNCV for Gaussian and Poisson data", {
  skip_if_not_installed("mgcv", minimum_version = "1.9.4")
  n <- 80L
  x <- data.frame(x = seq(0, 1, length.out = n))
  x$y <- sin(2 * pi * x$x) + 0.2 * cos(seq_len(n) * 1.7)
  folds <- custom_folds(x, as.list(seq_len(n)), lapply(seq_len(n), function(i) {
    seq.int(max(1L, i - 1L), min(n, i + 1L))
  }))
  nei <- as_nei(folds)
  normal <- mgcv::gam(y ~ s(x, k = 6), data = x, method = "NCV", nei = nei)
  quadratic <- mgcv::gam(y ~ s(x, k = 6), data = x, method = "QNCV", nei = nei)
  expect_equal(normal$sp, quadratic$sp, tolerance = 1e-6)
  expect_equal(as.numeric(normal$gcv.ubre), as.numeric(quadratic$gcv.ubre), tolerance = 1e-6)
  x$y <- as.integer(round(exp(1 + sin(2 * pi * x$x))))
  for (method in c("NCV", "QNCV")) {
    fit <- mgcv::gam(y ~ s(x, k = 6),
      data = x, family = stats::poisson(),
      method = method, nei = nei, na.action = stats::na.fail
    )
    expect_true(is.finite(fit$gcv.ubre))
    expect_true(all(is.finite(stats::coef(fit))))
    expect_true(fit$converged)
  }
})

test_that("sampled bam NCV receives complete neighbourhoods for covariance", {
  skip_if_not_installed("mgcv", minimum_version = "1.9.4")
  n <- 80L
  x <- data.frame(x = seq(0, 1, length.out = n))
  x$y <- sin(2 * pi * x$x) + 0.2 * cos(seq_len(n) * 1.7)
  folds <- custom_folds(x, as.list(seq_len(n)), lapply(seq_len(n), function(i) {
    seq.int(max(1L, i - 1L), min(n, i + 1L))
  }))
  sample <- seq.int(1L, n, by = 3L)
  full <- as_nei(folds)
  nei <- as_nei(folds, sample = sample)
  expect_identical(nei[c("a", "ma", "d", "md")], full)
  expect_length(nei$ma, n)
  fit <- mgcv::bam(y ~ s(x, k = 8),
    data = x, method = "NCV", nei = nei,
    discrete = TRUE, nthreads = 1L, na.action = stats::na.fail
  )
  expect_identical(fit$method, "NCV")
  expect_true(is.finite(fit$gcv.ubre))
  expect_true(all(is.finite(fit$Vp)))
  expect_length(stats::fitted(fit), n)
})
