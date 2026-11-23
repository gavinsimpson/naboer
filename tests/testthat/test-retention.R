test_that("folds retain membership and column names without storing observations", {
  data <- tibble::tibble(id = letters[1:6], time = c(1L, 5L, 3L, 2L, 6L, 4L), y = 1:6)
  window <- temporal_window(time_span(1))
  folds <- neighbourhoods(data, index = time, id = id, temporal = window)
  other <- neighbourhoods(
    transform(data, y = rev(y)), index = time, id = id, temporal = window
  )

  expect_named(folds$ctx, c("n", "key", "index", "by", "id"), ignore.order = TRUE)
  expect_identical(folds, other)
  expect_identical(fold_data(folds, data), data)
  expect_error(fold_data(folds), class = "naboer_error_data")
  expect_error(buffer_folds(folds, temporal = window), class = "naboer_error_data")
})

test_that("rsample export requires observations", {
  skip_if_not_installed("rsample")
  expect_error(as_rset(loo_folds(tibble::tibble(y = 1:5))),
    class = "naboer_error_data"
  )
})

test_that("phylogenetic fold specifications do not capture caller observations", {
  skip_if_not_installed("ape")
  construct <- function(y) {
    data <- tibble::tibble(taxon = letters[1:4], y = y)
    tree <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):2);")
    window <- phylogenetic_window(tree, taxon, radius = 2)
    list(
      neighbourhoods(data, phylogenetic = window),
      buffer_folds(loo_folds(data), data, phylogenetic = window)
    )
  }

  expect_identical(serialize(construct(1:4), NULL), serialize(construct(5:8), NULL))
})

test_that("mgcv fits with data-free neighbourhoods and original model data", {
  time <- seq_len(48L)
  data <- tibble::tibble(time = time, y = sin(time / 5) + 0.1 * cos(time * 2))
  folds <- neighbourhoods(data, index = time, temporal = temporal_window(time_span(1)))
  fit <- mgcv::gam(y ~ s(time, k = 6), data = data, method = "NCV",
    nei = as_nei(folds, data = data), na.action = stats::na.fail
  )

  expect_identical(fit$method, "NCV")
  expect_equal(stats::nobs(fit), nrow(data))
  expect_true(all(is.finite(fit$sp)))
})
