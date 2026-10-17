test_that("nested groups and domains stay intact", {
  x <- tibble::tibble(
    region = rep(c("east", "west"), each = 4),
    site = rep(c("a", "a", "b", "b"), 2), y = 1:8
  )
  folds <- group_folds(x, c(region, site))
  expect_equal(assessment_indices(folds), list(1:2, 3:4, 5:6, 7:8))
  domains <- group_folds(x, site, .by = region)
  expect_equal(assessment_indices(domains), assessment_indices(folds))
  folds <- group_folds(x, c(region, site), v = 2, seed = 42)
  expect_equal(sort(lengths(assessment_indices(folds))), c(4L, 4L))
  bad <- x
  bad$site[[1]] <- NA_character_
  expect_error(group_folds(bad, site), "missing", class = "naboer_error")
})

test_that("environmental clustering scales predictors and preserves RNG state", {
  x <- tibble::tibble(
    a = c(0, 0.1, 0.2, 10, 10.1, 10.2),
    b = c(0, 1, 2, 100, 101, 102), constant = 2
  )
  set.seed(18)
  old_seed <- .Random.seed
  expect_message(
    folds <- environmental_folds(x, c(a, b, constant), v = 2, seed = 9),
    "constant"
  )
  expect_identical(.Random.seed, old_seed)
  expect_equal(sort(lengths(assessment_indices(folds))), c(3L, 3L))
  expect_equal(sort(assessment_indices(folds)[[1]]), 1:3)
  expect_equal(folds$spec$excluded[[1]], "constant")
  expect_equal(folds$spec$clusters, 2)
  raw <- environmental_folds(x, c(a, b), v = 2, standardise = FALSE, seed = 9)
  expect_equal(assessment_indices(raw), assessment_indices(folds))
  exact <- environmental_folds(x, a, v = nrow(x), seed = 1)
  expect_equal(assessment_indices(exact), as.list(seq_len(nrow(x))))
})

test_that("environmental clustering reports numerically unrepresentable scales", {
  x <- tibble::tibble(a = c(1, 2, 3, 4) * 1e-300)
  expect_error(environmental_folds(x, a, v = 2), "Rescale extreme")
})

test_that("phylogenetic radius agrees with full-tree distances after pruning", {
  skip_if_not_installed("ape")
  set.seed(731)
  tree <- ape::rtree(9)
  selected <- tree$tip.label[c(1, 3, 5, 7, 9)]
  x <- tibble::tibble(
    taxon = sample(rep(selected, c(3, 2, 4, 1, 3))),
    domain = rep(c("north", "south"), length.out = 13)
  )
  distances <- ape::cophenetic.phylo(tree)
  radius <- stats::median(distances[upper.tri(distances)])
  folds <- neighbourhoods(x,
    phylogenetic = phylogenetic_window(tree, taxon, radius),
    .by = domain
  )
  expected <- lapply(seq_len(nrow(x)), function(i) {
    unname(which(x$domain == x$domain[[i]] & distances[x$taxon[[i]], x$taxon] <= radius))
  })
  expect_equal(omission_indices(folds), expected)
})

test_that("invalid tree topology is rejected before ape algorithms or traversal", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- tibble::tibble(taxon = tree$tip.label)
  missing_edges <- tree
  missing_edges$edge <- NULL
  expect_error(phylogenetic_window(missing_edges, taxon, 1), "valid.*edges")
  duplicate_parent <- tree
  duplicate_parent$edge[2, 2] <- duplicate_parent$edge[1, 2]
  expect_error(phylogenetic_folds(x, duplicate_parent, taxon, v = 2), "parented edges")
  invalid_count <- tree
  invalid_count$Nnode <- NA_integer_
  expect_error(phylogenetic_window(invalid_count, taxon, 1), "internal-node count")
  # Retains all children once and root numbering, but disconnects a cycle of
  # two internal nodes from the root. Clade traversal must never enter it.
  cyclic <- tree
  ntip <- length(tree$tip.label)
  root <- ntip + 1L
  nodes <- setdiff(unique(tree$edge[, 1]), root)
  cyclic$edge[cyclic$edge[, 2] == nodes[[1]], 1] <- nodes[[2]]
  cyclic$edge[cyclic$edge[, 2] == nodes[[2]], 1] <- nodes[[1]]
  expect_error(phylogenetic_folds(x, cyclic, taxon,
    clades = list(1L, 2L, 3L, 4L)
  ), "connected and acyclic")
  reordered <- tree
  reordered$edge <- tree$edge[nrow(tree$edge):1, , drop = FALSE]
  reordered$edge.length <- tree$edge.length[nrow(tree$edge):1]
  attr(reordered, "order") <- NULL
  expect_equal(
    assessment_indices(phylogenetic_folds(x, reordered, taxon, v = 2)),
    assessment_indices(phylogenetic_folds(x, tree, taxon, v = 2))
  )
})

test_that("environmental strata and independent clustering are supported", {
  x <- tibble::tibble(
    a = rep(c(0, 0.1, 10, 10.1), 2),
    region = rep(c("a", "b"), each = 4),
    habitat = rep(c("low", "low", "high", "high"), 2)
  )
  folds <- environmental_folds(x, a, v = 2, .by = region, seed = 11)
  expect_equal(length(assessment_indices(folds)), 4L)
  expect_true(all(vapply(assessment_indices(folds), function(i) {
    length(unique(x$region[i])) == 1L
  }, logical(1))))
  strata <- environmental_folds(x, strata = habitat, .by = region)
  expect_equal(assessment_indices(strata), list(1:2, 3:4, 5:6, 7:8))
  expect_equal(length(assessment_indices(environmental_folds(
    x,
    strata = habitat, v = 2, .by = region, seed = 2
  ))), 2L)
})

test_that("environmental invalid inputs fail clearly", {
  x <- tibble::tibble(a = 1:4, b = letters[1:4], constant = 1, bad = c(1, 2, NA, 4))
  expect_error(environmental_folds(x), "exactly one")
  expect_error(environmental_folds(x, a, v = 2, strata = b), "exactly one")
  expect_error(environmental_folds(x, b, v = 2), "numeric")
  expect_error(environmental_folds(x, bad, v = 2), "missing")
  expect_error(environmental_folds(x, a, v = 2, standardise = NA), "TRUE or FALSE")
  expect_error(environmental_folds(x, a), "at least 2")
  expect_error(environmental_folds(x, a, v = 1), "at least 2")
  expect_error(environmental_folds(x, a, v = 5), "distinct")
  expect_message(
    expect_error(environmental_folds(x, constant, v = 2), "constant"),
    "constant"
  )
  x$a[[1]] <- Inf
  expect_error(environmental_folds(x, a, v = 2), "finite")
})

test_that("phylogenetic neighbourhoods reuse taxa and respect domains", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):2);")
  x <- tibble::tibble(
    taxon = c("a", "a", "b", "c", "d", "a"),
    domain = c(1, 1, 1, 1, 1, 2)
  )
  ctx <- prepare_data(x, by = rlang::quo(domain))
  members <- phylogenetic_membership(ctx, phylogenetic_window(tree, taxon, 2))
  expect_equal(index_get(members, 1), 1:3)
  expect_equal(index_get(members, 2), 1:3)
  expect_equal(index_get(members, 4), 4:5)
  expect_equal(index_get(members, 6), 6L)
  expect_s3_class(members, "naboer_shared")
  zero <- phylogenetic_membership(ctx, phylogenetic_window(tree, taxon, 0))
  expect_equal(index_get(zero, 1), 1:2)
  polytomy <- ape::read.tree(text = "(a:1,b:1,c:2,d:3);")
  poly <- phylogenetic_membership(ctx, phylogenetic_window(polytomy, taxon, 2))
  expect_equal(index_get(poly, 1), 1:3)
  expect_equal(index_get(poly, 4), 4L)
})

test_that("patristic clustering uses complete linkage and unique taxa", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):2);")
  x <- tibble::tibble(taxon = c("a", "a", "b", "c", "d"))
  folds <- phylogenetic_folds(x, tree, taxon, v = 2)
  expect_equal(assessment_indices(folds), list(1:3, 4:5))
  expect_equal(
    assessment_indices(phylogenetic_folds(x, tree, taxon, height = 2)),
    list(1:3, 4:5)
  )
  expect_equal(
    assessment_indices(phylogenetic_folds(x, tree, taxon, height = 0)),
    list(1:2, 3L, 4L, 5L)
  )
  unused <- ape::read.tree(text = "(((a:1,b:1):2,(c:1,d:1):2):1,unused:1);")
  expect_equal(
    assessment_indices(phylogenetic_folds(x, unused, taxon, v = 2)),
    list(1:3, 4:5)
  )
})

test_that("clade holdouts require complete nonoverlapping coverage", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):2,(c:1,d:1):2);")
  x <- tibble::tibble(taxon = c("a", "b", "c", "d"))
  folds <- phylogenetic_folds(x, tree, taxon, clades = list(c("a", "b"), c("c", "d")))
  expect_equal(assessment_indices(folds), list(1:2, 3:4))
  a_node <- ape::getMRCA(tree, c("a", "b"))
  c_node <- ape::getMRCA(tree, c("c", "d"))
  expect_equal(assessment_indices(phylogenetic_folds(
    x, tree, taxon,
    clades = list(a_node, c_node)
  )), list(1:2, 3:4))
  tips <- as.list(tree$tip.label)
  expect_equal(
    assessment_indices(phylogenetic_folds(x, tree, taxon, clades = tips)),
    as.list(1:4)
  )
  expect_error(phylogenetic_folds(x, tree, taxon, clades = list(c("a", "b"))), "cover")
  expect_error(phylogenetic_folds(x, tree, taxon,
    clades = list(c("a", "b"), "a", c("c", "d"))
  ), "overlap")
  expect_error(phylogenetic_folds(x, tree, taxon,
    clades = list(c("a", "c"), c("b", "d"))
  ), "complete clades")
  expect_error(phylogenetic_folds(x, ape::unroot(tree), taxon, clades = tips), "rooted")
  expect_error(phylogenetic_folds(x, tree, taxon, clades = list()), "nonempty")
  expect_error(phylogenetic_folds(x, tree, taxon, clades = list(100L)), "valid node")
})

test_that("phylogenetic validation and memory guards are explicit", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:0,b:0):2,c:1);")
  x <- tibble::tibble(taxon = c("a", "b", "c"))
  expect_equal(
    assessment_indices(phylogenetic_folds(x, tree, taxon, height = 0)),
    list(1:2, 3L)
  )
  expect_error(phylogenetic_window(tree, taxon, -1), "nonnegative")
  expect_error(phylogenetic_window(NULL, taxon, 1), "phylo")
  bad <- tree
  bad$edge.length <- NULL
  expect_error(phylogenetic_window(bad, taxon, 1), "branch")
  bad <- tree
  bad$edge.length[[1]] <- -1
  expect_error(phylogenetic_window(bad, taxon, 1), "nonnegative")
  bad <- tree
  bad$tip.label[[2]] <- "a"
  expect_error(phylogenetic_window(bad, taxon, 1), "unique")
  x$taxon[[1]] <- "unmatched"
  expect_error(phylogenetic_folds(x, tree, taxon, v = 2), "match")
  x$taxon[[1]] <- NA_character_
  expect_error(phylogenetic_folds(x, tree, taxon, v = 2), "missing")
  x$taxon[[1]] <- "a"
  expect_error(phylogenetic_folds(x, tree, taxon), "exactly one")
  expect_error(phylogenetic_folds(x, tree, taxon, v = 2, height = 1), "exactly one")
  expect_error(phylogenetic_folds(x, tree, taxon, v = 4), "number of observed")
  expect_error(phylogenetic_folds(x, tree, taxon, height = -1), "nonnegative")
  expect_error(phylogenetic_folds(x, tree, taxon, clades = list("a"), v = 2), "clades.*alone")
  x$taxon <- "a"
  expect_error(phylogenetic_folds(x, tree, taxon, height = 0), "at least two")
  ctx <- prepare_data(x)
  zero <- phylogenetic_membership(ctx, phylogenetic_window(tree, taxon, 0))
  expect_equal(index_get(zero, 1), 1:3)
  old <- options(naboer.max_bytes = 1)
  on.exit(options(old))
  expect_error(taxa_distances(tree, c("a", "b")), "exceed", class = "naboer_error_allocation")
})
