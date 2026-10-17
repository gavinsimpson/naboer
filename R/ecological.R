#' Hold out observational groups
#'
#' Groups are identified by the combination of the selected columns. Select
#' all relevant levels, for example `group = c(region, site)`, for nested groups.
#' Groups never cross independent domains defined by `.by`.
#'
#' @param data An observation data frame.
#' @param group Columns defining whole groups, using tidy selection.
#' @param v Number of folds. `NULL` leaves each group out in turn.
#' @param balance Balance observation counts when allocating groups to folds.
#' @param seed Optional local random seed.
#' @param .by Columns defining independent domains.
#' @param id Optional unique observation identifier column.
#' @return A `naboer_folds` object.
#' @export
group_folds <- function(data, group, v = NULL, balance = TRUE, seed = NULL,
                        .by = NULL, id = NULL) {
  ctx <- prepare_data(data, by = rlang::enquo(.by), id = rlang::enquo(id))
  cols <- select_cols(ctx$data, rlang::enquo(group), "group", allow_empty = FALSE)
  check_ecological_columns(ctx$data, cols, "group")
  block <- group_ids(ctx$data, cols)
  group_folds_from_blocks(ctx, block, v, balance, seed,
    spec = list(strategy = "group", group = cols)
  )
}

#' Construct environmental cross-validation folds
#'
#' Select numeric predictors for k-means clustering, or supply existing strata.
#' Clustering is performed separately in each independent domain. Standardising
#' predictors gives each nonconstant variable equal scale. Constant predictors
#' are reported and excluded within each domain.
#' For clustering, `v` is the number of clusters per domain, and each cluster
#' forms one holdout. Clusters can have unequal sizes. `balance` only affects
#' allocation of preexisting `strata`, not k-means clustering.
#'
#' @inheritParams group_folds
#' @param predictors Numeric predictor columns, using tidy selection.
#' @param strata Existing environmental strata columns, using tidy selection.
#' @param standardise Standardise numeric predictors before clustering.
#' @param v Number of k-means clusters, or folds of whole existing strata.
#' @export
environmental_folds <- function(data, predictors = NULL, v = NULL,
                                strata = NULL, standardise = TRUE,
                                balance = TRUE, seed = NULL,
                                .by = NULL, id = NULL) {
  ctx <- prepare_data(data, by = rlang::enquo(.by), id = rlang::enquo(id))
  predictors <- select_cols(ctx$data, rlang::enquo(predictors), "predictors")
  strata <- select_cols(ctx$data, rlang::enquo(strata), "strata")
  if ((length(predictors) == 0L) == (length(strata) == 0L)) {
    abort_naboer("Supply exactly one of `predictors` and `strata`.")
  }
  if (length(strata)) {
    check_ecological_columns(ctx$data, strata, "strata")
    return(group_folds_from_blocks(
      ctx, group_ids(ctx$data, strata), v, balance, seed,
      spec = list(strategy = "environmental", strata = strata)
    ))
  }
  if (!is.logical(standardise) || length(standardise) != 1L || is.na(standardise)) {
    abort_naboer("`standardise` must be TRUE or FALSE.")
  }
  check_ecological_columns(ctx$data, predictors, "predictors")
  if (!all(vapply(ctx$data[predictors], is.numeric, logical(1)))) {
    abort_naboer("Environmental predictors must all be numeric.")
  }
  values <- as.matrix(ctx$data[predictors])
  if (any(!is.finite(values))) {
    abort_naboer("Environmental predictors must contain only finite values.")
  }
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || !is.finite(v) ||
    v < 2 || v > .Machine$integer.max || v != floor(v)) {
    abort_naboer("For predictor clustering, `v` must be an integer of at least 2.")
  }
  block <- integer(ctx$n)
  excluded <- list()
  with_seed(seed, {
    for (domain in unique(ctx$domains)) {
      rows <- which(ctx$domains == domain)
      x <- values[rows, , drop = FALSE]
      keep <- vapply(seq_len(ncol(x)), function(j) {
        length(unique(x[, j])) > 1L
      }, logical(1))
      excluded[[as.character(domain)]] <- predictors[!keep]
      if (any(!keep)) {
        rlang::inform(paste0(
          "Excluding constant environmental predictors in domain ", domain,
          ": ", paste(predictors[!keep], collapse = ", "), "."
        ))
      }
      if (!any(keep)) {
        abort_naboer("Every environmental predictor is constant in a domain.")
      }
      x <- x[, keep, drop = FALSE]
      if (standardise) x <- scale(x)
      if (any(!is.finite(x))) {
        abort_naboer(paste0(
          "Standardising environmental predictors produced nonfinite values. ",
          "Rescale extreme predictor values before clustering."
        ))
      }
      distinct <- nrow(unique(x))
      if (distinct < v) {
        abort_naboer("Each domain must have at least `v` distinct predictor rows.")
      }
      # Hartigan-Wong requires fewer centres than rows. When every distinct
      # predictor vector is a centre, its assignment is exact without fitting.
      block[rows] <- if (distinct == v) {
        group_ids(tibble::as_tibble(x), colnames(x))
      } else {
        stats::kmeans(x, centers = v, nstart = 10L)$cluster
      }
    }
  })
  group_folds_from_blocks(
    ctx, block,
    spec = list(
      strategy = "environmental", predictors = predictors,
      standardise = standardise, excluded = excluded, clusters = v, seed = seed
    )
  )
}

#' Specify a phylogenetic neighbourhood
#'
#' Distances use branch lengths between observed taxa, not distances between
#' replicated observation rows. Unused tree tips are allowed.
#'
#' @param tree An `ape` phylogenetic tree with finite, nonnegative branch lengths.
#' @param taxon Column containing tree tip labels, using tidy selection.
#' @param radius Inclusive patristic distance threshold.
#' @return A phylogenetic window specification for [neighbourhoods()].
#' @export
phylogenetic_window <- function(tree, taxon, radius) {
  validate_tree(tree)
  if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) ||
    !is.finite(radius) || radius < 0) {
    abort_naboer("`radius` must be one finite nonnegative branch-length distance.")
  }
  structure(list(tree = tree, taxon = rlang::enquo(taxon), radius = radius),
    class = "naboer_phylogenetic_window"
  )
}

#' Construct phylogenetic cross-validation folds
#'
#' Supply `v` or `height` to cut complete-linkage clustering of patristic
#' distances. Alternatively, provide explicitly selected nonoverlapping clades
#' of a rooted tree. A clade is a tree node number or a character vector of all
#' descendant tip labels. Selected clades must cover all observed taxa.
#'
#' @inheritParams group_folds
#' @inheritParams phylogenetic_window
#' @param v Number of patristic-distance clusters.
#' @param height Maximum complete-linkage cluster height.
#' @param clades A list of nodes or complete sets of descendant tip labels.
#' @export
phylogenetic_folds <- function(data, tree, taxon, v = NULL, height = NULL,
                               clades = NULL, .by = NULL, id = NULL) {
  ctx <- prepare_data(data, by = rlang::enquo(.by), id = rlang::enquo(id))
  taxa <- observed_taxa(ctx, tree, rlang::enquo(taxon))
  observed <- unique(taxa)
  if (!is.null(clades)) {
    if (!is.null(v) || !is.null(height)) {
      abort_naboer("Use `clades` alone, or choose exactly one of `v` and `height`.")
    }
    if (!ape::is.rooted(tree)) {
      abort_naboer("Clade holdouts require a rooted tree.")
    }
    if (!is.list(clades) || !length(clades)) {
      abort_naboer("`clades` must be a nonempty list of nodes or tip-label vectors.")
    }
    tips <- lapply(clades, resolve_clade, tree = tree)
    joined <- unlist(tips, use.names = FALSE)
    if (anyDuplicated(joined)) abort_naboer("Selected clades must not overlap.")
    if (!all(observed %in% joined)) {
      abort_naboer("Selected clades must cover every observed taxon.")
    }
    assignments <- rep(seq_along(tips), lengths(tips))
    block <- assignments[match(taxa, joined)]
  } else {
    if (is.null(v) == is.null(height)) {
      abort_naboer("Choose exactly one of `v` and `height`, or supply `clades`.")
    }
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || is.na(v) ||
      !is.finite(v) || v < 2 || v > length(observed) || v != floor(v))) {
      abort_naboer("`v` must be an integer between 2 and the number of observed taxa.")
    }
    if (!is.null(height) && (!is.numeric(height) || length(height) != 1L ||
      is.na(height) || !is.finite(height) || height < 0)) {
      abort_naboer("`height` must be one finite nonnegative distance.")
    }
    if (length(observed) < 2L) {
      abort_naboer("Phylogenetic clustering requires at least two observed taxa.")
    }
    distances <- taxa_distances(tree, observed)
    clusters <- stats::hclust(stats::as.dist(distances), method = "complete")
    assignments <- stats::cutree(clusters, k = v, h = height)
    block <- unname(assignments[match(taxa, observed)])
  }
  group_folds_from_blocks(
    ctx, block,
    spec = list(
      strategy = "phylogenetic", taxon = rlang::as_label(rlang::enquo(taxon)),
      v = v, height = height, clades = clades
    )
  )
}

phylogenetic_membership <- function(ctx, window) {
  if (!inherits(window, "naboer_phylogenetic_window")) {
    abort_naboer("`phylogenetic` must be created with phylogenetic_window().")
  }
  taxa <- observed_taxa(ctx, window$tree, window$taxon)
  observed <- unique(taxa)
  distances <- taxa_distances(window$tree, observed)
  taxa_id <- match(taxa, observed)
  # Search among taxa once, then reuse each taxon's set for replicated rows.
  lookup <- group_ids(
    tibble::tibble(domain = ctx$domains, taxon = taxa_id),
    c("domain", "taxon")
  )
  membership <- vector("list", max(lookup))
  for (domain in unique(ctx$domains)) {
    domain_rows <- which(ctx$domains == domain)
    by_taxon <- split(domain_rows, taxa_id[domain_rows])
    for (j in names(by_taxon)) {
      near <- which(distances[as.integer(j), ] <= window$radius)
      neighbours <- sort.int(unlist(
        by_taxon[intersect(
          names(by_taxon),
          as.character(near)
        )],
        use.names = FALSE
      ))
      membership[[lookup[[by_taxon[[j]][[1L]]]]]] <- neighbours
    }
  }
  shared_index(list_index(membership), lookup)
}

check_ecological_columns <- function(data, cols, arg) {
  if (any(vapply(data[cols], anyNA, logical(1)))) {
    abort_naboer(paste0("Columns selected by `", arg, "` must not contain missing values."))
  }
}

validate_tree <- function(tree) {
  need_package("ape")
  if (!inherits(tree, "phylo") || !is.character(tree$tip.label) || !length(tree$tip.label) ||
    anyNA(tree$tip.label) || anyDuplicated(tree$tip.label)) {
    abort_naboer("`tree` must be a phylo object with unique nonmissing tip labels.")
  }
  # ape algorithms rely on a connected, acyclic tree numbered with tips first
  # and root n_tip + 1. Validate topology before calling them or traversing
  # descendant clades: malformed cycles can otherwise hang or mislabel taxa.
  n_tip <- length(tree$tip.label)
  if (!is.numeric(tree$Nnode) || length(tree$Nnode) != 1L ||
    is.na(tree$Nnode) || !is.finite(tree$Nnode) || tree$Nnode < 1 ||
    tree$Nnode != floor(tree$Nnode) || tree$Nnode > .Machine$integer.max - n_tip) {
    abort_naboer("The tree must have a valid internal-node count.")
  }
  total <- n_tip + tree$Nnode
  edge <- tree$edge
  if (!is.matrix(edge) || ncol(edge) != 2L || !is.numeric(edge) ||
    anyNA(edge) || any(!is.finite(edge)) || any(edge != floor(edge)) ||
    any(edge < 1 | edge > total) || nrow(edge) != total - 1L ||
    any(edge[, 1L] <= n_tip) || anyDuplicated(edge[, 2L]) ||
    !setequal(edge[, 2L], setdiff(seq_len(total), n_tip + 1L))) {
    abort_naboer("The tree must have valid, uniquely parented edges and root n_tip + 1.")
  }
  parents <- integer(total)
  parents[edge[, 2L]] <- edge[, 1L]
  state <- integer(total)
  state[[n_tip + 1L]] <- 2L
  for (start in seq_len(total)) {
    node <- start
    while (state[[node]] == 0L) {
      state[[node]] <- 1L
      node <- parents[[node]]
    }
    if (state[[node]] == 1L) {
      abort_naboer("Tree edges must be connected and acyclic.")
    }
    node <- start
    while (state[[node]] == 1L) {
      state[[node]] <- 2L
      node <- parents[[node]]
    }
  }
  if (!is.numeric(tree$edge.length) || length(tree$edge.length) != nrow(edge) ||
    any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)) {
    abort_naboer("The tree must have finite, nonnegative lengths for every branch.")
  }
  invisible(tree)
}

observed_taxa <- function(ctx, tree, taxon) {
  validate_tree(tree)
  cols <- select_cols(ctx$data, taxon, "taxon", allow_empty = FALSE)
  if (length(cols) != 1L) abort_naboer("Select exactly one `taxon` column.")
  check_ecological_columns(ctx$data, cols, "taxon")
  taxa <- as.character(ctx$data[[cols]])
  if (any(!taxa %in% tree$tip.label)) {
    abort_naboer("Every observed taxon must match a tree tip label.")
  }
  taxa
}

taxa_distances <- function(tree, taxa) {
  # Prune before computing the dense taxon matrix: unused tips can be numerous.
  n <- length(taxa)
  if (8 * as.double(n)^2 > getOption("naboer.max_bytes", 1024^3)) {
    abort_naboer("Taxon distances exceed `getOption('naboer.max_bytes')`; reduce the tree or raise the limit.",
      class = "allocation"
    )
  }
  if (n == 1L) {
    return(matrix(0, 1L, 1L, dimnames = list(taxa, taxa)))
  }
  pruned <- ape::keep.tip(tree, taxa)
  ape::cophenetic.phylo(pruned)[taxa, taxa, drop = FALSE]
}

descendant_tips <- function(tree, node) {
  ntip <- length(tree$tip.label)
  if (node <= ntip) {
    return(tree$tip.label[[node]])
  }
  frontier <- node
  found <- integer()
  while (length(frontier)) {
    children <- tree$edge[tree$edge[, 1L] %in% frontier, 2L]
    found <- c(found, children[children <= ntip])
    frontier <- children[children > ntip]
  }
  tree$tip.label[found]
}

resolve_clade <- function(clade, tree) {
  if (is.numeric(clade) && length(clade) == 1L && !is.na(clade) &&
    is.finite(clade) && clade == floor(clade) && clade >= 1L &&
    clade <= length(tree$tip.label) + tree$Nnode) {
    return(descendant_tips(tree, clade))
  }
  if (!is.character(clade) || !length(clade) || anyNA(clade) ||
    anyDuplicated(clade) || any(!clade %in% tree$tip.label)) {
    abort_naboer("Each clade must be a valid node or a nonempty set of unique tree tips.")
  }
  if (length(clade) == 1L) {
    return(clade)
  }
  node <- ape::getMRCA(tree, clade)
  if (is.null(node) || !setequal(descendant_tips(tree, node), clade)) {
    abort_naboer("Tip-label sets supplied in `clades` must be complete clades.")
  }
  clade
}
