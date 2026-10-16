# Internal compressed membership stores. Ranges index a shared permutation;
# packed sets use a single values vector and cumulative offsets.
point_index <- function(n) range_index(seq_len(n), seq_len(n))

range_index <- function(start, end, order = NULL) {
  structure(
    list(
      start = as.integer(start), end = as.integer(end), order = order,
      inverse = if (!is.null(order)) order(order)
    ),
    class = c("naboer_ranges", "naboer_index")
  )
}

list_index <- function(x) {
  if (inherits(x, "naboer_index")) {
    return(x)
  }
  if (!is.list(x)) abort_naboer("Membership must be a list of integer vectors.", "indices")
  # Avoid allocating a second list while testing whether ranges suffice.
  contiguous <- vapply(x, function(z) {
    !length(z) ||
      (is.numeric(z) && !anyNA(z) && all(diff(z) == 1))
  }, logical(1))
  if (all(contiguous)) {
    start <- vapply(x, function(z) if (length(z)) z[[1L]] else 1, numeric(1))
    end <- vapply(x, function(z) if (length(z)) z[[length(z)]] else 0, numeric(1))
    return(range_index(start, end))
  }
  lengths <- lengths(x)
  structure(
    list(
      values = as.integer(unlist(x, use.names = FALSE)),
      ends = cumsum(as.double(lengths))
    ),
    class = c("naboer_packed", "naboer_index")
  )
}

index_nfold <- function(x) {
  if (inherits(x, "naboer_ranges")) {
    length(x$start)
  } else if (inherits(x, "naboer_shared")) {
    length(x$lookup)
  } else if (inherits(x, "naboer_temporal_complement")) {
    length(x$start)
  } else if (inherits(x, "naboer_combined")) {
    index_nfold(x$indices[[1L]])
  } else if (inherits(x, "naboer_multiranges")) {
    length(x$ends)
  } else {
    length(x$ends)
  }
}

index_sizes <- function(x) {
  if (inherits(x, "naboer_ranges")) {
    pmax(0, as.double(x$end) - x$start + 1)
  } else if (inherits(x, "naboer_shared")) {
    index_sizes(x$members)[x$lookup]
  } else if (inherits(x, "naboer_temporal_complement")) {
    as.double(x$domain_end) - x$domain_start + 1 - pmax(0, x$end - x$start + 1)
  } else if (inherits(x, "naboer_combined")) {
    vapply(seq_len(index_nfold(x)), function(i) length(index_get(x, i)), numeric(1))
  } else if (inherits(x, "naboer_multiranges")) {
    x$sizes
  } else {
    diff(c(0, x$ends))
  }
}

index_get <- function(x, i) {
  if (inherits(x, "naboer_multiranges")) {
    first <- if (i == 1L) 1L else x$ends[[i - 1L]] + 1L
    last <- x$ends[[i]]
    if (last < first) {
      return(integer())
    }
    ranges <- seq.int(first, last)
    pos <- unlist(Map(seq.int, x$start[ranges], x$end[ranges]), use.names = FALSE)
    return(sort(as.integer(x$order[pos])))
  }
  if (inherits(x, "naboer_shared")) {
    return(index_get(x$members, x$lookup[[i]]))
  }
  if (inherits(x, "naboer_temporal_complement")) {
    left <- if (x$start[[i]] > x$domain_start[[i]]) seq.int(x$domain_start[[i]], x$start[[i]] - 1L) else integer()
    right <- if (x$end[[i]] < x$domain_end[[i]]) seq.int(x$end[[i]] + 1L, x$domain_end[[i]]) else integer()
    return(as.integer(x$order[c(left, right)]))
  }
  if (inherits(x, "naboer_combined")) {
    if (x$combine == "union") {
      return(sort(unique(unlist(lapply(x$indices, index_get, i = i), use.names = FALSE))))
    }
    sizes <- vapply(x$indices, function(z) index_size_one(z, i), numeric(1))
    first <- which.min(sizes)
    out <- index_get(x$indices[[first]], i)
    for (j in setdiff(seq_along(x$indices), first)) out <- out[index_contains(x$indices[[j]], i, out)]
    return(sort(out))
  }
  if (inherits(x, "naboer_ranges")) {
    if (x$end[[i]] < x$start[[i]]) {
      return(integer())
    }
    z <- seq.int(x$start[[i]], x$end[[i]])
    if (!is.null(x$order)) z <- x$order[z]
    return(as.integer(z))
  }
  start <- if (i == 1L) 1 else x$ends[[i - 1L]] + 1
  end <- x$ends[[i]]
  if (end < start) {
    return(integer())
  }
  x$values[seq.int(start, end)]
}

shared_index <- function(members, lookup) {
  structure(list(members = list_index(members), lookup = as.integer(lookup)),
    class = c("naboer_shared", "naboer_index")
  )
}

combined_index <- function(indices, combine) {
  structure(list(indices = indices, combine = combine), class = c("naboer_combined", "naboer_index"))
}

index_size_one <- function(x, i) {
  if (inherits(x, "naboer_multiranges")) {
    return(x$sizes[[i]])
  }
  if (inherits(x, "naboer_ranges")) {
    return(max(0, x$end[[i]] - x$start[[i]] + 1))
  }
  if (inherits(x, "naboer_shared")) {
    return(index_size_one(x$members, x$lookup[[i]]))
  }
  if (inherits(x, "naboer_packed")) {
    return(x$ends[[i]] - if (i == 1L) 0 else x$ends[[i - 1L]])
  }
  length(index_get(x, i))
}

index_contains <- function(x, i, candidates) {
  if (inherits(x, "naboer_ranges")) {
    position <- if (is.null(x$order)) candidates else x$inverse[candidates]
    return(position >= x$start[[i]] & position <= x$end[[i]])
  }
  if (inherits(x, "naboer_shared") && !is.null(x$group_members)) {
    groups <- index_get(x$group_members, x$fold_groups[[i]])
    return(x$observation_groups[candidates] %in% groups & x$observation_domains[candidates] == x$observation_domains[[i]])
  }
  candidates %in% index_get(x, i)
}

pack_index <- function(x) {
  sizes <- index_sizes(x)
  total <- sum(sizes)
  values <- integer(total)
  # Fill directly rather than expanding a million intermediate list elements.
  pos <- 0
  for (i in seq_len(index_nfold(x))) {
    k <- sizes[[i]]
    if (k) values[seq.int(pos + 1, pos + k)] <- sort(index_get(x, i))
    pos <- pos + k
  }
  list(values = values, ends = as.integer(cumsum(sizes)))
}

# Validate compressed relations directly, without expanding their memberships.
validate_compact_membership <- function(x) {
  a <- x$assessment
  o <- x$omission
  n <- x$ctx$n
  focal <- inherits(a, "naboer_ranges") && is.null(a$order) &&
    identical(a$start, seq_len(n)) && identical(a$end, a$start)
  if (focal && inherits(o, "naboer_shared")) {
    members <- o$members
    rows_by_group <- split(seq_len(n), o$lookup)
    inverse <- if (!is.null(members$order)) order(members$order) else NULL
    for (j in names(rows_by_group)) {
      i <- as.integer(j)
      rows <- rows_by_group[[j]]
      if (inherits(members, "naboer_multiranges")) {
        first <- if (i == 1L) 1L else members$ends[[i - 1L]] + 1L
        last <- members$ends[[i]]
        ranges <- seq.int(first, last)
        position <- inverse[rows]
        hit <- findInterval(position, members$start[ranges])
        if (any(hit == 0L) || any(position > members$end[ranges][hit])) abort_naboer("Assessment must be included in omission.", "indices")
      } else if (!all(rows %in% index_get(members, i))) abort_naboer("Assessment must be included in omission.", "indices")
    }
    return(TRUE)
  }
  if (focal && inherits(o, "naboer_ranges")) {
    position <- if (is.null(o$order)) seq_len(n) else order(o$order)
    if (any(position < o$start | position > o$end)) abort_naboer("Assessment must be included in omission.", "indices")
    return(TRUE)
  }
  if (focal && inherits(o, "naboer_multiranges")) {
    counts <- diff(c(0, o$ends))
    fold <- rep.int(seq_len(n), counts)
    position <- order(o$order)[fold]
    covers <- position >= o$start & position <= o$end
    if (any(tabulate(fold[covers], n) == 0L)) abort_naboer("Assessment must be included in omission.", "indices")
    return(TRUE)
  }
  if (inherits(a, "naboer_ranges") && inherits(o, "naboer_temporal_complement")) {
    if (!identical(a$order, o$order) ||
      any(a$start < o$domain_start | a$end > o$domain_end) ||
      any(a$start <= o$end & a$end >= o$start)) {
      abort_naboer("Forecast assessment must be omitted within its domain.", "indices")
    }
    return(TRUE)
  }
  FALSE
}

validate_index <- function(x, n) {
  fail <- function() abort_naboer("Invalid compact membership structure.", "indices")
  whole <- function(z, lower = 0, upper = .Machine$integer.max) {
    is.numeric(z) && !anyNA(z) && all(is.finite(z)) &&
      all(z == floor(z) & z >= lower & z <= upper)
  }
  permutation <- function(z) {
    if (!whole(z, 1, n) || length(z) != n || anyDuplicated(z)) fail()
  }
  endpoints <- function(ends, total) {
    if (!whole(ends) || any(diff(c(0, ends)) < 0) ||
      !length(ends) || tail(ends, 1) != total) {
      fail()
    }
  }
  ranges <- function(z) {
    if (!whole(z$start, 1, n) || !whole(z$end, 0, n) ||
      length(z$start) != length(z$end) || any(z$end < z$start - 1)) {
      fail()
    }
    if (!is.null(z$order)) permutation(z$order)
  }
  if (!inherits(x, "naboer_index")) fail()
  if (inherits(x, "naboer_ranges")) {
    ranges(x)
    if (!is.null(x$order) && !identical(x$inverse, order(x$order))) fail()
  } else if (inherits(x, "naboer_packed")) {
    if (!whole(x$values, 1, n)) fail()
    endpoints(x$ends, length(x$values))
    for (i in which(diff(c(0, x$ends)) > 1)) {
      if (anyDuplicated(index_get(x, i))) fail()
    }
  } else if (inherits(x, "naboer_shared")) {
    validate_index(x$members, n)
    if (!whole(x$lookup, 1, index_nfold(x$members))) fail()
  } else if (inherits(x, "naboer_multiranges")) {
    ranges(x)
    permutation(x$order)
    endpoints(x$ends, length(x$start))
    if (!whole(x$sizes) || length(x$sizes) != length(x$ends)) fail()
    counts <- diff(c(0, x$ends))
    fold <- rep.int(seq_along(x$ends), counts)
    adjacent <- which(diff(fold) == 0L)
    if (any(x$start[adjacent + 1L] <= x$end[adjacent])) fail()
    cumulative <- c(0, cumsum(as.double(x$end) - x$start + 1))
    if (!identical(as.double(x$sizes), diff(c(0, cumulative[x$ends + 1])))) fail()
  } else if (inherits(x, "naboer_temporal_complement")) {
    ranges(x)
    permutation(x$order)
    if (!whole(x$domain_start, 1, n) || !whole(x$domain_end, 1, n) ||
      length(x$domain_start) != length(x$start) || length(x$domain_end) != length(x$start) ||
      any(x$domain_start > x$start | x$domain_end < x$end)) {
      fail()
    }
  } else if (inherits(x, "naboer_combined")) {
    if (!is.list(x$indices) || length(x$indices) < 2L ||
      !identical(x$combine, "union") && !identical(x$combine, "intersection")) {
      fail()
    }
    for (z in x$indices) validate_index(z, n)
    if (length(unique(vapply(x$indices, index_nfold, integer(1)))) != 1L) fail()
  } else {
    fail()
  }
  invisible(x)
}
