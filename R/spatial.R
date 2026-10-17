#' Specify a spatial neighbourhood
#'
#' Distances are measured between the original geometries, including polygons.
#' The current sf spherical/ellipsoidal setting is respected and never changed.
#' Distances use stored coordinates without implicit precision rounding. To use
#' rounded coordinates, transform them explicitly before constructing folds.
#' @param radius A nonnegative scalar with physical distance units, for example
#'   `units::set_units(5, "km")`.
#' @return A spatial window specification.
#' @export
spatial_window <- function(radius) {
  need_package("units")
  if (!inherits(radius, "units") || length(radius) != 1L ||
    !is.finite(as.numeric(radius)) || as.numeric(radius) < 0) {
    abort_naboer("`radius` must be one finite, nonnegative distance with units.")
  }
  tryCatch(units::set_units(radius, "m"), error = function(e) {
    abort_naboer("`radius` must have physical distance units.")
  })
  structure(list(radius = radius), class = "naboer_spatial_window")
}

spatial_geometry <- function(ctx) {
  need_package("sf")
  cols <- names(ctx$data)[vapply(ctx$data, inherits, logical(1), "sfc")]
  active <- attr(ctx$data, "sf_column", exact = TRUE)
  if (length(active) == 1L && active %in% cols) cols <- active
  if (length(cols) != 1L) {
    abort_naboer("Spatial operations require one active sf geometry column.", "geometry")
  }
  geom <- ctx$data[[cols]]
  # GEOS predicates may apply sf precision rounding whereas st_distance() uses
  # stored coordinates. Use one coordinate interpretation throughout queries;
  # the canonical data retain their original geometry and precision metadata.
  geom <- sf::st_set_precision(geom, 0)
  if (is.na(sf::st_crs(geom))) {
    abort_naboer("Spatial operations require a declared CRS; set it with sf::st_set_crs().", "geometry")
  }
  if (any(sf::st_is_empty(geom)) || anyNA(sf::st_is_valid(geom)) ||
    !all(sf::st_is_valid(geom))) {
    abort_naboer("Geometry must be nonempty and valid; repair it explicitly before constructing folds.", "geometry")
  }
  geom
}

spatial_engine <- function(geom) {
  if (!sf::st_is_longlat(geom)) "GEOS" else if (sf::sf_use_s2()) "s2" else "ellipsoidal"
}

# Query unique geometries first, then expand to observation rows only when needed.
spatial_membership <- function(ctx, window) {
  if (!inherits(window, "naboer_spatial_window")) {
    abort_naboer("`spatial` must be created with spatial_window().")
  }
  geom <- spatial_geometry(ctx)
  unique_geom <- geom[!duplicated(geom)]
  lookup <- match(geom, unique_geom)
  membership <- vector("list", length(unique_geom))
  batch_size <- getOption("naboer.spatial_batch_size", 1000L)
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
    !is.finite(batch_size) || batch_size < 1 ||
    batch_size != floor(batch_size) || batch_size > .Machine$integer.max) {
    abort_naboer("`naboer.spatial_batch_size` must be a positive integer within the R index range.")
  }
  if (!sf::st_is_longlat(geom)) {
    radius <- as.numeric(units::set_units(window$radius, sf::st_crs(geom)$ud_unit,
      mode = "standard"
    ))
    boxes <- if (radius == 0) {
      unique_geom
    } else {
      lapply(unique_geom, function(g) {
        b <- sf::st_bbox(g)
        xx <- c(b[[1L]] - radius, b[[3L]] + radius)
        yy <- c(b[[2L]] - radius, b[[4L]] + radius)
        sf::st_polygon(list(cbind(xx[c(1, 2, 2, 1, 1)], yy[c(1, 1, 2, 2, 1)])))
      })
    }
    if (!inherits(boxes, "sfc")) boxes <- sf::st_sfc(boxes, crs = sf::st_crs(geom))
  }
  for (start in seq.int(1L, length(unique_geom), by = as.integer(batch_size))) {
    rows <- seq.int(start, min(length(unique_geom), start + batch_size - 1L))
    if (sf::st_is_longlat(geom)) {
      hits <- sf::st_is_within_distance(unique_geom[rows], unique_geom,
        dist = window$radius, sparse = TRUE
      )
    } else {
      hits <- sf::st_intersects(boxes[rows], unique_geom, sparse = TRUE)
      hits <- lapply(seq_along(rows), function(j) {
        candidates <- hits[[j]]
        distances <- sf::st_distance(unique_geom[rep(rows[[j]], length(candidates))],
          unique_geom[candidates],
          by_element = TRUE
        )
        candidates[as.numeric(distances) <= radius]
      })
    }
    membership[rows] <- lapply(hits, as.integer)
  }
  # Distinct domains may share a geometry, but never an omission set.
  pair <- group_ids(tibble::tibble(geometry = lookup, domain = ctx$domains), c("geometry", "domain"))
  group <- as.integer(pair)
  first <- match(seq_len(max(pair)), group)
  group_order <- order(group, seq_len(ctx$n))
  counts <- tabulate(group)
  group_ends <- cumsum(counts)
  group_starts <- c(1L, head(group_ends, -1L) + 1L)
  domain_groups <- split(seq_along(first), ctx$domains[first])
  targets <- lapply(first, function(i) {
    candidates <- domain_groups[[as.character(ctx$domains[[i]])]]
    candidates[lookup[first[candidates]] %in% membership[[lookup[[i]]]]]
  })
  target_sizes <- lengths(targets)
  targets <- as.integer(unlist(targets, use.names = FALSE))
  starts <- group_starts[targets]
  stops <- group_ends[targets]
  ends <- cumsum(as.double(target_sizes))
  accumulated <- c(0, cumsum(as.double(stops) - starts + 1))
  members <- structure(
    list(
      start = as.integer(starts), end = as.integer(stops),
      ends = ends, sizes = diff(c(0, accumulated[ends + 1])),
      order = group_order
    ),
    class = c("naboer_multiranges", "naboer_index")
  )
  out <- shared_index(members, group)
  out$group_members <- list_index(membership)
  out$observation_groups <- lookup
  out$fold_groups <- lookup
  out$domains <- ctx$domains
  out$observation_domains <- ctx$domains
  attr(out, "spatial_engine") <- spatial_engine(geom)
  out
}

#' Construct spatial cross-validation folds
#'
#' Blocks are defined independently of buffers and fold counts. Shared boundaries
#' are assigned to the first region in deterministic grid/region order. Features
#' intersecting several blocks are also assigned to the first intersected block;
#' their original geometries are preserved. Geographic grids use a local
#' azimuthal equidistant projection centred on the data bounding box unless
#' `grid_crs` is supplied. Custom regions retain their own geometry.
#' @param data A data frame, tibble, sf object, tsibble, or cubble.
#' @param size A positive scalar with physical distance units, defining grid width.
#' @param shape Square or hexagonal grid cells.
#' @param regions Optional sf/sfc polygon regions instead of a generated grid.
#' @param assignment Optional tidy selection of a column containing integer region
#'   row numbers. Required to resolve overlapping custom regions.
#' @param origin Optional numeric x/y grid anchor in the grid CRS. The grid
#'   extends in both directions from the anchor: square cell boundaries pass
#'   through it, and a hexagonal cell is centred on it. The default is
#'   the coordinate-wise floor of the bounding-box minimum divided by cell width.
#' @param grid_crs Optional projected CRS for generating grid cells.
#' @param v Optional number of folds to assign whole blocks to.
#' @param balance Balance observation counts when assigning whole blocks.
#' @param seed Optional local random seed.
#' @param key,index,.by,id Tidy selections of series keys, temporal index,
#'   independent domains, and observation identifiers.
#' @return A `naboer_folds` object.
#' @export
spatial_folds <- function(data, size = NULL, shape = c("square", "hexagonal"),
                          regions = NULL, assignment = NULL, origin = NULL,
                          grid_crs = NULL, v = NULL, balance = TRUE, seed = NULL,
                          key = NULL, index = NULL, .by = NULL, id = NULL) {
  ctx <- prepare_data(data,
    key = rlang::enquo(key), index = rlang::enquo(index),
    by = rlang::enquo(.by), id = rlang::enquo(id)
  )
  shape <- match.arg(shape)
  block <- spatial_blocks(
    ctx, size, shape, regions, rlang::enquo(assignment),
    origin, grid_crs
  )
  group_folds_from_blocks(ctx, block,
    v = v, balance = balance, seed = seed,
    spec = list(
      strategy = "spatial", size = size,
      shape = shape, origin = attr(block, "origin"),
      regions = regions, assignment = attr(block, "assignment"),
      grid_crs = attr(block, "grid_crs"),
      spatial_engine = attr(block, "spatial_engine")
    )
  )
}

spatial_blocks <- function(ctx, size = NULL, shape = "square", regions = NULL,
                           assignment = rlang::quo(NULL), origin = NULL,
                           grid_crs = NULL) {
  geom <- spatial_geometry(ctx)
  assigned <- select_cols(ctx$data, assignment, "assignment")
  if (length(assigned) > 1L) abort_naboer("`assignment` must select one column.")
  if (!is.null(regions)) {
    if (!inherits(regions, c("sf", "sfc"))) {
      abort_naboer("`regions` must be an sf or sfc polygon object.")
    }
    regions <- sf::st_set_precision(sf::st_geometry(regions), 0)
    if (is.na(sf::st_crs(regions)) || !length(regions) ||
      any(sf::st_is_empty(regions)) || !all(sf::st_is_valid(regions)) ||
      !all(sf::st_geometry_type(regions) %in% c("POLYGON", "MULTIPOLYGON"))) {
      abort_naboer("`regions` must contain nonempty valid polygons with a CRS.", "geometry")
    }
    regions <- sf::st_transform(regions, sf::st_crs(geom))
    if (!length(assigned)) {
      if (sf::st_is_longlat(regions) && sf::sf_use_s2()) {
        hits <- sf::st_intersects(regions, sparse = TRUE, model = "closed")
        overlapping <- FALSE
        for (i in seq_along(hits)) {
          others <- hits[[i]][hits[[i]] > i]
          if (!length(others)) next
          intersection <- sf::st_intersection(regions[i], regions[others])
          if (any(as.numeric(sf::st_area(intersection)) > 0)) {
            overlapping <- TRUE
            break
          }
        }
      } else {
        intersections <- suppressMessages(sf::st_relate(regions, pattern = "T********", sparse = TRUE))
        overlapping <- any(lengths(intersections) > 1L)
      }
      if (overlapping) {
        abort_naboer("Custom regions overlap in their interiors; supply an explicit `assignment` column.", "geometry")
      }
    }
  } else {
    if (length(assigned)) abort_naboer("`assignment` requires custom `regions`.")
    width <- spatial_window(size)$radius
    if (as.numeric(width) <= 0) abort_naboer("Grid `size` must be positive.")
    if (!is.null(grid_crs)) {
      if (isTRUE(sf::st_is_longlat(sf::st_crs(grid_crs)))) {
        abort_naboer("`grid_crs` must be projected.")
      }
      geom <- sf::st_transform(geom, grid_crs)
    }
    if (sf::st_is_longlat(geom)) {
      bb <- sf::st_bbox(geom)
      if (bb[[3]] - bb[[1]] > 180) {
        abort_naboer("This geographic extent spans more than 180 degrees of longitude; supply an appropriate projected `grid_crs` to resolve the grid across the antimeridian.", "geometry")
      }
      projection <- sprintf(
        "+proj=aeqd +lat_0=%.10f +lon_0=%.10f +datum=WGS84 +units=m",
        mean(bb[c(2, 4)]), mean(bb[c(1, 3)])
      )
      geom <- sf::st_transform(geom, projection)
    }
    if (sf::st_is_longlat(geom)) abort_naboer("`grid_crs` must be projected.")
    width <- as.numeric(units::set_units(width, sf::st_crs(geom)$ud_unit,
      mode = "standard"
    ))
    if (is.null(origin)) {
      origin <- floor(as.numeric(sf::st_bbox(geom)[c(1, 2)]) / width) * width
    }
    if (!is.numeric(origin) || length(origin) != 2L || any(!is.finite(origin))) {
      abort_naboer("`origin` must contain finite x and y grid coordinates.")
    }
    bb <- sf::st_bbox(geom)
    if (shape == "hexagonal") {
      regions <- spatial_hex_grid(geom, width, origin)
    } else {
      # Include the lower adjacent cell so boundary assignments do not depend
      # on which observation happens to define the data extent.
      offset <- origin + (floor((as.numeric(bb[c(1, 2)]) - origin) / width) - 1) * width
      bb[c(3, 4)] <- bb[c(3, 4)] + width * 1e-8
      cells <- prod(ceiling((as.numeric(bb[c(3, 4)]) - offset) / width))
      spatial_check_grid_size(cells)
      regions <- sf::st_make_grid(sf::st_as_sfc(bb), cellsize = width, offset = offset)
    }
  }
  if (length(assigned)) {
    block <- ctx$data[[assigned]]
    if (!is.numeric(block) || anyNA(block) || any(!is.finite(block)) ||
      any(block != floor(block)) ||
      any(block < 1L | block > length(regions))) {
      abort_naboer("`assignment` must contain valid integer region row numbers.")
    }
    matched <- sf::st_intersects(geom, regions, sparse = TRUE)
    if (any(!vapply(seq_along(block), function(i) block[[i]] %in% matched[[i]], logical(1)))) {
      abort_naboer("Each assigned region must intersect its observation geometry.")
    }
  } else {
    args <- list(geom, regions, sparse = TRUE)
    if (sf::st_is_longlat(geom) && sf::sf_use_s2()) args$model <- "closed"
    matched <- do.call(sf::st_intersects, args)
    if (any(lengths(matched) == 0L)) abort_naboer("Some observations do not intersect any spatial region.")
    block <- vapply(matched, min, integer(1))
  }
  block <- group_ids(tibble::tibble(domain = ctx$domains, block = block), c("domain", "block"))
  attr(block, "origin") <- origin
  attr(block, "grid_crs") <- sf::st_crs(geom)
  attr(block, "assignment") <- assigned
  attr(block, "spatial_engine") <- if (sf::st_is_longlat(geom) && sf::sf_use_s2()) "s2" else "GEOS"
  block
}

#' Construct spatiotemporal block folds
#'
#' @inheritParams spatial_folds
#' @param spatial_size Spatial grid width with physical distance units.
#' @param temporal_size A time span specifying temporal block width.
#' @param mode Cross space and time blocks, or hold out space/time blocks alone.
#' @param spatial_origin,temporal_origin Origins for the respective block systems.
#' @param within Temporal grouping, `"key"` or `"domain"`; the default is
#'   `"domain"` so temporal blocks align across sites.
#' @return A `naboer_folds` object.
#' @export
spatiotemporal_folds <- function(data, spatial_size = NULL, temporal_size = NULL,
                                 mode = c("crossed", "space", "time"),
                                 shape = c("square", "hexagonal"), regions = NULL,
                                 assignment = NULL, spatial_origin = NULL,
                                 temporal_origin = NULL, grid_crs = NULL,
                                 within = "domain", v = NULL, balance = TRUE,
                                 seed = NULL, key = NULL, index = NULL,
                                 .by = NULL, id = NULL) {
  mode <- match.arg(mode)
  shape <- match.arg(shape)
  ctx <- prepare_data(data,
    key = rlang::enquo(key), index = rlang::enquo(index),
    by = rlang::enquo(.by), id = rlang::enquo(id)
  )
  space <- if (mode != "time") {
    spatial_blocks(
      ctx, spatial_size, shape, regions, rlang::enquo(assignment),
      spatial_origin, grid_crs
    )
  }
  time <- if (mode != "space") temporal_blocks(ctx, temporal_size, temporal_origin, within)
  blocks <- switch(mode,
    space = space,
    time = time,
    crossed = group_ids(tibble::tibble(space = space, time = time), c("space", "time"))
  )
  group_folds_from_blocks(ctx, blocks,
    v = v, balance = balance, seed = seed,
    spec = list(
      strategy = "spatiotemporal", mode = mode,
      spatial_size = spatial_size,
      temporal_size = temporal_size, shape = shape, regions = regions,
      assignment = if (!is.null(space)) attr(space, "assignment") else character(),
      spatial_origin = if (!is.null(space)) attr(space, "origin"),
      temporal_origin = temporal_origin, within = within,
      grid_crs = if (!is.null(space)) attr(space, "grid_crs"),
      spatial_engine = if (!is.null(space)) attr(space, "spatial_engine")
    )
  )
}

# Public cubble accessors retain the original geometry and the long-form order.
cubble_data <- function(data) {
  need_package("cubble")
  key <- tsibble::key_vars(data)
  index <- cubble::index_var(data)
  nested <- inherits(data, "spatial_cubble_df")
  spatial <- tibble::as_tibble(cubble::spatial(data))
  if (!length(key) || !all(key %in% names(spatial)) || anyDuplicated(spatial[key])) {
    abort_naboer("Cubble spatial keys must uniquely identify spatial records.", "alignment")
  }
  if (nested) {
    lengths <- vapply(data$ts, nrow, integer(1))
    source <- tibble::tibble(
      .source_row = rep(seq_len(nrow(data)), lengths),
      .source_subrow = sequence(lengths)
    )
    long <- tibble::as_tibble(cubble::face_temporal(data))
  } else {
    long <- tibble::as_tibble(data)
    site <- vctrs::vec_match(long[key], spatial[key])
    if (anyNA(site)) abort_naboer("Cubble observations have unmatched spatial keys.")
    source <- tibble::tibble(
      .source_row = seq_len(nrow(long)),
      .source_subrow = NA_integer_, .spatial_row = site
    )
  }
  site <- vctrs::vec_match(long[key], spatial[key])
  for (column in setdiff(names(spatial), names(long))) {
    long[[column]] <- spatial[[column]][site]
  }
  if (!any(vapply(long, inherits, logical(1), "sfc"))) {
    need_package("sf")
    coords <- cubble::coords(data)
    if (length(coords) == 2L && all(coords %in% names(long))) {
      original_crs <- tryCatch(sf::st_crs(data), error = function(e) sf::NA_crs_)
      if (is.na(original_crs)) {
        original_crs <- tryCatch(sf::st_crs(cubble::spatial(data)),
          error = function(e) sf::NA_crs_
        )
      }
      converted <- sf::st_as_sf(long,
        coords = coords, crs = original_crs,
        remove = FALSE, na.fail = FALSE
      )
      long <- tibble::as_tibble(converted)
    }
  }
  active <- attr(long, "sf_column", exact = TRUE)
  if (is.null(active)) active <- attr(spatial, "sf_column", exact = TRUE)
  if (length(active) == 1L && active %in% names(long) && inherits(long[[active]], "sfc")) {
    attr(long, "sf_column") <- active
  }
  list(data = long, key = key, index = index, source = source)
}

# Joint spatial/time intersection without expanding a domain-wide time window
# for each focal row. Query each neighbouring station's sorted time index with
# a vector of focal bounds, retaining each result as a range.
spatiotemporal_membership <- function(ctx, spatial, window, within = "domain") {
  if (!inherits(window, "naboer_temporal_window")) {
    abort_naboer("`temporal` must be constructed with temporal_window().")
  }
  if (window$before$type == "count") {
    return(combined_index(
      list(spatial, temporal_membership(ctx, window, within)),
      "intersection"
    ))
  }
  x <- temporal_check_index(ctx)
  domains <- temporal_domain_ids(ctx, within)
  geometry <- spatial$observation_groups
  groups <- group_ids(
    tibble::tibble(domain = domains, geometry = geometry),
    c("domain", "geometry")
  )
  order <- order(groups, x, seq_len(ctx$n))
  positions <- split(seq_len(ctx$n), groups[order])
  first <- vapply(positions, function(z) order[z[[1L]]], integer(1))
  first_geometry <- geometry[first]
  first_domain <- domains[first]
  # Link geometry neighbours to temporal groups within the same domain/key.
  domain_groups <- split(seq_along(first), first_domain)
  lower <- temporal_numeric(temporal_shift(x, window$before, -1))
  upper <- temporal_numeric(temporal_shift(x, window$after, 1))
  times <- temporal_numeric(x)
  starts <- stops <- focal <- list()
  edge <- 0L
  for (group in seq_along(positions)) {
    focal_rows <- order[positions[[group]]]
    neighbours <- index_get(spatial$group_members, first_geometry[[group]])
    candidates <- domain_groups[[as.character(first_domain[[group]])]]
    target_groups <- candidates[first_geometry[candidates] %in% neighbours]
    for (target in target_groups) {
      target_positions <- positions[[target]]
      target_values <- times[order[target_positions]]
      lo <- findInterval(lower[focal_rows], target_values, left.open = TRUE) + 1L
      hi <- findInterval(upper[focal_rows], target_values)
      keep <- hi >= lo
      if (!any(keep)) next
      edge <- edge + 1L
      focal[[edge]] <- focal_rows[keep]
      starts[[edge]] <- as.integer(target_positions[[1L]] + lo[keep] - 1L)
      stops[[edge]] <- as.integer(target_positions[[1L]] + hi[keep] - 1L)
    }
  }
  focal <- unlist(focal, use.names = FALSE)
  starts <- unlist(starts, use.names = FALSE)
  stops <- unlist(stops, use.names = FALSE)
  range_order <- order(focal, starts)
  starts <- starts[range_order]
  stops <- stops[range_order]
  counts <- tabulate(focal, nbins = ctx$n)
  ends <- cumsum(as.double(counts))
  accumulated <- c(0, cumsum(as.double(stops) - starts + 1))
  sizes <- diff(c(0, accumulated[ends + 1]))
  structure(
    list(
      start = starts, end = stops, ends = ends,
      sizes = sizes, order = order
    ),
    class = c("naboer_multiranges", "naboer_index")
  )
}

# A global pointy-top lattice. Each shared vertex is computed from identical
# integer lattice coordinates, avoiding gaps from separately added offsets.
spatial_hex_grid <- function(geometry, width, origin) {
  bb <- sf::st_bbox(geometry)
  side <- width / sqrt(3)
  ystep <- 1.5 * side
  first <- ceiling((bb[[2]] - origin[[2]] - side) / ystep)
  last <- floor((bb[[4]] - origin[[2]] + side) / ystep)
  spatial_check_grid_size(last - first + 1)
  row <- seq.int(first, last)
  parity <- row %% 2
  left <- ceiling((bb[[1]] - origin[[1]]) / width - 0.5 - parity / 2)
  right <- floor((bb[[3]] - origin[[1]]) / width + 0.5 - parity / 2)
  counts <- pmax(0, right - left + 1)
  spatial_check_grid_size(sum(counts))
  # Rows are ordered south to north, and columns west to east. This ordering is
  # independent of the extent and resolves shared-boundary ties reproducibly.
  lattice_y <- rep.int(3 * row, counts)
  lattice_x <- unlist(lapply(seq_along(row), function(i) {
    if (!counts[[i]]) {
      return(numeric())
    }
    2 * seq.int(left[[i]], right[[i]]) + parity[[i]]
  }), use.names = FALSE)
  dx <- c(1, 0, -1, -1, 0, 1, 1)
  dy <- c(1, 2, 1, -1, -2, -1, 1)
  cells <- lapply(seq_along(lattice_x), function(i) {
    sf::st_polygon(list(cbind(
      origin[[1]] + (lattice_x[[i]] + dx) * width / 2,
      origin[[2]] + (lattice_y[[i]] + dy) * side / 2
    )))
  })
  sf::st_sfc(cells, crs = sf::st_crs(geometry))
}

spatial_check_grid_size <- function(n) {
  limit <- getOption("naboer.max_bytes", 1024^3)
  if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) || limit < 0) {
    abort_naboer("`naboer.max_bytes` must be nonnegative or Inf.")
  }
  if (!is.finite(n) || n > .Machine$integer.max) {
    abort_naboer("The spatial grid exceeds R's index range; increase `size` or reduce the extent.", "overflow")
  }
  # Conservative allowance for each polygon, its coordinates, and work arrays.
  if (n * 1024 > limit) {
    abort_naboer("The spatial grid exceeds `naboer.max_bytes`; increase `size`, reduce the extent, or explicitly raise the limit.", "memory")
  }
  invisible(n)
}
