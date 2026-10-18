# Rebuild the vector logo and local pkgdown favicon assets from the repo root:
# Rscript tools/build-branding.R
# Development-only dependencies: magick, sf, xml2, jsonlite and systemfonts.
# The approved raster reference is retained under tools/branding. Houses and
# foliage are contour-traced; the foreground sign and lettering are redrawn.
# All displayed artwork is vector geometry, with no embedded raster image.

required <- c("magick", "sf", "xml2", "jsonlite", "systemfonts")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install development dependencies: ", paste(missing, collapse = ", "))
stopifnot(file.exists("DESCRIPTION"), file.exists("tools/branding/ramsay-street-reference.png"))
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("pkgdown/favicon", recursive = TRUE, showWarnings = FALSE)

teal <- "#003F48"
cream <- "#FFF0CF"
coral <- "#F67C56"
sage <- "#84AAA0"

svg_document <- function(body, view_box, title, description) {
  paste0('<?xml version="1.0" encoding="UTF-8"?>\n',
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="', view_box,
    '" role="img" aria-labelledby="title desc">\n',
    '<title id="title">', title, '</title>\n',
    '<desc id="desc">', description, '</desc>\n', body, '\n</svg>\n')
}

write_svg <- function(body, file, view_box = "0 0 256 256", title, description) {
  writeLines(svg_document(body, view_box, title, description), file, sep = "")
  invisible(xml2::read_xml(file))
}

read_svg <- function(svg, size) {
  # Use ImageMagick's SVG delegate; magick::image_read_svg() additionally needs
  # the separate rsvg R package, which this workflow does not require.
  # The delegate paints unused canvas white. None of the artwork uses pure
  # white, so remove that canvas before downsampling, retaining alpha edges.
  magick::image_resize(magick::image_transparent(
    magick::image_read(svg, density = "192"), "white"),
    sprintf("%dx%d!", size, size))
}

render_png <- function(svg, file, size) {
  img <- read_svg(svg, size)
  magick::image_write(img, file, format = "png")
  invisible(img)
}

# Trace a softened, reduced palette, rather than the reference's grain and
# anti-aliasing noise. Each colour consists of closed contours with holes.
reference <- magick::image_read("tools/branding/ramsay-street-reference.png")
reference <- magick::image_resize(reference, "1000x1000!")
reference <- magick::image_blur(reference, radius = 0, sigma = 0.65)
reference <- magick::image_quantize(reference, max = 40, colorspace = "sRGB", dither = FALSE)
pixels <- magick::image_data(reference, channels = "rgb")
colours <- matrix(grDevices::rgb(as.integer(pixels[1, , ]),
  as.integer(pixels[2, , ]), as.integer(pixels[3, , ]), maxColorValue = 255), 1000, 1000)
palette <- names(sort(table(colours), decreasing = TRUE))
paths <- vapply(palette, function(colour) {
  mask <- matrix(0, 1002, 1002)
  mask[2:1001, 2:1001] <- colours == colour
  contours <- grDevices::contourLines(x = seq(-0.5, 1000.5),
    y = seq(-0.5, 1000.5), z = mask, levels = 0.5)
  parts <- vapply(contours, function(contour) {
    xy <- cbind(contour$x, contour$y)
    area <- abs(sum(xy[-nrow(xy), 1] * xy[-1, 2] -
      xy[-1, 1] * xy[-nrow(xy), 2])) / 2
    if (area < 3) return("")
    line <- sf::st_simplify(sf::st_sfc(sf::st_linestring(xy)), dTolerance = 0.9)
    xy <- sf::st_coordinates(line)[, 1:2, drop = FALSE] * 1.254
    xy <- xy[-nrow(xy), , drop = FALSE]
    midpoint <- (xy + xy[c(2:nrow(xy), 1), , drop = FALSE]) / 2
    paste0("M", sprintf("%.2f,%.2f", midpoint[nrow(midpoint), 1], midpoint[nrow(midpoint), 2]),
      paste(sprintf("Q%.2f,%.2f %.2f,%.2f", xy[, 1], xy[, 2],
        midpoint[, 1], midpoint[, 2]), collapse = ""), "Z")
  }, character(1))
  paste0('<path fill="', colour, '" stroke="', colour,
    '" stroke-width="1.5" stroke-linejoin="round" fill-rule="evenodd" d="',
    paste0(parts, collapse = ""), '"/>')
}, character(1))

# Clip the original off-white exterior away while retaining the teal hex edge.
hex <- "M626 12L1171 326L1171 927L626 1242L81 927L81 326Z"
logo <- paste0('<defs><clipPath id="hex"><path d="', hex,
  '"/></clipPath><clipPath id="scene"><rect x="110" y="200" width="1031" height="639"/>',
  '</clipPath></defs>\n<g clip-path="url(#hex)">\n',
  '<path fill="#FFF3D7" d="', hex, '"/>\n<g clip-path="url(#scene)">',
  paste(paths, collapse = "\n"), '\n</g></g>')
# Glyph outlines make all lettering independent of fonts at display time.
# The Avenir families are needed only when rebuilding the artwork.
outlined_letters <- function(text, x, y, width, height, family = "Avenir") {
  font <- systemfonts::match_fonts(family, weight = "heavy")
  glyphs <- systemfonts::glyph_info(text, path = font$path, index = font$index, size = 100)
  outlines <- systemfonts::glyph_outline(glyphs$index, path = font$path,
    index = font$index, size = 100, tolerance = 0.1)
  offset <- c(0, cumsum(glyphs$x_advance[-nrow(glyphs)]))
  outlines$x <- outlines$x + offset[outlines$glyph]
  outlines$y <- -outlines$y
  outlines$x <- x + (outlines$x - min(outlines$x)) / diff(range(outlines$x)) * width
  outlines$y <- y + (outlines$y - min(outlines$y)) / diff(range(outlines$y)) * height
  contours <- split(outlines, interaction(outlines$glyph, outlines$contour, drop = TRUE))
  d <- vapply(contours, function(z) paste0("M",
    paste(sprintf("%.2f,%.2f", z$x, z$y), collapse = "L"), "Z"), character(1))
  paste0('<path fill="', cream, '" fill-rule="evenodd" d="', paste0(d, collapse = ""), '"/>')
}

# Redraw the foreground elements so the badge edge, sign and type remain crisp
# at print sizes; retain the approved scene's composition and perspective.
foreground <- paste0('<g clip-path="url(#hex)">',
  '<path d="', hex, '" fill="none" stroke="', teal, '" stroke-width="60"/>',
  '<path fill="#FFF3D7" d="M414 220L1097 370V610L414 460Z"/>',
  '<g transform="matrix(1 .22 0 1 420 247)">',
  '<rect width="665" height="202" rx="17" fill="#0B5046" stroke="', teal, '" stroke-width="14"/>',
  '<rect x="16" y="16" width="633" height="170" rx="11" fill="none" stroke="', cream, '" stroke-width="9"/>',
  outlined_letters("Ramsay St.", 43, 38, 578, 130, "Avenir Next Condensed"),
  '</g>',
  '<path fill="', teal, '" d="M81 839H1171V927L626 1242L81 927Z"/>',
  outlined_letters("naboer", 253, 859, 777, 164), '</g>')
write_svg(paste0(logo, foreground), "man/figures/logo.svg", "0 0 1254 1254",
  "naboer: Ramsay Street neighbourhood",
  "A green Ramsay St. sign in perspective above suburban houses, with naboer below, inside a teal hexagon. The exterior is transparent.")
render_png("man/figures/logo.svg", "man/figures/logo.png", 1254)
render_png("man/figures/logo.svg", "tools/branding/logo-preview.png", 320)

icon_background <- paste0('<rect width="256" height="256" rx="48" fill="', teal, '"/>')
sign_start <- paste0(icon_background,
  '<rect x="47" y="65" width="18" height="166" rx="6" fill="', sage, '"/>',
  '<g transform="matrix(1 .15 0 1 0 -18)">',
  '<rect x="32" y="46" width="204" height="132" rx="15" fill="', teal,
  '" stroke="', cream, '" stroke-width="9"/>')
letter_n <- paste0('<path fill="', cream, '" d="M88 78H111V89',
  'C118 79 129 75 141 75C163 75 177 90 177 114V156H152V116',
  'C152 104 146 98 136 98C121 98 113 108 113 122V156H88Z"/>')
n_sign <- paste0(sign_start, letter_n, '</g>')
rs_sign <- paste0(sign_start, outlined_letters("RS", 73, 77, 116, 80), '</g>')
house <- paste0(icon_background,
  '<path fill="', coral, '" d="M40 121L128 39L216 121H193V216H63V121Z"/>',
  '<path fill="', cream, '" d="M110 148H146V216H110Z"/>')
write_svg(n_sign, "pkgdown/favicon/favicon.svg", title = "naboer street-sign icon",
  description = "A cream lowercase n on an angled teal street sign.")
write_svg(rs_sign, "tools/branding/favicon-rs.svg", title = "Ramsay Street monogram icon",
  description = "Cream RS lettering on an angled teal street sign.")
write_svg(house, "tools/branding/favicon-house.svg", title = "naboer house icon",
  description = "A coral house with a cream doorway on teal.")

for (size in c(16, 32, 48, 96)) {
  render_png("pkgdown/favicon/favicon.svg",
    sprintf("pkgdown/favicon/favicon-%dx%d.png", size, size), size)
}
render_png("pkgdown/favicon/favicon.svg", "pkgdown/favicon/apple-touch-icon.png", 180)
for (size in c(192, 512)) {
  render_png("pkgdown/favicon/favicon.svg",
    sprintf("pkgdown/favicon/web-app-manifest-%dx%d.png", size, size), size)
}
ico <- magick::image_join(lapply(c(16, 32, 48), function(size)
  magick::image_read(sprintf("pkgdown/favicon/favicon-%dx%d.png", size, size))))
magick::image_write(ico, "pkgdown/favicon/favicon.ico", format = "ico")
jsonlite::write_json(list(name = "naboer", short_name = "naboer",
  icons = lapply(c(192, 512), function(size) list(
    src = sprintf("web-app-manifest-%dx%d.png", size, size),
    sizes = sprintf("%dx%d", size, size), type = "image/png", purpose = "any")),
  start_url = "./", scope = "./", display = "standalone",
  theme_color = teal, background_color = cream),
  "pkgdown/favicon/site.webmanifest", auto_unbox = TRUE, pretty = TRUE)

# A comparison at actual 16/32/48 px sizes, plus enlarged views, on light and
# dark backgrounds. Upscaling uses point filtering to expose tiny-size detail.
icon_files <- c("pkgdown/favicon/favicon.svg", "tools/branding/favicon-rs.svg",
  "tools/branding/favicon-house.svg")
labels <- c("n street sign (default)", "RS street sign", "Coral house")
rows <- lapply(seq_along(icon_files), function(i) {
  row <- magick::image_blank(900, 240, "#FFF9ED")
  row <- magick::image_annotate(row, labels[i], location = "+24+26",
    size = 22, color = teal, font = "Helvetica")
  large <- read_svg(icon_files[i], 128)
  row <- magick::image_composite(row, large, offset = "+30+76")
  for (j in seq_along(c(16, 32, 48))) {
    size <- c(16, 32, 48)[j]
    small <- read_svg(icon_files[i], size)
    left <- 230 + (j - 1) * 190
    row <- magick::image_annotate(row, paste0(size, " px"),
      location = sprintf("+%d+80", left), size = 17, color = teal, font = "Helvetica")
    row <- magick::image_composite(row, small, offset = sprintf("+%d+110", left))
    dark <- magick::image_blank(92, 70, "#15242B")
    dark <- magick::image_composite(dark, small, offset = sprintf("+12+%d", (70 - size) %/% 2))
    row <- magick::image_composite(row, dark, offset = sprintf("+%d+164", left))
  }
  row
})
magick::image_write(magick::image_append(magick::image_join(rows), stack = TRUE),
  "tools/branding/favicon-comparison.png")
message("Built vector logo, PNG and local pkgdown favicons.")
