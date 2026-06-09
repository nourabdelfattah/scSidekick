# =============================================================================
# scSidekick Color Palette System
# Two palettes: Nour18 (original 18-color qualitative) and
#               Nour20 (20-color full-spectrum, great for large cell-type sets)
# Both plug into the same Nour_pal / scale_color_Nour / scale_fill_Nour system.
# =============================================================================

# -----------------------------------------------------------------------------
# Nour18 - original 18-color qualitative palette
# -----------------------------------------------------------------------------

#' Nour18 color palette
#'
#' A named vector of 18 hand-curated colors spanning blues, purples, pinks,
#' reds, oranges, yellows, and greens. Suitable for coloring up to 18 cell
#' types or groups.
#'
#' @format A named character vector of length 18.
#' @export
Nour18 <- c(
  petrol      = "#003f5c",
  lightpetrol = "#008796",
  darkblue    = "#2f4b7c",
  blue        = "#227ed4",
  darkpurple  = "#665191",
  purple      = "#9e71c7",
  darkpink    = "#a05195",
  pink        = "#e082c1",
  fuchsia     = "#d45087",
  salmon      = "#f57689",
  darkred     = "#bd3e3e",
  lightred    = "#f95d6a",
  pumpkin     = "#e07b18",
  darkorange  = "#ff7c43",
  lightorange = "#ffa600",
  yellow      = "#f5d256",
  green       = "#42a665",
  grass       = "#8db032"
)

#' Access colors from the Nour18 palette by name
#'
#' @param ... Color names to retrieve. If empty, returns all 18 colors.
#' @return A named character vector of hex color codes.
#' @export
Nour_cols <- function(...) {
  cols <- c(...)
  if (is.null(cols)) return(Nour18)
  Nour18[cols]
}

# -----------------------------------------------------------------------------
# Nour20 - 20-color full-spectrum palette (blues → teals → greens →
#           yellows → oranges → reds → purples)
# -----------------------------------------------------------------------------

#' Nour20 full-spectrum color palette
#'
#' A named vector of 20 colors sweeping the full visible spectrum from deep
#' navy through teal, green, gold, orange, red, and ending at indigo. Ideal
#' for single-cell datasets with many cell types (up to 20).
#'
#' @format A named character vector of length 20.
#' @export
Nour20 <- c(
  navyblue   = "#0D47A1",
  steelblue  = "#3F81BD",
  cobalt     = "#0F6FC6",
  cerulean   = "#009DD9",
  deepteal   = "#1B587C",
  teal       = "#4BACC6",
  cyan       = "#0BD0D9",
  mint       = "#10CF9B",
  lime       = "#7CCA62",
  chartreuse = "#A5C249",
  gold       = "#FFD000",
  amber      = "#F79646",
  burnt      = "#E65100",
  coral      = "#FD817E",
  rose       = "#FD625E",
  crimson    = "#C00000",
  lavender   = "#8872C4",
  violet     = "#9B57D3",
  plum       = "#9030A0",
  indigo     = "#300890"
)

#' Access colors from the Nour20 palette by name
#'
#' @param ... Color names to retrieve. If empty, returns all 20 colors.
#' @return A named character vector of hex color codes.
#' @export
Nour_cols20 <- function(...) {
  cols <- c(...)
  if (is.null(cols)) return(Nour20)
  Nour20[cols]
}

# -----------------------------------------------------------------------------
# Named palette lists - used by Nour_pal()
# -----------------------------------------------------------------------------

#' Named list of scSidekick color palettes
#'
#' Palettes available:
#' \describe{
#'   \item{main}{8 colors from Nour18: petrol, darkblue, darkpurple, darkpink, fuchsia, lightred, darkorange, lightorange}
#'   \item{cool}{8 cool-toned colors from Nour18}
#'   \item{hot}{8 warm-toned colors from Nour18}
#'   \item{all}{All 18 colors from Nour18}
#'   \item{spectrum}{All 20 colors from Nour20}
#'   \item{blues}{5 blue-to-teal colors from Nour20}
#'   \item{warm}{6 gold-to-crimson colors from Nour20}
#'   \item{purples}{4 lavender-to-indigo colors from Nour20}
#' }
#'
#' @export
Nour_palettes <- list(
  # Nour18-based
  main    = Nour_cols("petrol", "darkblue", "darkpurple", "darkpink",
                      "fuchsia", "lightred", "darkorange", "lightorange"),
  cool    = Nour_cols("petrol", "lightpetrol", "darkblue", "blue",
                      "darkpurple", "purple", "green", "grass"),
  hot     = Nour_cols("darkred", "lightred", "pumpkin", "darkorange",
                      "lightorange", "yellow", "fuchsia", "salmon"),
  all     = Nour_cols(),

  # Nour20-based
  spectrum = Nour_cols20(),
  blues    = Nour_cols20("navyblue", "steelblue", "cobalt", "cerulean", "deepteal"),
  warm     = Nour_cols20("gold", "amber", "burnt", "coral", "rose", "crimson"),
  purples  = Nour_cols20("lavender", "violet", "plum", "indigo")
)

# -----------------------------------------------------------------------------
# Palette interpolation and ggplot2 scale constructors
# -----------------------------------------------------------------------------

#' Return an interpolated scSidekick color palette function
#'
#' @param palette Name of palette in `Nour_palettes`. Default `"main"`.
#' @param reverse Logical. Reverse the palette direction. Default `FALSE`.
#' @param ... Additional arguments passed to [grDevices::colorRampPalette()].
#' @return A palette function that takes an integer `n` and returns `n` colors.
#' @export
Nour_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- Nour_palettes[[palette]]
  if (is.null(pal)) stop("Palette '", palette, "' not found. Choose from: ",
                         paste(names(Nour_palettes), collapse = ", "))
  if (reverse) pal <- rev(pal)
  grDevices::colorRampPalette(pal, ...)
}

#' ggplot2 color scale using scSidekick palettes
#'
#' @param palette Name of palette in `Nour_palettes`. Default `"main"`.
#' @param discrete Logical. Use a discrete scale? Default `TRUE`.
#' @param reverse Logical. Reverse the palette? Default `FALSE`.
#' @param ... Additional arguments passed to [ggplot2::discrete_scale()] or
#'   [ggplot2::scale_color_gradientn()].
#' @return A ggplot2 scale object.
#' @export
scale_color_Nour <- function(palette = "main", discrete = TRUE,
                              reverse = FALSE, ...) {
  pal <- Nour_pal(palette = palette, reverse = reverse)
  if (discrete) {
    ggplot2::discrete_scale("colour", paste0("Nour_", palette), palette = pal, ...)
  } else {
    ggplot2::scale_color_gradientn(colours = pal(256), ...)
  }
}

#' ggplot2 fill scale using scSidekick palettes
#'
#' @param palette Name of palette in `Nour_palettes`. Default `"main"`.
#' @param discrete Logical. Use a discrete scale? Default `TRUE`.
#' @param reverse Logical. Reverse the palette? Default `FALSE`.
#' @param ... Additional arguments passed to [ggplot2::discrete_scale()] or
#'   [ggplot2::scale_fill_gradientn()].
#' @return A ggplot2 scale object.
#' @export
scale_fill_Nour <- function(palette = "main", discrete = TRUE,
                             reverse = FALSE, ...) {
  pal <- Nour_pal(palette = palette, reverse = reverse)
  if (discrete) {
    ggplot2::discrete_scale("fill", paste0("Nour_", palette), palette = pal, ...)
  } else {
    ggplot2::scale_fill_gradientn(colours = pal(256), ...)
  }
}


# =============================================================================
# SelectColors - generate a named color vector from any vector or factor
# =============================================================================

#' Generate a named color vector from a variable
#'
#' Accepts any vector (factor or not), resolves its levels, and returns a named
#' character vector of colors - one per level - ready to pass to
#' \code{scale_color_manual()}, \code{PrepObject(custom_colors = ...)}, or any
#' scSidekick \code{colors} argument.
#'
#' If \code{x} is not already a factor it is converted automatically and a
#' message is printed.  Supply \code{levels} to control the order; any levels
#' present in the data but not listed are appended at the end.
#'
#' @param x A vector (factor, character, or numeric).
#' @param palette Color source.  One of:
#'   \itemize{
#'     \item A \strong{Nour_pal name} string (e.g. \code{"all"}, \code{"spectrum"}).
#'     \item A \strong{character vector} of hex colors - interpolated to the
#'       required length via \code{colorRampPalette}.
#'     \item A \strong{palette function} that accepts an integer \code{n} and
#'       returns \code{n} colors (e.g. \code{scales::hue_pal()}).
#'   }
#'   Default \code{"all"}.
#' @param levels Optional character vector specifying level order.  Levels
#'   present in \code{x} but absent here are appended at the end.
#' @param reverse Logical.  Reverse the color order.  Default \code{FALSE}.
#'
#' @return A named character vector of hex colors, one per level.
#' @export
#'
#' @examples
#' \dontrun{
#' # From a Seurat metadata column
#' cols <- SelectColors(SeuratObj$seurat_clusters)
#'
#' # Custom level order + warm palette
#' cols <- SelectColors(SeuratObj$Group,
#'   levels  = c("Control", "Treated"),
#'   palette = "warm")
#'
#' # Any color vector as palette
#' cols <- SelectColors(SeuratObj$CellType,
#'   palette = c("#003f5c", "#bc5090", "#ffa600"))
#' }
SelectColors <- function(x,
                          palette = "all",
                          levels  = NULL,
                          reverse = FALSE) {

  # ── 1. Factor check / conversion ─────────────────────────────────────────
  if (!is.factor(x)) {
    x <- factor(x)   # R's default: alphabetical order, no unique() abuse
    message("scSidekick: Variable was converted to a factor (levels set alphabetically). ",
            "If you need a different order, use `levels = c(\"level1\", \"level2\", ...)`.")
  }

  # ── 2. Apply / update levels if requested ────────────────────────────────
  if (!is.null(levels)) {
    # Warn about any requested levels not actually present in the data
    phantom <- setdiff(levels, base::levels(x))
    if (length(phantom))
      warning("scSidekick: The following levels were not found in x and will be ignored: ",
              paste(phantom, collapse = ", "))
    levels <- intersect(levels, base::levels(x))        # keep only real ones
    rest   <- setdiff(base::levels(x), levels)          # unlisted → append
    x      <- factor(x, levels = c(levels, rest))
  }

  lvls <- base::levels(x)
  n    <- length(lvls)

  # ── 3. Resolve palette → color vector of length n ────────────────────────
  cols <- if (is.function(palette)) {
    out <- palette(n)
    if (reverse) rev(out) else out

  } else if (is.character(palette) && length(palette) == 1L &&
             palette %in% names(Nour_palettes)) {
    Nour_pal(palette, reverse = reverse)(n)

  } else if (is.character(palette) && length(palette) > 1L) {
    if (reverse) palette <- rev(palette)
    grDevices::colorRampPalette(palette)(n)

  } else {
    stop("'palette' must be a Nour_pal name (\"all\", \"spectrum\", etc.), ",
         "a character vector of colors, or a palette function.\n",
         "Available Nour_pal names: ",
         paste(names(Nour_palettes), collapse = ", "))
  }

  stats::setNames(cols, lvls)
}
