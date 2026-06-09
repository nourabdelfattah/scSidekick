# =============================================================================
# scSidekick - PrepObject  (prep_object.R)
#
# Exported:
#   PrepObject()   - assign colors + factor levels + analysis defaults
#                    once; all downstream scSidekick functions read from there
#   ShowColors()   - swatch plot of all stored color assignments
#   GetColors()    - retrieve stored color vector for one or all variables
#
# Internal helpers (used by all other scSidekick functions):
#   .nk_colors()   - look up colors for a metadata variable; auto-generates
#                    via Nour_pal if not stored
#   .nk_setting()  - retrieve any scalar setting (group.by, split.by, etc.)
# =============================================================================


# -----------------------------------------------------------------------------
# .nk_colors()
#
# Look up the stored color vector for a metadata variable.
# Falls back to automatic Nour_pal assignment when:
#   • PrepObject has not been run, OR
#   • the variable is not in the stored color map
#
# Palette auto-selection: ≤ 8 levels → "all" (warm/compact),  > 8 → "spectrum"
# -----------------------------------------------------------------------------
.nk_colors <- function(seurat_object, variable_name,
                        palette   = "auto",
                        reverse   = FALSE) {
  if (is.null(variable_name)) return(NULL)

  # ── Check stored color map ─────────────────────────────────────────────────
  stored <- tryCatch(
    seurat_object@misc$nk_settings$colors[[variable_name]],
    error = function(e) NULL
  )
  if (!is.null(stored) && length(stored) > 0) return(stored)

  # ── Auto-generate from Nour_pal ────────────────────────────────────────────
  if (!variable_name %in% colnames(seurat_object@meta.data)) return(NULL)

  col_data <- seurat_object@meta.data[[variable_name]]
  lvls <- if (is.factor(col_data)) levels(col_data) else
    sort(unique(as.character(col_data)))
  n <- length(lvls)

  pal_name <- if (identical(palette, "auto")) {
    if (n <= 8) "all" else "spectrum"
  } else palette

  cols <- scSidekick::Nour_pal(pal_name)(n)
  if (isTRUE(reverse)) cols <- rev(cols)
  stats::setNames(cols, lvls)
}


# -----------------------------------------------------------------------------
# .nk_setting()
# Retrieve a scalar setting stored by PrepObject (e.g. "group.by", "split.by",
# "object_name", "subset_name").  Returns NULL silently if not found.
# -----------------------------------------------------------------------------
.nk_setting <- function(seurat_object, setting_name) {
  tryCatch(
    seurat_object@misc$nk_settings[[setting_name]],
    error = function(e) NULL
  )
}


# -----------------------------------------------------------------------------
# .nk_autosave()
#
# Returns TRUE (auto-save enabled) unless PrepObject was called with
# AutoSavePlots = FALSE.  Used by every plotting function before walking up
# the stored output_dir: if FALSE the walk-up is suppressed so no files are
# written unless the caller explicitly passes output_dir in the function call.
# -----------------------------------------------------------------------------
.nk_autosave <- function(seurat_object) {
  val <- tryCatch(
    seurat_object@misc$nk_settings$autosave_plots,
    error = function(e) NULL
  )
  # Default TRUE when the setting has never been stored
  if (is.null(val)) TRUE else isTRUE(val)
}


# =============================================================================
# PrepObject
# =============================================================================

#' Prepare a Seurat Object for Uniform Downstream Analysis
#'
#' @description
#' Assigns color palettes and factor-level orderings to one or more metadata
#' variables and stores them in \code{seurat_object@misc$nk_settings}.  All
#' downstream scSidekick functions look there first for colors and default
#' parameters; they only fall back to automatic palette generation when a
#' variable is not stored.
#'
#' \strong{Additive by default.} Calling \code{PrepObject} again for a new
#' variable (e.g. after cell-type assignment) does not overwrite existing
#' settings - it only adds or updates the variables you name.  Set
#' \code{force = TRUE} to regenerate a variable that is already stored.
#'
#' \strong{Auto-palette selection:}
#' \itemize{
#'   \item \code{"auto"}: \code{Nour_pal("all")} for variables with \eqn{\leq 8}
#'     levels (warm, distinct), \code{Nour_pal("spectrum")} for \eqn{> 8} levels
#'     (broader color range).
#'   \item Specify per-variable, either \strong{by name}
#'     \code{palettes = c(Cluster = "all", Sample = "spectrum")} or
#'     \strong{by position} (one entry per \code{variables}, same order)
#'     \code{palettes = c("spectrum", "all", "all")}.
#' }
#'
#' @param seurat_object A Seurat object.
#' @param variables Character vector of metadata column names to color and
#'   prepare.  Can be a single name; re-run \code{PrepObject} at any time to
#'   add more.
#' @param palettes Palette name(s).  One of:
#'   \itemize{
#'     \item A \strong{single string} (\code{"auto"}, \code{"all"},
#'       \code{"spectrum"}, \dots) applied to every variable.
#'     \item A \strong{named} character vector keyed by variable
#'       (\code{c(Cluster = "all", Sample = "spectrum")}); names not listed
#'       fall back to \code{"auto"}.
#'     \item An \strong{unnamed} character vector of the same length as
#'       \code{variables}, matched \emph{by position}
#'       (\code{c("spectrum", "all", "all")}).
#'   }
#'   Default \code{"auto"}.
#' @param reverse Logical, named logical vector, or an unnamed logical vector
#'   matched by position to \code{variables}.  \code{TRUE} reverses the palette
#'   for all variables; a named vector reverses selectively
#'   (e.g. \code{c(Cluster = TRUE, Sample = FALSE)}).  Default \code{FALSE}.
#' @param custom_colors Named list of named character vectors that override
#'   automatic palette generation entirely for specific variables.
#'   E.g. \code{list(Sample = c(Control = "#3B4CC0", NRF1_mRNA = "#B40426"))}.
#' @param var_levels Named list of character vectors that set the \emph{factor
#'   level ordering} for each variable, controlling sort order in heatmap
#'   columns, violin plots, and ggplot2 facets.
#'   E.g. \code{list(Sample = c("Control", "NRF1_mRNA"))}.
#' @param group.by Character.  Store the preferred primary grouping column.
#'   Downstream scSidekick functions use this as the default \code{group.by} when
#'   none is specified in the call.  Pass \code{NA} to explicitly clear a
#'   previously stored value.
#' @param split.by Character.  Store the preferred secondary split column.
#'   Downstream functions use this as the default \code{split.by}.  Pass
#'   \code{NA} to explicitly disable the stored default (even if one exists).
#' @param output_dir Character.  Store the preferred output directory for
#'   saving PDFs and plots.  All downstream scSidekick functions walk up to
#'   this value when their own \code{output_dir} argument is \code{NULL}.
#'   Pass \code{NA} to explicitly clear a stored value.
#' @param object_name Character.  Label used to prefix output file names.
#' @param subset_name Character.  Optional subset label appended to output
#'   file names (e.g. \code{"NoNeurons"}).
#' @param AutoSavePlots Logical or \code{NULL}.  Controls whether plotting
#'   functions automatically save PDFs to the stored \code{output_dir} when no
#'   \code{output_dir} is passed explicitly in the function call.
#'   \itemize{
#'     \item \code{TRUE} (default when not set): plotting functions walk up and
#'       use the stored \code{output_dir} automatically — existing behaviour.
#'     \item \code{FALSE}: the walk-up is suppressed; functions return the plot
#'       object without saving unless the caller explicitly passes
#'       \code{output_dir = "..."} in the plotting call.
#'   }
#'   Pass \code{NULL} to leave the current setting unchanged.
#' @param force Logical.  If \code{FALSE} (default), variables that already
#'   have stored colors are left unchanged.  If \code{TRUE}, existing colors
#'   are regenerated for every variable listed in \code{variables}.
#'
#' @return The \code{seurat_object} with \code{@misc$nk_settings} updated.
#'   Return the object and re-assign: \code{seurat_object <- PrepObject(...)}.
#'
#' @seealso \code{\link{ShowColors}}, \code{\link{GetColors}}
#'
#' @export
PrepObject <- function(
    seurat_object,
    variables,
    palettes      = "auto",
    reverse       = FALSE,
    custom_colors = list(),
    var_levels    = list(),
    group.by      = NULL,
    split.by      = NULL,
    output_dir    = NULL,
    object_name   = NULL,
    subset_name   = NULL,
    AutoSavePlots = NULL,
    force         = FALSE
) {

  # ── 0. Validate ─────────────────────────────────────────────────────────────
  if (!inherits(seurat_object, "Seurat"))
    stop("'seurat_object' must be a Seurat object.")
  if (missing(variables) || length(variables) == 0)
    stop("'variables' must be a non-empty character vector of metadata column names.")

  # Check all variables exist in metadata
  bad_vars <- setdiff(variables, colnames(seurat_object@meta.data))
  if (length(bad_vars) > 0)
    stop("Variable(s) not found in seurat_object@meta.data:\n  ",
         paste(bad_vars, collapse = ", "))

  # ── 1. Initialise nk_settings if absent ────────────────────────────────────
  if (is.null(seurat_object@misc$nk_settings))
    seurat_object@misc$nk_settings <- list(colors = list(), levels = list())
  if (is.null(seurat_object@misc$nk_settings$colors))
    seurat_object@misc$nk_settings$colors <- list()
  if (is.null(seurat_object@misc$nk_settings$levels))
    seurat_object@misc$nk_settings$levels <- list()

  # ── 2. Resolve palette + reverse per variable ───────────────────────────────
  # Per-variable lookup priority:
  #   1. NAMED vector/list  -> look up by variable name
  #   2. POSITIONAL vector  -> length == length(variables): map by position
  #      (so palettes = c("spectrum","all","all") follows `variables` order)
  #   3. SINGLE value       -> applied to every variable
  #   4. otherwise          -> "auto"
  .pal_for <- function(var) {
    if (is.character(palettes) && !is.null(names(palettes)))
      palettes[[var]] %||% "auto"
    else if (length(palettes) == length(variables))
      palettes[[match(var, variables)]]
    else if (length(palettes) == 1)
      palettes
    else "auto"
  }
  .rev_for <- function(var) {
    if (is.logical(reverse) && length(reverse) == 1) return(isTRUE(reverse))
    if (!is.null(names(reverse)) && var %in% names(reverse))
      return(isTRUE(reverse[[var]]))
    if (length(reverse) == length(variables))
      return(isTRUE(reverse[[match(var, variables)]]))
    FALSE
  }
  # Simple null-coalescing helper (internal only)
  # ── 3. Assign colors per variable ──────────────────────────────────────────
  for (var in variables) {
    col_data <- seurat_object@meta.data[[var]]

    # Factor levels: use var_levels if supplied, else existing factor or sorted unique
    # (computed BEFORE the skip check so we can detect new levels)
    if (!is.null(var_levels[[var]])) {
      lvls <- var_levels[[var]]
    } else if (is.factor(col_data)) {
      lvls <- levels(col_data)
    } else {
      lvls <- sort(unique(as.character(col_data)))
    }
    n <- length(lvls)

    stored_cols    <- seurat_object@misc$nk_settings$colors[[var]]
    already_stored <- !is.null(stored_cols) && length(stored_cols) > 0
    new_lvls       <- setdiff(lvls, names(stored_cols))

    if (already_stored && length(new_lvls) == 0 && !isTRUE(force)) {
      message("  ", var, ": already stored (use force = TRUE to regenerate). Skipping.")
      next
    }
    if (already_stored && length(new_lvls) > 0) {
      message("  ", var, ": new level(s) detected - regenerating colors: ",
              paste(new_lvls, collapse = ", "))
    }

    # Custom colors take priority
    if (!is.null(custom_colors[[var]])) {
      cols <- custom_colors[[var]]
      # Ensure all levels are covered; fill missing with grey
      missing_lvls <- setdiff(lvls, names(cols))
      if (length(missing_lvls) > 0) {
        greys <- grDevices::grey(seq(0.4, 0.8, length.out = length(missing_lvls)))
        cols  <- c(cols, stats::setNames(greys, missing_lvls))
      }
      cols <- cols[lvls]  # ensure correct order
    } else {
      # Automatic palette
      pal_name <- .pal_for(var)
      if (identical(pal_name, "auto")) pal_name <- if (n <= 8) "all" else "spectrum"
      rev_pal  <- .rev_for(var)
      cols     <- scSidekick::Nour_pal(pal_name)(n)
      if (rev_pal) cols <- rev(cols)
      cols     <- stats::setNames(cols, lvls)
    }

    seurat_object@misc$nk_settings$colors[[var]] <- cols

    # Set factor levels on the metadata column
    seurat_object@meta.data[[var]] <- factor(
      as.character(seurat_object@meta.data[[var]]),
      levels = lvls
    )
    # Propagate to nk_settings$levels for reference
    seurat_object@misc$nk_settings$levels[[var]] <- lvls

    message("  ", var, ": ", n, " levels - ",
            if (!is.null(custom_colors[[var]])) "custom colors"
            else paste0("Nour_pal('",
                        if (identical(.pal_for(var), "auto"))
                          if (n <= 8) "all" else "spectrum"
                        else .pal_for(var),
                        "')", if (.rev_for(var)) " (reversed)" else ""))
  }

  # ── 4. Store scalar analysis defaults (only update when explicitly supplied) ─
  # NA = explicitly clear; NULL = do not touch
  if (!is.null(group.by))
    seurat_object@misc$nk_settings$group.by <-
      if (isTRUE(is.na(group.by))) NULL else as.character(group.by)

  if (!is.null(split.by))
    seurat_object@misc$nk_settings$split.by <-
      if (isTRUE(is.na(split.by))) NULL else as.character(split.by)

  if (!is.null(output_dir))
    seurat_object@misc$nk_settings$output_dir <-
      if (isTRUE(is.na(output_dir))) NULL else as.character(output_dir)

  if (!is.null(object_name))
    seurat_object@misc$nk_settings$object_name <- as.character(object_name)

  if (!is.null(subset_name))
    seurat_object@misc$nk_settings$subset_name <- as.character(subset_name)

  if (!is.null(AutoSavePlots))
    seurat_object@misc$nk_settings$autosave_plots <- isTRUE(AutoSavePlots)

  # ── 5. Summary ───────────────────────────────────────────────────────────────
  cfg <- seurat_object@misc$nk_settings
  message("\nPrepObject complete.")
  message("  Variables with stored colors: ",
          paste(names(cfg$colors), collapse = ", "))
  if (!is.null(cfg$group.by))
    message("  Default group.by:  ", cfg$group.by)
  if (!is.null(cfg$split.by))
    message("  Default split.by:  ", cfg$split.by)
  if (!is.null(cfg$output_dir))
    message("  Default output_dir: ", cfg$output_dir)
  if (!is.null(cfg$object_name))
    message("  object_name: ", cfg$object_name)
  if (!is.null(cfg$subset_name) && nchar(cfg$subset_name) > 0)
    message("  subset_name: ", cfg$subset_name)
  if (!is.null(cfg$autosave_plots))
    message("  AutoSavePlots: ", cfg$autosave_plots)

  invisible(seurat_object)
}


# =============================================================================
# ShowColors
# =============================================================================

#' Display All Stored Color Assignments
#'
#' @description
#' Plots a swatch grid of every variable stored by \code{\link{PrepObject}}.
#' One row per variable, one tile per level, labeled with the level name.
#' Useful for verifying color assignments before running analyses.
#'
#' @param seurat_object A Seurat object (must have run \code{PrepObject} first).
#'
#' @return A \code{ggplot2} object (invisibly). The plot is printed to the
#'   active device.
#'
#' @export
ShowColors <- function(seurat_object) {
  cfg <- tryCatch(seurat_object@misc$nk_settings$colors, error = function(e) NULL)
  if (is.null(cfg) || length(cfg) == 0) {
    message("No colors stored. Run PrepObject() first.")
    return(invisible(NULL))
  }

  # Build a long data frame: variable | level | hex_color | level_index
  rows <- lapply(names(cfg), function(var) {
    cols <- cfg[[var]]
    data.frame(
      variable    = var,
      level       = names(cols),
      color       = unname(cols),
      level_index = seq_along(cols),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df$variable <- factor(df$variable, levels = rev(names(cfg)))  # top-to-bottom order

  p <- ggplot2::ggplot(df, ggplot2::aes(x = level_index, y = variable,
                                         fill = color)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4, width = 0.9, height = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = level),
                       size = 2.8, fontface = "bold", color = "black") +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::labs(title = "scSidekick stored color assignments",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_blank(),
      axis.ticks.x     = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_text(face = "bold", colour = "black",
                                                size = 9),
      panel.grid       = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                size = 10),
      plot.margin      = ggplot2::unit(c(0.3, 0.5, 0.3, 0.5), "cm")
    )

  print(p)
  invisible(p)
}


# =============================================================================
# GetColors
# =============================================================================

#' Retrieve Stored Color Vectors from a Prepared Seurat Object
#'
#' @description
#' Returns the named color vector for one metadata variable, or the complete
#' color list if no variable is specified.
#'
#' @param seurat_object A Seurat object.
#' @param variable Character.  Name of the metadata variable to retrieve.
#'   \code{NULL} (default) returns the full list for all stored variables.
#'
#' @return A named character vector (one variable) or a named list (all
#'   variables), or \code{NULL} if nothing is stored.
#'
#' @export
GetColors <- function(seurat_object, variable = NULL) {
  cfg <- tryCatch(seurat_object@misc$nk_settings$colors, error = function(e) NULL)
  if (is.null(cfg)) {
    message("No colors stored. Run PrepObject() first.")
    return(invisible(NULL))
  }
  if (is.null(variable)) cfg else cfg[[variable]]
}
