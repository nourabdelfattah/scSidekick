# =============================================================================
# scSidekick — Group Heatmap  (group_heatmap.R)
#
# Exported:
#   GroupHeatmap() — ComplexHeatmap of mean z-scored expression per group,
#                    with optional dot overlay (% expressing), row split by
#                    feature category, column split / faceting by secondary
#                    metadata, and textbox row labels.
#
# Key design choices:
#   • group.by is a CHARACTER VECTOR: each unique combination of its levels
#     becomes one heatmap column.  Pass c("CellType", "Treatment") to get
#     every CellType × Treatment pair automatically.
#   • column_split_by FACETS: it creates separate column panels.  Every
#     group.by combination appears inside each split panel — it does NOT
#     try to assign each group level to a single split category.
#   • add_bg draws a white disc behind the dot, not a full-cell rectangle.
# =============================================================================


#' Group-Level Expression Heatmap with Optional Dot Overlay
#'
#' @description
#' Computes mean expression per combination of \code{group.by} levels for the
#' supplied \code{features}, scales across groups (z-score by default), and
#' draws a \code{ComplexHeatmap} with:
#' \itemize{
#'   \item \strong{Colour fill} = scaled mean expression
#'   \item \strong{Dot overlay} (optional) — dots sized by \% of cells
#'     expressing each gene within each column group
#'   \item \strong{Row split + textbox labels} — rows grouped and annotated
#'     by \code{feature_split} categories
#'   \item \strong{Column faceting} — each level of \code{column_split_by}
#'     forms a panel; within every panel \emph{all} \code{group.by}
#'     combinations are shown
#'   \item \strong{Column annotation} — one colour bar per \code{group.by}
#'     variable (analogous to ComplexHeatmap's multi-bar top annotation)
#' }
#'
#' @section Multiple group.by variables:
#' \code{group.by} may be a character vector of length > 1.  The heatmap
#' columns are the \emph{observed} unique combinations of all listed
#' variables, ordered by the first variable, then the second, etc.
#' The top annotation automatically shows one colour bar per variable.
#'
#' @section column_split_by behaviour:
#' When \code{column_split_by} is supplied the function computes expression
#' for \emph{every} (group.by × column_split_by) combination found in the
#' object.  This means a CellType that appears in both "Vehicle" and "Drug"
#' produces \emph{two} columns (one per treatment), placed in their
#' respective panels.  This is the correct faceting behaviour.
#'
#' @param seurat_object A Seurat object.
#' @param group.by Character vector of metadata column(s) whose unique
#'   combinations become the heatmap columns (e.g.
#'   \code{c("CellTypes", "Treatment")}).
#' @param features Character vector of gene names to plot.
#' @param feature_split Named character or factor vector (names = gene
#'   symbols, values = row-group labels).  When supplied, rows are split
#'   and labelled.
#' @param column_split_by Optional metadata column used to \emph{facet}
#'   (visually split) the columns into panels.  May be one of the
#'   \code{group.by} variables or a separate one.
#' @param assay Assay to use. Default \code{"RNA"}.
#' @param layer Layer for expression extraction. Default \code{"data"}.
#' @param scale_method \code{"zscore"} (default), \code{"minmax"}, or
#'   \code{"none"}.
#' @param add_dot Logical.  Overlay dots sized by \% expressing. Default
#'   \code{TRUE}.
#' @param add_bg Logical.  Draw a white disc behind each dot so it pops
#'   against dark expression colours.  Default \code{TRUE}.
#' @param dot_size Base dot size as a \code{grid::unit} object. Default
#'   \code{grid::unit(4, "mm")}.
#' @param pct_cutoff Minimum fraction expressing (0–1) to draw a dot.
#'   Default \code{0}.
#' @param palette RColorBrewer diverging palette for the expression colour
#'   scale.  Default \code{"RdBu"}.
#' @param color_scale Optional pre-built \code{circlize::colorRamp2} object.
#'   Overrides \code{palette}.
#' @param scale_limits \code{c(min, max)} for the colour scale. Default
#'   \code{c(-2, 2)}.
#' @param group_colors Colours for \code{group.by} levels.  Accepts:
#'   \itemize{
#'     \item A \strong{named list} where each element is a named colour
#'       vector for one \code{group.by} variable, e.g.
#'       \code{list(CellType = ct_cols, Treatment = tx_cols)}.
#'     \item A \strong{named vector} when \code{group.by} has length 1
#'       (backward-compatible).
#'     \item \code{NULL} (default) — auto-assigned from PrepObject or
#'       \code{Nour_pal}.
#'   }
#' @param column_split_colors Named colour vector for \code{column_split_by}
#'   levels.
#' @param feature_split_colors Named colour vector for \code{feature_split}
#'   categories (cosmetic; not currently used in cell_fun).
#' @param textbox_labels Logical.  Use \code{anno_textbox} when available
#'   (\pkg{ComplexHeatmap} ≥ 2.13).  Default \code{TRUE}.
#' @param cluster_rows Logical.  Cluster rows within each split.  Default
#'   \code{TRUE}.
#' @param cluster_row_slices Logical.  Cluster row slices.  Default
#'   \code{FALSE}.
#' @param cluster_columns Logical.  Cluster columns within each split.
#'   Default \code{FALSE}.
#' @param heatmap_params Named list of extra arguments forwarded to
#'   \code{ComplexHeatmap::Heatmap()}.
#' @param auto_draw Logical.  Render automatically after building. Default
#'   \code{TRUE}.
#' @param legend_side Where to place the legend.  Default \code{"right"}.
#' @param width Heatmap body width as a \code{grid::unit}. Default
#'   \code{grid::unit(5, "in")}.
#' @param height Heatmap body height.  Default \code{grid::unit(7, "in")}.
#' @param output_dir Directory to save a PDF.  \code{NULL} = no file.
#' @param object_name Prefix for output file names.
#' @param subset_name Optional subset label.
#'
#' @return Invisibly, a named list:
#' \describe{
#'   \item{\code{heatmap}}{The \code{Heatmap} object.}
#'   \item{\code{expression_matrix}}{Scaled mean-expression matrix.}
#'   \item{\code{percent_matrix}}{Fraction expressing per column group.}
#'   \item{\code{column_groups}}{Data frame of all column combinations.}
#'   \item{\code{feature_split}}{The \code{feature_split} vector used.}
#'   \item{\code{dot_legend}}{The dot-size \code{Legend} object or
#'     \code{NULL}.}
#' }
#'
#' @export
GroupHeatmap <- function(
    seurat_object,
    group.by,
    features,
    feature_split         = NULL,
    column_split_by       = NULL,
    assay                 = "RNA",
    layer                 = "data",
    scale_method          = "zscore",
    add_dot               = TRUE,
    add_bg                = TRUE,
    dot_size              = grid::unit(4, "mm"),
    pct_cutoff            = 0,
    palette               = "RdBu",
    color_scale           = NULL,
    scale_limits          = c(-2, 2),
    group_colors          = NULL,
    column_split_colors   = NULL,
    feature_split_colors  = NULL,
    textbox_labels        = TRUE,
    cluster_rows          = TRUE,
    cluster_row_slices    = FALSE,
    cluster_columns       = FALSE,
    heatmap_params        = list(
      show_row_names      = TRUE,
      row_names_side      = "right",
      row_names_gp        = grid::gpar(fontsize = 8),
      show_row_dend       = FALSE,
      show_column_names   = TRUE,
      column_names_side   = "top",
      row_title           = NULL,
      column_title        = NULL
    ),
    auto_draw             = TRUE,
    legend_side           = "right",
    width                 = grid::unit(5, "in"),
    height                = grid::unit(7, "in"),
    output_dir            = NULL,
    object_name           = "",
    subset_name           = ""
) {

  # ── Walk-up PrepObject defaults ────────────────────────────────────────────
  output_dir  <- output_dir %||% .nk_setting(seurat_object, "output_dir")
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  # ── 0. Validate ─────────────────────────────────────────────────────────────
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("ComplexHeatmap is required. Install: BiocManager::install('ComplexHeatmap')")
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("circlize is required. Install: BiocManager::install('circlize')")

  scale_method <- match.arg(scale_method, c("zscore", "minmax", "none"))
  meta <- seurat_object@meta.data

  # Validate all metadata columns
  all_check_cols <- unique(c(group.by, if (!is.null(column_split_by)) column_split_by))
  for (col in all_check_cols) {
    if (!col %in% colnames(meta))
      stop("'", col, "' not found in seurat_object@meta.data.")
  }

  # ── 1. Feature validation ────────────────────────────────────────────────────
  avail_genes <- tryCatch(
    rownames(SeuratObject::GetAssay(seurat_object, assay = assay)),
    error = function(e) {
      tryCatch(rownames(seurat_object[[assay]]), error = function(e2) rownames(seurat_object))
    }
  )
  features_use <- intersect(features, avail_genes)
  if (length(features_use) == 0L)
    stop("None of the supplied features found in assay '", assay, "'.")
  if (length(features_use) < length(features))
    message("  ", length(features) - length(features_use),
            " feature(s) not found in assay '", assay, "' — skipped.")

  # Align feature_split to available features (preserves split order)
  if (!is.null(feature_split)) {
    feature_split <- feature_split[names(feature_split) %in% features_use]
    features_use  <- names(feature_split)
  }

  # ── 2. Build column combination table ────────────────────────────────────────
  # All metadata variables that define columns (group.by + optional split)
  all_group_vars <- unique(c(group.by,
                              if (!is.null(column_split_by)) column_split_by))

  # Helper: ordered levels for a metadata variable
  .lvls <- function(v) {
    x <- meta[[v]]
    if (is.factor(x)) levels(x) else sort(unique(as.character(x)))
  }

  # Convert all relevant meta columns to character for consistent matching
  meta_chr <- as.data.frame(
    lapply(meta[, all_group_vars, drop = FALSE], as.character),
    stringsAsFactors = FALSE
  )
  rownames(meta_chr) <- rownames(meta)

  # Unique observed combinations
  combo_df <- unique(meta_chr[, all_group_vars, drop = FALSE])
  rownames(combo_df) <- NULL

  # Sort order: column_split_by first (outer), then group.by vars (inner)
  sort_vars <- c(
    if (!is.null(column_split_by)) column_split_by,
    group.by
  )
  sort_vars <- unique(sort_vars)   # deduplicate in case split is also in group.by
  for (v in sort_vars) {
    lvl_ord      <- .lvls(v)
    combo_df[[v]] <- factor(combo_df[[v]], levels = lvl_ord)
  }
  combo_df <- combo_df[do.call(order, as.list(combo_df[, sort_vars, drop = FALSE])),
                        , drop = FALSE]
  combo_df[] <- lapply(combo_df, as.character)   # back to character
  rownames(combo_df) <- NULL

  n_combos <- nrow(combo_df)
  message("Building GroupHeatmap: ", n_combos, " column group(s) × ",
          length(features_use), " features...")

  # ── 3. Extract expression data ───────────────────────────────────────────────
  mat_raw <- .get_layer_data(seurat_object, assay = assay, layer = layer,
                              features = features_use)
  mat_raw  <- mat_raw[features_use, , drop = FALSE]

  # Align cells: meta and mat_raw may differ if some cells were filtered
  common_cells <- intersect(rownames(meta_chr), colnames(mat_raw))
  meta_chr <- meta_chr[common_cells, , drop = FALSE]
  mat_raw  <- mat_raw[ , common_cells, drop = FALSE]

  # For each combination row, find matching cells and compute stats
  get_cells <- function(k) {
    idx <- rep(TRUE, nrow(meta_chr))
    for (v in all_group_vars)
      idx <- idx & (meta_chr[[v]] == combo_df[k, v])
    which(idx)
  }

  mean_list <- lapply(seq_len(n_combos), function(k) {
    cells <- get_cells(k)
    if (length(cells) == 0L) rep(NA_real_, length(features_use))
    else Matrix::rowMeans(mat_raw[, cells, drop = FALSE], na.rm = TRUE)
  })
  pct_list <- lapply(seq_len(n_combos), function(k) {
    cells <- get_cells(k)
    if (length(cells) == 0L) rep(0, length(features_use))
    else Matrix::rowSums(mat_raw[, cells, drop = FALSE] > 0) / length(cells)
  })

  mean_mat <- do.call(cbind, mean_list)
  pct_mat  <- do.call(cbind, pct_list)
  rownames(mean_mat) <- features_use
  rownames(pct_mat)  <- features_use

  # ── 4. Column identifiers ────────────────────────────────────────────────────
  # Internal unique IDs (used as matrix/annotation rownames) — always unique
  col_ids <- apply(combo_df[, all_group_vars, drop = FALSE], 1,
                   function(r) paste(r, collapse = " │ "))  # " | " separator
  if (anyDuplicated(col_ids))
    col_ids <- paste0(col_ids, " [", seq_len(n_combos), "]")

  # Display labels shown on heatmap axis (only group.by part, not split)
  col_display <- if (length(group.by) == 1L) {
    combo_df[[group.by]]
  } else {
    apply(combo_df[, group.by, drop = FALSE], 1, paste, collapse = "\n")
  }

  colnames(mean_mat) <- col_ids
  colnames(pct_mat)  <- col_ids

  # ── 5. Scale expression ──────────────────────────────────────────────────────
  expr_mat <- switch(scale_method,
    zscore = {
      sc <- t(scale(t(mean_mat)))
      sc[is.nan(sc)] <- 0
      sc
    },
    minmax = {
      rng   <- t(apply(mean_mat, 1, range, na.rm = TRUE))
      denom <- rng[, 2] - rng[, 1]
      denom[denom == 0] <- 1
      (mean_mat - rng[, 1]) / denom
    },
    none = mean_mat
  )

  # ── 6. Colour scale ──────────────────────────────────────────────────────────
  if (!is.null(color_scale)) {
    col_fun <- color_scale
  } else {
    pal_colors <- tryCatch(
      rev(RColorBrewer::brewer.pal(9, palette)),
      error = function(e) {
        message("scSidekick: Could not load palette '", palette,
                "'. Falling back to blue-white-red.")
        c("#007dd1", "#FFFFFF", "#ab3000")
      }
    )
    col_fun <- circlize::colorRamp2(
      seq(scale_limits[1], scale_limits[2], length.out = length(pal_colors)),
      pal_colors
    )
  }

  # ── 7. Column split factor ───────────────────────────────────────────────────
  col_split <- NULL
  if (!is.null(column_split_by)) {
    present   <- sort(unique(combo_df[[column_split_by]]))
    lvl_order <- .lvls(column_split_by)
    lvl_order <- lvl_order[lvl_order %in% present]
    col_split <- factor(combo_df[[column_split_by]], levels = lvl_order)
  }

  # ── 8. Column annotation ─────────────────────────────────────────────────────
  annot_df     <- data.frame(row.names = col_ids)
  annot_colors <- list()

  for (v in group.by) {
    annot_df[[v]] <- factor(combo_df[[v]], levels = .lvls(v))

    # Resolve per-variable colors
    v_cols <- if (is.list(group_colors)) {
      group_colors[[v]]               # user supplied named list
    } else if (!is.null(group_colors) && length(group.by) == 1L) {
      group_colors                    # backward compat: single named vector
    } else {
      NULL
    }
    v_cols <- v_cols %||% .nk_colors(seurat_object, v)
    if (!is.null(v_cols)) annot_colors[[v]] <- v_cols
  }

  # Add column_split_by as a separate annotation bar if it's not already in group.by
  if (!is.null(column_split_by) && !column_split_by %in% group.by) {
    annot_df[[column_split_by]] <- col_split
    cs_cols <- column_split_colors %||% .nk_colors(seurat_object, column_split_by)
    if (!is.null(cs_cols)) annot_colors[[column_split_by]] <- cs_cols
  }

  top_annot <- ComplexHeatmap::columnAnnotation(
    df                   = annot_df,
    col                  = annot_colors,
    show_annotation_name = TRUE,
    show_legend          = TRUE
  )

  # ── 9. Row (feature) annotation ──────────────────────────────────────────────
  left_annot <- NULL
  if (!is.null(feature_split) && isTRUE(textbox_labels)) {

    has_textbox <- tryCatch(
      utils::packageVersion("ComplexHeatmap") >= "2.13.0",
      error = function(e) FALSE
    )

    if (has_textbox) {
      text_list <- stats::setNames(
        as.list(as.character(unique(feature_split))),
        as.character(unique(feature_split))
      )
      tb <- ComplexHeatmap::anno_textbox(
        feature_split, text_list,
        side          = "left",
        background_gp = grid::gpar(fill = "#fffff9", color = "black"),
        gp            = grid::gpar(fontsize = 8)
      )
      left_annot <- ComplexHeatmap::rowAnnotation(
        textbox              = tb,
        show_annotation_name = FALSE
      )
    } else {
      message(
        "scSidekick: anno_textbox requires ComplexHeatmap >= 2.13.0 (Bioconductor 3.17). ",
        "Detected version ", utils::packageVersion("ComplexHeatmap"), ". ",
        "Falling back to simple row labels. ",
        "Upgrade: BiocManager::install('ComplexHeatmap')"
      )
      left_annot <- ComplexHeatmap::rowAnnotation(
        Category = ComplexHeatmap::anno_text(
          as.character(feature_split),
          gp = grid::gpar(fontsize = 8)
        ),
        show_annotation_name = FALSE
      )
    }
  }

  # ── 10. Dot overlay (cell_fun) ────────────────────────────────────────────────
  # Capture in closure to avoid environment lookup issues
  .pct_mat  <- pct_mat
  .pct_cut  <- pct_cutoff
  .dot_size <- dot_size
  .add_bg   <- add_bg

  cell_fn <- if (isTRUE(add_dot)) {
    function(j, i, x, y, w, h, fill) {
      perc <- ComplexHeatmap::pindex(.pct_mat, i, j)
      if (!is.na(perc) && perc >= .pct_cut) {
        if (.add_bg) {
          # White halo disc — slightly larger than the actual dot so it
          # creates a visible ring that lifts the dot off dark backgrounds.
          # We add a fixed 0.8 mm "halo ring" on top of the scaled dot size.
          grid::grid.points(
            x, y, pch = 21,
            size = .dot_size * perc + grid::unit(0.8, "mm"),
            gp   = grid::gpar(fill = "white", col = NA)
          )
        }
        grid::grid.points(
          x, y, pch = 21,
          size = .dot_size * perc,
          gp   = grid::gpar(col = "black", lwd = 0.5, fill = fill)
        )
      }
    }
  } else NULL

  # ── 11. Dot legend ────────────────────────────────────────────────────────────
  dot_lgd <- if (isTRUE(add_dot)) {
    pct_breaks <- seq(0.2, 1, length.out = 5)
    ComplexHeatmap::Legend(
      labels      = paste0(round(pct_breaks * 100), "%"),
      title       = "Percent\nExpressing",
      type        = "points",
      pch         = 21,
      size        = dot_size * pct_breaks,
      grid_height = dot_size * pct_breaks * 0.8,
      grid_width  = dot_size,
      legend_gp   = grid::gpar(fill = "grey40"),
      border      = FALSE,
      background  = "transparent",
      direction   = "vertical"
    )
  } else NULL

  # ── 12. Build Heatmap object ─────────────────────────────────────────────────
  message("  Rendering heatmap...")

  default_ht_args <- list(
    matrix               = as.matrix(expr_mat),
    name                 = scale_method,
    col                  = col_fun,
    top_annotation       = top_annot,
    left_annotation      = left_annot,
    row_split            = if (!is.null(feature_split)) feature_split else NULL,
    column_split         = col_split,
    column_labels        = col_display,    # display-only labels (no split prefix)
    cluster_rows         = cluster_rows,
    cluster_row_slices   = cluster_row_slices,
    cluster_columns      = cluster_columns,
    cell_fun             = cell_fn,
    width                = width,
    height               = height,
    show_row_names       = TRUE,
    row_names_side       = "right",
    row_names_gp         = grid::gpar(fontsize = 8),
    show_row_dend        = FALSE,
    show_column_names    = TRUE,
    column_names_side    = "top",
    row_title            = NULL,
    column_title         = NULL,
    heatmap_legend_param = list(
      title     = paste0("Mean\n", scale_method),
      direction = "vertical"
    )
  )

  ht_args <- modifyList(default_ht_args, heatmap_params)
  ht      <- do.call(ComplexHeatmap::Heatmap, ht_args)

  # ── 13. Draw and / or save ───────────────────────────────────────────────────
  lgd_list <- Filter(Negate(is.null), list(dot_lgd))

  .do_draw <- function() {
    ComplexHeatmap::draw(ht,
                         heatmap_legend_side = legend_side,
                         heatmap_legend_list = lgd_list)
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    pfx      <- paste(c(object_name, subset_name)[nchar(c(object_name, subset_name)) > 0],
                      collapse = "_")
    pdf_path <- file.path(output_dir,
                           paste0(pfx, if (nchar(pfx) > 0) " " else "",
                                  "GroupHeatmap.pdf"))
    w_in <- tryCatch(as.numeric(grid::convertUnit(width,  "in")), error = function(e) 7)
    h_in <- tryCatch(as.numeric(grid::convertUnit(height, "in")), error = function(e) 9)
    grDevices::pdf(pdf_path, width = w_in + 3, height = h_in + 1)
    .do_draw()
    grDevices::dev.off()
    message("scSidekick: Saved to ", pdf_path)
  }

  if (isTRUE(auto_draw)) .do_draw()

  # ── 14. Return ───────────────────────────────────────────────────────────────
  invisible(list(
    heatmap           = ht,
    expression_matrix = expr_mat,
    percent_matrix    = pct_mat,
    column_groups     = combo_df,       # data.frame of all column combinations
    feature_split     = feature_split,
    dot_legend        = dot_lgd
  ))
}
