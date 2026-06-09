# =============================================================================
# scSidekick - Composition Visualizations  (viz_composition.R)
#
# Exported:
#   PlotComposition()  - stacked / percent / dodged bar charts showing
#                        sample-cluster composition in both directions
#   PlotPieUMAP()      - UMAP scatter with pie charts at cluster centroids
#                        showing split.by composition (uses ggforce, stable)
# =============================================================================


# =============================================================================
# PlotComposition
# =============================================================================

#' Cell Composition Bar Charts
#'
#' @description
#' Visualizes the composition between two metadata variables (e.g. cluster and
#' sample) as stacked, 100\%-stacked, or dodged bar charts.  By default two
#' complementary panels are produced:
#' \enumerate{
#'   \item \emph{Group → Split}: "What percentage of each cluster comes from
#'     each sample?"
#'   \item \emph{Split → Group}: "How are clusters distributed within each
#'     sample?"
#' }
#'
#' @param seurat_object A Seurat object.
#' @param group.by Primary grouping variable (e.g. cluster, cell type).
#' @param split.by Secondary variable (e.g. sample, condition).
#' @param plot_type One of \code{"percent"} (100\%-stacked, default),
#'   \code{"stacked"} (raw cell counts), or \code{"dodged"} (side-by-side bars).
#' @param both_panels Logical.  If \code{TRUE} (default) both group→split and
#'   split→group panels are shown side by side.  \code{FALSE} returns only the
#'   group→split panel.
#' @param show_counts Logical.  Overlay cell-count labels on bars.  Default
#'   \code{FALSE}.
#' @param group_colors Named character vector of colors for \code{group.by}
#'   levels. Auto-assigned from PrepObject or \code{Nour_pal} if \code{NULL}.
#' @param split_colors Named character vector of colors for \code{split.by}
#'   levels. Auto-assigned from PrepObject or \code{Nour_pal} if \code{NULL}.
#' @param output_dir Directory to save a PDF. \code{NULL} = no save.
#' @param object_name Label prefix for output file names.
#' @param subset_name Optional subset label.
#'
#' @return A \code{ggplot2} or \code{patchwork} object.
#'
#' @export
PlotComposition <- function(
    seurat_object,
    group.by,
    split.by,
    plot_type    = "percent",
    both_panels  = TRUE,
    show_counts  = FALSE,
    group_colors = NULL,
    split_colors = NULL,
    output_dir   = NULL,
    object_name  = "",
    subset_name  = ""
) {
  plot_type <- match.arg(plot_type, c("percent", "stacked", "dodged"))
  meta      <- seurat_object@meta.data

  for (col in c(group.by, split.by)) {
    if (!col %in% colnames(meta))
      stop("'", col, "' not found in seurat_object@meta.data.")
  }

  # ── Factor levels (respect PrepObject order) ─────────────────────────────
  .lvls <- function(col_name) {
    col <- meta[[col_name]]
    if (is.factor(col)) levels(col) else sort(unique(as.character(col)))
  }
  group_levels <- .lvls(group.by)
  split_levels <- .lvls(split.by)

  meta[[group.by]] <- factor(as.character(meta[[group.by]]), levels = group_levels)
  meta[[split.by]] <- factor(as.character(meta[[split.by]]), levels = split_levels)

  # ── Colors ────────────────────────────────────────────────────────────────
  n_grp <- length(group_levels)
  n_spl <- length(split_levels)
  if (is.null(group_colors))
    group_colors <- .nk_colors(seurat_object, group.by) %||%
      stats::setNames(Nour_pal(if (n_grp <= 8) "all" else "spectrum")(n_grp), group_levels)
  if (is.null(split_colors))
    split_colors <- .nk_colors(seurat_object, split.by) %||%
      stats::setNames(Nour_pal("spectrum")(n_spl), split_levels)

  # ── Build contingency table → long data frame ─────────────────────────────
  ct <- as.data.frame(table(meta[[group.by]], meta[[split.by]]),
                       stringsAsFactors = FALSE)
  colnames(ct) <- c("group", "split", "n")
  ct$group <- factor(ct$group, levels = group_levels)
  ct$split <- factor(ct$split, levels = split_levels)

  # Percentages
  ct <- ct |>
    dplyr::group_by(group) |>
    dplyr::mutate(pct_within_group = n / sum(n) * 100) |>
    dplyr::ungroup() |>
    dplyr::group_by(split) |>
    dplyr::mutate(pct_within_split = n / sum(n) * 100) |>
    dplyr::ungroup()

  # Y-axis variable
  y_var    <- switch(plot_type,
    percent = "pct_within_group",
    stacked = "n",
    dodged  = "pct_within_group"
  )
  y_label  <- switch(plot_type,
    percent = "Percentage of cells (%)",
    stacked = "Cell count",
    dodged  = "Percentage of cells (%)"
  )
  pos_geom <- if (plot_type == "dodged") "dodge" else "stack"

  # ── Shared theme ──────────────────────────────────────────────────────────
  .comp_theme <- ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(color = "black", linewidth = 0.5),
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1,
                                                   face = "bold", color = "black"),
      axis.text.y        = ggplot2::element_text(color = "black"),
      axis.title         = ggplot2::element_text(face = "bold"),
      strip.background   = ggplot2::element_rect(fill = "white", color = "black"),
      strip.text         = ggplot2::element_text(face = "bold"),
      legend.title       = ggplot2::element_text(face = "bold"),
      plot.title         = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.margin        = ggplot2::unit(c(0.3, 0.5, 0.2, 0.3), "cm")
    )

  # ── Panel 1: group.by on x, fill by split.by ─────────────────────────────
  p1 <- ggplot2::ggplot(ct,
    ggplot2::aes(x    = group,
                 y    = .data[[y_var]],
                 fill = split)) +
    ggplot2::geom_col(position = pos_geom, width = 0.8, color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = split_colors, name = split.by) +
    ggplot2::labs(title = paste("Composition by", group.by),
                  x = group.by, y = y_label) +
    .comp_theme

  if (plot_type == "percent")
    p1 <- p1 + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02)),
                                             limits = c(0, 101))

  if (show_counts)
    p1 <- p1 + ggplot2::geom_text(ggplot2::aes(label = n),
                                    position = if (plot_type == "dodged")
                                      ggplot2::position_dodge(0.8)
                                    else ggplot2::position_stack(vjust = 0.5),
                                    size = 2.5, color = "white", fontface = "bold")

  if (!isTRUE(both_panels)) {
    .save_plot(p1, output_dir, object_name, subset_name,
               filename = "composition.pdf",
               w = max(3, n_grp * 0.5 + 2), h = 5)
    if (!is.null(output_dir)) {
      pfx  <- paste(c(object_name, subset_name)[nchar(c(object_name, subset_name)) > 0], collapse = "_")
      .write_legend_sidecar(
        file.path(output_dir, paste0(pfx, if (nchar(pfx)) " " else "", "composition.pdf")),
        paste0(
          plot_type, " bar chart showing the composition of ", group.by,
          " levels colored by ", split.by, ". ",
          "The x-axis shows each ", group.by, " level; bar fill represents ",
          split.by, " proportions."
        )
      )
    }
    print(p1); return(invisible(p1))
  }

  # ── Panel 2: split.by on x, fill by group.by ─────────────────────────────
  ct2 <- ct
  ct2$x_var <- ct2$split
  ct2$fill_var <- ct2$group
  ct2$y2 <- ct2$pct_within_split

  p2 <- ggplot2::ggplot(ct2,
    ggplot2::aes(x    = x_var,
                 y    = if (plot_type == "stacked") n else pct_within_split,
                 fill = fill_var)) +
    ggplot2::geom_col(position = pos_geom, width = 0.8, color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = group_colors, name = group.by) +
    ggplot2::labs(title = paste("Composition by", split.by),
                  x = split.by, y = y_label) +
    .comp_theme

  if (plot_type == "percent")
    p2 <- p2 + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02)),
                                             limits = c(0, 101))

  out_plot <- patchwork::wrap_plots(p1, p2, nrow = 1)

  .save_plot(out_plot, output_dir, object_name, subset_name,
             filename = "composition.pdf",
             w = max(6, (n_grp + n_spl) * 0.45 + 4), h = 5)
  if (!is.null(output_dir)) {
    pfx  <- paste(c(object_name, subset_name)[nchar(c(object_name, subset_name)) > 0], collapse = "_")
    .write_legend_sidecar(
      file.path(output_dir, paste0(pfx, if (nchar(pfx)) " " else "", "composition.pdf")),
      paste0(
        "Two-panel ", plot_type, " bar chart showing the reciprocal composition of ",
        group.by, " and ", split.by, ". ",
        "Left panel: each ", group.by, " level on the x-axis, bars filled by ",
        split.by, " proportions. ",
        "Right panel: each ", split.by, " level on the x-axis, bars filled by ",
        group.by, " proportions."
      )
    )
  }
  print(out_plot)
  invisible(out_plot)
}


# ── Shared save helper ───────────────────────────────────────────────────────
.save_plot <- function(p, output_dir, object_name, subset_name, filename, w, h) {
  if (is.null(output_dir)) return(invisible(NULL))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pfx  <- paste(c(object_name, subset_name)[nchar(c(object_name, subset_name)) > 0],
                collapse = "_")
  path <- file.path(output_dir, paste0(pfx, if (nchar(pfx)) " " else "", filename))
  pdf(path, width = w, height = h)
  print(p)
  grDevices::dev.off()
  message("Saved: ", basename(path))
}


# =============================================================================
# PlotPieUMAP
# =============================================================================

#' UMAP with Pie Charts at Cluster Centroids
#'
#' @description
#' Overlays proportional pie charts on a UMAP (or any 2D reduction) - one pie
#' per level of \code{group.by} - where each slice shows the fraction of cells
#' from each level of \code{split.by} (e.g. sample composition per cluster).
#'
#' Uses \code{ggforce::geom_arc_bar} for stable pie rendering that works
#' correctly with ggplot2 >= 3.4 (unlike the old \code{coord_polar} approach).
#'
#' @param seurat_object A Seurat object.
#' @param group.by Metadata column defining the groups/clusters. One pie chart
#'   is placed at each group's centroid in the reduction space.
#' @param split.by Metadata column whose level proportions define the pie
#'   slices (e.g. \code{"Sample"} shows how much of each cluster is each
#'   sample).
#' @param reduction Dimensionality reduction to use. Default \code{"umap"}.
#' @param pie_scale Scaling multiplier for pie radius (1 = auto-sized).
#'   Increase for larger pies, decrease if pies overlap. Default \code{1}.
#' @param bg_color Background scatter point color. Default \code{"lightgray"}.
#' @param bg_alpha Background scatter point transparency. Default \code{0.3}.
#' @param pt.size Background scatter point size. Default \code{0.4}.
#' @param split_colors Named character vector of colors for \code{split.by}
#'   levels. Auto-assigned from PrepObject or \code{Nour_pal} if \code{NULL}.
#' @param label Logical. Add cluster labels at centroids. Default \code{TRUE}.
#' @param label_size Font size for cluster labels. Default \code{3}.
#' @param output_dir Directory to save a PDF. \code{NULL} = no save.
#' @param object_name Label prefix for output file names.
#' @param subset_name Optional subset label.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
PlotPieUMAP <- function(
    seurat_object,
    group.by,
    split.by,
    reduction    = "umap",
    pie_scale    = 1,
    bg_color     = "lightgray",
    bg_alpha     = 0.3,
    pt.size      = 0.4,
    split_colors = NULL,
    label        = TRUE,
    label_size   = 3,
    output_dir   = NULL,
    object_name  = "",
    subset_name  = ""
) {
  if (!requireNamespace("ggforce", quietly = TRUE))
    stop("Package 'ggforce' is required for PlotPieUMAP.\n",
         "Install with: install.packages('ggforce')")

  meta <- seurat_object@meta.data
  for (col in c(group.by, split.by)) {
    if (!col %in% colnames(meta))
      stop("'", col, "' not found in seurat_object@meta.data.")
  }

  # ── Extract 2D coordinates ────────────────────────────────────────────────
  red_key <- tolower(reduction)
  if (!red_key %in% names(seurat_object@reductions))
    stop("Reduction '", reduction, "' not found in seurat_object@reductions.")

  coords <- as.data.frame(
    SeuratObject::Embeddings(seurat_object, reduction = red_key)[, 1:2]
  )
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  coords[[group.by]] <- as.character(meta[[group.by]])
  coords[[split.by]] <- as.character(meta[[split.by]])

  # ── Factor levels ─────────────────────────────────────────────────────────
  .lvls <- function(col_name) {
    col <- meta[[col_name]]
    if (is.factor(col)) levels(col) else sort(unique(as.character(col)))
  }
  group_levels <- .lvls(group.by)
  split_levels <- .lvls(split.by)

  # ── Colors ────────────────────────────────────────────────────────────────
  n_spl <- length(split_levels)
  if (is.null(split_colors))
    split_colors <- .nk_colors(seurat_object, split.by) %||%
      stats::setNames(Nour_pal("spectrum")(n_spl), split_levels)

  # ── Auto pie radius ───────────────────────────────────────────────────────
  x_range  <- diff(range(coords$UMAP_1))
  y_range  <- diff(range(coords$UMAP_2))
  umap_sc  <- min(x_range, y_range)
  n_grp    <- length(group_levels)
  auto_r   <- (umap_sc / (2 * sqrt(n_grp))) * 0.55 * pie_scale

  # ── Compute arc segments per group ────────────────────────────────────────
  arc_rows <- lapply(group_levels, function(grp) {
    cells      <- coords[[group.by]] == grp
    cx         <- mean(coords$UMAP_1[cells])
    cy         <- mean(coords$UMAP_2[cells])
    comp       <- table(factor(coords[[split.by]][cells], levels = split_levels))
    fracs      <- as.numeric(comp) / sum(comp)
    end_ang    <- cumsum(fracs) * 2 * pi - pi / 2   # start from 12 o'clock
    start_ang  <- c(-pi / 2, end_ang[-length(end_ang)])

    data.frame(
      x0        = cx,
      y0        = cy,
      r0        = 0,
      r         = auto_r,
      start     = start_ang,
      end       = end_ang,
      fill      = split_levels,
      group_lbl = grp,
      stringsAsFactors = FALSE
    )
  })
  arc_df <- do.call(rbind, arc_rows)
  arc_df$fill <- factor(arc_df$fill, levels = split_levels)

  # Centroids for labels
  centroid_df <- do.call(rbind, lapply(group_levels, function(grp) {
    cells <- coords[[group.by]] == grp
    data.frame(UMAP_1 = mean(coords$UMAP_1[cells]),
               UMAP_2 = mean(coords$UMAP_2[cells]),
               label  = grp)
  }))

  # ── Plot ─────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot() +
    # Background scatter
    ggplot2::geom_point(
      data    = coords,
      ggplot2::aes(x = UMAP_1, y = UMAP_2),
      color   = bg_color,
      alpha   = bg_alpha,
      size    = pt.size,
      stroke  = 0
    ) +
    # Pie charts
    ggforce::geom_arc_bar(
      data = arc_df,
      ggplot2::aes(x0 = x0, y0 = y0, r0 = r0, r = r,
                   start = start, end = end, fill = fill),
      color = "white", linewidth = 0.3
    ) +
    ggplot2::scale_fill_manual(values = split_colors, name = split.by) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste(split.by, "composition per", group.by),
      x = paste0(toupper(red_key), "_1"),
      y = paste0(toupper(red_key), "_2")
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid    = ggplot2::element_blank(),
      axis.title    = ggplot2::element_text(face = "bold", color = "black"),
      axis.text     = ggplot2::element_text(color = "black"),
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.title  = ggplot2::element_text(face = "bold"),
      plot.margin   = ggplot2::unit(c(0.3, 0.3, 0.3, 0.3), "cm")
    )

  # Cluster labels
  if (isTRUE(label)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data         = centroid_df,
        ggplot2::aes(x = UMAP_1, y = UMAP_2, label = label),
        size         = label_size,
        fontface     = "bold",
        color        = "black",
        box.padding  = 0.2,
        max.overlaps = 30,
        segment.size = 0.2,
        show.legend  = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        data    = centroid_df,
        ggplot2::aes(x = UMAP_1, y = UMAP_2, label = label),
        size    = label_size,
        fontface = "bold",
        color   = "black"
      )
    }
  }

  .save_plot(p, output_dir, object_name, subset_name,
             filename = paste0("PieUMAP_", group.by, "_by_", split.by, ".pdf"),
             w = 7, h = 6)
  if (!is.null(output_dir)) {
    pfx      <- paste(c(object_name, subset_name)[nchar(c(object_name, subset_name)) > 0], collapse = "_")
    pie_file <- paste0("PieUMAP_", group.by, "_by_", split.by, ".pdf")
    .write_legend_sidecar(
      file.path(output_dir, paste0(pfx, if (nchar(pfx)) " " else "", pie_file)),
      paste0(
        reduction, " scatter plot with proportional pie charts placed at each ",
        group.by, " centroid. Each pie shows the fraction of cells from each ",
        split.by, " level within that ", group.by, " group. ",
        "Background points are colored ", bg_color, "; pie slices are colored by ",
        split.by, "."
      )
    )
  }
  print(p)
  invisible(p)
}


# =============================================================================
# PlotAlluvial
# =============================================================================

#' Parallel Sets (Alluvial) Plot
#'
#' @description
#' Creates a parallel sets (alluvial) chart - a clean alternative to SCP's
#' alluvials that preserves your own color palettes and theme.
#'
#' Supports **two or more axes** via the \code{axes} parameter.  Passing three
#' or more columns (e.g. \code{axes = c("BraakStage", "ThalPhase", "CERAD")})
#' produces a multi-axis diagram ideal for multiparametric comparisons.
#'
#' Use \code{donor.by} to aggregate to the **donor / patient level** before
#' computing counts - streams then show donor flows rather than cell flows.
#' This is the recommended approach for neuropathology scores, clinical
#' metadata, or any variable that is defined per-individual rather than
#' per-cell.
#'
#' When \code{facet.by} is supplied, one alluvial panel is drawn per level and
#' assembled with \pkg{patchwork}.  Plot margins are auto-expanded based on the
#' longest axis label so text is never clipped.
#'
#' @param data A Seurat object **or** a plain \code{data.frame}/\code{tibble}
#'   containing the metadata columns.
#' @param axes Character vector of **2 or more** column names to use as axes,
#'   from left to right (e.g. \code{c("BraakStage", "ThalPhase", "CERAD")}).
#'   Stream color is determined by the first axis.
#' @param group.by Character. Left-axis column (retained for backwards
#'   compatibility).  Ignored when \code{axes} is supplied.
#' @param split.by Character. Right-axis column (retained for backwards
#'   compatibility).  Ignored when \code{axes} is supplied.
#' @param donor.by Character or \code{NULL}.  Metadata column identifying the
#'   donor / patient (e.g. \code{"Donor.ID"}).  When supplied, metadata is
#'   first deduplicated to one row per donor so stream widths reflect
#'   \emph{donor counts} rather than cell counts - the right choice for
#'   patient-level pathological scores (BraakStage, Thal phase, CERAD, etc.).
#'   Default \code{NULL}.
#' @param facet.by Character or \code{NULL}.  Metadata column to facet by -
#'   one alluvial panel per level (e.g. \code{"Diagnosis"}).  Default
#'   \code{NULL}.
#' @param ncol_facet Integer or \code{NULL}.  Number of panel columns when
#'   \code{facet.by} is used.  Auto-computed (max 3) if \code{NULL}.
#' @param alpha Numeric. Stream transparency. Default \code{0.75}.
#' @param axis_width Numeric. Width of the axis rectangles.  The stream width
#'   is \code{1.5 × axis_width}.  Default \code{0.1}.
#' @param label Logical. Annotate each axis segment with its level name.
#'   Default \code{TRUE}.
#' @param label_size Numeric. Label font size. Default \code{4}.
#' @param label_nudge Numeric. Horizontal gap between axis bar and label.
#'   Default \code{0.1}.
#' @param group_colors Named character vector of colors for the first axis
#'   (streams + left axis).  Auto-assigned if \code{NULL}.
#' @param split_colors Named character vector of colors for the last axis.
#'   Auto-assigned if \code{NULL}.
#' @param axes_colors Named \code{list} of color vectors, keyed by column
#'   name (e.g. \code{list(BraakStage = BraakColors, ThalPhase = ThalColors)}).
#'   Takes precedence over \code{group_colors}/\code{split_colors}.
#' @param output_dir Directory to save a PDF.  \code{NULL} = no save.
#' @param object_name Label prefix for the output filename.
#' @param pdf_width,pdf_height PDF dimensions per panel in inches.  Default
#'   \code{7 x 5.5}.
#'
#' @return A \code{ggplot2} or \code{patchwork} object (invisibly).
#' @export
PlotAlluvial <- function(
    data,
    axes         = NULL,
    group.by     = NULL,
    split.by     = NULL,
    donor.by     = NULL,
    facet.by     = NULL,
    ncol_facet   = NULL,
    alpha        = 0.75,
    axis_width   = 0.1,
    label        = TRUE,
    label_size   = 4,
    label_nudge  = 0.1,
    group_colors = NULL,
    split_colors = NULL,
    axes_colors  = list(),
    output_dir   = NULL,
    object_name  = "",
    pdf_width    = 7,
    pdf_height   = 5.5
) {
  if (!requireNamespace("ggforce", quietly = TRUE))
    stop("Package 'ggforce' is required for PlotAlluvial.\n",
         "Install with: install.packages('ggforce')")

  # ── Resolve axes ───────────────────────────────────────────────────────────
  if (is.null(axes)) {
    if (is.null(group.by) || is.null(split.by))
      stop("Provide 'axes' (a vector of 2+ column names) OR both 'group.by' and 'split.by'.")
    axes <- c(group.by, split.by)
  }
  if (length(axes) < 2L) stop("'axes' must specify at least 2 columns.")
  n_axes <- length(axes)

  # ── Extract metadata ──────────────────────────────────────────────────────
  meta <- if (inherits(data, "Seurat")) data@meta.data else as.data.frame(data)
  for (col in c(axes, donor.by, facet.by))
    if (!is.null(col) && !col %in% colnames(meta))
      stop("'", col, "' not found in metadata.")

  # Drop rows with NA in any axis, donor, or facet column
  drop_cols <- c(axes, donor.by, facet.by)
  drop_cols <- drop_cols[lengths(list(drop_cols)) > 0 & !sapply(drop_cols, is.null)]
  meta <- meta[stats::complete.cases(meta[, drop_cols, drop = FALSE]), , drop = FALSE]

  # ── Donor-level aggregation ────────────────────────────────────────────────
  if (!is.null(donor.by)) {
    meta <- meta[!duplicated(meta[[donor.by]]), , drop = FALSE]
    message("scSidekick: Aggregated to ", nrow(meta), " unique donors via '", donor.by, "'.")
  }
  unit_word <- if (!is.null(donor.by)) "donor" else "cell"

  # ── Factor levels (from full data - consistent across facets) ─────────────
  .lvls <- function(col_name) {
    col <- meta[[col_name]]
    if (is.factor(col)) levels(col) else sort(unique(as.character(col)))
  }
  axes_levels <- stats::setNames(lapply(axes, .lvls), axes)

  # ── Color resolution (per axis) ───────────────────────────────────────────
  colors_per_axis <- lapply(seq_along(axes), function(i) {
    ax   <- axes[i]
    lvls <- axes_levels[[ax]]
    n    <- length(lvls)

    # Priority: axes_colors > group_colors (axis 1) > split_colors (last axis) > auto
    if (!is.null(axes_colors[[ax]]))  return(axes_colors[[ax]][lvls])
    if (i == 1L     && !is.null(group_colors)) return(group_colors[lvls])
    if (i == n_axes && !is.null(split_colors)) return(split_colors[lvls])
    auto <- if (inherits(data, "Seurat")) .nk_colors(data, ax) else NULL
    if (!is.null(auto)) return(auto[lvls])
    stats::setNames(Nour_pal(if (n <= 8) "all" else "spectrum")(n), lvls)
  })
  names(colors_per_axis) <- axes
  # unlist(unname(...)) is the ONLY form that gives bare level names ("AD", "0")
  # rather than list-path-prefixed names ("Dementia.AD.AD", "Braak.0").
  # Both unlist(use.names=TRUE) and do.call(c, ...) prefix with the list element
  # name when the list itself has names - breaking scale_fill_manual lookups.
  # Stripping outer names with unname() first lets unlist keep just the inner names.
  # When the same value appears in multiple axes, the first axis's color wins
  # (consistent with ggforce treating identical values as the same node).
  all_colors <- unlist(unname(colors_per_axis))

  # ── Inner panel builder (one metadata slice → one ggplot) ─────────────────
  .one_panel <- function(meta_sub, panel_title = NULL) {
    # Levels present in this slice (may be fewer than full-data levels)
    axes_lvls_sub <- lapply(axes, function(ax) {
      meta_sub[[ax]] <- factor(as.character(meta_sub[[ax]]), levels = axes_levels[[ax]])
      intersect(axes_levels[[ax]], unique(as.character(meta_sub[[ax]])))
    })
    names(axes_lvls_sub) <- axes
    for (ax in axes)
      meta_sub[[ax]] <- factor(as.character(meta_sub[[ax]]), levels = axes_lvls_sub[[ax]])

    # Aggregate - group by ALL axes at once
    df <- meta_sub |>
      dplyr::group_by(dplyr::across(dplyr::all_of(axes))) |>
      dplyr::tally() |>
      dplyr::ungroup() |>
      as.data.frame()

    if (nrow(df) == 0L) return(NULL)

    # Build ggforce data (supports 2+ axes via 1:n_axes)
    pdata <- ggforce::gather_set_data(df, seq_len(n_axes))
    pdata$x <- factor(pdata$x, levels = unique(pdata$x))
    # y levels: axes[1] levels first, then axes[2], ..., axes[n_axes].
    # unique() is required: if the same value (e.g. "0") appears in multiple axes,
    # unlist() produces duplicates which factor() rejects.  ggforce treats identical
    # values across axes as the SAME node - visually and semantically fine when the
    # user's axes are independent ordinal scores (Braak, Thal, CERAD).
    pdata$y <- factor(pdata$y,
                      levels = unique(unlist(axes_lvls_sub, use.names = FALSE)))

    # ── Label positioning: left for axis-1, right for axis-N, center middle ──
    # stat='parallel_sets_axes' generates rows in axis order, level order within
    dlabels <- do.call(rbind, lapply(seq_along(axes), function(i) {
      n_lvl <- length(axes_lvls_sub[[axes[i]]])
      if (i == 1L)      data.frame(hjust = rep(1,   n_lvl), nudge_x = rep(-label_nudge, n_lvl))
      else if (i == n_axes) data.frame(hjust = rep(0,   n_lvl), nudge_x = rep( label_nudge, n_lvl))
      else                  data.frame(hjust = rep(0.5, n_lvl), nudge_x = rep(0,            n_lvl))
    }))

    # Auto-pad: left margin by first axis, right by last axis
    left_pad  <- max(8, max(nchar(axes_lvls_sub[[axes[1]]]))      * label_size * 0.5 + label_nudge * 20)
    right_pad <- max(8, max(nchar(axes_lvls_sub[[axes[n_axes]]])) * label_size * 0.5 + label_nudge * 20)

    p <- ggplot2::ggplot(pdata, ggplot2::aes(x, id = id, split = y, value = n)) +
      ggforce::geom_parallel_sets(
        ggplot2::aes(fill = .data[[axes[1]]]),
        alpha = alpha, axis.width = axis_width * 1.5
      ) +
      ggforce::geom_parallel_sets_axes(
        ggplot2::aes(fill = y),
        color = "black", axis.width = axis_width
      ) +
      ggplot2::scale_fill_manual(values = all_colors) +
      ggplot2::scale_x_discrete(labels = axes) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position  = "none",
        axis.title       = ggplot2::element_blank(),
        axis.text.x      = ggplot2::element_text(face = "bold", color = "black", size = 13),
        axis.text.y      = ggplot2::element_blank(),
        axis.ticks       = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border     = ggplot2::element_blank(),
        plot.margin      = ggplot2::unit(c(5, right_pad, 5, left_pad), "mm")
      )

    if (!is.null(panel_title))
      p <- p +
        ggplot2::labs(title = panel_title) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12))

    if (isTRUE(label))
      p <- p + ggplot2::geom_text(
        ggplot2::aes(y = n, split = y),
        stat     = ggforce::StatParallelSetsAxes,
        fontface = "bold",
        size     = label_size,
        hjust    = dlabels$hjust,
        nudge_x  = dlabels$nudge_x
      )
    p
  }

  # ── Build plot(s) ─────────────────────────────────────────────────────────
  if (!is.null(facet.by)) {
    facet_lvls <- if (is.factor(meta[[facet.by]])) levels(meta[[facet.by]])
                  else sort(unique(as.character(meta[[facet.by]])))

    plot_list <- Filter(Negate(is.null), lapply(facet_lvls, function(fv)
      .one_panel(meta[as.character(meta[[facet.by]]) == fv, , drop = FALSE],
                 panel_title = fv)
    ))

    ncol_f   <- ncol_facet %||% min(length(plot_list), 3L)
    nrow_f   <- ceiling(length(plot_list) / ncol_f)
    combined <- patchwork::wrap_plots(plot_list, ncol = ncol_f)
    pdf_w    <- pdf_width  * ncol_f
    pdf_h    <- pdf_height * nrow_f
  } else {
    combined <- .one_panel(meta)
    pdf_w    <- pdf_width
    pdf_h    <- pdf_height
  }

  # ── Filename ──────────────────────────────────────────────────────────────
  parts <- c(
    if (nchar(object_name) > 0) object_name,
    paste(axes, collapse = "_"),
    if (!is.null(donor.by)) "Donors",
    if (!is.null(facet.by)) paste0("by_", facet.by),
    "Alluvial"
  )
  fname <- paste(parts, collapse = "_")

  # ── Save + .legend sidecar ────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))
    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(combined)
    grDevices::dev.off()
    message("scSidekick: Saved ", basename(fpath),
            " (", round(pdf_w, 1), " × ", round(pdf_h, 1), " in)")

    .write_legend_sidecar(fpath, paste0(
      "Parallel sets (alluvial) diagram showing the flow of ",
      unit_word, "s across ", paste(axes, collapse = " → "),
      ". Stream width is proportional to ", unit_word, " count. ",
      "Streams are colored by ", axes[1], " level.",
      if (n_axes > 2L)
        paste0(" Middle axes: ", paste(axes[-c(1L, n_axes)], collapse = ", "), ".")
      else "",
      if (!is.null(donor.by))
        paste0(" Aggregated to one row per unique ", donor.by, " before counting.")
      else "",
      if (!is.null(facet.by))
        paste0(" Separate panels shown for each level of ", facet.by, ".")
      else "",
      if (nchar(object_name) > 0) paste0(" Dataset: ", object_name, ".") else ""
    ))
  }

  print(combined)
  invisible(combined)
}


# =============================================================================
# PlotRose
# =============================================================================

#' Rose (Polar Bar) Composition Chart
#'
#' @description
#' Visualizes the composition of \code{stat.by} within each level of
#' \code{group.by} as a polar bar (rose) chart.  Each petal of the rose
#' corresponds to a \code{group.by} level; the petal is split by color into
#' \code{stat.by} proportions or counts.
#'
#' When \code{facet.by} is supplied, native \code{facet_wrap} is used so each
#' panel shows its own rose.  Plot margins are auto-expanded based on the
#' longest \code{group.by} label.
#'
#' @param data A Seurat object or a plain \code{data.frame}/\code{tibble}.
#' @param stat.by Character. Variable whose levels become the fill color
#'   (e.g. \code{"Sample"} or \code{"CERAD"}).
#' @param group.by Character. Variable that defines the x-axis petals before
#'   polar transformation (e.g. \code{"Assignment"} or \code{"BraakStage"}).
#' @param donor.by Character or \code{NULL}.  Column identifying the donor /
#'   patient.  When supplied, metadata is first deduplicated to one row per
#'   donor so petal sizes reflect \emph{donor counts} - appropriate for
#'   patient-level neuropathological or clinical variables.  Default
#'   \code{NULL}.
#' @param facet.by Character or \code{NULL}.  Metadata column to facet by -
#'   one rose panel per level.  Default \code{NULL}.
#' @param ncol_facet Integer or \code{NULL}.  Number of facet columns when
#'   \code{facet.by} is used.  Auto-computed (max 3) if \code{NULL}.
#' @param stat_type Character.  \code{"percent"} (default) plots the fraction
#'   of \code{stat.by} within each \code{group.by} level;
#'   \code{"count"} plots raw counts.
#' @param colors Named character vector for \code{stat.by} levels.
#'   Auto-assigned if \code{NULL}.
#' @param alpha Numeric. Fill transparency. Default \code{0.85}.
#' @param label Logical. Show proportion/count labels on each segment.
#'   Default \code{TRUE}.
#' @param label_size Numeric. Label text size. Default \code{2.5}.
#' @param min_label_pct Numeric. Minimum percentage below which labels are
#'   suppressed (avoids clutter in tiny slices). Default \code{2}.
#' @param keep_empty Logical. Include \code{group.by × stat.by} combinations
#'   with zero observations. Default \code{FALSE}.
#' @param output_dir Directory to save a PDF.  \code{NULL} = no save.
#' @param object_name Label prefix for the output filename.
#' @param pdf_width,pdf_height PDF dimensions per panel in inches.  Scales
#'   with number of facet panels.  Default \code{7 x 7}.
#'
#' @return A \code{ggplot2} object (invisibly).
#' @export
PlotRose <- function(
    data,
    stat.by,
    group.by,
    donor.by       = NULL,
    facet.by       = NULL,
    ncol_facet     = NULL,
    stat_type      = "percent",
    colors         = NULL,
    alpha          = 0.85,
    label          = TRUE,
    label_size     = 2.5,
    min_label_pct  = 2,
    keep_empty     = FALSE,
    output_dir     = NULL,
    object_name    = "",
    pdf_width      = 7,
    pdf_height     = 7
) {
  stat_type <- match.arg(stat_type, c("percent", "count"))

  # ── Metadata ──────────────────────────────────────────────────────────────
  meta <- if (inherits(data, "Seurat")) data@meta.data else as.data.frame(data)
  for (col in c(stat.by, group.by, donor.by, facet.by))
    if (!is.null(col) && !col %in% colnames(meta))
      stop("'", col, "' not found in metadata.")

  keep_cols <- c(stat.by, group.by, donor.by, facet.by)
  meta <- meta[stats::complete.cases(meta[, keep_cols[!sapply(keep_cols, is.null)],
                                           drop = FALSE]), , drop = FALSE]

  # ── Donor-level aggregation ────────────────────────────────────────────────
  if (!is.null(donor.by)) {
    meta <- meta[!duplicated(meta[[donor.by]]), , drop = FALSE]
    message("scSidekick: Aggregated to ", nrow(meta), " unique donors via '", donor.by, "'.")
  }
  unit_word <- if (!is.null(donor.by)) "donor" else "cell"

  .lvls <- function(col_name) {
    col <- meta[[col_name]]
    if (is.factor(col)) levels(col) else sort(unique(as.character(col)))
  }
  stat_levels  <- .lvls(stat.by)
  group_levels <- .lvls(group.by)
  n_stat <- length(stat_levels)

  # ── Colors ────────────────────────────────────────────────────────────────
  if (is.null(colors)) {
    colors <- if (inherits(data, "Seurat")) .nk_colors(data, stat.by) else NULL
    if (is.null(colors))
      colors <- stats::setNames(
        Nour_pal(if (n_stat <= 8) "all" else "spectrum")(n_stat), stat_levels)
  }

  meta[[stat.by]]  <- factor(as.character(meta[[stat.by]]),  levels = stat_levels)
  meta[[group.by]] <- factor(as.character(meta[[group.by]]), levels = group_levels)

  # ── Compute counts + proportions ──────────────────────────────────────────
  # Include facet.by in the grouping if provided
  agg_vars <- c(if (!is.null(facet.by)) facet.by, group.by, stat.by)
  pct_vars <- c(if (!is.null(facet.by)) facet.by, group.by)

  ct <- meta |>
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_vars))) |>
    dplyr::tally() |>
    dplyr::ungroup()

  if (isTRUE(keep_empty)) {
    combo_list <- c(
      if (!is.null(facet.by)) list(.lvls(facet.by)) else list(),
      list(factor(group_levels, levels = group_levels),
           factor(stat_levels,  levels = stat_levels))
    )
    names(combo_list) <- agg_vars
    all_combos <- do.call(expand.grid, c(combo_list, list(stringsAsFactors = FALSE)))
    ct <- dplyr::left_join(all_combos, ct, by = agg_vars) |>
      dplyr::mutate(n = ifelse(is.na(n), 0L, n))
  }

  ct <- ct |>
    dplyr::group_by(dplyr::across(dplyr::all_of(pct_vars))) |>
    dplyr::mutate(pct = n / sum(n) * 100) |>
    dplyr::ungroup()

  y_col <- if (stat_type == "percent") "pct" else "n"
  y_lab <- if (stat_type == "percent") "Percentage (%)" else "Cell count"

  # ── Auto-pad margins for outermost axis.text.x labels ─────────────────────
  # polar coord: labels radiate from the ring; generous uniform padding prevents clipping
  rose_pad <- max(10, max(nchar(group_levels)) * label_size * 0.4)

  # ── Plot ──────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(ct, ggplot2::aes(
      x    = .data[[group.by]],
      y    = .data[[y_col]],
      fill = .data[[stat.by]]
    )) +
    ggplot2::geom_bar(
      stat = "identity", position = "stack",
      color = "white", linewidth = 0.2, alpha = alpha
    ) +
    ggplot2::scale_fill_manual(values = colors, name = stat.by) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::labs(
      title = paste(stat.by, "composition by", group.by),
      x = NULL, y = y_lab
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      legend.title  = ggplot2::element_text(face = "bold"),
      legend.text   = ggplot2::element_text(size = 8),
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.x   = ggplot2::element_text(face = "bold", color = "black", size = 8),
      strip.text    = ggplot2::element_text(face = "bold", size = 10),
      plot.margin   = ggplot2::unit(rep(rose_pad, 4), "mm")
    )

  # ── Faceting ──────────────────────────────────────────────────────────────
  n_facets <- if (!is.null(facet.by)) length(unique(as.character(meta[[facet.by]]))) else 1L
  ncol_f   <- ncol_facet %||% min(n_facets, 3L)

  if (!is.null(facet.by))
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet.by]]), ncol = ncol_f)

  # ── Labels on segments ────────────────────────────────────────────────────
  if (isTRUE(label)) {
    grp_vars_lab <- c(if (!is.null(facet.by)) facet.by, group.by)
    label_df <- ct |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_vars_lab))) |>
      dplyr::arrange(dplyr::desc(.data[[stat.by]]), .by_group = TRUE) |>
      dplyr::mutate(
        y_cum = cumsum(.data[[y_col]]),
        y_mid = y_cum - .data[[y_col]] / 2,
        lab   = if (stat_type == "percent")
          paste0(round(pct, 1), "%") else as.character(n)
      ) |>
      dplyr::ungroup() |>
      dplyr::filter(pct >= min_label_pct)

    p <- p + ggplot2::geom_text(
      data        = label_df,
      ggplot2::aes(x = .data[[group.by]], y = y_mid, label = lab),
      size        = label_size,
      color      = "black",
      fontface    = "bold",
      inherit.aes = FALSE
    )
  }

  # ── PDF dimensions ─────────────────────────────────────────────────────────
  nrow_f <- ceiling(n_facets / ncol_f)
  pdf_w  <- pdf_width  * ncol_f
  pdf_h  <- pdf_height * nrow_f

  # ── Filename ──────────────────────────────────────────────────────────────
  parts <- c(
    if (nchar(object_name) > 0) object_name,
    stat.by, group.by,
    if (!is.null(donor.by)) "Donors",
    if (!is.null(facet.by)) paste0("by_", facet.by),
    "Rose"
  )
  fname <- paste(parts, collapse = "_")

  # ── Save + .legend sidecar ────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))
    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(p)
    grDevices::dev.off()
    message("scSidekick: Saved ", basename(fpath),
            " (", round(pdf_w, 1), " × ", round(pdf_h, 1), " in)")

    .write_legend_sidecar(fpath, paste0(
      if (stat_type == "percent") "Percentage" else "Count",
      " rose (polar bar) chart showing ", stat.by, " composition within each ",
      group.by, " level. Each petal represents one ", group.by,
      " level, stacked and filled by ", stat.by, " proportion. ",
      "Values reflect ", unit_word, " counts.",
      if (!is.null(donor.by))
        paste0(" Aggregated to one row per unique ", donor.by, " before counting.")
      else "",
      if (!is.null(facet.by))
        paste0(" Separate panels shown for each level of ", facet.by, ".")
      else "",
      if (nchar(object_name) > 0) paste0(" Dataset: ", object_name, ".") else ""
    ))
  }

  print(p)
  invisible(p)
}


# =============================================================================
# PlotChord
# =============================================================================

#' Chord Diagram Between Two Categorical Variables
#'
#' @description
#' Draws a chord diagram showing the connectivity between two metadata
#' variables (e.g. cell-type assignment and sample).  Each arc represents the
#' number of cells shared between a level of \code{var1} and a level of
#' \code{var2}.  Uses \pkg{circlize}.
#'
#' Because \pkg{circlize} draws to the active graphics device (base R
#' graphics), this function returns the contingency matrix invisibly rather
#' than a \code{ggplot2} object.
#'
#' When \code{facet.by} is supplied, chord diagrams are arranged
#' \strong{side-by-side in a grid} on a single PDF page using
#' \code{par(fig = ..., new = TRUE)}.  The PDF width/height scale
#' automatically with the number of panels.  The filename includes
#' \code{"by_{facet.by}"}.
#'
#' \code{canvas_padding} expands the circlize canvas so long sector labels are
#' never clipped - increase it if you still see clipping.
#'
#' @param data A Seurat object or a plain \code{data.frame}/\code{tibble}.
#' @param var1 Character. First variable (e.g. \code{"BraakStage"}).  Sectors
#'   for \code{var1} levels appear first in the diagram.
#' @param var2 Character. Second variable (e.g. \code{"ThalPhase"}).
#' @param donor.by Character or \code{NULL}.  Column identifying the donor /
#'   patient.  When supplied, metadata is deduplicated to one row per donor so
#'   arc widths reflect \emph{donor counts} - the right choice for any
#'   patient-level pathological or clinical variable.  Default \code{NULL}.
#' @param facet.by Character or \code{NULL}.  Metadata column to facet by -
#'   one chord diagram per level, arranged side-by-side in a grid.  Default
#'   \code{NULL}.
#' @param ncol_facet Integer or \code{NULL}.  Number of columns in the chord
#'   grid when \code{facet.by} is used.  Auto-computed (max 3) if \code{NULL}.
#' @param canvas_padding Numeric.  Extra space added around the chord circle
#'   (in canvas units) to prevent sector labels from being clipped.  The
#'   canvas runs from \code{-1 - padding} to \code{1 + padding} on each axis.
#'   Default \code{0.3}.  Increase to e.g. \code{0.5} for very long labels.
#' @param var1_colors Named character vector for \code{var1} levels.
#'   Auto-assigned if \code{NULL}.
#' @param var2_colors Named character vector for \code{var2} levels.
#'   Auto-assigned if \code{NULL}.
#' @param alpha Numeric. Chord fill opacity (0-1).  Default \code{0.8}.
#' @param label_cex Numeric. Sector label font size multiplier.  Default
#'   \code{0.8}.
#' @param directional Logical.  Draw directional arrows from \code{var1} to
#'   \code{var2}.  Default \code{FALSE}.
#' @param output_dir Directory to save a PDF.  \code{NULL} = no save.
#' @param object_name Label prefix for the output filename.
#' @param pdf_width,pdf_height PDF dimensions \emph{per panel} in inches.
#'   Total PDF width = \code{pdf_width × ncol}; total height = \code{pdf_height
#'   × nrow}.  Default \code{7 x 7}.
#'
#' @return The contingency matrix, or a named list of matrices when
#'   \code{facet.by} is used (invisibly).
#' @export
PlotChord <- function(
    data,
    var1,
    var2,
    donor.by        = NULL,
    facet.by        = NULL,
    ncol_facet      = NULL,
    canvas_padding  = 0.3,
    var1_colors     = NULL,
    var2_colors     = NULL,
    alpha           = 0.8,
    label_cex       = 0.8,
    directional     = FALSE,
    output_dir      = NULL,
    object_name     = "",
    pdf_width       = 7,
    pdf_height      = 7
) {
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Package 'circlize' is required for PlotChord.\n",
         "Install with: install.packages('circlize')")

  # ── Metadata ──────────────────────────────────────────────────────────────
  meta <- if (inherits(data, "Seurat")) data@meta.data else as.data.frame(data)
  for (col in c(var1, var2, facet.by))
    if (!is.null(col) && !col %in% colnames(meta))
      stop("'", col, "' not found in metadata.")

  keep_cols <- c(var1, var2, donor.by, facet.by)
  meta <- meta[stats::complete.cases(
    meta[, keep_cols[!sapply(keep_cols, is.null)], drop = FALSE]), , drop = FALSE]

  # ── Donor-level aggregation ────────────────────────────────────────────────
  if (!is.null(donor.by)) {
    meta <- meta[!duplicated(meta[[donor.by]]), , drop = FALSE]
    message("scSidekick: Aggregated to ", nrow(meta), " unique donors via '", donor.by, "'.")
  }
  unit_word <- if (!is.null(donor.by)) "donor" else "cell"

  .lvls <- function(col_name) {
    col <- meta[[col_name]]
    if (is.factor(col)) levels(col) else sort(unique(as.character(col)))
  }
  v1_levels <- .lvls(var1)
  v2_levels <- .lvls(var2)
  n_v1 <- length(v1_levels)
  n_v2 <- length(v2_levels)

  # ── Colors ────────────────────────────────────────────────────────────────
  if (is.null(var1_colors)) {
    var1_colors <- if (inherits(data, "Seurat")) .nk_colors(data, var1) else NULL
    if (is.null(var1_colors))
      var1_colors <- stats::setNames(
        Nour_pal(if (n_v1 <= 8) "all" else "spectrum")(n_v1), v1_levels)
  }
  if (is.null(var2_colors)) {
    var2_colors <- if (inherits(data, "Seurat")) .nk_colors(data, var2) else NULL
    if (is.null(var2_colors))
      var2_colors <- stats::setNames(Nour_pal("spectrum")(n_v2), v2_levels)
  }
  all_colors <- c(var1_colors[v1_levels], var2_colors[v2_levels])

  # Canvas limits - expanding beyond ±1 gives labels room on all sides
  clim <- c(-1 - canvas_padding, 1 + canvas_padding)

  # ── Inner draw function for one matrix ────────────────────────────────────
  # Called once per mfrow cell; par(mfrow) is set by .draw_grid before any
  # circlize calls so each sub-panel gets its own properly-sized plot region.
  .draw_one <- function(mat_one, panel_title = NULL) {
    # Drop zero-sum rows/columns: they create zero-width sectors that circlize
    # cannot render, triggering "not enough space for cells" errors.
    mat_one <- mat_one[rowSums(mat_one) > 0, colSums(mat_one) > 0, drop = FALSE]
    if (nrow(mat_one) == 0L || ncol(mat_one) == 0L) {
      graphics::plot.new()
      graphics::plot.window(xlim = clim, ylim = clim)
      lbl <- if (!is.null(panel_title) && nzchar(panel_title))
        paste0(panel_title, "\n(no data)") else "(no data)"
      graphics::text(0, 0, lbl, cex = 1.0, adj = c(0.5, 0.5))
      return(invisible(NULL))
    }
    circlize::circos.clear()
    circlize::circos.par(canvas.xlim = clim, canvas.ylim = clim)
    circlize::chordDiagram(
      mat_one,
      grid.col          = all_colors,
      transparency      = 1 - alpha,
      directional       = if (directional) 1L else 0L,
      direction.type    = if (directional) c("diffHeight", "arrows") else "diffHeight",
      annotationTrack   = "grid",
      preAllocateTracks = list(track.height = 0.1)
    )
    circlize::circos.trackPlotRegion(
      track.index = 1,
      panel.fun   = function(x, y) {
        sector.name <- circlize::get.cell.meta.data("sector.index")
        xlim        <- circlize::get.cell.meta.data("xlim")
        ylim        <- circlize::get.cell.meta.data("ylim")
        circlize::circos.text(
          mean(xlim), ylim[1] + 0.1, sector.name,
          facing = "clockwise", niceFacing = TRUE,
          adj = c(0, 0.5), cex = label_cex
        )
      },
      bg.border = NA
    )
    circlize::circos.clear()
    # Panel title in the top margin (mtext works in base-R space after circos.clear)
    if (!is.null(panel_title) && nzchar(panel_title))
      graphics::mtext(panel_title, side = 3, line = 0.3, cex = 1.0, font = 2)
  }

  # ── Build contingency matrices ─────────────────────────────────────────────
  if (!is.null(facet.by)) {
    facet_lvls <- if (is.factor(meta[[facet.by]])) levels(meta[[facet.by]])
                  else sort(unique(as.character(meta[[facet.by]])))
    mat_list <- stats::setNames(lapply(facet_lvls, function(fv) {
      meta_f <- meta[as.character(meta[[facet.by]]) == fv, , drop = FALSE]
      table(factor(meta_f[[var1]], levels = v1_levels),
            factor(meta_f[[var2]], levels = v2_levels))
    }), facet_lvls)
  } else {
    mat_list <- list(table(factor(meta[[var1]], levels = v1_levels),
                           factor(meta[[var2]], levels = v2_levels)))
  }

  # ── Grid layout helpers ────────────────────────────────────────────────────
  n_plots <- length(mat_list)
  ncol_f  <- if (n_plots == 1L) 1L else (ncol_facet %||% min(n_plots, 3L))
  nrow_f  <- ceiling(n_plots / ncol_f)

  # Draw all diagrams in a grid using par(mfrow=...).
  # par(fig=..., new=TRUE) was the previous approach but it doesn't reset
  # margins per cell, so circlize inherits the full default margins in a
  # shrunken sub-region and runs out of space even for simple diagrams.
  # par(mfrow) handles per-cell margin accounting correctly.
  .draw_grid <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    # Small margins: circlize needs room for the circle, not axis labels.
    # 1.5 lines on top leaves space for the panel title drawn via mtext().
    graphics::par(mfrow = c(nrow_f, ncol_f), mar = c(1, 1, 1.5, 1))

    for (i in seq_len(n_plots)) {
      .draw_one(mat_list[[i]],
                panel_title = if (!is.null(facet.by)) names(mat_list)[i] else NULL)
    }
    # Fill any trailing empty cells so the mfrow grid is complete
    remainder <- ncol_f * nrow_f - n_plots
    for (i in seq_len(remainder)) graphics::plot.new()
  }

  # ── Filename ──────────────────────────────────────────────────────────────
  parts <- c(
    if (nchar(object_name) > 0) object_name,
    var1, var2,
    if (!is.null(donor.by)) "Donors",
    if (!is.null(facet.by)) paste0("by_", facet.by),
    "Chord"
  )
  fname <- paste(parts, collapse = "_")

  # ── Save + .legend sidecar ────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))

    # Total PDF = per-panel dims × grid shape
    grDevices::pdf(fpath,
                   width  = pdf_width  * ncol_f,
                   height = pdf_height * nrow_f)
    .draw_grid()
    grDevices::dev.off()

    message("scSidekick: Saved ", basename(fpath),
            " (", pdf_width * ncol_f, " × ", pdf_height * nrow_f, " in",
            if (n_plots > 1L) paste0(", ", ncol_f, "×", nrow_f, " grid") else "",
            ")")

    .write_legend_sidecar(fpath, paste0(
      "Chord diagram showing the distribution of ", unit_word, "s between ",
      var1, " and ", var2, ". Arc width is proportional to ", unit_word, " count.",
      if (!is.null(donor.by))
        paste0(" Aggregated to one row per unique ", donor.by, " before counting.")
      else "",
      if (!is.null(facet.by))
        paste0(" Separate diagrams arranged in a ", ncol_f, "-column grid,",
               " one per level of ", facet.by, ".")
      else "",
      if (nchar(object_name) > 0) paste0(" Dataset: ", object_name, ".") else ""
    ))
  }

  # ── Draw to active device ─────────────────────────────────────────────────
  .draw_grid()

  if (!is.null(facet.by)) invisible(mat_list) else invisible(mat_list[[1L]])
}
