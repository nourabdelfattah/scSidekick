# =============================================================================
# scSidekick — Composition Visualizations  (viz_composition.R)
#
# Exported:
#   PlotComposition()  — stacked / percent / dodged bar charts showing
#                        sample-cluster composition in both directions
#   PlotPieUMAP()      — UMAP scatter with pie charts at cluster centroids
#                        showing split.by composition (uses ggforce, stable)
# =============================================================================


# =============================================================================
# PlotComposition
# =============================================================================

#' Cell Composition Bar Charts
#'
#' @description
#' Visualises the composition between two metadata variables (e.g. cluster and
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
#' @param group_colors Named character vector of colours for \code{group.by}
#'   levels. Auto-assigned from PrepObject or \code{Nour_pal} if \code{NULL}.
#' @param split_colors Named character vector of colours for \code{split.by}
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
      panel.border       = ggplot2::element_rect(colour = "black", linewidth = 0.5),
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1,
                                                   face = "bold", colour = "black"),
      axis.text.y        = ggplot2::element_text(colour = "black"),
      axis.title         = ggplot2::element_text(face = "bold"),
      strip.background   = ggplot2::element_rect(fill = "white", colour = "black"),
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
                                             limits = c(0, 100))

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
                                             limits = c(0, 100))

  out_plot <- patchwork::wrap_plots(p1, p2, nrow = 1)

  .save_plot(out_plot, output_dir, object_name, subset_name,
             filename = "composition.pdf",
             w = max(6, (n_grp + n_spl) * 0.45 + 4), h = 5)
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
#' Overlays proportional pie charts on a UMAP (or any 2D reduction) — one pie
#' per level of \code{group.by} — where each slice shows the fraction of cells
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
#' @param bg_color Background scatter point colour. Default \code{"lightgrey"}.
#' @param bg_alpha Background scatter point transparency. Default \code{0.3}.
#' @param pt.size Background scatter point size. Default \code{0.4}.
#' @param split_colors Named character vector of colours for \code{split.by}
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
    bg_color     = "lightgrey",
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
      axis.title    = ggplot2::element_text(face = "bold", colour = "black"),
      axis.text     = ggplot2::element_text(colour = "black"),
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
  print(p)
  invisible(p)
}
