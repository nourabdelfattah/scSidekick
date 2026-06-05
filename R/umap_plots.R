# =============================================================================
# scSidekick UMAP / Dim-plot visualization  (umap_plots.R)
#
# PlotDimPlots()       — unified UMAP grid with optional composition trend
# PlotGridUMAP()       — deprecated thin wrapper → PlotDimPlots
# PlotTrendAndUMAP()   — deprecated thin wrapper → PlotDimPlots(show_trend=TRUE)
#
# Layout:
#   • UMAP: per-row ggarrange grids (no empty panels; shorter rows padded
#     with theme_void() blanks that have no visible border or fill)
#   • Trend: one ggplot per row group, stacked with ggarrange (nrow=n_rows)
#   • BOTH stacks use ggarrange(nrow=n_rows, ncol=1) → each row gets the
#     same fraction of height → corresponding rows are automatically aligned
#   • Legend placed below the combined figure
# =============================================================================

# Internal: case-insensitive Embeddings lookup
.resolve_emb <- function(seurat_object, reduction) {
  avail <- SeuratObject::Reductions(seurat_object)
  hit   <- avail[tolower(avail) == tolower(reduction)]
  if (length(hit) == 0)
    stop("Reduction '", reduction, "' not found. Available: ",
         paste(avail, collapse = ", "))
  SeuratObject::Embeddings(seurat_object, hit[1])
}

# Internal: invisible blank filler panel (no border, no background)
.blank_panel <- function() ggplot2::ggplot() + ggplot2::theme_void()


# =============================================================================
# PlotDimPlots  — main function
# =============================================================================

#' UMAP grid with optional composition trend panel
#'
#' Produces a UMAP grid where **columns** correspond to \code{split.by}
#' levels and **rows** (optionally) to \code{row.by} levels.  No empty
#' panels: samples missing from a row are omitted rather than shown as blank
#' bordered boxes.  The optional composition trend panel
#' (\code{show_trend = TRUE}) is placed to the left and stacked with the same
#' \code{ggarrange} row structure so row heights always align.
#'
#' @param seurat_object A Seurat object.
#' @param group.by Character. Metadata column for point colour.
#' @param split.by Character. Metadata column for columns (samples/conditions).
#' @param row.by Character or \code{NULL}. Metadata column for rows.
#' @param show_trend Logical. Add a stacked composition trend panel to the
#'   left. Default \code{FALSE}.
#' @param show_labels Logical. Overlay numbered cluster labels. Default
#'   \code{FALSE}.
#' @param label_by Character or \code{NULL}. Column for label centroids.
#'   \code{NULL} → \code{group.by}.
#' @param colors Named colour vector. \code{NULL} auto-resolves from
#'   PrepObject or \code{Nour_pal}.
#' @param reduction Character. Reduction name (case-insensitive).
#' @param pt.size Numeric. Point size. Default \code{0.1}.
#' @param label.size Numeric. Label text size. Default \code{3}.
#' @param bar.width Numeric. Bar width in the trend plot. Default \code{0.6}.
#' @param legendnrow Integer. Legend rows. Default \code{2}.
#' @param legendtitle Character or \code{NA}. Legend title.
#' @param trend_width Numeric. Relative width of the trend panel vs. one UMAP
#'   column. Default \code{1}.
#' @param split.col.levels Character vector. Explicit ordering for
#'   \code{split.by} levels.
#' @param output_dir Character or \code{NULL}. Directory to save a PDF.
#' @param object_name Character. Output file-name prefix.
#' @param subset_name Character. Optional second prefix.
#' @param ... Ignored (forward compatibility).
#'
#' @return A \code{ggpubr} combined plot.
#' @export
PlotDimPlots <- function(seurat_object,
                          group.by         = NULL,
                          split.by         = NULL,
                          row.by           = NULL,
                          show_trend       = FALSE,
                          show_labels      = FALSE,
                          label_by         = NULL,
                          colors           = NULL,
                          reduction        = "umap",
                          pt.size          = 0.1,
                          label.size       = 3,
                          bar.width        = 0.6,
                          legendnrow       = 2,
                          legendtitle      = NA,
                          trend_width      = 1,
                          split.col.levels = NULL,
                          output_dir       = NULL,
                          object_name      = "",
                          subset_name      = "",
                          ...) {

  # ── Walk-up PrepObject defaults ────────────────────────────────────────────
  output_dir  <- output_dir %||% .nk_setting(seurat_object, "output_dir")
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""
  if (missing(group.by))
    group.by <- .nk_setting(seurat_object, "group.by") %||%
      stop("'group.by' is required. Provide it directly or store via PrepObject().")
  if (missing(split.by))
    split.by <- .nk_setting(seurat_object, "split.by")

  # ── Multiple group.by: recurse once per variable ──────────────────────────
  if (length(group.by) > 1L) {
    results <- lapply(stats::setNames(group.by, group.by), function(gb) {
      PlotDimPlots(
        seurat_object    = seurat_object,
        group.by         = gb,
        split.by         = split.by,
        row.by           = row.by,
        show_trend       = show_trend,
        show_labels      = show_labels,
        label_by         = label_by,
        colors           = colors,
        reduction        = reduction,
        pt.size          = pt.size,
        label.size       = label.size,
        bar.width        = bar.width,
        legendnrow       = legendnrow,
        legendtitle      = legendtitle,
        trend_width      = trend_width,
        split.col.levels = split.col.levels,
        output_dir       = output_dir,
        object_name      = object_name,
        subset_name      = subset_name
      )
    })
    return(invisible(results))
  }

  # ── Embeddings (case-insensitive) ──────────────────────────────────────────
  raw_emb  <- .resolve_emb(seurat_object, reduction)
  emb      <- as.data.frame(raw_emb[, 1:2, drop = FALSE])
  dim_labs <- colnames(raw_emb)[1:2]
  colnames(emb) <- c("Dim1", "Dim2")

  # ── Metadata ───────────────────────────────────────────────────────────────
  has_split <- !is.null(split.by)
  meta      <- seurat_object@meta.data
  req_cols  <- unique(c(group.by, if (has_split) split.by, row.by))
  bad       <- setdiff(req_cols, colnames(meta))
  if (length(bad))
    stop("Metadata column(s) not found: ", paste(bad, collapse = ", "))

  dat <- cbind(emb, meta[rownames(emb), req_cols, drop = FALSE])

  grp_lvls  <- if (is.factor(meta[[group.by]])) levels(meta[[group.by]])
               else sort(unique(as.character(meta[[group.by]])))
  dat$Group <- factor(as.character(dat[[group.by]]), levels = grp_lvls)

  if (has_split) {
    col_lvls     <- if (!is.null(split.col.levels)) split.col.levels
                    else if (is.factor(meta[[split.by]])) levels(meta[[split.by]])
                    else sort(unique(as.character(meta[[split.by]])))
    dat$ColSplit <- factor(as.character(dat[[split.by]]), levels = col_lvls)
  } else {
    col_lvls     <- group.by   # single panel; title = variable name
    dat$ColSplit <- factor(group.by, levels = col_lvls)
  }

  has_row <- !is.null(row.by)
  if (has_row) {
    row_lvls     <- if (is.factor(meta[[row.by]])) levels(meta[[row.by]])
                    else sort(unique(as.character(meta[[row.by]])))
    dat$RowSplit <- factor(as.character(dat[[row.by]]), levels = row_lvls)
  } else {
    row_lvls <- "All"   # single synthetic row for no-row case
  }

  n_rows <- length(row_lvls)

  # ── Colors ─────────────────────────────────────────────────────────────────
  colors <- colors %||% .nk_colors(seurat_object, group.by) %||%
    stats::setNames(
      Nour_pal(if (length(grp_lvls) <= 8) "all" else "spectrum")(length(grp_lvls)),
      grp_lvls
    )
  leg_title     <- if (is.na(legendtitle)) group.by else legendtitle
  cluster_nums  <- seq_along(grp_lvls)
  legend_labels <- stats::setNames(paste0(cluster_nums, ": ", grp_lvls), grp_lvls)
  num_lookup    <- stats::setNames(cluster_nums, grp_lvls)

  # ── Centroids (only when show_labels = TRUE) ───────────────────────────────
  centroid_df <- NULL
  if (show_labels) {
    c_grp <- if (has_row) dplyr::group_by(dat, RowSplit, ColSplit, Group)
             else         dplyr::group_by(dat, ColSplit, Group)
    centroid_df <- dplyr::summarise(c_grp,
      Dim1 = stats::median(Dim1), Dim2 = stats::median(Dim2), .groups = "drop")
    centroid_df <- dplyr::mutate(centroid_df,
      LabelNum = num_lookup[as.character(Group)])
    centroid_df <- centroid_df[!is.na(centroid_df$LabelNum), , drop = FALSE]
  }

  # ── Per-row column map (only columns that actually exist in each row) ───────
  if (has_row) {
    row_col_map <- lapply(stats::setNames(row_lvls, row_lvls), function(rl) {
      present <- as.character(dat$ColSplit[as.character(dat$RowSplit) == rl])
      col_lvls[col_lvls %in% unique(present)]
    })
  } else {
    row_col_map <- list(All = col_lvls)
  }
  max_cols <- max(vapply(row_col_map, length, integer(1)))

  # ── Internal: one UMAP panel ───────────────────────────────────────────────
  .make_umap <- function(panel_dat, title_txt, centroid_sub = NULL,
                          show_legend_here = FALSE) {
    p <- ggplot2::ggplot(panel_dat,
                         ggplot2::aes(x = Dim1, y = Dim2, color = Group)) +
      ggplot2::geom_point(size = pt.size, alpha = 0.7) +
      ggplot2::scale_color_manual(
        values = colors,
        labels = if (show_labels) legend_labels else grp_lvls,
        drop   = FALSE
      ) +
      ggplot2::labs(x = dim_labs[1], y = dim_labs[2]) +
      ggplot2::ggtitle(title_txt) +
      ggplot2::theme_classic() +
      ggplot2::scale_y_continuous(breaks = NULL) +
      ggplot2::scale_x_continuous(breaks = NULL) +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                 size = 13),
        panel.border     = ggplot2::element_rect(color = "darkgray", fill = NA),
        axis.line        = ggplot2::element_line(color = "darkgray")
      )

    if (!is.null(centroid_sub) && nrow(centroid_sub) > 0)
      p <- p + ggrepel::geom_text_repel(
        data         = centroid_sub,
        ggplot2::aes(x = Dim1, y = Dim2, label = LabelNum),
        color        = "black", size = label.size, fontface = "bold",
        max.overlaps = Inf, segment.size = 0.3,
        bg.color     = "white", bg.r = 0.1, inherit.aes = FALSE
      )

    if (!show_legend_here) {
      p <- p + Seurat::NoLegend()
    } else {
      p <- p +
        ggplot2::guides(colour = ggplot2::guide_legend(
          override.aes   = list(size = 5),
          nrow           = legendnrow,
          title.position = "top",
          title          = leg_title
        )) +
        ggplot2::theme(legend.position = "bottom",
                       legend.direction = "horizontal")
    }
    p
  }

  # ── Build UMAP row grids ───────────────────────────────────────────────────
  umap_row_grids <- lapply(row_lvls, function(rl) {
    if (has_row) {
      row_dat <- dat[as.character(dat$RowSplit) == rl, , drop = FALSE]
    } else {
      row_dat <- dat
    }
    cols_here <- row_col_map[[rl]]

    panels <- lapply(cols_here, function(cv) {
      pd  <- row_dat[as.character(row_dat$ColSplit) == cv, , drop = FALSE]
      ctr <- if (!is.null(centroid_df)) {
        if (has_row)
          centroid_df[centroid_df$ColSplit == cv & centroid_df$RowSplit == rl, ]
        else
          centroid_df[centroid_df$ColSplit == cv, ]
      } else NULL
      .make_umap(pd, cv, ctr, show_legend_here = FALSE)
    })

    # Pad to max_cols with invisible blank panels
    while (length(panels) < max_cols)
      panels[[length(panels) + 1L]] <- .blank_panel()

    row_g <- ggpubr::ggarrange(plotlist = panels, nrow = 1, ncol = max_cols)

    # Row label on right (only when row.by is used)
    if (has_row)
      row_g <- ggpubr::annotate_figure(
        row_g,
        right = ggpubr::text_grob(rl, rot = 270, face = "bold", size = 14)
      )
    row_g
  })

  # Stack row grids vertically
  umap_stack <- ggpubr::ggarrange(plotlist = umap_row_grids,
                                   nrow = n_rows, ncol = 1)

  # Shared legend (extracted from a full-data dummy panel)
  lgd_plot   <- .make_umap(dat, "", centroid_sub = NULL, show_legend_here = TRUE)
  shared_lgd <- ggpubr::get_legend(lgd_plot)

  # ── Build trend row plots (optional) ──────────────────────────────────────
  if (show_trend) {

    # Helper: build bar + area data for one row group with LOCAL col numbering
    .trend_data <- function(row_dat, rc_lvls) {
      num_map <- stats::setNames(seq_along(rc_lvls), as.character(rc_lvls))

      if (has_row) {
        bd <- dplyr::group_by(row_dat, ColSplit, Group, .drop = FALSE)
      } else {
        bd <- dplyr::group_by(row_dat, ColSplit, Group, .drop = FALSE)
      }
      bd <- dplyr::summarise(bd, Count = dplyr::n(), .groups = "drop")
      bd <- dplyr::group_by(bd, ColSplit)
      bd <- dplyr::mutate(bd,
        Fraction = Count / pmax(sum(Count), 1L),
        ColNum   = num_map[as.character(ColSplit)]
      )

      ad <- dplyr::bind_rows(
        dplyr::mutate(bd, X_Position = ColNum - bar.width / 2),
        dplyr::mutate(bd, X_Position = ColNum + bar.width / 2)
      )
      ad <- dplyr::arrange(ad, ColNum, X_Position)
      list(bar = bd, area = ad, num_map = num_map, rc_lvls = rc_lvls)
    }

    # Helper: build one trend ggplot from pre-built data
    .make_trend <- function(td, show_y_label = TRUE) {
      ggplot2::ggplot() +
        ggplot2::geom_area(
          data    = td$area,
          ggplot2::aes(x = X_Position, y = Fraction, fill = Group),
          alpha   = 0.5, position = "stack", color = "darkgray"
        ) +
        ggplot2::geom_bar(
          data    = td$bar,
          ggplot2::aes(x = ColNum, y = Fraction, fill = Group),
          stat = "identity", position = "stack",
          width = bar.width, alpha = 1, color = "black", linewidth = 0.2
        ) +
        ggplot2::scale_fill_manual(values = colors, labels = legend_labels) +
        ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
        ggplot2::scale_x_continuous(
          breaks = td$num_map, labels = td$rc_lvls, expand = c(0, 0)
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(y = if (show_y_label) "Percentage" else NULL, x = NULL) +
        ggplot2::theme(
          axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1,
                                                  size = 11, color = "black"),
          axis.text.y     = ggplot2::element_text(color = "black"),
          legend.position = "none",
          plot.margin     = ggplot2::margin(r = 10)
        )
    }

    trend_row_plots <- lapply(seq_along(row_lvls), function(i) {
      rl <- row_lvls[i]
      if (has_row) {
        row_dat <- dat[as.character(dat$RowSplit) == rl, , drop = FALSE]
      } else {
        row_dat <- dat
      }
      td <- .trend_data(row_dat, row_col_map[[rl]])
      .make_trend(td, show_y_label = (i == 1L))
    })

    # Stack trend rows with same nrow as UMAP — rows will align
    trend_stack <- ggpubr::ggarrange(plotlist = trend_row_plots,
                                      nrow = n_rows, ncol = 1,
                                      align = "v")

    # Combine trend + UMAP side by side (equal n_rows → rows align)
    inner <- ggpubr::ggarrange(
      trend_stack,
      umap_stack,
      nrow   = 1,
      ncol   = 2,
      widths = c(trend_width, max_cols)
    )
  } else {
    inner <- umap_stack
  }

  # Add shared legend below the combined figure
  result <- ggpubr::ggarrange(
    inner,
    shared_lgd,
    nrow    = 2,
    ncol    = 1,
    heights = c(1, 0.12)
  )

  # ── Save or return ──────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    parts <- c(
      if (nchar(object_name) > 0) object_name,
      if (nchar(subset_name) > 0) subset_name,
      group.by,
      "DimPlots"
    )
    fname <- gsub("[^A-Za-z0-9._-]", "_", paste(parts, collapse = "_"))
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))
    pdf_w <- max_cols * 3 + (if (show_trend) trend_width * 2 + 1 else 0) + 1
    pdf_h <- n_rows * 3 + 1.5   # +1.5 for legend
    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(result)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath)
    return(invisible(result))
  }

  result
}


# =============================================================================
# Backward-compatible wrappers
# =============================================================================

#' @describeIn PlotDimPlots Deprecated — use \code{PlotDimPlots(show_trend=FALSE)}.
#' @export
PlotGridUMAP <- function(seurat_object,
                          group.by,
                          split.by  = NULL,
                          row.by    = NULL,
                          colors    = NULL,
                          reduction = "umap",
                          pt.size   = 0.1,
                          col.split = NULL,
                          row.split = NULL) {
  if (is.null(split.by) && !is.null(col.split)) split.by <- col.split
  if (is.null(row.by)   && !is.null(row.split)) row.by   <- row.split
  if (is.null(split.by)) stop("'split.by' (column variable) is required.")
  PlotDimPlots(seurat_object,
               group.by   = group.by,
               split.by   = split.by,
               row.by     = row.by,
               colors     = colors,
               reduction  = reduction,
               pt.size    = pt.size,
               show_trend = FALSE)
}

#' @describeIn PlotDimPlots Deprecated — use \code{PlotDimPlots(show_trend=TRUE)}.
#' @export
PlotTrendAndUMAP <- function(seurat_object,
                              group.by,
                              split.by         = NULL,
                              row.by           = NULL,
                              colors           = NULL,
                              reduction        = "umap",
                              pt.size          = 0.1,
                              bar.width        = 0.6,
                              label.size       = 3,
                              legendnrow       = 2,
                              legendtitle      = NA,
                              trend_width      = 1,
                              split.col.levels = NULL,
                              split.cols       = NULL,
                              split.rows       = NULL) {
  if (is.null(split.by) && !is.null(split.cols)) split.by <- split.cols
  if (is.null(row.by)   && !is.null(split.rows)) row.by   <- split.rows
  if (is.null(split.by)) stop("'split.by' (column variable) is required.")
  PlotDimPlots(seurat_object,
               group.by         = group.by,
               split.by         = split.by,
               row.by           = row.by,
               show_trend       = TRUE,
               colors           = colors,
               reduction        = reduction,
               pt.size          = pt.size,
               bar.width        = bar.width,
               legendnrow       = legendnrow,
               legendtitle      = legendtitle,
               trend_width      = trend_width,
               split.col.levels = split.col.levels)
}
