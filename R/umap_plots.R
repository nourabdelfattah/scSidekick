# =============================================================================
# scSidekick UMAP / Dim-plot visualization  (umap_plots.R)
#
# PlotDimPlots()       - unified UMAP grid with optional composition trend
# PlotGridUMAP()       - deprecated thin wrapper → PlotDimPlots
# PlotTrendAndUMAP()   - deprecated thin wrapper → PlotDimPlots(show_trend=TRUE)
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
# PlotDimPlots  - main function
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
#' @param group.by Character. Metadata column for point color.
#' @param split.by Character. Metadata column for columns (samples/conditions).
#' @param row.by Character or \code{NULL}. Metadata column for rows.
#' @param number_labels Logical. When \code{TRUE}, prefix each legend entry
#'   and centroid label with a zero-padded index (e.g., "01. CellType").
#'   Useful when there are many groups and alphanumeric ordering aids
#'   interpretation. Default \code{FALSE}.
#' @param show_trend Logical. Add a stacked composition trend panel to the
#'   left. Default \code{FALSE}.
#' @param show_labels Logical. Overlay numbered cluster labels. Default
#'   \code{FALSE}.
#' @param label.by Character or \code{NULL}. Column for label centroids.
#'   \code{NULL} → \code{group.by}.
#' @param colors Named color vector. \code{NULL} auto-resolves from
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
#' @param width Numeric or \code{NULL}. Saved PDF width in inches. \code{NULL}
#'   (default) auto-sizes to fit the panel grid and the legend.
#' @param height Numeric or \code{NULL}. Saved PDF height in inches. \code{NULL}
#'   (default) auto-sizes from the number of rows plus the legend strip.
#' @param output_dir Character or \code{NULL}. Directory to save a PDF.
#' @param object_name Character. Output file-name prefix.
#' @param subset_name Character. Optional second prefix.
#' @param file_name Character or \code{NULL}. Base name (no extension) for the
#'   saved PDF. \code{NULL} (default) auto-deduces from \code{object_name},
#'   \code{subset_name}, and \code{group.by}.
#' @param join_plots Logical. When \code{group.by} has more than one variable,
#'   combine the per-variable plots into a \emph{single} figure / PDF instead of
#'   one file each. Default \code{FALSE}.
#' @param join_nrow,join_ncol Integer or \code{NULL}. Grid layout for
#'   \code{join_plots}. \code{NULL} auto-arranges into a near-square grid.
#' @param ... Ignored (forward compatibility).
#'
#' @return A \code{patchwork} combined plot.  With multiple \code{group.by} and
#'   \code{join_plots = FALSE}, a named list of plots (one per variable).
#' @export
PlotDimPlots <- function(seurat_object,
                          group.by         = NULL,
                          split.by         = NULL,
                          row.by           = NULL,
                          number_labels    = FALSE,
                          show_trend       = FALSE,
                          show_labels      = FALSE,
                          label.by         = NULL,
                          colors           = NULL,
                          reduction        = "umap",
                          pt.size          = 0.1,
                          label.size       = 3,
                          bar.width        = 0.6,
                          legendnrow       = 2,
                          legendtitle      = NA,
                          trend_width      = 1,
                          split.col.levels = NULL,
                          width            = NULL,
                          height           = NULL,
                          join_plots       = FALSE,
                          join_nrow        = NULL,
                          join_ncol        = NULL,
                          output_dir       = NULL,
                          object_name      = "",
                          subset_name      = "",
                          file_name        = NULL,
                          ...) {

  # ── Walk-up PrepObject defaults ────────────────────────────────────────────
  # NA_character_ is a sentinel passed by the outer join call meaning
  # "explicitly no save — do not walk up". Convert it to NULL before use.
  if (identical(output_dir, NA_character_)) output_dir <- NULL else
    output_dir <- output_dir %||%
      if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""
  group.by <- group.by %||% .nk_setting(seurat_object, "group.by") %||%
    stop("'group.by' is required. Provide it directly or store via PrepObject().")
  # Use %||% (NULL-check) rather than missing() so an explicit split.by value
  # is always respected regardless of whether the call comes from the console,
  # a lapply, or the multi-group recursive path.
  split.by <- split.by %||% .nk_setting(seurat_object, "split.by")

  # ── Multiple group.by ──────────────────────────────────────────────────────
  if (length(group.by) > 1L) {

    # Build one plot per variable. In join mode they are NOT saved individually
    # (output_dir = NULL) and width/height apply to the COMBINED figure, not
    # each panel; otherwise each variable saves its own PDF as before.
    sub_dir <- if (isTRUE(join_plots)) NA_character_ else output_dir

    plots <- lapply(stats::setNames(group.by, group.by), function(gb) {
      PlotDimPlots(
        seurat_object    = seurat_object,
        group.by         = gb,
        split.by         = split.by,
        row.by           = row.by,
        number_labels    = number_labels,
        show_trend       = show_trend,
        show_labels      = show_labels,
        label.by         = label.by,
        colors           = colors,
        reduction        = reduction,
        pt.size          = pt.size,
        label.size       = label.size,
        bar.width        = bar.width,
        legendnrow       = legendnrow,
        legendtitle      = legendtitle,
        trend_width      = trend_width,
        split.col.levels = split.col.levels,
        width            = if (isTRUE(join_plots)) NULL else width,
        height           = if (isTRUE(join_plots)) NULL else height,
        output_dir       = sub_dir,
        object_name      = object_name,
        subset_name      = subset_name,
        file_name        = if (isTRUE(join_plots)) NULL else file_name
      )
    })

    if (!isTRUE(join_plots)) return(invisible(plots))

    # ── Join: tile every per-variable plot into one figure ──────────────────
    n <- length(plots)
    if (is.null(join_ncol) && is.null(join_nrow)) {
      ncol_j <- ceiling(sqrt(n)); nrow_j <- ceiling(n / ncol_j)
    } else if (is.null(join_ncol)) {
      nrow_j <- join_nrow; ncol_j <- ceiling(n / nrow_j)
    } else {
      ncol_j <- join_ncol; nrow_j <- ceiling(n / ncol_j)
    }

    combined <- patchwork::wrap_plots(plots, nrow = nrow_j, ncol = ncol_j)

    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      base <- if (!is.null(file_name) && nzchar(file_name)) file_name else
        paste(c(if (nchar(object_name) > 0) object_name,
                if (nchar(subset_name) > 0) subset_name,
                paste(group.by, collapse = "-"), "DimPlots"), collapse = "_")
      fname <- gsub("[^A-Za-z0-9._-]", "_", base)
      fpath <- file.path(output_dir, paste0(fname, ".pdf"))

      # Per-plot footprint (same split.by/row.by for all) -> tile it.
      one_dim <- attr(plots[[1]], "nk_pdf_dims") %||% c(8, 6)
      pdf_w <- width  %||% (ncol_j * one_dim[1])
      pdf_h <- height %||% (nrow_j * one_dim[2])

      grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
      print(combined)
      grDevices::dev.off()
      message("scSidekick: Saved to ", fpath,
              " (", round(pdf_w, 1), " x ", round(pdf_h, 1), " in)")
      .write_legend_sidecar(fpath, paste0(
        toupper(reduction), " plots of ",
        paste(group.by, collapse = ", "),
        if (!is.null(split.by)) paste0(", split by ", split.by) else "",
        ", arranged in a ", nrow_j, " x ", ncol_j, " grid",
        if (nchar(object_name) > 0) paste0(" for ", object_name) else "", "."
      ))
      return(invisible(combined))
    }
    return(combined)
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
  leg_title    <- if (is.na(legendtitle)) group.by else legendtitle
  num_pad      <- sprintf("%02d", seq_along(grp_lvls))
  # number_labels: legend shows "01. Full Name"; centroid shows "01"
  # default:       legend shows "1: Full Name";  centroid shows "1"
  legend_labels <- if (number_labels) {
    stats::setNames(paste0(num_pad, ". ", grp_lvls), grp_lvls)
  } else {
    stats::setNames(paste0(seq_along(grp_lvls), ": ", grp_lvls), grp_lvls)
  }
  num_lookup <- if (number_labels) {
    stats::setNames(num_pad, grp_lvls)        # "01", "02", … on centroids
  } else {
    stats::setNames(seq_along(grp_lvls), grp_lvls)
  }

  # ── Centroids (only when show_labels = TRUE) ───────────────────────────────
  centroid_df <- NULL
  if (show_labels) {
    c_grp <- if (has_row) dplyr::group_by(dat, RowSplit, ColSplit, Group)
             else         dplyr::group_by(dat, ColSplit, Group)
    centroid_df <- dplyr::summarize(c_grp,
      Dim1 = stats::median(Dim1), Dim2 = stats::median(Dim2), .groups = "drop")
    centroid_df <- dplyr::mutate(centroid_df,
      LabelNum = num_lookup[as.character(Group)])
    centroid_df <- centroid_df[!is.na(centroid_df$LabelNum), , drop = FALSE]
  }

  # ── Per-row column map ────────────────────────────────────────────────────
  # Always show ALL split-by columns in every row so that e.g. a Sex × AD
  # grid has both AD levels in both Sex rows even when one combination is
  # absent from the data (renders as an empty panel with the column title).
  if (has_row) {
    row_col_map <- stats::setNames(rep(list(col_lvls), length(row_lvls)), row_lvls)
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
        labels = if (show_labels || number_labels) legend_labels else grp_lvls,
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
        # axis.line is blank: panel.border already provides the frame, and
        # having both creates double-lines at edges that visually clip points.
        axis.line        = ggplot2::element_blank()
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
        ggplot2::guides(color = ggplot2::guide_legend(
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

  # ── Helpers: trend data + plot ────────────────────────────────────────────
  # (defined here so they are available whether or not show_trend is TRUE)

  # Build bar + area data for one row group with LOCAL column numbering
  .trend_data <- function(row_dat, rc_lvls) {
    num_map <- stats::setNames(seq_along(rc_lvls), as.character(rc_lvls))
    bd <- dplyr::group_by(row_dat, ColSplit, Group, .drop = FALSE)
    bd <- dplyr::summarize(bd, Count = dplyr::n(), .groups = "drop")
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

  # Build one trend ggplot from pre-built data
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
        # b = 30 gives rotated x-axis labels room; patchwork will add matching
        # top whitespace so the trend panel aligns with UMAP panel areas.
        plot.margin     = ggplot2::margin(r = 10, b = 30)
      )
  }

  # ── Build individual UMAP panels per row ──────────────────────────────────
  # Keeping them as individual ggplots (not pre-assembled with ggarrange) lets
  # patchwork align their panel areas with the trend panel in the same row.
  umap_panels_by_row <- stats::setNames(
    lapply(row_lvls, function(rl) {
      row_dat <- if (has_row) dat[as.character(dat$RowSplit) == rl, , drop = FALSE]
                 else dat
      cols_here <- row_col_map[[rl]]

      panels <- lapply(cols_here, function(cv) {
        pd  <- row_dat[as.character(row_dat$ColSplit) == cv, , drop = FALSE]
        ctr <- if (!is.null(centroid_df)) {
          if (has_row) centroid_df[centroid_df$ColSplit == cv &
                                   centroid_df$RowSplit == rl, ]
          else         centroid_df[centroid_df$ColSplit == cv, ]
        } else NULL
        .make_umap(pd, cv, ctr, show_legend_here = FALSE)
      })

      # Pad shorter rows with invisible blank panels
      while (length(panels) < max_cols)
        panels[[length(panels) + 1L]] <- .blank_panel()

      panels
    }),
    row_lvls
  )

  # ── Build trend row plots (optional) ──────────────────────────────────────
  trend_row_plots <- if (show_trend) {
    lapply(seq_along(row_lvls), function(i) {
      rl      <- row_lvls[i]
      row_dat <- if (has_row) dat[as.character(dat$RowSplit) == rl, , drop = FALSE]
                 else dat
      td <- .trend_data(row_dat, row_col_map[[rl]])
      .make_trend(td, show_y_label = (i == 1L))
    })
  } else NULL

  # ── Combine per row using patchwork (guarantees panel-area alignment) ─────
  # patchwork reads each ggplot's gtable and adds matching whitespace so that
  # data rectangles line up, even when plots have different title/axis heights.
  .make_row_combo <- function(i) {
    rl          <- row_lvls[i]
    umap_panels <- umap_panels_by_row[[rl]]

    # Flatten all content panels into one list (trend + UMAPs)
    content <- if (show_trend)
      c(list(trend_row_plots[[i]]), umap_panels)
    else
      umap_panels

    if (has_row) {
      # Add row label as a flat element at the same patchwork level — avoids
      # ggpubr::annotate_figure() which only annotates the rightmost panel of a
      # wrap_plots result and destroys the multi-column structure.
      lbl <- patchwork::wrap_elements(
        full = grid::textGrob(rl, rot = 270,
                              gp = grid::gpar(fontface = "bold", fontsize = 14))
      )
      panels <- c(content, list(lbl))
      widths <- if (show_trend)
        c(trend_width, rep(1, max_cols), 0.1)
      else
        c(rep(1, max_cols), 0.1)
      patchwork::wrap_plots(panels, nrow = 1L, widths = widths)
    } else {
      widths <- if (show_trend) c(trend_width, rep(1L, max_cols)) else NULL
      patchwork::wrap_plots(content, nrow = 1L, widths = widths)
    }
  }

  row_combos <- lapply(seq_along(row_lvls), .make_row_combo)

  # Stack rows (single row: return as-is; multi-row: wrap with equal heights)
  inner <- if (n_rows == 1L) row_combos[[1L]]
           else patchwork::wrap_plots(row_combos, ncol = 1L,
                                      heights = rep(1L, n_rows))

  # ── Shared legend ──────────────────────────────────────────────────────────
  lgd_plot   <- .make_umap(dat, "", centroid_sub = NULL, show_legend_here = TRUE)
  shared_lgd <- ggpubr::get_legend(lgd_plot)

  # Estimate the space the legend needs from the number of levels and the
  # longest label, so it is never squeezed into a fixed fraction (the old
  # ggarrange heights = c(1, 0.12) cropped large legends).
  leg_labels    <- if (show_labels || number_labels) legend_labels else grp_lvls
  max_lbl_chars <- max(nchar(as.character(leg_labels)), 1L)
  per_row       <- ceiling(length(grp_lvls) / max(legendnrow, 1L))
  entry_w_in    <- 0.45 + max_lbl_chars * 0.075          # key + label width
  legend_w_in   <- per_row * entry_w_in + 1.0            # + legend title
  legend_h_in   <- 0.45 + legendnrow * 0.30              # title + rows

  # ── Final assembly via patchwork ──────────────────────────────────────────
  # patchwork stacks the grid over the legend, giving the legend a FIXED inch
  # height (a real "null" + fixed-unit layout) instead of a squeezable fraction.
  result <- patchwork::wrap_plots(
    inner,
    patchwork::wrap_elements(full = shared_lgd),
    ncol    = 1,
    heights = grid::unit.c(grid::unit(1, "null"), grid::unit(legend_h_in, "in"))
  )

  # ── Auto PDF size (computed always so join mode can tile by it) ────────────
  grid_w_in <- max_cols * 3 + (if (show_trend) trend_width * 2 + 1 else 0) + 1
  grid_h_in <- n_rows * 3
  pdf_w <- width  %||% max(grid_w_in, legend_w_in)
  pdf_h <- height %||% (grid_h_in + legend_h_in + 0.5)
  attr(result, "nk_pdf_dims") <- c(pdf_w, pdf_h)

  # ── Save or return ──────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # File name: file_name (verbatim) > object_subset_group.by_DimPlots
    if (!is.null(file_name) && nzchar(file_name)) {
      base <- file_name
    } else {
      parts <- c(
        if (nchar(object_name) > 0) object_name,
        if (nchar(subset_name) > 0) subset_name,
        group.by, "DimPlots"
      )
      base <- paste(parts, collapse = "_")
    }
    fname <- gsub("[^A-Za-z0-9._-]", "_", base)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))

    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(result)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath,
            " (", round(pdf_w, 1), " x ", round(pdf_h, 1), " in)")
    .write_legend_sidecar(fpath, paste0(
      toupper(reduction), " plot colored by ", group.by,
      if (!is.null(split.by)) paste0(", split into one panel per ", split.by,
                                     " level") else "",
      if (has_row) paste0(", with rows by ", row.by) else "",
      if (show_trend) "; the left panel shows the stacked composition trend" else "",
      ". A single shared color legend is shown below the panels",
      if (show_labels || number_labels)
        "; clusters are marked with numbered labels keyed to the legend" else "",
      if (nchar(object_name) > 0) paste0(". Dataset: ", object_name) else "", "."
    ))
    return(invisible(result))
  }

  result
}


# =============================================================================
# Backward-compatible wrappers
# =============================================================================

#' @describeIn PlotDimPlots Deprecated - use \code{PlotDimPlots(show_trend=FALSE)}.
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

#' @describeIn PlotDimPlots Deprecated - use \code{PlotDimPlots(show_trend=TRUE)}.
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
                              legendmargin     = NULL,  # accepted for back-compat; not forwarded
                              plotwidths       = NULL,  # c(trend, umap): trend_width <- plotwidths[1]
                              split.col.levels = NULL,
                              split.cols       = NULL,
                              split.rows       = NULL) {

  # ── Map legacy parameter names ─────────────────────────────────────────────
  if (is.null(split.by) && !is.null(split.cols)) split.by <- split.cols

  # split.rows = "all" (or NULL) means no row grouping — never treat "all" as
  # a metadata column name (which was the previous bug causing a crash).
  if (is.null(row.by) && !is.null(split.rows) && split.rows != "all")
    row.by <- split.rows

  if (is.null(split.by)) stop("'split.by' (column variable) is required.")

  # plotwidths = c(trend_frac, umap_frac): use the first value as trend_width
  # when the caller hasn't already set trend_width explicitly.
  if (!is.null(plotwidths) && trend_width == 1L)
    trend_width <- plotwidths[1L]

  PlotDimPlots(seurat_object,
               group.by         = group.by,
               split.by         = split.by,
               row.by           = row.by,
               show_trend       = TRUE,
               colors           = colors,
               reduction        = reduction,
               pt.size          = pt.size,
               bar.width        = bar.width,
               label.size       = label.size,
               legendnrow       = legendnrow,
               legendtitle      = legendtitle,
               trend_width      = trend_width,
               split.col.levels = split.col.levels)
}
