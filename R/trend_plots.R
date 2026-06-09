# =============================================================================
# scSidekick trend plots
#
# PlotTrendLabeled - stacked area/bar composition trend faceted by a split
#   variable, with ggrepel labels on the right side of the last time point.
#   Pure ggplot2 reimplementation - no SCP dependency.
# =============================================================================

#' Labeled stacked composition trend plot
#'
#' Computes cell-type proportions across ordered levels of a grouping variable
#' (e.g., time points or conditions), displays them as a stacked area + bar
#' chart, and places [ggrepel::geom_text_repel()] labels at the rightmost
#' position so each cell type is identified without a legend.
#'
#' Panels are faceted by `split.by` and combined with
#' [patchwork::wrap_plots()]. One panel per `split.by` level is produced.
#'
#' @param seurat_object A Seurat object.
#' @param stat.by Character. Metadata column whose proportions are plotted
#'   (e.g., cell type assignments). This becomes the fill color.
#' @param group.by Character. Ordered x-axis variable (e.g., `"TimePoint"`).
#' @param split.by Character. Metadata column that defines independent panels.
#' @param colors Named character vector mapping `stat.by` levels to colors.
#' @param group.levels Character vector. Explicit ordering of `group.by`
#'   levels. `NULL` uses factor or alphabetical order.
#' @param label.side Character. Where to place labels: `"right"` (default,
#'   last x position) or `"both"` (first and last).
#' @param bar.width Numeric. Width of stacked bars. Default `0.6`.
#' @param label.size Numeric. ggrepel text size. Default `3.5`.
#' @param label.nudge Numeric. Horizontal nudge for labels (in x-axis units).
#'   Default `0.3`.
#' @param label.margin Numeric. Right plot margin (pt) to accommodate labels.
#'   Default `80`.
#' @param segment.size Numeric. Connector line width. Default `0.4`.
#' @param max.overlaps Integer. Passed to [ggrepel::geom_text_repel()].
#'   Default `100`.
#' @param show.legend Logical. Show a fill legend? Default `FALSE` (labels
#'   replace the legend).
#' @param border.size Numeric. Panel border line width. Default `0.5`.
#' @param base_size Numeric. Base font size. Default `12`.
#' @param output_dir Character or `NULL`. Directory to save the PDF. Walks up
#'   from `PrepObject()` settings when `NULL`.
#' @param file_name Character or `NULL`. Base filename (no extension). Auto-
#'   deduced from `object_name`, `subset_name`, and column names when `NULL`.
#' @param object_name Character. Prefix for the auto-deduced filename.
#' @param subset_name Character. Second prefix component.
#' @param nrow,ncol Integer or `NULL`. Number of rows/columns for the panel
#'   grid. `NULL` (default) auto-arranges into a single row.
#' @param width,height Numeric. PDF dimensions in inches. Auto-sized when `NULL`.
#'
#' @return A patchwork object (one panel per `split.by` level). Returns
#'   invisibly when saved to disk.
#' @export
PlotTrendLabeled <- function(seurat_object,
                              stat.by,
                              group.by,
                              split.by,
                              colors        = NULL,
                              group.levels  = NULL,
                              label.side    = "right",
                              bar.width     = 0.6,
                              label.size    = 3,
                              label.nudge   = 0.1,
                              label.margin  = 120,
                              segment.size  = 0.4,
                              max.overlaps  = 100,
                              hjust         = 1,
                              show.legend   = FALSE,
                              border.size   = 0.5,
                              base_size     = 12,
                              nrow          = NULL,
                              ncol          = NULL,
                              output_dir    = NULL,
                              file_name     = NULL,
                              object_name   = NULL,
                              subset_name   = NULL,
                              width         = NULL,
                              height        = NULL) {

  if (!label.side %in% c("right", "both"))
    stop("label.side must be 'right' or 'both'.")

  # ── Walk up PrepObject settings ───────────────────────────────────────────
  if (is.null(output_dir))  output_dir  <- .nk_setting(seurat_object, "output_dir")
  if (is.null(object_name)) object_name <- .nk_setting(seurat_object, "object_name") %||% ""
  if (is.null(subset_name)) subset_name <- .nk_setting(seurat_object, "subset_name") %||% ""

  # Resolve colors from PrepObject or auto-generate
  if (is.null(colors)) {
    colors <- .nk_colors(seurat_object, stat.by)
    if (is.null(colors)) {
      stat_lvls <- if (is.factor(seurat_object@meta.data[[stat.by]]))
        levels(seurat_object@meta.data[[stat.by]])
      else sort(unique(as.character(seurat_object@meta.data[[stat.by]])))
      n <- length(stat_lvls)
      colors <- stats::setNames(
        Nour_pal(if (n <= 8) "all" else "spectrum")(n), stat_lvls)
    }
  }

  # 1. Fetch data
  dat <- Seurat::FetchData(seurat_object, vars = c(stat.by, group.by, split.by))
  colnames(dat) <- c("StatBy", "GroupBy", "SplitBy")

  # 2. Factor levels
  if (!is.null(group.levels)) {
    dat$GroupBy <- factor(dat$GroupBy, levels = group.levels)
  } else if (!is.factor(dat$GroupBy)) {
    dat$GroupBy <- factor(dat$GroupBy)
  }
  if (!is.factor(dat$StatBy)) dat$StatBy <- factor(dat$StatBy)

  grp_levels  <- levels(dat$GroupBy)
  stat_levels <- levels(dat$StatBy)
  n_x         <- length(grp_levels)

  # 3. Compute proportions
  prop_data <- dat |>
    dplyr::group_by(SplitBy, GroupBy, StatBy) |>
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
    tidyr::complete(SplitBy, GroupBy, StatBy, fill = list(Count = 0)) |>
    dplyr::group_by(SplitBy, GroupBy) |>
    dplyr::mutate(
      Proportion = Count / sum(Count),
      GroupNum   = as.numeric(GroupBy)
    ) |>
    dplyr::ungroup()

  # 4. Area ribbon data (flat-top effect)
  area_data <- dplyr::bind_rows(
    dplyr::mutate(prop_data, X_pos = GroupNum - bar.width / 2),
    dplyr::mutate(prop_data, X_pos = GroupNum + bar.width / 2)
  ) |> dplyr::arrange(SplitBy, GroupNum, X_pos)

  # 5. Label data - cumulative centroid at the rightmost (and optionally
  #    leftmost) x position
  label_x_vals <- if (label.side == "both") c(grp_levels[1], grp_levels[n_x]) else
    grp_levels[n_x]

  label_data <- prop_data |>
    dplyr::filter(as.character(GroupBy) %in% label_x_vals) |>
    dplyr::arrange(SplitBy, GroupBy, dplyr::desc(StatBy)) |>
    dplyr::group_by(SplitBy, GroupBy) |>
    dplyr::mutate(
      CumProp  = cumsum(Proportion),
      y_center = CumProp - 0.5 * Proportion
    ) |>
    dplyr::ungroup()

  # 6. Build one plot per split level
  split_levels <- unique(dat$SplitBy)
  plot_list    <- lapply(split_levels, function(sp) {

    pd  <- dplyr::filter(prop_data, SplitBy == sp)
    ad  <- dplyr::filter(area_data,  SplitBy == sp)
    ld  <- dplyr::filter(label_data, SplitBy == sp)

    p <- ggplot2::ggplot() +
      ggplot2::geom_area(
        data     = ad,
        ggplot2::aes(x = X_pos, y = Proportion, fill = StatBy),
        alpha    = 0.45, position = "stack", color = "darkgray",
        linewidth = 0.2
      ) +
      ggplot2::geom_bar(
        data     = pd,
        ggplot2::aes(x = GroupNum, y = Proportion, fill = StatBy),
        stat     = "identity", position = "stack",
        width    = bar.width, alpha = 1,
        color    = "black", linewidth = 0.2
      ) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        expand = c(0, 0)
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_len(n_x),
        labels = grp_levels,
        expand = ggplot2::expansion(add = c(0.3, 0.3))
      ) +
      ggplot2::labs(title = sp, y = "Proportion", x = NULL) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1,
                                                color = "black"),
        axis.text.y     = ggplot2::element_text(color = "black"),
        legend.position = if (show.legend) "right" else "none",
        panel.border    = ggplot2::element_rect(fill = NA, color = "black",
                                               linewidth = border.size),
        plot.margin     = ggplot2::margin(t = 10, r = label.margin,
                                         b = 10, l = 10)
      ) +
      ggplot2::coord_cartesian(clip = "off")

    # ggrepel labels at rightmost (and optionally leftmost) position.
    # Key: anchor x = GroupNum (actual bar position), nudge_x shifts slightly
    # rightward, then xlim = c(n_x * 5, NA) forces the FINAL text position
    # completely outside the panel (matches the proven xlim = c(20, NA) for
    # n_x=4 from the original reference script). coord_cartesian(clip="off")
    # lets labels render in the right-margin physical space.
    p <- p + ggrepel::geom_text_repel(
      data              = ld,
      ggplot2::aes(x = GroupNum, y = y_center,
                   label = StatBy, color = StatBy),
      direction         = "y",
      nudge_x           = label.nudge,
      hjust             = hjust,
      xlim              = c(n_x * 5, NA),
      segment.size      = segment.size,
      segment.color     = "black",
      segment.curvature = 0,
      max.overlaps      = max.overlaps,
      size              = label.size,
      box.padding       = 1,
      show.legend       = FALSE,
      inherit.aes       = FALSE
    ) +
      ggplot2::scale_color_manual(values = colors)

    p
  })

  names(plot_list) <- split_levels

  # ── Compute layout grid ───────────────────────────────────────────────────
  n_panels <- length(split_levels)
  if (!is.null(nrow) && !is.null(ncol)) {
    grid_nrow <- nrow; grid_ncol <- ncol
  } else if (!is.null(nrow)) {
    grid_nrow <- nrow; grid_ncol <- ceiling(n_panels / nrow)
  } else if (!is.null(ncol)) {
    grid_ncol <- ncol; grid_nrow <- ceiling(n_panels / ncol)
  } else {
    grid_nrow <- 1L;   grid_ncol <- n_panels
  }

  result <- patchwork::wrap_plots(plot_list, nrow = grid_nrow, ncol = grid_ncol) +
    patchwork::plot_layout(axes = "collect")

  # ── Save or return ────────────────────────────────────────────────────────
  if (.nk_autosave(seurat_object) && !is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    base <- if (!is.null(file_name) && nzchar(file_name)) file_name else {
      parts <- c(
        if (nchar(object_name) > 0) object_name,
        if (nchar(subset_name) > 0) subset_name,
        stat.by, group.by, "TrendLabeled"
      )
      paste(parts, collapse = "_")
    }
    fname <- gsub("[^A-Za-z0-9._-]", "_", base)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))

    pdf_w <- width  %||% (grid_ncol * 5)
    pdf_h <- height %||% (grid_nrow * 5)

    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(result)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath,
            " (", round(pdf_w, 1), " x ", round(pdf_h, 1), " in)")
    .write_legend_sidecar(fpath, paste0(
      "Stacked composition trend of ", stat.by, " across ", group.by,
      ", split by ", split.by, "."
    ))
    return(invisible(result))
  }

  result
}
