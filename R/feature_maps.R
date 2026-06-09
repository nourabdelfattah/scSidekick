# =============================================================================
# scSidekick single-cell feature map generator
#
# GenerateFeatureMaps - batch FeaturePlot across genes, split by a metadata
#   variable, with shared per-gene color limits, a common legend, and
#   optional PDF output. Supports two layout modes:
#     "auto"     - auto-computes rows × cols from number of split levels
#     "metadata" - arranges panels in a row × column grid defined by two
#                  metadata variables
#
#   Optional add_boxplot = TRUE appends a per-row violin/box plot showing
#   expression across split.by groups with Wilcoxon significance brackets.
# =============================================================================

# ---------------------------------------------------------------------------
# Internal helper - build one expression violin/boxplot panel
# ---------------------------------------------------------------------------
.make_expr_boxplot <- function(df,
                                gene,
                                group_name,
                                plot_type,
                                alpha,
                                ref_group,
                                comparisons,
                                group_colors,
                                label_format,
                                hide_ns,
                                y_lim,
                                subtitle = NULL) {

  lvls   <- levels(df$group)
  n_lvls <- length(lvls)

  # ---- Colours - same resolution as PlotFeature ----
  fill_vals <- if (!is.null(group_colors) &&
                    all(lvls %in% names(group_colors))) {
    group_colors[lvls]
  } else {
    SelectColors(factor(lvls, levels = lvls), palette = "all")
  }

  # ---- Base plot ----
  p <- ggplot2::ggplot(df,
    ggplot2::aes(x = group, y = expr, fill = group)) +
    ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
    ggplot2::scale_y_continuous(
      limits = y_lim,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      x     = group_name,
      y     = paste0(gene, "\n(log-normalized)"),
      title = subtitle
    ) +
    theme_NourMin() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5, size = 11),
      plot.title   = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
      strip.text   = ggplot2::element_text(face = "bold"),
      plot.margin  = ggplot2::margin(t = 5, r = 5, b = 5, l = 10, unit = "mm")
    )

  if (plot_type %in% c("violin", "both")) {
    p <- p + ggplot2::geom_violin(trim = TRUE, scale = "width",
                                   alpha = alpha, color = NA)
  }
  if (plot_type %in% c("box", "both")) {
    bw <- if (plot_type == "both") 0.12 else 0.45
    if (plot_type == "both") {
      p <- p + ggplot2::geom_boxplot(
        width = bw, outlier.shape = NA,
        fill  = "white", alpha = alpha, linewidth = 0.35, color = "grey30"
      )
    } else {
      p <- p + ggplot2::geom_boxplot(
        width = bw, outlier.shape = NA,
        alpha = alpha, linewidth = 0.35, color = "grey30"
      )
    }
  }

  # ---- % expressed annotation at base of each group ----
  pct_df <- dplyr::summarise(
    dplyr::group_by(df, group),
    pct = round(mean(expr > 0) * 100, 1),
    .groups = "drop"
  )
  y_pct <- y_lim[1] + diff(y_lim) * 0.01
  for (lvl in lvls) {
    pct_val <- pct_df$pct[pct_df$group == lvl]
    if (length(pct_val)) {
      p <- p + ggplot2::annotate(
        "text", x = lvl, y = y_pct,
        label = paste0(pct_val, "%"),
        size = 2.2, color = "grey50", hjust = 0.5, vjust = 0
      )
    }
  }

  # ---- Determine comparison pairs ----
  all_pairs <- if (!is.null(comparisons)) {
    comparisons
  } else if (!is.null(ref_group) && ref_group %in% lvls) {
    lapply(setdiff(lvls, ref_group), function(g) c(ref_group, g))
  } else {
    utils::combn(lvls, 2, simplify = FALSE)
  }

  # ---- Compute Wilcoxon p-values ----
  stat_rows <- lapply(all_pairs, function(pair) {
    if (!all(pair %in% lvls)) return(NULL)
    g1 <- df$expr[df$group == pair[1]]
    g2 <- df$expr[df$group == pair[2]]
    if (length(g1) < 3 || length(g2) < 3) return(NULL)
    wt <- tryCatch(stats::wilcox.test(g1, g2), error = function(e) NULL)
    if (is.null(wt)) return(NULL)
    data.frame(group1 = pair[1], group2 = pair[2],
               p = wt$p.value, stringsAsFactors = FALSE)
  })
  stat_tbl <- do.call(rbind, Filter(Negate(is.null), stat_rows))

  if (!is.null(stat_tbl) && nrow(stat_tbl) > 0) {

    # BH-adjusted p-values
    stat_tbl$p.adj <- stats::p.adjust(stat_tbl$p, method = "BH")

    # Format labels
    if (label_format == "p.format") {
      stat_tbl$ann <- ifelse(
        stat_tbl$p.adj < 0.001,
        sprintf("p=%.1e", stat_tbl$p.adj),
        sprintf("p=%.3f", stat_tbl$p.adj)
      )
    } else {
      stat_tbl$ann <- dplyr::case_when(
        stat_tbl$p.adj < 0.0001 ~ "****",
        stat_tbl$p.adj < 0.001  ~ "***",
        stat_tbl$p.adj < 0.01   ~ "**",
        stat_tbl$p.adj < 0.05   ~ "*",
        TRUE                     ~ "ns"
      )
    }

    stat_tbl$is_ns <- stat_tbl$ann == "ns"

    # x positions → sort by span (short first, avoids crossings)
    x_pos <- stats::setNames(seq_along(lvls), lvls)
    stat_tbl$span <- abs(x_pos[stat_tbl$group2] - x_pos[stat_tbl$group1])
    stat_tbl <- stat_tbl[order(stat_tbl$is_ns, stat_tbl$span), ]

    # y positions - evenly spaced in the bracket zone (top 18% of axis)
    n_p   <- nrow(stat_tbl)
    bz_lo <- y_lim[1] + diff(y_lim) * 0.80
    bz_hi <- y_lim[1] + diff(y_lim) * 0.97
    stat_tbl$y_pos <- seq(bz_lo, bz_hi, length.out = max(n_p, 1))

    tsz  <- if (label_format == "p.format") 2.2 else 3.0
    sig_tbl <- stat_tbl[!stat_tbl$is_ns, ]
    ns_tbl  <- stat_tbl[ stat_tbl$is_ns, ]

    # Significant brackets - black
    if (nrow(sig_tbl) > 0) {
      p <- p + ggsignif::geom_signif(
        comparisons = lapply(seq_len(nrow(sig_tbl)),
          function(j) c(sig_tbl$group1[j], sig_tbl$group2[j])),
        annotations = sig_tbl$ann,
        y_position  = sig_tbl$y_pos,
        tip_length  = 0.01,
        vjust       = 0.4,
        size        = 0.4,
        textsize    = tsz,
        color       = "black",
        test        = NULL
      )
    }

    # ns brackets - grey (shown unless hide_ns = TRUE)
    if (!hide_ns && nrow(ns_tbl) > 0) {
      p <- p + ggsignif::geom_signif(
        comparisons = lapply(seq_len(nrow(ns_tbl)),
          function(j) c(ns_tbl$group1[j], ns_tbl$group2[j])),
        annotations = ns_tbl$ann,
        y_position  = ns_tbl$y_pos,
        tip_length  = 0.01,
        vjust       = 0.4,
        size        = 0.4,
        textsize    = tsz,
        color       = "grey60",
        test        = NULL
      )
    }
  }

  p
}


# ---------------------------------------------------------------------------
# Internal helper - add cluster/cell-type labels to a UMAP ggplot panel.
#
# centroid_df: data.frame with columns x, y, label, color
# repel:       logical - use ggrepel (TRUE) or geom_text (FALSE)
#
# Uses I(color) inside aes() to pass literal hex colors without interfering
# with the existing continuous color scale on the base FeaturePlot.
# ---------------------------------------------------------------------------
.add_umap_labels <- function(p, centroid_df, label.size, repel) {
  if (is.null(centroid_df) || nrow(centroid_df) == 0) return(p)

  if (isTRUE(repel) && requireNamespace("ggrepel", quietly = TRUE)) {
    p + ggrepel::geom_text_repel(
      data         = centroid_df,
      ggplot2::aes(x = x, y = y, label = label, colour = I(color)),
      size         = label.size,
      fontface     = "bold",
      box.padding  = 0.3,
      max.overlaps = Inf,
      bg.color     = "white",
      bg.r         = 0.1,
      segment.size = 0.25,
      min.segment.length = 0,
      show.legend  = FALSE
    )
  } else {
    p + ggplot2::geom_text(
      data        = centroid_df,
      ggplot2::aes(x = x, y = y, label = label, colour = I(color)),
      size        = label.size,
      fontface    = "bold",
      show.legend = FALSE
    )
  }
}


# Safe print wrapper - catches ggrepel render-time errors and returns a
# helpful message instead of crashing. Returns TRUE if print succeeded.
.safe_print_repel <- function(plot_obj) {
  tryCatch({
    print(plot_obj)
    invisible(TRUE)
  }, error = function(e) {
    message("Repel caused an error due to size constraints, ",
            "increase output size or remove repel.")
    invisible(FALSE)
  })
}


#' Batch single-cell feature plots with shared color limits
#'
#' Loops over `features`, generates a [Seurat::FeaturePlot()] split by
#' `split.by`, applies a shared per-gene color gradient (0 → max expression),
#' strips legends from individual panels, and attaches a single shared legend.
#'
#' When `add_boxplot = TRUE` a violin/box plot column is appended to the right
#' of the UMAP grid.  In `"metadata"` layout each row gets its own boxplot
#' (cells from that `row.by` level only); in `"auto"` layout a single boxplot
#' covers all cells.  Pairwise Wilcoxon tests are computed with BH correction
#' and displayed as brackets: significant comparisons in black, `ns` in grey.
#'
#' When `output_dir` is provided the figure is saved as a PDF and nothing is
#' returned. When `output_dir = NULL` the assembled plot is returned invisibly.
#'
#' @param seurat_object A Seurat object.
#' @param assay Character. assay to use for expression values. Default `"RNA"`.
#' @param reduction Character. Dimensional reduction for the UMAP. Default
#'   `"umap"`.
#' @param features Character vector of gene names to plot.
#' @param layout_method Character. `"auto"` (default) or `"metadata"`. In
#'   `"metadata"` mode `row.by` must be specified.
#' @param split.by Character. Metadata column to split UMAP panels by.
#' @param row.by Character or `NULL`. In `"metadata"` layout, the metadata
#'   column whose levels define the rows of the grid.
#' @param colors Character vector. Color gradient from low to high expression.
#' @param output_dir Character or `NULL`. Directory to write PDFs. If `NULL`,
#'   the plot is returned instead of saved.
#' @param object_name Character. Label appended to PDF filenames. Default `""`.
#' @param subset_name Character. Second label appended to PDF filenames.
#'   Default `""`.
#' @param order Logical. Plot higher-expression cells on top. Default `TRUE`.
#' @param label Logical. Overlay cell-type / cluster labels on UMAP panels.
#'   Default `FALSE`.
#' @param label.by Character. Metadata column to use for labels.  If `NULL`
#'   (default) and `label = TRUE`, labels are drawn using `Seurat::Idents()`.
#'   Accepts **any** metadata column - avoids Seurat's ident-order problem and
#'   the error that occurs when a column is literally named `"idents"`.
#' @param label_colors Named character vector mapping each label value to a
#'   hex colour (e.g. `c("Excitatory" = "#E41A1C", "Inhibitory" = "#377EB8")`).
#'   Labels not present in the vector default to `"black"`.  `NULL` = all black.
#' @param repel Logical. Use `ggrepel::geom_text_repel` to avoid overlapping
#'   labels (default `FALSE`).  If repel fails due to output-size constraints a
#'   message is printed and the panel is skipped rather than crashing.
#' @param label.size Numeric. Label font size. Default `2`.
#' @param pt.size Numeric. Point size. Default `0.001`.
#' @param add_boxplot Logical. Append an expression violin/box plot column.
#'   Default `FALSE`. Requires `split.by` to be set.
#' @param plot_type Character. Type of expression plot: `"violin"` (default),
#'   `"box"`, or `"both"` (violin with a thin white box overlay).
#' @param ref_group Character or `NULL`. If set, each group is compared only
#'   to this reference level. `NULL` (default) performs all pairwise tests.
#' @param comparisons List of character vectors or `NULL`. Custom pairs to
#'   test, e.g. `list(c("NCI","AD"), c("MCI","AD"))`. Overrides `ref_group`.
#' @param group_colors Named character vector or `NULL`. Fill colours for each
#'   `split.by` level. `NULL` auto-assigns from `Nour18`.
#' @param label_format Character. `"p.signif"` (default, shows `*`/`**`/`ns`)
#'   or `"p.format"` (shows `p = 0.023`).
#' @param hide_ns Logical. If `TRUE`, omit `ns` brackets entirely. Default
#'   `FALSE` (shows `ns` in grey).
#'
#' @return Invisibly returns the last assembled plot when `output_dir = NULL`;
#'   otherwise writes PDFs and returns `NULL` invisibly.
#' @export
GenerateFeatureMaps <- function(seurat_object,
                                 assay         = "RNA",
                                 reduction     = "umap",
                                 features,
                                 layout_method = "auto",
                                 split.by      = NULL,
                                 row.by       = NULL,
                                 colors        = c("#053061", "#2166AC",
                                                   "#D1E5F0", "#FDDBC7",
                                                   "#F4A582", "#D6604D",
                                                   "#B2182B", "#67001F"),
                                 output_dir       = NULL,
                                 object_name      = "",
                                 subset_name   = "",
                                 order         = TRUE,
                                 label         = FALSE,
                                 label.by      = NULL,
                                 label_colors  = NULL,
                                 repel         = FALSE,
                                 label.size    = 2,
                                 pt.size       = 0.001,
                                 # boxplot params
                                 add_boxplot   = FALSE,
                                 plot_type     = "violin",
                                 alpha         = 0.7,
                                 ref_group     = NULL,
                                 comparisons   = NULL,
                                 group_colors  = NULL,
                                 label_format  = "p.signif",
                                 hide_ns       = FALSE) {

  if (!layout_method %in% c("auto", "metadata"))
    stop("layout_method must be 'auto' or 'metadata'.")
  if (layout_method == "metadata" && is.null(row.by))
    stop("'row.by' must be specified when layout_method = 'metadata'.")
  if (add_boxplot && is.null(split.by))
    warning("add_boxplot = TRUE has no effect when split.by = NULL.")
  if (!plot_type %in% c("violin", "box", "both"))
    stop("plot_type must be 'violin', 'box', or 'both'.")

  # Walk up to PrepObject-stored defaults when not explicitly supplied
  output_dir  <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  # ── Multiple split.by: recurse once per variable ───────────────────────────
  if (!is.null(split.by) && length(split.by) > 1L) {
    results <- lapply(stats::setNames(split.by, split.by), function(sb) {
      GenerateFeatureMaps(
        seurat_object = seurat_object,
        assay         = assay,
        reduction     = reduction,
        features      = features,
        layout_method = layout_method,
        split.by      = sb,
        row.by       = row.by,
        colors        = colors,
        output_dir    = output_dir,
        object_name   = object_name,
        subset_name   = subset_name,
        order         = order,
        label         = label,
        label.by      = label.by,
        label_colors  = label_colors,
        repel         = repel,
        label.size    = label.size,
        pt.size       = pt.size,
        add_boxplot   = add_boxplot,
        plot_type     = plot_type,
        ref_group     = ref_group,
        comparisons   = comparisons,
        group_colors  = group_colors,
        label_format  = label_format,
        hide_ns       = hide_ns
      )
    })
    return(invisible(results))
  }

  # Resolve categorical colors from PrepObject
  # group_colors: fills for split.by levels in violin/box panels
  if (is.null(group_colors) && !is.null(split.by))
    group_colors <- .nk_colors(seurat_object, split.by)
  # label_colors: UMAP label text colors; resolve when label.by is set
  if (is.null(label_colors) && isTRUE(label) && !is.null(label.by))
    label_colors <- .nk_colors(seurat_object, label.by)

  Seurat::DefaultAssay(seurat_object) <- assay

  # ---- Pre-compute group levels (used for layout AND boxplots) ----
  if (!is.null(split.by)) {
    col_levels <- if (is.factor(seurat_object@meta.data[[split.by]]))
      levels(seurat_object@meta.data[[split.by]])
    else sort(unique(as.character(seurat_object@meta.data[[split.by]])))
  } else {
    col_levels <- NULL
  }

  # ---- Metadata layout: build temporary combined split variable ----
  if (layout_method == "metadata") {
    row_levels <- if (is.factor(seurat_object@meta.data[[row.by]]))
      levels(seurat_object@meta.data[[row.by]])
    else sort(unique(as.character(seurat_object@meta.data[[row.by]])))

    # Use @meta.data directly - compatible with Seurat v3 and v5; avoids the
    # seurat_object[[col, drop=TRUE]] syntax which behaves differently in v5.
    temp_split       <- "Temp_Grid_Split"
    seurat_object@meta.data[[temp_split]] <- paste(
      as.character(seurat_object@meta.data[[row.by]]),
      as.character(seurat_object@meta.data[[split.by]]),
      sep = "_"
    )
    active_split <- temp_split
    on.exit({
      if (temp_split %in% colnames(seurat_object@meta.data))
        seurat_object@meta.data[[temp_split]] <- NULL
    }, add = TRUE)
  } else {
    row_levels   <- NULL
    active_split <- split.by
  }

  # ---- Estimate number of pairs (for y-axis headroom) ----
  .n_pairs <- function(lvls) {
    if (!is.null(comparisons)) length(comparisons)
    else if (!is.null(ref_group)) max(length(lvls) - 1L, 1L)
    else max(choose(length(lvls), 2L), 1L)
  }

  # ---- Compute label centroids (once, before gene loop) ----
  # Centroids are computed from the full embedding regardless of split.by so
  # labels appear at biologically correct positions in every panel.
  centroid_df <- NULL
  if (isTRUE(label)) {
    # Determine which metadata column drives labels
    lby <- if (!is.null(label.by) && label.by %in% colnames(seurat_object@meta.data)) {
      label.by
    } else {
      # Fall back to Seurat Idents - stored as a temporary column to avoid
      # the 'idents' column-name conflict
      seurat_object@meta.data$.nk_ident_lbl <- as.character(Seurat::Idents(seurat_object))
      ".nk_ident_lbl"
    }

    red_key <- tolower(reduction)
    if (red_key %in% names(seurat_object@reductions)) {
      coords <- as.data.frame(
        SeuratObject::Embeddings(seurat_object, reduction = red_key)[, 1:2]
      )
      colnames(coords) <- c("x", "y")
      coords$lbl <- as.character(seurat_object@meta.data[[lby]])
      coords <- coords[!is.na(coords$lbl), ]

      centroid_df <- do.call(rbind, lapply(
        split(coords, coords$lbl),
        function(g) data.frame(
          label = g$lbl[1],
          x     = stats::median(g$x),
          y     = stats::median(g$y),
          stringsAsFactors = FALSE
        )
      ))

      # Assign colors: named vector → match by label name; unmatched → black
      centroid_df$color <- "black"
      if (!is.null(label_colors) && length(label_colors) > 0) {
        mc <- label_colors[centroid_df$label]
        ok <- !is.na(mc)
        centroid_df$color[ok] <- mc[ok]
      }
    }
  }

  last_plot <- NULL

  for (i in seq_along(features)) {
    gene      <- features[i]
    box_panel <- NULL   # reset each iteration so stale value never bleeds through

    # Extract full expression vector (needed for both colour scale and boxplot)
    gene_expr <- .get_layer_data(seurat_object, assay = assay, layer = "data",
                                  features = gene)[gene, ]
    max_exp   <- max(gene_expr, na.rm = TRUE)
    color_lim <- round(max(max_exp, 0.01), digits = 2)

    # Shared y-axis limits for boxplots:
    # data range (0 → max_exp) + headroom zone for brackets (top ~22%)
    if (add_boxplot && !is.null(split.by)) {
      n_p     <- .n_pairs(col_levels)
      headroom <- max_exp * (0.05 + n_p * 0.09)
      y_lim_box <- c(0, max_exp + headroom)
    }

    # ---- FeaturePlot ----
    # Always label = FALSE - labels are added manually below via .add_umap_labels
    basic <- Seurat::FeaturePlot(
      seurat_object, features = gene, reduction = reduction,
      split.by = active_split,
      order = order, label = FALSE,
      pt.size = pt.size
    ) & ggplot2::scale_color_gradientn(
      colors = colors, limits = c(0, color_lim), na.value = "red"
    )

    legend <- ggpubr::get_legend(
      basic[[1]] + theme_NourMin() + Seurat::NoAxes() +
        ggplot2::guides(colour = ggplot2::guide_colourbar(
          title          = gene,
          title.position = "top",
          title.theme    = ggplot2::element_text(size = 10)
        ))
    )

    n_plots <- length(basic)

    # ====================================================================
    # UMAP grid assembly
    # ====================================================================
    if (layout_method == "metadata") {

      # max_row_len is always length(col_levels): every row gets one slot per
      # column-level, whether a panel exists or not.  Empty slots are filled
      # with theme_void() placeholders RIGHT WHERE THEY BELONG, so split.by
      # levels that are absent in a given row.by level (e.g. "Reference" only
      # having Male samples) stay in the correct column instead of shifting left.
      max_row_len <- length(col_levels)

      row_plot_lists <- list()

      for (r in row_levels) {
        row_plots   <- list()
        has_any_panel <- FALSE

        for (co in col_levels) {
          combo        <- paste(r, co, sep = "_")
          panel_found  <- FALSE

          for (pp in seq_len(n_plots)) {
            if (!is.null(basic[[pp]]$labels$title) &&
                basic[[pp]]$labels$title == combo) {
              panel_p <- basic[[pp]] +
                theme_NourMin() + Seurat::NoLegend() + Seurat::NoAxes() +
                ggplot2::ggtitle(co) +
                ggplot2::theme(plot.title = ggplot2::element_text(
                  hjust = 0.5, face = "bold", size = 11))
              panel_p <- .add_umap_labels(panel_p, centroid_df, label.size, repel)
              row_plots[[length(row_plots) + 1]] <- panel_p
              panel_found   <- TRUE
              has_any_panel <- TRUE
              break
            }
          }

          if (!panel_found) {
            # Insert empty placeholder - preserves column alignment when a
            # split.by level is absent for this row.by level
            row_plots[[length(row_plots) + 1]] <-
              ggplot2::ggplot() + ggplot2::theme_void()
          }
        }

        row_plot_lists[[r]] <- if (has_any_panel) row_plots else list()
      }

      list_of_rows <- list()
      valid_row_levels <- character(0)

      for (r in row_levels) {
        rp <- row_plot_lists[[r]]
        if (length(rp) == 0) next
        valid_row_levels <- c(valid_row_levels, r)
        # No end-padding needed - every row already has max_row_len slots
        single_row <- ggpubr::ggarrange(plotlist = rp,
                                         nrow = 1, ncol = max_row_len)
        list_of_rows[[length(list_of_rows) + 1]] <- ggpubr::annotate_figure(
          single_row,
          left = ggpubr::text_grob(r, rot = 90, face = "bold", size = 12)
        )
      }

      n_rows     <- length(list_of_rows)
      n_cols     <- max_row_len
      inner_grid <- ggpubr::ggarrange(plotlist = list_of_rows,
                                       nrow = n_rows, ncol = 1)

      # ---- Boxplot column: one panel per valid row.by level ----
      if (add_boxplot && !is.null(split.by)) {
        box_plots <- lapply(valid_row_levels, function(r) {
          cell_idx <- which(
            as.character(seurat_object@meta.data[[row.by]]) == r
          )
          if (length(cell_idx) < 3)
            return(ggplot2::ggplot() + ggplot2::theme_void())

          df_r <- data.frame(
            expr  = as.numeric(gene_expr[cell_idx]),
            group = factor(
              as.character(seurat_object@meta.data[[split.by]][cell_idx]),
              levels = col_levels
            ),
            stringsAsFactors = FALSE
          )
          df_r <- df_r[!is.na(df_r$group), ]

          .make_expr_boxplot(
            df           = df_r,
            gene         = gene,
            group_name   = split.by,
            plot_type    = plot_type,
            alpha        = alpha,
            ref_group    = ref_group,
            comparisons  = comparisons,
            group_colors = group_colors,
            label_format = label_format,
            hide_ns      = hide_ns,
            y_lim        = y_lim_box,
            subtitle     = r
          )
        })
        box_panel <- ggpubr::ggarrange(plotlist = box_plots,
                                        nrow = n_rows, ncol = 1,
                                        align = "v")
      }

    } else {
      # ---- Auto layout ----
      n_rows    <- if (n_plots < 5) 1L else if (n_plots < 10) 2L else
        ceiling(n_plots / 4L)
      n_cols    <- ceiling(n_plots / n_rows)
      sub_list  <- lapply(seq_len(n_plots), function(pp) {
        panel_p <- basic[[pp]] + theme_NourMin() + Seurat::NoLegend() +
          Seurat::NoAxes() +
          ggplot2::theme(plot.title = ggplot2::element_text(
            hjust = 0.5, face = "bold", size = 11))
        .add_umap_labels(panel_p, centroid_df, label.size, repel)
      })
      inner_grid <- ggpubr::ggarrange(plotlist = sub_list,
                                       nrow = n_rows, ncol = n_cols)

      # ---- Single boxplot: all cells grouped by split.by ----
      if (add_boxplot && !is.null(split.by)) {
        df_all <- data.frame(
          expr  = as.numeric(gene_expr),
          group = factor(
            as.character(seurat_object@meta.data[[split.by]]),
            levels = col_levels
          ),
          stringsAsFactors = FALSE
        )
        df_all    <- df_all[!is.na(df_all$group), ]
        box_panel <- .make_expr_boxplot(
          df           = df_all,
          gene         = gene,
          group_name   = split.by,
          plot_type    = plot_type,
          alpha        = alpha,
          ref_group    = ref_group,
          comparisons  = comparisons,
          group_colors = group_colors,
          label_format = label_format,
          hide_ns      = hide_ns,
          y_lim        = y_lim_box
        )
      }
    }

    # ====================================================================
    # Assemble final figure
    # ====================================================================
    spacer  <- ggpubr::ggparagraph(text = " ", size = 0)
    box_w   <- max(1.3, n_cols * 0.35)    # boxplot column width (proportional)

    if (add_boxplot && !is.null(split.by) && !is.null(box_panel)) {
      CustomPlot <- ggpubr::ggarrange(
        spacer, inner_grid, box_panel, legend,
        nrow = 1, ncol = 4,
        widths = c(0.03, n_cols, box_w, 0.5)
      )
    } else {
      CustomPlot <- ggpubr::ggarrange(
        spacer, inner_grid, legend,
        nrow = 1, ncol = 3,
        widths = c(0.03, n_cols, 0.5)
      )
    }

    # ---- Legend text for InspectPlot / PPTX sidecar ----
    split_desc <- if (!is.null(split.by)) paste0(", split by ", split.by) else ""
    row_desc   <- if (!is.null(row.by))  paste0(", rows by ", row.by)   else ""
    box_desc   <- if (add_boxplot && !is.null(split.by)) {
      cmp_desc <- if (!is.null(comparisons))
        paste0("custom pairs: ",
               paste(sapply(comparisons, paste, collapse = " vs "),
                     collapse = "; "))
      else if (!is.null(ref_group))
        paste0("each group vs ", ref_group)
      else "all pairwise"
      paste0(" ",
             switch(plot_type,
               violin = "Violin plots",
               box    = "Box plots",
               both   = "Violin-box plots",
               "Plots"),
             " on the right show expression per ",
             split.by, " group (", cmp_desc,
             "; Wilcoxon rank-sum, BH-adjusted; significant comparisons in ",
             "black, ns in grey).",
             if (plot_type %in% c("box", "both"))
               " Box plot elements: centre line = median; box limits = 25th-75th percentile (IQR); whiskers extend to the furthest observation within 1.5x IQR from the box; outliers beyond this range are not shown."
             else "")
    } else ""

    leg_txt <- sprintf(
      paste0("UMAP feature plot showing expression of %s in %s%s%s. ",
             "Colour scale: low (dark blue) to high (dark red), capped at ",
             "the per-gene maximum (%.2f). Each panel represents one %s level; ",
             "legend shows the shared colour bar.%s"),
      gene,
      if (nchar(object_name) > 0) object_name else "the dataset",
      split_desc, row_desc,
      color_lim,
      if (!is.null(split.by)) split.by else "group",
      box_desc
    )
    attr(CustomPlot, "legend_text") <- leg_txt

    # ---- Save or return ----
    if (!is.null(output_dir)) {
      box_extra <- if (add_boxplot && !is.null(split.by)) box_w else 0
      pdf_h     <- n_rows * 3
      pdf_w     <- (n_cols * 3) + 2 + box_extra
      fname     <- file.path(
        output_dir,
        paste0(gene,
               if (!is.null(split.by)) paste0("_", split.by) else "",
               " featuremap ",
               object_name,
               if (nchar(subset_name) > 0) paste0(" ", subset_name) else "",
               ".pdf")
      )
      grDevices::pdf(fname, width = pdf_w, height = pdf_h)
      if (isTRUE(repel) && !is.null(centroid_df)) {
        .safe_print_repel(CustomPlot)
      } else {
        print(CustomPlot)
      }
      grDevices::dev.off()
      .write_legend_sidecar(fname, leg_txt)
    } else {
      # output_dir = NULL: always print so the user sees the plot even if they
      # forgot to capture the return value in a variable.
      last_plot <- CustomPlot
      if (isTRUE(repel) && !is.null(centroid_df)) {
        .safe_print_repel(CustomPlot)
      } else {
        print(CustomPlot)
      }
    }

    message(sprintf("Plotted %s (%d of %d) - layout: %dx%d%s",
                    gene, i, length(features), n_rows, n_cols,
                    if (add_boxplot && !is.null(split.by)) " + boxplot" else ""))
  }

  invisible(last_plot)
}

# Short-name alias
#' @describeIn GenerateFeatureMaps Preferred short alias.
#' @export
PlotFeaturePlots <- GenerateFeatureMaps
