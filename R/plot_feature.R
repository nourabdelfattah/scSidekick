# =============================================================================
# plot_feature.R
#
# PlotFeature() - violin / box / bar plots for gene expression and numeric
#                 metadata variables, with full faceting, statistics, and
#                 PrepObject integration.
#
# Replaces SCP::FeatureStatPlot for scSidekick workflows.
# The violin panel used by GenerateFeatureMaps (add_boxplot = TRUE) also
# delegates to the shared internal builder .build_feature_panel().
# =============================================================================


# =============================================================================
# Internal: build one feature panel (one ggplot for one feature)
# =============================================================================
.build_feature_panel <- function(
    df,              # data.frame with columns: Value, Group, (Split), (Row)
    feature_name,    # character - displayed as y-axis label
    is_gene,         # logical - TRUE = add "(log-normalized)" to y label
    group.by,        # character - x-axis column name (for axis label)
    split.by,        # character or NULL
    row.by,          # character or NULL
    plot_type,       # "violin" / "box" / "both" / "bar"
    add_points,
    point_size,
    alpha,
    add_stats,
    comparisons,
    ref_group,
    hide_ns,
    stat_alpha,      # significance threshold (default 0.05)
    label_format,
    fill_colors      # named character vector
) {

  y_lab <- if (is_gene) paste0(feature_name, "\n(log-normalized)") else feature_name

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = Group, y = Value, fill = Group))

  # ── Geom layers ──────────────────────────────────────────────────────────
  if (plot_type %in% c("violin", "both"))
    p <- p + ggplot2::geom_violin(trim   = TRUE, scale = "width",
                                   alpha  = alpha, color = NA)

  if (plot_type %in% c("box", "both")) {
    bw <- if (plot_type == "both") 0.12 else 0.45
    # "both": white fill overlays the colored violin underneath
    # "box":  inherit fill from aes(fill = Group) so groups are colored
    if (plot_type == "both") {
      p <- p + ggplot2::geom_boxplot(
        width = bw, outlier.shape = NA,
        fill  = "white", alpha = 0.7, linewidth = 0.35, color = "gray30"
      )
    } else {
      p <- p + ggplot2::geom_boxplot(
        width = bw, outlier.shape = NA,
        alpha = 0.7, linewidth = 0.35, color = "gray30"
      )
    }
  }

  if (plot_type == "bar") {
    p <- p +
      ggplot2::stat_summary(fun = mean, geom = "bar", alpha = alpha) +
      ggplot2::stat_summary(
        fun.data  = ggplot2::mean_se, geom = "errorbar",
        width     = 0.25, linewidth = 0.5, color = "gray30"
      )
  }

  if (add_points)
    p <- p + ggplot2::geom_jitter(
      width = 0.15, size = point_size,
      alpha = 0.3, color = "black", show.legend = FALSE
    )

  # ── Statistics brackets ───────────────────────────────────────────────────
  if (add_stats && requireNamespace("ggsignif", quietly = TRUE)) {
    grp_lvls   <- levels(droplevels(df$Group))
    all_pairs  <- if (!is.null(comparisons)) {
      comparisons
    } else if (!is.null(ref_group) && ref_group %in% grp_lvls) {
      lapply(setdiff(grp_lvls, ref_group), function(g) c(ref_group, g))
    } else {
      utils::combn(grp_lvls, 2L, simplify = FALSE)
    }

    if (length(all_pairs) > 0L) {
      sig_args <- list(
        comparisons      = all_pairs,
        test             = "wilcox.test",
        test.args        = list(exact = FALSE),
        step_increase    = 0.08,
        tip_length       = 0.01,
        size             = 0.4,
        vjust            = 0.3,
        textsize         = if (label_format == "stars") 3.5 else 2.5
      )
      if (label_format == "stars") {
        sig_args$map_signif_level <- c("***" = stat_alpha / 50,
                                       "**"  = stat_alpha / 5,
                                       "*"   = stat_alpha,
                                       "ns"  = 1)
      } else {
        sig_args$map_signif_level <- FALSE
      }
      if (hide_ns)
        sig_args$map_signif_level <- c("***" = stat_alpha / 50,
                                       "**"  = stat_alpha / 5,
                                       "*"   = stat_alpha)

      p <- p + do.call(ggsignif::geom_signif, sig_args) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0.05, 0.35))
        )
    }
  }

  # ── Scales and labels ─────────────────────────────────────────────────────
  p <- p +
    ggplot2::scale_fill_manual(values = fill_colors, guide = "none") +
    ggplot2::labs(x = group.by, y = y_lab) +
    theme_NourMin() +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.x  = ggplot2::element_text(size = 11),
      axis.title.y  = ggplot2::element_text(angle = 90, vjust = 0.5, size = 11),
      strip.text    = ggplot2::element_text(face = "bold"),
      panel.spacing = ggplot2::unit(0.3, "lines"),
      plot.margin   = ggplot2::margin(t = 5, r = 5, b = 5, l = 10, unit = "mm")
    )

  # ── Faceting ──────────────────────────────────────────────────────────────
  if (!is.null(split.by) && !is.null(row.by)) {
    p <- p + ggplot2::facet_grid(
      rows   = ggplot2::vars(Row),
      cols   = ggplot2::vars(Split),
      scales = "free_y"
    )
  } else if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(~ Split)
  } else if (!is.null(row.by)) {
    p <- p + ggplot2::facet_wrap(~ Row, ncol = 1L)
  }

  p
}


# =============================================================================
# PlotFeature - public function
# =============================================================================

#' Violin / box / bar plots for gene expression and numeric metadata
#'
#' Plots one panel per feature - either a gene expression value pulled from
#' the specified assay layer or a numeric metadata column - grouped by a
#' categorical metadata variable on the x-axis.  Multiple features are
#' assembled with \pkg{patchwork}.
#'
#' Feature type is auto-detected: gene names found in the assay are extracted
#' as log-normalized expression; numeric metadata columns are used directly.
#' Passing a **categorical** metadata column raises a helpful error pointing
#' to the appropriate scSidekick visualization function instead.
#'
#' @param data A Seurat object or a plain data frame.
#' @param features Character vector of gene names and/or numeric metadata
#'   column names.  Mixed vectors are supported.
#' @param group.by Character.  Categorical metadata column for the x-axis
#'   (groups / conditions).
#' @param split.by Character or \code{NULL}.  Creates column facets within
#'   each feature panel - one panel per level of \code{split.by}.
#' @param row.by Character or \code{NULL}.  Creates row facets within each
#'   feature panel.  Combined with \code{split.by} produces a two-way
#'   \code{facet_grid}.
#' @param plot_type One of \code{"violin"} (default), \code{"box"},
#'   \code{"both"} (violin with a thin white box overlay), or \code{"bar"}
#'   (mean ± SE bar chart).
#' @param add_points Logical.  Overlay jittered individual points.
#'   Default \code{FALSE}.
#' @param point_size Numeric.  Size of jittered points.  Default \code{0.3}.
#' @param alpha Numeric.  Fill transparency for violins / bars.  Default
#'   \code{0.7}.
#' @param add_stats Logical.  Add Wilcoxon significance brackets via
#'   \pkg{ggsignif}.  Default \code{FALSE}.
#' @param comparisons List of length-2 character vectors specifying pairs to
#'   compare.  \code{NULL} (default) tests all pairwise combinations, or all
#'   vs \code{ref_group} when that is set.
#' @param ref_group Character or \code{NULL}.  Reference group for
#'   comparisons (all others tested against it).
#' @param hide_ns Logical.  Suppress \code{"ns"} brackets.  Default
#'   \code{FALSE}.
#' @param stat_alpha Numeric.  Significance threshold used for bracket labeling
#'   and for determining star levels (\code{*} at \code{stat_alpha}, \code{**}
#'   at \code{stat_alpha / 5}, \code{***} at \code{stat_alpha / 50}).  Default
#'   \code{0.05}.
#' @param label_format One of \code{"stars"} (default: \code{*},
#'   \code{**}, \code{***}) or \code{"p.format"} (numeric p-value shown).
#' @param ncol Integer or \code{NULL}.  Columns in the patchwork assembly.
#'   \code{NULL} = \code{min(n_features, 3)}.
#' @param downsample Integer or \code{NULL}.  If supplied, downsample to at
#'   most this many cells total (balanced across \code{group.by} levels) before
#'   plotting.  \code{NULL} (default) uses all cells.
#' @param exclude Named list of values to \strong{exclude} before plotting.
#'   Each name is a column name and each value is a character vector of levels
#'   to drop.  Can target \code{group.by}, \code{split.by}, \code{row.by}, or
#'   any other metadata column.
#'   Example: \code{list(Cognitive.Status = "Reference", Sex = "Unknown")}.
#' @param colors Named character vector of fill colors for \code{group.by}
#'   levels.  \code{NULL} auto-resolves from \code{PrepObject} (if a Seurat
#'   object is passed) or \code{SelectColors()}.
#' @param assay Character.  Seurat assay to pull gene expression from.
#'   Default \code{"RNA"}.
#' @param layer Character.  Assay layer / slot.  Default \code{"data"}
#'   (log-normalized counts).
#' @param output_dir Character or \code{NULL}.  Save directory.  Walks up
#'   from \code{PrepObject} when a Seurat object is passed.
#' @param object_name Character.  Filename prefix.  Walks up from
#'   \code{PrepObject} when a Seurat object is passed.
#' @param pdf_width Numeric or \code{NULL}.  Override the auto-calculated PDF
#'   width in inches.  \code{NULL} (default) sizes automatically from the
#'   number of groups, splits, and features.
#' @param pdf_height Numeric or \code{NULL}.  Override the auto-calculated PDF
#'   height in inches.  \code{NULL} (default) sizes automatically from the
#'   number of row facets and features.
#'
#' @return When one feature is requested, a single \code{ggplot2} object.
#'   When multiple features are requested, a \code{patchwork} combined plot.
#'   The plot is also saved as a PDF when \code{output_dir} is available.
#' @export
#'
#' @examples
#' \dontrun{
#' # Gene expression violin plot
#' PlotFeature(SeuratObj, features = c("CD3E", "CD79A", "LYZ"),
#'   group.by = "CellType")
#'
#' # Numeric metadata with split + stats, excluding a reference group
#' PlotFeature(SeuratObj, features = "Age.at.Death",
#'   group.by  = "Dementia.AD",
#'   split.by  = "Sex",
#'   exclude   = list(Cognitive.Status = "Reference"),
#'   add_stats = TRUE,
#'   plot_type = "box")
#'
#' # Two-way facet: split.by columns × row.by rows
#' PlotFeature(SeuratObj, features = c("CD3E", "LYZ"),
#'   group.by = "Group",
#'   split.by = "Timepoint",
#'   row.by   = "Sex")
#' }
PlotFeature <- function(data,
                         features,
                         group.by,
                         split.by     = NULL,
                         row.by       = NULL,
                         exclude      = NULL,
                         plot_type    = c("violin", "box", "both", "bar"),
                         add_points   = FALSE,
                         point_size   = 0.3,
                         alpha        = 0.7,
                         add_stats    = FALSE,
                         comparisons  = NULL,
                         ref_group    = NULL,
                         hide_ns      = FALSE,
                         stat_alpha   = 0.05,
                         label_format = c("stars", "p.format"),
                         ncol         = NULL,
                         downsample   = NULL,
                         colors       = NULL,
                         assay        = "RNA",
                         layer        = "data",
                         output_dir   = NULL,
                         object_name  = "",
                         pdf_width    = NULL,
                         pdf_height   = NULL) {

  # Accept c("violin") or "violin" - take first element before match.arg
  if (length(plot_type) > 1L) plot_type <- plot_type[1L]
  plot_type    <- match.arg(plot_type, c("violin", "box", "both", "bar"))
  label_format <- match.arg(label_format)

  # ── Walk-up PrepObject defaults ───────────────────────────────────────────
  if (inherits(data, "Seurat")) {
    output_dir  <- output_dir %||%
      if (.nk_autosave(data)) .nk_setting(data, "output_dir") else NULL
    object_name <- if (nchar(object_name) > 0) object_name else
      .nk_setting(data, "object_name") %||% ""
  }

  # ── Extract metadata ──────────────────────────────────────────────────────
  meta <- if (inherits(data, "Seurat")) data@meta.data else as.data.frame(data)

  # ── Validate group.by ─────────────────────────────────────────────────────
  if (!group.by %in% colnames(meta))
    stop("'group.by' column '", group.by, "' not found in metadata.")
  if (is.numeric(meta[[group.by]]))
    stop("'group.by' must be a categorical variable. '", group.by,
         "' appears to be numeric. Did you mean to use it as a feature instead?")

  # ── Validate split.by / row.by ────────────────────────────────────────────
  for (v in c(split.by, row.by)) {
    if (!is.null(v) && !v %in% colnames(meta))
      stop("Column '", v, "' not found in metadata.")
  }

  # ── Apply exclusions ──────────────────────────────────────────────────────
  if (!is.null(exclude)) {
    for (col in names(exclude)) {
      if (!col %in% colnames(meta)) {
        warning("Exclusion column '", col,
                "' not found in metadata - skipping.")
        next
      }
      meta <- meta[!as.character(meta[[col]]) %in%
                     as.character(exclude[[col]]), , drop = FALSE]
    }
    message("scSidekick PlotFeature: ",
            format(nrow(meta), big.mark = ","),
            " cells remain after exclusions.")
  }

  # ── Get available gene names (assay-agnostic) ─────────────────────────────
  gene_names <- if (inherits(data, "Seurat")) {
    tryCatch(
      rownames(SeuratObject::GetAssay(data, assay = assay)),
      error = function(e)
        tryCatch(rownames(data[[assay]]), error = function(e2) character(0))
    )
  } else character(0)

  # ── Classify each feature ─────────────────────────────────────────────────
  valid_features  <- character(0)
  feature_sources <- list()   # "gene" or "meta"

  for (feat in features) {
    if (feat %in% colnames(meta)) {
      if (!is.numeric(meta[[feat]])) {
        stop(
          "'", feat, "' is a categorical metadata column.\n",
          "PlotFeature() is for numeric variables and gene expression.\n",
          "For categorical variables use:\n",
          "  • PlotDimPlots()   - UMAP colored by category\n",
          "  • GroupHeatmap()  - average expression heatmap\n",
          "  • SplitDotPlot()  - dot plot by marker gene set\n",
          "  • PlotMetaSummary() - patient-level distribution bars"
        )
      }
      valid_features <- c(valid_features, feat)
      feature_sources[[feat]] <- "meta"

    } else if (feat %in% gene_names) {
      valid_features <- c(valid_features, feat)
      feature_sources[[feat]] <- "gene"

    } else {
      warning("Feature '", feat, "' not found in metadata or assay '",
              assay, "' - skipping.")
    }
  }

  if (length(valid_features) == 0L)
    stop("No valid features found. Check spelling and assay name (assay = \"",
         assay, "\").")

  # ── Resolve group.by levels and colors ────────────────────────────────────
  grp_vals <- as.character(meta[[group.by]])
  grp_lvls <- if (is.factor(meta[[group.by]])) levels(meta[[group.by]])
              else sort(unique(grp_vals))
  # Drop levels absent from data - critical after exclusions so ggsignif
  # doesn't try to compare against an empty (excluded) group.
  grp_lvls <- intersect(grp_lvls, unique(grp_vals))

  fill_colors <- colors
  if (is.null(fill_colors) && inherits(data, "Seurat"))
    fill_colors <- tryCatch(.nk_colors(data, group.by), error = function(e) NULL)
  if (is.null(fill_colors))
    fill_colors <- SelectColors(
      factor(grp_vals, levels = grp_lvls), palette = "all"
    )

  # ── Resolve cells for gene features (once, before the loop) ─────────────
  # Keep the lazy matrix reference so each gene is extracted with a single
  # BPCells streaming scan (same pattern as PercentageFeatureSet - fast).
  has_gene_feat <- any(vapply(feature_sources, `==`, logical(1), "gene"))

  mat_lazy      <- NULL   # lazy BPCells / Seurat matrix reference
  gene_cell_idx <- NULL   # integer column indices into mat_lazy

  gene_cells <- if (has_gene_feat) {
    # Grab the lazy matrix - colnames() is metadata-only for BPCells
    mat_lazy <- tryCatch(
      SeuratObject::LayerData(data, assay = assay, layer = layer),
      error = function(e1) tryCatch(
        Seurat::GetAssayData(data, assay = assay, layer = layer),
        error = function(e2)
          # Seurat v3/v4 fallback: `slot=` is the correct argument for old objects
          # but triggers a lifecycle warning in SeuratObject 5+.  Suppress it here
          # since the warning is a false positive for genuinely old objects.
          suppressWarnings(
            Seurat::GetAssayData(data, assay = assay, slot = layer)
          )
      )
    )
    assay_cell_names <- colnames(mat_lazy)
    common <- intersect(assay_cell_names, rownames(meta))

    if (length(common) == 0L)
      stop("No cells overlap between assay '", assay,
           "' and the Seurat metadata. Check the assay name.")
    if (length(common) < nrow(meta))
      message("scSidekick: Assay '", assay, "' has ",
              format(length(common), big.mark = ","),
              " cells (sketch / subset of ",
              format(nrow(meta), big.mark = ","), " total). ",
              "Plotting those cells only.")

    # Optional downsample (off by default - preserves rare populations)
    if (!is.null(downsample) && length(common) > downsample) {
      grp_for_ds <- as.character(meta[common, group.by])
      lvls_ds    <- unique(grp_for_ds)
      per_grp    <- max(1L, ceiling(downsample / length(lvls_ds)))
      common     <- unlist(lapply(lvls_ds, function(g) {
        cells_g <- common[grp_for_ds == g]
        if (length(cells_g) <= per_grp) cells_g else sample(cells_g, per_grp)
      }), use.names = FALSE)
      message("scSidekick: Downsampled to ",
              format(length(common), big.mark = ","), " cells.")
    }

    # Pre-compute integer indices once - faster than name lookup per gene
    gene_cell_idx <- match(common, assay_cell_names)
    common
  } else rownames(meta)

  gene_meta <- meta[gene_cells, , drop = FALSE]

  # ── Cached split/row factor levels ────────────────────────────────────────
  split_lvls <- if (!is.null(split.by))
    if (is.factor(meta[[split.by]])) levels(meta[[split.by]])
    else sort(unique(as.character(meta[[split.by]])))
  else NULL

  row_lvls <- if (!is.null(row.by))
    if (is.factor(meta[[row.by]])) levels(meta[[row.by]])
    else sort(unique(as.character(meta[[row.by]])))
  else NULL

  # ── Build one panel per feature ───────────────────────────────────────────
  plot_list <- lapply(valid_features, function(feat) {

    # ── Extract values using BPCells-efficient streaming ─────────────────────
    if (feature_sources[[feat]] == "meta") {
      vals     <- as.numeric(meta[[feat]])
      use_meta <- meta
    } else {
      # as.matrix() first: converts BPCells S4 to a plain R object.
      # as.numeric() then: flattens to a vector regardless of whether
      # as.matrix() returned a matrix, named vector, or sparse matrix.
      vals     <- as.numeric(as.matrix(
        mat_lazy[feat, gene_cell_idx, drop = FALSE]
      ))
      use_meta <- gene_meta
    }

    df <- data.frame(
      Value = vals,
      Group = factor(as.character(use_meta[[group.by]]), levels = grp_lvls),
      stringsAsFactors = FALSE
    )
    if (!is.null(split.by))
      df$Split <- factor(as.character(use_meta[[split.by]]), levels = split_lvls)
    if (!is.null(row.by))
      df$Row   <- factor(as.character(use_meta[[row.by]]),   levels = row_lvls)
    df <- df[!is.na(df$Value), , drop = FALSE]

    .build_feature_panel(
      df           = df,
      feature_name = feat,
      is_gene      = feature_sources[[feat]] == "gene",
      group.by     = group.by,
      split.by     = split.by,
      row.by       = row.by,
      plot_type    = plot_type,
      add_points   = add_points,
      point_size   = point_size,
      alpha        = alpha,
      add_stats    = add_stats,
      comparisons  = comparisons,
      ref_group    = ref_group,
      hide_ns      = hide_ns,
      stat_alpha   = stat_alpha,
      label_format = label_format,
      fill_colors  = fill_colors
    )
  })
  names(plot_list) <- valid_features

  # ── Assemble with patchwork ───────────────────────────────────────────────
  n_feat  <- length(plot_list)
  ncol    <- ncol %||% min(n_feat, 3L)
  combined <- if (n_feat == 1L) plot_list[[1L]]
              else patchwork::wrap_plots(plot_list, ncol = ncol)

  # ── Auto-save PDF + .legend sidecar ──────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    n_groups <- length(grp_lvls)
    n_splits <- if (!is.null(split.by)) length(split_lvls) else 1L
    n_rows_f <- if (!is.null(row.by))   length(row_lvls)   else 1L
    nrow_pw  <- ceiling(n_feat / ncol)

    # Width: each feature panel has n_splits facet columns,
    #        each column shows n_groups bars.
    #        0.45 in/bar, minimum 1.8 in per facet column.
    facet_col_w <- max(n_groups * 0.45, 1.8)
    panel_w     <- facet_col_w * n_splits + 0.8   # +0.8 for y-axis + strip

    # Height: each feature panel has n_rows_f facet rows,
    #         each row needs room for the violin + rotated x-labels.
    #         2.8 in/row, minimum 2.5 in per facet row.
    facet_row_h <- max(2.8, 2.5)
    panel_h     <- facet_row_h * n_rows_f + 1.0   # +1.0 for x-labels + strip

    # Total patchwork size - cap at 50 in to keep PDFs openable
    # User-supplied pdf_width / pdf_height override the auto-calculation
    pdf_w <- pdf_width  %||% min(panel_w * ncol + 0.5, 50)
    pdf_h <- pdf_height %||% min(panel_h * nrow_pw + 0.5, 50)

    excl_tags <- if (!is.null(exclude) && length(exclude) > 0L)
      unlist(lapply(names(exclude), function(col)
        paste0("no_", paste(exclude[[col]], collapse = "_"))))
    else NULL

    parts <- c(
      if (nchar(object_name) > 0) object_name,
      paste(valid_features, collapse = "_"),
      group.by,
      split.by,
      row.by,
      excl_tags,
      plot_type,
      "PlotFeature"
    )
    fname <- gsub("[^A-Za-z0-9._-]", "_", paste(parts, collapse = "_"))
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))

    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(combined)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath,
            " (", round(pdf_w, 1), " × ", round(pdf_h, 1), " in)")

    .write_legend_sidecar(fpath, paste0(
      switch(plot_type,
        violin = "Violin plot",  box  = "Box plot",
        both   = "Violin-box plot", bar = "Bar chart (mean ± SE)"),
      " of ", paste(valid_features, collapse = ", "),
      " grouped by ", group.by,
      if (!is.null(split.by)) paste0(", split by ", split.by) else "",
      if (!is.null(row.by))   paste0(", rows by ", row.by)   else "",
      ". Gene expression values are log-normalized counts.",
      if (plot_type %in% c("box", "both"))
        " Box plot elements: center line = median; box limits = 25th-75th percentile (IQR); whiskers extend to the furthest observation within 1.5x IQR from the box; outliers beyond this range are not shown."
      else "",
      if (add_stats)
        " Statistical comparisons: Wilcoxon rank-sum test (brackets show significance)."
      else "",
      if (!is.null(exclude) && length(exclude) > 0)
        paste0(" The following groups were excluded before plotting: ",
               paste(mapply(function(col, vals)
                 paste0(col, " = ", paste(vals, collapse = ", ")),
                 names(exclude), exclude), collapse = "; "), ".")
      else "",
      if (!is.null(object_name) && nchar(object_name) > 0)
        paste0(" Dataset: ", object_name, ".") else ""
    ))
  }

  combined
}
