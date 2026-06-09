# =============================================================================
# scSidekick - PlotMultiFeature  (viz_multifeature.R)
#
# Modernised from the classic Seurat v2 MultiFeaturePlot helper.
# Supports Seurat v3/v5 and BPCells via .get_layer_data().
# =============================================================================

#' Plot multiple features on a single embedding
#'
#' Overlays several genes onto one UMAP / tSNE / PCA plot, coloring cells by
#' which gene they express above a threshold.  Three modes are available:
#'
#' \describe{
#'   \item{\strong{"gene"} (default)}{Each cell is colored by the last gene
#'     in \code{features} that it expresses above threshold (later genes
#'     overwrite earlier ones for multi-expressing cells).  Cells below all
#'     thresholds are drawn in \code{null_color}.}
#'   \item{\strong{"intersection"}}{Cells that express at least
#'     \code{pct_cutoff} × n_genes genes are highlighted as co-expressing.
#'     All others are \code{null_color}.}
#'   \item{\strong{"scale"}}{Continuous color gradient showing expression
#'     level of a single gene (\code{features[1]} only).}
#' }
#'
#' @param seurat_object A Seurat object (v3 or v5).
#' @param features Character vector of gene names to visualise.
#' @param colors Named or unnamed character vector of colors (one per gene).
#'   Auto-filled from \code{Nour_pal("all")} for missing entries. If
#'   \code{NULL}, the full \code{Nour_pal("all")} palette is used.
#' @param null_color Color for cells that express no gene above threshold.
#'   Default \code{"grey70"}.
#' @param threshold Numeric scalar, or named numeric vector (names = gene
#'   names), giving the minimum expression value for a cell to be considered
#'   positive.  Default \code{0.0125}.
#' @param plot_mode One of \code{"gene"} (default), \code{"intersection"},
#'   or \code{"scale"}.
#' @param pct_cutoff For \code{plot_mode = "intersection"}: the fraction of
#'   the gene list a cell must express to count as co-expressing.  Default
#'   \code{0.3} (i.e. >= 30% of genes).
#' @param intersection_label Label shown in the legend for co-expressing cells.
#'   Default \code{"Co-expressing"}.
#' @param assay Assay to pull expression from.  Default \code{"RNA"}.
#' @param layer Layer / slot to pull (Seurat v5 layer name or v3 slot name).
#'   Default \code{"data"}.
#' @param reduction Dimensional reduction to use for coordinates.  Default
#'   \code{"umap"}.
#' @param pt.size Point size.  Default \code{0.5}.
#' @param alpha Point transparency \code{[0, 1]}.  Default \code{0.9}.
#' @param shape ggplot2 point shape integer.  Default \code{16} (filled
#'   circle).
#' @param title Plot title.  Default \code{NULL}.
#' @param subtitle Plot subtitle.  Default \code{NULL}.
#' @param legend_title Legend title.  Default \code{"Gene"}.
#' @param show_legend Logical; show color legend?  Default \code{TRUE}.
#' @param output_dir Directory to save a PDF.  If \code{NULL} (default) the
#'   plot is printed to the active graphics device.
#' @param object_name Optional prefix used when building the output file name.
#' @param subset_name Optional second prefix for the output file name.
#'
#' @return A \code{ggplot} object (invisibly when \code{output_dir} is set).
#'
#' @export
PlotMultiFeature <- function(
    seurat_object,
    features,
    colors             = NULL,
    null_color         = "grey70",
    threshold          = 0.0125,
    plot_mode          = c("gene", "intersection", "scale"),
    pct_cutoff         = 0.3,
    intersection_label = "Co-expressing",
    assay              = "RNA",
    layer              = "data",
    reduction          = "umap",
    pt.size            = 0.5,
    alpha              = 0.9,
    shape              = 16,
    title              = NULL,
    subtitle           = NULL,
    legend_title       = "Gene",
    show_legend        = TRUE,
    output_dir         = NULL,
    object_name        = "",
    subset_name        = ""
) {

  plot_mode <- match.arg(plot_mode)

  # ── input validation ────────────────────────────────────────────────────
  if (!is.character(features) || length(features) == 0L)
    stop("'features' must be a non-empty character vector of gene names.")

  if (plot_mode == "scale" && length(features) > 1L) {
    message("scSidekick: plot_mode = 'scale' uses only the first feature ('",
            features[1L], "').  Pass a single gene to suppress this message.")
    features <- features[1L]
  }

  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("SeuratObject is required. Install with: install.packages('SeuratObject')")

  # ── embedding coordinates ───────────────────────────────────────────────
  emb <- tryCatch(
    SeuratObject::Embeddings(seurat_object, reduction = reduction),
    error = function(e)
      stop("Reduction '", reduction, "' not found in seurat_object. ",
           "Available reductions: ",
           paste(SeuratObject::Reductions(seurat_object), collapse = ", "))
  )
  frame    <- as.data.frame(emb[, 1:2, drop = FALSE])
  dim_labs <- colnames(emb)[1:2]
  colnames(frame) <- c("dim1", "dim2")

  # ── validate features against the assay ────────────────────────────────
  all_genes <- tryCatch({
    rownames(seurat_object[[assay]])          # works for both v3 & v5 assay objects
  }, error = function(e) {
    tryCatch(
      rownames(Seurat::GetAssayData(seurat_object, assay = assay, slot = "data")),
      error = function(e2) NULL
    )
  })

  if (!is.null(all_genes)) {
    missing_f <- setdiff(features, all_genes)
    if (length(missing_f) > 0L) {
      message("scSidekick: Features not found in assay '", assay,
              "', skipping: ", paste(missing_f, collapse = ", "))
      features <- intersect(features, all_genes)
    }
    if (length(features) == 0L)
      stop("No valid features found in assay '", assay, "'.")
  }

  # ── expression matrix (features × cells) ───────────────────────────────
  expr_mat <- .get_layer_data(seurat_object,
                              assay    = assay,
                              layer    = layer,
                              features = features)

  # Align to embedding cells by name; fall back to positional alignment
  common_cells <- intersect(rownames(frame), colnames(expr_mat))
  if (length(common_cells) == 0L) {
    message("scSidekick: Cell names don't overlap between the embedding and the ",
            "expression matrix - aligning by position.")
    n_use    <- min(nrow(frame), ncol(expr_mat))
    frame    <- frame[seq_len(n_use), , drop = FALSE]
    expr_mat <- expr_mat[, seq_len(n_use), drop = FALSE]
  } else {
    frame    <- frame[common_cells, , drop = FALSE]
    expr_mat <- expr_mat[, common_cells, drop = FALSE]
  }

  # ── per-gene thresholds ─────────────────────────────────────────────────
  if (length(threshold) == 1L) {
    thresh_vec <- setNames(rep(as.numeric(threshold), length(features)), features)
  } else {
    thresh_vec <- as.numeric(threshold[features])
    thresh_vec[is.na(thresh_vec)] <- as.numeric(threshold[[1L]])
    names(thresh_vec) <- features
  }

  # ── colors ──────────────────────────────────────────────────────────────
  if (is.null(colors)) {
    colors <- stats::setNames(Nour_pal("spectrum")(length(features)), features)
  } else if (is.null(names(colors))) {
    colors <- stats::setNames(rep_len(colors, length(features)), features)
  } else {
    missing_c <- setdiff(features, names(colors))
    if (length(missing_c) > 0L) {
      extra_cols <- stats::setNames(Nour_pal("spectrum")(length(missing_c)), missing_c)
      colors     <- c(colors, extra_cols)
    }
  }

  # ── plot construction ───────────────────────────────────────────────────

  if (plot_mode == "scale") {
    # ---- continuous gradient (single gene) --------------------------------
    gene          <- features[1L]
    frame$expr_val <- as.numeric(expr_mat[gene, ])
    above          <- frame$expr_val >= thresh_vec[[gene]]
    f_null         <- frame[!above, , drop = FALSE]
    f_expr         <- frame[ above, , drop = FALSE]
    high_col       <- unname(colors[gene])

    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = dim1, y = dim2)) +
      ggplot2::geom_point(data   = f_null,
                          colour = null_color, size  = pt.size,
                          shape  = shape,      alpha = alpha) +
      ggplot2::geom_point(data        = f_expr,
                          ggplot2::aes(colour = expr_val),
                          size        = pt.size, shape = shape,
                          alpha       = alpha,   show.legend = show_legend) +
      ggplot2::scale_colour_gradient(low  = null_color,
                                     high = high_col,
                                     name = gene) +
      ggplot2::labs(title    = title,
                    subtitle = subtitle,
                    x        = dim_labs[1],
                    y        = dim_labs[2]) +
      theme_NourMin()

  } else if (plot_mode == "intersection") {
    # ---- highlight co-expressing cells ------------------------------------
    # Build binary (ncells × ngenes) matrix
    bin_mat <- vapply(features,
                      function(g) as.numeric(expr_mat[g, ]) >= thresh_vec[[g]],
                      logical(nrow(frame)))
    if (!is.matrix(bin_mat))
      bin_mat <- matrix(bin_mat, ncol = 1L, dimnames = list(NULL, features))

    gene_sums  <- rowSums(bin_mat)
    required_n <- ceiling(length(features) * pct_cutoff)

    frame$Plot.Status <- factor(
      ifelse(gene_sums >= required_n, intersection_label, "None"),
      levels = c(intersection_label, "None")
    )
    f_null <- frame[frame$Plot.Status == "None", , drop = FALSE]
    f_pos  <- frame[frame$Plot.Status != "None", , drop = FALSE]

    int_colors <- setNames(unname(colors[1L]), intersection_label)

    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = dim1, y = dim2)) +
      ggplot2::geom_point(data   = f_null,
                          colour = null_color, size  = pt.size,
                          shape  = shape,      alpha = alpha) +
      ggplot2::geom_point(data        = f_pos,
                          ggplot2::aes(colour = Plot.Status),
                          size        = pt.size, shape = shape,
                          alpha       = alpha,   show.legend = show_legend) +
      ggplot2::scale_colour_manual(values = int_colors) +
      ggplot2::labs(title    = title,
                    subtitle = subtitle,
                    x        = dim_labs[1],
                    y        = dim_labs[2]) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 4L),
        title        = legend_title)) +
      theme_NourMin()

  } else {
    # ---- "gene" mode: color by last expressed gene ------------------------
    # Iterating features in order means the last gene overwrites earlier
    # assignments for multi-expressing cells - intentional "last wins" logic
    # that makes cells unique to the rarest gene most visible.
    frame$Plot.Status <- "None"
    for (g in features) {
      expressed               <- as.numeric(expr_mat[g, ]) >= thresh_vec[[g]]
      frame$Plot.Status[expressed] <- g
    }
    frame$Plot.Status <- factor(frame$Plot.Status,
                                levels = c(features, "None"))

    f_null <- frame[frame$Plot.Status == "None", , drop = FALSE]
    f_pos  <- frame[frame$Plot.Status != "None", , drop = FALSE]

    gene_colors <- colors[features]
    na_idx <- is.na(gene_colors)
    if (any(na_idx)) {
      fill_pal            <- Nour_pal("all")
      gene_colors[na_idx] <- fill_pal[seq_len(sum(na_idx))]
    }
    names(gene_colors) <- features

    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = dim1, y = dim2)) +
      ggplot2::geom_point(data   = f_null,
                          colour = null_color, size  = pt.size,
                          shape  = shape,      alpha = alpha) +
      ggplot2::geom_point(data        = f_pos,
                          ggplot2::aes(colour = Plot.Status),
                          size        = pt.size, shape = shape,
                          alpha       = alpha,   show.legend = show_legend) +
      ggplot2::scale_colour_manual(values = gene_colors, drop = FALSE) +
      ggplot2::labs(title    = title,
                    subtitle = subtitle,
                    x        = dim_labs[1],
                    y        = dim_labs[2]) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 4L),
        title        = legend_title)) +
      theme_NourMin()
  }

  # Hide legend if requested (belt-and-suspenders with show.legend above)
  if (!show_legend)
    p <- p + ggplot2::theme(legend.position = "none")

  # ── output ──────────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    gene_tag <- paste(features[seq_len(min(3L, length(features)))],
                      collapse = "_")
    parts <- c(
      if (nchar(object_name) > 0L) object_name,
      if (nchar(subset_name) > 0L) subset_name,
      "MultiFeature",
      gene_tag
    )
    fname <- gsub("[^A-Za-z0-9._-]", "_", paste(parts, collapse = "_"))
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))
    ggplot2::ggsave(fpath, plot = p, width = 6, height = 5)
    message("scSidekick: Plot saved to ", fpath)
    return(invisible(p))
  }

  print(p)
  invisible(p)
}
