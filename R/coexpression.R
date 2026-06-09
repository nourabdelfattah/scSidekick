# =============================================================================
# scSidekick gene coexpression / correlation utilities
#
# CorrelateGene   - correlate one gene against all others within each cell type
#                   at single-cell, pseudobulk (donor), or group-average level.
#                   Optional compare_by: run separately within each group level
#                   (e.g. Sex) so you can compare coexpression patterns across
#                   groups within cell types.
# SplitByGene     - add a cell-level metadata column splitting cells into
#                   high / low (or expressed / zero) groups per cell type,
#                   ready to plug into RunGSEA_pseudobulk or FindMarkers
# PlotCorrelation - per-cell-type volcano plots assembled with wrap_plots()
# =============================================================================

# ── Internal helpers ──────────────────────────────────────────────────────────

.cor_pval <- function(r, n) {
  t_stat <- r * sqrt(n - 2) / sqrt(pmax(1 - r^2, .Machine$double.eps))
  2 * stats::pt(-abs(t_stat), df = n - 2)
}

.get_expr_mat <- function(seurat_object, assay, layer) {
  tryCatch(
    Seurat::GetAssayData(seurat_object, assay = assay, layer = layer),
    error = function(e)
      Seurat::GetAssayData(seurat_object, assay = assay, slot = layer)
  )
}

# Core computation shared by all three levels; returns a tidy data frame or NULL.
.cor_one_result <- function(goi_vec, gene_mat, gene_names, ct, method, n) {
  if (length(goi_vec) < 3 || ncol(gene_mat) == 0) return(NULL)
  r_vec    <- as.numeric(cor(goi_vec, gene_mat, method = method))
  p_vec    <- .cor_pval(r_vec, n)
  padj_vec <- stats::p.adjust(p_vec, method = "BH")
  data.frame(gene = gene_names, cell_type = ct,
             r = r_vec, p = p_vec, padj = padj_vec, n = n,
             stringsAsFactors = FALSE)
}

# ── Level helpers ─────────────────────────────────────────────────────────────

.cor_gene_sc <- function(seurat_object, gene, group_by, assay, layer,
                          method, expressed_only, min_pct, min_cells) {
  mat        <- .get_expr_mat(seurat_object, assay, layer)
  cell_types <- sort(unique(as.character(seurat_object@meta.data[[group_by]])))

  results <- lapply(cell_types, function(ct) {
    cells    <- rownames(seurat_object@meta.data)[
      seurat_object@meta.data[[group_by]] == ct]
    sub      <- mat[, cells, drop = FALSE]
    goi_vals <- as.numeric(as.matrix(sub[gene, , drop = FALSE]))

    if (expressed_only) {
      keep <- which(goi_vals > 0)
      if (length(keep) < min_cells) {
        message("  Skipping ", ct, ": only ", length(keep),
                " cells express ", gene)
        return(NULL)
      }
      sub      <- sub[, keep, drop = FALSE]
      goi_vals <- goi_vals[keep]
    } else if (ncol(sub) < min_cells) return(NULL)

    pct_expr   <- Matrix::rowSums(sub > 0) / ncol(sub) * 100
    keep_genes <- names(pct_expr)[pct_expr >= min_pct & names(pct_expr) != gene]
    if (!length(keep_genes)) return(NULL)

    gene_mat <- as.matrix(t(sub[keep_genes, , drop = FALSE]))
    .cor_one_result(goi_vals, gene_mat, keep_genes, ct, method, ncol(sub))
  })

  do.call(rbind, Filter(Negate(is.null), results))
}

.cor_gene_pseudobulk <- function(seurat_object, gene, group_by, sample_column,
                                  assay, method, min_cells, min_donors) {
  for (pkg in c("edgeR"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for level = 'pseudobulk'.")

  seurat_object[[".pb_cor"]] <- paste(
    seurat_object@meta.data[[sample_column]],
    seurat_object@meta.data[[group_by]], sep = "||")

  agg <- tryCatch(
    Seurat::AggregateExpression(seurat_object, assays = assay,
                                group.by = ".pb_cor", layer = "counts")[[assay]],
    error = function(e)
      Seurat::AggregateExpression(seurat_object, assays = assay,
                                  group.by = ".pb_cor", slot  = "counts")[[assay]]
  )
  agg  <- as.matrix(agg)
  samp <- data.frame(.pb       = colnames(agg),
                     sample    = sub("\\|\\|.*",  "", colnames(agg)),
                     cell_type = sub(".*\\|\\|",  "", colnames(agg)),
                     stringsAsFactors = FALSE)

  y      <- edgeR::DGEList(round(agg))
  y      <- edgeR::calcNormFactors(y)
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)

  results <- lapply(sort(unique(samp$cell_type)), function(ct) {
    ks <- samp$.pb[samp$cell_type == ct]
    if (length(ks) < min_donors) {
      message("  Skipping ", ct, ": only ", length(ks), " donors"); return(NULL)
    }
    sub <- logcpm[, ks, drop = FALSE]
    if (!gene %in% rownames(sub)) return(NULL)
    other <- setdiff(rownames(sub), gene)
    .cor_one_result(as.numeric(as.matrix(sub[gene, , drop = FALSE])),
                    t(sub[other, , drop = FALSE]),
                    other, ct, method, length(ks))
  })

  do.call(rbind, Filter(Negate(is.null), results))
}

.cor_gene_group <- function(seurat_object, gene, group_by, group_column,
                             assay, layer, method, min_groups) {
  mat   <- .get_expr_mat(seurat_object, assay, layer)
  md    <- seurat_object@meta.data
  combo <- paste(as.character(md[[group_column]]),
                 as.character(md[[group_by]]), sep = "||")
  combos <- unique(combo)

  avg_mat <- do.call(cbind, lapply(combos, function(co)
    Matrix::rowMeans(mat[, which(combo == co), drop = FALSE])))
  colnames(avg_mat) <- combos
  avg_mat <- as.matrix(avg_mat)

  results <- lapply(sort(unique(as.character(md[[group_by]]))), function(ct) {
    ct_cols <- grep(paste0("\\|\\|", ct, "$"), combos, value = TRUE)
    if (length(ct_cols) < min_groups) {
      message("  Skipping ", ct, ": only ", length(ct_cols), " groups")
      return(NULL)
    }
    sub <- avg_mat[, ct_cols, drop = FALSE]
    if (!gene %in% rownames(sub)) return(NULL)
    other <- setdiff(rownames(sub), gene)
    .cor_one_result(as.numeric(as.matrix(sub[gene, , drop = FALSE])),
                    t(sub[other, , drop = FALSE]),
                    other, ct, method, length(ct_cols))
  })

  do.call(rbind, Filter(Negate(is.null), results))
}

# =============================================================================
# CorrelateGene
# =============================================================================

#' Correlate one gene against all others within each cell type
#'
#' Computes Spearman (or Pearson) correlation between a gene of interest (GOI)
#' and every other gene, separately within each cell type. Three resolution
#' levels are supported:
#'
#' * **`"sc"`**: cell level. Restricts to cells where GOI > 0
#'   (`expressed_only = TRUE`) by default to avoid co-absence artifacts from
#'   dropout. Exploratory; validate hits at pseudobulk level.
#' * **`"pseudobulk"`**: donor level. Aggregates raw counts per donor x
#'   cell type, converts to edgeR log-CPM, then runs Spearman across donors.
#'   Statistically most valid.
#' * **`"group"`**: averages log-normalized expression within each
#'   `group_column x cell_type` combination, then correlates across groups.
#'
#' When `compare_by` is supplied, the analysis is run separately within each
#' level of that variable (e.g. `compare_by = "Sex"` produces independent
#' correlation tables for males and females). Results gain a `compare_group`
#' column and [PlotCorrelation()] will lay the groups out side by side.
#'
#' @param seurat_object A Seurat object.
#' @param gene Character. The gene of interest (GOI).
#' @param group_by Character. Cell-type identity column.
#' @param level Character. `"sc"` (default), `"pseudobulk"`, or `"group"`.
#' @param compare_by Character or `NULL`. Optional metadata column whose levels
#'   define independent comparison groups (e.g. `"Sex"`, `"Diagnosis"`).
#'   When supplied, correlation is computed separately within each level so you
#'   can ask which genes are correlated with GOI in group A but not group B.
#' @param sample_column Character. Donor/sample ID. Required for `"pseudobulk"`.
#' @param group_column Character. Condition column. Required for `"group"`.
#' @param assay Character or `NULL`. Assay to use.
#' @param layer Character. Layer/slot for SC and group levels. Default `"data"`.
#' @param method Character. `"spearman"` (default) or `"pearson"`.
#' @param expressed_only Logical. SC level: restrict to GOI > 0 cells? Default
#'   `TRUE`.
#' @param min_pct Numeric. SC level: minimum % of (GOI-expressing) cells
#'   that must express a tested gene. Default `10`.
#' @param min_cells Integer. SC level: minimum cells per cell type. Default
#'   `20`.
#' @param min_donors Integer. Pseudobulk level: minimum donors. Default `5`.
#' @param min_groups Integer. Group level: minimum groups. Default `3`.
#'
#' @return A data frame: `gene`, `cell_type`, `r`, `p`, `padj`, `n`, `goi`,
#'   `level`, `method`, and (when `compare_by` is set) `compare_group`.
#'   Sorted by cell type then descending `|r|`.
#' @seealso [SplitByGene()], [PlotCorrelation()]
#' @export
CorrelateGene <- function(seurat_object,
                           gene,
                           group_by,
                           level          = c("sc", "pseudobulk", "group"),
                           compare_by     = NULL,
                           sample_column  = NULL,
                           group_column   = NULL,
                           assay          = NULL,
                           layer          = "data",
                           method         = c("spearman", "pearson"),
                           expressed_only = TRUE,
                           min_pct        = 10,
                           min_cells      = 20L,
                           min_donors     = 5L,
                           min_groups     = 3L) {

  level  <- match.arg(level)
  method <- match.arg(method)
  assay  <- assay %||% Seurat::DefaultAssay(seurat_object)

  if (!gene %in% rownames(seurat_object))
    stop("Gene '", gene, "' not found.")
  if (!group_by %in% colnames(seurat_object@meta.data))
    stop("'", group_by, "' not found in meta.data.")
  if (!is.null(compare_by) && !compare_by %in% colnames(seurat_object@meta.data))
    stop("'", compare_by, "' not found in meta.data.")

  .run_one <- function(seu) {
    switch(level,
      sc = .cor_gene_sc(seu, gene, group_by, assay, layer,
                         method, expressed_only, min_pct, min_cells),
      pseudobulk = {
        if (is.null(sample_column))
          stop("'sample_column' is required for level = 'pseudobulk'.")
        .cor_gene_pseudobulk(seu, gene, group_by, sample_column,
                              assay, method, min_cells, min_donors)
      },
      group = {
        if (is.null(group_column))
          stop("'group_column' is required for level = 'group'.")
        .cor_gene_group(seu, gene, group_by, group_column,
                         assay, layer, method, min_groups)
      }
    )
  }

  message("CorrelateGene: '", gene, "'  level = '", level, "'",
          if (!is.null(compare_by)) paste0("  compare_by = '", compare_by, "'") else "")

  if (!is.null(compare_by)) {
    grp_levels <- sort(unique(stats::na.omit(
      as.character(seurat_object@meta.data[[compare_by]]))))
    res_list <- lapply(grp_levels, function(grp) {
      cells <- rownames(seurat_object@meta.data)[
        !is.na(seurat_object@meta.data[[compare_by]]) &
        seurat_object@meta.data[[compare_by]] == grp]
      message("  Group: ", grp, " (", length(cells), " cells)")
      r <- .run_one(seurat_object[, cells])
      if (!is.null(r) && nrow(r)) r$compare_group <- grp
      r
    })
    res <- do.call(rbind, Filter(Negate(is.null), res_list))
  } else {
    res <- .run_one(seurat_object)
  }

  if (is.null(res) || !nrow(res)) {
    warning("CorrelateGene: no results returned.")
    return(data.frame(gene = character(), cell_type = character(),
                      r = numeric(), p = numeric(), padj = numeric(),
                      n = integer()))
  }

  res$goi    <- gene
  res$level  <- level
  res$method <- method
  res <- res[order(res$cell_type, -abs(res$r)), ]
  rownames(res) <- NULL

  # Summary
  grp_col <- if (!is.null(compare_by)) "compare_group" else NULL
  by_cols  <- c("cell_type", grp_col)
  splits   <- split(res, res[, by_cols, drop = FALSE])
  message("Significant genes (padj < 0.05) per cell type",
          if (!is.null(compare_by)) " x compare group" else "", ":")
  for (nm in names(splits))
    message("  ", nm, ": ", sum(splits[[nm]]$padj < 0.05, na.rm = TRUE))

  res
}

# =============================================================================
# SplitByGene
# =============================================================================

#' Split cells into expression groups per cell type
#'
#' Adds a metadata column labeling each cell as `"high"` / `"low"` (or
#' `"expressed"` / `"zero"`) for a given gene. The split is computed
#' independently within each cell type. The column can be passed directly as
#' `group_column` in [RunGSEA_pseudobulk()] or used with `Seurat::FindMarkers`.
#'
#' **Split modes**
#' \describe{
#'   \item{`"zero_vs_expressed"`}{GOI = 0 → `"zero"`; GOI > 0 →
#'     `"expressed"`. Cleanest for pseudobulk: maps to donor-level proportion.}
#'   \item{`"quantile"`}{Bottom `quantile_range[1]` → `"low"`; top
#'     `1 - quantile_range[2]` → `"high"`; middle → `NA` (excluded).}
#'   \item{`"expressed_quantile"`}{Same as `"quantile"` but among GOI > 0 cells
#'     only. Zeros → `NA`. Avoids dropout in the low group but can confound
#'     with library size - regress out `nCount_RNA` downstream.}
#' }
#'
#' @param seurat_object A Seurat object.
#' @param gene Character. Gene to split on.
#' @param group_by Character. Cell-type identity column.
#' @param split Character. `"zero_vs_expressed"` (default), `"quantile"`, or
#'   `"expressed_quantile"`.
#' @param quantile_range Numeric length-2. Quantile bounds for the two
#'   threshold modes. Default `c(0.25, 0.75)`.
#' @param assay Character or `NULL`. Assay to use.
#' @param layer Character. Layer/slot. Default `"data"`.
#' @param col_name Character or `NULL`. New metadata column name. Defaults to
#'   `"<gene>_split"`.
#'
#' @return Seurat object with a new metadata column. Prints a count table.
#' @seealso [CorrelateGene()], [RunGSEA_pseudobulk()]
#' @export
SplitByGene <- function(seurat_object,
                         gene,
                         group_by,
                         split          = c("zero_vs_expressed", "quantile",
                                            "expressed_quantile"),
                         quantile_range = c(0.25, 0.75),
                         assay          = NULL,
                         layer          = "data",
                         col_name       = NULL) {

  split    <- match.arg(split)
  assay    <- assay %||% Seurat::DefaultAssay(seurat_object)
  col_name <- col_name %||% paste0(make.names(gene), "_split")

  if (!gene %in% rownames(seurat_object))
    stop("Gene '", gene, "' not found.")
  if (!group_by %in% colnames(seurat_object@meta.data))
    stop("'", group_by, "' not found in meta.data.")
  if (length(quantile_range) != 2 ||
      any(quantile_range < 0 | quantile_range > 1) ||
      quantile_range[1] >= quantile_range[2])
    stop("quantile_range must be two values in (0,1) with [1] < [2].")

  mat     <- .get_expr_mat(seurat_object, assay, layer)
  goi_all <- stats::setNames(
    as.numeric(as.matrix(mat[gene, , drop = FALSE])),
    colnames(mat))

  cell_types <- as.character(seurat_object@meta.data[[group_by]])
  split_vec  <- stats::setNames(rep(NA_character_, nrow(seurat_object@meta.data)),
                                rownames(seurat_object@meta.data))

  for (ct in unique(cell_types)) {
    cells <- rownames(seurat_object@meta.data)[cell_types == ct]
    vals  <- goi_all[cells]

    lbl <- switch(split,

      zero_vs_expressed = ifelse(vals > 0, "expressed", "zero"),

      quantile = {
        lo  <- stats::quantile(vals, quantile_range[1])
        hi  <- stats::quantile(vals, quantile_range[2])
        l   <- rep(NA_character_, length(vals))
        l[vals <= lo] <- "low"
        l[vals >= hi] <- "high"
        l
      },

      expressed_quantile = {
        expr_vals <- vals[vals > 0]
        l         <- rep(NA_character_, length(vals))
        if (length(expr_vals) >= 10) {
          lo <- stats::quantile(expr_vals, quantile_range[1])
          hi <- stats::quantile(expr_vals, quantile_range[2])
          l[vals > 0 & vals <= lo] <- "low"
          l[vals >= hi]            <- "high"
        } else {
          message("  Skipping ", ct, ": only ", length(expr_vals),
                  " expressing cells")
        }
        l
      }
    )

    split_vec[cells] <- lbl
  }

  seurat_object@meta.data[[col_name]] <- split_vec
  tab <- table(split_vec, cell_types, useNA = "no",
               dnn = c(col_name, group_by))
  message("SplitByGene: '", col_name, "'  (split = '", split, "')")
  print(tab)
  seurat_object
}

# =============================================================================
# PlotCorrelation
# =============================================================================

#' Volcano plots of CorrelateGene results
#'
#' Draws one panel per cell type (and per comparison group when `compare_by`
#' was used in [CorrelateGene()]) assembled with [patchwork::wrap_plots()].
#' Each panel shows the top `n_pos` most positively and top `n_neg` most
#' negatively correlated genes so that both directions are always represented.
#'
#' @param cor_result Data frame returned by [CorrelateGene()].
#' @param cell_type Character vector or `NULL`. Cell types to plot (all by
#'   default).
#' @param n_pos Integer. Top positive genes to show per panel. Default `15`.
#' @param n_neg Integer. Top negative genes to show per panel. Default `15`.
#' @param label_n Integer. Significant genes to label per direction per panel.
#'   Default `5` (5 positive + 5 negative).
#' @param padj_cutoff Numeric. Significance threshold. Default `0.05`.
#' @param r_cutoff Numeric or `NULL`. Optional ±r dotted reference lines.
#' @param ncol Integer or `NULL`. Columns in the assembled figure. `NULL`
#'   auto-selects: `n_compare_groups` when `compare_by` was used, else 3.
#' @param output_dir Character or `NULL`. Directory to save a PDF.
#' @param file_name Character or `NULL`. PDF base name.
#'
#' @return A patchwork object (invisibly when saved).
#' @seealso [CorrelateGene()]
#' @export
PlotCorrelation <- function(cor_result,
                             cell_type   = NULL,
                             n_pos       = 15L,
                             n_neg       = 15L,
                             label_n     = 5L,
                             padj_cutoff = 0.05,
                             r_cutoff    = NULL,
                             ncol        = NULL,
                             output_dir  = NULL,
                             file_name   = NULL) {

  df <- cor_result
  if (!is.null(cell_type))
    df <- df[df$cell_type %in% cell_type, , drop = FALSE]
  if (!nrow(df)) stop("No results for the selected cell type(s).")
  df <- df[!is.na(df$r) & !is.na(df$padj), ]

  has_compare <- "compare_group" %in% colnames(df)
  panel_cols  <- if (has_compare) c("cell_type", "compare_group") else "cell_type"

  # Build panel key data frame (one row per plot).
  # Coerce to plain data.frame so nrow/[[ behave consistently regardless of
  # whether df is a tibble, data.table, or grouped_df.
  panel_keys <- as.data.frame(
    unique(df[, panel_cols, drop = FALSE]),
    stringsAsFactors = FALSE)
  if (has_compare) {
    panel_keys <- panel_keys[
      order(panel_keys$cell_type, panel_keys$compare_group), , drop = FALSE]
  } else {
    panel_keys <- panel_keys[order(panel_keys$cell_type), , drop = FALSE]
  }
  rownames(panel_keys) <- NULL
  n_panels <- nrow(panel_keys)

  # Select top n_pos positive + top n_neg negative per panel
  df_plot <- do.call(rbind, lapply(seq_len(n_panels), function(i) {
    d <- df
    for (col in panel_cols) d <- d[d[[col]] == panel_keys[[col]][i], ]
    pos <- d[d$r >  0, ][order(-d[d$r >  0, "r"]), ][seq_len(min(n_pos, sum(d$r >  0))), ]
    neg <- d[d$r <= 0, ][order( d[d$r <= 0, "r"]), ][seq_len(min(n_neg, sum(d$r <= 0))), ]
    rbind(pos, neg)
  }))

  df_plot$significant <- df_plot$padj < padj_cutoff
  df_plot$direction   <- ifelse(df_plot$r > 0, "positive", "negative")
  df_plot$log10padj   <- -log10(pmax(df_plot$padj, .Machine$double.eps))

  # One volcano per panel key
  plots <- lapply(seq_len(n_panels), function(i) {
    d <- df_plot
    for (col in panel_cols) d <- d[d[[col]] == panel_keys[[col]][i], ]

    # Labels: top label_n significant per direction
    sig     <- d[d$significant, ]
    lbl_pos <- sig[sig$r >  0, ][order(-sig[sig$r >  0, "r"]), ][seq_len(min(label_n, sum(sig$r >  0))), ]
    lbl_neg <- sig[sig$r <= 0, ][order( sig[sig$r <= 0, "r"]), ][seq_len(min(label_n, sum(sig$r <= 0))), ]
    lbl_df  <- rbind(lbl_pos, lbl_neg)

    ptitle <- panel_keys[["cell_type"]][i]
    if (has_compare)
      ptitle <- paste0(ptitle, "\n", panel_keys[["compare_group"]][i])

    p <- ggplot2::ggplot(d, ggplot2::aes(x = r, y = log10padj,
                                          color = direction)) +
      ggplot2::geom_point(ggplot2::aes(size  = abs(r),
                                        alpha = significant)) +
      ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.25),
                                   guide  = "none") +
      ggplot2::scale_color_manual(
        values = c(positive = "#B5232E", negative = "#2C5F8A"),
        name   = "Direction") +
      ggplot2::scale_size_continuous(range = c(0.5, 3), guide = "none") +
      ggrepel::geom_text_repel(
        data         = lbl_df,
        ggplot2::aes(label = gene),
        size         = 2.8,
        max.overlaps = 20,
        show.legend  = FALSE
      ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                           color = "gray60", linewidth = 0.4) +
      ggplot2::geom_hline(yintercept = -log10(padj_cutoff),
                           linetype = "dashed", color = "gray60",
                           linewidth = 0.4) +
      ggplot2::labs(
        x     = paste0(unique(d$method)[1], " r"),
        y     = expression(-log[10](p[adj])),
        title = ptitle
      ) +
      theme_NourMin() +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(size = 9, face = "bold",
                                             lineheight = 1.1),
        axis.title  = ggplot2::element_text(size = 8),
        axis.text   = ggplot2::element_text(size = 7),
        # Extra bottom + right margin so rotated/long axis text isn't clipped
        plot.margin = ggplot2::margin(4, 10, 10, 4, "mm")
      )

    if (!is.null(r_cutoff))
      p <- p +
        ggplot2::geom_vline(xintercept =  r_cutoff, linetype = "dotted",
                             color = "gray40", linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = -r_cutoff, linetype = "dotted",
                             color = "gray40", linewidth = 0.4)
    p
  })

  # Layout: when compare_by was used, arrange as cell_types x compare_groups
  n_panels <- length(plots)
  n_col_out <- ncol %||% (if (has_compare)
    length(unique(df_plot$compare_group)) else min(3L, n_panels))

  combined <- patchwork::wrap_plots(plots, ncol = n_col_out,
                                     guides = "collect") +
    patchwork::plot_annotation(
      title = paste0("Coexpression with ", unique(df_plot$goi)[1],
                     " - ", unique(df_plot$level)[1], " level"),
      theme = ggplot2::theme(
        plot.title  = ggplot2::element_text(size = 12, face = "bold"),
        plot.margin = ggplot2::margin(4, 4, 4, 4, "mm")
      )
    )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fname <- if (!is.null(file_name) && nzchar(file_name)) file_name else {
      parts <- c(
        "CorrelateGene",
        unique(df_plot$goi)[1],
        unique(df_plot$level)[1],
        unique(df_plot$method)[1],
        if (has_compare) paste(sort(unique(df_plot$compare_group)), collapse = "_vs_") else NULL
      )
      paste(parts, collapse = "_")
    }
    fpath <- file.path(output_dir,
                       paste0(gsub("[^A-Za-z0-9._-]", "_", fname), ".pdf"))
    n_row_out <- ceiling(n_panels / n_col_out)
    ggplot2::ggsave(fpath, combined,
                    width     = n_col_out * 5.0,
                    height    = n_row_out * 4.8 + 0.8,
                    limitsize = FALSE)
    message("scSidekick: Saved to ", fpath)

    if (exists(".write_legend_sidecar")) {
      goi    <- unique(df_plot$goi)[1]
      lvl    <- unique(df_plot$level)[1]
      mth    <- unique(df_plot$method)[1]
      .write_legend_sidecar(fpath, paste0(
        "Volcano plots of gene-gene ", mth, " correlation with ", goi,
        " at the ", lvl, " level",
        if (has_compare)
          paste0(", stratified by comparison groups: ",
                 paste(sort(unique(df_plot$compare_group)), collapse = " vs "))
        else "",
        ". Each panel = one cell type",
        if (has_compare) " x comparison group" else "",
        ". X-axis: ", mth, " r. Y-axis: -log10(BH-adjusted p). ",
        "Red = positively correlated with ", goi,
        "; blue = negatively correlated. ",
        "Each panel shows the top ", n_pos,
        " positive and top ", n_neg, " negative genes by |r|. ",
        "Labelled = top ", label_n, " significant genes per direction. ",
        "Dashed lines: r = 0 and padj = ", padj_cutoff, "."
      ))
    }
    return(invisible(combined))
  }

  combined
}
