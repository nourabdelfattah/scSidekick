# =============================================================================
# scSidekick gene expression filters
#
# FilterGenesByPct  - keep genes expressed above a % threshold across groups
# PlotPctHeatmap    - heatmap of per-group expression percentages
# =============================================================================

#' Filter genes by percent expression across groups
#'
#' Returns the subset of `gene_list` where the fraction of cells expressing
#' each gene (above `expr_threshold`) meets `pct_threshold` in either any or
#' all groups. Works with Assay, Assay5, and BPCells-backed objects.
#'
#' @param seurat_object A Seurat object.
#' @param gene_list Character vector of gene names to test.
#' @param group_by Character. Metadata column defining the groups.
#' @param assay Character or `NULL`. Assay to pull expression from. Defaults
#'   to `DefaultAssay()`.
#' @param layer Character. Layer/slot to use (`"data"`, `"counts"`, etc.).
#'   Default `"data"`.
#' @param expr_threshold Numeric. A cell is considered "expressing" the gene
#'   if its value is strictly above this threshold. Default `0`.
#' @param pct_threshold Numeric. Minimum percent (0–100) of cells in a group
#'   that must express the gene. Default `20`.
#' @param mode Character. `"any_group"` (default) keeps a gene if it passes
#'   the threshold in at least one group; `"all_groups"` requires every group
#'   to pass.
#'
#' @return A named list:
#' \describe{
#'   \item{`genes`}{Character vector of genes that passed the filter.}
#'   \item{`pct_table`}{Data frame (genes × groups) of percent-expressed
#'     values, with an extra `gene` column. Pass to [PlotPctHeatmap()].}
#' }
#' @seealso [PlotPctHeatmap()]
#' @export
FilterGenesByPct <- function(seurat_object,
                              gene_list,
                              group_by,
                              assay          = NULL,
                              layer          = "data",
                              expr_threshold = 0,
                              pct_threshold  = 20,
                              mode           = c("any_group", "all_groups")) {

  mode  <- match.arg(mode)
  assay <- assay %||% Seurat::DefaultAssay(seurat_object)

  if (!group_by %in% colnames(seurat_object@meta.data))
    stop("'", group_by, "' not found in meta.data.")

  # GetAssayData uses `layer` in Seurat v5, `slot` in v3/v4
  mat <- tryCatch(
    Seurat::GetAssayData(seurat_object, assay = assay, layer = layer),
    error = function(e)
      Seurat::GetAssayData(seurat_object, assay = assay, slot = layer)
  )

  present_genes <- intersect(gene_list, rownames(mat))
  missing_genes <- setdiff(gene_list, rownames(mat))
  if (length(missing_genes))
    warning(length(missing_genes), " gene(s) not found in assay '", assay,
            "': ", paste(missing_genes, collapse = ", "))
  if (!length(present_genes))
    return(list(genes = character(0), pct_table = NULL))

  # Subset to the genes we care about first, then standardise to dgCMatrix.
  # This keeps BPCells-backed objects manageable — the gene slice is small.
  mat_sub <- mat[present_genes, , drop = FALSE]
  if (!inherits(mat_sub, c("dgCMatrix", "matrix")))
    mat_sub <- methods::as(mat_sub, "dgCMatrix")

  groups         <- as.character(seurat_object@meta.data[[group_by]])
  group_levels   <- unique(groups)
  cells_by_group <- lapply(group_levels, function(g) which(groups == g))
  names(cells_by_group) <- group_levels

  # Vectorised over all genes at once per group — avoids per-gene R loop
  pct_mat <- do.call(cbind, lapply(cells_by_group, function(idx) {
    grp <- mat_sub[, idx, drop = FALSE]
    Matrix::rowSums(grp > expr_threshold) / length(idx) * 100
  }))
  colnames(pct_mat) <- group_levels

  pct_df       <- as.data.frame(pct_mat)
  pct_df$gene  <- rownames(pct_df)

  keep <- if (mode == "all_groups") {
    rownames(pct_df)[apply(pct_mat, 1, function(x) all(x > pct_threshold))]
  } else {
    rownames(pct_df)[apply(pct_mat, 1, function(x) any(x > pct_threshold))]
  }

  message("FilterGenesByPct: ", length(keep), " / ", length(present_genes),
          " genes passed (mode = '", mode, "', pct_threshold = ", pct_threshold, "%)")

  list(genes = keep, pct_table = pct_df)
}

#' Heatmap of per-group gene expression percentages
#'
#' Visualizes the `pct_table` returned by [FilterGenesByPct()] as a tile
#' heatmap — genes on the y-axis, groups on the x-axis, color encoding the
#' percent of cells expressing each gene.
#'
#' @param pct_table Data frame as returned by [FilterGenesByPct()] (gene names
#'   as rownames or in a `gene` column; remaining columns are group names).
#' @param title Character. Plot title. Default `"% cells expressing each gene"`.
#' @param output_dir Character or `NULL`. Directory to save a PDF. `NULL`
#'   returns the plot.
#' @param file_name Character or `NULL`. Base name (no extension) for the
#'   saved PDF. Defaults to `"PctHeatmap"`.
#'
#' @return A ggplot2 object (invisibly when saved).
#' @seealso [FilterGenesByPct()]
#' @export
PlotPctHeatmap <- function(pct_table,
                            title      = "% cells expressing each gene",
                            output_dir = NULL,
                            file_name  = NULL) {

  df <- pct_table
  if (!"gene" %in% colnames(df))
    df$gene <- rownames(df)

  group_cols <- setdiff(colnames(df), "gene")
  if (!length(group_cols))
    stop("pct_table has no group columns.")

  df_long <- tidyr::pivot_longer(df, cols = dplyr::all_of(group_cols),
                                  names_to = "group", values_to = "pct")

  p <- ggplot2::ggplot(df_long,
                       ggplot2::aes(x = group, y = gene, fill = pct)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(name = "% expressed", option = "plasma",
                                   begin = 0, end = 1) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_NourMin() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8)
    )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fname <- if (!is.null(file_name) && nzchar(file_name)) file_name else "PctHeatmap"
    fpath <- file.path(output_dir, paste0(gsub("[^A-Za-z0-9._-]", "_", fname), ".pdf"))
    n_genes  <- length(unique(df_long$gene))
    n_groups <- length(unique(df_long$group))
    ggplot2::ggsave(fpath, p,
                    width     = max(4, n_groups * 0.7 + 2.5),
                    height    = max(4, n_genes  * 0.25 + 2.0),
                    limitsize = FALSE)
    message("scSidekick: Saved to ", fpath)
    return(invisible(p))
  }

  p
}
