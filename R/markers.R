# =============================================================================
# scSidekick - Marker Gene Discovery  (markers.R)
#
# Exported:
#   RunWilcoxAUC()  - Seurat-v5 / BPCells-compatible wrapper around
#                     presto::wilcoxauc that bypasses presto's internal
#                     Seurat extraction (which fails on lazy / multi-layer
#                     objects).  Extracts the matrix via .get_layer_data()
#                     first, then calls presto with a plain dgCMatrix so
#                     S3 dispatch lands on wilcoxauc.default cleanly -
#                     no ::: required.
# =============================================================================


#' Fast Wilcoxon Rank-Sum + AUC Marker Test (BPCells / Seurat-v5 Compatible)
#'
#' @description
#' A drop-in replacement for \code{presto::wilcoxauc(seurat_object, ...)} that
#' works reliably with **BPCells-backed** and **Seurat v5 multi-layer** objects.
#'
#' The standard \code{presto::wilcoxauc.Seurat} method extracts the expression
#' matrix internally in a way that can error with lazy matrices or v5 layer
#' syntax.  \code{RunWilcoxAUC} bypasses that path: it pulls the matrix through
#' \code{.get_layer_data()} (which handles BPCells lazy subsetting, Seurat v3
#' \code{slot=} vs v5 \code{layer=}, and sparse → dgCMatrix coercion) and then
#' calls \code{presto::wilcoxauc(matrix, y)} so S3 dispatch lands on
#' \code{wilcoxauc.default} directly - no internal \code{:::} access needed.
#'
#' When \code{split.by} is supplied the test is run independently within each
#' level of that column (e.g. each cell type), and the results are combined into
#' a single data frame with an extra column named after \code{split.by}.
#'
#' @param seurat_object A Seurat object.
#' @param group.by Character. Metadata column whose levels define the
#'   comparison groups (e.g. \code{"Condition"}, \code{"Assignment"}).
#'   When \code{NULL}, \code{\link[Seurat]{Idents}} is used.
#' @param split.by Character or \code{NULL}. Metadata column to split on before
#'   running the test.  The test is run independently within each level
#'   (e.g. \code{split.by = "CellType"} finds condition markers per cell type).
#'   The output gains a column named after \code{split.by} identifying each
#'   level.  Levels with fewer than 2 groups after filtering are skipped with a
#'   message.  Default \code{NULL} (no split).
#' @param groups_use Character vector or \code{NULL}. Restrict the test to
#'   a subset of \code{group.by} levels.  Applied within each split level when
#'   \code{split.by} is set.  \code{NULL} = all levels.
#' @param assay Character. Assay to extract. Default \code{"RNA"}.
#' @param layer Character. Layer / slot to extract. Default \code{"data"}
#'   (log-normalized counts, recommended for AUC statistics).
#'   \code{"counts"} uses raw counts; \code{"scale.data"} uses scaled values
#'   (note: negative values make AUC less interpretable but the test is still
#'   valid).
#' @param caffeinate Logical. When \code{TRUE}, prevents the machine from
#'   sleeping during the run (macOS only; uses \code{caffeinate -i}).
#'   Default \code{FALSE}.
#' @param ... Additional arguments forwarded to \code{\link[presto]{wilcoxauc}}.
#'
#' @return A \code{data.frame} with one row per (feature × group) combination
#'   containing \code{auc}, \code{pval}, \code{padj}, \code{logFC},
#'   \code{pct_in}, \code{pct_out} - the standard \pkg{presto} output.
#'   When \code{split.by} is set, an additional column named after
#'   \code{split.by} identifies which split level each row came from.
#'
#' @seealso \code{\link[presto]{wilcoxauc}}
#'
#' @examples
#' \dontrun{
#' # Basic usage - use active identity
#' markers <- RunWilcoxAUC(SeuratObj)
#'
#' # Specify a metadata column
#' markers <- RunWilcoxAUC(SeuratObj, group.by = "Condition")
#'
#' # Find condition markers within each cell type
#' markers <- RunWilcoxAUC(SeuratObj,
#'                          group.by = "Condition",
#'                          split.by = "CellType")
#'
#' # Filter results for one cell type
#' library(dplyr)
#' markers |> filter(CellType == "Microglia", padj < 0.05, auc > 0.6)
#'
#' # Find cluster markers within each sample
#' markers <- RunWilcoxAUC(SeuratObj,
#'                          group.by = "Cluster",
#'                          split.by = "Sample")
#' }
#'
#' @export
RunWilcoxAUC <- function(
    seurat_object,
    group.by   = NULL,
    split.by   = NULL,
    groups_use = NULL,
    assay      = "RNA",
    layer      = "data",
    caffeinate = FALSE,
    ...
) {
  if (caffeinate) { .caff <- .nk_caffeinate(); on.exit(.nk_decaffeinate(.caff), add = TRUE) }
  if (!requireNamespace("presto", quietly = TRUE))
    stop("Package 'presto' is required.\n",
         "Install with: remotes::install_github('immunogenomics/presto')")

  meta <- seurat_object@meta.data

  # ── Validate split.by ─────────────────────────────────────────────────────
  if (!is.null(split.by) && !split.by %in% colnames(meta))
    stop("'", split.by, "' not found in seurat_object@meta.data.")

  # ── Extract expression matrix ─────────────────────────────────────────────
  mat <- .get_layer_data(seurat_object, assay = assay, layer = layer)

  # ── Build cell → group label vector ───────────────────────────────────────
  if (is.null(group.by)) {
    idents <- SeuratObject::Idents(seurat_object)
    y <- stats::setNames(as.character(idents), names(idents))
  } else {
    if (!group.by %in% colnames(meta))
      stop("'", group.by, "' not found in seurat_object@meta.data.")
    y <- stats::setNames(as.character(meta[[group.by]]), rownames(meta))
  }

  # ── Align cells ───────────────────────────────────────────────────────────
  common <- intersect(colnames(mat), names(y))
  if (length(common) == 0L)
    stop("No cells overlap between expression matrix columns and group labels. ",
         "Check that 'assay' and 'layer' point to the right slot.")
  if (length(common) < ncol(mat))
    message("scSidekick: ", ncol(mat) - length(common),
            " cell(s) in the matrix have no group label and will be excluded.")

  mat <- mat[, common, drop = FALSE]
  y   <- y[common]

  # ── Split path ────────────────────────────────────────────────────────────
  if (!is.null(split.by)) {
    split_vals   <- stats::setNames(as.character(meta[[split.by]]), rownames(meta))
    split_vals   <- split_vals[common]
    split_levels <- unique(split_vals[!is.na(split_vals)])

    results <- lapply(split_levels, function(lvl) {
      keep_s  <- !is.na(split_vals) & split_vals == lvl
      mat_s   <- mat[, keep_s, drop = FALSE]
      y_s     <- y[keep_s]

      if (!is.null(groups_use)) {
        keep_g <- y_s %in% groups_use
        mat_s  <- mat_s[, keep_g, drop = FALSE]
        y_s    <- y_s[keep_g]
      }

      n_grp <- length(unique(y_s))
      if (ncol(mat_s) == 0L || n_grp < 2L) {
        message("scSidekick: Skipping ", split.by, " = '", lvl,
                "' (", ncol(mat_s), " cell(s), ", n_grp, " group(s) — need ≥ 2).")
        return(NULL)
      }

      message("scSidekick: [", lvl, "] wilcoxauc — ",
              nrow(mat_s), " features × ", ncol(mat_s), " cells, ",
              n_grp, " groups",
              if (!is.null(groups_use)) paste0(" [", paste(groups_use, collapse = ", "), "]") else "",
              ".")

      res <- presto::wilcoxauc(X = mat_s, y = y_s, ...)
      res[[split.by]] <- lvl
      res
    })

    out <- do.call(rbind, Filter(Negate(is.null), results))
    # Move split column to front for readability
    out <- out[, c(split.by, setdiff(colnames(out), split.by)), drop = FALSE]
    return(out)
  }

  # ── No-split path ─────────────────────────────────────────────────────────
  if (!is.null(groups_use)) {
    keep <- y %in% groups_use
    if (!any(keep))
      stop("None of 'groups_use' values found in '",
           if (is.null(group.by)) "Idents()" else group.by, "'.\n",
           "Available: ", paste(unique(y), collapse = ", "))
    mat <- mat[, keep, drop = FALSE]
    y   <- y[keep]
  }

  message("scSidekick: Running wilcoxauc — ",
          nrow(mat), " features × ", ncol(mat), " cells, ",
          length(unique(y)), " groups",
          if (!is.null(groups_use)) paste0(" [", paste(groups_use, collapse = ", "), "]") else "",
          ".")

  presto::wilcoxauc(X = mat, y = y, ...)
}
