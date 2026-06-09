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
#' @param seurat_object A Seurat object.
#' @param group.by Character. Metadata column whose levels define the
#'   comparison groups (e.g. \code{"Cluster"}, \code{"Assignment"}).
#'   When \code{NULL}, \code{\link[Seurat]{Idents}} is used.
#' @param groups_use Character vector or \code{NULL}. Restrict the test to
#'   a subset of \code{group.by} levels.  \code{NULL} = all levels.
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
#'
#' @seealso \code{\link[presto]{wilcoxauc}}
#'
#' @examples
#' \dontrun{
#' # Basic usage - use active identity
#' markers <- RunWilcoxAUC(SeuratObj)
#'
#' # Specify a metadata column
#' markers <- RunWilcoxAUC(SeuratObj, group.by = "Cluster")
#'
#' # Only compare specific groups
#' markers <- RunWilcoxAUC(SeuratObj, group.by = "Cluster",
#'                          groups_use = c("Microglia", "Astrocytes"))
#'
#' # Top markers per cluster
#' library(dplyr)
#' top_markers <- markers |>
#'   filter(padj < 0.05, logFC > 0.5, auc > 0.6) |>
#'   group.by(group) |>
#'   slice_max(auc, n = 20)
#' }
#'
#' @export
RunWilcoxAUC <- function(
    seurat_object,
    group.by   = NULL,
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

  # ── Extract expression matrix ─────────────────────────────────────────────
  # .get_layer_data handles BPCells lazy subset, v3 slot= vs v5 layer=,
  # and coerces the result to a dgCMatrix.
  mat <- .get_layer_data(seurat_object, assay = assay, layer = layer)

  # ── Build cell → group label vector ───────────────────────────────────────
  if (is.null(group.by)) {
    idents <- SeuratObject::Idents(seurat_object)
    y <- stats::setNames(as.character(idents), names(idents))
  } else {
    meta <- seurat_object@meta.data
    if (!group.by %in% colnames(meta))
      stop("'", group.by, "' not found in seurat_object@meta.data.")
    y <- stats::setNames(as.character(meta[[group.by]]), rownames(meta))
  }

  # ── Align cells ───────────────────────────────────────────────────────────
  # Matrix columns and label names must agree.  Subset to the intersection
  # (handles the case where the assay has fewer cells than the full metadata).
  common <- intersect(colnames(mat), names(y))
  if (length(common) == 0L)
    stop("No cells overlap between expression matrix columns and group labels. ",
         "Check that 'assay' and 'layer' point to the right slot.")
  if (length(common) < ncol(mat))
    message("scSidekick: ", ncol(mat) - length(common),
            " cell(s) in the matrix have no group label and will be excluded.")

  mat <- mat[, common, drop = FALSE]
  y   <- y[common]

  # ── Filter to requested groups ────────────────────────────────────────────
  if (!is.null(groups_use)) {
    keep <- y %in% groups_use
    if (!any(keep))
      stop("None of 'groups_use' values found in '",
           if (is.null(group.by)) "Idents()" else group.by, "'.\n",
           "Available: ", paste(unique(y), collapse = ", "))
    mat <- mat[, keep, drop = FALSE]
    y   <- y[keep]
  }

  message("scSidekick: Running wilcoxauc - ",
          nrow(mat), " features × ", ncol(mat), " cells, ",
          length(unique(y)), " groups",
          if (!is.null(groups_use)) paste0(" [", paste(groups_use, collapse = ", "), "]") else "",
          ".")

  # ── Call presto ───────────────────────────────────────────────────────────
  # Passing a dgCMatrix dispatches to wilcoxauc.default (or wilcoxauc.dgCMatrix)
  # via S3 - bypasses the wilcoxauc.Seurat path entirely.  No ::: needed.
  presto::wilcoxauc(X = mat, y = y, ...)
}
