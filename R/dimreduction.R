# =============================================================================
# Dimensionality reduction helpers
# =============================================================================

#' Determine the optimal number of dimensions for downstream analysis
#'
#' Uses two complementary heuristics to choose how many dimensions from a
#' linear reduction (PCA, LSI, etc.) to pass to UMAP, neighbor graphs, and
#' clustering.
#'
#' **Heuristic 1 — cumulative variance plateau**: finds the first dimension
#' where cumulative variance explained exceeds `threshold_cum_pct` *and* the
#' per-dimension contribution has already dropped below `threshold_pct`. This
#' marks the point of diminishing returns.
#'
#' **Heuristic 2 — elbow**: finds the last dimension where the drop in
#' per-dimension variance between successive dimensions is still larger than
#' `threshold_diff_pct`. Beyond this point the curve has flattened.
#'
#' The function returns the **minimum** of the two, which is the more
#' conservative choice.
#'
#' @param seurat_object A Seurat object with the target reduction already
#'   computed (e.g. via `RunPCA()` or `RunLSI()`).
#' @param threshold_pct Numeric. A selected dimension must contribute less than
#'   this percentage of total variance. Default `5`.
#' @param reduction Character. Name of the reduction slot to inspect. Default
#'   `"pca"`.
#' @param threshold_cum_pct Numeric. Cumulative variance (%) that must be
#'   reached before Heuristic 1 can trigger. Default `90`.
#' @param threshold_diff_pct Numeric. Minimum drop in per-dimension variance
#'   between consecutive dimensions for Heuristic 2. Default `0.1`.
#' @param elbow_plot Logical. Return an annotated elbow plot? Default `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{n_dims}{Recommended number of dimensions (integer).}
#'   \item{cum_thresh}{Dimension selected by Heuristic 1.}
#'   \item{elbow_thresh}{Dimension selected by Heuristic 2.}
#'   \item{plot}{A ggplot2 elbow plot, or `NULL` when `elbow_plot = FALSE`.}
#' }
#'
#' @export
Determine_nDims <- function(seurat_object,
                             threshold_pct      = 5,
                             reduction          = "pca",
                             threshold_cum_pct  = 90,
                             threshold_diff_pct = 0.1,
                             elbow_plot         = TRUE) {

  stdev   <- seurat_object[[reduction]]@stdev
  pct_var <- stdev / sum(stdev) * 100
  cum_var <- cumsum(pct_var)

  # Heuristic 1 — cumulative plateau:
  # first dim where cumulative variance > threshold_cum_pct
  # AND per-dim contribution has already fallen below threshold_pct
  h1 <- which(cum_var > threshold_cum_pct & pct_var < threshold_pct)[1L]

  # Heuristic 2 — elbow:
  # last dim where the step-down in variance between consecutive dims
  # is still larger than threshold_diff_pct (i.e. the curve is still dropping)
  h2 <- sort(which(abs(diff(pct_var)) > threshold_diff_pct),
             decreasing = TRUE)[1L] + 1L

  n_dims <- min(h1, h2)

  message(sprintf("Heuristic 1 (cumulative plateau) : dim %d", h1))
  message(sprintf("Heuristic 2 (elbow)              : dim %d", h2))
  message(sprintf("Selected (conservative minimum)  : dim %d", n_dims))

  plot_obj <- NULL
  if (elbow_plot) {
    plot_df <- data.frame(
      pct_var = pct_var,
      cum_var = cum_var,
      dim     = seq_along(pct_var)
    )

    plot_obj <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(cum_var, pct_var, label = dim, color = dim <= n_dims)
    ) +
      ggplot2::geom_text(size = 4, fontface = "bold") +
      ggplot2::labs(
        x = "Cumulative variance (%)",
        y = "Variance per dimension (%)"
      ) +
      ggplot2::geom_vline(
        xintercept = threshold_cum_pct,
        linetype = "dashed", color = "black", alpha = 0.5
      ) +
      ggplot2::annotate(
        "text",
        x = threshold_cum_pct, y = max(pct_var),
        label = sprintf("%d%% threshold → dim %d", threshold_cum_pct, h1),
        angle = 90, hjust = 1, vjust = -0.3, size = 3.5, color = "black"
      ) +
      ggplot2::geom_vline(
        xintercept = plot_df$cum_var[n_dims],
        linetype = "dotdash", color = "black", alpha = 0.75
      ) +
      ggplot2::annotate(
        "text",
        x = plot_df$cum_var[n_dims], y = max(pct_var),
        label = sprintf("Selected: %d dims", n_dims),
        angle = 90, hjust = 1, vjust = -0.3, size = 3.5, color = "black"
      ) +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "dodgerblue", "FALSE" = "darkorange"),
        guide  = "none"
      ) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(axis.text = ggplot2::element_text(color = "black"))
  }

  list(n_dims       = n_dims,
       cum_thresh   = h1,
       elbow_thresh = h2,
       plot         = plot_obj)
}
