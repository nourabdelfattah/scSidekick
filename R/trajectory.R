# =============================================================================
# scSidekick pseudotime / trajectory wrappers  (trajectory.R)
#
# RunSlingshot()       - fits Slingshot directly on a Seurat reduction (any
#                        reduction, including Harmony / other batch-corrected
#                        embeddings) + a cluster metadata column. Never
#                        converts to SingleCellExperiment - the expression
#                        matrix is never touched, so there is no risk of
#                        re-normalizing already-normalized data, and nothing
#                        needs to be stripped from the object first (no assay,
#                        no version, no BPCells/sketch compatibility issue).
#
# SlingshotLineages()  - human-readable summary of what RunSlingshot() found:
#                        one row per lineage, with its cluster path and cell
#                        count, so picking a lineage means reading a path
#                        instead of guessing from a number.
#
# PlotPseudotime()     - engine-agnostic pseudotime distribution panel: takes
#                        any numeric meta.data column (a Slingshot lineage, a
#                        Monocle3 pseudotime, or anything else), so the same
#                        function serves every trajectory method instead of
#                        duplicating density/ridge plotting code per engine.
#
# PlotTrajectory()     - ggplot2-only trajectory overlay on any DimPlot.
#                        style = "straight" draws the cluster-level MST as
#                        straight segments between cluster centroids (the
#                        clean, DDRTree-like look that is easier to read than
#                        Slingshot's snake curves or Monocle3's dense graph in
#                        complex settings) - built from slingMST(), so it is
#                        guaranteed to match the fitted lineages rather than
#                        being a separately hand-built minimum spanning tree.
#                        style = "curves" draws Slingshot's own smoothed
#                        principal curves instead, when the plot reduction is
#                        the same one Slingshot was fit on.
#
# Both follow scSidekick conventions: results are stored on the object
# (@meta.data for per-cell pseudotime, @misc$slingshot for the fitted model)
# so they can be recalled anytime without recomputation; missing packages
# produce an actionable stop() rather than a raw error; and a run always
# writes analysis_params.json (via .write_subdir_params()) when output_dir
# is available, so the method is transparently documented alongside the
# figures that use it.
# =============================================================================


# Resolve which cells to fit on for a fast-preview run.
# Prefers the object's own "sketch" assay cells (an already-principled
# representative subsample computed by SketchData()) over an arbitrary
# downsample, and never calls subset() on the full Seurat object - only the
# already-extracted embedding matrix and cluster vector are subset, so this
# stays fast regardless of object size.
.slingshot_preview_cells <- function(seurat_object, use_sketch, downsample, cluster_vec) {
  all_cells <- colnames(seurat_object)

  if (isTRUE(use_sketch) && "sketch" %in% SeuratObject::Assays(seurat_object)) {
    sketch_cells <- tryCatch(SeuratObject::Cells(seurat_object[["sketch"]]),
                             error = function(e) NULL)
    if (!is.null(sketch_cells) && length(sketch_cells) > 0 &&
        length(sketch_cells) < length(all_cells)) {
      message("scSidekick: 'sketch' assay found - using its ",
              format(length(sketch_cells), big.mark = ","),
              " representative cells for a fast preview run. Set ",
              "use_sketch = FALSE to run on all ",
              format(length(all_cells), big.mark = ","), " cells.")
      return(sketch_cells)
    }
  }

  if (!is.null(downsample)) {
    by_cluster <- split(all_cells, cluster_vec)
    sampled <- unlist(lapply(by_cluster, function(cells) {
      if (length(cells) > downsample) sample(cells, downsample) else cells
    }), use.names = FALSE)
    message("scSidekick: downsampling to ", downsample,
            " cells per cluster for a fast preview run (",
            format(length(sampled), big.mark = ","), " of ",
            format(length(all_cells), big.mark = ","), " total). Set ",
            "downsample = NULL to run on all cells.")
    return(sampled)
  }

  all_cells
}


#' Fit Slingshot pseudotime directly on a Seurat reduction
#'
#' @description
#' Runs \code{\link[slingshot]{slingshot}} directly on a Seurat dimensional
#' reduction (PCA, Harmony, or any other integrated embedding) and a cluster
#' identity column - no \code{SingleCellExperiment} conversion involved.
#' Slingshot's own algorithm only ever needs a reduced-dimension matrix and
#' cluster labels (it builds a cluster-level minimum spanning tree, then fits
#' a principal curve per lineage), so bypassing the object conversion avoids
#' the two biggest practical problems with Seurat-wrapper packages for
#' trajectory inference:
#' \itemize{
#'   \item the expression matrix is never touched, so there is no
#'     renormalize-already-normalized-data risk downstream, and
#'   \item there is no assay/version stripping to do first - it does not
#'     matter whether the object has BPCells, sketch, or multiple v5 assays,
#'     since only the embedding matrix (always dense) and a metadata column
#'     are read.
#' }
#'
#' Results are written back onto \code{seurat_object} so they can be reused
#' without recomputation: per-cell pseudotime as \code{meta.data} columns
#' (\code{<prefix>_Pseudotime_1}, \code{_2}, ...) and the fitted model, cluster
#' MST, and lineage paths in \code{seurat_object@misc$slingshot}. Use
#' \code{\link{SlingshotLineages}} afterward to see which lineage number
#' corresponds to which cluster path.
#'
#' @param seurat_object A Seurat object with at least one dimensional
#'   reduction computed.
#' @param cluster.by Character. Metadata column of cluster/cell-type
#'   identities to build the lineage structure on.
#' @param reduction Character or \code{NULL}. Reduction to fit on. \code{NULL}
#'   (default) auto-detects: prefers \code{"harmony"} or any other
#'   integration-style reduction (name matching "integrated", "scvi", "mnn",
#'   "scanorama") over a plain \code{"pca"}, so batch-corrected embeddings are
#'   used by default when present.
#' @param dims Integer vector or \code{NULL}. Which components of
#'   \code{reduction} to use. \code{NULL} (default) uses the first
#'   \code{max_dims} available - see \code{max_dims}.
#' @param max_dims Integer. When \code{dims} is \code{NULL}, the number of
#'   leading components to use from \code{reduction}. Default \code{20L}.
#'   Principal-curve fitting cost scales with dimensionality (roughly 6x
#'   slower going from 2 to 50 dims at matched cell count in local testing),
#'   so a Harmony/PCA reduction with 30-50 computed components is capped
#'   rather than used in full. Ignored if \code{dims} is supplied explicitly,
#'   or if \code{reduction} has fewer than \code{max_dims} components (e.g. a
#'   2D UMAP).
#' @param start.clus Character or \code{NULL}. Root cluster (level of
#'   \code{cluster.by}) that lineages should originate from. Strongly
#'   recommended - without it, Slingshot's MST is unrooted and lineage
#'   direction is arbitrary.
#' @param end.clus Character vector or \code{NULL}. Optional terminal
#'   cluster(s) to constrain lineage endpoints.
#' @param approx_points Integer or \code{NULL}. Passed to
#'   \code{slingshot::slingshot()}; reduces the number of points used to fit
#'   principal curves, which speeds up large runs at a small cost to
#'   smoothness. \code{NULL} uses Slingshot's default (no approximation).
#' @param use_sketch Logical. If the object has a \code{"sketch"} assay
#'   (from \code{SketchData()}), fit on those cells only, as a fast preview.
#'   Default \code{TRUE}. Has no effect if no \code{"sketch"} assay exists.
#' @param downsample Integer or \code{NULL}. If supplied (and no usable
#'   \code{"sketch"} assay was found), cap each \code{cluster.by} level at
#'   this many cells before fitting, for a fast preview on large objects.
#'   Cells not included in the fit get \code{NA} pseudotime. Default
#'   \code{NULL} (use all cells).
#' @param prefix Character. Prefix for the pseudotime \code{meta.data}
#'   columns. Default \code{"Slingshot"}.
#' @param progress Logical. Show live "still fitting" heartbeat messages
#'   while Slingshot runs, so a long fit is never silent. Requires the
#'   \code{callr} package (runs the fit in a background process so R can
#'   poll and print while waiting); this also keeps the wait itself
#'   interruptible via Ctrl+C at any point, unlike being blocked directly
#'   inside Slingshot's own fitting code. Falls back to a plain blocking call
#'   with an upfront message only if \code{callr} is not installed. Default
#'   \code{TRUE}.
#' @param poll_interval Numeric. Seconds between heartbeat messages when
#'   \code{progress = TRUE}. Default \code{10}.
#' @param caffeinate Logical. Prevent the machine from sleeping for the
#'   duration of the fit (macOS/Windows only; see \code{\link{RunCellChat}}).
#'   Default \code{FALSE}.
#' @param output_dir Character or \code{NULL}. Directory to write
#'   \code{analysis_params.json} documenting this run. \code{NULL} (default)
#'   walks up to the \code{output_dir} stored by \code{\link{PrepObject}}
#'   (respecting \code{AutoSavePlots}); pass \code{NA} to suppress entirely.
#'
#' @return \code{seurat_object} with pseudotime columns added to
#'   \code{meta.data} and the fitted result in \code{@misc$slingshot}.
#' @export
RunSlingshot <- function(seurat_object,
                         cluster.by,
                         reduction     = NULL,
                         dims          = NULL,
                         max_dims      = 20L,
                         start.clus    = NULL,
                         end.clus      = NULL,
                         approx_points = NULL,
                         use_sketch    = TRUE,
                         downsample    = NULL,
                         prefix        = "Slingshot",
                         progress      = TRUE,
                         poll_interval = 10,
                         caffeinate    = FALSE,
                         output_dir    = NULL) {

  if (!requireNamespace("slingshot", quietly = TRUE))
    stop("Package 'slingshot' is required. Install with:\n",
         "  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
         "  BiocManager::install('slingshot')")

  if (caffeinate) { .caff <- .nk_caffeinate(); on.exit(.nk_decaffeinate(.caff), add = TRUE) }

  if (missing(cluster.by) || is.null(cluster.by))
    stop("'cluster.by' is required - name the metadata column of cluster/",
         "cell-type identities to build the trajectory on.")
  if (!cluster.by %in% colnames(seurat_object@meta.data))
    stop("'", cluster.by, "' not found in seurat_object@meta.data.")

  if (identical(output_dir, NA_character_) || identical(output_dir, NA)) {
    output_dir <- NULL
  } else if (missing(output_dir) || is.null(output_dir)) {
    output_dir <- if (.nk_autosave(seurat_object))
      .nk_setting(seurat_object, "output_dir") else NULL
  }

  # ── Resolve embedding (harmony-aware) and full cluster vector ──────────────
  emb <- .get_embedding(seurat_object, reduction = reduction, dims = dims)
  reduction_used <- attr(emb, "reduction")

  # Principal-curve fitting cost scales with dimensionality, not just cell
  # count (~6x slower going from 2 to 50 dims at matched n in local testing).
  # Harmony/PCA reductions commonly carry 30-50 components; cap to max_dims
  # by default unless the caller explicitly asked for specific dims - this is
  # also standard practice (most trajectory analyses use ~10-20 PCs, not all
  # computed components).
  if (is.null(dims) && ncol(emb) > max_dims) {
    message("scSidekick: '", reduction_used, "' has ", ncol(emb),
            " components; using the first ", max_dims, " for Slingshot ",
            "(principal-curve fitting cost scales with dimensionality). ",
            "Pass dims = 1:", ncol(emb), " to use them all, or a different ",
            "max_dims / dims to change the cutoff.")
    emb <- emb[, seq_len(max_dims), drop = FALSE]
  }
  dims_used <- dims %||% seq_len(ncol(emb))

  cl_full <- as.character(seurat_object@meta.data[[cluster.by]])
  names(cl_full) <- colnames(seurat_object)

  if (!is.null(start.clus) && !start.clus %in% cl_full)
    stop("start.clus = '", start.clus, "' is not a level of '", cluster.by,
         "'. Available: ", paste(sort(unique(cl_full)), collapse = ", "))
  if (!is.null(end.clus) && !all(end.clus %in% cl_full))
    stop("end.clus contains levels not present in '", cluster.by,
         "'. Available: ", paste(sort(unique(cl_full)), collapse = ", "))

  # ── Fast-preview cell selection (never subsets the Seurat object itself) ──
  use_cells <- .slingshot_preview_cells(seurat_object, use_sketch, downsample, cl_full)
  fast_preview <- length(use_cells) < ncol(seurat_object)

  emb_sub <- emb[use_cells, , drop = FALSE]
  cl_sub  <- cl_full[use_cells]

  # ── Fit Slingshot on the embedding + cluster labels only ───────────────────
  message("scSidekick: Fitting Slingshot on ", format(length(use_cells), big.mark = ","),
          " cells x ", ncol(emb_sub), " dims",
          if (ncol(emb_sub) > 10)
            " - higher dimensionality means this can take a while." else ".")

  # Fitting is a single opaque, mostly-uninterruptible call inside slingshot's
  # compiled code - R cannot print anything else while blocked inside it. The
  # only way to get live "still running" messages (and to stay interruptible
  # via Ctrl+C in the meantime) is to run the fit in a separate process and
  # poll it from here, where a plain Sys.sleep() loop always responds to
  # interrupts immediately.
  use_bg <- isTRUE(progress) && requireNamespace("callr", quietly = TRUE)
  if (isTRUE(progress) && !use_bg)
    message("scSidekick: Install 'callr' (install.packages('callr')) for live progress ",
            "heartbeats and a safely-interruptible fit; proceeding without them.")

  if (use_bg) {
    rb <- callr::r_bg(
      func = function(data, clusterLabels, start.clus, end.clus, approx_points) {
        args <- list(data = data, clusterLabels = clusterLabels)
        if (!is.null(start.clus))    args$start.clus    <- start.clus
        if (!is.null(end.clus))      args$end.clus      <- end.clus
        if (!is.null(approx_points)) args$approx_points <- approx_points
        do.call(slingshot::slingshot, args)
      },
      args = list(data = emb_sub, clusterLabels = cl_sub, start.clus = start.clus,
                 end.clus = end.clus, approx_points = approx_points)
    )
    # Kills the background fit if the user interrupts (or an error unwinds)
    # while we're polling, rather than leaving it running as an orphan.
    on.exit(if (rb$is_alive()) rb$kill(), add = TRUE)

    t0 <- Sys.time()
    while (rb$is_alive()) {
      Sys.sleep(poll_interval)
      if (rb$is_alive())
        message("scSidekick: ...still fitting (",
                round(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
                "s elapsed). Safe to interrupt (Ctrl+C/Esc) - the background ",
                "fit will be cleaned up.")
    }
    res <- rb$get_result()

  } else {
    sling_args <- list(data = emb_sub, clusterLabels = cl_sub)
    if (!is.null(start.clus))    sling_args$start.clus    <- start.clus
    if (!is.null(end.clus))      sling_args$end.clus      <- end.clus
    if (!is.null(approx_points)) sling_args$approx_points <- approx_points
    res <- do.call(slingshot::slingshot, sling_args)
  }

  pt       <- slingshot::slingPseudotime(res)
  lineages <- slingshot::slingLineages(res)
  mst_df   <- tryCatch(slingshot::slingMST(res, as.df = TRUE), error = function(e) NULL)
  n_lineages <- ncol(pt)
  n_cells_per_lineage <- stats::setNames(colSums(!is.na(pt)), colnames(pt))

  # ── Write flat pseudotime columns back to the FULL object (NA where unfit) ─
  for (i in seq_len(n_lineages)) {
    col_name <- paste0(prefix, "_Pseudotime_", i)
    full_col <- rep(NA_real_, ncol(seurat_object))
    names(full_col) <- colnames(seurat_object)
    full_col[rownames(pt)] <- pt[, i]
    seurat_object@meta.data[[col_name]] <- full_col[colnames(seurat_object)]
  }

  seurat_object@misc$slingshot <- list(
    result              = res,
    lineages            = lineages,
    mst_df              = mst_df,
    n_cells_per_lineage = n_cells_per_lineage,
    cluster.by          = cluster.by,
    reduction           = reduction_used,
    dims                = dims_used,
    start.clus          = start.clus,
    end.clus            = end.clus,
    prefix              = prefix,
    cells_used          = use_cells,
    fast_preview        = fast_preview,
    date                = format(Sys.Date())
  )

  # ── Console summary: cluster path per lineage, so numbers mean something ──
  message("scSidekick: Slingshot found ", n_lineages, " lineage",
          if (n_lineages != 1) "s" else "", ":")
  for (i in seq_len(n_lineages)) {
    ln <- paste0("Lineage", i)
    message("  ", ln, ": ", paste(lineages[[ln]], collapse = " -> "),
            "  (", format(n_cells_per_lineage[[ln]], big.mark = ","), " cells)")
  }
  message("scSidekick: Pseudotime stored as meta.data column",
          if (n_lineages != 1) "s" else "", " '", prefix, "_Pseudotime_1'",
          if (n_lineages != 1)
            paste0("..'", prefix, "_Pseudotime_", n_lineages, "'") else "",
          ". Model stored in seurat_object@misc$slingshot. Use ",
          "SlingshotLineages() to re-print this summary.")

  # ── Methods JSON (transparency: every method run documents itself) ────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    unit <- .nk_unit_label(seurat_object)
    methods_text <- paste0(
      "Cell lineages and pseudotime were inferred using Slingshot (Street et ",
      "al., 2018) on the ", reduction_used, " embedding (dimensions ",
      min(dims_used), "-", max(dims_used), ")",
      if (fast_preview) paste0(
        ", fit on a representative subsample of ",
        format(length(use_cells), big.mark = ","), " of ",
        format(ncol(seurat_object), big.mark = ","), " ", unit,
        " for a fast preview") else "",
      ". A cluster-level minimum spanning tree was built on '", cluster.by,
      "' identities",
      if (!is.null(start.clus)) paste0(", rooted at '", start.clus, "'") else "",
      ", yielding ", n_lineages, " lineage", if (n_lineages != 1) "s" else "",
      "; per-cell pseudotime was estimated along each lineage's principal curve."
    )
    # .write_subdir_params() only ever merges with a PARENT json (one found
    # by walking up from dirname(output_dir)) - it never reads a json already
    # sitting inside output_dir itself. When output_dir is the same folder
    # another method already wrote to, that local file has to be read and
    # merged here first, or it gets overwritten with only Slingshot's fields.
    local_json <- file.path(output_dir, "analysis_params.json")
    prev <- if (requireNamespace("jsonlite", quietly = TRUE) && file.exists(local_json))
      tryCatch(jsonlite::read_json(local_json, simplifyVector = TRUE),
               error = function(e) list())
    else list()

    .write_subdir_params(output_dir, utils::modifyList(prev, list(
      date                    = format(Sys.Date()),
      slingshot_reduction     = reduction_used,
      slingshot_dims          = paste0(min(dims_used), "-", max(dims_used)),
      slingshot_cluster_by    = cluster.by,
      slingshot_start_clus    = start.clus,
      slingshot_n_lineages    = n_lineages,
      slingshot_n_cells_used  = length(use_cells),
      slingshot_fast_preview  = fast_preview,
      slingshot_methods_text  = methods_text
    )))
  }

  seurat_object
}


#' Summarize Slingshot lineages found by RunSlingshot()
#'
#' @description
#' Returns one row per lineage discovered by \code{\link{RunSlingshot}}: its
#' cluster path (root to terminal cluster) and the number of cells with
#' non-\code{NA} pseudotime along it. Meant to answer "which lineage number is
#' actually the one I care about" without re-deriving it from the model object.
#'
#' @param seurat_object A Seurat object previously processed by
#'   \code{\link{RunSlingshot}}.
#'
#' @return A data frame with columns \code{Lineage}, \code{Path}, and
#'   \code{n_cells}.
#' @export
SlingshotLineages <- function(seurat_object) {
  sl <- seurat_object@misc$slingshot
  if (is.null(sl))
    stop("No Slingshot results found on this object. Run RunSlingshot() first.")

  ln_names <- names(sl$lineages)
  data.frame(
    Lineage   = ln_names,
    Path      = vapply(sl$lineages, paste, character(1), collapse = " -> "),
    n_cells   = as.integer(sl$n_cells_per_lineage[ln_names]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}


# Build cluster-centroid straight-line segments from slingMST(as.df=TRUE),
# reusing the fitted MST's cluster connectivity (which clusters connect, and
# in what order per lineage) but positioning nodes at this PLOT reduction's
# cluster centroids - so the topology always matches what was fit, even when
# plotting on a different reduction (e.g. fit on harmony, plot on UMAP).
.slingshot_mst_segments <- function(mst_df, centers, lineages_keep) {
  mst_df$Lineage <- as.character(mst_df$Lineage)
  if (!is.null(lineages_keep)) mst_df <- mst_df[mst_df$Lineage %in% lineages_keep, ]

  segs <- lapply(split(mst_df, mst_df$Lineage), function(df) {
    df <- df[order(df$Order), ]
    if (nrow(df) < 2) return(NULL)
    from <- centers$cluster_id[match(df$Cluster[-nrow(df)], centers$cluster)]
    to   <- centers$cluster_id[match(df$Cluster[-1],        centers$cluster)]
    data.frame(
      x       = centers$x[from], y    = centers$y[from],
      xend    = centers$x[to],   yend = centers$y[to],
      Lineage = unique(df$Lineage)
    )
  })
  do.call(rbind, segs)
}


#' Overlay a Slingshot trajectory on a DimPlot
#'
#' @description
#' Draws a standard cluster DimPlot with the Slingshot trajectory overlaid,
#' in one of two styles:
#' \itemize{
#'   \item \code{style = "straight"} (default) - straight segments between
#'     cluster centroids, following the cluster-level minimum spanning tree
#'     Slingshot fit (\code{slingshot::slingMST()}). This is the clean,
#'     DDRTree-like look that stays legible with many clusters, and works on
#'     \emph{any} plot reduction since only the MST's cluster connectivity is
#'     reused - the node positions come from this plot's own reduction.
#'   \item \code{style = "curves"} - Slingshot's smoothed per-lineage
#'     principal curves (\code{slingshot::slingCurves()}). These live in
#'     whatever reduction Slingshot was \emph{fit} on, so this style requires
#'     \code{reduction} to match \code{seurat_object@misc$slingshot$reduction}
#'     (e.g. fit and plot both on \code{"umap"}); otherwise the curve
#'     coordinates would not correspond to the plotted axes, and this
#'     function stops rather than draw a misleading overlay.
#' }
#'
#' @param seurat_object A Seurat object previously processed by
#'   \code{\link{RunSlingshot}}.
#' @param reduction Character. Reduction to plot on. Default \code{"umap"}.
#' @param group.by Character or \code{NULL}. Metadata column for point color.
#'   \code{NULL} (default) uses the \code{cluster.by} stored by
#'   \code{\link{RunSlingshot}}.
#' @param lineages Integer/character vector or \code{NULL}. Which lineages to
#'   draw (e.g. \code{c(1, 2)} or \code{c("Lineage1", "Lineage2")}).
#'   \code{NULL} (default) draws all.
#' @param style \code{"straight"} or \code{"curves"}. See Description.
#' @param colors Named color vector. \code{NULL} auto-resolves from
#'   \code{\link{PrepObject}} or \code{Nour_pal}.
#' @param pt.size Numeric. Point size. Default \code{0.4}.
#' @param line.color Character. Trajectory line color. Default \code{"black"}.
#' @param line.width Numeric. Trajectory line width. Default \code{0.8}.
#' @param label_lineages Logical. When more than one lineage is drawn, map
#'   line type (solid/dashed/dotted/...) to lineage and add a "Lineage"
#'   legend, so overlapping trajectory lines can be told apart. Color is left
#'   alone (still \code{line.color}) so this doesn't collide with the cluster
#'   color legend. Default \code{TRUE}; has no visible effect (and no legend)
#'   when only one lineage is drawn.
#' @param arrow Logical. Draw arrowheads on straight segments (ignored for
#'   \code{style = "curves"}, which is unordered per point). Default
#'   \code{TRUE}.
#' @param label Logical. Label cluster centroids. Default \code{TRUE}.
#' @param legend.side Character. Legend position, matching
#'   \code{\link{PlotDimPlots}}: \code{"bottom"} (default), \code{"top"},
#'   \code{"left"}, or \code{"right"}.
#' @param legend_nrow Integer. Rows for the cluster legend when
#'   \code{legend.side} is \code{"bottom"}/\code{"top"}. Default \code{2},
#'   matching \code{PlotDimPlots}'s default - keeping this the same makes the
#'   two functions' output comparable side by side.
#' @param output_dir Character or \code{NULL}. Directory to save a PDF.
#'   \code{NULL} (default) walks up to the \code{output_dir} stored by
#'   \code{\link{PrepObject}} (respecting \code{AutoSavePlots}); pass
#'   \code{NA} to suppress entirely.
#' @param object_name,subset_name Character. Output file-name prefixes.
#' @param file_name Character or \code{NULL}. Base name (no extension) for
#'   the saved PDF. \code{NULL} auto-deduces.
#' @param width,height Numeric or \code{NULL}. Saved PDF dimensions in
#'   inches. \code{NULL} (default) uses a fixed 6 x 5.5 in panel.
#'
#' @return A \code{ggplot} object (invisibly, if saved).
#' @export
PlotTrajectory <- function(seurat_object,
                           reduction      = "umap",
                           group.by       = NULL,
                           lineages       = NULL,
                           style          = c("straight", "curves"),
                           colors         = NULL,
                           pt.size        = 0.4,
                           line.color     = "black",
                           line.width     = 0.8,
                           label_lineages = TRUE,
                           arrow          = TRUE,
                           label          = TRUE,
                           legend.side    = "bottom",
                           legend_nrow    = 2,
                           output_dir     = NULL,
                           object_name    = "",
                           subset_name    = "",
                           file_name      = NULL,
                           width          = NULL,
                           height         = NULL) {

  style <- match.arg(style)

  sl <- seurat_object@misc$slingshot
  if (is.null(sl))
    stop("No Slingshot results found on this object. Run RunSlingshot() first.")

  group.by <- group.by %||% sl$cluster.by

  if (identical(output_dir, NA_character_) || identical(output_dir, NA)) {
    output_dir <- NULL
  } else if (missing(output_dir) || is.null(output_dir)) {
    output_dir <- if (.nk_autosave(seurat_object))
      .nk_setting(seurat_object, "output_dir") else NULL
  }
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  # Normalize requested lineage names to "Lineage<n>" strings
  lineages_keep <- if (is.null(lineages)) NULL else {
    ifelse(grepl("^Lineage", lineages), lineages, paste0("Lineage", lineages))
  }

  emb <- .resolve_emb(seurat_object, reduction)

  # Respect factor levels already set by PrepObject() (or the user) rather
  # than flattening to character, which would silently drop back to
  # alphabetical ordering - same convention as PlotDimPlots/PlotComposition.
  group_col <- seurat_object@meta.data[[group.by]]
  group_lvls <- if (is.factor(group_col)) levels(group_col) else
    sort(unique(as.character(group_col)))
  plot_df <- data.frame(x = emb[, 1], y = emb[, 2],
                        cluster = factor(as.character(group_col), levels = group_lvls))

  centers <- stats::aggregate(cbind(x, y) ~ cluster, data = plot_df, FUN = mean)
  centers$cluster_id <- seq_len(nrow(centers))

  colors <- colors %||% .nk_colors(seurat_object, group.by)

  # Same visual grammar as PlotDimPlots: theme_classic(), no numeric tick
  # labels (UMAP/PCA coordinates are arbitrary), bordered panel with
  # axis.line blanked to avoid the double-line clipping that comes from
  # having both a panel border and axis lines, and axis titles taken
  # verbatim from the embedding's own column names (not a re-cased guess),
  # so output sits comfortably next to PlotDimPlots in the same figure.
  lgd_pos <- if (legend.side %in% c("bottom", "top")) "bottom" else "right"
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = cluster), size = pt.size) +
    ggplot2::labs(x = colnames(emb)[1], y = colnames(emb)[2],
                 color = group.by, title = "Slingshot trajectory") +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes   = list(size = 5),
      nrow           = legend_nrow,
      title.position = "top",
      title          = group.by
    )) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
      panel.border     = ggplot2::element_rect(color = "darkgray", fill = NA),
      axis.line        = ggplot2::element_blank(),
      legend.position  = lgd_pos,
      legend.direction = "horizontal"
    )
  if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors)

  if (style == "straight") {
    if (is.null(sl$mst_df))
      stop("No cluster MST stored for this object (slingMST() failed at fit ",
           "time). Re-run RunSlingshot().")
    segs <- .slingshot_mst_segments(sl$mst_df, centers, lineages_keep)
    if (!is.null(segs) && nrow(segs) > 0) {
      # Map line TYPE (not color) to lineage when there's more than one, so
      # overlapping trajectory lines can be told apart without a second color
      # scale colliding with the cluster color legend.
      if (isTRUE(label_lineages) && length(unique(segs$Lineage)) > 1) {
        p <- p + ggplot2::geom_segment(
          data = segs,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linetype = Lineage),
          inherit.aes = FALSE, color = line.color, linewidth = line.width,
          arrow = if (isTRUE(arrow))
            grid::arrow(length = grid::unit(0.12, "inches")) else NULL
        ) + ggplot2::guides(linetype = ggplot2::guide_legend(title = "Lineage"))
      } else {
        p <- p + ggplot2::geom_segment(
          data = segs,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          inherit.aes = FALSE, color = line.color, linewidth = line.width,
          arrow = if (isTRUE(arrow))
            grid::arrow(length = grid::unit(0.12, "inches")) else NULL
        )
      }
    }
  } else {
    fit_reduction <- sl$reduction
    hit <- SeuratObject::Reductions(seurat_object)[
      tolower(SeuratObject::Reductions(seurat_object)) == tolower(reduction)]
    if (length(hit) == 0 || !identical(tolower(hit[1]), tolower(fit_reduction)))
      stop("style = 'curves' requires reduction = '", fit_reduction,
           "' (the reduction Slingshot was fit on), not '", reduction,
           "'. Use style = 'straight' to overlay the trajectory on a ",
           "different reduction, or re-run RunSlingshot(reduction = '",
           reduction, "').")

    crv <- slingshot::slingCurves(sl$result)
    if (!is.null(lineages_keep))
      crv <- crv[paste0("Lineage", seq_along(crv)) %in% lineages_keep]

    curve_df <- do.call(rbind, lapply(seq_along(crv), function(i) {
      s <- crv[[i]]$s[crv[[i]]$ord, , drop = FALSE]
      data.frame(x = s[, 1], y = s[, 2], Lineage = paste0("Lineage", i))
    }))
    if (isTRUE(label_lineages) && length(unique(curve_df$Lineage)) > 1) {
      p <- p + ggplot2::geom_path(
        data = curve_df,
        ggplot2::aes(x = x, y = y, group = Lineage, linetype = Lineage),
        inherit.aes = FALSE, color = line.color, linewidth = line.width
      ) + ggplot2::guides(linetype = ggplot2::guide_legend(title = "Lineage"))
    } else {
      p <- p + ggplot2::geom_path(
        data = curve_df,
        ggplot2::aes(x = x, y = y, group = Lineage),
        inherit.aes = FALSE, color = line.color, linewidth = line.width
      )
    }
  }

  if (isTRUE(label)) {
    p <- p + ggrepel::geom_label_repel(
      data = centers,
      ggplot2::aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE, fontface = "bold", show.legend = FALSE
    )
  }

  # ── Save or return ──────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    base <- if (!is.null(file_name) && nzchar(file_name)) file_name else {
      paste(c(if (nchar(object_name) > 0) object_name,
              if (nchar(subset_name) > 0) subset_name,
              "SlingshotTrajectory", style), collapse = "_")
    }
    fname <- gsub("[^A-Za-z0-9._-]", "_", base)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))
    n_ln <- length(unique(if (is.null(lineages_keep)) names(sl$lineages) else lineages_keep))

    # Auto-size must budget space for BOTH legends (cluster color + Lineage
    # linetype, when shown) - previously a flat 6x5.5 assumed no legend at
    # all, so with enough cluster levels the legend (and even a cluster or
    # two) ran off the saved PDF's edge.
    cluster_lgd <- .cat_legend_dims(length(group_lvls), group_lvls, "horizontal", legend_nrow)
    lineage_lgd_w <- 0; lineage_lgd_h <- 0
    if (isTRUE(label_lineages) && n_ln > 1) {
      ln_dims <- .cat_legend_dims(n_ln, paste0("Lineage", seq_len(n_ln)), "horizontal", 1)
      lineage_lgd_w <- ln_dims[["w"]]; lineage_lgd_h <- ln_dims[["h"]]
    }
    panel_w <- 6; panel_h <- 5.5
    if (legend.side %in% c("bottom", "top")) {
      pdf_w <- width  %||% max(panel_w, cluster_lgd[["w"]] + lineage_lgd_w)
      pdf_h <- height %||% (panel_h + max(cluster_lgd[["h"]], lineage_lgd_h))
    } else {
      pdf_w <- width  %||% (panel_w + max(cluster_lgd[["w"]], lineage_lgd_w))
      pdf_h <- height %||% max(panel_h, cluster_lgd[["h"]] + lineage_lgd_h)
    }

    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(p)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath,
            " (", round(pdf_w, 1), " x ", round(pdf_h, 1), " in)")

    ctx <- .nk_legend_context(seurat_object)
    .write_legend_sidecar(fpath, paste0(
      .nk_obs_sentence(ctx, reduction = reduction, group.by = group.by),
      ". Slingshot trajectory overlaid (", style, " style) for ", n_ln,
      " lineage", if (n_ln != 1) "s" else "", " rooted at '",
      sl$start.clus %||% "unspecified", "', fit on the ", sl$reduction,
      " embedding. Cluster centroid labels shown.",
      if (isTRUE(label_lineages) && n_ln > 1)
        " Line type distinguishes lineages (see legend)." else ""
    ))
    return(invisible(p))
  }

  p
}


#' Plot the distribution of pseudotime across clusters
#'
#' @description
#' Engine-agnostic pseudotime panel: takes any numeric \code{meta.data}
#' column, so the same function plots a Slingshot lineage
#' (\code{"Slingshot_Pseudotime_1"}), a Monocle3 pseudotime
#' (\code{"Monocle3_Pseudotime"}), or any other pseudotime-like column,
#' without re-implementing the same density/ridge plotting code per method.
#'
#' Produces a density + rug plot of pseudotime per \code{group.by} level, with
#' peak labels, and optionally a ridge plot of the same values split/filled by
#' a second variable (e.g. treatment response).
#'
#' @param seurat_object A Seurat object.
#' @param pseudotime.by Character. Numeric \code{meta.data} column holding a
#'   pseudotime value (e.g. from \code{\link{RunSlingshot}}).
#' @param group.by Character or \code{NULL}. Cluster/cell-type column for
#'   color and the y-axis of the ridge plot. \code{NULL} (default) uses the
#'   \code{cluster.by} stored by \code{\link{RunSlingshot}} when present.
#' @param split.by Character or \code{NULL}. Secondary variable (e.g.
#'   \code{"Response"}, \code{"Sample"}) used to fill the ridge plot.
#'   \code{NULL} (default) skips the ridge panel.
#' @param colors Named color vector for \code{group.by}. \code{NULL}
#'   auto-resolves from \code{\link{PrepObject}} or \code{Nour_pal}.
#' @param split_colors Named color vector for \code{split.by}. \code{NULL}
#'   auto-resolves the same way.
#' @param output_dir Character or \code{NULL}. Directory to save a PDF.
#'   \code{NULL} (default) walks up to the \code{output_dir} stored by
#'   \code{\link{PrepObject}} (respecting \code{AutoSavePlots}); pass
#'   \code{NA} to suppress entirely.
#' @param object_name,subset_name Character. Output file-name prefixes.
#' @param file_name Character or \code{NULL}. Base name (no extension) for
#'   the saved PDF. \code{NULL} auto-deduces from \code{pseudotime.by}.
#' @param width,height Numeric or \code{NULL}. Saved PDF dimensions in
#'   inches. \code{NULL} (default) auto-sizes to 1 or 2 panels.
#'
#' @return A \code{patchwork} (or single \code{ggplot}) object.
#' @export
PlotPseudotime <- function(seurat_object,
                           pseudotime.by,
                           group.by     = NULL,
                           split.by     = NULL,
                           colors       = NULL,
                           split_colors = NULL,
                           output_dir   = NULL,
                           object_name  = "",
                           subset_name  = "",
                           file_name    = NULL,
                           width        = NULL,
                           height       = NULL) {

  if (!pseudotime.by %in% colnames(seurat_object@meta.data))
    stop("'", pseudotime.by, "' not found in seurat_object@meta.data.")

  group.by <- group.by %||% seurat_object@misc$slingshot$cluster.by %||%
    stop("'group.by' is required (no stored Slingshot cluster.by found). ",
         "Provide the cluster/cell-type metadata column directly.")

  if (identical(output_dir, NA_character_) || identical(output_dir, NA)) {
    output_dir <- NULL
  } else if (missing(output_dir) || is.null(output_dir)) {
    output_dir <- if (.nk_autosave(seurat_object))
      .nk_setting(seurat_object, "output_dir") else NULL
  }
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  meta <- seurat_object@meta.data

  # Respect factor levels already set by PrepObject() (or the user) rather
  # than flattening to character, which would silently drop back to
  # alphabetical ordering - same convention as PlotDimPlots/PlotComposition.
  group_col  <- meta[[group.by]]
  group_lvls <- if (is.factor(group_col)) levels(group_col) else
    sort(unique(as.character(group_col)))

  df <- data.frame(pseudo  = meta[[pseudotime.by]],
                   cluster = factor(as.character(group_col), levels = group_lvls))
  df <- df[!is.na(df$pseudo), ]

  colors <- colors %||% .nk_colors(seurat_object, group.by)

  # ── Peak labels: one per cluster with enough points for a density estimate ─
  label_data <- do.call(rbind, lapply(split(df, df$cluster), function(sub) {
    if (nrow(sub) < 2) return(NULL)
    d <- stats::density(sub$pseudo)
    data.frame(cluster = sub$cluster[1],
              x_peak  = d$x[which.max(d$y)],
              y_peak  = max(d$y))
  }))

  # theme_NourMin() is the package's canonical theme for statistical/
  # distribution plots (see coexpression.R, pseudobulk.R, meta_summary.R); it
  # blanks axis titles by default, so both are restored here since the
  # pseudotime column name and "Density"/group.by labels are the point.
  p_density <- ggplot2::ggplot(df, ggplot2::aes(x = pseudo, fill = cluster)) +
    ggplot2::geom_rug() +
    ggplot2::geom_density(alpha = 0.9) +
    ggplot2::labs(x = pseudotime.by, y = "Density", fill = group.by) +
    theme_NourMin() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 11),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5, size = 11),
      # theme_NourMin()'s default plot.margin has a NEGATIVE left value
      # (tuned for its default axis.title = element_blank()); since the axis
      # title is restored above, the margin must be too, or the rotated
      # y-axis title clips off the left edge of the page.
      plot.margin  = ggplot2::margin(t = 5, r = 5, b = 5, l = 10, unit = "mm")
    )
  if (!is.null(colors)) p_density <- p_density + ggplot2::scale_fill_manual(values = colors)
  if (!is.null(label_data) && nrow(label_data) > 0) {
    p_density <- p_density + ggrepel::geom_label_repel(
      data = label_data,
      ggplot2::aes(x = x_peak, y = y_peak, label = cluster, fill = cluster),
      inherit.aes = FALSE, color = "white", fontface = "bold",
      nudge_y = 0.05, segment.color = "black", show.legend = FALSE
    )
  }

  p_ridge <- NULL
  if (!is.null(split.by)) {
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      warning("Package 'ggridges' is not installed - skipping the ridge panel. ",
              "Install with install.packages('ggridges') to enable it.")
    } else {
      split_col  <- meta[[split.by]]
      split_lvls <- if (is.factor(split_col)) levels(split_col) else
        sort(unique(as.character(split_col)))

      df2 <- data.frame(pseudo  = meta[[pseudotime.by]],
                        cluster = factor(as.character(group_col), levels = group_lvls),
                        split   = factor(as.character(split_col), levels = split_lvls))
      df2 <- df2[!is.na(df2$pseudo), ]
      split_colors <- split_colors %||% .nk_colors(seurat_object, split.by)

      p_ridge <- ggplot2::ggplot(df2, ggplot2::aes(x = pseudo, y = cluster, fill = split)) +
        ggridges::geom_density_ridges2(alpha = 0.7) +
        ggplot2::labs(x = pseudotime.by, y = group.by, fill = split.by) +
        theme_NourMin() +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 11),
          axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5, size = 11),
          plot.margin  = ggplot2::margin(t = 5, r = 5, b = 5, l = 10, unit = "mm")
        )
      if (!is.null(split_colors))
        p_ridge <- p_ridge + ggplot2::scale_fill_manual(values = split_colors)
    }
  }

  result <- if (!is.null(p_ridge))
    patchwork::wrap_plots(list(p_density, p_ridge), ncol = 2)
  else p_density

  # ── Save or return ──────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    base <- if (!is.null(file_name) && nzchar(file_name)) file_name else {
      paste(c(if (nchar(object_name) > 0) object_name,
              if (nchar(subset_name) > 0) subset_name,
              pseudotime.by), collapse = "_")
    }
    fname <- gsub("[^A-Za-z0-9._-]", "_", base)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))
    pdf_w <- width %||% (if (!is.null(p_ridge)) 11 else 5.5)
    pdf_h <- height %||% 5

    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(result)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath,
            " (", round(pdf_w, 1), " x ", round(pdf_h, 1), " in)")

    ctx <- .nk_legend_context(seurat_object)
    .write_legend_sidecar(fpath, paste0(
      .nk_obs_sentence(ctx, group.by = group.by,
                      split.by = if (!is.null(p_ridge)) split.by else NULL),
      ". Distribution of '", pseudotime.by, "' shown as a density + rug plot",
      if (!is.null(p_ridge))
        paste0(", alongside a ridge plot of the same values by ", group.by,
              " filled by ", split.by) else "",
      "; peak labels mark each ", group.by, " level's modal pseudotime."
    ))
    return(invisible(result))
  }

  result
}


#' Plot gene/feature expression trends across pseudotime
#'
#' @description
#' Smoothed expression (or any numeric feature) vs. pseudotime, one facet
#' panel per feature, optionally colored by a grouping variable (e.g.
#' \code{"Sample"}) to compare trends across replicates or conditions.
#' Equivalent to a \code{geom_smooth()} + \code{facet_wrap(~Gene)} script run
#' by hand, but feature extraction goes through the same assay/version-safe
#' path as the rest of scSidekick (\code{.nk_feature_matrix()}, which is what
#' every other feature-plotting function in the package uses): genes are
#' pulled via the v5/BPCells/sketch-safe \code{.get_layer_data()}, and numeric
#' \code{meta.data} columns (module scores, cNMF usages, QC metrics, other
#' pseudotime columns, ...) are used directly - both are accepted
#' interchangeably in \code{features}.
#'
#' @param seurat_object A Seurat object, typically after
#'   \code{\link{RunSlingshot}}.
#' @param pseudotime.by Character. Numeric \code{meta.data} column for the
#'   x-axis (e.g. \code{"Slingshot_Pseudotime_1"}).
#' @param features Character vector. Gene names and/or numeric
#'   \code{meta.data} columns. Each gets its own facet panel.
#' @param assay Character or \code{NULL}. Assay to pull gene expression from.
#'   \code{NULL} (default) uses \code{DefaultAssay(seurat_object)}.
#' @param layer Character. Assay layer/slot. Default \code{"data"}
#'   (log-normalized counts).
#' @param min_expr Numeric or \code{NULL}. Expression floor applied before
#'   smoothing: observations with expression \eqn{\le} \code{min_expr} are
#'   dropped so the trend reflects the expressing population rather than being
#'   pulled toward zero by dropout. \code{min_expr = 0} removes exact zeros.
#'   \code{NULL} (default) keeps all cells.
#' @param donor.by Character or \code{NULL}. When a donor/patient column is
#'   given, cells are collapsed to \strong{one point per donor} (mean x and
#'   mean y) before fitting, so both the trend and points reflect independent
#'   patients rather than cells. This is the honest way to relate a donor-level
#'   variable (e.g. age, a clinical score) to a per-cell value, and avoids
#'   drawing thousands of overlapping cell points. \code{show_points} defaults
#'   to \code{TRUE} when this is set. \code{NULL} (default) keeps cell-level.
#' @param color.by Character or \code{NULL}. Metadata column to color the
#'   smoothed trend lines by (e.g. \code{"Sample"}, to see whether replicates
#'   agree). \code{NULL} (default) draws a single trend line per feature.
#' @param group.by Character or \code{NULL}. Metadata column to restrict
#'   cells by (e.g. plot the trend within one cluster only).
#' @param groups Character vector or \code{NULL}. Which level(s) of
#'   \code{group.by} to keep. \code{NULL} (default, when \code{group.by} is
#'   set) keeps all levels.
#' @param method Character. Smoothing method passed to
#'   \code{ggplot2::geom_smooth()}. \code{"loess"} (default) fits a flexible,
#'   non-parametric curve - appropriate when the relationship's shape is
#'   unknown or non-monotonic. Use \code{"lm"} for a straight-line (linear)
#'   fit when the relationship looks linear/monotonic (e.g. after inspecting
#'   it with \code{add_cor = TRUE}, or from prior knowledge of the variables);
#'   \code{"glm"} is also accepted for a generalized linear fit.
#' @param se Logical. Show the smoother's standard-error ribbon. Default
#'   \code{FALSE}.
#' @param add_cor Logical. Annotate each facet with a correlation coefficient
#'   and p-value (via \code{ggpubr::stat_cor()}), computed directly from the
#'   plotted x/y (donor-level when \code{donor.by} is set) - this is the
#'   direct way to test the relationship shown by the trend line without
#'   routing through \code{\link{CorrelateGene}()}/\code{\link{PlotCorrelation}()},
#'   which are built around a gene-vs-many-genes screen keyed on a Seurat
#'   object rather than an arbitrary x/y pair. When \code{color.by} is set,
#'   one annotation is shown per color group. Default \code{FALSE}.
#' @param cor_method Character. \code{"spearman"} (default; rank-based,
#'   appropriate for monotonic but not necessarily linear relationships) or
#'   \code{"pearson"} (linear correlation). Only used when \code{add_cor =
#'   TRUE}.
#' @param show_points Logical. Overlay the raw observations under the trend.
#'   Besides showing individual points, this trains the y-axis on the full data
#'   range - \code{geom_smooth()} alone auto-scales the axis to the fitted line,
#'   which can look "clipped" relative to the spread of the data. Default
#'   \code{FALSE}.
#' @param point_size,point_alpha Size and transparency of the overlaid points
#'   when \code{show_points = TRUE}. Defaults \code{0.4} and \code{0.25}.
#' @param scales Character. Passed to \code{facet_wrap()}. Default
#'   \code{"free_y"} (expression ranges differ across genes).
#' @param ncol Integer or \code{NULL}. Facet columns. \code{NULL} auto-picks.
#' @param colors Named color vector for \code{color.by}. \code{NULL}
#'   auto-resolves from \code{\link{PrepObject}} or \code{Nour_pal}.
#' @param output_dir Character or \code{NULL}. Directory to save a PDF.
#'   \code{NULL} (default) walks up to the \code{output_dir} stored by
#'   \code{\link{PrepObject}} (respecting \code{AutoSavePlots}); pass
#'   \code{NA} to suppress entirely.
#' @param object_name,subset_name Character. Output file-name prefixes.
#' @param file_name Character or \code{NULL}. Base name (no extension) for
#'   the saved PDF. \code{NULL} auto-deduces.
#' @param width,height Numeric or \code{NULL}. Saved PDF dimensions in
#'   inches. \code{NULL} (default) auto-sizes from the number of facets.
#'
#' @return A \code{ggplot} object.
#' @export
PlotFeatureTrend <- function(seurat_object,
                             pseudotime.by,
                             features,
                             assay       = NULL,
                             layer       = "data",
                             min_expr    = NULL,
                             donor.by    = NULL,
                             color.by    = NULL,
                             group.by    = NULL,
                             groups      = NULL,
                             method      = "loess",
                             se          = FALSE,
                             add_cor     = FALSE,
                             cor_method  = c("spearman", "pearson"),
                             show_points = FALSE,
                             point_size  = 0.4,
                             point_alpha = 0.25,
                             scales      = "free_y",
                             ncol        = NULL,
                             colors      = NULL,
                             output_dir  = NULL,
                             object_name = "",
                             subset_name = "",
                             file_name   = NULL,
                             width       = NULL,
                             height      = NULL) {

  cor_method <- match.arg(cor_method)

  if (!pseudotime.by %in% colnames(seurat_object@meta.data))
    stop("'", pseudotime.by, "' not found in seurat_object@meta.data.")

  # Patient-level aggregation yields few points, so show them by default (unless
  # the caller explicitly set show_points).
  if (!is.null(donor.by) && missing(show_points)) show_points <- TRUE

  if (isTRUE(add_cor) && !requireNamespace("ggpubr", quietly = TRUE))
    stop("Package 'ggpubr' is required for add_cor = TRUE.")

  if (identical(output_dir, NA_character_) || identical(output_dir, NA)) {
    output_dir <- NULL
  } else if (missing(output_dir) || is.null(output_dir)) {
    output_dir <- if (.nk_autosave(seurat_object))
      .nk_setting(seurat_object, "output_dir") else NULL
  }
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  fm <- .nk_feature_matrix(seurat_object, features, assay = assay, layer = layer)
  if (is.null(fm$mat))
    stop("None of the requested features were found as genes in assay '",
         assay %||% SeuratObject::DefaultAssay(seurat_object),
         "' or as numeric meta.data columns: ", paste(features, collapse = ", "))
  if (length(fm$missing) > 0)
    warning("Feature(s) not found, skipping: ", paste(fm$missing, collapse = ", "))

  meta <- seurat_object@meta.data
  wide <- data.frame(pseudo = meta[[pseudotime.by]], t(fm$mat),
                     check.names = FALSE, stringsAsFactors = FALSE)

  if (!is.null(color.by)) {
    color_col  <- meta[[color.by]]
    color_lvls <- if (is.factor(color_col)) levels(color_col) else
      sort(unique(as.character(color_col)))
    wide$.color <- factor(as.character(color_col), levels = color_lvls)
  }
  if (!is.null(group.by)) {
    wide$.group <- as.character(meta[[group.by]])
    keep <- if (!is.null(groups)) wide$.group %in% groups else TRUE
    wide <- wide[keep, ]
  }
  if (!is.null(donor.by)) {
    if (!donor.by %in% colnames(meta))
      stop("'donor.by' column '", donor.by, "' not found in meta.data.")
    wide$.donor <- as.character(meta[rownames(wide), donor.by])
  }
  wide <- wide[!is.na(wide$pseudo), ]

  long <- tidyr::pivot_longer(wide, cols = fm$found,
                              names_to = "Feature", values_to = "Expression")
  long$Feature <- factor(long$Feature, levels = fm$found)

  # Drop non-expressing (dropout) observations before fitting the trend so the
  # smoother reflects the expressing population rather than being pulled toward
  # zero by the dropout pile-up. `min_expr = 0` removes exact zeros; a higher
  # value applies a stricter expression floor. Applied per feature (long form).
  if (!is.null(min_expr)) {
    n_before <- nrow(long)
    long     <- long[!is.na(long$Expression) & long$Expression > min_expr, ]
    if (nrow(long) == 0)
      stop("No observations remain after applying min_expr = ", min_expr,
           "; lower the threshold.")
    message("scSidekick: min_expr = ", min_expr, " kept ", nrow(long), " of ",
            n_before, " (feature, cell) observations for the trend.")
  }

  # Patient-level aggregation: collapse cells to one point per donor (mean x and
  # mean y) so the trend and points reflect independent patients rather than
  # cells - removing the cells-per-patient weighting that confounds a cell-level
  # trend of a donor-level variable. Both the smoother and the (optional) points
  # are then drawn at donor level.
  if (!is.null(donor.by)) {
    grp_cols <- c(".donor", "Feature",
                  if (!is.null(color.by)) ".color",
                  if (!is.null(group.by)) ".group")
    long <- long |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) |>
      dplyr::summarize(pseudo     = mean(pseudo,     na.rm = TRUE),
                       Expression = mean(Expression, na.rm = TRUE),
                       .groups = "drop")
    long$Feature <- factor(as.character(long$Feature), levels = fm$found)
    message("scSidekick: aggregated to ",
            length(unique(long$.donor)), " donors (", donor.by, ") per feature.")
  }

  aes_args <- list(x = quote(pseudo), y = quote(Expression))
  if (!is.null(color.by)) aes_args$color <- quote(.color)

  p <- ggplot2::ggplot(long, do.call(ggplot2::aes, aes_args))
  # Optional raw points UNDER the trend. Besides showing individual
  # observations, this trains the y-scale on the actual data range so the axis
  # spans the full spread of values rather than only the fitted-line range
  # (geom_smooth alone auto-scales to the trend, which looks "clipped").
  if (isTRUE(show_points))
    p <- p + ggplot2::geom_point(size = point_size, alpha = point_alpha)
  p <- p +
    ggplot2::geom_smooth(method = method, se = se) +
    ggplot2::facet_wrap(~Feature, scales = scales, ncol = ncol)
  # Correlation annotation (rho/r + p), computed directly from the plotted x/y
  # (pseudotime vs. feature - donor-level when donor.by is set) rather than
  # through CorrelateGene()/PlotCorrelation(), which are built around gene x
  # many-genes screens keyed on a Seurat object, not an arbitrary x/y pair.
  # ggpubr::stat_cor() reads the same mapping as geom_smooth(), so it respects
  # facet_wrap(~Feature) and color.by automatically (one annotation per color
  # group, stacked, when color.by is set).
  if (isTRUE(add_cor))
    p <- p + ggpubr::stat_cor(method = cor_method, show.legend = FALSE)
  p <- p +
    ggplot2::labs(x = pseudotime.by, y = "Expression", color = color.by) +
    theme_NourMin() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 11),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5, size = 11),
      strip.text   = ggplot2::element_text(face = "bold"),
      plot.margin  = ggplot2::margin(t = 5, r = 5, b = 5, l = 10, unit = "mm")
    )

  if (!is.null(color.by)) {
    colors <- colors %||% .nk_colors(seurat_object, color.by)
    if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors)
  }

  # ── Save or return ──────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    base <- if (!is.null(file_name) && nzchar(file_name)) file_name else {
      paste(c(if (nchar(object_name) > 0) object_name,
              if (nchar(subset_name) > 0) subset_name,
              pseudotime.by, "FeatureTrend"), collapse = "_")
    }
    fname <- gsub("[^A-Za-z0-9._-]", "_", base)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))
    n_feat <- length(fm$found)
    facet_ncol <- ncol %||% min(n_feat, 3)
    facet_nrow <- ceiling(n_feat / facet_ncol)

    # Auto-size must budget space for the color.by legend (default right-hand
    # position) - it was previously left out entirely, so long labels/many
    # levels ran off the edge of the saved PDF instead of just being cramped.
    lgd_w <- 0; lgd_h <- 0
    if (!is.null(color.by)) {
      lgd_dims <- .cat_legend_dims(length(levels(wide$.color)), levels(wide$.color),
                                   "vertical", 1)
      lgd_w <- lgd_dims[["w"]]; lgd_h <- lgd_dims[["h"]]
    }
    pdf_w <- width  %||% (facet_ncol * 3.5 + lgd_w)
    pdf_h <- height %||% max(facet_nrow * 3.2, lgd_h)

    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(p)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath,
            " (", round(pdf_w, 1), " x ", round(pdf_h, 1), " in)")

    ctx        <- .nk_legend_context(seurat_object)
    method_lbl <- switch(method,
      lm    = "Linear (lm)",
      loess = "Loess-smoothed",
      glm   = "GLM-fitted",
      paste0(method, "-fitted"))
    unit_lbl   <- if (!is.null(donor.by)) "donors" else ctx$unit

    .write_legend_sidecar(fpath, paste0(
      .nk_obs_sentence(ctx, group.by = color.by),
      ". ", method_lbl, " trend of ", n_feat, " feature",
      if (n_feat != 1) "s" else "", " (", paste(fm$found, collapse = ", "),
      ") vs. '", pseudotime.by, "'",
      if (!is.null(group.by)) paste0(
        ", restricted to '", group.by, "' = ",
        paste(groups %||% "all levels", collapse = ", ")) else "",
      if (!is.null(color.by)) paste0(", colored by ", color.by) else "",
      if (!is.null(donor.by)) paste0(
        ". Cells were aggregated to one point per ", donor.by,
        " (mean of both axes) before fitting, so the trend reflects ",
        length(unique(long$.donor)), " independent ", unit_lbl,
        " rather than being weighted by cells per ", donor.by)
      else "",
      if (isTRUE(show_points)) paste0(
        ". Individual ", unit_lbl, " are shown as points beneath the trend line")
      else "",
      if (!is.null(min_expr)) paste0(
        ". Observations with expression <= ", min_expr,
        " were excluded before smoothing to avoid the dropout pile-up at zero")
      else "",
      if (isTRUE(add_cor)) paste0(
        ". ", if (cor_method == "spearman") "Spearman's rho" else "Pearson's r",
        " and its p-value are annotated on each facet",
        if (!is.null(color.by)) paste0(" (one per ", color.by, " level)") else "",
        ", computed directly from the plotted values")
      else "", "."
    ))
    return(invisible(p))
  }

  p
}


#' Density or ridge plot for any gene or numeric metadata feature
#'
#' @description
#' A more general version of Seurat's \code{RidgePlot()} - one density
#' (or ridge) per \code{group.by} level, for any single feature - but with
#' the two things \code{RidgePlot()} does not support: a \code{split.by} fill
#' \emph{within} each ridge, and a \code{facet.by} panel split. That combination
#' answers questions like "one ridge per Sample, one panel per Response" to
#' check whether replicate samples show the same distribution shape within
#' each condition.
#'
#' \code{feature} is resolved the same assay/version-safe way as every other
#' feature-plotting function in the package (\code{.nk_feature_matrix()}):
#' gene names go through \code{.get_layer_data()} (v5/BPCells/sketch-safe),
#' numeric \code{meta.data} columns (a pseudotime, a module score, a QC
#' metric, ...) are used directly.
#'
#' @param seurat_object A Seurat object.
#' @param feature Character vector. One or more gene names and/or numeric
#'   \code{meta.data} columns. When several are supplied, one PDF is written
#'   per feature (the file name always includes the feature).
#' @param group.by Character. Categorical \code{meta.data} column - one
#'   density/ridge per level.
#' @param split.by Character or \code{NULL}. Categorical \code{meta.data}
#'   column used as a secondary fill \emph{within} each \code{group.by}
#'   ridge. Only used when \code{style = "ridge"}.
#' @param facet.by Character or \code{NULL}. Categorical \code{meta.data}
#'   column to split the whole plot into panels via \code{facet_wrap()}
#'   (e.g. \code{group.by = "Sample"}, \code{facet.by = "Response"} to check
#'   whether replicates overlap within each response group).
#' @param style \code{"ridge"} (default, via \code{ggridges}) or
#'   \code{"density"} (overlaid \code{geom_density()} curves, fill =
#'   \code{group.by}; \code{split.by} is ignored for this style).
#' @param assay Character or \code{NULL}. Assay to pull gene expression from.
#'   \code{NULL} (default) uses \code{DefaultAssay(seurat_object)}.
#' @param layer Character. Assay layer/slot. Default \code{"data"}.
#' @param min_expr Numeric or \code{NULL}. Expression floor: cells with
#'   expression \eqn{\le} \code{min_expr} are excluded so the distribution
#'   reflects the expressing population rather than a spike at zero from
#'   dropout. \code{min_expr = 0} removes exact zeros. \code{NULL} (default)
#'   keeps all cells.
#' @param colors Named color vector for the fill variable
#'   (\code{split.by} if set, else \code{group.by}). \code{NULL}
#'   auto-resolves from \code{\link{PrepObject}} or \code{Nour_pal}.
#' @param scales Character. Passed to \code{facet_wrap()} when
#'   \code{facet.by} is set. Default \code{"fixed"} (keeps axes directly
#'   comparable across panels - usually what you want when checking whether
#'   distributions overlap).
#' @param ncol Integer or \code{NULL}. Facet columns. \code{NULL} auto-picks.
#' @param output_dir Character or \code{NULL}. Directory to save a PDF.
#'   \code{NULL} (default) walks up to the \code{output_dir} stored by
#'   \code{\link{PrepObject}} (respecting \code{AutoSavePlots}); pass
#'   \code{NA} to suppress entirely.
#' @param object_name,subset_name Character. Output file-name prefixes.
#' @param file_name Character or \code{NULL}. Base name (no extension) for
#'   the saved PDF. \code{NULL} auto-deduces a name that includes the feature,
#'   \code{group.by} and \code{split.by}. When \code{file_name} is given and
#'   several features are requested, the feature is appended so each still
#'   gets its own file.
#' @param width,height Numeric or \code{NULL}. Saved PDF dimensions in
#'   inches. \code{NULL} (default) auto-sizes.
#'
#' @return A \code{ggplot} object for a single feature, or (for several
#'   features) a named list of \code{ggplot} objects, returned invisibly when
#'   saved to disk.
#' @export
PlotFeatureDistribution <- function(seurat_object,
                                    feature,
                                    group.by,
                                    split.by    = NULL,
                                    facet.by    = NULL,
                                    style       = c("ridge", "density"),
                                    assay       = NULL,
                                    layer       = "data",
                                    min_expr    = NULL,
                                    colors      = NULL,
                                    scales      = "fixed",
                                    ncol        = NULL,
                                    output_dir  = NULL,
                                    object_name = "",
                                    subset_name = "",
                                    file_name   = NULL,
                                    width       = NULL,
                                    height      = NULL) {

  style <- match.arg(style)
  meta  <- seurat_object@meta.data

  if (is.numeric(meta[[group.by]]))
    stop("'group.by' must be a categorical variable. '", group.by, "' is numeric.")

  if (identical(output_dir, NA_character_) || identical(output_dir, NA)) {
    output_dir <- NULL
  } else if (missing(output_dir) || is.null(output_dir)) {
    output_dir <- if (.nk_autosave(seurat_object))
      .nk_setting(seurat_object, "output_dir") else NULL
  }
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  # ── Shared setup (independent of which feature is plotted) ──────────────────
  group_col  <- meta[[group.by]]
  group_lvls <- if (is.factor(group_col)) levels(group_col) else
    sort(unique(as.character(group_col)))

  use_split_fill <- !is.null(split.by) && style == "ridge"
  if (!is.null(split.by) && style == "density")
    warning("split.by is ignored for style = 'density' (would overlay too ",
            "many colors); use style = 'ridge', or pass this variable as ",
            "facet.by instead.")

  split_lvls <- if (use_split_fill) {
    sc <- meta[[split.by]]
    if (is.factor(sc)) levels(sc) else sort(unique(as.character(sc)))
  } else NULL
  facet_lvls <- if (!is.null(facet.by)) {
    fc <- meta[[facet.by]]
    if (is.factor(fc)) levels(fc) else sort(unique(as.character(fc)))
  } else NULL

  fill_var   <- if (use_split_fill) "split" else "group"
  fill_label <- if (use_split_fill) split.by else group.by
  colors <- colors %||%
    .nk_colors(seurat_object, if (use_split_fill) split.by else group.by)

  if (style == "ridge" && !requireNamespace("ggridges", quietly = TRUE))
    stop("Package 'ggridges' is required for style = 'ridge'. Install with ",
         "install.packages('ggridges'), or use style = 'density'.")

  features <- feature          # `feature` may name one OR several genes/columns
  multi    <- length(features) > 1L

  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE,
                                       showWarnings = FALSE)

  # ── Build one distribution plot for a single feature ────────────────────────
  .build_one <- function(feat) {
    if (feat %in% colnames(meta) && !is.numeric(meta[[feat]]))
      stop("'", feat, "' is a categorical metadata column.\n",
           "PlotFeatureDistribution() is for numeric variables and gene ",
           "expression. For categorical variables use PlotDimPlots(), ",
           "PlotComposition(), or SplitDotPlot().")

    vals <- .nk_feature_values(seurat_object, feat, assay = assay, layer = layer)
    df <- data.frame(value = as.numeric(vals[rownames(meta)]),
                     group = factor(as.character(group_col), levels = group_lvls))
    if (use_split_fill)
      df$split <- factor(as.character(meta[[split.by]]), levels = split_lvls)
    if (!is.null(facet.by))
      df$facet <- factor(as.character(meta[[facet.by]]), levels = facet_lvls)
    df <- df[!is.na(df$value), ]

    # Drop non-expressing (dropout) cells so the distribution reflects the
    # expressing population rather than a spike at zero. `min_expr = 0` removes
    # exact zeros; a higher value applies a stricter expression floor.
    if (!is.null(min_expr)) {
      df <- df[df$value > min_expr, ]
      if (nrow(df) == 0)
        stop("'", feat, "': no cells remain after applying min_expr = ",
             min_expr, "; lower the threshold.")
    }

    p <- if (style == "ridge") {
      ggplot2::ggplot(df, ggplot2::aes(x = value, y = group,
                                       fill = .data[[fill_var]])) +
        ggridges::geom_density_ridges2(alpha = 0.7)
    } else {
      ggplot2::ggplot(df, ggplot2::aes(x = value, fill = .data[[fill_var]])) +
        ggplot2::geom_density(alpha = 0.7)
    }
    p <- p +
      ggplot2::labs(x = feat, y = if (style == "ridge") group.by else "Density",
                   fill = fill_label) +
      theme_NourMin() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 11),
        axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5, size = 11),
        strip.text   = ggplot2::element_text(face = "bold"),
        plot.margin  = ggplot2::margin(t = 5, r = 5, b = 5, l = 10, unit = "mm")
      )
    if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
    if (!is.null(facet.by))
      p <- p + ggplot2::facet_wrap(~facet, scales = scales, ncol = ncol)
    list(plot = p, df = df)
  }

  # ── Save one feature's plot to its own PDF ──────────────────────────────────
  .save_one <- function(feat, res) {
    p <- res$plot; df <- res$df
    # File name includes group.by and split.by so different groupings/splits
    # never overwrite each other. With an explicit file_name and several
    # features, the gene is appended so each still gets its own file.
    base <- if (!is.null(file_name) && nzchar(file_name)) {
      if (multi) paste(file_name, feat, sep = "_") else file_name
    } else {
      paste(c(if (nchar(object_name) > 0) object_name,
              if (nchar(subset_name) > 0) subset_name,
              feat, group.by,
              if (!is.null(split.by)) split.by,
              style),
            collapse = "_")
    }
    fname <- gsub("[^A-Za-z0-9._-]", "_", base)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))
    n_facets <- if (!is.null(facet.by)) length(levels(df$facet)) else 1

    # Auto-size must budget space for: (1) the fill legend (default right-hand
    # position), and (2) for ridge style, the group.by category labels on the
    # y-axis, which can be wide with long sample/cluster names.
    fill_lvls <- levels(df[[fill_var]])
    lgd_dims  <- .cat_legend_dims(length(fill_lvls), fill_lvls, "vertical", 1)
    lgd_w     <- lgd_dims[["w"]]; lgd_h <- lgd_dims[["h"]]
    y_axis_w  <- if (style == "ridge")
      0.3 + max(nchar(as.character(levels(df$group)))) * 0.07 else 0
    facet_w   <- if (!is.null(facet.by)) min(n_facets, 4) * 3.2 else 4.5
    pdf_w <- width  %||% (facet_w + y_axis_w + lgd_w)
    pdf_h <- height %||% max(if (style == "ridge") 5 else 4, lgd_h)

    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(p)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath,
            " (", round(pdf_w, 1), " x ", round(pdf_h, 1), " in)")

    ctx <- .nk_legend_context(seurat_object)
    .write_legend_sidecar(fpath, paste0(
      .nk_obs_sentence(ctx, group.by = group.by, split.by = facet.by), ". ",
      if (style == "ridge") "Ridge" else "Density", " plot of '", feat, "'",
      if (use_split_fill) paste0(", filled by ", split.by) else "",
      if (!is.null(facet.by)) paste0(", faceted by ", facet.by) else "",
      if (!is.null(min_expr)) paste0(
        ". Cells with expression <= ", min_expr,
        " were excluded to avoid the dropout pile-up at zero")
      else "", "."
    ))
  }

  # ── Loop over features: one PDF per gene ────────────────────────────────────
  plots <- list()
  for (feat in features) {
    res <- .build_one(feat)
    if (!is.null(output_dir)) .save_one(feat, res)
    plots[[feat]] <- res$plot
  }

  if (!is.null(output_dir)) return(invisible(if (multi) plots else plots[[1]]))
  if (!multi) return(plots[[1]])
  invisible(plots)
}
