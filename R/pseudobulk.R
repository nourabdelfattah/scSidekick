# =============================================================================
# pseudobulk.R
#
# ComputePseudobulk() - aggregate cell-level expression to donor-level means
#                       and cache the result in seurat_object@misc$pseudobulk.
#
# PlotPseudoBulk()    - box + dot plots at pseudobulk (donor) resolution.
#                       Dots are one per donor (or donor x celltype, etc.).
#                       Pulls from cache when available; computes lazily if not.
# =============================================================================


# Internal helper: make a character key from selected columns of a data frame.
# Used for matching pseudobulk rows across cache updates.
.pb_key <- function(df, cols)
  do.call(paste, c(lapply(cols, function(g) as.character(df[[g]])),
                   list(sep = "\t")))


# =============================================================================
# ComputePseudobulk
# =============================================================================

#' Compute pseudobulk (donor-level mean) expression
#'
#' Aggregates cell-level log-normalized expression to the level of unique
#' combinations of metadata columns (e.g. donor, or donor x cell type).
#' The result is stored in \code{seurat_object@misc$pseudobulk} under a key
#' derived from the sorted \code{group.by} columns so that
#' \code{\link{PlotPseudoBulk}} can retrieve it instantly on subsequent calls.
#'
#' Aggregation uses a single sparse matrix multiplication
#' (\code{genes x cells} . \code{cells x groups}), making it efficient
#' for BPCells / on-disk assays: the matrix is streamed through once.
#'
#' When called again with the same \code{group.by} but different (or
#' additional) genes, only the newly requested genes are computed and appended
#' to the existing cache - already-cached genes are never recomputed.
#'
#' @param seurat_object A Seurat object.
#' @param group.by Character vector.  Metadata columns that together define
#'   one pseudobulk sample.  For donor-level plots pass
#'   \code{c("Donor.ID")}.  For donor x cell-type plots pass
#'   \code{c("Donor.ID", "CellType")}.
#' @param genes One of:
#'   \describe{
#'     \item{\code{"variable"} (default)}{Use the variable features stored by
#'       \code{FindVariableFeatures()}.}
#'     \item{\code{NULL}}{All genes in the assay (can be slow).}
#'     \item{Character vector}{Specific gene names.}
#'   }
#' @param assay Character.  Seurat assay.  Default \code{"RNA"}.
#' @param layer Character.  Assay layer / slot.  Default \code{"data"}
#'   (log-normalized).
#' @param min_cells Integer.  Drop pseudobulk groups with fewer than this many
#'   cells.  Default \code{0L} (keep all groups).
#' @param verbose Logical.  Print progress messages.  Default \code{TRUE}.
#'
#' @return The \code{seurat_object} with the pseudobulk data frame stored in
#'   \code{seurat_object@misc$pseudobulk[[cache_key]]}.  The data frame has
#'   one row per unique group combination and columns:
#'   \code{group.by columns}, \code{n_cells}, then one column per gene.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Per-donor pseudobulk (variable genes)
#' SeuratObj <- ComputePseudobulk(SeuratObj,
#'   group.by = "Donor.ID")
#'
#' # Per-donor x per-cell-type (specific genes)
#' SeuratObj <- ComputePseudobulk(SeuratObj,
#'   group.by = c("Donor.ID", "CellType"),
#'   genes    = c("STAT3", "APOE", "CD3E"))
#' }
ComputePseudobulk <- function(seurat_object,
                               group.by,
                               genes     = "variable",
                               assay     = "RNA",
                               layer     = "data",
                               min_cells = 0L,
                               verbose   = TRUE) {

  meta <- seurat_object@meta.data

  # ── Validate group.by ─────────────────────────────────────────────────────
  missing_cols <- setdiff(group.by, colnames(meta))
  if (length(missing_cols) > 0L)
    stop("group.by column(s) not found in metadata: ",
         paste(missing_cols, collapse = ", "))

  # ── Determine gene set ────────────────────────────────────────────────────
  if (identical(genes, "variable")) {
    genes_to_use <- Seurat::VariableFeatures(seurat_object)
    if (length(genes_to_use) == 0L)
      stop("No variable features found. Run FindVariableFeatures() first, ",
           "or supply genes = c('GENE1', 'GENE2') explicitly.")
    if (verbose)
      message("scSidekick ComputePseudobulk: Using ",
              length(genes_to_use), " variable features.")
  } else if (is.null(genes)) {
    genes_to_use <- tryCatch(
      rownames(SeuratObject::GetAssay(seurat_object, assay = assay)),
      error = function(e) rownames(seurat_object[[assay]])
    )
  } else {
    genes_to_use <- as.character(genes)
  }

  # ── Get expression matrix (lazy - safe for BPCells) ───────────────────────
  mat <- tryCatch(
    SeuratObject::LayerData(seurat_object, assay = assay, layer = layer),
    error = function(e1) tryCatch(
      Seurat::GetAssayData(seurat_object, assay = assay, layer = layer),
      error = function(e2)
        Seurat::GetAssayData(seurat_object, assay = assay, slot = layer)
    )
  )

  cells       <- intersect(colnames(mat), rownames(meta))
  genes_found <- intersect(genes_to_use, rownames(mat))

  if (length(cells) == 0L)
    stop("No cells overlap between assay '", assay, "' and metadata.")
  if (length(genes_found) == 0L)
    stop("None of the requested genes found in assay '", assay, "'.")
  if (length(genes_found) < length(genes_to_use) && verbose)
    message("scSidekick: ",
            length(genes_to_use) - length(genes_found),
            " gene(s) not found in assay '", assay, "' - skipped.")

  # ── Check cache: only compute genes not already stored ────────────────────
  cache_key <- paste(sort(group.by), collapse = "+")
  existing  <- seurat_object@misc$pseudobulk[[cache_key]]

  if (!is.null(existing)) {
    existing_genes   <- setdiff(colnames(existing), c(group.by, "n_cells"))
    genes_to_compute <- setdiff(genes_found, existing_genes)

    if (length(genes_to_compute) == 0L) {
      if (verbose)
        message("scSidekick ComputePseudobulk: All ", length(genes_found),
                " gene(s) already in cache for '", cache_key,
                "'. Nothing to do.")
      return(seurat_object)
    }
    if (verbose)
      message("scSidekick ComputePseudobulk: Adding ",
              length(genes_to_compute), " new gene(s) to existing cache '",
              cache_key, "'.")
  } else {
    genes_to_compute <- genes_found
  }

  # ── Subset matrix to requested genes and available cells ──────────────────
  mat_sub  <- mat[genes_to_compute, cells, drop = FALSE]
  meta_sub <- meta[cells, group.by, drop = FALSE]

  # ── Assign each cell to a group index ────────────────────────────────────
  # Get the unique group combinations directly from the metadata subset -
  # no string parsing needed.
  meta_keys  <- .pb_key(meta_sub, group.by)
  uniq_keys  <- unique(meta_keys)
  n_grps     <- length(uniq_keys)
  grp_idx    <- match(meta_keys, uniq_keys)   # integer index per cell
  n_cells_v  <- length(cells)

  # Recover group data frame from the first occurrence of each key
  first_rows <- match(uniq_keys, meta_keys)
  group_df   <- meta_sub[first_rows, group.by, drop = FALSE]
  rownames(group_df) <- NULL

  # Group sizes via tabulate (much faster than building a sparse matrix)
  grp_sizes <- as.integer(tabulate(grp_idx, nbins = n_grps))

  if (verbose)
    message("scSidekick ComputePseudobulk: ",
            format(n_cells_v, big.mark = ","), " cells  ->  ",
            n_grps, " pseudobulk groups x ",
            length(genes_to_compute), " genes.")

  # ── Compute group means ───────────────────────────────────────────────────
  # Strategy depends on gene count:
  #
  # FEW genes (≤ 50): one sequential BPCells scan to materialize genes x cells,
  #   then pure-R group means.  Same access pattern as PercentageFeatureSet -
  #   fast even for BPCells on network drives.
  #
  # MANY genes (> 50): sparse indicator matrix multiplication.  Avoids loading
  #   a huge dense matrix into memory; relies on BPCells streaming for %*%.
  #
  .FEW_GENES <- 50L

  if (length(genes_to_compute) <= .FEW_GENES) {
    # Materialise via BPCells: as.numeric() on a raw S4 BPCells object fails,
    # so go through as.matrix() first (triggers BPCells' dense conversion),
    # then as.numeric() flattens to a plain R vector, then matrix() rebuilds
    # with guaranteed 2D shape - handles single-gene dim-drop and all assay
    # types (BPCells, dgCMatrix, standard R matrix).
    expr_mat <- matrix(
      as.numeric(as.matrix(mat_sub)),
      nrow = length(genes_to_compute),
      ncol = n_cells_v
    )
    rownames(expr_mat) <- genes_to_compute

    # Group means in R memory - fast, no large matrix construction needed
    mean_mat <- vapply(seq_len(n_grps), function(k) {
      cols <- which(grp_idx == k)
      if (length(cols) == 0L) return(rep(0, length(genes_to_compute)))
      if (length(cols) == 1L) return(as.numeric(expr_mat[, cols]))
      rowMeans(expr_mat[, cols, drop = FALSE])
    }, numeric(length(genes_to_compute)))
    rownames(mean_mat) <- genes_to_compute

  } else {
    # Large gene set: sparse indicator matrix multiplication
    indicator <- Matrix::sparseMatrix(
      i    = seq_len(n_cells_v),
      j    = grp_idx,
      x    = 1,
      dims = c(n_cells_v, n_grps)
    )
    sum_mat  <- mat_sub %*% indicator
    mean_mat <- sweep(as.matrix(sum_mat), 2L, grp_sizes, "/")
  }

  # ── Restore factor levels where applicable ────────────────────────────────
  for (g in group.by) {
    if (is.factor(meta[[g]]))
      group_df[[g]] <- factor(group_df[[g]], levels = levels(meta[[g]]))
  }
  group_df$n_cells <- grp_sizes

  # ── Apply min_cells filter ────────────────────────────────────────────────
  if (min_cells > 0L) {
    keep  <- group_df$n_cells >= min_cells
    n_drop <- sum(!keep)
    if (n_drop > 0L && verbose)
      message("scSidekick: Dropped ", n_drop,
              " pseudobulk group(s) with < ", min_cells, " cells.")
    group_df  <- group_df[keep,  , drop = FALSE]
    mean_mat  <- mean_mat[, keep,  drop = FALSE]
  }

  # ── Assemble data frame: group columns + gene expression columns ──────────
  pb_df <- cbind(group_df,
                 as.data.frame(t(mean_mat), stringsAsFactors = FALSE))
  rownames(pb_df) <- NULL

  # ── Merge with existing cache (incremental update) ────────────────────────
  if (!is.null(existing)) {
    new_keys <- .pb_key(pb_df, group.by)
    old_keys <- .pb_key(existing, group.by)
    idx      <- match(old_keys, new_keys)

    if (anyNA(idx))
      warning("scSidekick: Some groups in the existing cache were not found ",
              "in the new computation. The cache may be stale - consider ",
              "rerunning ComputePseudobulk() to rebuild from scratch.")

    extra_cols <- pb_df[idx, genes_to_compute, drop = FALSE]
    rownames(extra_cols) <- NULL
    pb_df <- cbind(existing, extra_cols)
  }

  # ── Store in misc ─────────────────────────────────────────────────────────
  if (is.null(seurat_object@misc$pseudobulk))
    seurat_object@misc$pseudobulk <- list()
  seurat_object@misc$pseudobulk[[cache_key]] <- pb_df

  if (verbose)
    message("scSidekick ComputePseudobulk: Done. Access via:\n",
            "  seurat_object@misc$pseudobulk[[\"", cache_key, "\"]]")

  seurat_object
}


# =============================================================================
# PlotPseudoBulk
# =============================================================================

#' Pseudobulk box + dot plots for gene expression
#'
#' Plots mean gene expression per donor (or per donor x cell type, etc.) as
#' box plots with individual donor dots overlaid - the standard pseudobulk
#' visualization that avoids pseudo-replication.
#'
#' Each dot represents one unique \code{donor.by} value (e.g. one patient).
#' Boxes show median, IQR, and 1.5x IQR whiskers.  Dots are black.
#'
#' On first call, pseudobulk means are computed on the fly and a hint is
#' printed suggesting \code{\link{ComputePseudobulk}} for BPCells / large
#' datasets.  Pre-warm the cache with \code{ComputePseudobulk()} to make
#' subsequent plots instant.
#'
#' @param seurat_object A Seurat object.
#' @param features Character vector.  Gene names to plot.
#' @param group.by Character.  Categorical metadata column shown on the x-axis
#'   (e.g. \code{"Cognitive.Status"} or \code{"Condition"}).
#' @param donor.by Character.  Metadata column that identifies one biological
#'   replicate (e.g. \code{"Donor.ID"}).  One dot per unique value.
#' @param split.by Character or \code{NULL}.  Column facets within each panel
#'   (e.g. \code{"CellType"}).
#' @param row.by Character or \code{NULL}.  Row facets within each panel.
#'   Combined with \code{split.by} gives a two-way grid.
#' @param exclude Named list of values to drop before plotting.  Each name is
#'   a column and each value is a character vector of levels to remove.
#'   Example: \code{list(Cognitive.Status = "Reference")}.
#' @param min_cells Integer.  Drop pseudobulk groups (donor x group
#'   combinations) with fewer than this many cells.  Default \code{0L}.
#' @param show_n Logical.  Annotate each x-axis group with the number of
#'   donors (or pseudobulk samples) in that group after all filtering.
#'   Shown as \code{n=X} below the x-axis ticks.  Default \code{TRUE}.
#' @param add_stats Logical.  Add Wilcoxon significance brackets via
#'   \pkg{ggsignif}.  Default \code{FALSE}.
#' @param comparisons List of length-2 character vectors specifying pairs.
#'   \code{NULL} (default) tests all pairs, or all vs \code{ref_group}.
#' @param ref_group Character or \code{NULL}.  Reference group for comparisons.
#' @param hide_ns Logical.  Suppress \code{"ns"} brackets.  Default
#'   \code{FALSE}.
#' @param label_format One of \code{"stars"} (default) or \code{"p.format"}.
#' @param point_size Numeric.  Size of donor dots.  Default \code{1.5}.
#' @param jitter_width Numeric.  Jitter width for donor dots.  Default
#'   \code{0.15}.
#' @param box_width Numeric.  Width of the box geom.  Default \code{0.5}.
#' @param alpha Numeric.  Fill transparency of boxes.  Default \code{0.7}.
#' @param ncol Integer or \code{NULL}.  Columns in the patchwork assembly.
#'   \code{NULL} = \code{min(n_features, 3)}.
#' @param colors Named character vector of fill colors for \code{group.by}
#'   levels.  \code{NULL} auto-resolves from \code{PrepObject} or
#'   \code{SelectColors()}.
#' @param assay Character.  Seurat assay.  Default \code{"RNA"}.
#' @param layer Character.  Assay layer.  Default \code{"data"}.
#' @param output_dir Character or \code{NULL}.  Save directory.  Walks up from
#'   \code{PrepObject}.
#' @param object_name Character.  Filename prefix.  Walks up from
#'   \code{PrepObject}.
#' @param pdf_width Numeric or \code{NULL}.  Override auto-calculated PDF width.
#' @param pdf_height Numeric or \code{NULL}.  Override auto-calculated PDF height.
#'
#' @return A \code{ggplot2} object (single feature) or a \code{patchwork}
#'   combined plot (multiple features).  Saved as PDF when \code{output_dir}
#'   is available.
#' @export
#'
#' @examples
#' \dontrun{
#' # Pre-warm cache (recommended for BPCells / large datasets)
#' SeuratObj <- ComputePseudobulk(SeuratObj,
#'   group.by = c("Donor.ID", "Cognitive.Status"))
#'
#' # Per-donor boxplot, exclude reference group, add stats
#' PlotPseudoBulk(SeuratObj,
#'   features  = c("STAT3", "APOE"),
#'   group.by  = "Cognitive.Status",
#'   donor.by  = "Donor.ID",
#'   exclude   = list(Cognitive.Status = "Reference"),
#'   add_stats = TRUE)
#'
#' # Per-donor x per-cell-type (split into column facets)
#' SeuratObj <- ComputePseudobulk(SeuratObj,
#'   group.by = c("Donor.ID", "Cognitive.Status", "CellType"))
#'
#' PlotPseudoBulk(SeuratObj,
#'   features   = "APOE",
#'   group.by   = "Cognitive.Status",
#'   donor.by   = "Donor.ID",
#'   split.by   = "CellType",
#'   min_cells  = 10,
#'   add_stats  = TRUE,
#'   ref_group  = "Control")
#' }
PlotPseudoBulk <- function(seurat_object,
                            features,
                            group.by,
                            donor.by,
                            split.by     = NULL,
                            row.by       = NULL,
                            exclude      = NULL,
                            min_cells    = 0L,
                            show_n       = TRUE,
                            add_stats    = FALSE,
                            comparisons  = NULL,
                            ref_group    = NULL,
                            hide_ns      = FALSE,
                            stat_alpha   = 0.05,
                            label_format = c("stars", "p.format"),
                            point_size   = 1.5,
                            jitter_width = 0.15,
                            box_width    = 0.5,
                            alpha        = 0.7,
                            ncol         = NULL,
                            colors       = NULL,
                            assay        = "RNA",
                            layer        = "data",
                            output_dir   = NULL,
                            object_name  = "",
                            pdf_width    = NULL,
                            pdf_height   = NULL) {

  label_format <- match.arg(label_format)

  # ── Walk-up PrepObject defaults ───────────────────────────────────────────
  output_dir  <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  meta <- seurat_object@meta.data

  # ── Validate columns ──────────────────────────────────────────────────────
  all_meta_cols <- unique(c(donor.by, group.by, split.by, row.by))
  missing_cols  <- setdiff(all_meta_cols, colnames(meta))
  if (length(missing_cols) > 0L)
    stop("Column(s) not found in metadata: ",
         paste(missing_cols, collapse = ", "))

  # ── Full pseudobulk grouping (donor + all facet / x-axis columns) ─────────
  pb_group_by <- unique(c(donor.by, group.by, split.by, row.by))
  cache_key   <- paste(sort(pb_group_by), collapse = "+")

  # ── Check cache FIRST - before any expensive assay operations ───────────────
  # For BPCells objects, Features() / GetAssay() can be slow even for rownames.
  # If all requested features are already cached we skip assay access entirely.
  cache    <- seurat_object@misc$pseudobulk[[cache_key]]
  in_cache <- if (!is.null(cache)) intersect(features, colnames(cache))
              else character(0)
  need_compute <- setdiff(features, in_cache)

  message("scSidekick PlotPseudoBulk: cache key = \"", cache_key, "\" | ",
          length(in_cache), " gene(s) found in cache, ",
          length(need_compute), " need computation.")

  # ── Validate gene features (only when something is missing from cache) ─────
  valid_features <- in_cache   # fast path: already have everything

  if (length(need_compute) > 0L) {
    # Diagnose key mismatch before the slow Features() call
    stored_keys <- names(seurat_object@misc$pseudobulk)
    if (length(stored_keys) > 0L && !cache_key %in% stored_keys) {
      message(
        "scSidekick PlotPseudoBulk: Cache miss - key mismatch.\n",
        "  Expected key : \"", cache_key, "\"\n",
        "  Stored key(s): ", paste0('"', stored_keys, '"', collapse = ", "), "\n",
        "  Fix: re-run ComputePseudobulk() with exactly:\n",
        "    group.by = c(", paste0('"', pb_group_by, '"', collapse = ", "), ")"
      )
    }

    # SeuratObject::Features() is the canonical v5 / BPCells-safe call.
    all_gene_names <- tryCatch(
      SeuratObject::Features(seurat_object, assay = assay),
      error = function(e1) tryCatch(
        rownames(SeuratObject::GetAssay(seurat_object, assay = assay)),
        error = function(e2) tryCatch({
          lm <- SeuratObject::LayerData(seurat_object, assay = assay,
                                        layer = layer)
          rownames(lm)
        }, error = function(e3) character(0))
      )
    )
    valid_missing <- intersect(need_compute, all_gene_names)
    not_found     <- setdiff(need_compute, all_gene_names)

    if (length(not_found) > 0L)
      warning("Feature(s) not found in assay '", assay, "': ",
              paste(not_found, collapse = ", "), " - skipped.")

    valid_features <- c(in_cache, valid_missing)
  }

  if (length(valid_features) == 0L)
    stop("No valid gene features found. Check spelling and assay name ",
         "(assay = \"", assay, "\").")

  # ── Compute missing genes on the fly if cache incomplete ──────────────────
  missing_genes <- setdiff(valid_features, in_cache)

  if (length(missing_genes) > 0L) {
    message(
      "scSidekick PlotPseudoBulk: ",
      paste0('"', missing_genes, '"', collapse = ", "),
      " not found in cache (not a variable feature, or cache was built ",
      "without it).\n",
      "  Computing on the fly - this is slow for BPCells / large datasets.\n",
      "  To pre-cache these gene(s) permanently run:\n",
      "  SeuratObj <- ComputePseudobulk(SeuratObj,\n",
      "    group.by = c(", paste0('"', pb_group_by, '"', collapse = ", "), "),\n",
      "    genes    = c(", paste0('"', missing_genes, '"', collapse = ", "), "))"
    )
    tmp <- ComputePseudobulk(
      seurat_object = seurat_object,
      group.by      = pb_group_by,
      genes         = missing_genes,
      assay         = assay,
      layer         = layer,
      min_cells     = 0L,
      verbose       = FALSE
    )
    cache <- tmp@misc$pseudobulk[[cache_key]]
  }

  # ── Apply min_cells filter ────────────────────────────────────────────────
  if (min_cells > 0L && "n_cells" %in% colnames(cache)) {
    keep  <- cache$n_cells >= min_cells
    n_drop <- sum(!keep)
    if (n_drop > 0L)
      message("scSidekick: Dropped ", n_drop,
              " pseudobulk group(s) with < ", min_cells, " cells.")
    cache <- cache[keep, , drop = FALSE]
  }

  # ── Apply exclusions ──────────────────────────────────────────────────────
  if (!is.null(exclude)) {
    n_before_excl <- nrow(cache)
    for (col in names(exclude)) {
      if (!col %in% colnames(cache)) {
        warning("Exclusion column '", col,
                "' not found in pseudobulk data - skipping.")
        next
      }
      present_vals <- unique(as.character(cache[[col]]))
      matched      <- intersect(as.character(exclude[[col]]), present_vals)
      unmatched    <- setdiff(as.character(exclude[[col]]), present_vals)

      if (length(unmatched) > 0L)
        warning("scSidekick PlotPseudoBulk: exclude value(s) not found in '",
                col, "': ", paste0('"', unmatched, '"', collapse = ", "), "\n",
                "  Available values: ",
                paste0('"', sort(present_vals), '"', collapse = ", "))

      if (length(matched) > 0L) {
        n_before_col <- nrow(cache)
        cache <- cache[
          !as.character(cache[[col]]) %in% matched,
          , drop = FALSE
        ]
        message("scSidekick: Excluded ", n_before_col - nrow(cache),
                " pseudobulk group(s) where ", col, " = ",
                paste0('"', matched, '"', collapse = ", "))
      }
    }
    message("scSidekick PlotPseudoBulk: ",
            nrow(cache), " pseudobulk samples remain (",
            n_before_excl - nrow(cache), " excluded).")
  }

  # ── Resolve group.by levels and fill colors ───────────────────────────────
  grp_vals <- as.character(cache[[group.by]])
  grp_lvls <- if (is.factor(cache[[group.by]])) levels(cache[[group.by]])
              else sort(unique(grp_vals))
  # Drop levels with no remaining data after exclusions
  grp_lvls <- intersect(grp_lvls, unique(grp_vals))

  fill_colors <- colors
  if (is.null(fill_colors))
    fill_colors <- tryCatch(.nk_colors(seurat_object, group.by),
                            error = function(e) NULL)
  if (is.null(fill_colors))
    fill_colors <- SelectColors(factor(grp_vals, levels = grp_lvls),
                                palette = "all")

  # ── Cached facet levels ───────────────────────────────────────────────────
  split_lvls <- if (!is.null(split.by))
    if (is.factor(cache[[split.by]])) levels(cache[[split.by]])
    else sort(unique(as.character(cache[[split.by]])))
  else NULL

  row_lvls <- if (!is.null(row.by))
    if (is.factor(cache[[row.by]])) levels(cache[[row.by]])
    else sort(unique(as.character(cache[[row.by]])))
  else NULL

  # ── Build one panel per feature ───────────────────────────────────────────
  plot_list <- lapply(valid_features, function(feat) {

    df <- data.frame(
      Value = as.numeric(cache[[feat]]),
      Group = factor(as.character(cache[[group.by]]), levels = grp_lvls),
      Donor = as.character(cache[[donor.by]]),
      stringsAsFactors = FALSE
    )
    if (!is.null(split.by))
      df$Split <- factor(as.character(cache[[split.by]]), levels = split_lvls)
    if (!is.null(row.by))
      df$Row   <- factor(as.character(cache[[row.by]]),   levels = row_lvls)
    df <- df[!is.na(df$Value), , drop = FALSE]

    y_lab <- paste0(feat, "\n(mean log-norm. per ", donor.by, ")")

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Value, fill = Group))

    # Box (group-colored fill, no outlier points - jitter covers them)
    p <- p + ggplot2::geom_boxplot(
      width         = box_width,
      outlier.shape = NA,
      alpha         = alpha,
      linewidth     = 0.4,
      color         = "grey20"
    )

    # Individual donor dots - always black
    p <- p + ggplot2::geom_jitter(
      width       = jitter_width,
      size        = point_size,
      color       = "black",
      alpha       = 0.8,
      show.legend = FALSE
    )

    # ── Statistics brackets ──────────────────────────────────────────────────
    if (add_stats && requireNamespace("ggsignif", quietly = TRUE)) {
      grp_present <- levels(droplevels(df$Group))

      all_pairs <- if (!is.null(comparisons)) {
        comparisons
      } else if (!is.null(ref_group) && ref_group %in% grp_present) {
        lapply(setdiff(grp_present, ref_group), function(g) c(ref_group, g))
      } else {
        utils::combn(grp_present, 2L, simplify = FALSE)
      }

      if (length(all_pairs) > 0L) {
        sig_args <- list(
          comparisons   = all_pairs,
          test          = "wilcox.test",
          test.args     = list(exact = FALSE),
          step_increase = 0.08,
          tip_length    = 0.01,
          size          = 0.4,
          vjust         = 0.3,
          textsize      = if (label_format == "stars") 3.5 else 2.5
        )
        sig_args$map_signif_level <-
          if (label_format == "stars")
            c("***" = stat_alpha / 50, "**" = stat_alpha / 5,
              "*"   = stat_alpha,      "ns" = 1)
          else FALSE
        if (hide_ns)
          sig_args$map_signif_level <-
            c("***" = stat_alpha / 50, "**" = stat_alpha / 5,
              "*"   = stat_alpha)

        p <- p + do.call(ggsignif::geom_signif, sig_args) +
          ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.05, 0.35))
          )
      }
    }

    # ── N annotations (number of donors per group per facet) ────────────────
    if (show_n) {
      n_cols  <- c("Group",
                   if (!is.null(split.by)) "Split",
                   if (!is.null(row.by))   "Row")
      n_df    <- df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(n_cols))) |>
        dplyr::summarise(n = dplyr::n(), .groups = "drop")
      n_df$n_label <- paste0("n=", n_df$n)

      p <- p + ggplot2::geom_text(
        data            = n_df,
        ggplot2::aes(x = Group, y = -Inf, label = n_label),
        vjust           = 1.5,
        size            = 2.8,
        color           = "grey40",
        inherit.aes     = FALSE
      ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(plot.margin = ggplot2::margin(t = 5, r = 5, b = 14,
                                                      l = 10, unit = "mm"))
    }

    # ── Scales, labels, theme ────────────────────────────────────────────────
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

    # ── Facets ───────────────────────────────────────────────────────────────
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
  })
  names(plot_list) <- valid_features

  # ── Assemble with patchwork ───────────────────────────────────────────────
  n_feat   <- length(plot_list)
  ncol     <- ncol %||% min(n_feat, 3L)
  combined <- if (n_feat == 1L) plot_list[[1L]]
              else patchwork::wrap_plots(plot_list, ncol = ncol)

  # ── Auto-save PDF + .legend sidecar ──────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    n_groups <- length(grp_lvls)
    n_splits <- if (!is.null(split.by)) length(split_lvls) else 1L
    n_rows_f <- if (!is.null(row.by))   length(row_lvls)   else 1L
    nrow_pw  <- ceiling(n_feat / ncol)

    facet_col_w <- max(n_groups * 0.55, 2.0)
    panel_w     <- facet_col_w * n_splits + 0.8
    panel_h     <- max(2.8, 2.5) * n_rows_f + 1.0

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
      "PseudoBulk"
    )
    fname <- gsub("[^A-Za-z0-9._-]", "_", paste(parts, collapse = "_"))
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))

    message("scSidekick PlotPseudoBulk: Rendering plot to PDF (",
            round(pdf_w, 1), " x ", round(pdf_h, 1), " in)...")
    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(combined)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath,
            " (", round(pdf_w, 1), " x ", round(pdf_h, 1), " in)")

    n_donors <- length(unique(cache[[donor.by]]))
    .write_legend_sidecar(fpath, paste0(
      "Pseudobulk box plot of ", paste(valid_features, collapse = ", "), ". ",
      "Each dot represents the mean log-normalized expression for one ",
      donor.by, " (n = ", n_donors, " donors), ",
      "grouped by ", group.by,
      if (!is.null(split.by)) paste0(", split by ", split.by) else "",
      if (!is.null(row.by))   paste0(", rows by ", row.by)   else "",
      ". Box plot elements: centre line = median; box limits = 25th-75th percentile (IQR); whiskers extend to the furthest observation within 1.5x IQR from the box; outliers beyond this range are not shown.",
      if (min_cells > 0L)
        paste0(" Groups with fewer than ", min_cells,
               " cells per donor were excluded.")
      else "",
      if (add_stats)
        " Statistical comparisons: Wilcoxon rank-sum test on donor-level means."
      else "",
      if (!is.null(exclude) && length(exclude) > 0L)
        paste0(" Excluded: ",
               paste(mapply(function(col, vals)
                 paste0(col, " = ", paste(vals, collapse = ", ")),
                 names(exclude), exclude), collapse = "; "), ".")
      else "",
      if (nchar(object_name) > 0) paste0(" Dataset: ", object_name, ".") else ""
    ))
  }

  combined
}
