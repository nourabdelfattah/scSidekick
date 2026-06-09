# =============================================================================
# scSidekick dot plot functions
#
# SplitDotPlot  - dot plot split by a metadata variable, with numbered labels
#                 and optional hierarchical clustering of cell groups
# SplitDotPlot2 - extends SplitDotPlot with gene-axis clustering and a
#                 minimum-expression filter
# FastDotPlot   - like SplitDotPlot2 but accepts a gene name pattern (regex)
#                 and can slice the gene dendrogram into k labelled groups
# FastDotPlot2  - FastDotPlot + optional EnrichR enrichment header aligned
#                 above each gene-cluster panel
#
# Shared internal helpers
# .build_id_map()   - create zero-padded numeric labels for cell groups
# .clamp()          - clamp a numeric vector to [lo, hi]
# .dotplot_theme()  - shared ggplot2 theme for all dot plots
# .cleanup_meta()   - on.exit handler to remove the temporary interaction col
# =============================================================================

# --------------------------------------------------------------------------- #
# Internal helpers                                                             #
# --------------------------------------------------------------------------- #

# Sorts cell groups respecting existing factor levels, returns zero-padded map
.build_id_map <- function(plot_data, raw_group_col) {
  existing_groups <- unique(plot_data$Y_Axis_Group)
  if (is.factor(raw_group_col)) {
    sorted <- levels(raw_group_col)[levels(raw_group_col) %in% existing_groups]
    remaining <- setdiff(existing_groups, sorted)
    if (length(remaining)) sorted <- c(sorted, sort(remaining))
  } else {
    sorted <- stringr::str_sort(existing_groups, numeric = TRUE)
  }
  id_map <- stats::setNames(sprintf("%02d", seq_along(sorted)), sorted)
  list(sorted = sorted, id_map = id_map)
}

# Clamp numeric vector to [lo, hi]
.clamp <- function(x, lo, hi) pmax(pmin(x, hi), lo)

# Shared theme used by all dot-plot functions
.dotplot_theme <- function() {
  ggplot2::theme(
    panel.grid.major  = ggplot2::element_blank(),
    panel.grid.minor  = ggplot2::element_blank(),
    panel.background  = ggplot2::element_blank(),
    panel.border      = ggplot2::element_rect(fill = NA, colour = "grey90",
                                              linewidth = 1),
    strip.background  = ggplot2::element_rect(fill = "grey95", colour = "grey90"),
    axis.text.x       = ggplot2::element_text(angle = 45, hjust = 1,
                                              vjust = 1, size = 10),
    axis.text.y       = ggplot2::element_text(size = 10),
    axis.line         = ggplot2::element_line(colour = "grey90",
                                             linewidth = ggplot2::rel(1)),
    strip.text.y      = ggplot2::element_text(face = "bold", size = 10,
                                             margin = ggplot2::margin(l = 20)),
    strip.text.x      = ggplot2::element_text(face = "bold", size = 10),
    panel.spacing     = ggplot2::unit(1, "lines")
  )
}

# Helper: remove temp column from Seurat metadata when function exits
.make_cleanup <- function(seurat_object, col) {
  function() {
    if (col %in% colnames(seurat_object@meta.data))
      seurat_object@meta.data[[col]] <<- NULL
  }
}

# Helper: build a save path for a dot plot.
#   file_name : explicit base name from the caller (highest priority; NULL = ignore)
#   auto_base : a deduced base (e.g. the markers_df variable name or gene pattern)
#   suffix    : function tag appended when no explicit file_name (e.g. "SplitDotPlot")
# When file_name is given it is used verbatim (object/subset prefixes are skipped),
# so the user gets exactly the name they asked for.
.dotplot_path <- function(output_dir, seurat_object, file_name, auto_base, suffix) {
  if (!is.null(file_name) && nzchar(file_name)) {
    base <- file_name
  } else {
    obj   <- .nk_setting(seurat_object, "object_name") %||% ""
    sub   <- .nk_setting(seurat_object, "subset_name") %||% ""
    # auto_base is used as the descriptive middle (df name or pattern) when it is
    # informative; otherwise fall back to the function suffix alone.
    mid   <- if (!is.null(auto_base) && nzchar(auto_base)) auto_base else NULL
    parts <- c(if (nchar(obj) > 0) obj, if (nchar(sub) > 0) sub, mid, suffix)
    base  <- paste(parts, collapse = "_")
  }
  fname <- gsub("[^A-Za-z0-9._-]", "_", base)
  file.path(output_dir, paste0(fname, ".pdf"))
}

# Helper: is a deparsed argument a usable variable name (vs an inline expression)?
# Rejects multi-line deparse, function calls, $ / [ indexing, and over-long names.
.usable_obj_name <- function(x) {
  length(x) == 1L && nzchar(x) && nchar(x) <= 60L &&
    !grepl("[()$\\[\\]]|, |::", x) &&
    grepl("^[A-Za-z.][A-Za-z0-9._]*$", x)
}

# --------------------------------------------------------------------------- #
# SplitDotPlot                                                                 #
# --------------------------------------------------------------------------- #

#' Dot plot split by a metadata variable with numbered group labels
#'
#' Generates a dot plot where dot color encodes scaled average expression and
#' dot size encodes percent of cells expressing each gene. Cell groups
#' (y-axis) are numbered and can be hierarchically clustered. Genes (x-axis)
#' are grouped by a gene-metadata column and faceted accordingly.
#'
#' A temporary metadata column (`Temp_Interaction_Var`) is added to the
#' Seurat object during computation and automatically removed on exit.
#'
#' @param seurat_object A Seurat object.
#' @param markers_df A data frame containing at least two columns: gene names
#'   and gene group labels.
#' @param gene_column Character. Column in `markers_df` with gene names.
#'   Default `"Genes"`.
#' @param gene_group_column Character. Column in `markers_df` with gene group
#'   labels (used for x-axis faceting). Default `"CellType"`.
#' @param group.by Character. Metadata column for the y-axis (cell groups).
#'   Default `"GlobalAssignment"`.
#' @param split.by Character or `NULL`. Metadata column to split rows by.
#'   `NULL` shows all groups together. Default `"CancerType"`.
#' @param ClusterFeatures Logical. Hierarchically cluster cell groups
#'   (y-axis)? Default `FALSE`.
#' @param col.min,col.max Numeric. Color scale clamp limits for scaled average
#'   expression. Defaults `0` and `2.5`.
#' @param scale.min,scale.max Numeric. Size scale clamp limits for percent
#'   expressed. Defaults `5` and `100`.
#' @param dot.scale Numeric. Maximum dot radius. Default `6`.
#' @param cols Ignored (kept for API compatibility). Color palette is viridis
#'   plasma.
#'
#' @param output_dir Character or \code{NULL}. Directory to save a PDF.
#'   \code{NULL} (default) returns the plot without saving.  When
#'   \code{NULL}, the function walks up to the \code{output_dir} stored by
#'   \code{\link{PrepObject}} — unless \code{AutoSavePlots = FALSE} was set
#'   there.
#' @param file_name Character or \code{NULL}. Base name (no extension) for the
#'   saved PDF. \code{NULL} (default) auto-deduces the name from the
#'   \code{markers_df} variable name, falling back to
#'   \code{object_name_group.by_SplitDotPlot}.
#' @return A ggplot2 object (invisibly when saved to disk).
#' @export
SplitDotPlot <- function(seurat_object,
                          markers_df,
                          gene_column       = "Genes",
                          gene_group_column = "CellType",
                          group.by          = "seurat_clusters",
                          split.by          = NULL,
                          # deprecated short aliases
                          gene_col          = NULL,
                          gene_group_col    = NULL,
                          ClusterFeatures = FALSE,
                          col.min        = 0,
                          col.max        = 2.5,
                          scale.min      = 5,
                          scale.max      = 100,
                          dot.scale      = 6,
                          cols           = "RdBu",
                          output_dir     = NULL,
                          file_name      = NULL) {

  # Deduce a base name from the markers_df variable (e.g. SplitDotPlot(o, my_markers)
  # -> "my_markers"); only used when file_name is not supplied.
  df_name <- deparse(substitute(markers_df))
  if (!.usable_obj_name(df_name)) df_name <- NULL

  # Deprecated short aliases
  if (!is.null(gene_col)       && identical(gene_column,       "Genes"))    gene_column       <- gene_col
  if (!is.null(gene_group_col) && identical(gene_group_column, "CellType")) gene_group_column <- gene_group_col

  # Walk up output_dir from PrepObject when not explicitly supplied
  output_dir <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL

  # 1. Filter to genes present in the object
  valid_genes <- markers_df[[gene_column]][markers_df[[gene_column]] %in%
                                          rownames(seurat_object)]
  markers_df  <- markers_df[markers_df[[gene_column]] %in% valid_genes, ]
  if (!is.factor(markers_df[[gene_group_column]]))
    markers_df[[gene_group_column]] <- factor(markers_df[[gene_group_column]])

  # 2. Build interaction grouping variable
  has_split        <- !is.null(split.by) && nzchar(split.by)
  interaction_var  <- "Temp_Interaction_Var"
  on.exit(.make_cleanup(seurat_object, interaction_var)(), add = TRUE)

  seurat_object[[interaction_var]] <- if (has_split) {
    paste(seurat_object@meta.data[[group.by]], seurat_object@meta.data[[split.by]], sep = "___")
  } else {
    seurat_object@meta.data[[group.by]]
  }

  # 3. Extract DotPlot data via Seurat
  p_base    <- Seurat::DotPlot(seurat_object, features = unique(valid_genes),
                               group.by = interaction_var, scale = TRUE)
  plot_data <- p_base$data

  # 4. Split interaction ID back into group + row-split
  if (has_split) {
    split_ids              <- strsplit(as.character(plot_data$id), "___")
    plot_data$Y_Axis_Group <- sapply(split_ids, `[`, 1)
    plot_data$Row_Split    <- sapply(split_ids, `[`, 2)
  } else {
    plot_data$Y_Axis_Group <- as.character(plot_data$id)
    plot_data$Row_Split    <- "All"
  }

  # 5. Build numeric ID map and apply optional hierarchical clustering
  key        <- .build_id_map(plot_data, seurat_object@meta.data[[group.by]])
  id_map     <- key$id_map
  plot_data$Y_Index <- id_map[as.character(plot_data$Y_Axis_Group)]

  if (ClusterFeatures) {
    mat <- plot_data |>
      dplyr::select(Y_Axis_Group, features.plot, avg.exp.scaled) |>
      tidyr::pivot_wider(names_from  = features.plot,
                         values_from = avg.exp.scaled,
                         values_fn   = mean) |>
      tibble::column_to_rownames("Y_Axis_Group")
    mat[is.na(mat)] <- 0
    hc           <- stats::hclust(stats::dist(mat), method = "ward.D2")
    final_levels <- rownames(mat)[hc$order]
  } else {
    final_levels <- key$sorted
  }

  plot_data$Y_Axis_Group <- factor(plot_data$Y_Axis_Group, levels = final_levels)
  plot_data$Y_Label_Full <- paste0(id_map[final_levels][
    match(as.character(plot_data$Y_Axis_Group), final_levels)
  ], " | ", plot_data$Y_Axis_Group)
  label_order            <- paste0(id_map[final_levels], " | ", final_levels)
  plot_data$Y_Label_Full <- factor(plot_data$Y_Label_Full, levels = label_order)

  # 6. Clamp expression and percent values
  plot_data$avg.exp.scaled <- .clamp(plot_data$avg.exp.scaled, col.min, col.max)
  plot_data$pct.exp        <- .clamp(plot_data$pct.exp, scale.min, scale.max)

  # 7. Merge gene metadata
  colnames(markers_df)[colnames(markers_df) == gene_column]       <- "features.plot"
  colnames(markers_df)[colnames(markers_df) == gene_group_column] <- "Gene_Group"
  plot_data <- dplyr::left_join(
    plot_data,
    markers_df[, c("features.plot", "Gene_Group")],
    by = "features.plot"
  )

  # 8. Suggest PDF dimensions
  n_groups  <- length(unique(plot_data$Y_Axis_Group))
  n_splits  <- length(unique(plot_data$Row_Split))
  n_genes   <- length(unique(plot_data$features.plot))
  n_panels  <- length(unique(plot_data$Gene_Group))
  sug_h     <- 2 + (n_groups * n_splits * 0.2)
  sug_w     <- 3 + (n_genes * 0.2) + (n_panels * 0.15)
  message(sprintf("[SplitDotPlot] Suggested: pdf(file, width=%.1f, height=%.1f)",
                  sug_w, sug_h))

  # 9. Mini-index labels at the right margin
  text_data <- unique(plot_data[, c("Y_Label_Full", "Y_Index", "Row_Split")])

  # 10. Build plot
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = features.plot, y = Y_Label_Full)) +
    ggplot2::geom_point(ggplot2::aes(size = pct.exp, color = avg.exp.scaled)) +
    ggplot2::scale_color_viridis_c(option = "plasma", begin = 0, end = 1) +
    ggplot2::scale_size(range = c(0, dot.scale),
                        limits = c(scale.min, scale.max)) +
    ggplot2::facet_grid(rows  = dplyr::vars(Row_Split),
                        cols  = dplyr::vars(Gene_Group),
                        scales = "free", space = "free") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL,
                  color = "Avg Exp", size = "% Expressed") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::geom_text(
      data        = text_data,
      ggplot2::aes(y = Y_Label_Full, label = Y_Index),
      x           = Inf, hjust = 0, nudge_x = 0.2,
      size        = 3, color = "darkgray", inherit.aes = FALSE
    ) +
    .dotplot_theme()

  # ── Save or return ──────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fpath <- .dotplot_path(output_dir, seurat_object, file_name,
                           auto_base = df_name %||% group.by,
                           suffix    = "SplitDotPlot")
    grDevices::pdf(fpath, width = sug_w, height = sug_h)
    print(p)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath)
    return(invisible(p))
  }

  p
}

# --------------------------------------------------------------------------- #
# SplitDotPlot2                                                                #
# --------------------------------------------------------------------------- #

#' Dot plot with both group and gene-axis hierarchical clustering
#'
#' Extends [SplitDotPlot()] with two additional features:
#' \itemize{
#'   \item **`ClusterGenes`**: hierarchically clusters genes on the x-axis.
#'   \item **`min.pct.exp`**: drops genes expressed in fewer than this percent
#'         of cells in every group.
#' }
#'
#' @inheritParams SplitDotPlot
#' @param ClusterGroups Logical. Hierarchically cluster cell groups (y-axis)?
#'   Default `FALSE`.
#' @param ClusterGenes Logical. Hierarchically cluster genes (x-axis)?
#'   Default `TRUE`.
#' @param min.pct.exp Numeric. Minimum percent expression a gene must reach
#'   in at least one group to be kept. Default `15`.
#'
#' @return A ggplot2 object (invisibly when saved to disk).
#' @export
SplitDotPlot2 <- function(seurat_object,
                           markers_df,
                           gene_column        = "Genes",
                           gene_group_column  = "CellType",
                           group.by           = "seurat_clusters",
                           split.by           = NULL,
                           gene_col           = NULL,
                           gene_group_col     = NULL,
                           ClusterGroups   = FALSE,
                           ClusterGenes    = TRUE,
                           min.pct.exp     = 15,
                           col.min         = 0,
                           col.max         = 2.5,
                           scale.min       = 5,
                           scale.max       = 100,
                           dot.scale       = 6,
                           cols            = "RdBu",
                           output_dir      = NULL,
                           file_name       = NULL) {

  # Deduce a base name from the markers_df variable when file_name is absent
  df_name <- deparse(substitute(markers_df))
  if (!.usable_obj_name(df_name)) df_name <- NULL

  # Walk up output_dir from PrepObject when not explicitly supplied
  output_dir <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL

  valid_genes <- markers_df[[gene_column]][markers_df[[gene_column]] %in%
                                          rownames(seurat_object)]
  markers_df  <- markers_df[markers_df[[gene_column]] %in% valid_genes, ]
  if (!is.factor(markers_df[[gene_group_column]]))
    markers_df[[gene_group_column]] <- factor(markers_df[[gene_group_column]])

  has_split       <- !is.null(split.by) && nzchar(split.by)
  interaction_var <- "Temp_Interaction_Var"
  on.exit(.make_cleanup(seurat_object, interaction_var)(), add = TRUE)

  seurat_object[[interaction_var]] <- if (has_split) {
    paste(seurat_object@meta.data[[group.by]], seurat_object@meta.data[[split.by]], sep = "___")
  } else {
    seurat_object@meta.data[[group.by]]
  }

  p_base    <- Seurat::DotPlot(seurat_object, features = unique(valid_genes),
                               group.by = interaction_var, scale = TRUE)
  plot_data <- p_base$data

  # Filter low-expressed genes
  if (min.pct.exp > 0) {
    keep <- plot_data |>
      dplyr::group_by(features.plot) |>
      dplyr::summarise(max_pct = max(pct.exp), .groups = "drop") |>
      dplyr::filter(max_pct >= min.pct.exp) |>
      dplyr::pull(features.plot) |>
      as.character()
    plot_data   <- dplyr::filter(plot_data, features.plot %in% keep)
    valid_genes <- valid_genes[valid_genes %in% keep]
    if (length(valid_genes) == 0)
      stop("No genes passed the min.pct.exp threshold of ", min.pct.exp, "%.")
  }

  if (has_split) {
    split_ids              <- strsplit(as.character(plot_data$id), "___")
    plot_data$Y_Axis_Group <- sapply(split_ids, `[`, 1)
    plot_data$Row_Split    <- sapply(split_ids, `[`, 2)
  } else {
    plot_data$Y_Axis_Group <- as.character(plot_data$id)
    plot_data$Row_Split    <- "All"
  }

  key       <- .build_id_map(plot_data, seurat_object@meta.data[[group.by]])
  id_map    <- key$id_map
  plot_data$Y_Index <- id_map[as.character(plot_data$Y_Axis_Group)]

  # Y-axis ordering
  if (ClusterGroups) {
    mat <- plot_data |>
      dplyr::select(Y_Axis_Group, features.plot, avg.exp.scaled) |>
      tidyr::pivot_wider(names_from = features.plot,
                         values_from = avg.exp.scaled, values_fn = mean) |>
      tibble::column_to_rownames("Y_Axis_Group")
    mat[is.na(mat)] <- 0
    hc           <- stats::hclust(stats::dist(mat), method = "ward.D2")
    final_levels <- rownames(mat)[hc$order]
  } else {
    final_levels <- key$sorted
  }

  plot_data$Y_Axis_Group <- factor(plot_data$Y_Axis_Group, levels = final_levels)
  plot_data$Y_Label_Full <- paste0(id_map[final_levels][
    match(as.character(plot_data$Y_Axis_Group), final_levels)
  ], " | ", plot_data$Y_Axis_Group)
  label_order            <- paste0(id_map[final_levels], " | ", final_levels)
  plot_data$Y_Label_Full <- factor(plot_data$Y_Label_Full, levels = label_order)

  # X-axis ordering
  if (ClusterGenes) {
    gmat <- plot_data |>
      dplyr::select(id, features.plot, avg.exp.scaled) |>
      tidyr::pivot_wider(names_from = id,
                         values_from = avg.exp.scaled, values_fn = mean) |>
      tibble::column_to_rownames("features.plot")
    gmat[is.na(gmat)] <- 0
    ghc         <- stats::hclust(stats::dist(gmat), method = "ward.D2")
    gene_levels <- rownames(gmat)[ghc$order]
  } else {
    gene_levels <- unique(valid_genes)
  }
  plot_data$features.plot <- factor(plot_data$features.plot, levels = gene_levels)

  plot_data$avg.exp.scaled <- .clamp(plot_data$avg.exp.scaled, col.min, col.max)
  plot_data$pct.exp        <- .clamp(plot_data$pct.exp, scale.min, scale.max)

  colnames(markers_df)[colnames(markers_df) == gene_column]       <- "features.plot"
  colnames(markers_df)[colnames(markers_df) == gene_group_column] <- "Gene_Group"
  plot_data <- dplyr::left_join(
    plot_data,
    markers_df[, c("features.plot", "Gene_Group")],
    by = "features.plot"
  )

  n_groups <- length(unique(plot_data$Y_Axis_Group))
  n_splits <- length(unique(plot_data$Row_Split))
  n_genes  <- length(unique(plot_data$features.plot))
  n_panels <- length(unique(plot_data$Gene_Group))
  sug_h    <- 2 + (n_groups * n_splits * 0.2)
  sug_w    <- 3 + (n_genes * 0.2) + (n_panels * 0.15)
  message(sprintf("[SplitDotPlot2] Suggested: pdf(file, width=%.1f, height=%.1f)",
                  sug_w, sug_h))

  text_data <- unique(plot_data[, c("Y_Label_Full", "Y_Index", "Row_Split")])

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = features.plot, y = Y_Label_Full)) +
    ggplot2::geom_point(ggplot2::aes(size = pct.exp, color = avg.exp.scaled)) +
    ggplot2::scale_color_viridis_c(option = "plasma", begin = 0, end = 1) +
    ggplot2::scale_size(range = c(0, dot.scale),
                        limits = c(scale.min, scale.max)) +
    ggplot2::facet_grid(rows = dplyr::vars(Row_Split),
                        cols = dplyr::vars(Gene_Group),
                        scales = "free", space = "free") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL,
                  color = "Avg Exp", size = "% Expressed") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::geom_text(
      data = text_data,
      ggplot2::aes(y = Y_Label_Full, label = Y_Index),
      x = Inf, hjust = 0, nudge_x = 0.2,
      size = 3, color = "darkgray", inherit.aes = FALSE
    ) +
    .dotplot_theme()

  # ── Save or return ──────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fpath <- .dotplot_path(output_dir, seurat_object, file_name,
                           auto_base = df_name %||% group.by,
                           suffix    = "SplitDotPlot2")
    grDevices::pdf(fpath, width = sug_w, height = sug_h)
    print(p)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath)
    return(invisible(p))
  }

  p
}

# --------------------------------------------------------------------------- #
# FastDotPlot                                                                  #
# --------------------------------------------------------------------------- #

#' Fast dot plot with regex gene selection and k-group gene clustering
#'
#' Like [SplitDotPlot2()] but accepts genes via a character vector **or** a
#' regex pattern, and can slice the gene dendrogram into `k_genes` labelled
#' panels using `cutree()`.
#'
#' @param seurat_object A Seurat object.
#' @param assay Character. Assay to use. Default `"RNA"`.
#' @param features Character vector of gene names, or `NULL` if using
#'   `pattern`.
#' @param pattern Character regex. All matching genes in `rownames(seurat_object)` are
#'   used. Ignored if `features` is provided.
#' @param group.by Character. Metadata column for y-axis. Default
#'   `"GlobalAssignment"`.
#' @param split.by Character or `NULL`. Metadata column for row facets.
#' @param ClusterGroups Logical. Cluster y-axis groups? Default `FALSE`.
#' @param ClusterGenes Logical. Cluster genes on x-axis? Default `TRUE`.
#' @param k_genes Integer. Number of gene clusters (panels) to cut the
#'   dendrogram into. `1` = no cutting (single panel). Default `1`.
#' @param min.pct.exp Numeric. Minimum percent expression threshold.
#'   Default `15`.
#' @param col.min,col.max Numeric. Color clamp limits. Defaults `0`, `2.5`.
#' @param scale.min,scale.max Numeric. Size clamp limits. Defaults `5`, `100`.
#' @param dot.scale Numeric. Max dot radius. Default `6`.
#'
#' @inheritParams SplitDotPlot
#' @param file_name Character or \code{NULL}. Base name (no extension) for the
#'   saved PDF. \code{NULL} (default) auto-deduces from the \code{pattern}
#'   (when used), falling back to \code{object_name_group.by_FastDotPlot}.
#' @return A ggplot2 object (invisibly when saved to disk).
#' @export
FastDotPlot <- function(seurat_object,
                         assay        = "RNA",
                         features     = NULL,
                         pattern      = NULL,
                         group.by     = "seurat_clusters",
                         split.by     = NULL,
                         ClusterGroups = FALSE,
                         ClusterGenes  = TRUE,
                         k_genes       = 1,
                         min.pct.exp   = 15,
                         col.min       = 0,
                         col.max       = 2.5,
                         scale.min     = 5,
                         scale.max     = 100,
                         dot.scale     = 6,
                         output_dir    = NULL,
                         file_name     = NULL) {

  # Deduce a base name from the regex pattern when file_name is absent
  auto_base <- if (!is.null(pattern)) paste0("pattern_", pattern) else group.by

  # Walk up output_dir from PrepObject when not explicitly supplied
  output_dir <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL

  Seurat::DefaultAssay(seurat_object) <- assay

  valid_genes <- if (!is.null(pattern)) {
    grep(pattern, rownames(seurat_object), value = TRUE)
  } else if (!is.null(features)) {
    features[features %in% rownames(seurat_object)]
  } else {
    stop("Provide either 'features' or 'pattern'.")
  }
  if (length(valid_genes) == 0) stop("No valid genes found.")

  has_split       <- !is.null(split.by) && nzchar(split.by)
  interaction_var <- "Temp_Interaction_Var"
  on.exit(.make_cleanup(seurat_object, interaction_var)(), add = TRUE)

  seurat_object[[interaction_var]] <- if (has_split) {
    paste(seurat_object@meta.data[[group.by]], seurat_object@meta.data[[split.by]], sep = "___")
  } else {
    seurat_object@meta.data[[group.by]]
  }

  p_base    <- Seurat::DotPlot(seurat_object, features = unique(valid_genes),
                               group.by = interaction_var, scale = TRUE)
  plot_data <- p_base$data

  if (min.pct.exp > 0) {
    keep <- plot_data |>
      dplyr::group_by(features.plot) |>
      dplyr::summarise(max_pct = max(pct.exp), .groups = "drop") |>
      dplyr::filter(max_pct >= min.pct.exp) |>
      dplyr::pull(features.plot) |>
      as.character()
    plot_data   <- dplyr::filter(plot_data, features.plot %in% keep)
    valid_genes <- valid_genes[valid_genes %in% keep]
    if (length(valid_genes) == 0)
      stop("No genes passed the min.pct.exp threshold of ", min.pct.exp, "%.")
  }

  if (has_split) {
    split_ids              <- strsplit(as.character(plot_data$id), "___")
    plot_data$Y_Axis_Group <- sapply(split_ids, `[`, 1)
    plot_data$Row_Split    <- sapply(split_ids, `[`, 2)
  } else {
    plot_data$Y_Axis_Group <- as.character(plot_data$id)
    plot_data$Row_Split    <- "All"
  }

  key       <- .build_id_map(plot_data, seurat_object@meta.data[[group.by]])
  id_map    <- key$id_map
  plot_data$Y_Index <- id_map[as.character(plot_data$Y_Axis_Group)]

  if (ClusterGroups) {
    mat <- plot_data |>
      dplyr::select(Y_Axis_Group, features.plot, avg.exp.scaled) |>
      tidyr::pivot_wider(names_from = features.plot,
                         values_from = avg.exp.scaled, values_fn = mean) |>
      tibble::column_to_rownames("Y_Axis_Group")
    mat[is.na(mat)] <- 0
    hc           <- stats::hclust(stats::dist(mat), method = "ward.D2")
    final_levels <- rownames(mat)[hc$order]
  } else {
    final_levels <- key$sorted
  }

  plot_data$Y_Axis_Group <- factor(plot_data$Y_Axis_Group, levels = final_levels)
  plot_data$Y_Label_Full <- paste0(id_map[final_levels][
    match(as.character(plot_data$Y_Axis_Group), final_levels)
  ], " | ", plot_data$Y_Axis_Group)
  label_order <- paste0(id_map[final_levels], " | ", final_levels)
  plot_data$Y_Label_Full <- factor(plot_data$Y_Label_Full, levels = label_order)

  # Gene axis ordering and optional k-tree slicing
  if (ClusterGenes) {
    gmat <- plot_data |>
      dplyr::select(id, features.plot, avg.exp.scaled) |>
      tidyr::pivot_wider(names_from = id,
                         values_from = avg.exp.scaled, values_fn = mean) |>
      tibble::column_to_rownames("features.plot")
    gmat[is.na(gmat)] <- 0
    ghc         <- stats::hclust(stats::dist(gmat), method = "ward.D2")
    gene_levels <- rownames(gmat)[ghc$order]

    if (k_genes > 1) {
      gene_clusters   <- stats::cutree(ghc, k = k_genes)
      ordered_clusters <- unique(gene_clusters[gene_levels])
      remap           <- stats::setNames(seq_len(k_genes), ordered_clusters)
      cluster_map     <- data.frame(
        features.plot = names(gene_clusters),
        Gene_Cluster  = factor(
          paste0("Pattern ", remap[as.character(gene_clusters)]),
          levels = paste0("Pattern ", seq_len(k_genes))
        )
      )
      plot_data <- dplyr::left_join(plot_data, cluster_map,
                                    by = "features.plot")
    } else {
      plot_data$Gene_Cluster <- factor("All")
    }
  } else {
    gene_levels            <- unique(valid_genes)
    plot_data$Gene_Cluster <- factor("All")
  }
  plot_data$features.plot <- factor(plot_data$features.plot, levels = gene_levels)

  plot_data$avg.exp.scaled <- .clamp(plot_data$avg.exp.scaled, col.min, col.max)
  plot_data$pct.exp        <- .clamp(plot_data$pct.exp, scale.min, scale.max)

  n_groups <- length(unique(plot_data$Y_Axis_Group))
  n_splits <- length(unique(plot_data$Row_Split))
  n_genes  <- length(unique(plot_data$features.plot))
  sug_h    <- 2 + (n_groups * n_splits * 0.2)
  sug_w    <- 3 + (n_genes * 0.2) + (ifelse(k_genes > 1, k_genes * 0.15, 0))
  message(sprintf("[FastDotPlot] Suggested: pdf(file, width=%.1f, height=%.1f)",
                  sug_w, sug_h))

  text_data        <- unique(plot_data[, c("Y_Label_Full", "Y_Index", "Row_Split")])
  has_gene_clusters <- ClusterGenes && k_genes > 1

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = features.plot, y = Y_Label_Full)) +
    ggplot2::geom_point(ggplot2::aes(size = pct.exp, color = avg.exp.scaled)) +
    ggplot2::scale_color_viridis_c(option = "plasma", begin = 0, end = 1) +
    ggplot2::scale_size(range = c(0, dot.scale),
                        limits = c(scale.min, scale.max)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL,
                  color = "Avg Exp", size = "% Expressed") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::geom_text(
      data = text_data,
      ggplot2::aes(y = Y_Label_Full, label = Y_Index),
      x = Inf, hjust = 0, nudge_x = 0.2,
      size = 3, color = "darkgray", inherit.aes = FALSE
    ) +
    .dotplot_theme()

  if (has_split && has_gene_clusters) {
    p <- p + ggplot2::facet_grid(rows = dplyr::vars(Row_Split),
                                 cols = dplyr::vars(Gene_Cluster),
                                 scales = "free", space = "free")
  } else if (has_split) {
    p <- p + ggplot2::facet_grid(rows = dplyr::vars(Row_Split),
                                 scales = "free_y", space = "free_y")
  } else if (has_gene_clusters) {
    p <- p + ggplot2::facet_grid(cols = dplyr::vars(Gene_Cluster),
                                 scales = "free_x", space = "free_x")
  }

  # ── Save or return ──────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fpath <- .dotplot_path(output_dir, seurat_object, file_name,
                           auto_base = auto_base,
                           suffix    = "FastDotPlot")
    grDevices::pdf(fpath, width = sug_w, height = sug_h)
    print(p)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath)
    return(invisible(p))
  }

  p
}

# --------------------------------------------------------------------------- #
# FastDotPlot2                                                                 #
# --------------------------------------------------------------------------- #

#' FastDotPlot with optional EnrichR enrichment header
#'
#' Extends [FastDotPlot()] with an optional enrichment annotation panel
#' aligned above the dot plot columns. When `RunEnrichR = TRUE`, the top
#' `TopN_Enrich` significant terms from the specified EnrichR database are
#' displayed as a text header above each gene-cluster panel.
#'
#' **Bug fixes vs. original script:**
#' \itemize{
#'   \item Gene-cluster strip labels are shown even when enrichment is
#'         disabled or fails (they were previously always hidden).
#'   \item The filtered gene list and per-pattern gene assignments are always
#'         returned in the output list regardless of enrichment status.
#' }
#'
#' @inheritParams FastDotPlot
#' @param RunEnrichR Logical. Run EnrichR enrichment? Default `TRUE`.
#' @param EnrichR_DB Character. EnrichR database name. Default
#'   `"CellMarker_2024"`.
#' @param TopN_Enrich Integer. Number of top terms to show per pattern.
#'   Default `3`.
#'
#' @return A named list:
#' \describe{
#'   \item{plot}{The combined (enrichment + dot plot) or standalone dot plot.}
#'   \item{dotplot}{The standalone dot plot ggplot2 object.}
#'   \item{enrich_plot}{The enrichment header ggplot2 object, or `NULL`.}
#'   \item{patterns}{Named list mapping each pattern label to its gene vector.}
#'   \item{filtered_genes}{Data frame with columns `features.plot` and
#'         `Gene_Cluster` for every gene that passed `min.pct.exp`.}
#'   \item{enrich_table}{Raw EnrichR results data frame, or `NULL`.}
#' }
#' @export
FastDotPlot2 <- function(seurat_object,
                          assay         = "RNA",
                          features      = NULL,
                          pattern       = NULL,
                          group.by      = "seurat_clusters",
                          split.by      = NULL,
                          ClusterGroups  = FALSE,
                          ClusterGenes   = TRUE,
                          k_genes        = 1,
                          min.pct.exp    = 15,
                          col.min        = 0,
                          col.max        = 2.5,
                          scale.min      = 5,
                          scale.max      = 100,
                          dot.scale      = 6,
                          RunEnrichR     = TRUE,
                          EnrichR_DB     = "CellMarker_2024",
                          TopN_Enrich    = 3,
                          output_dir     = NULL,
                          file_name      = NULL) {

  # Deduce a base name from the regex pattern when file_name is absent
  auto_base <- if (!is.null(pattern)) paste0("pattern_", pattern) else group.by

  # Walk up output_dir from PrepObject when not explicitly supplied
  output_dir <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL

  Seurat::DefaultAssay(seurat_object) <- assay

  valid_genes <- if (!is.null(pattern)) {
    grep(pattern, rownames(seurat_object), value = TRUE)
  } else if (!is.null(features)) {
    features[features %in% rownames(seurat_object)]
  } else {
    stop("Provide either 'features' or 'pattern'.")
  }
  if (length(valid_genes) == 0) stop("No valid genes found.")

  has_split       <- !is.null(split.by) && nzchar(split.by)
  interaction_var <- "Temp_Interaction_Var"
  on.exit(.make_cleanup(seurat_object, interaction_var)(), add = TRUE)

  seurat_object[[interaction_var]] <- if (has_split) {
    paste(seurat_object@meta.data[[group.by]], seurat_object@meta.data[[split.by]], sep = "___")
  } else {
    seurat_object@meta.data[[group.by]]
  }

  p_base    <- Seurat::DotPlot(seurat_object, features = unique(valid_genes),
                               group.by = interaction_var, scale = TRUE)
  plot_data <- p_base$data

  if (min.pct.exp > 0) {
    keep <- plot_data |>
      dplyr::group_by(features.plot) |>
      dplyr::summarise(max_pct = max(pct.exp), .groups = "drop") |>
      dplyr::filter(max_pct >= min.pct.exp) |>
      dplyr::pull(features.plot) |>
      as.character()
    plot_data   <- dplyr::filter(plot_data, features.plot %in% keep)
    valid_genes <- valid_genes[valid_genes %in% keep]
    if (length(valid_genes) == 0)
      stop("No genes passed the min.pct.exp threshold of ", min.pct.exp, "%.")
  }

  if (has_split) {
    split_ids              <- strsplit(as.character(plot_data$id), "___")
    plot_data$Y_Axis_Group <- sapply(split_ids, `[`, 1)
    plot_data$Row_Split    <- sapply(split_ids, `[`, 2)
  } else {
    plot_data$Y_Axis_Group <- as.character(plot_data$id)
    plot_data$Row_Split    <- "All"
  }

  key       <- .build_id_map(plot_data, seurat_object@meta.data[[group.by]])
  id_map    <- key$id_map
  plot_data$Y_Index <- id_map[as.character(plot_data$Y_Axis_Group)]

  if (ClusterGroups) {
    mat <- plot_data |>
      dplyr::select(Y_Axis_Group, features.plot, avg.exp.scaled) |>
      tidyr::pivot_wider(names_from = features.plot,
                         values_from = avg.exp.scaled, values_fn = mean) |>
      tibble::column_to_rownames("Y_Axis_Group")
    mat[is.na(mat)] <- 0
    hc           <- stats::hclust(stats::dist(mat), method = "ward.D2")
    final_levels <- rownames(mat)[hc$order]
  } else {
    final_levels <- key$sorted
  }

  plot_data$Y_Axis_Group <- factor(plot_data$Y_Axis_Group, levels = final_levels)
  plot_data$Y_Label_Full <- paste0(id_map[final_levels][
    match(as.character(plot_data$Y_Axis_Group), final_levels)
  ], " | ", plot_data$Y_Axis_Group)
  label_order <- paste0(id_map[final_levels], " | ", final_levels)
  plot_data$Y_Label_Full <- factor(plot_data$Y_Label_Full, levels = label_order)

  # Gene ordering and k-tree slicing
  pattern_list <- list()
  if (ClusterGenes) {
    gmat <- plot_data |>
      dplyr::select(id, features.plot, avg.exp.scaled) |>
      tidyr::pivot_wider(names_from = id,
                         values_from = avg.exp.scaled, values_fn = mean) |>
      tibble::column_to_rownames("features.plot")
    gmat[is.na(gmat)] <- 0
    ghc         <- stats::hclust(stats::dist(gmat), method = "ward.D2")
    gene_levels <- rownames(gmat)[ghc$order]

    if (k_genes > 1) {
      gene_clusters    <- stats::cutree(ghc, k = k_genes)
      ordered_clusters <- unique(gene_clusters[gene_levels])
      remap            <- stats::setNames(seq_len(k_genes), ordered_clusters)
      cluster_map      <- data.frame(
        features.plot = names(gene_clusters),
        Gene_Cluster  = factor(
          paste0("Pattern ", remap[as.character(gene_clusters)]),
          levels = paste0("Pattern ", seq_len(k_genes))
        )
      )
      plot_data    <- dplyr::left_join(plot_data, cluster_map,
                                       by = "features.plot")
      raw_split    <- split(as.character(cluster_map$features.plot),
                            cluster_map$Gene_Cluster)
      pattern_list <- lapply(raw_split, function(g)
        g[order(match(g, gene_levels))])
    } else {
      plot_data$Gene_Cluster <- factor("All")
      pattern_list           <- list(All = gene_levels)
    }
  } else {
    gene_levels            <- unique(valid_genes)
    plot_data$Gene_Cluster <- factor("All")
    pattern_list           <- list(Unclustered_Genes = gene_levels)
  }
  plot_data$features.plot <- factor(plot_data$features.plot, levels = gene_levels)

  # Always export filtered gene table regardless of enrichment
  filtered_genes <- unique(plot_data[, c("features.plot", "Gene_Cluster")])

  plot_data$avg.exp.scaled <- .clamp(plot_data$avg.exp.scaled, col.min, col.max)
  plot_data$pct.exp        <- .clamp(plot_data$pct.exp, scale.min, scale.max)

  n_groups <- length(unique(plot_data$Y_Axis_Group))
  n_splits <- length(unique(plot_data$Row_Split))
  n_genes  <- length(unique(plot_data$features.plot))
  sug_h    <- 2 + (n_groups * n_splits * 0.2)
  sug_w    <- 3 + (n_genes * 0.2) + (ifelse(k_genes > 1, k_genes * 0.15, 0))
  message(sprintf("[FastDotPlot2] Suggested: pdf(file, width=%.1f, height=%.1f)",
                  sug_w, sug_h))

  text_data         <- unique(plot_data[, c("Y_Label_Full", "Y_Index", "Row_Split")])
  has_gene_clusters <- ClusterGenes && k_genes > 1

  # Build dot plot - strip labels shown when no enrichment header is present;
  # hidden only when the enrichment header will cover them.
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = features.plot, y = Y_Label_Full)) +
    ggplot2::geom_point(ggplot2::aes(size = pct.exp, color = avg.exp.scaled)) +
    ggplot2::scale_color_viridis_c(option = "plasma", begin = 0, end = 1) +
    ggplot2::scale_size(range = c(0, dot.scale),
                        limits = c(scale.min, scale.max)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL,
                  color = "Avg Exp", size = "% Expressed") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::geom_text(
      data = text_data,
      ggplot2::aes(y = Y_Label_Full, label = Y_Index),
      x = Inf, hjust = 0, nudge_x = 0.2,
      size = 3, color = "darkgray", inherit.aes = FALSE
    ) +
    .dotplot_theme()

  if (has_split && has_gene_clusters) {
    p <- p + ggplot2::facet_grid(rows = dplyr::vars(Row_Split),
                                 cols = dplyr::vars(Gene_Cluster),
                                 scales = "free", space = "free")
  } else if (has_split) {
    p <- p + ggplot2::facet_grid(rows = dplyr::vars(Row_Split),
                                 scales = "free_y", space = "free_y")
  } else if (has_gene_clusters) {
    p <- p + ggplot2::facet_grid(cols = dplyr::vars(Gene_Cluster),
                                 scales = "free_x", space = "free_x")
  }

  # EnrichR enrichment header
  enrich_res_df <- NULL
  enrich_plot   <- NULL
  combined_plot <- p   # default: standalone dot plot

  if (RunEnrichR && length(pattern_list) > 0) {
    if (!requireNamespace("enrichR", quietly = TRUE))
      stop("Package 'enrichR' is required for RunEnrichR = TRUE.")

    message("Running EnrichR (", EnrichR_DB, ") for ", length(pattern_list),
            " gene pattern(s)...")

    enrich_list <- lapply(names(pattern_list), function(pat_name) {
      genes <- pattern_list[[pat_name]]
      if (length(genes) == 0) return(NULL)
      res <- enrichR::enrichr(genes, databases = EnrichR_DB)
      if (!is.null(res[[EnrichR_DB]]) && nrow(res[[EnrichR_DB]]) > 0) {
        df          <- res[[EnrichR_DB]]
        df$Pattern  <- pat_name
        return(df)
      }
      NULL
    })
    enrich_res_df <- dplyr::bind_rows(enrich_list)

    if (!is.null(enrich_res_df) && nrow(enrich_res_df) > 0) {

      # Midpoints for centering text within each facet panel
      pattern_widths <- plot_data |>
        dplyr::group_by(Gene_Cluster) |>
        dplyr::summarise(n_genes = dplyr::n_distinct(features.plot),
                         .groups = "drop") |>
        dplyr::mutate(mid_x = (n_genes + 1) / 2)

      plot_df <- enrich_res_df |>
        dplyr::mutate(Log_Combined_Score = log(Combined.Score)) |>
        dplyr::filter(Adjusted.P.value < 0.05) |>
        dplyr::group_by(Pattern) |>
        dplyr::arrange(dplyr::desc(Log_Combined_Score)) |>
        dplyr::slice_head(n = TopN_Enrich) |>
        dplyr::mutate(
          Rank       = dplyr::row_number(),
          Term_Clean = stringr::str_trunc(Term, width = 35),
          Label      = paste0(Rank, ". ", Term_Clean,
                              " (", round(Log_Combined_Score, 1), ")"),
          Gene_Cluster = factor(Pattern,
                                levels = levels(plot_data$Gene_Cluster))
        ) |>
        dplyr::ungroup() |>
        dplyr::left_join(pattern_widths, by = "Gene_Cluster")

      if (nrow(plot_df) > 0) {
        # When the enrichment header is drawn, hide the dotplot strip labels
        # to avoid double-labelling the same column group
        p <- p + ggplot2::theme(
          strip.background.x = ggplot2::element_blank(),
          strip.text.x       = ggplot2::element_blank()
        )

        enrich_plot <- ggplot2::ggplot() +
          ggplot2::geom_blank(data = plot_data,
                              ggplot2::aes(x = features.plot, y = -1)) +
          ggplot2::geom_text(
            data      = plot_df,
            ggplot2::aes(x = mid_x, y = -Rank, label = Label,
                         color = Gene_Cluster),
            hjust     = 0.5, size = 3.5, fontface = "bold",
            show.legend = FALSE
          ) +
          ggplot2::facet_grid(cols  = dplyr::vars(Gene_Cluster),
                              scales = "free_x", space = "free_x") +
          ggplot2::scale_color_viridis_d(option = "plasma",
                                         begin = 0.2, end = 0.8) +
          ggplot2::scale_y_continuous(
            limits = c(-TopN_Enrich - 0.5, -0.5)
          ) +
          ggplot2::labs(title = EnrichR_DB) +
          ggplot2::theme_void() +
          ggplot2::theme(
            plot.title   = ggplot2::element_text(face = "bold", size = 16,
                                                 hjust = 0.5,
                                                 margin = ggplot2::margin(b = 10)),
            strip.background = ggplot2::element_rect(fill = "grey95",
                                                     colour = "grey90",
                                                     linewidth = 1),
            strip.text   = ggplot2::element_text(face = "bold", size = 11,
                                                 margin = ggplot2::margin(t = 6, b = 6)),
            panel.border = ggplot2::element_rect(fill = NA, colour = "grey90",
                                                 linewidth = 1),
            panel.spacing = ggplot2::unit(1, "lines"),
            plot.margin  = ggplot2::margin(b = 2)
          )

        combined_plot <- enrich_plot / p +
          patchwork::plot_layout(heights = c(0.6, 3))
      }
    }
  }

  # ── Save or return ──────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fpath <- .dotplot_path(output_dir, seurat_object, file_name,
                           auto_base = auto_base,
                           suffix    = "FastDotPlot2")
    grDevices::pdf(fpath, width = sug_w, height = sug_h)
    print(combined_plot)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath)
  } else {
    print(combined_plot)
  }

  invisible(list(
    plot           = combined_plot,
    dotplot        = p,
    enrich_plot    = enrich_plot,
    patterns       = pattern_list,
    filtered_genes = filtered_genes,
    enrich_table   = enrich_res_df
  ))
}
