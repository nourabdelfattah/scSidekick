# =============================================================================
# scSidekick utility functions
#
# Includes the package-wide null-coalescing operator %||% used throughout.
#
# ExtractMethods — reads seurat_object@commands to auto-generate a methods paragraph
#                  summarising preprocessing, integration, clustering, and
#                  dimensionality reduction steps.
#
# InspectPlot    — inspects a ggplot2 object and returns a plain-text
#                  description of what was plotted (geom types, aesthetic
#                  mappings, faceting, scales). For scSidekick functions the
#                  description is attached as attr(plot, "legend_text").
# =============================================================================

# Package-wide null-coalescing operator (not exported — internal only)
# Usage: a %||% b  →  returns a if non-NULL, else b
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Extract analysis parameters and generate a methods text from a Seurat object
#'
#' Reads the command history stored in `seurat_object@commands` to reconstruct the
#' parameters used for normalization, variable feature selection, scaling,
#' PCA, integration (Harmony / FastMNN / Seurat CCA), UMAP, and clustering.
#' Returns both a structured summary list and a ready-to-paste methods
#' paragraph.
#'
#' @param seurat_object A Seurat object with commands recorded (i.e., the standard
#'   Seurat workflow functions were called directly on the object, not via
#'   custom wrappers that bypass command logging).
#' @param cite_seurat Logical. Append a Seurat citation reminder at the end
#'   of the methods paragraph. Default `TRUE`.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{`summary`}{Named list of extracted parameter values.}
#'   \item{`methods_text`}{Character string — a draft methods paragraph
#'     ready to paste into a manuscript.}
#' }
#' @export
ExtractMethods <- function(seurat_object, cite_seurat = TRUE) {

  cmds  <- names(seurat_object@commands)
  get_p <- function(cmd, param) {
    if (cmd %in% cmds) seurat_object@commands[[cmd]]@params[[param]] else NULL
  }

  # ---- Basic counts ----
  n_cells <- ncol(seurat_object)
  n_genes <- nrow(seurat_object)
  assay   <- Seurat::DefaultAssay(seurat_object)

  # ---- Normalization ----
  norm_method  <- get_p("NormalizeData", "normalization.method") %||% "LogNormalize"
  scale_factor <- get_p("NormalizeData", "scale.factor")         %||% 10000

  # ---- Variable features ----
  n_var_feat   <- get_p("FindVariableFeatures", "nfeatures")     %||% 2000
  var_method   <- get_p("FindVariableFeatures", "selection.method") %||% "vst"

  # ---- Scaling / regression ----
  regress_vars <- get_p("ScaleData", "vars.to.regress")

  # ---- PCA ----
  n_pcs_computed <- get_p("RunPCA", "npcs") %||%
    if ("pca" %in% names(seurat_object@reductions)) ncol(seurat_object@reductions$pca) else NA

  # ---- Integration ----
  integration_method <- NULL
  integration_params <- list()
  if ("RunHarmony" %in% cmds) {
    integration_method  <- "Harmony"
    integration_params  <- list(
      group_by = get_p("RunHarmony", "group.by.vars"),
      theta    = get_p("RunHarmony", "theta")
    )
  } else if ("IntegrateLayers" %in% cmds) {
    integration_method  <- "Seurat v5 IntegrateLayers"
    integration_params  <- list(
      method = get_p("IntegrateLayers", "method")
    )
  } else if ("RunFastMNN" %in% cmds) {
    integration_method  <- "FastMNN"
  }

  # ---- UMAP ----
  umap_dims     <- get_p("RunUMAP", "dims")
  umap_nn       <- get_p("RunUMAP", "n.neighbors") %||% 30
  umap_dist     <- get_p("RunUMAP", "min.dist")    %||% 0.3
  umap_red_used <- get_p("RunUMAP", "reduction")   %||% "pca"

  n_dims <- if (!is.null(umap_dims)) max(umap_dims) else NA

  # ---- Clustering ----
  clust_res  <- get_p("FindClusters", "resolution")
  clust_alg  <- switch(as.character(get_p("FindClusters", "algorithm") %||% 1),
                       "1" = "Louvain",
                       "2" = "Louvain with multilevel refinement",
                       "3" = "SLM",
                       "4" = "Leiden",
                       "Louvain")
  n_clusters <- if ("seurat_clusters" %in% colnames(seurat_object@meta.data))
    length(levels(seurat_object@meta.data$seurat_clusters)) else NA

  # ---- Available reductions ----
  reductions <- names(seurat_object@reductions)

  # ---- Assemble structured summary ----
  summary_list <- list(
    n_cells            = n_cells,
    n_genes            = n_genes,
    assay              = assay,
    normalization      = list(method = norm_method, scale_factor = scale_factor),
    variable_features  = list(n = n_var_feat, method = var_method),
    regression_vars    = regress_vars,
    pca_npcs           = n_pcs_computed,
    integration        = list(method = integration_method,
                              params = integration_params),
    umap               = list(dims       = n_dims,
                              reduction  = umap_red_used,
                              n.neighbors = umap_nn,
                              min.dist   = umap_dist),
    clustering         = list(algorithm  = clust_alg,
                              resolution = clust_res,
                              n_clusters = n_clusters),
    reductions_present = reductions
  )

  # ---- Build methods paragraph ----
  lines <- character(0)

  lines <- c(lines, sprintf(
    "Single-cell RNA-seq data (%s cells, %s genes; assay: %s) were analyzed using Seurat.",
    format(n_cells, big.mark = ","),
    format(n_genes, big.mark = ","),
    assay
  ))

  lines <- c(lines, sprintf(
    "Gene expression counts were normalized using %s (scale factor %s) and %s highly variable genes were identified using the %s method.",
    norm_method, format(scale_factor, big.mark = ","),
    format(n_var_feat, big.mark = ","), var_method
  ))

  if (!is.null(regress_vars) && length(regress_vars) > 0) {
    lines <- c(lines, sprintf(
      "Data were scaled with regression of: %s.",
      paste(regress_vars, collapse = ", ")
    ))
  }

  if (!is.na(n_pcs_computed))
    lines <- c(lines, sprintf("PCA was performed (%d PCs computed).", n_pcs_computed))

  if (!is.null(integration_method)) {
    int_detail <- switch(integration_method,
      "Harmony" = sprintf(
        "Batch correction was performed using Harmony (group.by.vars = '%s'%s).",
        paste(integration_params$group_by, collapse = ", "),
        if (!is.null(integration_params$theta))
          paste0(", theta = ", integration_params$theta) else ""
      ),
      sprintf("Data integration was performed using %s.", integration_method)
    )
    lines <- c(lines, int_detail)
  }

  if (!is.na(n_dims))
    lines <- c(lines, sprintf(
      "UMAP was computed on the first %d %s components (n.neighbors = %s, min.dist = %s).",
      n_dims, umap_red_used, umap_nn, umap_dist
    ))

  if (!is.null(clust_res) && !is.na(n_clusters))
    lines <- c(lines, sprintf(
      "Cells were clustered using the %s algorithm at resolution %.2f, yielding %d clusters.",
      clust_alg, clust_res, n_clusters
    ))

  if (cite_seurat)
    lines <- c(lines,
               "(Cite: Hao et al., Cell 2021 / Hao et al., Nature Methods 2024 as appropriate.)")

  methods_text <- paste(lines, collapse = " ")

  list(summary = summary_list, methods_text = methods_text)
}

# Null-coalescing operator (like %||% in rlang, but without the dependency)
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Round a number up to the next "nice" value
#'
#' Given a positive number `x`, returns the smallest element of
#' `10^floor(log10(x)) * nice` that is \eqn{\geq x}.  Useful for setting
#' clean axis limits or bin widths.
#'
#' @param x A single positive numeric value.
#' @param nice Numeric vector of multipliers to cycle through within each
#'   power-of-ten decade.  Default `1:10`.
#'
#' @return A single numeric value rounded up to the next nice number.
#' @export
roundUpNice <- function(x, nice = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) {
  if (length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}


#' Inspect a ggplot2 object and describe what it shows
#'
#' Extracts metadata from a ggplot2 object — geom types, key aesthetic
#' mappings, axis labels, faceting, and color/fill scale names — and returns
#' a plain-text description suitable for use as a figure legend draft.
#'
#' For plots created by scSidekick functions, a pre-formatted legend is stored
#' as `attr(plot, "legend_text")` and is returned directly when present.
#'
#' @param p A ggplot2 object.
#' @param title Character or `NULL`. Optional figure panel title to include
#'   in the legend text (e.g., `"Figure 2A"`).
#'
#' @return A character string — a draft legend sentence.
#' @export
InspectPlot <- function(p, title = NULL) {

  # ---- scSidekick composite objects (ggarrange / gtable) ----
  # GenerateFeatureMaps, PlotGridUMAP etc. return assembled grids that are
  # not pure ggplot objects. Check for the attached legend_text first.
  stored <- attr(p, "legend_text")
  if (!is.null(stored)) {
    out <- if (!is.null(title)) paste0(title, ". ", stored) else stored
    message(out)
    return(invisible(out))
  }

  if (!inherits(p, "gg"))
    stop("'p' is not a ggplot2 object and has no 'legend_text' attribute.\n",
         "scSidekick composite plots (ggarrange grids) carry their description in ",
         "attr(p, 'legend_text') — make sure the function returned the plot ",
         "object (out_dir = NULL) and that the legend_text attribute was set.")

  # ---- Geom types ----
  # Seurat uses GeomDrawGrob / GeomScattermore internally — these are opaque
  # renderers and not informative on their own. Map them to a readable label.
  .OPAQUE_GEOMS <- c("DrawGrob", "Scattermore", "CustomAnn",
                     "SpatialImage", "Blank")

  geom_raw <- vapply(p$layers, function(l) {
    sub("^Geom", "", class(l$geom)[1])
  }, character(1))

  is_opaque <- all(geom_raw %in% .OPAQUE_GEOMS) || length(geom_raw) == 0
  geom_str  <- if (is_opaque) NULL else paste(unique(geom_raw), collapse = " + ")

  # ---- Labels (always reliable, even for Seurat plots) ----
  lbs       <- p$labels
  x_lab     <- lbs$x     %||% ""
  y_lab     <- lbs$y     %||% ""
  col_lab   <- lbs$colour %||% lbs$color %||% lbs$fill %||% ""
  plot_title <- lbs$title   %||% lbs$subtitle %||% ""

  # ---- Detect reduction type from axis labels ----
  reduction <- dplyr::case_when(
    grepl("UMAP",        x_lab, ignore.case = TRUE) ~ "UMAP",
    grepl("tSNE|TSNE",  x_lab, ignore.case = TRUE) ~ "tSNE",
    grepl("^PC[_\\s]?1", x_lab, ignore.case = TRUE) ~ "PCA",
    grepl("Harmony",    x_lab, ignore.case = TRUE) ~ "Harmony",
    TRUE ~ ""
  )

  # ---- Detect plot type from labels when geom is opaque ----
  # A FeaturePlot colours by a continuous gene/score; a DimPlot by a discrete
  # identity. Seurat puts the feature name in p$labels$colour for FeaturePlot
  # and "Identity" / the group.by column for DimPlot.
  colour_scale_type <- tryCatch({
    sc <- Filter(function(s) any(s$aesthetics %in% c("colour","color","fill")),
                 p$scales$scales)
    if (length(sc) == 0) "unknown"
    else if (inherits(sc[[1]], "ScaleContinuous")) "continuous"
    else "discrete"
  }, error = function(e) "unknown")

  seurat_plot_type <-
    if (nchar(reduction) > 0 && colour_scale_type == "continuous")
      "feature plot"
    else if (nchar(reduction) > 0 && colour_scale_type == "discrete")
      "cell identity plot"
    else if (nchar(reduction) > 0)
      "embedding plot"
    else NULL

  # ---- Colour / fill scale limits ----
  scale_detail <- tryCatch({
    sc <- Filter(function(s) any(s$aesthetics %in% c("colour","color","fill")),
                 p$scales$scales)
    if (length(sc) == 0 || colour_scale_type != "continuous") return(NULL)
    lims <- sc[[1]]$limits
    if (!is.null(lims) && !any(is.na(lims)))
      sprintf("colour range %.2g – %.2g", lims[1], lims[2])
    else NULL
  }, error = function(e) NULL)

  # ---- Faceting ----
  facet_class <- class(p$facet)[1]
  facet_str   <- switch(facet_class,
    FacetGrid = {
      rows <- names(p$facet$params$rows)
      cols <- names(p$facet$params$cols)
      parts <- c(
        if (length(rows)) paste("rows:", paste(rows, collapse = ", ")),
        if (length(cols)) paste("cols:",  paste(cols, collapse = ", "))
      )
      if (length(parts)) paste("split by facet_grid(", paste(parts, collapse = "; "), ")")
      else ""
    },
    FacetWrap = {
      fvars <- names(p$facet$params$facets)
      if (length(fvars))
        paste("split by", paste(fvars, collapse = " × "))
      else ""
    },
    ""
  )

  # ---- Assemble the description ----
  # Priority: use Seurat-aware phrasing when possible; fall back to raw geom info.
  if (!is.null(seurat_plot_type)) {
    # e.g. "UMAP feature plot of APOE, coloured by expression (range 0–3.2)"
    feature_part <- if (nchar(col_lab) > 0 && col_lab != "value")
      paste0("of ", col_lab) else ""
    panel_part   <- if (nchar(plot_title) > 0 && plot_title != col_lab)
      paste0(" [panel: ", plot_title, "]") else ""
    colour_part  <- if (colour_scale_type == "continuous")
      paste0("coloured by expression",
             if (!is.null(scale_detail)) paste0(" (", scale_detail, ")") else "")
    else if (nchar(col_lab) > 0)
      paste0("coloured by ", col_lab)
    else ""

    bits <- Filter(nzchar, c(
      if (!is.null(title))       title,
      paste(reduction, seurat_plot_type, feature_part),
      panel_part,
      colour_part,
      facet_str
    ))
    out <- paste(bits, collapse = ". ")

  } else {
    # Generic non-Seurat path
    aes_map <- p$mapping
    aes_str <- if (length(aes_map)) {
      paste(vapply(names(aes_map), function(nm) {
        val <- tryCatch(rlang::as_label(aes_map[[nm]]),
                        error = function(e) "?")
        paste0(nm, " = ", val)
      }, character(1)), collapse = ", ")
    } else ""

    plot_type <- if (!is.null(geom_str)) paste(geom_str, "plot") else "plot"

    bits <- Filter(nzchar, c(
      if (!is.null(title))           title,
      plot_type,
      if (nchar(plot_title) > 0)     paste("—", plot_title),
      if (nchar(aes_str) > 0)        paste("showing", aes_str),
      if (nchar(x_lab) > 0 && x_lab != "x")   paste("x:", x_lab),
      if (nchar(y_lab) > 0 && y_lab != "y")   paste("y:", y_lab),
      if (nchar(col_lab) > 0)        paste("coloured by", col_lab),
      if (!is.null(scale_detail))    scale_detail,
      facet_str
    ))
    out <- paste(bits, collapse = ". ")
  }

  if (!nzchar(trimws(out))) out <- "Plot (no descriptive metadata found)"
  out <- sub("\\.*$", ".", out)   # ensure single trailing full stop

  message(out)
  invisible(out)
}


# =============================================================================
# .write_subdir_params()
# Internal helper used by RunGSEA and RunCellChat to write an
# analysis_params.json into their own output directory.
#
# Strategy:
#   1. Walk up the directory tree from output_dir to find the main
#      analysis_params.json (the one written by log_analysis_params()).
#   2. Load it as the base (contains n_cells, resolution, dataset name, …).
#   3. Merge with method-specific fields (passed as extra_params).
#      Method-specific keys take precedence; shared keys such as n_cells are
#      kept from the parent so they appear in the PPTX methods slide.
#   4. Write the merged JSON to file.path(output_dir, "analysis_params.json").
#
# When create_analysis_pptx() is later called on this sub-directory, the
# walk-up finds the LOCAL file first and gets all fields (base + method).
# =============================================================================
.write_subdir_params <- function(output_dir, extra_params) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("jsonlite not available — skipping analysis_params.json write.")
    return(invisible(NULL))
  }

  # Walk up from the PARENT directory (not output_dir itself, which is the
  # sub-folder we are about to write into).
  parent_json <- .find_params_json(dirname(normalizePath(output_dir,
                                                          mustWork = FALSE)))
  base_params <- if (!is.null(parent_json) && file.exists(parent_json)) {
    jsonlite::read_json(parent_json, simplifyVector = TRUE)
  } else {
    list()
  }

  # Merge: start from base, then let extra_params add or override.
  # Keys that already exist in base_params are kept UNLESS extra_params
  # explicitly provides a new value (e.g. a method-specific methods_text
  # replaces the generic one from log_analysis_params).
  merged <- base_params
  for (nm in names(extra_params))
    merged[[nm]] <- extra_params[[nm]]

  out_path <- file.path(output_dir, "analysis_params.json")
  jsonlite::write_json(merged, out_path, auto_unbox = TRUE, pretty = TRUE)
  message("Method params written to: ", out_path)
  invisible(out_path)
}


# =============================================================================
# .get_layer_data()
# Version-agnostic Seurat expression data extractor.
#
# Handles three common failure modes transparently:
#   • Seurat v5  — uses SeuratObject::LayerData() with layer= parameter
#   • Seurat v3  — falls back to Seurat::GetAssayData() with slot= parameter
#   • BPCells / lazy matrices (e.g. "RenameDims") — subsets requested rows
#     BEFORE materialising so memory stays manageable, then coerces result to
#     a proper dgCMatrix (or dense matrix as last resort)
#
# Arguments:
#   seurat_object       Seurat object
#   assay     Assay name; NULL = DefaultAssay(seurat_object)
#   layer     Layer / slot name (default "data" = log-normalised counts)
#   features  Optional character vector of gene names to extract.
#             When supplied the BPCells matrix is subsetted BEFORE coercion,
#             which avoids materialising the full expression matrix.
#             NULL = return all genes.
#
# Returns a dgCMatrix (sparse) or matrix.  Row order matches 'features' when
# supplied and present in the object.
# =============================================================================
.get_layer_data <- function(seurat_object, assay = NULL, layer = "data",
                             features = NULL) {
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(seurat_object)

  # Attempt extraction: v5 LayerData → v5 GetAssayData(layer=) → v3 GetAssayData(slot=)
  mat <- tryCatch(
    SeuratObject::LayerData(seurat_object, assay = assay, layer = layer),
    error = function(e1) tryCatch(
      Seurat::GetAssayData(seurat_object, assay = assay, layer = layer),
      error = function(e2) tryCatch(
        Seurat::GetAssayData(seurat_object, assay = assay, slot = layer),
        error = function(e3)
          stop("Cannot extract '", layer, "' from assay '", assay,
               "': ", conditionMessage(e3))
      )
    )
  )

  # Subset rows before materialising (efficient for BPCells: only reads requested genes)
  if (!is.null(features) && length(features) > 0) {
    avail <- intersect(features, rownames(mat))
    if (length(avail) > 0) mat <- mat[avail, , drop = FALSE]
  }

  # Materialise any lazy / BPCells class (e.g. RenameDims) to a proper sparse matrix
  if (!is.matrix(mat) &&
      !inherits(mat, c("dgCMatrix", "dgRMatrix", "dgTMatrix",
                       "dsCMatrix", "CsparseMatrix"))) {
    mat <- tryCatch(
      methods::as(mat, "dgCMatrix"),
      error = function(e) as.matrix(mat)
    )
  }

  mat
}
