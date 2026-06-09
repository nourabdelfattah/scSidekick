# =============================================================================
# scSidekick utility functions
#
# Includes the package-wide null-coalescing operator %||% used throughout.
#
# ExtractMethods - reads seurat_object@commands to auto-generate a methods paragraph
#                  summarising preprocessing, integration, clustering, and
#                  dimensionality reduction steps.
#
# InspectPlot    - inspects a ggplot2 object and returns a plain-text
#                  description of what was plotted (geom types, aesthetic
#                  mappings, faceting, scales). For scSidekick functions the
#                  description is attached as attr(plot, "legend_text").
# =============================================================================

# Package-wide null-coalescing operator (not exported - internal only)
# Usage: a %||% b  →  returns a if non-NULL, else b
`%||%` <- function(a, b) if (!is.null(a)) a else b


# ── Sleep prevention helpers ──────────────────────────────────────────────────
# .nk_caffeinate() / .nk_decaffeinate() keep the machine awake during long
# computations (CellChat, GSEA, annotation, etc.).
#
# macOS   — spawns `caffeinate -i -w <R PID>` (built-in since 10.8). The -w
#            flag auto-exits caffeinate when this R process exits, so cleanup
#            is guaranteed even if on.exit() is skipped due to a crash.
# Windows — uses processx (widely available; dep of devtools/pkgdown/usethis)
#            to run a hidden PowerShell loop calling SetThreadExecutionState
#            every 30 s. Falls back silently if processx is absent.
# Linux   — no-op (servers don't sleep; desktop users can set a global
#            inhibitor via systemd-inhibit externally if needed).
#
# All failure paths return NULL. .nk_decaffeinate(NULL) is always a no-op,
# so these helpers can never break the calling function.

.nk_caffeinate <- function() {
  os <- Sys.info()[["sysname"]]

  if (os == "Darwin") {
    # system2(..., wait=FALSE) forks without creating a pipe — no COW memory lock
    # that would starve subsequent large allocations (e.g. ULM matrix operations).
    # -w <pid>: caffeinate auto-exits when this R process exits, so no PID needed.
    ok <- tryCatch({
      system2("caffeinate", args = c("-i", "-w", as.character(Sys.getpid())),
              wait = FALSE, stdout = FALSE, stderr = FALSE)
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok) {
      message("scSidekick: caffeinate active — machine will stay awake.")
      return(list(os = "Darwin", pid = NA_integer_))
    }

  } else if (os == "Windows") {
    if (!requireNamespace("processx", quietly = TRUE)) return(NULL)
    ps1 <- tempfile(fileext = ".ps1")
    writeLines(c(
      'Add-Type -MemberDefinition \'[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);\' -Name SleepHelper -Namespace Win32 -ErrorAction SilentlyContinue',
      'while ($true) {',
      '  [Win32.SleepHelper]::SetThreadExecutionState(0x80000003)',
      '  Start-Sleep -Seconds 30',
      '}'
    ), ps1)
    proc <- try(
      processx::process$new(
        "powershell",
        c("-ExecutionPolicy", "Bypass", "-NonInteractive", "-File", ps1),
        windows_verbatim_args = FALSE,
        stdout = NULL, stderr = NULL
      ),
      silent = TRUE
    )
    if (!inherits(proc, "try-error") && proc$is_alive()) {
      message("scSidekick: Windows sleep prevention active.")
      return(list(os = "Windows", proc = proc, ps1 = ps1))
    }
    try(unlink(ps1), silent = TRUE)
  }

  NULL  # Linux or any failure — silent no-op
}

.nk_decaffeinate <- function(handle) {
  if (is.null(handle)) return(invisible(NULL))

  if (handle$os == "Darwin") {
    try(tools::pskill(handle$pid), silent = TRUE)

  } else if (handle$os == "Windows") {
    try(handle$proc$kill(), silent = TRUE)
    # Reset ES so Windows can sleep again
    ps_reset <- paste0(
      'Add-Type -MemberDefinition \'[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);\' ',
      '-Name SleepHelper -Namespace Win32 -ErrorAction SilentlyContinue; ',
      '[Win32.SleepHelper]::SetThreadExecutionState(0x80000000)'
    )
    try(system(
      paste("powershell -NonInteractive -ExecutionPolicy Bypass -Command",
            shQuote(ps_reset)),
      ignore.stdout = TRUE, ignore.stderr = TRUE
    ), silent = TRUE)
    try(unlink(handle$ps1), silent = TRUE)
  }

  invisible(NULL)
}


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
#' @param modality Character or `NULL`. Override the sequencing modality label
#'   used in the methods text. Typical values: `"scRNA-seq"`, `"snRNA-seq"`.
#'   `NULL` (default) auto-detects from `seurat_object@misc$scSidekick_params`
#'   or the `"modality"` metadata column.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{`summary`}{Named list of extracted parameter values.}
#'   \item{`methods_text`}{Character string - a draft methods paragraph
#'     ready to paste into a manuscript.}
#' }
#' @export
ExtractMethods <- function(seurat_object, cite_seurat = TRUE, modality = NULL) {

  cmds  <- names(seurat_object@commands)
  get_p <- function(cmd, param) {
    if (cmd %in% cmds) seurat_object@commands[[cmd]]@params[[param]] else NULL
  }

  # ---- Modality (sc vs sn) ----
  # Priority: explicit argument > @misc$scSidekick_params$modality
  #           > @misc$scSidekick_params$assay_type (set by log_analysis_params)
  #           > modality metadata column > default "scRNA-seq"
  modality <- modality %||%
    seurat_object@misc$scSidekick_params$modality %||%
    seurat_object@misc$scSidekick_params$assay_type %||%
    if ("modality" %in% colnames(seurat_object@meta.data))
      unique(seurat_object@meta.data$modality)[1]
    else "scRNA-seq"
  is_sn       <- grepl("sn", modality, ignore.case = TRUE)
  unit_label  <- if (is_sn) "nuclei" else "cells"
  method_term <- if (is_sn) "Single-nucleus RNA-seq (snRNA-seq)"
                 else       "Single-cell RNA-seq (scRNA-seq)"

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
    modality           = modality,
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
    "%s data (%s %s, %s genes; assay: %s) were analyzed using Seurat.",
    method_term,
    format(n_cells, big.mark = ","),
    unit_label,
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
      "%s were clustered using the %s algorithm at resolution %.2f, yielding %d clusters.",
      if (is_sn) "Nuclei" else "Cells",
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


# =============================================================================
# RenameClusters - rename numeric cluster IDs to prefixed, zero-padded labels
# =============================================================================

#' Rename numeric cluster IDs to prefixed labels
#'
#' Converts 0-based or 1-based integer cluster IDs (as produced by Seurat's
#' \code{FindClusters}) into zero-padded labels with a custom prefix, e.g.
#' cluster 0 → \code{"C01"}, cluster 1 → \code{"C02"}.
#'
#' @param x A vector of numeric or character cluster IDs (e.g. \code{0, 1, 2, ...}).
#' @param prefix Character prefix to prepend. Default \code{"C"}.
#' @param start Integer. The output number assigned to the lowest cluster.
#'   \code{start = 1} (default) maps cluster 0 → \code{"C01"}.
#'   \code{start = 0} maps cluster 0 → \code{"C00"}.
#'
#' @return A factor with levels in cluster order.
#' @export
#'
#' @examples
#' \dontrun{
#' SeuratObj$Cluster <- RenameClusters(SeuratObj$seurat_clusters, prefix = "C", start = 1)
#' # cluster 0 → "C01", cluster 1 → "C02", ...
#' }
RenameClusters <- function(x, prefix = "C", start = 1L) {
  x_num <- suppressWarnings(as.integer(as.character(x)))
  if (anyNA(x_num))
    stop("'x' must contain numeric cluster IDs. Non-numeric values found: ",
         paste(unique(as.character(x)[is.na(x_num)]), collapse = ", "))

  sorted_ids  <- sort(unique(x_num))
  n_clusters  <- length(sorted_ids)
  max_out_num <- start + n_clusters - 1L
  pad_width   <- nchar(as.character(max_out_num))

  out_labels  <- paste0(prefix,
                         sprintf(paste0("%0", pad_width, "d"),
                                 seq(start, start + n_clusters - 1L)))
  names(out_labels) <- as.character(sorted_ids)

  factor(out_labels[as.character(x_num)], levels = out_labels)
}


# =============================================================================
# NumberLabels - prepend/append sequential numbers to factor levels
#               or create group × id interaction labels
# =============================================================================

#' Prepend or append sequential numbers to factor levels
#'
#' Two modes:
#' \describe{
#'   \item{Mode 1 (no \code{id_by})}{Numbers the levels sequentially:
#'     \code{"01.B cell"}, \code{"02.CD4 T"}, …  Useful for shortening long
#'     names on UMAP centroid labels while keeping the full name in the legend.}
#'   \item{Mode 2 (with \code{id_by})}{Creates a \emph{group × id} interaction
#'     label.  Each unique value of \code{id_by} within a group gets its own
#'     sequential number: \code{"Male_1"}, \code{"Male_2"}, \code{"Female_1"}, …}
#' }
#'
#' If \code{x} is not already a factor it is converted automatically (with a
#' message).  Run \code{SetLevels()} first if you need a specific level order.
#'
#' @param x A vector (factor, character, or numeric).
#' @param id_by Optional vector of the same length as \code{x}.  When supplied,
#'   enables Mode 2: one sequential number per unique \code{id_by} value within
#'   each level of \code{x}.
#' @param position \code{"before"} (default) places the number before the label;
#'   \code{"after"} places it after.
#' @param sep Separator between the number and the label.  Default \code{"."}.
#' @param pad_width Integer.  Width for zero-padding.  \code{NULL} (default)
#'   auto-detects based on the maximum number.
#'
#' @return A factor with relabeled levels.
#' @export
#'
#' @examples
#' \dontrun{
#' # Mode 1: number long SingleR labels for UMAP
#' SeuratObj$SingleRAssignment <- NumberLabels(SeuratObj$SingleRAssignment)
#' # levels: "01.B cell", "02.CD4 T", ...
#'
#' # Mode 2: one number per patient within each Sex group
#' SeuratObj$SexByPatient <- NumberLabels(SeuratObj$Sex, id_by = SeuratObj$patient)
#' # levels: "Female_1", "Female_2", "Male_1", "Male_2", ...
#'
#' # Append number after, custom separator
#' NumberLabels(SeuratObj$CellType, position = "after", sep = " - ")
#' }
NumberLabels <- function(x,
                          id_by     = NULL,
                          position  = c("before", "after"),
                          sep       = ".",
                          pad_width = NULL) {
  position <- match.arg(position)

  if (!is.factor(x)) {
    x <- factor(x)
    message("scSidekick: Variable was converted to a factor. ",
            "Run SetLevels() first if you need a specific level order.")
  }

  lvls <- levels(x)

  if (is.null(id_by)) {
    # ── Mode 1: sequential numbering of levels ───────────────────────────────
    n     <- length(lvls)
    width <- pad_width %||% nchar(as.character(n))
    nums  <- sprintf(paste0("%0", width, "d"), seq_len(n))

    new_lvls  <- if (position == "before") paste0(nums, sep, lvls)
                 else                       paste0(lvls, sep, nums)
    label_map <- stats::setNames(new_lvls, lvls)
    return(factor(as.character(label_map[as.character(x)]), levels = new_lvls))
  }

  # ── Mode 2: group × id interaction ───────────────────────────────────────
  if (length(id_by) != length(x))
    stop("'id_by' must be the same length as 'x'.")

  grp_char   <- as.character(x)
  id_char    <- as.character(id_by)
  labels_out <- character(length(x))
  all_new_lvls <- character(0)

  for (grp in lvls) {
    in_grp     <- grp_char == grp
    unique_ids <- unique(id_char[in_grp])
    n_ids      <- length(unique_ids)
    width      <- pad_width %||% max(1L, nchar(as.character(n_ids)))
    nums       <- stats::setNames(
      sprintf(paste0("%0", width, "d"), seq_len(n_ids)),
      unique_ids
    )

    num_mapped          <- nums[id_char[in_grp]]
    labels_out[in_grp]  <- if (position == "before") paste0(num_mapped, sep, grp)
                            else                       paste0(grp, sep, num_mapped)

    grp_lvls     <- if (position == "before") paste0(nums, sep, grp)
                    else                       paste0(grp, sep, nums)
    all_new_lvls <- c(all_new_lvls, sort(grp_lvls))
  }

  factor(labels_out, levels = all_new_lvls)
}


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


# =============================================================================
# SetLevels - relevel a factor (or any vector) without repeating its name
# =============================================================================

#' Relevel a factor or vector concisely
#'
#' A drop-in replacement for the verbose
#' \code{factor(x, levels = c(...))} pattern.  Converts non-factors
#' automatically (with a message), places any listed levels first in the
#' specified order, and appends all remaining levels at the end - so you only
#' need to list the levels you care about.
#'
#' Designed for pipe- and assignment-friendly use:
#' \preformatted{
#' SeuratObj$Group <- SetLevels(SeuratObj$Group, c("Control", "Treated"))
#' }
#'
#' @param x A vector (factor, character, or numeric).
#' @param levels Character vector of desired levels, in order.  Levels present
#'   in \code{x} but not listed here are appended at the end unchanged.
#' @param append_rest Logical.  If \code{TRUE} (default), unlisted levels are
#'   appended after the specified ones.  If \code{FALSE}, unlisted levels are
#'   dropped (values become \code{NA}).
#'
#' @return A factor with the requested level ordering.
#' @export
#'
#' @examples
#' \dontrun{
#' # Reorder two key levels; remaining levels appended automatically
#' SeuratObj$Group <- SetLevels(SeuratObj$Group, c("Control", "Treated"))
#'
#' # Full explicit ordering
#' SeuratObj$Dementia <- SetLevels(SeuratObj$Dementia,
#'   c("Reference.Control", "No dementia.Control", "Dementia.Control",
#'     "No dementia.AD",    "Dementia.AD"))
#' }
SetLevels <- function(x, levels, append_rest = TRUE) {

  if (!is.factor(x)) {
    x <- factor(x)
    message("scSidekick: Variable was converted to a factor. ",
            "If you need a different order, use `levels = c(\"level1\", \"level2\", ...)`.")
  }

  existing <- base::levels(x)

  # Warn about any requested level not actually in the data
  phantom <- setdiff(levels, existing)
  if (length(phantom))
    warning("scSidekick: The following levels are not present in x and will be ignored: ",
            paste(phantom, collapse = ", "))

  levels <- intersect(levels, existing)   # keep only real ones, in requested order

  final_levels <- if (append_rest) {
    c(levels, setdiff(existing, levels))  # listed first, rest appended
  } else {
    levels                                # strict - unlisted become NA
  }

  factor(x, levels = final_levels)
}


#' Inspect a ggplot2 object and describe what it shows
#'
#' Extracts metadata from a ggplot2 object - geom types, key aesthetic
#' mappings, axis labels, faceting, and color/fill scale names - and returns
#' a plain-text description suitable for use as a figure legend draft.
#'
#' For plots created by scSidekick functions, a pre-formatted legend is stored
#' as `attr(plot, "legend_text")` and is returned directly when present.
#'
#' @param p A ggplot2 object.
#' @param title Character or `NULL`. Optional figure panel title to include
#'   in the legend text (e.g., `"Figure 2A"`).
#'
#' @return A character string - a draft legend sentence.
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
         "attr(p, 'legend_text') - make sure the function returned the plot ",
         "object (out_dir = NULL) and that the legend_text attribute was set.")

  # ---- Geom types ----
  # Seurat uses GeomDrawGrob / GeomScattermore internally - these are opaque
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
  col_lab   <- lbs$color %||% lbs$color %||% lbs$fill %||% ""
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
  # A FeaturePlot colors by a continuous gene/score; a DimPlot by a discrete
  # identity. Seurat puts the feature name in p$labels$color for FeaturePlot
  # and "Identity" / the group.by column for DimPlot.
  color_scale_type <- tryCatch({
    sc <- Filter(function(s) any(s$aesthetics %in% c("color","color","fill")),
                 p$scales$scales)
    if (length(sc) == 0) "unknown"
    else if (inherits(sc[[1]], "ScaleContinuous")) "continuous"
    else "discrete"
  }, error = function(e) "unknown")

  seurat_plot_type <-
    if (nchar(reduction) > 0 && color_scale_type == "continuous")
      "feature plot"
    else if (nchar(reduction) > 0 && color_scale_type == "discrete")
      "cell identity plot"
    else if (nchar(reduction) > 0)
      "embedding plot"
    else NULL

  # ---- Color / fill scale limits ----
  scale_detail <- tryCatch({
    sc <- Filter(function(s) any(s$aesthetics %in% c("color","color","fill")),
                 p$scales$scales)
    if (length(sc) == 0 || color_scale_type != "continuous") return(NULL)
    lims <- sc[[1]]$limits
    if (!is.null(lims) && !any(is.na(lims)))
      sprintf("color range %.2g - %.2g", lims[1], lims[2])
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
    # e.g. "UMAP feature plot of APOE, colored by expression (range 0-3.2)"
    feature_part <- if (nchar(col_lab) > 0 && col_lab != "value")
      paste0("of ", col_lab) else ""
    panel_part   <- if (nchar(plot_title) > 0 && plot_title != col_lab)
      paste0(" [panel: ", plot_title, "]") else ""
    color_part  <- if (color_scale_type == "continuous")
      paste0("colored by expression",
             if (!is.null(scale_detail)) paste0(" (", scale_detail, ")") else "")
    else if (nchar(col_lab) > 0)
      paste0("colored by ", col_lab)
    else ""

    bits <- Filter(nzchar, c(
      if (!is.null(title))       title,
      paste(reduction, seurat_plot_type, feature_part),
      panel_part,
      color_part,
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
      if (nchar(plot_title) > 0)     paste("-", plot_title),
      if (nchar(aes_str) > 0)        paste("showing", aes_str),
      if (nchar(x_lab) > 0 && x_lab != "x")   paste("x:", x_lab),
      if (nchar(y_lab) > 0 && y_lab != "y")   paste("y:", y_lab),
      if (nchar(col_lab) > 0)        paste("colored by", col_lab),
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
    message("jsonlite not available - skipping analysis_params.json write.")
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
#   • Seurat v5  - uses SeuratObject::LayerData() with layer= parameter
#   • Seurat v3  - falls back to Seurat::GetAssayData() with slot= parameter
#   • BPCells / lazy matrices (e.g. "RenameDims") - subsets requested rows
#     BEFORE materialising so memory stays manageable, then coerces result to
#     a proper dgCMatrix (or dense matrix as last resort)
#
# Arguments:
#   seurat_object       Seurat object
#   assay     Assay name; NULL = DefaultAssay(seurat_object)
#   layer     Layer / slot name (default "data" = log-normalized counts)
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


# =============================================================================
# .nk_feature_matrix()
#
# Resolve a set of "features" to a features x cells numeric matrix, pulling each
# feature from the assay (genes) OR from a numeric meta.data column (e.g. cNMF
# usages, module scores, QC metrics). Lets plotting functions accept gene names
# and continuous metadata interchangeably.
#
# Returns list(mat = features x cells matrix of the found features in the
# requested order, found = character, missing = character). mat is NULL when
# nothing is found.
# =============================================================================
.nk_feature_matrix <- function(seurat_object, features, assay = NULL,
                                layer = "data") {
  assay    <- assay %||% SeuratObject::DefaultAssay(seurat_object)
  features <- unique(as.character(features))
  meta     <- seurat_object@meta.data

  avail_genes <- tryCatch(rownames(seurat_object[[assay]]),
                          error = function(e)
                            tryCatch(rownames(seurat_object),
                                     error = function(e2) character(0)))

  is_meta_num <- function(f) f %in% colnames(meta) && is.numeric(meta[[f]])

  gene_feats <- features[features %in% avail_genes]
  meta_feats <- features[!(features %in% avail_genes) &
                           vapply(features, is_meta_num, logical(1))]
  found      <- features[features %in% c(gene_feats, meta_feats)]
  missing    <- setdiff(features, found)

  if (length(found) == 0L)
    return(list(mat = NULL, found = character(0), missing = missing))

  parts <- list()
  if (length(gene_feats) > 0L) {
    gm <- .get_layer_data(seurat_object, assay = assay, layer = layer,
                          features = gene_feats)
    parts$genes <- as.matrix(gm[gene_feats, , drop = FALSE])
  }
  if (length(meta_feats) > 0L) {
    mm <- t(as.matrix(meta[, meta_feats, drop = FALSE]))
    rownames(mm) <- meta_feats
    colnames(mm) <- rownames(meta)
    parts$meta <- mm
  }

  common <- Reduce(intersect, lapply(parts, colnames))
  mat <- do.call(rbind, lapply(parts, function(m) m[, common, drop = FALSE]))
  mat <- mat[found, , drop = FALSE]
  list(mat = mat, found = found, missing = missing)
}

# Single-feature convenience: named numeric vector (named by cell).
.nk_feature_values <- function(seurat_object, feature, assay = NULL,
                               layer = "data") {
  fm <- .nk_feature_matrix(seurat_object, feature, assay = assay, layer = layer)
  if (is.null(fm$mat))
    stop("Feature '", feature, "' not found in assay '",
         assay %||% SeuratObject::DefaultAssay(seurat_object),
         "' or as a numeric meta.data column.")
  stats::setNames(as.numeric(fm$mat[1, ]), colnames(fm$mat))
}


# =============================================================================
# .msigdbr_get()
#
# Version-agnostic wrapper around msigdbr::msigdbr().  msigdbr 10.0.0 renamed
# the arguments  category -> collection  and  subcategory -> subcollection,
# and split some MSigDB subcollections (notably KEGG: the old "CP:KEGG" became
# "CP:KEGG_LEGACY" / "CP:KEGG_MEDICUS").  This helper lets the rest of the
# package keep using the historical category / subcategory vocabulary and the
# old "CP:KEGG" code, and translates to whatever the installed msigdbr expects.
#
#   category    : top-level collection code, e.g. "H", "C2", "C5"
#   subcategory : subcollection code, e.g. "CP:KEGG", "CP:REACTOME" (or NULL)
#   species     : e.g. "Homo sapiens", "Mus musculus"
#
# Returns the msigdbr data frame.  The output always carries gs_name and
# gene_symbol columns, which is all downstream fgsea code relies on.
# =============================================================================
.msigdbr_get <- function(species, category, subcategory = NULL) {
  new_api <- "collection" %in% names(formals(msigdbr::msigdbr))

  # Try a list of candidate subcategory codes in order; first that returns rows
  # (without erroring) wins.  Covers the KEGG split in MSigDB 2023+.
  sub_candidates <- if (is.null(subcategory)) {
    list(NULL)
  } else if (identical(subcategory, "CP:KEGG")) {
    # Old objects: try legacy KEGG first, then MEDICUS, then the original code.
    list("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", "CP:KEGG")
  } else {
    list(subcategory)
  }

  # Build the argument list dynamically, omitting NULL collection/subcollection
  # (NULL category = "all collections"), and using the right arg names per API.
  call_msigdbr <- function(sub) {
    args <- list(species = species)
    if (new_api) {
      if (!is.null(category)) args$collection    <- category
      if (!is.null(sub))      args$subcollection <- sub
    } else {
      if (!is.null(category)) args$category    <- category
      if (!is.null(sub))      args$subcategory <- sub
    }
    do.call(msigdbr::msigdbr, args)
  }

  last_err <- NULL
  for (sub in sub_candidates) {
    res <- tryCatch(call_msigdbr(sub), error = function(e) { last_err <<- e; NULL })
    if (!is.null(res) && nrow(res) > 0) return(res)
  }

  stop("msigdbr returned no gene sets for collection '", category, "'",
       if (!is.null(subcategory)) paste0(" / subcollection '", subcategory, "'"),
       ". ", if (!is.null(last_err)) conditionMessage(last_err) else "",
       call. = FALSE)
}


# =============================================================================
# .as_inches()
# Coerce a size argument to inches.  Accepts a grid::unit (converted) OR a bare
# numeric (interpreted as inches).  Falls back to `default_in` when conversion
# fails - this fixes the bug where passing width = 10 (a number, not a unit)
# silently reverted the saved PDF to the hard-coded default size.
# =============================================================================
.as_inches <- function(x, default_in) {
  if (inherits(x, "unit")) {
    out <- tryCatch(as.numeric(grid::convertUnit(x, "in")),
                    error = function(e) NA_real_)
    if (is.na(out)) default_in else out
  } else if (is.numeric(x) && length(x) == 1L && is.finite(x)) {
    x
  } else {
    default_in
  }
}


# =============================================================================
# .heatmap_pdf_dims()
#
# Compute PDF width/height (inches) for a ComplexHeatmap so that long row and
# column names are never clipped.  Sizes the device as:
#
#   body  +  row-name label space (one side)  +  legend  +  margins
#   body  +  column-name label space (top/bottom)  +  title  +  annotations
#
# Character width is estimated as 0.6 * fontsize(pt) / 72 inches per character.
# Column-name space depends on rotation: vertical (90 deg) consumes height equal
# to the label length; horizontal (0 deg) consumes none vertically; 45 deg
# splits the difference via sin/cos.
#
#   body_w_in, body_h_in : heatmap BODY size in inches
#   row_names, col_names : the actual label vectors (for max nchar)
#   row_fontsize, col_fontsize : label font sizes in pt
#   row_names_side  : "left" | "right"   (which side carries row labels)
#   column_names_side : "top" | "bottom"
#   column_names_rot  : label rotation in degrees (0, 45, 90 typical)
#   legend_in    : horizontal allowance for the heatmap legend
#   extra_right_in : extra right allowance (e.g. a dot-size legend)
#   title_in     : vertical allowance for the column title
#   n_top_anno   : number of top annotation bars (each ~0.25 in)
#   margin_in    : base margin on every side
#   max_in       : hard cap per dimension (device sanity limit)
#
# Returns list(width, height) in inches.
# =============================================================================
.heatmap_pdf_dims <- function(body_w_in, body_h_in,
                              row_names, col_names,
                              row_fontsize     = 8,
                              col_fontsize     = 10,
                              row_names_side   = "right",
                              column_names_side = "bottom",
                              column_names_rot = 90,
                              legend_in        = 1.5,
                              extra_right_in   = 0,
                              extra_left_in    = 0,
                              title_in         = 0.4,
                              n_top_anno       = 0,
                              n_col_split      = 0,
                              margin_in        = 0.4,
                              max_in           = 200) {

  row_char_w <- 0.6 * row_fontsize / 72   # inches per character (row labels)
  col_char_w <- 0.6 * col_fontsize / 72   # inches per character (column labels)

  max_rn <- if (length(row_names)) max(nchar(as.character(row_names)), na.rm = TRUE) else 0
  max_cn <- if (length(col_names)) max(nchar(as.character(col_names)), na.rm = TRUE) else 0
  if (!is.finite(max_rn)) max_rn <- 0
  if (!is.finite(max_cn)) max_cn <- 0

  # Row-name label space (horizontal text on one side)
  rn_space <- max_rn * row_char_w + 0.1

  # Column-name space: how much VERTICAL room the (possibly rotated) labels need,
  # plus a small HORIZONTAL overhang for the outermost label.
  rot          <- column_names_rot %% 180
  vert_frac    <- if (rot >= 80 && rot <= 100) 1 else if (rot == 0) 0 else sin(rot * pi / 180)
  horiz_frac   <- if (rot == 0) 1 else cos(rot * pi / 180)
  cn_space_v   <- max_cn * col_char_w * vert_frac + 0.1
  cn_overhang  <- max_cn * col_char_w * horiz_frac * 0.5

  # Column-split gap: ComplexHeatmap inserts ~5mm between each split panel
  col_split_gap_in <- max(0L, n_col_split - 1L) * 0.2

  # Width
  left_in  <- margin_in + extra_left_in +
              (if (identical(row_names_side, "left"))  rn_space else 0)
  right_in <- margin_in + (if (identical(row_names_side, "right")) rn_space else 0) +
              legend_in + extra_right_in
  width_in <- body_w_in + left_in + right_in + cn_overhang + col_split_gap_in

  # Height: column-split titles add ~0.3 in each on top
  col_split_title_in <- if (n_col_split > 0) 0.3 else 0
  top_in    <- margin_in + title_in + col_split_title_in + n_top_anno * 0.25 +
               (if (identical(column_names_side, "top"))    cn_space_v else 0)
  bottom_in <- margin_in +
               (if (identical(column_names_side, "bottom")) cn_space_v else 0)
  height_in <- body_h_in + top_in + bottom_in

  list(width  = min(max_in, max(3, width_in)),
       height = min(max_in, max(3, height_in)))
}
