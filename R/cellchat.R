# =============================================================================
# scSidekick CellChat wrappers
#
# RunCellChat     - run the full CellChat pipeline for each level of a
#                   grouping variable, saving all standard visualizations
#                   (circle, chord, bubble, pathway, communication patterns)
#                   and the CellChat object as an RDS file.
#
# CompareCellChat - compare a named list of CellChat objects across all shared
#                   pathways, generating per-pathway circle, chord, and heatmap
#                   comparison PDFs.
#
# Both functions eliminate all setwd() calls, use explicit path construction,
# and accept species = "mouse" | "human" to select the right CellChatDB.
# =============================================================================

#' Run the full CellChat analysis pipeline per group
#'
#' For each level of `group.col` in the Seurat object, subsets the data,
#' builds a CellChat object, runs the standard pipeline
#' (overexpressed genes → interactions → communication probability →
#' pathway probability → aggregate network → communication patterns),
#' and writes visualizations + the CellChat RDS file to `output_dir`.
#'
#' @param seurat_object A Seurat object. Ignored when `redo_plots = TRUE`.
#' @param cell.by Character. Metadata column for cell identities (y-axis of
#'   chord/circle diagrams, e.g., `"subAssignment"`). `identity_column` is
#'   accepted as a deprecated alias.
#' @param group.by Character. Metadata column whose levels each get their own
#'   CellChat run (e.g., `"Group"`). `group.col` is accepted as a deprecated alias.
#' @param groups Character vector or `NULL`. Specific levels of `group.by` to
#'   process. `NULL` uses all levels.
#' @param species Character. Determines which CellChatDB ligand-receptor
#'   database and PPI network are loaded. One of:
#'   \itemize{
#'     \item `"human"` (default) - uses `CellChatDB.human` and `PPI.human`
#'     \item `"mouse"` - uses `CellChatDB.mouse` and `PPI.mouse`
#'   }
#' @param colors Named character vector mapping `cell.by` levels to colors.
#'   `NULL` uses CellChat's default palette. Tip: use the same color vector as
#'   your Seurat DimPlots so all figures are consistent.
#' @param assay_type Character or `NULL`. Seurat assay to use for ligand-receptor
#'   expression. `NULL` (default) auto-detects: uses `"RNA"` for Seurat v4 objects
#'   and `"RNA"` / `"Spatial"` for v5 objects depending on available assays.
#' @param caffeinate Logical. If `TRUE`, prevents the Mac from sleeping during
#'   the run using the `caffeinate` system command. Default `FALSE`.
#' @param output_dir Character. Root directory for all output. Per-group
#'   sub-folders are created automatically.
#' @param robj_dir Character or `NULL`. Directory for RDS checkpoint files.
#'   Defaults to `file.path(output_dir, "CellChatObjects")`.
#' @param min.cells Integer. Minimum number of cells a cell type must have in
#'   a given condition for its interactions to be retained by
#'   [CellChat::filterCommunication()]. Increase to reduce noise from rare
#'   populations. Default `0` (keep all).
#' @param thresh.p Numeric. P-value threshold used by
#'   [CellChat::identifyOverExpressedGenes()] and
#'   [CellChat::computeCommunProbPathway()]. Lowering this (e.g. `0.01`) makes
#'   the analysis more conservative; raising it (e.g. `0.1`) recovers more
#'   interactions. Default `0.05`.
#'
#'   **Note:** this is the pipeline significance threshold. The `thresh`
#'   parameter in [CompareCellChat()] independently controls which L-R
#'   interactions are *shown* in comparison plots.
#' @param run.patterns Logical. Run communication pattern analysis for
#'   outgoing and incoming signals across k values in `pattern_k_range`?
#'   Adds substantially to run time. Default `TRUE`.
#' @param pattern_k_range Integer vector. Range of latent pattern numbers (k)
#'   tested by [CellChat::identifyCommunicationPatterns()]. Each k×direction
#'   combination produces one PDF. Default `3:10`. Reduce (e.g. `3:6`) if
#'   run time is a concern.
#' @param vertex.receiver Integer vector. Indices of receiver cell types for
#'   the hierarchy layout in per-pathway plots. Default `seq(1, 4)` (first
#'   four cell types). Adjust to match your biological receiver populations.
#' @param save.rds Logical. Save each CellChat object as an RDS file?
#'   Default `TRUE`.
#' @param resume Logical. Skip groups that have already completed and recover
#'   groups interrupted mid-run. Default `FALSE`. Requires `save.rds = TRUE`
#'   and `output_dir` to be set. Two checkpoints are used:
#'   \itemize{
#'     \item **Pipeline checkpoint** (`{grp}CellChat_ckpt.rds`): written
#'       immediately after `netAnalysis_computeCentrality` (before any
#'       plotting). If a run is interrupted during plotting or pattern
#'       analysis, this file lets the next run skip the expensive pipeline
#'       steps and jump straight to visualizations.
#'     \item **Final RDS** (`{grp}CellChat.rds`): written at the end of a
#'       successful group run and replaces the checkpoint. If this file
#'       exists, that group is considered fully complete and skipped.
#'   }
#' @param redo_plots Logical. Re-generate all per-group visualization PDFs from
#'   existing RDS files without re-running the communication probability
#'   pipeline. Default `FALSE`. When `TRUE`:
#'   \itemize{
#'     \item For each group, `{grp}CellChat.rds` is loaded from `robj_dir`.
#'     \item All plotting code (circle, pathway, bubble, scatter, patterns) is
#'       re-executed using the loaded object.
#'     \item The pipeline (`computeCommunProb`, `aggregateNet`, etc.) is
#'       completely skipped - the `seurat_object` argument is ignored.
#'     \item Use this after [RenameCellTypeInCC()] + `saveRDS()` to regenerate
#'       figures with corrected cell-type labels.
#'   }
#'
#' @return A named list of CellChat objects (one per group level processed).
#' @export
RunCellChat <- function(seurat_object       = NULL,
                         cell.by,
                         group.by,
                         groups          = NULL,
                         # deprecated aliases
                         identity_column = NULL,
                         group.col       = NULL,
                         species         = "human",
                         colors          = NULL,
                         output_dir,
                         robj_dir        = NULL,
                         min.cells       = 0L,
                         thresh.p        = 0.05,
                         run.patterns    = TRUE,
                         pattern_k_range = 3:10,
                         vertex.receiver = seq(1, 4),
                         save.rds        = TRUE,
                         resume          = FALSE,
                         redo_plots      = FALSE,
                         assay_type      = NULL,
                         caffeinate      = FALSE) {

  if (caffeinate) { .caff <- .nk_caffeinate(); on.exit(.nk_decaffeinate(.caff), add = TRUE) }

  # Deprecated aliases
  if (!is.null(identity_column)) cell.by  <- identity_column
  if (!is.null(group.col))       group.by <- group.col

  if (!requireNamespace("CellChat", quietly = TRUE))
    stop("Package 'CellChat' is required. Install the latest version from GitHub:\n",
         "  devtools::install_github('jinworks/CellChat')\n",
         "scSidekick requires CellChat >= 2.1.0 (jinworks fork). ",
         "The older sqjin/CellChat is not supported.")

  if (!species %in% c("human", "mouse"))
    stop("species must be 'human' or 'mouse'.")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(robj_dir)) robj_dir <- file.path(output_dir, "CellChatObjects")
  dir.create(robj_dir, recursive = TRUE, showWarnings = FALSE)

  # Select species-specific resources
  CCdb  <- if (species == "mouse") CellChat::CellChatDB.mouse else
    CellChat::CellChatDB.human
  PPI   <- if (species == "mouse") CellChat::PPI.mouse else
    CellChat::PPI.human

  # redo_plots = TRUE: seurat_object is not used - groups come from robj_dir
  if (redo_plots) {
    if (is.null(groups)) {
      rds_found  <- list.files(robj_dir, pattern = "CellChat\\.rds$")
      grp_levels <- sub("CellChat\\.rds$", "", rds_found)
      if (length(grp_levels) == 0)
        stop("redo_plots = TRUE but no *CellChat.rds files found in robj_dir: ",
             robj_dir)
    } else {
      grp_levels <- groups
    }
  } else {
    if (is.null(seurat_object))
      stop("seurat_object is required when redo_plots = FALSE.")
    cell_info  <- seurat_object@meta.data
    grp_levels <- if (!is.null(groups)) groups else
      levels(factor(cell_info[[group.by]]))
  }

  # ── Resolve assay_type: explicit arg > parent params > default ──────────────
  effective_at <- assay_type
  if (is.null(effective_at)) {
    pj <- .find_params_json(dirname(normalizePath(output_dir, mustWork = FALSE)))
    if (!is.null(pj)) {
      pp <- tryCatch(jsonlite::read_json(pj, simplifyVector = TRUE), error = function(e) list())
      effective_at <- pp$assay_type
    }
  }
  effective_at <- effective_at %||% "scRNAseq"
  at_long <- switch(effective_at,
    snRNAseq   = "single-nucleus RNA-seq",
    scRNAseq   = "single-cell RNA-seq",
    scATACseq  = "single-cell ATAC-seq",
    scMultiome = "single-cell multiome",
    Spatial    = "spatial transcriptomics",
    Visium     = "Visium spatial transcriptomics",
    VisiumHD   = "Visium HD spatial transcriptomics",
    Xenium     = "Xenium in situ transcriptomics",
    effective_at   # passthrough
  )
  cell_word <- switch(effective_at,
    snRNAseq = "nuclei", Spatial = "spots", Visium = "spots",
    VisiumHD = "bins", "cells")

  # ── Cell count for this specific CellChat object ─────────────────────────────
  # The parent analysis_params.json n_cells_final reflects the FULL Seurat
  # object (e.g. 500k cells). RunCellChat may receive a smaller subset (e.g.
  # 50k sketched cells), so we compute n_cells_final from the actual seurat_object
  # and pass it explicitly to override the parent value.
  # In redo_plots mode seurat_object is NULL; fall back to the value already
  # stored in the local params JSON from the original run.
  n_cells_cc <- if (!is.null(seurat_object)) {
    ncol(seurat_object)
  } else {
    local_json <- file.path(output_dir, "analysis_params.json")
    if (file.exists(local_json)) {
      lp <- tryCatch(jsonlite::read_json(local_json, simplifyVector = TRUE),
                     error = function(e) list())
      lp$n_cells_final   # preserves value from original run
    } else NULL
  }
  nc_fmt <- if (!is.null(n_cells_cc))
    format(as.integer(n_cells_cc), big.mark = ",") else NULL

  # ── Write method params so create_analysis_pptx() finds them in output_dir ──
  # n_cells_final is set from seurat_object - overrides whatever the parent params
  # say so the CellChat PPTX reports the correct analyzed cell count.
  .write_subdir_params(output_dir, list(
    date                = format(Sys.Date()),
    assay_type          = effective_at,
    n_cells_final       = n_cells_cc,          # ← this object, not the parent
    cellchat_n_groups   = length(grp_levels),
    cellchat_species    = species,
    cellchat_ident_col  = cell.by,
    cellchat_group_col  = group.by,
    cellchat_groups     = as.list(grp_levels),
    cellchat_thresh_p   = thresh.p,
    cellchat_min_cells  = min.cells,
    cellchat_db         = paste0("CellChatDB.", species),
    methods_text        = paste0(
      "Cell-cell communication was inferred using CellChat (Jin et al., ",
      "Nature Communications, 2021) applied to ", species, " ",
      at_long, " data",
      if (!is.null(nc_fmt))
        paste0(" comprising ", nc_fmt, " ", cell_word) else "",
      " across ", length(grp_levels), " condition(s) (",
      paste(grp_levels, collapse = "; "), "). ",
      "Cell-type identities were defined by the '", cell.by,
      "' metadata column. ",
      "Ligand-receptor interactions were drawn from CellChatDB.", species,
      ". Overexpressed ligands and receptors were identified per cell type ",
      "using a Wilcoxon rank-sum test (p < ", thresh.p, "). ",
      "Communication probability was computed with a mass-action signaling ",
      "model; probabilities were aggregated to the pathway level and filtered ",
      "to remove interactions involving fewer than ", min.cells, " ", cell_word, ". ",
      "Signaling roles were characterized by network centrality analysis ",
      "(netAnalysis_computeCentrality). ",
      if (run.patterns)
        paste0("Communication patterns were identified for k = ",
               min(pattern_k_range), " to ", max(pattern_k_range),
               " latent patterns using non-negative matrix factorization. ")
      else "",
      "All analyses were conducted in R using the scSidekick package wrapper."
    ),
    cellchat_run_patterns    = isTRUE(run.patterns),
    cellchat_pattern_k_range = if (isTRUE(run.patterns)) as.list(pattern_k_range) else NA
  ))

  results <- list()

  for (grp in grp_levels) {
    message("\n========== CellChat: ", grp, " ==========")
    out_grp  <- file.path(output_dir, grp)
    dir.create(out_grp, showWarnings = FALSE)

    rds_final <- file.path(robj_dir, paste0(grp, "CellChat.rds"))
    rds_ckpt  <- file.path(robj_dir, paste0(grp, "CellChat_ckpt.rds"))

    # ── redo_plots: load final RDS and jump straight to visualizations ──────
    # The pipeline (computeCommunProb etc.) is completely skipped.
    # Use this after RenameCellTypeInCC() + saveRDS() to regenerate all
    # per-group PDFs with corrected cell-type labels.
    if (redo_plots) {
      if (!file.exists(rds_final)) {
        message("  redo_plots: no RDS found for ", grp, " - skipping.")
        next
      }
      cc_check <- try(readRDS(rds_final), silent = TRUE)
      if (!inherits(cc_check, "CellChat")) {
        message("  redo_plots: RDS for ", grp,
                " is not a valid CellChat object - skipping.")
        next
      }
      message("  redo_plots: loaded RDS for ", grp,
              " - re-running visualizations only.")
      cc <- cc_check

    } else {
      # ── Normal path: resume checks then pipeline ──────────────────────────

      # Skip if final RDS already exists and is valid (resume mode)
      if (resume && save.rds && file.exists(rds_final)) {
        cc_check <- try(readRDS(rds_final), silent = TRUE)
        if (inherits(cc_check, "CellChat")) {
          message("  Resuming: final RDS found - skipping ", grp)
          results[[grp]] <- cc_check
          next
        }
        message("  Final RDS found but not a valid CellChat - re-running ", grp)
      }

      # Pipeline checkpoint: skip expensive steps if present
      pipeline_done <- FALSE
      if (resume && file.exists(rds_ckpt)) {
        cc_check <- try(readRDS(rds_ckpt), silent = TRUE)
        if (inherits(cc_check, "CellChat")) {
          message("  Resuming: pipeline checkpoint found for ", grp,
                  " - skipping to visualizations")
          cc            <- cc_check
          pipeline_done <- TRUE
        }
      }

      if (!pipeline_done) {
        # Subset to group
        Seurat::Idents(seurat_object) <- cell_info[[group.by]]
        cells_grp <- Seurat::WhichCells(seurat_object, idents = grp)
        sub_obj   <- subset(seurat_object, cells = cells_grp)
        Seurat::Idents(sub_obj) <- sub_obj@meta.data[[cell.by]]

        # Build CellChat object
        data_input <- .get_layer_data(sub_obj, assay = "RNA", layer = "data")
        labels     <- Seurat::Idents(sub_obj)
        meta_df    <- data.frame(group = labels, row.names = names(labels))
        cc         <- CellChat::createCellChat(object   = data_input,
                                               meta     = meta_df,
                                               group.by = "group")
        cc@DB      <- CCdb

        # Standard pipeline
        options(future.globals.maxSize = 891289600)
        cc <- CellChat::subsetData(cc)
        cc <- CellChat::identifyOverExpressedGenes(cc, only.pos = TRUE,
                                                   do.fast = TRUE,
                                                   min.cells = min.cells,
                                                   thresh.p  = thresh.p)
        cc <- CellChat::identifyOverExpressedInteractions(cc)
        cc <- CellChat::smoothData(cc, adj = PPI)
        cc <- CellChat::computeCommunProb(cc)
        cc <- CellChat::filterCommunication(cc, min.cells = min.cells)
        cc <- CellChat::computeCommunProbPathway(cc, thresh = thresh.p)
        cc <- CellChat::aggregateNet(cc)
        cc <- CellChat::netAnalysis_computeCentrality(cc, slot.name = "netP")

        # Pipeline checkpoint: save before any plotting
        if (save.rds)
          saveRDS(cc, rds_ckpt)
      }
    } # end normal path

    grp_size <- as.numeric(table(cc@idents))
    grp_cols <- if (!is.null(colors)) colors[levels(cc@meta$group)] else NULL

    # ---- Circle plots ----
    f_count <- file.path(out_grp, paste0(grp, " Number of Interactions.pdf"))
    grDevices::pdf(f_count, width = 6, height = 6)
    CellChat::netVisual_circle(cc@net$count, vertex.weight = grp_size,
                               color.use = grp_cols, weight.scale = TRUE,
                               label.edge = FALSE,
                               title.name = "Number of interactions")
    grDevices::dev.off()
    .write_legend_sidecar(f_count, paste0(
      "Circle plot showing the total number of predicted ligand-receptor ",
      "interactions between cell populations in ", grp, ". Each node represents ",
      "a cell type; edge width is proportional to the number of interactions ",
      "and node size reflects the number of cells in each population."
    ))

    f_weight <- file.path(out_grp, paste0(grp, " Interaction Weights.pdf"))
    grDevices::pdf(f_weight, width = 6, height = 6)
    CellChat::netVisual_circle(cc@net$weight, vertex.weight = grp_size,
                               color.use = grp_cols, weight.scale = TRUE,
                               label.edge = FALSE,
                               title.name = "Interaction weights/strength")
    grDevices::dev.off()
    .write_legend_sidecar(f_weight, paste0(
      "Circle plot showing the aggregate interaction strength (communication ",
      "probability) between cell populations in ", grp, ". Edge width reflects ",
      "the total signaling weight of predicted interactions; node size reflects ",
      "the number of cells in each population."
    ))

    # ---- Individual groups ----
    mat    <- cc@net$weight
    f_indv <- file.path(out_grp, paste0(grp, " Individual Groups.pdf"))
    grDevices::pdf(f_indv, width = 10, height = 6)
    graphics::par(mfrow = c(2, 4), xpd = TRUE)
    for (b in seq_len(nrow(mat))) {
      mat2        <- matrix(0, nrow(mat), ncol(mat), dimnames = dimnames(mat))
      mat2[b, ]   <- mat[b, ]
      CellChat::netVisual_circle(mat2, vertex.weight = grp_size,
                                 weight.scale = TRUE, color.use = grp_cols,
                                 edge.weight.max = max(mat),
                                 title.name = rownames(mat)[b])
    }
    grDevices::dev.off()
    .write_legend_sidecar(f_indv, paste0(
      "Panel of circle plots showing outgoing interaction weights from each ",
      "individual cell population in ", grp, ". Each sub-panel displays ",
      "interactions initiated by a single sender cell type, revealing dominant ",
      "communication hubs. Edge weight is scaled to the global maximum across ",
      "all panels to allow direct comparison."
    ))

    # ---- Per-pathway plots ----
    paths     <- cc@netP$pathways
    path_dir  <- file.path(out_grp, "Pathways")
    dir.create(path_dir, showWarnings = FALSE)
    graphics::par(mfrow = c(1, 1))

    for (pw in paths) {
      message("  Pathway: ", pw)
      # All plots are called *inside* the PDF device.
      # netVisual_aggregate("circle") and ("chord") use circlize/base graphics
      # and draw directly to the current device - they must NOT be called
      # before pdf() is opened or their output goes to the null device.
      # netVisual_aggregate("hierarchy") and contribution/scatter return ggplot
      # objects and need print(); heatmap returns a ComplexHeatmap and needs draw().
      f_pw <- file.path(path_dir,
                        paste0(grp, " ", pw, " cell-cell communication.pdf"))
      grDevices::pdf(f_pw, width = 10, height = 7, onefile = TRUE)
      # 1. Network role (draws directly with grid/base graphics)
      try(CellChat::netAnalysis_signalingRole_network(cc, signaling = pw,
                                                      width = 8, height = 2.5,
                                                      font.size = 10,
                                                      color.use = grp_cols),
          silent = TRUE)
      # 2. Hierarchy layout (ggplot - needs print)
      try(print(CellChat::netVisual_aggregate(cc, signaling = pw,
                                              vertex.receiver = vertex.receiver,
                                              layout = "hierarchy",
                                              color.use = grp_cols)), silent = TRUE)
      # 3. Circle layout (base/circos - draws directly, no print)
      try(CellChat::netVisual_aggregate(cc, signaling = pw,
                                        layout = "circle",
                                        color.use = grp_cols), silent = TRUE)
      # 4. Chord layout (base/circos - draws directly, no print)
      try(CellChat::netVisual_aggregate(cc, signaling = pw,
                                        layout = "chord",
                                        color.use = grp_cols), silent = TRUE)
      # 5. Heatmap (ComplexHeatmap - print calls draw())
      try(print(CellChat::netVisual_heatmap(cc, signaling = pw,
                                            color.heatmap = "Reds",
                                            color.use = grp_cols)), silent = TRUE)
      # 6. L-R pair contribution (ggplot)
      try(print(CellChat::netAnalysis_contribution(cc, signaling = pw)),
          silent = TRUE)
      # 7. Sender/receiver scatter (ggplot)
      try(print(CellChat::netAnalysis_signalingRole_scatter(cc, signaling = pw,
                                                            title = pw,
                                                            color.use = grp_cols)),
          silent = TRUE)
      grDevices::dev.off()
      .write_legend_sidecar(f_pw, paste0(
        "Cell-cell communication networks for the ", pw, " signaling pathway ",
        "in ", grp, ". Panels show hierarchy, circle, and chord diagram layouts ",
        "of predicted ligand-receptor interactions, a heatmap of pairwise ",
        "communication probability, L-R pair contribution scores, and a scatter ",
        "plot of sender/receiver signaling roles."
      ))
    }

    # ---- Bubble plots per identity ----
    ids <- levels(cc@idents)
    for (j in seq_along(ids)) {
      bp <- try(CellChat::netVisual_bubble(cc, sources.use = j,
                                           remove.isolate = FALSE), silent = TRUE)
      f_bub <- file.path(out_grp, paste0(grp, " ", ids[j], " bubble plot.pdf"))
      grDevices::pdf(f_bub, pointsize = 8, width = 6, height = 16, onefile = TRUE)
      try(print(bp), silent = TRUE)
      grDevices::dev.off()
      .write_legend_sidecar(f_bub, paste0(
        "Bubble plot of predicted ligand-receptor interactions with ",
        ids[j], " as the source cell population in ", grp, ". ",
        "Each bubble represents one L-R pair; bubble size reflects ",
        "communication probability and color indicates interaction significance. ",
        "All target cell populations are shown on the y-axis."
      ))
    }

    # ---- Dominant sender/receiver ----
    gg1    <- CellChat::netAnalysis_signalingRole_scatter(cc)
    f_role <- file.path(out_grp, paste0(grp, " dominant sender receiver.pdf"))
    grDevices::pdf(f_role, width = 5, height = 5)
    print(gg1)
    grDevices::dev.off()
    .write_legend_sidecar(f_role, paste0(
      "Scatter plot of signaling roles across cell populations in ", grp, ". ",
      "Each point represents a cell type positioned by its outgoing (x-axis) ",
      "and incoming (y-axis) interaction strength, aggregated across all ",
      "signaling pathways. Cell types in the upper-right quadrant act as both ",
      "dominant senders and receivers of intercellular signals."
    ))

    # ---- Communication patterns ----
    if (run.patterns) {
      if (!requireNamespace("ggalluvial", quietly = TRUE))
        message("  'ggalluvial' not installed - skipping pattern river plots.")

      for (xx in pattern_k_range) {
        for (direction in c("outgoing", "incoming")) {
          message("  Patterns k=", xx, " (", direction, ")")
          # Use a temp variable so cc is NEVER overwritten with a try-error.
          # If the assignment were `cc <- try(...)` and it failed, cc would
          # become a try-error, subsequent iterations would also fail, and
          # the corrupted cc would be saveRDS()'d - breaking CompareCellChat.
          cc_pat <- try(CellChat::identifyCommunicationPatterns(
            cc, pattern = direction, k = xx, heatmap.show = FALSE
          ), silent = TRUE)
          if (inherits(cc_pat, "try-error")) next
          cc <- cc_pat   # update only on success

          colorsx <- if (!is.null(colors)) {
            clgrp <- if (direction == "outgoing")
              unique(cc@netP$pattern$outgoing$pattern$cell$CellGroup) else
              unique(cc@netP$pattern$incoming$pattern$cell$CellGroup)
            colors[levels(factor(clgrp))]
          } else NULL

          colors1 <- Nour_pal("main", reverse = TRUE)(xx)
          names(colors1) <- paste0("Pattern ", seq_len(xx))

          sig_lev <- try(unique(levels(
            if (direction == "outgoing")
              cc@netP$pattern$outgoing$pattern$signaling$Signaling else
              cc@netP$pattern$incoming$pattern$signaling$Signaling
          )), silent = TRUE)
          colors2 <- if (!inherits(sig_lev, "try-error") && length(sig_lev) > 0) {
            cv <- Nour_pal("cool", reverse = TRUE)(length(sig_lev))
            stats::setNames(cv, sig_lev)
          } else NULL

          f_pat <- file.path(out_grp,
                             paste0(grp, " k", xx, " ", direction,
                                    " communication patterns.pdf"))
          grDevices::pdf(f_pat, width = 15, height = 10, onefile = TRUE)
          try(CellChat::identifyCommunicationPatterns(
            cc, pattern = direction, k = xx, color.use = colorsx,
            color.heatmap = "RdYlBu", width = 12, height = 15
          ), silent = TRUE)
          try(print(CellChat::netAnalysis_river(cc, pattern = direction,
                                               color.use = colorsx,
                                               color.use.pattern = colors1,
                                               color.use.signaling = colors2)),
              silent = TRUE)
          try(print(CellChat::netAnalysis_dot(cc, pattern = direction,
                                             color.use = colorsx)),
              silent = TRUE)
          grDevices::dev.off()
          .write_legend_sidecar(f_pat, paste0(
            "Communication pattern analysis for ", direction, " signals ",
            "in ", grp, " (k = ", xx, " patterns). The heatmap shows the ",
            "contribution of each cell type to each latent communication ",
            "pattern; the river (alluvial) plot links patterns to the ",
            "signaling pathways they coordinate; the dot plot shows ",
            "pathway-level enrichment per pattern."
          ))
        }
      }
    }

    if (save.rds) {
      saveRDS(cc, rds_final)
      # Remove pipeline checkpoint now that the final RDS is written.
      # This keeps robj_dir clean and prevents stale checkpoints from
      # being mistakenly loaded on the next fresh run.
      if (file.exists(rds_ckpt)) unlink(rds_ckpt)
    }

    results[[grp]] <- cc
    message("  Done: ", grp)
  }

  results
}

# --------------------------------------------------------------------------- #
# CompareCellChat                                                               #
# --------------------------------------------------------------------------- #

#' Compare CellChat objects across groups for all shared pathways
#'
#' Takes a named list of CellChat objects (or paths to RDS files) and for
#' each pathway present in any object, generates a multi-panel comparison PDF
#' with circle, chord, and heatmap layouts.
#'
#' @param cellchat_list Named list. Each element is either a CellChat object
#'   or a character path to an RDS file. Names are used as panel titles.
#' @param output_dir Character. Directory for comparison PDFs.
#' @param groups.show Character vector or `NULL`. Cell identities to include.
#'   `NULL` uses all shared identities.
#' @param group.merged Named character vector or `NULL`. Maps identity names to
#'   broader categories for chord group aggregation. Names must match
#'   `groups.show`.
#' @param colors Named character vector mapping identity names to colors.
#' @param species Character. `"human"` or `"mouse"`. Default `"human"`.
#' @param thresh Numeric. P-value threshold applied to individual ligand-receptor
#'   interactions when drawing circle, chord, and heatmap plots. Default `0.05`
#'   (matches CellChat's standard significance filter). Set `thresh = 1` to
#'   display **all** predicted interactions regardless of statistical
#'   significance - this is essential for comparing pathways that are active
#'   in some conditions but did not reach pathway-level significance in others
#'   (i.e. the pathway is in `cc@net$prob` but not in `cc@netP$pathways`).
#'   Internally, all visualizations use the bundled `netVisual_aggregate2`-style
#'   functions which read `cc@net$prob` directly so the comparison is not
#'   limited to pathways that passed CellChat's aggregation threshold.
#' @param caffeinate Logical. If `TRUE`, prevents the Mac from sleeping during
#'   the run using the `caffeinate` system command. Default `FALSE`.
#'
#' @return A named list of the (possibly updated) CellChat objects.
#' @export
CompareCellChat <- function(cellchat_list,
                             output_dir,
                             groups.show   = NULL,
                             group.merged  = NULL,
                             colors        = NULL,
                             species       = "human",
                             thresh        = 0.05,
                             caffeinate    = FALSE) {

  if (caffeinate) { .caff <- .nk_caffeinate(); on.exit(.nk_decaffeinate(.caff), add = TRUE) }

  if (!requireNamespace("CellChat", quietly = TRUE))
    stop("Package 'CellChat' is required.")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Load RDS files if paths were provided
  cellchat_list <- lapply(cellchat_list, function(x) {
    if (is.character(x)) readRDS(x) else x
  })

  # Update objects to current CellChat version.
  # Keep the original on failure - never let a try-error replace a valid object.
  cellchat_list <- lapply(cellchat_list, function(cc) {
    if (!inherits(cc, "CellChat")) return(cc)
    cc_up <- try(CellChat::updateCellChat(cc), silent = TRUE)
    if (inherits(cc_up, "try-error")) cc else cc_up
  })

  # Compute centrality if not already done.
  # Same defensive pattern: on failure keep the current (valid) cc.
  cellchat_list <- lapply(cellchat_list, function(cc) {
    if (!inherits(cc, "CellChat")) return(cc)
    cc_up <- try(CellChat::netAnalysis_computeCentrality(cc, slot.name = "netP"),
                 silent = TRUE)
    if (inherits(cc_up, "try-error")) cc else cc_up
  })

  # Determine pathways to compare.
  # When thresh < 1 (default 0.05): union of cc@netP$pathways - same set
  #   CellChat considers significant.
  # When thresh >= 1 (show-all mode): also include pathways present in
  #   cc@LR$LRsig$pathway_name so interactions that exist at the L-R level
  #   but didn't reach pathway-level significance are still compared.
  all_pathways <- unique(unlist(lapply(cellchat_list, function(cc) {
    if (!inherits(cc, "CellChat")) return(character(0))
    pws <- cc@netP$pathways
    if (thresh >= 1 && !is.null(cc@LR$LRsig) &&
        "pathway_name" %in% colnames(cc@LR$LRsig))
      pws <- union(pws, cc@LR$LRsig$pathway_name[
        !is.na(cc@LR$LRsig$pathway_name)])
    pws
  })))

  # Subset to requested groups if specified
  if (!is.null(groups.show)) {
    cellchat_list <- lapply(cellchat_list, function(cc) {
      if (!inherits(cc, "CellChat")) return(cc)
      # Reorder idents and prob array
      cc@meta$group <- factor(cc@meta$group, levels = groups.show)
      cc@idents     <- cc@meta$group
      present       <- groups.show[groups.show %in% dimnames(cc@net$prob)[[1]]]
      cc@net$prob   <- cc@net$prob[present, present, , drop = FALSE]
      cc
    })
  }

  group_colors <- function(cc) {
    if (is.null(colors)) return(NULL)
    ids <- rownames(cc@net$prob)
    colors[ids[ids %in% names(colors)]]
  }

  # Determine the canonical cell-type set for heatmap padding.
  # netVisual_heatmap creates rows/cols from cc@netP$prob dimensions.
  # Objects where some cell types send/receive no signal in a given pathway
  # will have smaller netP$prob arrays, making Reduce("+", ht_list) crash
  # ("nrow of all heatmaps should be the same").
  # We compute the union here and pad each object's netP$prob per pathway
  # so every resulting heatmap has the same dimensions.
  all_cells_hm <- {
    cell_sets <- lapply(cellchat_list, function(cc) {
      if (!inherits(cc, "CellChat")) return(character(0))
      dimnames(cc@netP$prob)[[1]]
    })
    all_u <- unique(unlist(cell_sets))
    # Respect groups.show ordering when provided
    if (!is.null(groups.show))
      c(groups.show[groups.show %in% all_u], all_u[!all_u %in% groups.show])
    else
      all_u
  }

  for (pw in all_pathways) {
    message("Comparing pathway: ", pw)

    # Max weight across objects for comparable edge scaling
    weights <- vapply(cellchat_list, function(cc) {
      if (!inherits(cc, "CellChat")) return(NA_real_)
      w <- try(CellChat::getMaxWeight(list(cc), slot.name = "netP",
                                      attribute = pw), silent = TRUE)
      if (inherits(w, "try-error")) NA_real_ else as.numeric(w)
    }, numeric(1))
    w_max <- max(weights, na.rm = TRUE)
    if (!is.finite(w_max)) next

    n_obj   <- length(cellchat_list)
    obj_ids <- if (!is.null(groups.show)) groups.show else NULL

    f_cmp <- file.path(output_dir, paste0(pw, " signaling comparison.pdf"))
    grDevices::pdf(f_cmp, width = n_obj * 3 + 2, height = 12, onefile = TRUE)

    # All three rows use the bundled .cc_aggregate / .cc_chord_cell /
    # .cc_heatmap functions (adapted from Yun lab NetViasualHack.R).
    # Unlike standard CellChat functions, these read cc@net$prob (L-R pair
    # level) rather than cc@netP$prob (significance-filtered pathway level),
    # so pathways that have communication but did not reach pathway-level
    # significance are still shown - controlled by the `thresh` parameter.

    # Row 1: circle
    graphics::par(mfrow = c(1, n_obj), xpd = TRUE, mar = c(0, 0, 0, 0))
    for (i in seq_along(cellchat_list)) {
      cc <- cellchat_list[[i]]
      if (!inherits(cc, "CellChat")) next
      try(.cc_aggregate(
        cc, signaling = pw, layout = "circle", thresh = thresh,
        edge.weight.max = w_max, edge.width.max = 10,
        color.use       = group_colors(cc),
        signaling.name  = paste(names(cellchat_list)[i], pw)
      ), silent = TRUE)
    }

    # Row 2: chord (cell-type level)
    graphics::par(mfrow = c(1, n_obj), xpd = TRUE, mar = c(0, 0, 0, 0))
    for (i in seq_along(cellchat_list)) {
      cc <- cellchat_list[[i]]
      if (!inherits(cc, "CellChat")) next
      try(.cc_aggregate(
        cc, signaling = pw, layout = "chord", thresh = thresh,
        edge.weight.max = w_max, edge.width.max = 10,
        color.use       = group_colors(cc),
        signaling.name  = paste(names(cellchat_list)[i], pw)
      ), silent = TRUE)
    }

    # Row 2b (optional): merged chord - group.merged is a named character
    # vector mapping cell-type labels to broader category names, e.g.:
    #   c(Ex_L23_IT = "Excitatory", Inh_Sst = "Inhibitory", Astrocytes = "Glia")
    # This uses .cc_chord_cell which forwards `group` to circlize directly.
    if (!is.null(group.merged)) {
      graphics::par(mfrow = c(1, n_obj), xpd = TRUE, mar = c(0, 0, 0, 0))
      for (i in seq_along(cellchat_list)) {
        cc <- cellchat_list[[i]]
        if (!inherits(cc, "CellChat")) next
        try(.cc_chord_cell(
          cc, signaling = pw, thresh = thresh,
          group      = group.merged,
          color.use  = group_colors(cc),
          title.name = paste0(pw, " - ", names(cellchat_list)[i])
        ), silent = TRUE)
      }
    }

    # Row 3: heatmaps - .cc_heatmap reads netP$prob with a zero-matrix
    # fallback so conditions where the pathway is absent are shown as a
    # white/empty heatmap rather than being dropped from the combined panel.
    ht_list <- lapply(seq_along(cellchat_list), function(i) {
      cc  <- cellchat_list[[i]]
      nm  <- names(cellchat_list)[i]
      col <- group_colors(cc)
      if (!inherits(cc, "CellChat")) return(NULL)

      # Pad the cell-type set so all heatmaps have identical nrow
      current_cells <- tryCatch(dimnames(cc@netP$prob)[[1]], error = function(e) character(0))
      if (length(current_cells) > 0 && !identical(current_cells, all_cells_hm)) {
        n_pw   <- dim(cc@netP$prob)[3]
        pw_nm  <- dimnames(cc@netP$prob)[[3]]
        n      <- length(all_cells_hm)
        padded <- array(0, dim = c(n, n, n_pw),
                        dimnames = list(all_cells_hm, all_cells_hm, pw_nm))
        present <- intersect(current_cells, all_cells_hm)
        if (length(present) > 0)
          padded[present, present, ] <-
            cc@netP$prob[present, present, , drop = FALSE]
        cc@netP$prob <- padded   # local copy only
      }

      try(.cc_heatmap(
        cc, signaling = pw, color.heatmap = "Reds",
        title.name = paste(pw, nm),
        color.use  = col
      ), silent = TRUE)
    })
    ht_list <- Filter(function(x) !is.null(x) && !inherits(x, "try-error"),
                      ht_list)
    if (length(ht_list) > 0 && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      ht_combined <- Reduce(`+`, ht_list)
      try(ComplexHeatmap::draw(ht_combined,
                               ht_gap = grid::unit(0.5, "cm")), silent = TRUE)
    }

    grDevices::dev.off()
    try(grDevices::dev.off(), silent = TRUE)
    .write_legend_sidecar(f_cmp, paste0(
      "Cross-condition comparison of the ", pw, " signaling pathway across ",
      n_obj, " CellChat object(s) (p-value threshold: ", thresh, "). ",
      "Row 1: circle plots of communication probability. ",
      "Row 2: chord diagrams of cell-type interactions. ",
      if (!is.null(group.merged))
        "Row 2b: merged chord diagrams with cell types grouped by broader category. "
      else "",
      "Row 3: heatmaps of pairwise communication probability. ",
      "Visualizations read ligand-receptor pair probabilities directly (cc@net$prob) ",
      "so conditions where the pathway did not reach pathway-level significance ",
      "are still shown rather than being silently excluded."
    ))
  }

  message("CompareCellChat complete - ", length(all_pathways), " pathways processed.")
  invisible(cellchat_list)
}

# --------------------------------------------------------------------------- #
# RankCellChatPathways                                                          #
# --------------------------------------------------------------------------- #

#' Rank and visualize differentially active CellChat pathways across conditions
#'
#' Extracts per-pathway information flow and sender/receiver cell-type marginals
#' from a named list of CellChat objects, ranks pathways by differential
#' activity, and produces four outputs:
#' \enumerate{
#'   \item A **flow heatmap**: top `top_n` pathways × conditions, row-normalized
#'         so color reflects which condition has higher activity.
#'   \item **Sender and receiver tile plots**: faceted by top pathway, showing
#'         how much each cell type sends/receives per condition.
#'   \item A **top changed links heatmap**: the `top_links` sender→receiver pairs
#'         with the highest variance across conditions, showing exactly which
#'         cell-type interactions drive the pathway-level changes.
#'   \item A **CSV** of all underlying metrics (flow, sender, receiver, link
#'         tables).
#' }
#'
#' @param cellchat_list Named list of CellChat objects (one per condition).
#'   Names are used as condition labels throughout.
#' @param output_dir Character or `NULL`. Directory for PDFs and CSVs.
#'   `NULL` returns objects only.
#' @param top_n Integer. Number of top pathways (ranked by flow range across
#'   conditions) to include in visualizations. Default `20`.
#' @param top_links Integer. Number of top sender→receiver links (ranked by
#'   variance of flow across conditions) to show in the links heatmap.
#'   Default `15`.
#' @param rank_by Character. Metric used to rank pathways:
#'   `"range"` (max − min flow; default) or `"cv"` (coefficient of variation).
#' @param colors Named character vector. Condition colors (names = condition
#'   names). `NULL` uses ggplot2 defaults.
#' @param cell_colors Named character vector. Cell-type colors (names = cell
#'   type labels). Used for annotation strips in tile plots.
#'
#' @return Invisibly, a named list:
#'   \describe{
#'     \item{`flow_table`}{data.frame: pathway × condition flow + ranking cols.}
#'     \item{`sender_table`}{long data.frame: cell_type, condition, pathway, flow.}
#'     \item{`receiver_table`}{long data.frame: same for receivers.}
#'     \item{`links_table`}{long data.frame: sender, receiver, pathway, condition, flow.}
#'     \item{`flow_heatmap`}{pheatmap object.}
#'     \item{`sender_plot`}{ggplot2 object.}
#'     \item{`receiver_plot`}{ggplot2 object.}
#'     \item{`links_plot`}{ggplot2 object.}
#'   }
#' @export
RankCellChatPathways <- function(cellchat_list,
                                  output_dir  = NULL,
                                  top_n       = 20L,
                                  top_links   = 15L,
                                  rank_by     = "range",
                                  colors      = NULL,
                                  cell_colors = NULL) {

  if (!requireNamespace("CellChat", quietly = TRUE))
    stop("Package 'CellChat' is required.")

  stopifnot(is.list(cellchat_list), !is.null(names(cellchat_list)))
  if (!rank_by %in% c("range", "cv"))
    stop("rank_by must be 'range' or 'cv'.")

  if (!is.null(output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Load RDS files if paths were provided - same as CompareCellChat.
  cellchat_list <- lapply(cellchat_list, function(x) {
    if (is.character(x)) readRDS(x) else x
  })

  grp_names    <- names(cellchat_list)
  valid        <- vapply(cellchat_list, inherits, logical(1), "CellChat")
  if (!any(valid)) stop("No valid CellChat objects found in cellchat_list.")

  # ── Step 1: Extract flow, sender/receiver marginals, and link tables ────────
  message("Extracting information flow from ", sum(valid), " CellChat objects...")

  all_pathways <- unique(unlist(lapply(cellchat_list[valid], function(cc)
    dimnames(cc@netP$prob)[[3]]))  )
  all_cells    <- unique(unlist(lapply(cellchat_list[valid], function(cc)
    dimnames(cc@netP$prob)[[1]])))

  # flow_mat: pathway × condition
  flow_mat <- matrix(0, nrow = length(all_pathways), ncol = length(grp_names),
                     dimnames = list(all_pathways, grp_names))

  sender_rows   <- list()   # accumulate long-format rows
  receiver_rows <- list()
  link_rows     <- list()

  for (gi in seq_along(cellchat_list)) {
    grp <- grp_names[gi]
    cc  <- cellchat_list[[gi]]
    if (!inherits(cc, "CellChat")) next

    prob    <- cc@netP$prob       # senders × receivers × pathways
    cc_pws  <- dimnames(prob)[[3]]
    cc_cells <- dimnames(prob)[[1]]

    for (pw in cc_pws) {
      mat <- prob[, , pw, drop = TRUE]   # cell × cell matrix
      if (all(mat == 0)) next

      flow_mat[pw, gi] <- sum(mat)

      # Sender marginals (rowSums = how much each cell type sends)
      s_flow <- rowSums(mat)
      s_flow <- s_flow[s_flow > 0]
      if (length(s_flow) > 0)
        sender_rows[[length(sender_rows) + 1]] <- data.frame(
          cell_type = names(s_flow), condition = grp,
          pathway = pw, flow = s_flow, row.names = NULL,
          stringsAsFactors = FALSE)

      # Receiver marginals (colSums)
      r_flow <- colSums(mat)
      r_flow <- r_flow[r_flow > 0]
      if (length(r_flow) > 0)
        receiver_rows[[length(receiver_rows) + 1]] <- data.frame(
          cell_type = names(r_flow), condition = grp,
          pathway = pw, flow = r_flow, row.names = NULL,
          stringsAsFactors = FALSE)

      # All non-zero sender→receiver links
      mat_df          <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
      colnames(mat_df) <- c("sender", "receiver", "flow")
      mat_df$pathway   <- pw
      mat_df$condition <- grp
      mat_df           <- mat_df[mat_df$flow > 0, ]
      if (nrow(mat_df) > 0)
        link_rows[[length(link_rows) + 1]] <- mat_df
    }
  }

  sender_table   <- do.call(rbind, sender_rows)
  receiver_table <- do.call(rbind, receiver_rows)
  links_table    <- do.call(rbind, link_rows)
  rownames(sender_table) <- rownames(receiver_table) <-
    rownames(links_table) <- NULL

  # ── Step 2: Rank pathways ───────────────────────────────────────────────────
  rank_score <- if (rank_by == "cv") {
    apply(flow_mat, 1, function(r) {
      m <- mean(r); if (m == 0) 0 else stats::sd(r) / m
    })
  } else {
    apply(flow_mat, 1, function(r) max(r) - min(r))
  }

  # Drop pathways that are zero in all conditions
  rank_score     <- rank_score[apply(flow_mat, 1, max) > 0]
  rank_score     <- sort(rank_score, decreasing = TRUE)
  top_pathways   <- names(rank_score)[seq_len(min(top_n, length(rank_score)))]

  # Assemble flow table with ranking columns
  flow_table <- as.data.frame(flow_mat)
  flow_table$pathway     <- rownames(flow_mat)
  flow_table$range       <- apply(flow_mat, 1, function(r) max(r) - min(r))
  flow_table$cv          <- apply(flow_mat, 1, function(r) {
    m <- mean(r); if (m == 0) 0 else stats::sd(r) / m })
  flow_table$rank_score  <- rank_score[rownames(flow_table)]
  flow_table             <- flow_table[order(-flow_table$range), ]

  # ── Step 3: Flow heatmap ────────────────────────────────────────────────────
  message("Building flow heatmap...")
  sub_flow <- flow_mat[top_pathways, , drop = FALSE]

  # Row-normalize to [0, 1]: color shows "high in which condition"
  sub_norm <- t(apply(sub_flow, 1, function(r) {
    rng <- range(r)
    if (rng[2] == rng[1]) rep(0.5, length(r)) else (r - rng[1]) / diff(rng)
  }))

  # Build pheatmap annotation only when every condition name is present in
  # colors - pheatmap::pheatmap errors if annotation_colors has fewer names
  # than factor levels (e.g. when colors is GroupColors keyed by full group
  # names but the CellChat list is keyed by shorter condition labels).
  ann_col <- if (!is.null(colors) &&
                 all(colnames(sub_norm) %in% names(colors))) {
    ac       <- data.frame(Condition = colnames(sub_norm),
                           row.names = colnames(sub_norm))
    cond_col <- colors[colnames(sub_norm)]        # exact subset, no NAs
    list(df = ac, colors = list(Condition = cond_col))
  } else NULL

  flow_hm_path <- if (!is.null(output_dir))
    file.path(output_dir, "CellChat_pathway_flow_heatmap.pdf") else NULL

  flow_heatmap <- pheatmap::pheatmap(
    sub_norm,
    cluster_rows  = TRUE, cluster_cols = FALSE,
    color         = viridis::viridis(50),
    border_color  = NA,
    fontsize_row  = 8, fontsize_col  = 9,
    cellwidth     = 22, cellheight    = 10,
    angle_col     = 45,
    annotation_col = if (!is.null(ann_col)) ann_col$df else NULL,
    annotation_colors = if (!is.null(ann_col)) ann_col$colors else NULL,
    main = paste0("Top ", length(top_pathways), " pathways by flow ",
                  rank_by, "  (row-normalized to [0,1])"),
    filename = if (!is.null(flow_hm_path)) flow_hm_path else NA
  )

  if (!is.null(flow_hm_path))
    .write_legend_sidecar(flow_hm_path, paste0(
      "Heatmap of the top ", length(top_pathways), " CellChat signaling pathways ",
      "ranked by information-flow ", rank_by, " across ", length(grp_names),
      " conditions. Each value is the total communication probability summed ",
      "across all sender-receiver pairs for that pathway, row-normalized to ",
      "[0-1] so color encodes which condition has relatively higher activity. ",
      "Pathways are hierarchically clustered; conditions retain their input order."
    ))

  # ── Step 4: Sender and Receiver tile plots ──────────────────────────────────
  message("Building sender/receiver tile plots...")

  .make_tile_plot <- function(long_df, role = "Sender") {
    df <- long_df[long_df$pathway %in% top_pathways, , drop = FALSE]
    if (nrow(df) == 0) return(ggplot2::ggplot())

    # Normalize flow within each pathway so color = relative activity
    df <- do.call(rbind, lapply(split(df, df$pathway), function(x) {
      mx <- max(x$flow, na.rm = TRUE)
      x$flow_norm <- if (mx > 0) x$flow / mx else x$flow
      x
    }))

    # Factor condition to preserve input order
    df$condition <- factor(df$condition, levels = grp_names)

    # For legibility, order cell types by total flow (descending)
    ct_order <- names(sort(tapply(df$flow, df$cell_type, sum),
                           decreasing = TRUE))
    df$cell_type <- factor(df$cell_type, levels = rev(ct_order))

    # Pathway order = same as top_pathways ranking
    df$pathway <- factor(df$pathway, levels = rev(top_pathways))

    fill_label  <- paste0(role, " flow\n(pathway-normalized)")
    fill_colors <- if (!is.null(colors))
      ggplot2::scale_fill_manual(values = colors) else
      ggplot2::scale_fill_viridis_c(option = "plasma", name = fill_label)

    ggplot2::ggplot(df, ggplot2::aes(x = condition, y = cell_type,
                                     fill = flow_norm)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.3) +
      ggplot2::scale_fill_viridis_c(option = "plasma", name = fill_label,
                                    limits = c(0, 1)) +
      ggplot2::facet_wrap(~ pathway, ncol = 4) +
      ggplot2::labs(title = paste(role, "information flow - top pathways"),
                    x = NULL, y = NULL) +
      theme_NourMin() +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y  = ggplot2::element_text(size  = 7),
        strip.text   = ggplot2::element_text(size  = 7, face = "bold"),
        panel.spacing = ggplot2::unit(0.3, "lines")
      )
  }

  sender_plot   <- .make_tile_plot(sender_table,   "Sender")
  receiver_plot <- .make_tile_plot(receiver_table, "Receiver")

  n_pw  <- length(top_pathways)
  n_rows_tile <- ceiling(n_pw / 4)
  tile_h <- max(6, n_rows_tile * 3 + 1.5)   # scale height to number of rows

  for (info in list(
    list(plot = sender_plot,   suffix = "sender_tile",
         legend = "Sender information flow (row = cell type, column = condition, facet = pathway). Color encodes each cell type's outgoing communication probability for that pathway, normalized within each pathway panel to [0-1]."),
    list(plot = receiver_plot, suffix = "receiver_tile",
         legend = "Receiver information flow (row = cell type, column = condition, facet = pathway). Color encodes each cell type's incoming communication probability for that pathway, normalized within each pathway panel to [0-1].")
  )) {
    if (!is.null(output_dir)) {
      p_path <- file.path(output_dir,
                          paste0("CellChat_pathway_", info$suffix, ".pdf"))
      ggplot2::ggsave(info$plot, filename = p_path,
                      width = 16, height = tile_h, limitsize = FALSE)
      .write_legend_sidecar(p_path, info$legend)
    }
  }

  # ── Step 5: Top changed links heatmap ──────────────────────────────────────
  message("Building top changed links heatmap...")

  links_plot <- ggplot2::ggplot()   # fallback empty plot

  if (!is.null(links_table) && nrow(links_table) > 0) {

    links_table$link_label <- paste0(links_table$sender, " → ",
                                     links_table$receiver, "\n(",
                                     links_table$pathway, ")")

    # Variance of flow across conditions per link - rank by this
    lv <- tapply(links_table$flow, links_table$link_label, stats::var)
    lv[is.na(lv)] <- 0
    top_link_labels <- names(sort(lv, decreasing = TRUE))[
      seq_len(min(top_links, length(lv)))]

    sub_links <- links_table[links_table$link_label %in% top_link_labels, ]
    sub_links$condition <- factor(sub_links$condition, levels = grp_names)
    sub_links$link_label <- factor(sub_links$link_label,
                                   levels = top_link_labels)

    # Color the sender cell type if cell_colors provided
    sender_fill <- if (!is.null(cell_colors)) {
      senders <- unique(sub_links$sender)
      sc <- cell_colors[senders[senders %in% names(cell_colors)]]
      list(
        ggplot2::geom_tile(ggplot2::aes(fill = flow), color = "white",
                           linewidth = 0.3),
        ggplot2::scale_fill_viridis_c(option = "inferno", name = "Flow")
      )
    } else {
      list(
        ggplot2::geom_tile(ggplot2::aes(fill = flow), color = "white",
                           linewidth = 0.3),
        ggplot2::scale_fill_viridis_c(option = "inferno", name = "Flow")
      )
    }

    links_plot <- ggplot2::ggplot(
      sub_links, ggplot2::aes(x = condition, y = link_label)) +
      sender_fill +
      ggplot2::labs(
        title = paste0("Top ", min(top_links, length(lv)),
                       " sender→receiver links by flow variance"),
        x = NULL, y = NULL
      ) +
      theme_NourMin() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = ggplot2::element_text(size = 7)
      )

    if (!is.null(output_dir)) {
      lp_path <- file.path(output_dir,
                           "CellChat_top_changed_links_heatmap.pdf")
      link_h  <- max(4, min(top_links, length(lv)) * 0.45 + 2)
      ggplot2::ggsave(links_plot, filename = lp_path,
                      width = max(5, length(grp_names) * 1.2 + 3),
                      height = link_h, limitsize = FALSE)
      .write_legend_sidecar(lp_path, paste0(
        "Heatmap of the top ", min(top_links, length(lv)),
        " sender→receiver-pathway links ranked by variance of ",
        "communication probability across all conditions. Each row is a unique ",
        "cell-type pair and signaling pathway; each column is a condition. ",
        "Color = raw communication probability (viridis inferno scale). ",
        "High variance rows identify the specific intercellular links that ",
        "change most between conditions."
      ))
    }
  }

  # ── Step 6: Write CSVs ──────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    utils::write.csv(flow_table,     file.path(output_dir, "CellChat_flow_table.csv"))
    utils::write.csv(sender_table,   file.path(output_dir, "CellChat_sender_table.csv"),   row.names = FALSE)
    utils::write.csv(receiver_table, file.path(output_dir, "CellChat_receiver_table.csv"), row.names = FALSE)
    utils::write.csv(links_table,    file.path(output_dir, "CellChat_links_table.csv"),    row.names = FALSE)
    message("CSVs written to ", output_dir)
  }

  message("RankCellChatPathways complete.")
  invisible(list(
    flow_table     = flow_table,
    sender_table   = sender_table,
    receiver_table = receiver_table,
    links_table    = links_table,
    flow_heatmap   = flow_heatmap,
    sender_plot    = sender_plot,
    receiver_plot  = receiver_plot,
    links_plot     = links_plot
  ))
}

# --------------------------------------------------------------------------- #
# RenameCellTypeInCC                                                            #
# --------------------------------------------------------------------------- #

#' Rename a cell type label inside a CellChat object without re-running analysis
#'
#' Patches every slot that carries cell-type names (idents, meta, net arrays,
#' netP arrays, centrality scores, pattern tables, LR significant pairs) so
#' that downstream visualization functions and `CompareCellChat` /
#' `RankCellChatPathways` see the new name.
#'
#' The communication probability matrices are **not** recomputed - only labels
#' change. This is safe because the underlying values are invariant to what the
#' cell type is called.
#'
#' **Typical workflow after renaming:**
#' ```r
#' # 1. Patch and re-save all RDS files
#' rds_files <- list.files("CellChatObjects/", "CellChat\\.rds$", full.names = TRUE)
#' for (f in rds_files) {
#'   cc <- RenameCellTypeInCC(readRDS(f), "Immune_contam", "NK_TCells")
#'   saveRDS(cc, f)
#' }
#'
#' # 2. Re-run comparison figures (fast - pure plotting)
#' cclist <- setNames(lapply(rds_files, readRDS),
#'                    sub("CellChat\\.rds$", "", basename(rds_files)))
#' CompareCellChat(cclist, output_dir = ...)
#' RankCellChatPathways(cclist, output_dir = ...)
#'
#' # 3. Per-group RunCellChat PDFs: delete old PDFs and re-run RunCellChat with
#' #    resume = TRUE - the pipeline checkpoint means only plots are regenerated.
#' ```
#'
#' @param cc A CellChat object.
#' @param old_name Character. The cell type label to replace.
#' @param new_name Character. The replacement label.
#'
#' @return The modified CellChat object (the original is not changed in place -
#'   assign the result back).
#' @export
RenameCellTypeInCC <- function(cc, old_name, new_name) {

  if (!inherits(cc, "CellChat"))
    stop("cc must be a CellChat object.")
  if (!is.character(old_name) || length(old_name) != 1)
    stop("old_name must be a single character string.")
  if (!is.character(new_name) || length(new_name) != 1)
    stop("new_name must be a single character string.")

  # ── Helpers ─────────────────────────────────────────────────────────────────

  # Rename one level/value in a factor or character vector
  ren_vec <- function(x) {
    if (is.factor(x)) {
      levels(x)[levels(x) == old_name] <- new_name
    } else if (is.character(x)) {
      x[x == old_name] <- new_name
    }
    x
  }

  # Rename in the first two dimensions of any array/matrix dimnames
  ren_dim12 <- function(x) {
    dn <- dimnames(x)
    if (length(dn) >= 1 && !is.null(dn[[1]]))
      dn[[1]][dn[[1]] == old_name] <- new_name
    if (length(dn) >= 2 && !is.null(dn[[2]]))
      dn[[2]][dn[[2]] == old_name] <- new_name
    dimnames(x) <- dn
    x
  }

  # Rename in the names() of a named vector or list
  ren_names <- function(x) {
    if (!is.null(names(x)))
      names(x)[names(x) == old_name] <- new_name
    x
  }

  # Safely attempt a slot update; skip silently on error
  try_patch <- function(expr) tryCatch(expr, error = function(e) NULL)

  # ── Cell identity slots ──────────────────────────────────────────────────────
  try_patch(cc@idents      <- ren_vec(cc@idents))
  try_patch(cc@meta$group  <- ren_vec(cc@meta$group))

  # ── Network probability arrays (sender × receiver × LR or pathway) ──────────
  for (nm in c("prob", "pval", "count", "count.sp", "weight")) {
    try_patch({
      if (!is.null(cc@net[[nm]]))
        cc@net[[nm]] <- ren_dim12(cc@net[[nm]])
    })
  }
  try_patch({
    if (!is.null(cc@netP$prob))
      cc@netP$prob <- ren_dim12(cc@netP$prob)
  })

  # ── Centrality scores (named list → named vectors) ───────────────────────────
  try_patch({
    if (!is.null(cc@netP$centr))
      cc@netP$centr <- lapply(cc@netP$centr, ren_names)
  })

  # ── Communication pattern tables (from identifyCommunicationPatterns) ────────
  for (direction in c("outgoing", "incoming")) {
    try_patch({
      pat <- cc@netP$pattern[[direction]]$pattern$cell
      if (!is.null(pat) && "CellGroup" %in% names(pat)) {
        pat$CellGroup <- ren_vec(pat$CellGroup)
        cc@netP$pattern[[direction]]$pattern$cell <- pat
      }
    })
  }

  # ── LR significant pairs (source / target columns) ───────────────────────────
  for (col in c("source", "target")) {
    try_patch({
      if (!is.null(cc@LR$LRsig) && col %in% names(cc@LR$LRsig))
        cc@LR$LRsig[[col]] <- ren_vec(cc@LR$LRsig[[col]])
    })
  }

  # ── Aggregate communication data frame (if present) ──────────────────────────
  try_patch({
    if (!is.null(cc@communication)) {
      for (col in c("source", "target")) {
        if (col %in% names(cc@communication))
          cc@communication[[col]] <- ren_vec(cc@communication[[col]])
      }
    }
  })

  message("Renamed '", old_name, "' → '", new_name, "' in CellChat object.")
  cc
}
