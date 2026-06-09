# =============================================================================
# RunGSEA_pseudobulk - donor/sample-level (pseudobulk) GSEA + per-patient ssGSEA
# -----------------------------------------------------------------------------
# Sibling of RunGSEA(): instead of cell-level Wilcoxon (cells as replicates),
# aggregates to pseudobulk per sample x cell type (samples = replicates), fits
# limma-voom with a group-means design, extracts user contrasts, ranks genes by
# the moderated t statistic and runs fgsea. Also computes per-patient pathway
# activity (ssGSEA via GSVA) on the pseudobulk and tests it across groups.
# Emits the same artifact set as RunGSEA (per-cell-type NES heatmaps, lollipops,
# CSVs, .legend sidecars).
# =============================================================================

#' Pseudobulk GSEA with limma contrasts and per-patient ssGSEA
#'
#' @description
#' The statistically rigorous alternative to [RunGSEA()]. Instead of treating
#' individual cells as replicates, this function aggregates counts to the
#' donor/sample level (pseudobulk), fits a limma-voom linear model, and ranks
#' genes by the moderated t-statistic before running fgsea. This respects the
#' true unit of replication and avoids the inflated degrees of freedom that
#' arise from cell-level comparisons.
#'
#' Optionally also computes per-patient pathway activity scores (ssGSEA via
#' GSVA) and tests them across groups, producing boxplots and heatmaps
#' alongside the GSEA results.
#'
#' @section Three required columns:
#' Three metadata columns drive the analysis:
#' \describe{
#'   \item{`group.by`}{Cell type column. GSEA is run independently for each
#'     level (e.g. `"Assignment"`). Cells are aggregated within each level.}
#'   \item{`sample.by`}{Donor/sample ID column. This is the statistical unit.
#'     Each unique value becomes one pseudobulk replicate in the limma model
#'     (e.g. `"DonorID"`, `"SampleID"`).}
#'   \item{`contrast.by`}{Comparison group column. Levels define the groups
#'     being contrasted (e.g. `"Disease"` with levels `"AD"` and `"Control"`).
#'     Cells with `NA` in this column are dropped before aggregation.}
#' }
#'
#' @section Gene set input modes:
#' Supply gene sets using exactly one of the following (priority order):
#' \enumerate{
#'   \item `gene_sets` - a named list of character vectors, one vector per gene set.
#'   \item `deg_df` - a differential expression data frame; top DEGs per group
#'     are extracted and used as gene sets (see `deg_*` parameters).
#'   \item `pathway_sets` (default) - one or more MSigDB collections fetched via
#'     `msigdbr`. Does not require `gene_sets` or `deg_df`.
#' }
#'
#' @param seurat_object A Seurat object. Either `seurat_object` or a
#'   precomputed `pseudobulk` list must be supplied.
#' @param group.by Character or `NULL`. Metadata column holding cell-type
#'   labels. One independent GSEA is run per level (e.g. `"Assignment"`).
#'   Pass `NULL` when the object already contains only one cell type (e.g. a
#'   Microglia-only object); all cells are treated as a single group called
#'   `"All"` and a single GSEA is run.
#' @param sample.by Character. Metadata column identifying the biological
#'   replicate (donor/patient ID). This is the statistical unit; each unique
#'   value becomes one pseudobulk sample in the limma-voom model
#'   (e.g. `"DonorID"`).
#' @param contrast.by Character. Metadata column holding the comparison group
#'   label (e.g. `"Disease"` with levels `"AD"` and `"Control"`). Cells with
#'   `NA` in this column are excluded before aggregation.
#' @param split.by Character or `NULL`. Optional metadata column used to
#'   create interaction groups before fitting the model. When supplied, the
#'   contrast groups become `<contrast.by>.<split.by>` combinations
#'   (e.g. `split.by = "Sex"` turns `"AD"` and `"Control"` into
#'   `"AD.Female"`, `"AD.Male"`, `"Control.Female"`, `"Control.Male"`),
#'   allowing sex-stratified or condition-stratified contrasts. Default `NULL`.
#' @param subset.by Character or `NULL`. Optional metadata column used to run
#'   the full pipeline separately for each level. Each level gets its own
#'   subdirectory under `output_dir`. Use this to produce independent analyses
#'   per broad cell class, brain region, or any other grouping variable without
#'   calling the function multiple times. Default `NULL`.
#' @param run_label Character or `NULL`. Short label appended to output
#'   filenames and plot titles to distinguish this run from others saved to the
#'   same `output_dir`. Set automatically when `subset.by` is used; otherwise
#'   leave `NULL` unless you need a custom tag.
#' @param pseudobulk Optional precomputed pseudobulk. A named list with
#'   elements `$counts` (genes x samples matrix) and `$samples` (sample-level
#'   metadata data frame). Supply the value returned by a previous
#'   `RunGSEA_pseudobulk()` call to skip the aggregation step on reruns.
#'   Default `NULL` (aggregation is performed from `seurat_object`).
#' @param covariates Character vector of sample-level metadata columns to
#'   include as covariates in the limma-voom model (e.g. `c("Age", "PMI")`).
#'   Numeric covariates are median-imputed for `NA`; factor/character
#'   covariates receive an `"unknown"` level. Default `NULL`.
#' @param contrasts Named list of limma contrast expressions. Every element
#'   must be a **single string** built from model coefficient names of the form
#'   `make.names(paste0("group", level))`. A small helper makes this concise:
#'   \preformatted{
#'   cf <- function(l) make.names(paste0("group", l))
#'   contrasts <- list(
#'     AD_vs_Ctrl   = paste(cf("AD"),     "-", cf("Control")),
#'     Interaction  = paste0("(", cf("AD.F"), "-", cf("Ctrl.F"), ")",
#'                           " - ",
#'                           "(", cf("AD.M"), "-", cf("Ctrl.M"), ")")
#'   )
#'   }
#'   Positive NES = enrichment in the first term. When `NULL` a standard set
#'   is generated automatically from the levels of `contrast.by`. When
#'   `redo_plots = TRUE`, a user-supplied list takes priority over the list
#'   saved in the RDS. Default `NULL`.
#'
#' @param gene_sets Named list of character vectors. Each element is one gene
#'   set (e.g. `list(MySet = c("GENE1", "GENE2", ...))`). When supplied,
#'   `pathway_sets` and `deg_df` are ignored and `msigdbr` is not required.
#'   Default `NULL`.
#' @param deg_df Data frame of differential expression results. Top DEGs per
#'   group are extracted and used as gene sets. When supplied, `pathway_sets`
#'   is ignored and `msigdbr` is not required. Use the `deg_*` parameters
#'   below to map your data frame's column names. Default `NULL`.
#' @param deg_gene_column Character. Column in `deg_df` containing gene names.
#'   Default `"feature"`.
#' @param deg_group_column Character (length 1 or 2). Column(s) in `deg_df`
#'   defining the group each row belongs to. Two columns create a two-level
#'   label used for heatmap row splitting. Default `"group"`.
#' @param deg_fc_column Character. Column in `deg_df` containing the fold
#'   change used to rank genes within each group. Default `"logFC"`.
#' @param deg_padj_column Character or `NULL`. Column in `deg_df` containing
#'   adjusted p-values for pre-filtering. Set to `NULL` to skip filtering.
#'   Default `"padj"`.
#' @param deg_padj_cutoff Numeric. Adjusted p-value threshold applied to
#'   `deg_df` before selecting top genes. Default `0.05`.
#' @param deg_top_n Integer. Number of top genes per group (ranked by
#'   `deg_fc_column`, descending) to include in each gene set. Default `20`.
#'
#' @param pathway_sets Named list of MSigDB collections to test. Each element
#'   must have a `category` field and an optional `subcategory` field
#'   (see [msigdbr::msigdbr_collections()]). The element name becomes the
#'   label used in output filenames. Default: Hallmark, KEGG, Reactome,
#'   WikiPathways. Ignored when `gene_sets` or `deg_df` is supplied.
#' @param species Character. Species name passed to [msigdbr::msigdbr()].
#'   Use `"human"` or `"mouse"` (shorthand accepted) or any full name returned
#'   by `msigdbr::msigdbr_species()`. Default `"human"`.
#' @param output_dir Character. Root directory for all CSV, PDF, and RDS
#'   output. Per-contrast and per-database subdirectories are created
#'   automatically.
#' @param assay Character or `NULL`. Seurat assay containing raw counts to
#'   aggregate. `NULL` uses the default assay. For BPCells sketch workflows,
#'   always pass `"RNA"` explicitly since pseudobulk aggregation requires
#'   raw counts from all cells. Default `NULL`.
#' @param min.cells Integer. Minimum number of cells a sample must contribute
#'   to a given cell type to be included in the pseudobulk for that cell type.
#'   Samples below this threshold are excluded to avoid noisy pseudo-replicates.
#'   Default `10`.
#' @param min.samples.per.group Integer. Minimum number of samples required
#'   per group level within a cell type. Cell types where any group has fewer
#'   than this many samples are skipped entirely. Default `3`.
#' @param run.ssgsea Logical. Also compute per-patient pathway activity via
#'   single-sample GSEA (ssGSEA, using `GSVA::gsva()`)? Produces score
#'   heatmaps and contrast boxplots per cell type per database. Requires the
#'   `GSVA` Bioconductor package. Default `TRUE`.
#' @param ssgsea.use.padj Logical. When `TRUE` (default), pathway significance
#'   in ssGSEA is judged by BH-adjusted p-value (`padj < 0.05`). Set to
#'   `FALSE` to use the raw p-value instead -- useful when many contrasts make
#'   BH correction too aggressive and real signals drop out. Affects pathway
#'   selection for heatmaps, boxplot panel selection, and significance bracket
#'   thresholds. Default `TRUE`.
#' @param top_n_heatmap Integer. Number of pathways shown in the NES heatmap
#'   per contrast/database combination, selected by maximum |NES| across cell
#'   types. Default `30`.
#' @param top_n_ssgsea Integer. Number of pathways shown in the ssGSEA score
#'   heatmap per cell type, selected by significance or variance. Default `30`.
#' @param nes.cutoff Numeric. Minimum absolute NES required for a pathway to
#'   appear in the **summary** NES heatmap (pathways x contrasts). The candidate
#'   row set is first built as the union of the top `top_n_heatmap` pathways
#'   from each individual per-contrast heatmap; `nes.cutoff` then removes rows
#'   whose maximum |NES| across all contrasts falls below this threshold.
#'   Default `1.0`. Set to `0` to skip the cutoff and show all union rows.
#' @param group_colors Named character vector mapping `contrast.by` levels to
#'   colors. Used consistently across all boxplots and heatmap annotations.
#'   `NULL` auto-generates colors from the Nour palette. Default `NULL`.
#' @param heatmap_params Named list of additional arguments forwarded to
#'   [ComplexHeatmap::Heatmap()]. Any default set internally can be overridden
#'   here (e.g. `list(clustering_distance_rows = "pearson")`).
#' @param heatmap_colors A `circlize::colorRamp2` color function for all NES
#'   and ssGSEA score heatmaps. `NULL` (default) uses a 5-stop blue-white-red
#'   diverging palette centered at 0. Breaks are fixed at `c(-2, -1, 0, 1, 2)`
#'   for ssGSEA row-z-score heatmaps and auto-scaled to the data for NES
#'   heatmaps. Supply any `colorRamp2` object to override all heatmaps at once.
#'   Examples:
#'   \itemize{
#'     \item \strong{Classic RdBu:}
#'       `circlize::colorRamp2(c(-3, 0, 3), c("#2166ac", "white", "#b2182b"))`
#'     \item \strong{Purple-green (PRGn):}
#'       `circlize::colorRamp2(c(-3, 0, 3), RColorBrewer::brewer.pal(3, "PRGn"))`
#'     \item \strong{Viridis plasma (for one-sided scores):}
#'       `circlize::colorRamp2(seq(0, 3, length.out = 9), viridis::plasma(9))`
#'     \item \strong{Tighter scale for subtle signals:}
#'       `circlize::colorRamp2(c(-1.5, -0.75, 0, 0.75, 1.5), c("#007dd1", "#b3d9f5", "white", "#f5c08a", "#ab3000"))`
#'   }
#' @param save.rds Logical. Save the pseudobulk counts, sample metadata, GSEA
#'   results, and ssGSEA results as `RunGSEA_pseudobulk_results.rds` in
#'   `output_dir`. Required for `redo_plots = TRUE`. Default `TRUE`.
#' @param redo_plots Logical. Reload a previously saved RDS and regenerate
#'   all plots without rerunning the full pipeline. Default `FALSE`. When
#'   `TRUE`:
#'   \itemize{
#'     \item `RunGSEA_pseudobulk_results.rds` is loaded from `output_dir`.
#'     \item NES heatmaps and lollipop plots are regenerated from the saved
#'       GSEA data frame.
#'     \item ssGSEA score CSVs on disk are reread and used to regenerate score
#'       heatmaps and boxplots.
#'     \item `seurat_object`, `group.by`, `sample.by`, and `contrast.by` are
#'       not needed and can be omitted.
#'   }
#'   Useful after changing `group_colors`, `heatmap_params`, or any plot
#'   parameter without rerunning the hours-long pipeline.
#' @param caffeinate Logical. Prevent the system from sleeping during the run
#'   (macOS only, uses the built-in `caffeinate` utility). Recommended for
#'   overnight runs. Default `FALSE`.
#'
#' @return Invisibly returns a named list:
#' \describe{
#'   \item{`pseudobulk`}{Genes x samples aggregated count matrix.}
#'   \item{`samples`}{Sample-level metadata data frame used in the model.}
#'   \item{`gsea`}{Long-format data frame of all fgsea results with columns
#'     `pathway`, `NES`, `padj`, `cell_type`, `contrast`, and `db`.}
#'   \item{`ssgsea`}{ssGSEA contrast results (one row per pathway/cell type/
#'     contrast), or `NULL` if `run.ssgsea = FALSE`.}
#'   \item{`contrasts`}{The named contrast list actually used in the model
#'     (useful when contrasts were auto-generated).}
#' }
#'
#' @seealso [RunGSEA()] for cell-level GSEA, [RunSCssGSEA()] for single-cell
#'   ssGSEA scoring, [PlotGSEAEnrichment()] to visualize enrichment plots from
#'   the output directory.
#'
#' @export
RunGSEA_pseudobulk <- function(seurat_object       = NULL,
                               group.by            = NULL,
                               sample.by           = NULL,
                               contrast.by         = NULL,
                               split.by            = NULL,
                               subset.by           = NULL,
                               run_label           = NULL,
                               label.by            = NULL,
                               pseudobulk          = NULL,
                               covariates          = NULL,
                               contrasts           = NULL,
                               # ── Gene set input (pick one; pathway_sets used when all NULL) ──
                               gene_sets            = NULL,      # named list of gene vectors
                               deg_df               = NULL,      # DE results data.frame
                               deg_gene_column      = "feature",
                               deg_group_column     = "group",
                               deg_fc_column        = "logFC",
                               deg_padj_column      = "padj",
                               deg_padj_cutoff      = 0.05,
                               deg_top_n            = 20L,
                               # deprecated aliases
                               identity_column     = NULL,
                               sample_column       = NULL,
                               group_column        = NULL,
                               pathway_sets     = list(
                                 Hallmark = list(category = "H"),
                                 KEGG     = list(category = "C2", subcategory = "CP:KEGG"),
                                 Reactome = list(category = "C2", subcategory = "CP:REACTOME"),
                                 WP       = list(category = "C2", subcategory = "CP:WIKIPATHWAYS")),
                               species          = "human",
                               output_dir,
                               assay            = NULL,
                               min.cells        = 10L,
                               min.samples.per.group = 3L,
                               run.ssgsea       = TRUE,
                               ssgsea.use.padj  = TRUE,
                               top_n_heatmap    = 30L,
                               top_n_ssgsea     = 30L,
                               nes.cutoff       = 1.0,
                               group_colors     = NULL,
                               heatmap_params   = list(row_names_side      = "left",
                                                       show_row_dend       = FALSE,
                                                       row_names_max_width = grid::unit(15, "cm")),
                               heatmap_colors   = NULL,
                               save.rds         = TRUE,
                               redo_plots       = FALSE,
                               caffeinate       = FALSE) {

  if (caffeinate) { .caff <- .nk_caffeinate(); on.exit(.nk_decaffeinate(.caff), add = TRUE) }

  # Backwards-compatibility aliases
  group.by    <- group.by    %||% identity_column
  sample.by   <- sample.by   %||% sample_column
  contrast.by <- contrast.by %||% group_column

  # subset.by: run the full pipeline once per level, saving results in subdirectories
  if (!is.null(subset.by) && !is.null(seurat_object)) {
    if (!subset.by %in% colnames(seurat_object@meta.data))
      stop("'subset.by' column '", subset.by, "' not found in meta.data.")
    sub_levels <- sort(unique(as.character(
      seurat_object@meta.data[[subset.by]][!is.na(seurat_object@meta.data[[subset.by]])])))
    message("subset.by = '", subset.by, "': running ", length(sub_levels),
            " separate analyses (", paste(sub_levels, collapse = ", "), ")")
    results <- lapply(sub_levels, function(sl) {
      message("\n--- subset.by '", subset.by, "' = '", sl, "' ---")
      keep    <- rownames(seurat_object@meta.data)[
        !is.na(seurat_object@meta.data[[subset.by]]) &
        seurat_object@meta.data[[subset.by]] == sl]
      sub_obj <- subset(seurat_object, cells = keep)
      tryCatch(
        RunGSEA_pseudobulk(
          seurat_object         = sub_obj,
          group.by              = group.by,
          sample.by             = sample.by,
          contrast.by           = contrast.by,
          split.by              = split.by,
          subset.by             = NULL,
          run_label             = paste0(subset.by, ": ", sl),
          label.by              = label.by,
          pseudobulk            = NULL,
          covariates            = covariates,
          contrasts             = contrasts,
          gene_sets             = gene_sets,
          deg_df                = deg_df,
          deg_gene_column       = deg_gene_column,
          deg_group_column      = deg_group_column,
          deg_fc_column         = deg_fc_column,
          deg_padj_column       = deg_padj_column,
          deg_padj_cutoff       = deg_padj_cutoff,
          deg_top_n             = deg_top_n,
          pathway_sets          = pathway_sets,
          species               = species,
          output_dir            = file.path(output_dir, make.names(sl), make.names(contrast.by)),
          assay                 = assay,
          min.cells             = min.cells,
          min.samples.per.group = min.samples.per.group,
          run.ssgsea            = run.ssgsea,
          ssgsea.use.padj       = ssgsea.use.padj,
          top_n_heatmap         = top_n_heatmap,
          top_n_ssgsea          = top_n_ssgsea,
          nes.cutoff            = nes.cutoff,
          group_colors          = group_colors,
          heatmap_params        = heatmap_params,
          save.rds              = save.rds,
          redo_plots            = redo_plots,
          caffeinate            = FALSE
        ),
        error = function(e) {
          message("  ERROR in subset '", sl, "': ", conditionMessage(e))
          NULL
        }
      )
    })
    names(results) <- make.names(sub_levels)
    n_ok <- sum(!vapply(results, is.null, logical(1)))
    message("subset.by complete: ", n_ok, "/", length(sub_levels), " subsets succeeded.")
    return(invisible(results))
  }

  use_custom <- !is.null(gene_sets) || !is.null(deg_df)
  for (pkg in c("limma", "edgeR", "fgsea", "ggplot2", "Matrix"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
  if (!use_custom && !requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required for MSigDB gene sets. ",
         "Install it with install.packages('msigdbr'), or supply gene_sets / deg_df.")
  if (run.ssgsea && !requireNamespace("GSVA", quietly = TRUE)) {
    warning("GSVA not available; skipping ssGSEA."); run.ssgsea <- FALSE
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ── redo_plots: reload saved results and re-run only visualisations ────────
  if (redo_plots) {
    rds_path <- file.path(output_dir, "RunGSEA_pseudobulk_results.rds")
    if (!file.exists(rds_path))
      stop("redo_plots = TRUE but no saved results found at:\n  ", rds_path,
           "\nRun RunGSEA_pseudobulk() without redo_plots first.")

    message("redo_plots: loading saved results from ", basename(rds_path), " ...")
    saved <- readRDS(rds_path)
    gsea  <- saved$gsea
    ssg   <- saved$ssgsea
    samp  <- saved$samples
    # User-provided contrasts take priority. Only fall back to the saved list
    # when the caller did not supply their own, so that re-running with a new
    # contrast list does not require deleting the RDS first.
    if (is.null(contrasts)) contrasts <- saved$contrasts

    # Resolve group_colors from the saved sample metadata (same logic as
    # .pb_ssgsea) so heatmap and boxplot always use identical colors.
    all_grps <- levels(droplevels(factor(samp$group[!is.na(samp$group)])))
    if (is.null(group_colors)) {
      group_colors <- stats::setNames(Nour_pal("all")(length(all_grps)), all_grps)
    } else {
      missing_grps <- all_grps[!all_grps %in% names(group_colors)]
      if (length(missing_grps) > 0) {
        extra <- stats::setNames(
          Nour_pal("all")(length(all_grps))[seq_along(missing_grps)],
          missing_grps)
        group_colors <- c(group_colors, extra)
      }
    }

    # ── NES heatmaps + lollipops ────────────────────────────────────────────
    if (!is.null(gsea) && nrow(gsea) > 0) {
      for (cn in unique(gsea$contrast)) {
        for (db in unique(gsea$db)) {
          sub <- gsea[gsea$contrast == cn & gsea$db == db, , drop = FALSE]
          if (!nrow(sub)) next
          db_dir <- file.path(output_dir, cn, db)
          dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)

          message("  NES heatmap: ", cn, " / ", db)
          .pb_nes_heatmap(sub, cn, db, db_dir, top_n_heatmap, heatmap_params,
                          heatmap_colors = heatmap_colors, run_label = run_label)

          for (ct in unique(sub$cell_type)) {
            fg <- sub[sub$cell_type == ct, , drop = FALSE]
            .pb_lollipop(fg, db, ct, cn, db_dir, run_label = run_label)
          }
        }
      }
      # Summary NES heatmap: pathways x contrasts, one per cell type x db.
      # Only generated when there are at least 2 contrasts.
      if (length(unique(gsea$contrast)) > 1L) {
        for (ct in unique(gsea$cell_type)) {
          for (db in unique(gsea$db)) {
            message("  NES summary heatmap: ", ct, " / ", db)
            .pb_nes_summary_heatmap(
              gsea[gsea$db == db, , drop = FALSE],
              ct             = ct,
              db             = db,
              output_dir     = output_dir,
              top_n          = top_n_heatmap,
              nes.cutoff     = nes.cutoff,
              heatmap_params = heatmap_params,
              heatmap_colors = heatmap_colors,
              run_label      = run_label)
          }
        }
      }
    }

    # ── ssGSEA heatmaps + boxplots ───────────────────────────────────────────
    sig_col <- if (isTRUE(ssgsea.use.padj)) "padj" else "p"
    if (run.ssgsea && !is.null(samp)) {
      cell_types <- sort(unique(samp$cell_type))
      dbs        <- if (!is.null(gsea) && nrow(gsea) > 0)
        unique(gsea$db) else names(pathway_sets)

      for (ct in cell_types) {
        for (db in dbs) {
          db_dir   <- file.path(output_dir, "ssGSEA", db)
          csv_path <- file.path(db_dir,
                                paste0("ssGSEA ", db, " ", make.names(ct), " scores.csv"))
          if (!file.exists(csv_path)) next

          message("  ssGSEA heatmap: ", ct, " / ", db)
          sc_df <- tryCatch(
            utils::read.csv(csv_path, check.names = FALSE, row.names = NULL),
            error = function(e) NULL)
          if (is.null(sc_df)) next

          sc          <- as.matrix(sc_df[, -1, drop = FALSE])
          rownames(sc) <- sc_df[[1]]

          ks <- rownames(samp)[samp$cell_type == ct]
          s  <- samp[intersect(ks, colnames(sc)), , drop = FALSE]
          if (nrow(s) < 2) next
          s$group <- droplevels(s$group)

          # Determine shared pathway set (same logic as in .pb_ssgsea)
          res_ct <- if (!is.null(ssg) && nrow(ssg) > 0)
            ssg[ssg$cell_type == ct, , drop = FALSE] else NULL
          sig_label_long_r <- if (isTRUE(ssgsea.use.padj)) "BH-adjusted p < 0.05" else "raw p < 0.05"
          sig_label_r      <- if (isTRUE(ssgsea.use.padj)) "padj" else "p (raw)"
          row_sel_r        <- NULL
          if (!is.null(res_ct) && nrow(res_ct) > 0 &&
              any(!is.na(res_ct[[sig_col]]) & res_ct[[sig_col]] < 0.05)) {
            sig_sub  <- res_ct[!is.na(res_ct[[sig_col]]) & res_ct[[sig_col]] < 0.05, ]
            sig_sub  <- sig_sub[order(-abs(sig_sub$diff)), ]
            show_pws <- unique(sig_sub$pathway)[
              seq_len(min(top_n_ssgsea, length(unique(sig_sub$pathway))))]
            row_sel_r <- paste0("top ", length(show_pws),
                                " significantly enriched pathways (", sig_label_long_r, ")")
          } else {
            row_var  <- apply(sc, 1, stats::var, na.rm = TRUE)
            show_pws <- names(sort(row_var, decreasing = TRUE))[
              seq_len(min(top_n_ssgsea, nrow(sc)))]
            row_sel_r <- paste0("top ", length(show_pws),
                                " pathways by score variance (no ", sig_label_r, " < 0.05 found)")
          }

          .ssgsea_heatmap(sc, s, db, ct, db_dir,
                          pathways        = show_pws,
                          group_colors    = group_colors,
                          top_n           = top_n_ssgsea,
                          heatmap_params  = heatmap_params,
                          heatmap_colors  = heatmap_colors,
                          run_label       = run_label,
                          contrast_label  = contrast.by,
                          row_selection   = row_sel_r)

          if (!is.null(res_ct) && nrow(res_ct) > 0)
            .ssgsea_boxplot(sc, res_ct, s, db, ct, db_dir,
                            pathways     = show_pws,
                            group_colors = group_colors,
                            top_n        = top_n_ssgsea,
                            contrasts    = contrasts,
                            use.padj     = ssgsea.use.padj)
        }
      }
    }

    message("redo_plots complete.")
    return(invisible(list(
      pseudobulk = saved$pseudobulk,
      samples    = samp,
      gsea       = gsea,
      ssgsea     = ssg
    )))
  }
  # ── end redo_plots ──────────────────────────────────────────────────────────

  if (is.null(assay) && !is.null(seurat_object)) assay <- SeuratObject::DefaultAssay(seurat_object)

  ## ---- write method params for create_analysis_pptx ----
  # RunGSEA and RunGSEA_pseudobulk are fundamentally different methods:
  # RunGSEA uses cell-level Wilcoxon (AUC pre-ranking via presto), while
  # RunGSEA_pseudobulk aggregates to donor/sample pseudobulk and uses
  # limma-voom moderated t-statistics as the ranking statistic - matching
  # standard bulk RNA-seq GSEA practice and respecting sample as the
  # statistical unit.
  {
    db_list_str  <- paste(names(pathway_sets), collapse = ", ")
    contrast_str <- if (!is.null(contrasts))
      paste(names(contrasts), collapse = "; ")
    else
      "auto-detected from contrast.by levels"

    .write_subdir_params(output_dir, list(
      date                    = format(Sys.Date()),
      gsea_method             = "pseudobulk (limma-voom)",
      gsea_ident_col          = group.by,
      gsea_sample_col         = sample.by,
      gsea_group_col          = contrast.by,
      gsea_databases          = names(pathway_sets),
      gsea_contrasts          = contrast_str,
      gsea_min_cells          = min.cells,
      gsea_min_samples        = min.samples.per.group,
      gsea_species            = species,
      gsea_nes_cutoff         = nes.cutoff,
      gsea_top_n_heatmap      = top_n_heatmap,
      gsea_run_ssgsea         = run.ssgsea,
      gsea_top_n_ssgsea       = top_n_ssgsea,
      gsea_ssgsea_sig_col     = if (isTRUE(ssgsea.use.padj)) "padj" else "p",
      methods_text       = paste0(
        "Donor-level gene set enrichment analysis was performed using a ",
        "pseudobulk approach. Raw counts (assay: '", assay %||% "default",
        "') were summed per biological sample ('", sample.by, "') within ",
        "each cell-type subset ('", group.by, "'), retaining only ",
        "sample x cell-type combinations with >= ", min.cells,
        " cells and groups with >= ", min.samples.per.group,
        " samples. Pseudobulk libraries were normalized using edgeR ",
        "(TMM normalization) and precision weights were estimated with ",
        "limma-voom. A group-means design was fitted and the following ",
        "contrasts were tested: ", contrast_str,
        if (!is.null(covariates) && length(covariates) > 0)
          paste0(". Covariates included in the model: ",
                 paste(covariates, collapse = ", "))
        else "",
        ". Genes were ranked by the moderated t-statistic from limma and ",
        "gene set enrichment was assessed using fgsea (Korotkevich et al., ",
        "2021) against the following MSigDB collections: ", db_list_str,
        ". NES > 0 indicates enrichment in the first term of each contrast.",
        if (nes.cutoff > 0)
          paste0(" Pathways with max |NES| < ", nes.cutoff,
                 " across all contrasts were excluded from summary heatmaps.")
        else "",
        if (run.ssgsea)
          paste0(" Per-patient pathway activity was additionally quantified ",
                 "using single-sample GSEA (ssGSEA via GSVA) on the ",
                 "pseudobulk log-CPM expression matrix. Differential pathway ",
                 "activity was tested with limma. Significance threshold: ",
                 if (isTRUE(ssgsea.use.padj)) "BH-adjusted p < 0.05."
                 else "raw p < 0.05 (BH correction not applied).")
        else ""
      )
    ))
  }

  ## ---- gene sets (query once) — three modes: custom list > deg_df > MSigDB ----
  if (use_custom) {
    gs_res    <- .build_gene_sets_sc(
      gene_sets        = gene_sets,
      deg_df           = deg_df,
      deg_gene_column  = deg_gene_column,
      deg_group_column = deg_group_column,
      deg_fc_column    = deg_fc_column,
      deg_padj_column  = deg_padj_column,
      deg_padj_cutoff  = deg_padj_cutoff,
      deg_top_n        = deg_top_n
    )
    fgsea_dbs <- list(DEG = gs_res$gene_sets)
  } else {
    fgsea_dbs <- lapply(names(pathway_sets), function(db) {
      ps <- pathway_sets[[db]]
      mdf <- .msigdbr_get(species = species, category = ps$category,
                          subcategory = ps$subcategory)
      split(mdf$gene_symbol, mdf$gs_name)
    })
    names(fgsea_dbs) <- names(pathway_sets)
  }

  ## ---- pseudobulk: reuse precomputed, or aggregate (sum counts per sample x ident) ----
  original_contrast_levels <- NULL
  original_split_levels    <- NULL
  contrast_src_levels      <- NULL   # factor levels from contrast.by in Seurat metadata
  if (!is.null(pseudobulk)) {
    message("Reusing precomputed pseudobulk (skipping aggregation).")
    agg  <- as.matrix(pseudobulk$counts)
    samp <- pseudobulk$samples
    if (is.null(samp$.pb)) samp$.pb <- colnames(agg)
    if (is.null(samp$sample) || is.null(samp$cell_type)) {
      pb_parts <- strsplit(samp$.pb, "\\|\\|", fixed = FALSE)
      n_parts  <- lengths(pb_parts)
      if (all(n_parts >= 3L)) {
        # New 3-part format: sample || cell_type || group
        if (is.null(samp$sample))    samp$sample    <- vapply(pb_parts, `[`, character(1L), 1L)
        if (is.null(samp$cell_type)) samp$cell_type <- vapply(pb_parts, `[`, character(1L), 2L)
        if (is.null(samp$group))     samp$group     <- vapply(pb_parts, `[`, character(1L), 3L)
      } else {
        # Legacy 2-part format: sample || cell_type (group from metadata)
        if (is.null(samp$sample))    samp$sample    <- sub("\\|\\|.*", "", samp$.pb)
        if (is.null(samp$cell_type)) samp$cell_type <- sub(".*\\|\\|", "", samp$.pb)
        if (is.null(samp$group) && !is.null(contrast.by) && contrast.by %in% colnames(samp))
          samp$group <- samp[[contrast.by]]
      }
    }
    rownames(samp) <- samp$.pb; samp <- samp[colnames(agg), ]
  } else {
    if (is.null(seurat_object)) stop("Provide either `seurat_object` or a precomputed `pseudobulk`.")
    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(seurat_object)
    md <- seurat_object@meta.data

    # group.by: when NULL, treat all cells as a single cell type ("All").
    # This is the correct behavior for single-cell-type objects (e.g. Microglia
    # only) where the user wants to compare contrast.by groups without any
    # further cell-type splitting.
    if (is.null(group.by)) {
      md$.cell_type_all <- "All"
      seurat_object@meta.data$.cell_type_all <- "All"
      group.by <- ".cell_type_all"
      message("group.by not supplied - treating all cells as one cell type ('All').")
    } else if (!group.by %in% colnames(md)) {
      stop("'group.by' column '", group.by, "' not found in meta.data.")
    }

    # Save factor levels from contrast.by before split.by may rename the column.
    # Used later to re-level samp$group so heatmap columns follow the user's
    # preset factor order rather than alphabetical ordering.
    if (!is.null(contrast.by) && contrast.by %in% colnames(md) &&
        is.factor(md[[contrast.by]]))
      contrast_src_levels <- levels(md[[contrast.by]])

    # split.by: create an interaction column so contrasts operate on group x split combos
    if (!is.null(split.by)) {
      if (!split.by %in% colnames(md))
        stop("'split.by' column '", split.by, "' not found in meta.data.")
      original_contrast_levels <- sort(unique(as.character(md[[contrast.by]][!is.na(md[[contrast.by]])])))
      original_split_levels    <- sort(unique(as.character(md[[split.by]][!is.na(md[[split.by]])])))
      md$.contrast_split <- paste(md[[contrast.by]], md[[split.by]], sep = ".")
      contrast.by <- ".contrast_split"
      run_label <- if (!is.null(run_label))
        paste0(run_label, " | split by ", split.by)
      else
        paste0("split by ", split.by)
    }

    md <- md[!is.na(md[[contrast.by]]) & !is.na(md[[group.by]]) & !is.na(md[[sample.by]]), ]
    seurat_object <- subset(seurat_object, cells = rownames(md))
    # contrast.by is included in the .pb key so cell-level groupings (e.g.
    # STAT3_split, SplitByGene labels) produce correctly separated pseudobulks
    seurat_object$.pb <- paste(seurat_object@meta.data[[sample.by]],
                               seurat_object@meta.data[[group.by]],
                               seurat_object@meta.data[[contrast.by]], sep = "||")
    message("Aggregating pseudobulk (", assay, " counts) by ",
            sample.by, " x ", group.by, " x ", contrast.by, " ...")
    agg <- Seurat::AggregateExpression(seurat_object, assays = assay, group.by = ".pb", slot = "counts")[[assay]]
    agg <- as.matrix(agg)
    cc <- as.data.frame(table(seurat_object$.pb)); colnames(cc) <- c(".pb", "n_cells")
    samp <- data.frame(.pb = colnames(agg), stringsAsFactors = FALSE)
    # Parse 3-part key: sample || cell_type || contrast_group
    pb_parts       <- strsplit(samp$.pb, "\\|\\|", fixed = FALSE)
    samp$sample    <- vapply(pb_parts, `[`, character(1L), 1L)
    samp$cell_type <- vapply(pb_parts, `[`, character(1L), 2L)
    samp$group     <- vapply(pb_parts, `[`, character(1L), 3L)
    # Join donor-level covariates (group is now derived from the key, not metadata lookup)
    if (length(covariates) > 0L) {
      sm_cov <- md[!duplicated(md[[sample.by]]), c(sample.by, covariates), drop = FALSE]
      colnames(sm_cov)[1L] <- "sample"
      samp <- merge(samp, sm_cov, by = "sample", all.x = TRUE)
    }
    samp <- merge(samp, cc, by = ".pb", all.x = TRUE)
    rownames(samp) <- samp$.pb; samp <- samp[colnames(agg), ]
  }
  samp$group <- droplevels(factor(samp$group))
  # Re-apply the user's original factor level order (saved above from the
  # Seurat object's contrast.by column). droplevels(factor(...)) above would
  # otherwise impose alphabetical ordering, which breaks heatmap column order.
  if (!is.null(contrast_src_levels)) {
    cur <- levels(samp$group)
    samp$group <- factor(samp$group, levels = c(
      contrast_src_levels[contrast_src_levels %in% cur],
      setdiff(cur, contrast_src_levels)))
  }
  glevels <- levels(samp$group)

  # covariate cleaning
  for (cv in covariates) {
    x <- samp[[cv]]
    if (is.numeric(x)) { x[is.na(x)] <- stats::median(x, na.rm = TRUE); samp[[cv]] <- as.numeric(scale(x)) }
    else { x <- as.character(x); x[is.na(x)] <- "unknown"; samp[[cv]] <- factor(x) }
  }

  ## ---- contrasts ----
  if (is.null(contrasts)) {
    contrasts <- if (!is.null(original_contrast_levels))
      .auto_contrasts_split(original_contrast_levels, original_split_levels)
    else
      .auto_contrasts_simple(glevels)
  }
  # Translate simple c("A","B") pairs to limma coefficient strings
  contrasts <- .resolve_contrasts(contrasts)
  message("Contrasts: ", paste(names(contrasts), collapse = ", "))

  ## ---- per cell type: limma-voom + fgsea ----
  idents <- sort(unique(samp$cell_type))
  # Respect factor levels from source metadata so heatmap columns follow the
  # same order the user already set on the Seurat object.
  if (!is.null(seurat_object) && !is.null(group.by) &&
      group.by %in% colnames(seurat_object@meta.data) &&
      is.factor(seurat_object@meta.data[[group.by]])) {
    src_levels <- levels(seurat_object@meta.data[[group.by]])
    idents <- c(src_levels[src_levels %in% idents], setdiff(idents, src_levels))
  }
  gsea_all <- list()           # long NES table
  for (ct in idents) {
    keep_s <- rownames(samp)[samp$cell_type == ct & (is.na(samp$n_cells) | samp$n_cells >= min.cells)]
    s <- samp[keep_s, , drop = FALSE]; s$group <- droplevels(s$group)
    if (any(table(s$group) < min.samples.per.group) || nlevels(s$group) < 2) {
      message("  skip ", ct, " (insufficient samples per group)"); next
    }
    y <- edgeR::DGEList(counts = round(agg[, keep_s, drop = FALSE]))
    keepg <- edgeR::filterByExpr(y, group = s$group); y <- y[keepg, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    form <- stats::as.formula(paste("~ 0 + group", if (length(covariates)) paste("+", paste(covariates, collapse = " + ")) else ""))
    design <- stats::model.matrix(form, data = s); colnames(design) <- make.names(colnames(design))
    v <- limma::voom(y, design)
    fit <- limma::lmFit(v, design)
    # only keep contrasts whose coefficients all exist for this cell type
    valid_mask <- vapply(contrasts, function(ex) .coef_in(ex, colnames(design)), logical(1))
    valid      <- contrasts[valid_mask]
    dropped    <- names(contrasts)[!valid_mask]
    if (length(dropped)) {
      # Show which coefficients are actually missing so the user can fix their cf() calls
      missing_by <- vapply(dropped, function(nm) {
        toks <- unique(regmatches(contrasts[[nm]],
                                  gregexpr("group[A-Za-z0-9._]+", contrasts[[nm]])[[1]]))
        bad <- toks[!toks %in% colnames(design)]
        if (length(bad)) paste(bad, collapse = ", ") else "?"
      }, character(1))
      message("  dropped contrasts (coefficient not in design):")
      for (nm in dropped)
        message("    ", nm, "  [missing: ", missing_by[[nm]], "]")
      message("  available coefficients: ", paste(grep("^group", colnames(design), value = TRUE), collapse = ", "))
    }
    if (!length(valid)) { message("  skip ", ct, " (no valid contrasts)"); next }
    cm <- limma::makeContrasts(contrasts = unlist(valid), levels = design)
    colnames(cm) <- names(valid)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
    for (cn in colnames(cm)) {
      tt <- limma::topTable(fit2, coef = cn, number = Inf, sort.by = "none")
      ranks <- sort(stats::setNames(tt$t, rownames(tt)), decreasing = TRUE)
      for (db in names(fgsea_dbs)) {
        fg <- tryCatch(fgsea::fgsea(fgsea_dbs[[db]], ranks, minSize = 10, maxSize = 500),
                       error = function(e) NULL)
        if (is.null(fg) || !nrow(fg)) next
        fg <- as.data.frame(fg[, c("pathway","pval","padj","NES","size")])
        fg$cell_type <- ct; fg$contrast <- cn; fg$db <- db
        gsea_all[[paste(ct,cn,db)]] <- fg
        # per cell-type/contrast/db CSV + lollipop, mirroring RunGSEA layout
        db_dir <- file.path(output_dir, cn, db); dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)
        csv_label <- if (!is.null(run_label)) paste0(" [", make.names(run_label), "]") else ""
        utils::write.csv(fg[order(-fg$NES), ],
                         file.path(db_dir, paste0(db, " ", make.names(ct), " ", cn, csv_label, ".csv")), row.names = FALSE)
        .pb_lollipop(fg, db, ct, cn, db_dir, run_label = run_label)
      }
    }
    message("  ", ct, ": ", length(valid), " contrasts x ", length(fgsea_dbs), " DBs")
  }
  gsea <- if (length(gsea_all)) do.call(rbind, gsea_all) else data.frame()
  # Set cell_type and contrast as ordered factors so tapply, redo_plots, and
  # the summary heatmap all respect the user's intended ordering.
  if (nrow(gsea) > 0L) {
    gsea$cell_type <- factor(gsea$cell_type, levels = idents)
    contrast_order <- names(contrasts)[names(contrasts) %in% unique(gsea$contrast)]
    gsea$contrast  <- factor(gsea$contrast, levels = contrast_order)
  }

  ## ---- per (contrast, db) NES heatmaps: pathways x cell types ----
  if (nrow(gsea)) {
    for (cn in levels(gsea$contrast)) for (db in names(fgsea_dbs)) {
      sub <- gsea[gsea$contrast == cn & gsea$db == db, ]
      if (!nrow(sub)) next
      .pb_nes_heatmap(sub, cn, db, file.path(output_dir, cn, db), top_n_heatmap, heatmap_params,
                      heatmap_colors = heatmap_colors, run_label = run_label)
    }
    ## ---- summary NES heatmap: pathways x contrasts, per cell type per db ----
    if (length(unique(gsea$contrast)) > 1L) {
      for (ct in levels(gsea$cell_type)) for (db in names(fgsea_dbs)) {
        message("  NES summary heatmap: ", ct, " / ", db)
        .pb_nes_summary_heatmap(
          gsea[gsea$db == db, , drop = FALSE],
          ct             = ct,
          db             = db,
          output_dir     = output_dir,
          top_n          = top_n_heatmap,
          nes.cutoff     = nes.cutoff,
          heatmap_params = heatmap_params,
          heatmap_colors = heatmap_colors,
          run_label      = run_label)
      }
    }
  }

  ## ---- per-patient ssGSEA activity ----
  ssg <- NULL
  if (run.ssgsea) ssg <- .pb_ssgsea(agg, samp, fgsea_dbs, contrasts, covariates,
                                    output_dir, min.cells, min.samples.per.group,
                                    group_colors   = group_colors,
                                    top_n          = top_n_ssgsea,
                                    heatmap_params = heatmap_params,
                                    heatmap_colors = heatmap_colors,
                                    run_label      = run_label,
                                    use.padj       = ssgsea.use.padj)

  if (save.rds) saveRDS(
    list(pseudobulk = agg, samples = samp, gsea = gsea, ssgsea = ssg, contrasts = contrasts),
    file.path(output_dir, "RunGSEA_pseudobulk_results.rds"))
  message("RunGSEA_pseudobulk complete.")
  invisible(list(pseudobulk = agg, samples = samp, gsea = gsea, ssgsea = ssg))
}

# ---- helpers ---------------------------------------------------------------

#' @keywords internal
.auto_contrasts_simple <- function(glevels) {
  if (length(glevels) < 2) return(list())
  ct <- list()
  for (i in seq_len(length(glevels) - 1)) {
    for (j in (i + 1):length(glevels)) {
      la <- glevels[i]; lb <- glevels[j]
      ct[[paste0(la, "_vs_", lb)]] <- paste0(
        make.names(paste0("group", la)), " - ", make.names(paste0("group", lb)))
    }
  }
  ct
}

#' @keywords internal
.auto_contrasts_split <- function(contrast_levels, split_levels) {
  ct <- list()
  sep <- "."
  for (s in split_levels) {
    for (i in seq_len(length(contrast_levels) - 1)) {
      for (j in (i + 1):length(contrast_levels)) {
        la <- contrast_levels[i]; lb <- contrast_levels[j]
        ca <- make.names(paste0("group", la, sep, s))
        cb <- make.names(paste0("group", lb, sep, s))
        ct[[paste0(la, "_vs_", lb, "_", s)]] <- paste0(ca, " - ", cb)
      }
    }
  }
  # Interaction term for the 2×2 case
  if (length(contrast_levels) == 2 && length(split_levels) == 2) {
    la <- contrast_levels[1]; lb <- contrast_levels[2]
    s1 <- split_levels[1];    s2 <- split_levels[2]
    coef <- function(l, s) make.names(paste0("group", l, sep, s))
    ct[["Interaction"]] <- paste0(
      "(", coef(la, s1), " - ", coef(lb, s1), ") - ",
      "(", coef(la, s2), " - ", coef(lb, s2), ")")
  }
  ct
}

#' @keywords internal
# Each element of contrasts must be a single limma expression string built from
# coefficient names of the form make.names(paste0("group", level)).
# Use a helper like cf <- function(l) make.names(paste0("group", l)) to build
# coefficient names, then combine them into expressions:
#   contrasts = list(
#     AD_vs_Ctrl = paste(cf("AD"), "-", cf("Control")),
#     Interaction = paste0("(", cf("AD.F"), "-", cf("Ctrl.F"), ") - (",
#                                cf("AD.M"), "-", cf("Ctrl.M"), ")")
#   )
.resolve_contrasts <- function(contrasts) {
  lapply(seq_along(contrasts), function(i) {
    x   <- contrasts[[i]]
    nm  <- names(contrasts)[i]
    tag <- if (!is.null(nm) && nzchar(nm)) paste0("contrast '", nm, "'") else
             paste0("contrast [[", i, "]]")
    if (!is.character(x) || length(x) != 1L)
      stop(tag, " must be a single limma expression string. ",
           "Build coefficient names with make.names(paste0('group', level)), e.g.:\n",
           "  cf <- function(l) make.names(paste0('group', l))\n",
           "  list(my_contrast = paste(cf('GroupA'), '-', cf('GroupB')))")
    if (!grepl("group[A-Za-z0-9._]+", x))
      warning(tag, ": '", x, "' does not contain any group* coefficient tokens. ",
              "Did you mean make.names(paste0('group', '", x, "'))? ",
              "Passing through as-is.")
    x
  }) |> stats::setNames(names(contrasts))
}

#' @keywords internal
.coef_in <- function(expr, coefs) {
  toks <- unique(regmatches(expr, gregexpr("group[A-Za-z0-9._]+", expr))[[1]])
  all(toks %in% coefs)
}

#' @keywords internal
.pb_lollipop <- function(fg, db, ct, cn, db_dir, n = 15, run_label = NULL) {
  fg <- fg[fg$padj < 0.25, , drop = FALSE]; if (!nrow(fg)) return(invisible())
  fg <- fg[order(-abs(fg$NES)), ][seq_len(min(n, nrow(fg))), ]
  fg$pw <- gsub("_", " ", sub("^[A-Z]+_", "", fg$pathway))
  plot_title <- paste0(db, " ", ct, " : ", cn)
  if (!is.null(run_label)) plot_title <- paste0(plot_title, "\n", run_label)
  p <- ggplot2::ggplot(fg, ggplot2::aes(stats::reorder(pw, NES), NES, color = NES > 0)) +
    ggplot2::geom_segment(ggplot2::aes(xend = pw, y = 0, yend = NES)) +
    ggplot2::geom_point(size = 3) + ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#B5232E", `FALSE` = "#2C5F8A"), guide = "none") +
    ggplot2::labs(title = plot_title, x = NULL, y = "NES") +
    ggplot2::theme_minimal(base_size = 10)
  lbl_part <- if (!is.null(run_label)) paste0(" [", make.names(run_label), "]") else ""
  fp <- file.path(db_dir, paste0(db, " - ", make.names(ct), " ", cn, lbl_part, " lollipop.pdf"))
  ggplot2::ggsave(fp, p, width = 6, height = 4.5)
  if (exists(".write_legend_sidecar")) .write_legend_sidecar(fp, paste0(
    "Lollipop of top ", nrow(fg), " gene sets (", db, ") for ", ct, ", contrast ", cn,
    " (donor-level pseudobulk, limma-voom moderated t pre-ranking, fgsea). NES>0 = enriched in the first term of the contrast."))
}

#' @keywords internal
.pb_nes_heatmap <- function(sub, cn, db, db_dir, top_n, heatmap_params = list(),
                            heatmap_colors = NULL, run_label = NULL) {
  dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)

  # Build NES and padj matrices (pathways × cell types)
  m  <- as.matrix(tapply(sub$NES,  list(sub$pathway, sub$cell_type),
                         function(x) x[1]))
  pm <- as.matrix(tapply(sub$padj, list(sub$pathway, sub$cell_type),
                         function(x) x[1]))
  m[is.na(m)] <- 0

  # Keep top top_n rows by max |NES|
  top <- names(sort(apply(abs(m), 1, max), decreasing = TRUE))[
    seq_len(min(top_n, nrow(m)))]
  m  <- m[top,  , drop = FALSE]
  pm <- pm[top, , drop = FALSE]

  # Reorder columns to match cell_type factor levels when present.
  # tapply already respects factor order, but this guards the redo_plots path
  # and any case where cell_type arrives as plain character.
  if (is.factor(sub$cell_type)) {
    desired  <- levels(sub$cell_type)
    existing <- desired[desired %in% colnames(m)]
    if (length(existing) > 0) {
      m  <- m[, existing, drop = FALSE]
      pm <- pm[, existing, drop = FALSE]
    }
  }

  # Significance star matrix (used in cell_fun)
  stars_mat <- ifelse(
    is.na(pm), "",
    ifelse(pm < 0.001, "***",
    ifelse(pm < 0.01,  "**",
    ifelse(pm < 0.05,  "*",  ""))))

  # Clean up pathway names (remove database prefix, underscores)
  rn           <- gsub("_", " ", sub("^[A-Z]+_", "", rownames(m)))
  rownames(m)  <- rn
  rownames(stars_mat) <- rn

  # cell_fun draws the significance stars inside each cell
  sf <- stars_mat   # captured in closure
  cell_fn <- function(j, i, x, y, width, height, fill) {
    s <- sf[i, j]
    if (!is.null(s) && nzchar(s))
      grid::grid.text(s, x, y,
                      gp = grid::gpar(fontsize = 8, col = "black"))
  }

  lbl_part  <- if (!is.null(run_label)) paste0(" [", make.names(run_label), "]") else ""
  fp        <- file.path(db_dir, paste0("GSEA ", db, " ", cn, lbl_part, " (NES)-Heatmap.pdf"))
  ht_title  <- paste(db, cn, "(NES, pseudobulk)")
  if (!is.null(run_label)) ht_title <- paste0(ht_title, "\n", run_label)

  # Use .gsea_ht() for consistent sizing and color, with stars via cell_fun
  .gsea_ht(
    mat            = m,
    title          = ht_title,
    filepath       = fp,
    heatmap_params = utils::modifyList(
      list(cell_fun = cell_fn),   # add stars by default
      heatmap_params              # user can override cell_fun or anything else
    ),
    heatmap_colors = heatmap_colors
  )

  if (exists(".write_legend_sidecar")) .write_legend_sidecar(fp, paste0(
    "Heatmap of pseudobulk GSEA Normalized Enrichment Scores (", db,
    ") for contrast ", cn, " across cell types. ",
    "Rows = top ", nrow(m), " gene sets by max |NES| across cell types. ",
    "Method: donor-level pseudobulk, edgeR TMM normalization, ",
    "limma-voom moderated t-statistic pre-ranking, fgsea. ",
    "NES > 0 = enriched in the first term of the contrast. ",
    "Overlaid significance stars: * BH < 0.05, ** < 0.01, *** < 0.001."))
}

# ---------------------------------------------------------------------------
# .pb_nes_summary_heatmap  - pathways x contrasts summary (one per cell type)
# ---------------------------------------------------------------------------
# Rows: union of the top top_n pathways from each individual per-contrast
#       heatmap, filtered to max |NES| >= nes.cutoff across all contrasts.
# Cols: all contrast names present for this cell type.
# Saves to output_dir/{db}/GSEA {db} {ct} NES summary.pdf
#' @keywords internal
.pb_nes_summary_heatmap <- function(gsea_db, ct, db, output_dir,
                                     top_n          = 30L,
                                     nes.cutoff     = 1.0,
                                     heatmap_params = list(),
                                     heatmap_colors = NULL,
                                     run_label      = NULL) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) return(invisible())

  sub <- gsea_db[gsea_db$cell_type == ct, , drop = FALSE]
  if (!nrow(sub)) return(invisible())

  contrasts_present <- unique(as.character(sub$contrast))
  if (length(contrasts_present) < 2L) return(invisible())

  # Full pathway x contrast NES and padj matrices
  m  <- as.matrix(tapply(sub$NES,  list(sub$pathway, sub$contrast), function(x) x[1]))
  pm <- as.matrix(tapply(sub$padj, list(sub$pathway, sub$contrast), function(x) x[1]))
  m[is.na(m)] <- 0

  # Step 1: union of top_n pathways from each contrast column (mirrors the
  # individual per-contrast heatmap row selection exactly)
  union_pws <- unique(unlist(lapply(contrasts_present, function(cn) {
    if (!cn %in% colnames(m)) return(character(0))
    v <- abs(m[, cn])
    names(sort(v, decreasing = TRUE))[seq_len(min(top_n, length(v)))]
  })))

  m_sub  <- m[union_pws,  , drop = FALSE]
  pm_sub <- pm[union_pws, , drop = FALSE]

  # Step 2: NES cutoff - drop rows where no contrast reaches |NES| >= nes.cutoff
  if (nes.cutoff > 0) {
    keep   <- apply(abs(m_sub), 1, max, na.rm = TRUE) >= nes.cutoff
    m_sub  <- m_sub[keep,  , drop = FALSE]
    pm_sub <- pm_sub[keep, , drop = FALSE]
  }

  if (!nrow(m_sub)) {
    message("    NES summary: 0 pathways pass |NES| >= ", nes.cutoff,
            " for ", ct, " / ", db, " - skipping")
    return(invisible())
  }

  # Respect contrast order from factor levels if present (mirrors column order
  # in individual heatmaps); fall back to alphabetical
  if (is.factor(sub$contrast)) {
    desired  <- levels(sub$contrast)
    existing <- desired[desired %in% colnames(m_sub)]
    if (length(existing)) { m_sub <- m_sub[, existing, drop = FALSE]; pm_sub <- pm_sub[, existing, drop = FALSE] }
  }

  # Clean pathway names (same transform as individual heatmaps)
  rn <- gsub("_", " ", sub("^[A-Z]+_", "", rownames(m_sub)))
  rownames(m_sub)  <- rn
  rownames(pm_sub) <- rn

  # Stars matrix
  stars_mat <- ifelse(is.na(pm_sub), "",
               ifelse(pm_sub < 0.001, "***",
               ifelse(pm_sub < 0.01,  "**",
               ifelse(pm_sub < 0.05,  "*", ""))))
  rownames(stars_mat) <- rn

  sf <- stars_mat
  cell_fn <- function(j, i, x, y, width, height, fill) {
    s <- sf[i, j]
    if (!is.null(s) && nzchar(s))
      grid::grid.text(s, x, y, gp = grid::gpar(fontsize = 8, col = "black"))
  }

  ht_title  <- paste(db, ct, "- NES across contrasts")
  lbl_part  <- if (!is.null(run_label)) paste0(" [", make.names(run_label), "]") else ""
  if (!is.null(run_label)) ht_title <- paste0(ht_title, "\n", run_label)

  db_dir <- file.path(output_dir, db)
  dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)
  fp <- file.path(db_dir,
                  paste0("GSEA ", db, " ", make.names(ct), lbl_part, " NES summary.pdf"))

  .gsea_ht(
    mat            = m_sub,
    title          = ht_title,
    filepath       = fp,
    heatmap_params = utils::modifyList(
      list(cell_fun           = cell_fn,
           column_names_rot   = 45,
           row_names_max_width = grid::unit(15, "cm")),
      heatmap_params
    ),
    heatmap_colors = heatmap_colors
  )

  if (exists(".write_legend_sidecar")) .write_legend_sidecar(fp, paste0(
    "Summary NES heatmap for ", ct, " cells (", db, " database). ",
    "Rows = union of top ", top_n, " pathways per contrast",
    if (nes.cutoff > 0) paste0(", filtered to max |NES| >= ", nes.cutoff) else "", ". ",
    "Columns = all contrasts tested. ",
    "Stars: * BH < 0.05, ** < 0.01, *** < 0.001."))

  invisible(m_sub)
}

#' @keywords internal
.pb_ssgsea <- function(agg, samp, fgsea_dbs, contrasts, covariates, output_dir,
                       min.cells, min.smp,
                       group_colors   = NULL,
                       top_n          = 30L,
                       heatmap_params = list(),
                       heatmap_colors = NULL,
                       run_label      = NULL,
                       use.padj       = TRUE) {

  # Resolve group_colors once so the heatmap and boxplot always use identical
  # colors regardless of which helper is called first or whether the user
  # supplied a partial/NULL vector.
  all_grps <- levels(droplevels(factor(samp$group[!is.na(samp$group)])))
  if (is.null(group_colors)) {
    group_colors <- stats::setNames(Nour_pal("all")(length(all_grps)), all_grps)
  } else {
    # Fill any missing levels with Nour_pal so the vector is always complete
    missing_grps <- all_grps[!all_grps %in% names(group_colors)]
    if (length(missing_grps) > 0) {
      extra <- stats::setNames(
        Nour_pal("all")(length(all_grps))[seq_along(missing_grps)],
        missing_grps)
      group_colors <- c(group_colors, extra)
    }
  }

  out <- list()
  for (ct in sort(unique(samp$cell_type))) {
    ks <- rownames(samp)[samp$cell_type == ct &
                         (is.na(samp$n_cells) | samp$n_cells >= min.cells)]
    s <- samp[ks, , drop = FALSE]; s$group <- droplevels(s$group)
    if (any(table(s$group) < min.smp) || nlevels(s$group) < 2) next

    y      <- edgeR::DGEList(round(agg[, ks, drop = FALSE]))
    y      <- edgeR::calcNormFactors(y)
    logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)

    for (db in names(fgsea_dbs)) {
      sc <- tryCatch({
        gp <- tryCatch(GSVA::ssgseaParam(logcpm, fgsea_dbs[[db]]),
                       error = function(e) NULL)
        if (!is.null(gp)) GSVA::gsva(gp) else
          GSVA::gsva(logcpm, fgsea_dbs[[db]], method = "ssgsea", verbose = FALSE)
      }, error = function(e) NULL)
      if (is.null(sc)) next

      db_dir <- file.path(output_dir, "ssGSEA", db)
      dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)

      # ── CSV: raw scores ──────────────────────────────────────────────────
      lbl_part <- if (!is.null(run_label)) paste0(" [", make.names(run_label), "]") else ""
      utils::write.csv(
        data.frame(pathway = rownames(sc), sc, check.names = FALSE),
        file.path(db_dir, paste0("ssGSEA ", db, " ", make.names(ct), lbl_part, " scores.csv")),
        row.names = FALSE)

      # ── Limma contrasts first - results determine which pathways to show ─
      # Running contrasts before the heatmap lets both figures show the same
      # pathway set (significant pathways, or top-by-variance as fallback).
      res    <- NULL
      form   <- stats::as.formula(
        paste("~ 0 + group",
              if (length(covariates))
                paste("+", paste(covariates, collapse = " + "))
              else ""))
      design <- stats::model.matrix(form, data = s)
      colnames(design) <- make.names(colnames(design))
      valid_mask <- vapply(contrasts, function(ex) .coef_in(ex, colnames(design)), logical(1))
      valid      <- contrasts[valid_mask]
      dropped    <- names(contrasts)[!valid_mask]
      if (length(dropped)) {
        missing_by <- vapply(dropped, function(nm) {
          toks <- unique(regmatches(contrasts[[nm]],
                                    gregexpr("group[A-Za-z0-9._]+", contrasts[[nm]])[[1]]))
          bad <- toks[!toks %in% colnames(design)]
          if (length(bad)) paste(bad, collapse = ", ") else "?"
        }, character(1))
        message("  ssGSEA dropped contrasts (coefficient not in design):")
        for (nm in dropped)
          message("    ", nm, "  [missing: ", missing_by[[nm]], "]")
      }

      if (length(valid)) {
        cm  <- limma::makeContrasts(contrasts = unlist(valid), levels = design)
        colnames(cm) <- names(valid)
        fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(sc, design), cm))
        res <- do.call(rbind, lapply(colnames(cm), function(cn) {
          tt <- limma::topTable(fit, coef = cn, number = Inf, sort.by = "none")
          data.frame(pathway = rownames(tt), cell_type = ct, contrast = cn,
                     diff    = tt$logFC, p = tt$P.Value, padj = tt$adj.P.Val)
        }))
        out[[paste(ct, db)]] <- res
        utils::write.csv(
          res,
          file.path(db_dir, paste0("ssGSEA ", db, " ", make.names(ct), " contrasts.csv")),
          row.names = FALSE)
      }

      # ── Shared pathway set: significant > variance fallback ──────────────
      # Significant pathways (any contrast, p/padj < 0.05) take priority so
      # the heatmap and boxplot show exactly the same rows.
      # If nothing is significant, fall back to top N by variance so the
      # heatmap is still informative even without reaching significance.
      sig_col   <- if (isTRUE(use.padj)) "padj" else "p"
      sig_label <- if (isTRUE(use.padj)) "padj" else "p (raw)"
      sig_label_long <- if (isTRUE(use.padj)) "BH-adjusted p < 0.05" else "raw p < 0.05"
      row_sel   <- NULL
      if (!is.null(res) && any(!is.na(res[[sig_col]]) & res[[sig_col]] < 0.05)) {
        sig_sub  <- res[!is.na(res[[sig_col]]) & res[[sig_col]] < 0.05, ]
        sig_sub  <- sig_sub[order(-abs(sig_sub$diff)), ]
        show_pws <- unique(sig_sub$pathway)[
          seq_len(min(top_n, length(unique(sig_sub$pathway))))]
        row_sel  <- paste0("top ", length(show_pws),
                           " significantly enriched pathways (", sig_label_long, ")")
      } else {
        row_var  <- apply(sc, 1, stats::var, na.rm = TRUE)
        show_pws <- names(sort(row_var, decreasing = TRUE))[
          seq_len(min(top_n, nrow(sc)))]
        row_sel  <- paste0("top ", length(show_pws),
                           " pathways by score variance (no ", sig_label, " < 0.05 found)")
        if (!is.null(res))
          message("    No significant ssGSEA pathways for ", ct, "/", db,
                  " (", sig_label, " < 0.05) - heatmap shows top by variance.")
      }

      # ── Heatmap and boxplot - same pathway set ───────────────────────────
      .ssgsea_heatmap(sc, s, db, ct, db_dir,
                      pathways        = show_pws,
                      group_colors    = group_colors,
                      top_n           = top_n,
                      heatmap_params  = heatmap_params,
                      heatmap_colors  = heatmap_colors,
                      run_label       = run_label,
                      row_selection   = row_sel)

      if (!is.null(res))
        .ssgsea_boxplot(sc, res, s, db, ct, db_dir,
                        pathways     = show_pws,
                        group_colors = group_colors,
                        top_n        = top_n,
                        contrasts    = contrasts,
                        use.padj     = use.padj)
    }
  }
  if (length(out)) do.call(rbind, out) else NULL
}

# ---------------------------------------------------------------------------
# .ssgsea_heatmap  - pathways × samples heatmap for ssGSEA scores
# ---------------------------------------------------------------------------
#' @keywords internal
.ssgsea_heatmap <- function(sc, samp_meta, db, ct, db_dir,
                             pathways        = NULL,
                             group_colors    = NULL,
                             top_n           = 30L,
                             heatmap_params  = list(),
                             heatmap_colors  = NULL,
                             run_label       = NULL,
                             contrast_label  = NULL,
                             row_selection   = NULL) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) return(invisible())
  if (!requireNamespace("circlize",       quietly = TRUE)) return(invisible())

  # Pathway selection: use pre-computed shared set when provided so the
  # heatmap and boxplot show the same rows. Fall back to top N by variance.
  if (!is.null(pathways) && length(pathways) > 0) {
    top_pws <- pathways[pathways %in% rownames(sc)]
  } else {
    row_var <- apply(sc, 1, stats::var, na.rm = TRUE)
    top_pws <- names(sort(row_var, decreasing = TRUE))[seq_len(min(top_n, nrow(sc)))]
  }
  if (length(top_pws) == 0) return(invisible())
  m <- sc[top_pws, , drop = FALSE]

  # Row z-score: normalize each pathway across samples so that pathways with
  # different absolute score magnitudes are comparable on a shared color scale.
  # Without this, high-variance pathways dominate and subtle group differences
  # in other pathways are invisible.
  m <- t(scale(t(m)))               # scale() works column-wise, so transpose twice
  m[is.nan(m) | is.infinite(m)] <- 0  # guard: zero-variance rows become NaN

  # Truncate very long pathway names so the PDF stays a manageable width
  max_rn_chars <- 80L
  rn_long <- nchar(rownames(m)) > max_rn_chars
  if (any(rn_long))
    rownames(m)[rn_long] <- paste0(substr(rownames(m)[rn_long], 1L, max_rn_chars - 3L), "...")

  # Clean pathway names
  rownames(m) <- gsub("_", " ", sub("^[A-Z]+_", "", rownames(m)))

  # Order columns by group, respecting factor levels (which carry the user's
  # preset ordering). If samp_meta$group is a plain character, this falls back
  # to alphabetical — same behavior as before.
  samp_meta <- samp_meta[colnames(m), , drop = FALSE]
  col_order  <- order(samp_meta$group)
  m          <- m[, col_order, drop = FALSE]
  samp_meta  <- samp_meta[col_order, , drop = FALSE]

  # Build clean column labels: "group | sample" — strips the full .pb key
  # (which looks like "sample||cell_type||group") so the heatmap shows
  # meaningful labels without internal || delimiters.
  col_lbls <- paste(as.character(samp_meta$group),
                    as.character(samp_meta$sample), sep = " | ")

  # Color for group annotation - group_colors is always pre-resolved by the
  # caller (.pb_ssgsea) so it is never NULL here and always covers all levels.
  grp_lvls <- levels(droplevels(samp_meta$group))
  col_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = samp_meta$group,
    col   = list(Group = group_colors[grp_lvls[grp_lvls %in% names(group_colors)]]),
    show_annotation_name = TRUE,
    annotation_name_gp   = grid::gpar(fontsize = 8)
  )

  # After z-scoring the values are centered on 0; cap the color scale at ±2
  # (most z-scores fall within this range; outliers are still visible at the
  # maximum color but the scale isn't distorted by single extreme samples).
  col_fun <- if (is.null(heatmap_colors)) {
    circlize::colorRamp2(c(-2, -1, 0, 1, 2),
                         c("#007dd1", "#b3d9f5", "white", "#f5c08a", "#ab3000"))
  } else {
    heatmap_colors
  }

  # Auto-size — extra height for 90-degree rotated column labels
  n_rows    <- nrow(m);  n_cols <- ncol(m)
  cell_h_pt <- max(7,  min(14, 400 / max(n_rows, 1)))
  cell_w_pt <- max(5,  min(12, 250 / max(n_cols, 1)))   # samples can be many
  rn_fs     <- max(5,  min(9,  cell_h_pt * 0.75))
  # row_names_max_width is set to 20 cm below, so reserve that much in the PDF.
  rn_max_in   <- max(20 / 2.54,
                     max(nchar(rownames(m)), na.rm = TRUE) * rn_fs * 0.50 / 72)
  cn_rot_h_in <- max(nchar(col_lbls), na.rm = TRUE) * 7 * 0.45 / 72
  pdf_h <- min(40, max(3.5, n_rows * cell_h_pt / 72 + cn_rot_h_in + 1.5))
  pdf_w <- min(40, max(5.0, n_cols * cell_w_pt / 72 + rn_max_in   + 2.5))

  default_args <- list(
    m,
    name              = "Row z-score",
    col               = col_fun,
    top_annotation    = col_ann,
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    show_row_names      = TRUE,
    show_column_names   = TRUE,
    column_labels       = col_lbls,
    row_names_gp        = grid::gpar(fontsize = rn_fs),
    row_names_max_width = grid::unit(20, "cm"),
    column_names_gp   = grid::gpar(fontsize = 7),
    column_names_rot  = 90,
    column_title      = if (!is.null(run_label))
                          paste0("ssGSEA scores - ", db, " ", ct, "\n", run_label)
                        else
                          paste("ssGSEA scores -", db, ct),
    column_title_gp   = grid::gpar(fontsize = 10, fontface = "bold"),
    border            = FALSE,
    use_raster        = FALSE,
    width             = grid::unit(n_cols * cell_w_pt, "pt"),
    height            = grid::unit(n_rows * cell_h_pt, "pt"),
    heatmap_legend_param = list(
      title         = "ssGSEA",
      legend_height = grid::unit(30, "mm"),
      title_gp      = grid::gpar(fontsize = 9),
      labels_gp     = grid::gpar(fontsize = 8)
    )
  )

  ht_args <- utils::modifyList(default_args, heatmap_params)
  ht      <- do.call(ComplexHeatmap::Heatmap, ht_args)

  lbl_part  <- if (!is.null(run_label))     paste0(" [", make.names(run_label),    "]") else ""
  by_part   <- if (!is.null(contrast_label)) paste0(" by ", make.names(contrast_label)) else ""
  fp <- file.path(db_dir, paste0("ssGSEA ", db, " ", make.names(ct), by_part, lbl_part, " scores heatmap.pdf"))
  grDevices::pdf(fp, width = pdf_w, height = pdf_h)
  ComplexHeatmap::draw(ht, padding = grid::unit(c(5, 5, 5, 5), "mm"))
  grDevices::dev.off()

  row_desc <- if (!is.null(row_selection)) row_selection else
    paste0("top ", nrow(m), " pathways by score variance across donors")
  if (exists(".write_legend_sidecar")) .write_legend_sidecar(fp, paste0(
    "Heatmap of single-sample GSEA (ssGSEA) enrichment scores for ", ct,
    " cells, ", db, " gene-set database. ",
    "Rows: ", row_desc, "; ",
    "columns: individual donors ordered by group. ",
    "Blue = low/negative enrichment; red = high/positive enrichment. ",
    "Scores computed using GSVA::gsva() on log-CPM pseudobulk expression."
  ))
  invisible(ht)
}

# ---------------------------------------------------------------------------
# .ssgsea_boxplot  - per-contrast boxplot of top significant pathways
# ---------------------------------------------------------------------------
#' @keywords internal
.ssgsea_boxplot <- function(sc, res, samp_meta, db, ct, db_dir,
                             pathways     = NULL,
                             group_colors = NULL,
                             top_n        = 15L,
                             contrasts    = NULL,
                             use.padj     = TRUE) {
  sig_col <- if (isTRUE(use.padj)) "padj" else "p"
  # Pathway selection: use pre-computed shared set when provided so the
  # boxplot and heatmap show the same rows.
  if (!is.null(pathways) && length(pathways) > 0) {
    top_pws <- pathways[pathways %in% rownames(sc)]
  } else {
    # Fallback: top N significant pathways by |diff|
    sig <- res[!is.na(res[[sig_col]]) & res[[sig_col]] < 0.05, , drop = FALSE]
    if (!nrow(sig)) return(invisible())
    sig     <- sig[order(-abs(sig$diff)), ]
    top_pws <- unique(sig$pathway)[seq_len(min(top_n, length(unique(sig$pathway))))]
  }

  sc_sub <- sc[top_pws[top_pws %in% rownames(sc)], , drop = FALSE]
  if (!nrow(sc_sub)) return(invisible())

  # Long format for ggplot
  sc_df          <- as.data.frame(sc_sub)
  sc_df$pathway  <- rownames(sc_df)
  long_df        <- reshape2::melt(sc_df, id.vars = "pathway",
                                   variable.name = "sample", value.name = "score")
  long_df$sample <- as.character(long_df$sample)
  long_df        <- merge(long_df,
                          samp_meta[, "group", drop = FALSE],
                          by.x = "sample", by.y = 0)
  long_df$pathway_clean <- gsub("_", " ", sub("^[A-Z]+_", "", long_df$pathway))

  # Build lookup from clean name back to original pathway ID for res lookups.
  # pathway_clean may collapse multiple originals — take first match.
  pw_orig_map <- stats::setNames(
    long_df$pathway[match(unique(long_df$pathway_clean), long_df$pathway_clean)],
    unique(long_df$pathway_clean))

  pw_names <- unique(long_df$pathway_clean)
  n_pw     <- length(pw_names)
  n_col_p  <- min(5, ceiling(sqrt(n_pw)))
  n_row_p  <- ceiling(n_pw / n_col_p)

  has_signif <- requireNamespace("ggsignif", quietly = TRUE)

  # One plot per pathway, combined with wrap_plots() for clean independent margins
  plots <- lapply(pw_names, function(pw) {
    df_pw   <- long_df[long_df$pathway_clean == pw, ]
    df_pw$group <- factor(df_pw$group, levels = levels(samp_meta$group))
    p <- ggplot2::ggplot(df_pw,
                         ggplot2::aes(x = group, y = score, fill = group)) +
      ggplot2::geom_boxplot(outlier.size = 0.8, width = 0.6, linewidth = 0.4) +
      ggplot2::geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
      theme_NourMin() +
      ggplot2::theme(
        axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
        axis.title.y    = ggplot2::element_text(size = 7),
        plot.title      = ggplot2::element_text(size = 7, face = "bold"),
        legend.position = "none",
        plot.margin     = ggplot2::margin(4, 6, 6, 4, "mm")
      ) +
      ggplot2::labs(title = pw, x = NULL, y = "ssGSEA score")
    if (!is.null(group_colors))
      p <- p + ggplot2::scale_fill_manual(values = group_colors)

    # P-value brackets from limma results - only significant pairwise contrasts.
    # Contrast names follow the "la_vs_lb" pattern from .auto_contrasts_simple().
    if (has_signif && !is.null(res)) {
      pw_orig  <- pw_orig_map[[pw]]
      pw_res   <- res[!is.na(res$pathway) & res$pathway == pw_orig &
                      !is.na(res[[sig_col]]) & res[[sig_col]] < 0.05, , drop = FALSE]
      if (nrow(pw_res) > 0) {
        # Limit to top 5 most significant to avoid visual clutter
        pw_res <- pw_res[order(pw_res[[sig_col]]), , drop = FALSE]
        pw_res <- pw_res[seq_len(min(5L, nrow(pw_res))), , drop = FALSE]

        grps_present <- levels(droplevels(factor(df_pw$group)))
        pairs   <- list()
        annots  <- character()
        for (i in seq_len(nrow(pw_res))) {
          cn   <- as.character(pw_res$contrast[i])
          # Look up the contrast expression (e.g. "groupFemale.AD - groupMale.AD")
          # from the resolved contrasts list so we extract groups from the
          # expression, not by trying to parse the user's contrast name.
          expr <- if (!is.null(contrasts) && cn %in% names(contrasts))
                    contrasts[[cn]] else cn
          toks  <- unique(regmatches(expr,
                                     gregexpr("group[A-Za-z0-9._]+", expr))[[1]])
          # Only add a bracket for simple two-group contrasts; interactions
          # (4 unique tokens) don't reduce to a single pairwise comparison.
          if (length(toks) != 2L) next
          la <- sub("^group", "", toks[1])
          lb <- sub("^group", "", toks[2])
          if (!la %in% grps_present || !lb %in% grps_present) next
          matched <- TRUE
          pairs  <- c(pairs,  list(c(la, lb)))
          sig_v  <- pw_res[[sig_col]][i]
          annots <- c(annots,
                      if (sig_v < 0.001) "***"
                      else if (sig_v < 0.01) "**"
                      else "*")
        }
        if (length(pairs) > 0)
          p <- p + ggsignif::geom_signif(
            comparisons   = pairs,
            annotations   = annots,
            step_increase = 0.10,
            tip_length    = 0.02,
            textsize      = 2.5,
            vjust         = 0.5
          )
      }
    }
    p
  })

  combined <- patchwork::wrap_plots(plots, ncol = n_col_p) +
    patchwork::plot_annotation(
      title = paste("ssGSEA scores -", db, ct, "(significant pathways)"),
      theme = ggplot2::theme(
        plot.title  = ggplot2::element_text(size = 11, face = "bold"),
        plot.margin = ggplot2::margin(4, 4, 4, 4, "mm")
      )
    )

  fp <- file.path(db_dir,
                  paste0("ssGSEA ", db, " ", make.names(ct), " boxplots.pdf"))
  ggplot2::ggsave(fp, combined,
                  width    = max(6, n_col_p * 3.0),
                  height   = max(4, n_row_p * 3.2 + 1.0),
                  limitsize = FALSE)

  sig_label_long <- if (isTRUE(use.padj)) "BH-adjusted p < 0.05" else "raw p < 0.05"
  if (exists(".write_legend_sidecar")) .write_legend_sidecar(fp, paste0(
    "Boxplots of ssGSEA enrichment scores for the top ", n_pw,
    " significantly differentially enriched pathways in ", ct,
    " cells (", db, " database; ", sig_label_long, "). ",
    "Each panel shows score distributions per group across individual donors."
  ))
  invisible(combined)
}
