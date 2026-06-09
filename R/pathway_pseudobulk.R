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

#' Pseudobulk (per-sample) GSEA with limma contrasts + per-patient ssGSEA
#'
#' @param seurat_object A Seurat object with raw counts in `assay`. Either
#'   `seurat_object` or a precomputed `pseudobulk` object must be supplied.
#' @param identity_column Character. Metadata column holding cell-type identities.
#'   One independent GSEA is run for each level (e.g. `"Assignment"`).
#' @param sample_column Character. Metadata column identifying the biological
#'   replicate (donor / patient ID). This is the statistical unit - samples
#'   are the replicates in the limma-voom model. (e.g. `"Sample"`, `"Donor.ID"`).
#' @param group_column Character. Metadata column holding the comparison group
#'   label (e.g. `"Group"` where levels are `"Dementia.Female"`,
#'   `"Dementia.Male"`, `"Reference.Female"`, …). Cells with `NA` in this
#'   column are dropped before aggregation - use `NA` to exclude groups you
#'   do not want to compare (e.g. a Reference group that is not part of any
#'   contrast).
#' @param pseudobulk Optional precomputed pseudobulk. A list with elements
#'   `counts` (genes × samples matrix) and `samples` (sample-level metadata
#'   data.frame). Supply this to skip the aggregation step on subsequent runs.
#'   Returned invisibly by the function itself.
#' @param covariates Character vector of sample-level metadata columns to
#'   include as covariates in the limma-voom model (e.g. `c("Age", "PMI")`).
#'   Numeric covariates are median-imputed for missing values; factor/character
#'   covariates receive an `"unknown"` level.
#' @param contrasts Named list of limma contrast expressions written in terms
#'   of the group-means model coefficients `group<level>` (e.g.
#'   `list(FvsM_Dem = "groupDementia.Female - groupDementia.Male")`).
#'   Positive NES in the output means enrichment in the **first** term of the
#'   contrast. When `NULL` and `group_column` levels follow the `<Disease>.<Sex>`
#'   convention (2 × 2 design), a standard set of contrasts is built
#'   automatically: Female-vs-Male within each disease, Disease-vs-control
#'   within each sex, and the Sex × Disease interaction term.
#' @param pathway_sets Named list of MSigDB databases. Same structure as in
#'   [RunGSEA()] - see that function's documentation for a full list of valid
#'   collections. Default: Hallmark (`"H"`), KEGG, Reactome, WikiPathways.
#' @param species Character. `"human"` (default) or `"mouse"`. Passed to
#'   [msigdbr::msigdbr()]. See `msigdbr::msigdbr_species()` for all options.
#' @param output_dir Character. Root directory for CSV, PDF, and RDS output.
#'   Per-contrast and per-database sub-directories are created automatically.
#' @param assay Character or `NULL`. Seurat assay containing raw counts to
#'   aggregate. `NULL` uses `SeuratObject::DefaultAssay()`. For BPCells
#'   sketch workflows always specify `"RNA"` explicitly (do not use the sketch
#'   assay - pseudobulk aggregation requires raw counts for all cells).
#' @param min.cells Integer. Minimum cells a sample must contribute to a given
#'   cell-type to be included in the pseudobulk for that cell type. Samples
#'   below this threshold are dropped to avoid noisy pseudo-replicates.
#'   Default `10`.
#' @param min.samples.per.group Integer. Minimum samples per group level
#'   required within a cell type for that cell type to be tested. Cell types
#'   with too few replicates in any group are skipped. Default `3`.
#' @param run.ssgsea Logical. Also compute per-patient pathway activity
#'   scores using single-sample GSEA (ssGSEA via `GSVA::gsva()`)?
#'   Produces one boxplot PDF per contrast per database showing group-level
#'   pathway activity. Requires the `GSVA` Bioconductor package. Default
#'   `TRUE`.
#' @param top_n_heatmap Integer. Number of top pathways (by maximum |NES|
#'   across cell types) to display in the NES heatmap per contrast × database
#'   combination. Default `30`.
#' @param save.rds Logical. Save the pseudobulk counts, sample metadata, GSEA
#'   results, and ssGSEA results as a single RDS file
#'   (`RunGSEA_pseudobulk_results.rds`) in `output_dir`? Useful for
#'   resuming or re-plotting without re-running the full pipeline.
#'   Default `TRUE`.
#' @return Invisibly returns a named list:
#' \describe{
#'   \item{`pseudobulk`}{Genes × samples aggregated count matrix.}
#'   \item{`samples`}{Sample-level metadata data.frame.}
#'   \item{`gsea`}{Long-format data.frame of all fgsea results (pathway, NES,
#'     padj, cell_type, contrast, db).}
#'   \item{`ssgsea`}{ssGSEA results list (one element per cell type), or
#'     `NULL` if `run.ssgsea = FALSE`.}
#' }
#' @param redo_plots Logical. Re-generate all visualisations from previously
#'   saved results without re-running the pipeline. Default `FALSE`. When
#'   `TRUE`:
#'   \itemize{
#'     \item `RunGSEA_pseudobulk_results.rds` is loaded from `output_dir`.
#'     \item NES heatmaps and lollipop plots are regenerated from the saved
#'       `gsea` data frame.
#'     \item ssGSEA score CSVs already on disk are read back and used to
#'       regenerate the score heatmaps and contrast boxplots.
#'     \item `seurat_object`, `identity_column`, `sample_column`, and `group_column` are
#'       not needed and can be omitted.
#'   }
#'   Useful after changing `group_colors`, `heatmap_params`, or any
#'   visualisation parameter without rerunning the hours-long pipeline.
#' @export
RunGSEA_pseudobulk <- function(seurat_object       = NULL,
                               identity_column        = NULL,
                               sample_column       = NULL,
                               group_column        = NULL,
                               pseudobulk       = NULL,
                               covariates       = NULL,
                               contrasts        = NULL,
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
                               top_n_heatmap    = 30L,
                               top_n_ssgsea     = 30L,
                               group_colors     = NULL,
                               heatmap_params   = list(row_names_side      = "left",
                                                       show_row_dend       = FALSE,
                                                       row_names_max_width = grid::unit(15, "cm")),
                               save.rds         = TRUE,
                               redo_plots       = FALSE) {

  for (pkg in c("limma", "edgeR", "fgsea", "msigdbr", "ggplot2", "Matrix"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
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
    saved     <- readRDS(rds_path)
    gsea      <- saved$gsea
    ssg       <- saved$ssgsea
    samp      <- saved$samples
    contrasts <- saved$contrasts

    # Resolve group_colors from the saved sample metadata (same logic as
    # .pb_ssgsea) so heatmap and boxplot always use identical colours.
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
          .pb_nes_heatmap(sub, cn, db, db_dir, top_n_heatmap, heatmap_params)

          for (ct in unique(sub$cell_type)) {
            fg <- sub[sub$cell_type == ct, , drop = FALSE]
            .pb_lollipop(fg, db, ct, cn, db_dir)
          }
        }
      }
    }

    # ── ssGSEA heatmaps + boxplots ───────────────────────────────────────────
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
          if (!is.null(res_ct) && nrow(res_ct) > 0 &&
              any(!is.na(res_ct$padj) & res_ct$padj < 0.05)) {
            sig_sub  <- res_ct[!is.na(res_ct$padj) & res_ct$padj < 0.05, ]
            sig_sub  <- sig_sub[order(-abs(sig_sub$diff)), ]
            show_pws <- unique(sig_sub$pathway)[
              seq_len(min(top_n_ssgsea, length(unique(sig_sub$pathway))))]
          } else {
            row_var  <- apply(sc, 1, stats::var, na.rm = TRUE)
            show_pws <- names(sort(row_var, decreasing = TRUE))[
              seq_len(min(top_n_ssgsea, nrow(sc)))]
          }

          .ssgsea_heatmap(sc, s, db, ct, db_dir,
                          pathways       = show_pws,
                          group_colors   = group_colors,
                          top_n          = top_n_ssgsea,
                          heatmap_params = heatmap_params)

          if (!is.null(res_ct) && nrow(res_ct) > 0)
            .ssgsea_boxplot(sc, res_ct, s, db, ct, db_dir,
                            pathways     = show_pws,
                            group_colors = group_colors,
                            top_n        = top_n_ssgsea)
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
      "auto-detected from group_column levels"

    .write_subdir_params(output_dir, list(
      date               = format(Sys.Date()),
      gsea_method        = "pseudobulk (limma-voom)",
      gsea_ident_col     = identity_column,
      gsea_sample_col    = sample_column,
      gsea_group_col     = group_column,
      gsea_databases     = names(pathway_sets),
      gsea_contrasts     = contrast_str,
      gsea_min_cells     = min.cells,
      gsea_min_samples   = min.samples.per.group,
      gsea_species       = species,
      gsea_run_ssgsea    = run.ssgsea,
      methods_text       = paste0(
        "Donor-level gene set enrichment analysis was performed using a ",
        "pseudobulk approach. Raw counts (assay: '", assay %||% "default",
        "') were summed per biological sample ('", sample_column, "') within ",
        "each cell-type subset ('", identity_column, "'), retaining only ",
        "sample × cell-type combinations with ≥ ", min.cells,
        " cells and groups with ≥ ", min.samples.per.group,
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
        if (run.ssgsea)
          paste0(" Per-patient pathway activity was additionally quantified ",
                 "using single-sample GSEA (ssGSEA via GSVA) on the ",
                 "pseudobulk expression matrix.")
        else ""
      )
    ))
  }

  ## ---- gene sets (query once) ----
  fgsea_dbs <- lapply(names(pathway_sets), function(db) {
    ps <- pathway_sets[[db]]
    mdf <- .msigdbr_get(species = species, category = ps$category,
                        subcategory = ps$subcategory)
    split(mdf$gene_symbol, mdf$gs_name)
  })
  names(fgsea_dbs) <- names(pathway_sets)

  ## ---- pseudobulk: reuse precomputed, or aggregate (sum counts per sample x ident) ----
  if (!is.null(pseudobulk)) {
    message("Reusing precomputed pseudobulk (skipping aggregation).")
    agg  <- as.matrix(pseudobulk$counts)
    samp <- pseudobulk$samples
    if (is.null(samp$.pb)) samp$.pb <- colnames(agg)
    if (is.null(samp$sample))    samp$sample    <- sub("\\|\\|.*", "", samp$.pb)
    if (is.null(samp$cell_type)) samp$cell_type <- sub(".*\\|\\|", "", samp$.pb)
    if (is.null(samp$group) && !is.null(group_column) && group_column %in% colnames(samp))
      samp$group <- samp[[group_column]]
    rownames(samp) <- samp$.pb; samp <- samp[colnames(agg), ]
  } else {
    if (is.null(seurat_object)) stop("Provide either `seurat_object` or a precomputed `pseudobulk`.")
    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(seurat_object)
    md <- seurat_object@meta.data
    md <- md[!is.na(md[[group_column]]) & !is.na(md[[identity_column]]) & !is.na(md[[sample_column]]), ]
    seurat_object <- subset(seurat_object, cells = rownames(md))
    seurat_object$.pb <- paste(seurat_object@meta.data[[sample_column]],
                            seurat_object@meta.data[[identity_column]], sep = "||")
    message("Aggregating pseudobulk (", assay, " counts) by ", sample_column, " x ", identity_column, " ...")
    agg <- Seurat::AggregateExpression(seurat_object, assays = assay, group.by = ".pb", slot = "counts")[[assay]]
    agg <- as.matrix(agg)
    cc <- as.data.frame(table(seurat_object$.pb)); colnames(cc) <- c(".pb", "n_cells")
    samp <- data.frame(.pb = colnames(agg), stringsAsFactors = FALSE)
    samp$sample    <- sub("\\|\\|.*", "", samp$.pb)
    samp$cell_type <- sub(".*\\|\\|", "", samp$.pb)
    sm <- md[!duplicated(md[[sample_column]]), c(sample_column, group_column, covariates)]
    colnames(sm)[1:2] <- c("sample", "group")
    samp <- merge(samp, sm, by = "sample", all.x = TRUE)
    samp <- merge(samp, cc, by = ".pb", all.x = TRUE)
    rownames(samp) <- samp$.pb; samp <- samp[colnames(agg), ]
  }
  samp$group <- droplevels(factor(samp$group))
  glevels <- levels(samp$group)

  # covariate cleaning
  for (cv in covariates) {
    x <- samp[[cv]]
    if (is.numeric(x)) { x[is.na(x)] <- stats::median(x, na.rm = TRUE); samp[[cv]] <- as.numeric(scale(x)) }
    else { x <- as.character(x); x[is.na(x)] <- "unknown"; samp[[cv]] <- factor(x) }
  }

  ## ---- contrasts ----
  if (is.null(contrasts)) contrasts <- .auto_contrasts_2x2(glevels)
  message("Contrasts: ", paste(names(contrasts), collapse = ", "))

  ## ---- per cell type: limma-voom + fgsea ----
  idents <- sort(unique(samp$cell_type))
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
    valid <- contrasts[vapply(contrasts, function(ex) all(.coef_in(ex, colnames(design))), logical(1))]
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
        utils::write.csv(fg[order(-fg$NES), ],
                         file.path(db_dir, paste0(db, " ", make.names(ct), " ", cn, ".csv")), row.names = FALSE)
        .pb_lollipop(fg, db, ct, cn, db_dir)
      }
    }
    message("  ", ct, ": ", length(valid), " contrasts x ", length(fgsea_dbs), " DBs")
  }
  gsea <- if (length(gsea_all)) do.call(rbind, gsea_all) else data.frame()

  ## ---- per (contrast, db) NES heatmaps: pathways x cell types ----
  if (nrow(gsea)) for (cn in unique(gsea$contrast)) for (db in names(fgsea_dbs)) {
    sub <- gsea[gsea$contrast == cn & gsea$db == db, ]
    if (!nrow(sub)) next
    .pb_nes_heatmap(sub, cn, db, file.path(output_dir, cn, db), top_n_heatmap, heatmap_params)
  }

  ## ---- per-patient ssGSEA activity ----
  ssg <- NULL
  if (run.ssgsea) ssg <- .pb_ssgsea(agg, samp, fgsea_dbs, contrasts, covariates,
                                    output_dir, min.cells, min.samples.per.group,
                                    group_colors   = group_colors,
                                    top_n          = top_n_ssgsea,
                                    heatmap_params = heatmap_params)

  if (save.rds) saveRDS(
    list(pseudobulk = agg, samples = samp, gsea = gsea, ssgsea = ssg, contrasts = contrasts),
    file.path(output_dir, "RunGSEA_pseudobulk_results.rds"))
  message("RunGSEA_pseudobulk complete.")
  invisible(list(pseudobulk = agg, samples = samp, gsea = gsea, ssgsea = ssg))
}

# ---- helpers ---------------------------------------------------------------

#' @keywords internal
.auto_contrasts_2x2 <- function(glevels) {
  parts <- strsplit(glevels, "\\."); ok <- all(lengths(parts) == 2)
  if (!ok) stop("Could not auto-build contrasts; supply `contrasts=` explicitly. Levels: ",
                paste(glevels, collapse = ", "))
  disease <- vapply(parts, `[`, "", 1); sex <- vapply(parts, `[`, "", 2)
  ud <- unique(disease); us <- unique(sex)
  g <- function(d, s) { lv <- glevels[disease == d & sex == s]; paste0("group", make.names(lv)) }
  ct <- list()
  # female vs male within each disease (positive NES = higher in females)
  fem <- grep("F", us, value = TRUE, ignore.case = TRUE)[1]; mal <- setdiff(us, fem)[1]
  for (d in ud) ct[[paste0("FvsM_", d)]] <- paste(g(d, fem), "-", g(d, mal))
  # identify control (No/Healthy/Control) vs case so disease contrasts are case - control
  # (positive NES = higher in DISEASE)
  ctrl <- ud[grepl("^(no|healthy|control|ctrl|hc)", ud, ignore.case = TRUE)]
  ctrl <- if (length(ctrl) >= 1) ctrl[1] else ud[1]
  case <- setdiff(ud, ctrl)
  for (cs in case) for (s in us)
    ct[[paste0(cs, "_vs_", ctrl, "_", s)]] <- paste(g(cs, s), "-", g(ctrl, s))
  # interaction (needs 2x2): (case.F - case.M) - (control.F - control.M)
  # positive NES = female-direction strengthens in disease
  if (length(case) >= 1 && length(us) >= 2) {
    cs <- case[1]
    ct[["SexByDisease"]] <- paste0("(", g(cs, fem), " - ", g(cs, mal), ") - (",
                                   g(ctrl, fem), " - ", g(ctrl, mal), ")")
  }
  ct
}

#' @keywords internal
.coef_in <- function(expr, coefs) {
  toks <- unique(regmatches(expr, gregexpr("group[A-Za-z0-9._]+", expr))[[1]])
  all(toks %in% coefs)
}

#' @keywords internal
.pb_lollipop <- function(fg, db, ct, cn, db_dir, n = 15) {
  fg <- fg[fg$padj < 0.25, , drop = FALSE]; if (!nrow(fg)) return(invisible())
  fg <- fg[order(-abs(fg$NES)), ][seq_len(min(n, nrow(fg))), ]
  fg$pw <- gsub("_", " ", sub("^[A-Z]+_", "", fg$pathway))
  p <- ggplot2::ggplot(fg, ggplot2::aes(stats::reorder(pw, NES), NES, color = NES > 0)) +
    ggplot2::geom_segment(ggplot2::aes(xend = pw, y = 0, yend = NES)) +
    ggplot2::geom_point(size = 3) + ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#B5232E", `FALSE` = "#2C5F8A"), guide = "none") +
    ggplot2::labs(title = paste0(db, " ", ct, " : ", cn), x = NULL, y = "NES") +
    ggplot2::theme_minimal(base_size = 10)
  fp <- file.path(db_dir, paste0(db, " - ", make.names(ct), " ", cn, " lollipop.pdf"))
  ggplot2::ggsave(fp, p, width = 6, height = 4.5)
  if (exists(".write_legend_sidecar")) .write_legend_sidecar(fp, paste0(
    "Lollipop of top ", nrow(fg), " gene sets (", db, ") for ", ct, ", contrast ", cn,
    " (donor-level pseudobulk, limma-voom moderated t pre-ranking, fgsea). NES>0 = enriched in the first term of the contrast."))
}

#' @keywords internal
.pb_nes_heatmap <- function(sub, cn, db, db_dir, top_n, heatmap_params = list()) {
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

  fp <- file.path(db_dir, paste0("GSEA ", db, " ", cn, " (NES)-Heatmap.pdf"))

  # Use .gsea_ht() for consistent sizing and colour, with stars via cell_fun
  .gsea_ht(
    mat            = m,
    title          = paste(db, cn, "(NES, pseudobulk)"),
    filepath       = fp,
    heatmap_params = utils::modifyList(
      list(cell_fun = cell_fn),   # add stars by default
      heatmap_params              # user can override cell_fun or anything else
    )
  )

  if (exists(".write_legend_sidecar")) .write_legend_sidecar(fp, paste0(
    "Heatmap of pseudobulk GSEA Normalised Enrichment Scores (", db,
    ") for contrast ", cn, " across cell types. ",
    "Rows = top ", nrow(m), " gene sets by max |NES| across cell types. ",
    "Method: donor-level pseudobulk, edgeR TMM normalization, ",
    "limma-voom moderated t-statistic pre-ranking, fgsea. ",
    "NES > 0 = enriched in the first term of the contrast. ",
    "Overlaid significance stars: * BH < 0.05, ** < 0.01, *** < 0.001."))
}

#' @keywords internal
.pb_ssgsea <- function(agg, samp, fgsea_dbs, contrasts, covariates, output_dir,
                       min.cells, min.smp,
                       group_colors   = NULL,
                       top_n          = 30L,
                       heatmap_params = list()) {

  # Resolve group_colors once so the heatmap and boxplot always use identical
  # colours regardless of which helper is called first or whether the user
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
      utils::write.csv(
        data.frame(pathway = rownames(sc), sc, check.names = FALSE),
        file.path(db_dir, paste0("ssGSEA ", db, " ", make.names(ct), " scores.csv")),
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
      valid  <- contrasts[vapply(contrasts,
                                 function(ex) .coef_in(ex, colnames(design)),
                                 logical(1))]

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
      # Significant pathways (any contrast, BH < 0.05) take priority so that
      # the heatmap and boxplot show exactly the same rows.
      # If nothing is significant, fall back to top N by variance so the
      # heatmap is still informative even without reaching significance.
      if (!is.null(res) && any(!is.na(res$padj) & res$padj < 0.05)) {
        sig_sub  <- res[!is.na(res$padj) & res$padj < 0.05, ]
        sig_sub  <- sig_sub[order(-abs(sig_sub$diff)), ]
        show_pws <- unique(sig_sub$pathway)[
          seq_len(min(top_n, length(unique(sig_sub$pathway))))]
      } else {
        row_var  <- apply(sc, 1, stats::var, na.rm = TRUE)
        show_pws <- names(sort(row_var, decreasing = TRUE))[
          seq_len(min(top_n, nrow(sc)))]
        if (!is.null(res))
          message("    No significant ssGSEA pathways for ", ct, "/", db,
                  " - heatmap shows top by variance.")
      }

      # ── Heatmap and boxplot - same pathway set ───────────────────────────
      .ssgsea_heatmap(sc, s, db, ct, db_dir,
                      pathways       = show_pws,
                      group_colors   = group_colors,
                      top_n          = top_n,
                      heatmap_params = heatmap_params)

      if (!is.null(res))
        .ssgsea_boxplot(sc, res, s, db, ct, db_dir,
                        pathways     = show_pws,
                        group_colors = group_colors,
                        top_n        = top_n)
    }
  }
  if (length(out)) do.call(rbind, out) else NULL
}

# ---------------------------------------------------------------------------
# .ssgsea_heatmap  - pathways × samples heatmap for ssGSEA scores
# ---------------------------------------------------------------------------
#' @keywords internal
.ssgsea_heatmap <- function(sc, samp_meta, db, ct, db_dir,
                             pathways       = NULL,
                             group_colors   = NULL,
                             top_n          = 30L,
                             heatmap_params = list()) {
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
  # different absolute score magnitudes are comparable on a shared colour scale.
  # Without this, high-variance pathways dominate and subtle group differences
  # in other pathways are invisible.
  m <- t(scale(t(m)))               # scale() works column-wise, so transpose twice
  m[is.nan(m) | is.infinite(m)] <- 0  # guard: zero-variance rows become NaN

  # Clean pathway names
  rownames(m) <- gsub("_", " ", sub("^[A-Z]+_", "", rownames(m)))

  # Order columns by group
  samp_meta <- samp_meta[colnames(m), , drop = FALSE]
  col_order  <- order(as.character(samp_meta$group))
  m          <- m[, col_order, drop = FALSE]
  samp_meta  <- samp_meta[col_order, , drop = FALSE]

  # Colour for group annotation - group_colors is always pre-resolved by the
  # caller (.pb_ssgsea) so it is never NULL here and always covers all levels.
  grp_lvls <- levels(droplevels(samp_meta$group))
  col_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = samp_meta$group,
    col   = list(Group = group_colors[grp_lvls[grp_lvls %in% names(group_colors)]]),
    show_annotation_name = TRUE,
    annotation_name_gp   = grid::gpar(fontsize = 8)
  )

  # After z-scoring the values are centred on 0; cap the colour scale at ±2
  # (most z-scores fall within this range; outliers are still visible at the
  # maximum colour but the scale isn't distorted by single extreme samples).
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))

  # Auto-size
  n_rows    <- nrow(m);  n_cols <- ncol(m)
  cell_h_pt <- max(7,  min(14, 400 / max(n_rows, 1)))
  cell_w_pt <- max(5,  min(12, 250 / max(n_cols, 1)))   # samples can be many
  rn_fs     <- max(5,  min(9,  cell_h_pt * 0.75))
  rn_max_in <- max(nchar(rownames(m)), na.rm = TRUE) * rn_fs * 0.50 / 72
  pdf_h <- min(40, max(3.5, n_rows * cell_h_pt / 72 + 2.0))
  pdf_w <- min(40, max(5.0, n_cols * cell_w_pt / 72 + rn_max_in + 2.5))

  default_args <- list(
    m,
    name              = "Row z-score",
    col               = col_fun,
    top_annotation    = col_ann,
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = grid::gpar(fontsize = rn_fs),
    column_names_gp   = grid::gpar(fontsize = 7),
    column_names_rot  = 90,
    column_title      = paste("ssGSEA scores -", db, ct),
    column_title_gp   = grid::gpar(fontsize = 10, fontface = "bold"),
    border            = FALSE,
    use_raster        = n_rows > 100 || n_cols > 50,
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

  fp <- file.path(db_dir, paste0("ssGSEA ", db, " ", make.names(ct), " scores heatmap.pdf"))
  grDevices::pdf(fp, width = pdf_w, height = pdf_h)
  ComplexHeatmap::draw(ht, padding = grid::unit(c(5, 5, 5, 5), "mm"))
  grDevices::dev.off()

  if (exists(".write_legend_sidecar")) .write_legend_sidecar(fp, paste0(
    "Heatmap of single-sample GSEA (ssGSEA) enrichment scores for ", ct,
    " cells, ", db, " gene-set database. ",
    "Rows: top ", nrow(m), " pathways by score variance across donors; ",
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
                             top_n        = 15L) {
  # Pathway selection: use pre-computed shared set when provided so the
  # boxplot and heatmap show the same rows.
  if (!is.null(pathways) && length(pathways) > 0) {
    top_pws <- pathways[pathways %in% rownames(sc)]
  } else {
    # Fallback: top N significant pathways by |diff|
    sig <- res[!is.na(res$padj) & res$padj < 0.05, , drop = FALSE]
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

  n_pw     <- length(unique(long_df$pathway_clean))
  n_col_p  <- min(5, ceiling(sqrt(n_pw)))
  n_row_p  <- ceiling(n_pw / n_col_p)

  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = group, y = score, fill = group)) +
    ggplot2::geom_boxplot(outlier.size = 0.8, width = 0.6, linewidth = 0.4) +
    ggplot2::geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
    ggplot2::facet_wrap(~ pathway_clean, scales = "free_y", ncol = n_col_p) +
    theme_NourMin() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      strip.text   = ggplot2::element_text(size  = 7),
      legend.position = "none"
    ) +
    ggplot2::labs(title = paste("ssGSEA scores -", db, ct,
                                "(significant pathways)"),
                  x = NULL, y = "ssGSEA score")

  if (!is.null(group_colors))
    p <- p + ggplot2::scale_fill_manual(values = group_colors)

  fp <- file.path(db_dir,
                  paste0("ssGSEA ", db, " ", make.names(ct), " boxplots.pdf"))
  ggplot2::ggsave(fp, p,
                  width    = max(6, n_col_p * 3.0),
                  height   = max(4, n_row_p * 2.8 + 1.5),
                  limitsize = FALSE)

  if (exists(".write_legend_sidecar")) .write_legend_sidecar(fp, paste0(
    "Boxplots of ssGSEA enrichment scores for the top ", n_pw,
    " significantly differentially enriched pathways in ", ct,
    " cells (", db, " database; BH-adjusted p < 0.05). ",
    "Each panel shows score distributions per group across individual donors."
  ))
  invisible(p)
}
