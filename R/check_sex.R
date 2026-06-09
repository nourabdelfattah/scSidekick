# =============================================================================
# scSidekick — Sex verification
#
# CheckSex() verifies the sex annotation of samples/donors by comparing
# annotated sex against expression of sex-specific genes.
#
# Decision level: PSEUDOBULK (per sample), not per cell.
#
# Why pseudobulk?
#   Single-cell dropout is random. A female cell may fail to capture XIST
#   in any given sequencing run — that is noise, not a biology. But a truly
#   male sample will have ZERO XIST across all cells. Averaging across
#   hundreds to thousands of cells makes the sex signal robust: a female
#   sample with 60% XIST dropout still has a bulk mean ~0.5; a mislabeled
#   male sample has bulk mean = 0.
#
# Mismatch rule (both conditions required):
#   Expected sex gene completely absent at pseudobulk level
#   AND opposite-sex gene clearly present at pseudobulk level.
#   Neither condition alone triggers a flag.
#
# Two female-scoring modes (auto-selected):
#   PRIMARY  — XIST/Tsix expression (X-inactivation signal)
#   FALLBACK — X-paralog biallelic dosage (DDX3X, KDM6A, KDM5C)
#              used automatically when XIST/Tsix are absent from the panel
#              (common when lncRNAs are filtered out or panel is targeted)
#
# Per-cell scores (sex_score_f, sex_score_m, sex_ratio) are retained for
# UMAP visualisation. predicted_sex and sex_mismatch are sample-level values
# propagated identically to all cells within each sample.
# =============================================================================

# ── Gene registries ───────────────────────────────────────────────────────────

.nk_sex_genes <- list(
  human = list(
    female_primary = c("XIST", "TSIX"),
    female_xpar    = c("DDX3X", "KDM6A", "KDM5C"),
    male           = c("DDX3Y", "UTY", "KDM5D", "RPS4Y1")
  ),
  mouse = list(
    female_primary = c("Xist", "Tsix"),
    female_xpar    = c("Ddx3x", "Kdm6a", "Kdm5c"),
    male           = c("Ddx3y", "Uty", "Kdm5d", "Eif2s3y")
  ),
  zebrafish = list(       # sex determination differs; no universal defaults
    female_primary = character(0),
    female_xpar    = character(0),
    male           = character(0)
  )
)

# ── Internal: plots ───────────────────────────────────────────────────────────

.nk_sex_plots <- function(seurat_object, sex_col, sample_col,
                           all_sex_genes, output_dir) {

  genes_in_panel <- intersect(all_sex_genes, rownames(seurat_object))
  has_redux      <- length(seurat_object@reductions) > 0
  has_sex_col    <- !is.null(sex_col) && sex_col %in% colnames(seurat_object@meta.data)
  has_sample_col <- !is.null(sample_col) && sample_col %in% colnames(seurat_object@meta.data)

  # 1. DimPlot: predicted sex — columns split by annotated sex
  #    Coloring = predicted sex, facet columns = what the sample was annotated as.
  #    Mismatches are immediately visible as wrong-colored cells in a panel.
  if (has_redux && "predicted_sex" %in% colnames(seurat_object@meta.data)) {
    tryCatch(
      PlotDimPlots(
        seurat_object,
        group.by   = "predicted_sex",
        split.by   = if (has_sex_col) sex_col else NULL,
        colors     = c(Female       = "#F37388",
                       Male         = "#003F5C",
                       Ambiguous    = "#F9BF31",
                       Undetermined = "gray80"),
        output_dir = output_dir,
        file_name  = "sex_check_predicted",
        pt.size    = 0.3,
        legendnrow = 1
      ),
      error = function(e)
        warning("scSidekick: sex DimPlot failed: ", e$message)
    )
  }

  # 2. FeaturePlots: one PDF per sex gene, columns split by annotated sex
  #    Immediately shows whether each gene is expressed in the expected panel
  #    (e.g. XIST visible only in Female panel, DDX3Y only in Male panel).
  if (has_redux && length(genes_in_panel) > 0) {
    tryCatch(
      GenerateFeatureMaps(
        seurat_object,
        features    = genes_in_panel,
        split.by    = if (has_sex_col) sex_col else NULL,
        output_dir  = output_dir,
        object_name = "sex_check",
        pt.size     = 0.3
      ),
      error = function(e)
        warning("scSidekick: sex FeaturePlots failed: ", e$message)
    )
  }

  # 3. DotPlot: sex genes × all samples, faceted by annotated sex
  #    min.pct.exp = 0 is essential: sex genes absent in a sample (e.g. XIST in
  #    a male or mislabeled sample) must appear as empty dots — that zero IS the
  #    signal we are looking for. Filtering at min.pct.exp > 0 would remove them.
  if (has_sample_col && length(genes_in_panel) > 0) {
    tryCatch(
      FastDotPlot(
        seurat_object,
        features      = genes_in_panel,
        group.by      = sample_col,
        split.by      = if (has_sex_col) sex_col else NULL,
        ClusterGenes  = FALSE,
        ClusterGroups = FALSE,
        min.pct.exp   = 0,
        output_dir    = output_dir,
        file_name     = "sex_check_dotplot"
      ),
      error = function(e)
        warning("scSidekick: sex DotPlot failed: ", e$message)
    )
  }

  invisible(NULL)
}

# ── CheckSex ──────────────────────────────────────────────────────────────────

#' Verify sample sex annotation using sex-specific gene expression
#'
#' Compares each sample's annotated sex against pseudobulk expression of
#' sex-specific genes. Sex is called at the **sample level** (averaged across
#' all cells in the donor), not per cell. This avoids false positives from
#' single-cell dropout: a female cell may fail to capture XIST in any given
#' sequencing run, but a truly male sample will have zero XIST across all cells.
#'
#' A sample is flagged only when **both** of the following hold:
#' 1. The expected sex marker is **completely absent** at pseudobulk level
#'    (bulk mean ≤ `noise_floor`).
#' 2. The opposite-sex marker is **clearly present** at pseudobulk level
#'    (bulk mean > `noise_floor`).
#'
#' This strict two-condition rule avoids flags on:
#' - Dropout (expected marker absent but no opposite signal either).
#' - X–Y paralog multi-mapping artifacts: DDX3X/DDX3Y, KDM6A/UTY, and
#'   KDM5C/KDM5D share 70–85\% sequence identity, so single cells can show
#'   spurious Y-gene signal from misaligned X reads. At the pseudobulk level
#'   this artifact is tiny (scattered reads in a few cells) and stays well
#'   below `noise_floor`.
#'
#' **Female-scoring modes** (auto-selected):
#' - **Primary** (XIST-based): uses `XIST`/`Xist`, the X-inactivation signal.
#'   Most specific when available.
#' - **Fallback** (X-paralog dosage): uses `DDX3X`, `KDM6A`, `KDM5C` (or
#'   mouse equivalents). Triggered when XIST/Tsix are absent from the panel
#'   (common when lncRNAs are filtered or targeted panels are used). In this
#'   mode mismatch is called using sample-level prediction rather than the
#'   noise_floor threshold, because X-paralogs are expressed in males at lower
#'   (not zero) levels.
#'
#' When `cell_type_col` is supplied, a per-cell-type pseudobulk breakdown is
#' also computed. A sample flagged only in specific cell types (e.g., immune
#' cells) suggests BMT/chimerism; one flagged across all cell types suggests a
#' labeling error.
#'
#' @param seurat_object A Seurat object.
#' @param sex_col Character. Metadata column containing annotated sex.
#'   Accepted values: `"F"/"M"`, `"Female"/"Male"`, `"female"/"male"`,
#'   `"1"/"2"`. Default `"sex"`.
#' @param sample_col Character. Metadata column identifying donors/samples.
#'   Used for per-donor pseudobulk aggregation. `NULL` skips sample-level
#'   analysis and returns only per-cell scores for visualisation.
#'   Default `"Sample"`.
#' @param cell_type_col Character or `NULL`. When supplied, a per-cell-type
#'   pseudobulk table is computed in addition to the sample-level table —
#'   useful for distinguishing BMT/chimerism (cell-type-specific mismatch)
#'   from labeling errors (all cell types mismatched). Default `NULL`.
#' @param species Character. `"mouse"`, `"human"`, or `"zebrafish"`. Auto-
#'   detected from [AnnotateFeatures()] if previously run. Zebrafish has no
#'   universal sex-gene defaults — supply `female_genes` and `male_genes`.
#' @param female_genes Character vector or `NULL`. Override default female
#'   marker genes. When `NULL` the package defaults are used (XIST/Tsix with
#'   automatic X-paralog fallback).
#' @param male_genes Character vector or `NULL`. Override default male marker
#'   genes. When `NULL` the package defaults are used.
#' @param assay Character or `NULL`. Assay to pull expression from. Defaults
#'   to `DefaultAssay()`.
#' @param layer Character. Layer within the assay. Default `"data"`
#'   (log-normalized counts). Falls back to `"counts"` automatically if
#'   `"data"` is not present (BPCells objects). Only ~4–8 sex genes are
#'   materialised into memory — near-instantaneous even on million-cell objects.
#' @param noise_floor Numeric. Pseudobulk mean expression threshold that
#'   separates "completely absent" from "present" at the sample level.
#'   Default `0.05`. Because this is applied to the sample mean (averaged over
#'   hundreds–thousands of cells), alignment-artifact signal from X–Y paralog
#'   multi-mapping (typically 0.001–0.01 in pseudobulk space) is well below
#'   this threshold, while true sex-gene expression (0.2–2.0) is well above.
#'   Raise to `0.1` if you see "Undetermined" on low-depth samples; lower to
#'   `0.02` for very high-depth data.
#' @param output_dir Character or `NULL`. Directory for saved figures and CSV
#'   summary tables. Walks up from [PrepObject()] when `NULL`.
#' @param plot Logical. Generate and save figures. Default `TRUE`.
#'
#' @return The input Seurat object with the following metadata columns added:
#'   \describe{
#'     \item{`predicted_sex`}{Sample-level sex call (`"Female"` / `"Male"` /
#'       `"Ambiguous"` / `"Undetermined"`), propagated identically to all
#'       cells within each sample.}
#'     \item{`sex_score_f`}{Per-cell mean female-marker expression (for
#'       visualisation; not used for the mismatch decision).}
#'     \item{`sex_score_m`}{Per-cell mean male-marker expression.}
#'     \item{`sex_ratio`}{Per-cell `log2(female_score + 0.01) / (male_score +
#'       0.01)`. Positive = female signal, negative = male signal. Used for
#'       the continuous UMAP ratio plot.}
#'     \item{`sex_mismatch`}{Sample-level mismatch flag (`TRUE`/`FALSE`/`NA`),
#'       propagated to all cells in the sample. `TRUE` only when expected
#'       markers are completely absent **and** opposite markers are present.}
#'   }
#' @seealso [AnnotateFeatures()]
#' @export
CheckSex <- function(seurat_object,
                     sex_col       = "sex",
                     sample_col    = "Sample",
                     cell_type_col = NULL,
                     species       = NULL,
                     female_genes  = NULL,
                     male_genes    = NULL,
                     assay         = NULL,
                     layer         = "data",
                     noise_floor   = 0.05,
                     output_dir    = NULL,
                     plot          = TRUE) {

  # ── 0. Resolve assay and layer ─────────────────────────────────────────────
  assay <- assay %||% Seurat::DefaultAssay(seurat_object)

  available_layers <- tryCatch(
    Seurat::Layers(seurat_object[[assay]]),
    error = function(e) character(0)
  )
  if (length(available_layers) > 0 && !layer %in% available_layers) {
    fallback <- if ("counts" %in% available_layers) "counts" else available_layers[1]
    message("scSidekick: Layer '", layer, "' not found in assay '", assay,
            "' — using '", fallback, "' instead.")
    layer <- fallback
  }

  # ── 1. Resolve species ──────────────────────────────────────────────────────
  species <- species %||% seurat_object@misc$nk_annotation_species
  if (is.null(species))
    stop("'species' must be supplied (\"mouse\", \"human\", or \"zebrafish\"), ",
         "or stored via AnnotateFeatures().")
  species  <- match.arg(species, names(.nk_sex_genes))
  defaults <- .nk_sex_genes[[species]]
  panel    <- rownames(seurat_object)

  # ── 2. Resolve male genes ───────────────────────────────────────────────────
  male_g <- intersect(male_genes %||% defaults$male, panel)
  if (length(male_g) == 0) {
    if (length(defaults$male) == 0)
      stop("No default male sex genes defined for ", species,
           ".\n  Supply male_genes manually.")
    stop("None of the male sex genes (",
         paste(defaults$male, collapse = ", "),
         ") are present in this panel.\n  Supply male_genes manually.")
  }

  # ── 3. Resolve female genes (primary → X-paralog → error) ──────────────────
  female_mode <- "user"
  if (!is.null(female_genes)) {
    female_g <- intersect(female_genes, panel)
    if (length(female_g) == 0)
      stop("Supplied female_genes not found in panel.")
  } else {
    primary_avail <- intersect(defaults$female_primary, panel)
    if (length(primary_avail) > 0) {
      female_g    <- primary_avail
      female_mode <- "primary"
    } else {
      xpar_avail <- intersect(defaults$female_xpar, panel)
      if (length(xpar_avail) > 0) {
        female_g    <- xpar_avail
        female_mode <- "xpar"
        message(
          "scSidekick: ",
          paste(defaults$female_primary, collapse = "/"),
          " not found in panel — switching to X-paralog female scoring\n",
          "  using: ", paste(female_g, collapse = ", "), "\n",
          "  Female scores reflect biallelic X-paralog expression,",
          " not X-inactivation signal."
        )
      } else {
        if (length(defaults$female_primary) == 0)
          stop("No default female sex genes defined for ", species,
               ".\n  Supply female_genes manually.")
        stop(
          "No female sex genes found in this panel.\n",
          "  Primary (", paste(defaults$female_primary, collapse = ", "), ")",
          " and X-paralog fallbacks (",
          paste(defaults$female_xpar, collapse = ", "),
          ") are all absent.\n  Supply female_genes manually."
        )
      }
    }
  }

  all_sex_genes <- unique(c(female_g, male_g))
  message(
    "scSidekick: CheckSex — ",
    switch(female_mode,
      primary = "XIST-based female scoring",
      xpar    = "X-paralog female scoring (XIST absent)",
      user    = "user-supplied female genes"
    ), "\n",
    "  Female : ", paste(female_g, collapse = ", "), "\n",
    "  Male   : ", paste(male_g,   collapse = ", ")
  )

  # ── 4. Fetch expression (only ~4–8 genes materialised — BPCells-safe) ───────
  expr <- Seurat::FetchData(seurat_object,
                             vars  = all_sex_genes,
                             assay = assay,
                             layer = layer)

  f_cols <- intersect(female_g, colnames(expr))
  m_cols <- intersect(male_g,   colnames(expr))

  # Per-cell scores — kept for visualisation only, not used for mismatch calls
  female_score <- if (length(f_cols) > 0)
    rowMeans(expr[, f_cols, drop = FALSE]) else rep(0, nrow(expr))
  male_score   <- if (length(m_cols) > 0)
    rowMeans(expr[, m_cols, drop = FALSE]) else rep(0, nrow(expr))
  sex_ratio    <- log2((female_score + 0.01) / (male_score + 0.01))

  # ── 5. Sample-level pseudobulk sex prediction ──────────────────────────────
  meta <- seurat_object@meta.data

  # Normalize annotated sex labels
  annotated_norm <- if (!is.null(sex_col) && sex_col %in% colnames(meta)) {
    ann_raw <- as.character(meta[[sex_col]])
    dplyr::case_when(
      tolower(ann_raw) %in% c("f", "female", "2") ~ "Female",
      tolower(ann_raw) %in% c("m", "male",   "1") ~ "Male",
      TRUE ~ ann_raw
    )
  } else {
    rep(NA_character_, nrow(meta))
  }

  if (!is.null(sample_col) && sample_col %in% colnames(meta)) {

    sample_labels <- meta[[sample_col]]
    samples       <- unique(sample_labels)

    # ── per-sample pseudobulk ────────────────────────────────────────────────
    sample_tbl <- do.call(rbind, lapply(samples, function(s) {

      idx    <- which(sample_labels == s)
      bulk_f <- mean(female_score[idx])
      bulk_m <- mean(male_score[idx])

      pred_sex <- dplyr::case_when(
        bulk_f >  noise_floor & bulk_m <= noise_floor ~ "Female",
        bulk_m >  noise_floor & bulk_f <= noise_floor ~ "Male",
        bulk_f >  noise_floor & bulk_m >  noise_floor ~ "Ambiguous",
        TRUE ~ "Undetermined"
      )

      ann <- {
        tbl <- sort(table(annotated_norm[idx]), decreasing = TRUE)
        names(tbl)[1]
      }

      # Mismatch: BOTH conditions required
      #   – Expected marker completely absent (bulk ≤ noise_floor)
      #   – Opposite marker clearly present  (bulk  > noise_floor)
      # In X-paralog mode the male signal is not binary so we compare
      # predicted vs annotated sex instead of using a hard threshold.
      mismatch <- if (female_mode %in% c("primary", "user")) {
        (!is.na(ann) & ann == "Female" & bulk_f <= noise_floor & bulk_m > noise_floor) |
        (!is.na(ann) & ann == "Male"   & bulk_m <= noise_floor & bulk_f > noise_floor)
      } else {
        # X-paralog mode: use predicted sex comparison (male signal not zero in males)
        !is.na(ann) & !is.na(pred_sex) &
          pred_sex %in% c("Female", "Male") & ann != pred_sex
      }

      flag <- dplyr::case_when(
        mismatch                  ~ paste0("⚠  Annotated ", ann, " — predicted ", pred_sex),
        pred_sex == "Ambiguous"   ~ "ℹ  Both F and M markers detected (possible BMT/chimerism)",
        pred_sex == "Undetermined"~ "ℹ  Neither marker detected (low coverage?)",
        TRUE                      ~ "✓"
      )

      data.frame(
        Sample    = s,
        N_Cells   = length(idx),
        Annotated = ann,
        Predicted = pred_sex,
        Bulk_F    = round(bulk_f, 4),
        Bulk_M    = round(bulk_m, 4),
        Mismatch  = mismatch,
        Flag      = flag,
        stringsAsFactors = FALSE
      )
    }))

    # ── optional cell-type pseudobulk (BMT detection) ────────────────────────
    ct_tbl <- NULL
    if (!is.null(cell_type_col) && cell_type_col %in% colnames(meta)) {

      ct_tbl <- do.call(rbind, lapply(samples, function(s) {
        idx_s <- which(sample_labels == s)
        ann_s <- sample_tbl$Annotated[sample_tbl$Sample == s]
        cts   <- unique(meta[[cell_type_col]][idx_s])

        do.call(rbind, lapply(cts, function(ct) {
          idx    <- idx_s[meta[[cell_type_col]][idx_s] == ct]
          bulk_f <- mean(female_score[idx])
          bulk_m <- mean(male_score[idx])

          ct_mismatch <- if (female_mode %in% c("primary", "user")) {
            (!is.na(ann_s) & ann_s == "Female" & bulk_f <= noise_floor & bulk_m > noise_floor) |
            (!is.na(ann_s) & ann_s == "Male"   & bulk_m <= noise_floor & bulk_f > noise_floor)
          } else {
            ct_pred <- dplyr::case_when(
              bulk_f > noise_floor & bulk_m <= noise_floor ~ "Female",
              bulk_m > noise_floor & bulk_f <= noise_floor ~ "Male",
              bulk_f > noise_floor & bulk_m >  noise_floor ~ "Ambiguous",
              TRUE ~ "Undetermined"
            )
            !is.na(ann_s) & ct_pred %in% c("Female","Male") & ann_s != ct_pred
          }

          data.frame(
            Sample    = s,
            Cell_Type = ct,
            N_Cells   = length(idx),
            Annotated = ann_s,
            Bulk_F    = round(bulk_f, 4),
            Bulk_M    = round(bulk_m, 4),
            Mismatch  = ct_mismatch,
            stringsAsFactors = FALSE
          )
        }))
      }))
    }

    # ── propagate sample-level values to per-cell metadata ───────────────────
    cell_to_sample  <- match(sample_labels, sample_tbl$Sample)
    predicted_sex   <- sample_tbl$Predicted[cell_to_sample]
    sex_mismatch    <- sample_tbl$Mismatch[cell_to_sample]

  } else {
    # No sample_col: per-cell prediction only (for visualisation)
    predicted_sex <- dplyr::case_when(
      female_score >  noise_floor & male_score <= noise_floor ~ "Female",
      male_score   >  noise_floor & female_score <= noise_floor ~ "Male",
      female_score >  noise_floor & male_score   >  noise_floor ~ "Ambiguous",
      TRUE ~ "Undetermined"
    )
    sex_mismatch <- rep(NA, nrow(meta))
    sample_tbl   <- NULL
    ct_tbl       <- NULL
  }

  # ── 6. Add to Seurat metadata ───────────────────────────────────────────────
  seurat_object <- Seurat::AddMetaData(seurat_object,
    metadata = data.frame(
      predicted_sex = predicted_sex,
      sex_score_f   = round(female_score, 4),
      sex_score_m   = round(male_score,   4),
      sex_ratio     = round(sex_ratio,    4),
      sex_mismatch  = sex_mismatch,
      row.names     = colnames(seurat_object),
      stringsAsFactors = FALSE
    )
  )

  # ── 7. Print summary ─────────────────────────────────────────────────────────
  if (!is.null(sample_tbl)) {

    n_mismatch   <- sum(sample_tbl$Mismatch,   na.rm = TRUE)
    n_ambiguous  <- sum(sample_tbl$Predicted == "Ambiguous",    na.rm = TRUE)
    n_undeter    <- sum(sample_tbl$Predicted == "Undetermined", na.rm = TRUE)

    message("\n── Sex check summary ", paste(rep("─", 50), collapse = ""))
    message("  ", nrow(sample_tbl), " samples  |  ",
            "Female marker: ", paste(female_g, collapse=", "), "  |  ",
            "Male marker: ",   paste(male_g,   collapse=", "))

    if (n_mismatch > 0) {
      message("\n⚠  ", n_mismatch, " sample(s) with sex mismatch ",
              "(expected marker absent, opposite marker present):")
      mis <- sample_tbl[sample_tbl$Mismatch, ]
      print(mis[, c("Sample","N_Cells","Annotated","Predicted","Bulk_F","Bulk_M","Flag")],
            row.names = FALSE)
    } else {
      message("✓  All samples are consistent with their annotated sex.")
    }

    if (n_ambiguous > 0) {
      message("\nℹ  ", n_ambiguous,
              " sample(s) with both F and M markers detected (Ambiguous —",
              " possible BMT or chimerism; supply cell_type_col to investigate):")
      amb <- sample_tbl[sample_tbl$Predicted == "Ambiguous", ]
      print(amb[, c("Sample","N_Cells","Annotated","Predicted","Bulk_F","Bulk_M","Flag")],
            row.names = FALSE)
    }

    if (n_undeter > 0) {
      message("\nℹ  ", n_undeter,
              " sample(s) where neither sex marker was detected",
              " (possible low-coverage samples).")
    }

    message(paste(rep("─", 68), collapse = ""), "\n")

    # Cell-type breakdown when BMT suspected
    if (!is.null(ct_tbl) && any(ct_tbl$Mismatch, na.rm = TRUE)) {
      message("Cell-type pseudobulk breakdown (samples with any mismatch):")
      affected_samples <- unique(ct_tbl$Sample[ct_tbl$Mismatch])
      ct_show <- ct_tbl[ct_tbl$Sample %in% affected_samples, ]
      print(ct_show, row.names = FALSE)
      message()
    }

    # Save
    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(sample_tbl,
        file.path(output_dir, "sex_check_summary.csv"), row.names = FALSE)
      if (!is.null(ct_tbl))
        utils::write.csv(ct_tbl,
          file.path(output_dir, "sex_check_celltype.csv"), row.names = FALSE)
    }
  }

  # ── 8. Figures ──────────────────────────────────────────────────────────────
  if (plot) {
    output_dir <- output_dir %||%
      if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL
    .nk_sex_plots(seurat_object, sex_col, sample_col, all_sex_genes, output_dir)
  }

  seurat_object
}
