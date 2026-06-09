# =============================================================================
# scSidekick — Copy Number Variation inference
#
# RunInferCNV  — wrapper around the inferCNV package. Handles downsampling,
#                gene-order generation (via AnnotateFeatures / GetGeneOrder),
#                running the inference, post-processing the map_metadata file,
#                and producing standard figures. Returns the Seurat object with
#                new metadata columns added.
#
# RunCopyKAT   — wrapper around CopyKAT for fast, reference-optional CNV
#                calling. Useful as a complement or quick first-pass.
#
# Both functions accept caffeinate = TRUE (default FALSE) to prevent the
# laptop from sleeping during long runs.
# =============================================================================

# ── Internal helpers ──────────────────────────────────────────────────────────

.nk_default_chromosomes <- list(
  mouse     = as.character(1:19),
  human     = as.character(c(1:22, "X", "Y")),
  zebrafish = as.character(1:25)
)

# Auto-ensure gene-order annotations are present then return GetGeneOrder()
.nk_ensure_gene_order <- function(seurat_object, species, chromosomes,
                                   gene_biotype, annotation_file) {
  feat_meta <- seurat_object@assays[[Seurat::DefaultAssay(seurat_object)]]@meta.features
  if (!"nk_chr" %in% colnames(feat_meta)) {
    message("scSidekick: Running AnnotateFeatures() for gene order ...")
    seurat_object <- AnnotateFeatures(seurat_object,
      species         = species,
      annotation_file = annotation_file
    )
  }
  list(
    object     = seurat_object,
    gene_order = GetGeneOrder(seurat_object,
      chromosomes  = chromosomes,
      gene_biotype = gene_biotype
    )
  )
}

# Compute per-cell CNV score = mean((expression - 1)^2) from the inferCNV
# smoothed expression matrices.  Reference cells anchor the threshold:
#   threshold = mean(ref_scores) + cnv_k * sd(ref_scores)
# Falls back to proportion_cnv rowSums if the matrices are not found.
.nk_infercnv_postprocess <- function(seurat_object, out_dir,
                                      cnv_score_threshold, cnv_k) {

  # ── 1. Matrix-based score (preferred) ─────────────────────────────────────
  # inferCNV names the final file by pipeline step; match the last file found
  # so it works regardless of whether denoise/HMM steps ran.
  obs_files <- sort(list.files(out_dir, pattern = "observations\\.txt$",
                                full.names = TRUE))
  ref_files <- sort(list.files(out_dir, pattern = "references\\.txt$",
                                full.names = TRUE))

  used_matrix_method <- FALSE

  if (length(obs_files) > 0 && length(ref_files) > 0) {
    obs_path <- tail(obs_files, 1)
    ref_path <- tail(ref_files, 1)

    message("scSidekick: Computing CNV scores from\n  ", basename(obs_path),
            "\n  ", basename(ref_path))

    obs_mat <- tryCatch(
      as.matrix(read.delim(obs_path, row.names = 1, check.names = FALSE)),
      error = function(e) NULL
    )
    ref_mat <- tryCatch(
      as.matrix(read.delim(ref_path, row.names = 1, check.names = FALSE)),
      error = function(e) NULL
    )

    if (!is.null(obs_mat) && !is.null(ref_mat)) {
      # score = mean((expression - 1)^2) per cell  (Tirosh 2016 approach)
      obs_scores <- colMeans((obs_mat - 1)^2)
      ref_scores <- colMeans((ref_mat - 1)^2)

      if (identical(cnv_score_threshold, "auto")) {
        threshold <- mean(ref_scores) + cnv_k * sd(ref_scores)
        message("scSidekick: Auto threshold = mean(ref) + ", cnv_k, " × SD",
                "\n  ref mean = ", round(mean(ref_scores), 5),
                "  ref SD = ",   round(sd(ref_scores), 5),
                "  threshold = ", round(threshold, 5))
      } else {
        threshold <- as.numeric(cnv_score_threshold)
        message("scSidekick: Manual CNV score threshold: ", threshold)
      }

      all_scores <- c(obs_scores, ref_scores)
      score_df <- data.frame(
        CNV_score      = as.numeric(all_scores),
        CNV_Prediction = ifelse(all_scores > threshold, "Malignant", "Normal"),
        row.names      = names(all_scores),
        stringsAsFactors = FALSE
      )
      common <- intersect(rownames(score_df), colnames(seurat_object))
      seurat_object <- Seurat::AddMetaData(seurat_object,
        metadata = score_df[common, , drop = FALSE]
      )
      used_matrix_method <- TRUE
    }
  }

  # ── 2. Fallback: proportion_cnv rowSums ───────────────────────────────────
  if (!used_matrix_method) {
    warning("scSidekick: inferCNV matrix files not found — ",
            "falling back to proportion_cnv-based scoring.\n",
            "  This method uses an arbitrary threshold; ",
            "consider rerunning with output files intact.")

    meta_path <- file.path(out_dir, "map_metadata_from_infercnv.txt")
    if (file.exists(meta_path)) {
      mapmeta  <- read.delim(meta_path, row.names = 1,
                              stringsAsFactors = FALSE, check.names = FALSE)
      cnv_cols <- grep("^proportion_cnv_chr", colnames(mapmeta), value = TRUE)
      if (length(cnv_cols) > 0) {
        subcnv <- mapmeta[, cnv_cols, drop = FALSE]
        fb_threshold <- if (identical(cnv_score_threshold, "auto")) 4 else
          as.numeric(cnv_score_threshold)
        fb_score  <- log2(rowSums(subcnv, na.rm = TRUE) + 1)
        score_df  <- data.frame(
          CNV_score      = fb_score,
          CNV_Prediction = ifelse(fb_score > fb_threshold, "Malignant", "Normal"),
          row.names      = rownames(mapmeta),
          stringsAsFactors = FALSE
        )
        common <- intersect(rownames(score_df), colnames(seurat_object))
        seurat_object <- Seurat::AddMetaData(seurat_object,
          metadata = score_df[common, , drop = FALSE]
        )
      }
    }
  }

  # ── 3. Always add proportion_cnv columns + subcluster (useful for plots) ──
  meta_path <- file.path(out_dir, "map_metadata_from_infercnv.txt")
  if (file.exists(meta_path)) {
    mapmeta  <- read.delim(meta_path, row.names = 1,
                            stringsAsFactors = FALSE, check.names = FALSE)
    cnv_cols <- grep("^proportion_cnv_chr", colnames(mapmeta), value = TRUE)
    if (length(cnv_cols) > 0) {
      subcnv <- mapmeta[, cnv_cols, drop = FALSE]
      common <- intersect(rownames(subcnv), colnames(seurat_object))
      seurat_object <- Seurat::AddMetaData(seurat_object,
        metadata = subcnv[common, , drop = FALSE]
      )
    }
    if ("subcluster" %in% colnames(mapmeta)) {
      common <- intersect(rownames(mapmeta), colnames(seurat_object))
      seurat_object <- Seurat::AddMetaData(seurat_object,
        metadata = mapmeta[common, "subcluster", drop = FALSE]
      )
    }
  }

  seurat_object
}

# Standard CNV visualizations saved to output_dir
.nk_infercnv_plots <- function(seurat_object, cell_type_col, out_dir) {
  cnv_cols <- grep("^proportion_cnv_chr", colnames(seurat_object@meta.data), value = TRUE)
  if (length(cnv_cols) == 0) return(invisible(NULL))

  cell_colors <- .nk_colors(seurat_object, cell_type_col)
  pred_colors <- c(Malignant = "#ab3000", Normal = "#2f4b7c")

  # 1. Proportion CNV heatmap (ComplexHeatmap)
  tryCatch({
    submat  <- as.matrix(seurat_object@meta.data[, cnv_cols])
    colnames(submat) <- gsub("proportion_cnv_", "", colnames(submat))

    col_fun <- circlize::colorRamp2(
      c(0, 0.2, 0.3, 0.4, 0.6),
      c("white", "#feffe6", "#F9BF31", "#FF8B2B", "#ab3000")
    )
    ra <- ComplexHeatmap::rowAnnotation(
      df  = seurat_object@meta.data[, cell_type_col, drop = FALSE],
      col = if (!is.null(cell_colors)) list(setNames(list(cell_colors), cell_type_col)) else list()
    )
    ht <- ComplexHeatmap::Heatmap(
      submat,
      name              = "Proportion CNV",
      col               = col_fun,
      right_annotation  = ra,
      row_split         = seurat_object@meta.data[[cell_type_col]],
      show_row_names    = FALSE,
      show_row_dend     = FALSE,
      cluster_columns   = FALSE,
      row_title_rot     = 0,
      column_names_rot  = 45
    )
    grDevices::pdf(file.path(out_dir, "CNV_proportion_heatmap.pdf"), width = 12, height = 8)
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()
  }, error = function(e) warning("scSidekick: CNV heatmap failed: ", e$message))

  # 2. DimPlot — CNV prediction + cell type
  if ("CNV_Prediction" %in% colnames(seurat_object@meta.data) &&
      "reduction" %in% names(seurat_object@reductions)) {

    p <- Seurat::DimPlot(seurat_object, group.by = "CNV_Prediction",
                          cols = pred_colors, label = FALSE) +
         Seurat::DimPlot(seurat_object, group.by = cell_type_col,
                          cols = cell_colors, label = TRUE, label.size = 3) +
         patchwork::plot_layout(ncol = 2) &
         ggplot2::theme_classic(base_size = 11) &
         ggplot2::theme(legend.position = "bottom")

    grDevices::pdf(file.path(out_dir, "CNV_dimplot.pdf"), width = 12, height = 5)
    print(p)
    grDevices::dev.off()
  }

  # 3. FeaturePlot — CNV score
  if ("CNV_score" %in% colnames(seurat_object@meta.data)) {
    p_score <- Seurat::FeaturePlot(seurat_object, features = "CNV_score",
                                    order = TRUE) +
               ggplot2::scale_color_viridis_c(option = "plasma") +
               ggplot2::ggtitle("inferCNV score (log2 sum)")
    grDevices::pdf(file.path(out_dir, "CNV_score_featureplot.pdf"), width = 6, height = 5)
    print(p_score)
    grDevices::dev.off()
  }

  invisible(NULL)
}

# ── RunInferCNV ───────────────────────────────────────────────────────────────

#' Run inferCNV copy-number inference
#'
#' Wraps the full inferCNV workflow: gene-order table generation (via
#' [AnnotateFeatures()] / [GetGeneOrder()]), downsampling, object creation,
#' inference, post-processing of `map_metadata_from_infercnv.txt`, and
#' standard figure output. Returns the Seurat object with new metadata columns:
#' `CNV_score`, `CNV_Prediction`, `infercnv_subcluster`, and one
#' `proportion_cnv_chrN` column per chromosome.
#'
#' @param seurat_object A Seurat object.
#' @param cell_type_col Character. Metadata column containing cell-type labels
#'   used to split cells and define reference groups. Default `"Assignment"`.
#' @param ref_group_names Character vector. Cell-type labels treated as the
#'   normal (diploid) reference — typically immune or stromal populations.
#' @param output_dir Character or `NULL`. Directory for inferCNV output and
#'   saved figures. Walks up from [PrepObject()] settings when `NULL`.
#' @param species Character. `"mouse"`, `"human"`, or `"zebrafish"`. Used to
#'   auto-generate gene order via [AnnotateFeatures()] if not yet annotated.
#' @param chromosomes Character or numeric vector of chromosomes to analyze.
#'   `NULL` uses species-appropriate autosomes (mouse 1–19, human 1–22,
#'   zebrafish 1–25).
#' @param gene_biotype Character. Biotype filter for the gene order table.
#'   Default `"protein_coding"`.
#' @param annotation_file Character or `NULL`. Passed to [AnnotateFeatures()]
#'   when manual gene annotation is preferred over the automatic sources.
#' @param downsample Integer. Maximum cells per identity class. `NULL` uses all
#'   cells (not recommended for large objects). Default `200L`.
#' @param cutoff Numeric. inferCNV expression cutoff. Use `0.1` for 10x
#'   Genomics and `1` for Smart-seq2. Default `0.1`.
#' @param HMM Logical. Run inferCNV's HMM for subclonal CNV calls? Dramatically
#'   increases runtime. Default `FALSE`.
#' @param denoise Logical. Apply inferCNV denoising. Default `TRUE`.
#' @param scale_data Logical. Scale expression data before inference.
#'   Default `TRUE`.
#' @param cluster_by_groups Logical. Cluster cells within their annotation
#'   groups. Default `TRUE`.
#' @param num_threads Integer. Threads for inferCNV. Default `4L`.
#' @param cnv_score_threshold Numeric or `"auto"`. When `"auto"` (default),
#'   the threshold is set to `mean(ref_scores) + cnv_k × sd(ref_scores)` using
#'   the reference cells' CNV score distribution — this is the data-driven
#'   approach from Tirosh et al. 2016. Pass a numeric value to override.
#' @param cnv_k Numeric. Number of standard deviations above the reference mean
#'   used for the auto threshold. Higher values = more conservative (fewer
#'   malignant calls). Default `2`.
#' @param resume Logical. If `TRUE` and output files already exist, skip
#'   re-running inferCNV and go straight to post-processing. Default `FALSE`.
#' @param caffeinate Logical. Prevent the machine from sleeping during the run.
#'   Default `FALSE`.
#'
#' @return The input Seurat object with additional metadata columns.
#' @seealso [RunCopyKAT()], [AnnotateFeatures()], [GetGeneOrder()]
#' @export
RunInferCNV <- function(seurat_object,
                         cell_type_col       = "Assignment",
                         ref_group_names,
                         output_dir          = NULL,
                         species             = NULL,
                         chromosomes         = NULL,
                         gene_biotype        = "protein_coding",
                         annotation_file     = NULL,
                         downsample          = 200L,
                         cutoff              = 0.1,
                         HMM                 = FALSE,
                         denoise             = TRUE,
                         scale_data          = TRUE,
                         cluster_by_groups   = TRUE,
                         num_threads         = 4L,
                         cnv_score_threshold = "auto",
                         cnv_k               = 2,
                         resume              = FALSE,
                         caffeinate          = FALSE) {

  if (caffeinate) { .caff <- .nk_caffeinate(); on.exit(.nk_decaffeinate(.caff), add = TRUE) }

  if (!requireNamespace("infercnv", quietly = TRUE))
    stop("Package 'infercnv' is required.\n",
         "  BiocManager::install('infercnv')")

  # ── Resolve output directory ───────────────────────────────────────────────
  output_dir <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL
  if (is.null(output_dir))
    stop("'output_dir' must be supplied or stored via PrepObject(output_dir = ...).")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ── Resolve species ────────────────────────────────────────────────────────
  species <- species %||% seurat_object@misc$nk_annotation_species
  if (is.null(species))
    stop("'species' must be supplied (\"mouse\", \"human\", or \"zebrafish\").")
  species <- match.arg(species, names(.nk_default_chromosomes))

  chromosomes <- chromosomes %||% .nk_default_chromosomes[[species]]

  # ── Gene order table ───────────────────────────────────────────────────────
  ann_result   <- .nk_ensure_gene_order(seurat_object, species, chromosomes,
                                         gene_biotype, annotation_file)
  seurat_object <- ann_result$object
  gene_order    <- ann_result$gene_order

  # ── Downsample ─────────────────────────────────────────────────────────────
  if (!is.null(downsample)) {
    Seurat::Idents(seurat_object) <- seurat_object@meta.data[[cell_type_col]]
    seurat_object <- subset(seurat_object, downsample = as.integer(downsample))
    message("scSidekick: Downsampled to ", ncol(seurat_object), " cells (",
            downsample, " per group).")
  }

  # ── Align genes: keep only genes present in both object and gene order ─────
  shared_genes  <- intersect(rownames(seurat_object), rownames(gene_order))
  gene_order    <- gene_order[shared_genes, ]
  counts_matrix <- Seurat::GetAssayData(seurat_object, layer = "counts",
                                         assay = Seurat::DefaultAssay(seurat_object))
  counts_matrix <- counts_matrix[shared_genes, ]
  meta_sub      <- seurat_object@meta.data[, cell_type_col, drop = FALSE]

  message("scSidekick: ", length(shared_genes), " genes × ",
          ncol(counts_matrix), " cells → inferCNV.")

  # ── Skip inference if resuming ─────────────────────────────────────────────
  run_marker <- file.path(output_dir, "run.final.infercnv_obj")
  if (resume && file.exists(run_marker)) {
    message("scSidekick: Resume mode — skipping inferCNV run, loading existing output.")
  } else {
    # ── Create inferCNV object ───────────────────────────────────────────────
    infercnv_obj <- infercnv::CreateInfercnvObject(
      raw_counts_matrix = counts_matrix,
      annotations_file  = meta_sub,
      delim             = "\t",
      gene_order_file   = gene_order,
      ref_group_names   = ref_group_names
    )

    # inferCNV's add_to_seurat() expects Seurat v3 slot layout
    old_assay_opt <- getOption("Seurat.object.assay.version")
    options(Seurat.object.assay.version = "v3")
    on.exit(options(Seurat.object.assay.version = old_assay_opt), add = TRUE)

    # ── Run ─────────────────────────────────────────────────────────────────
    infercnv_obj <- infercnv::run(
      infercnv_obj,
      cutoff            = cutoff,
      out_dir           = output_dir,
      scale_data        = scale_data,
      cluster_by_groups = cluster_by_groups,
      denoise           = denoise,
      HMM               = HMM,
      num_threads       = as.integer(num_threads),
      output_format     = "pdf",
      reassignCNVs      = TRUE
    )
  }

  # ── Add inferCNV results back to Seurat object ────────────────────────────
  tryCatch(
    seurat_object <- infercnv::add_to_seurat(
      seurat_obj   = seurat_object,
      assay_name   = Seurat::DefaultAssay(seurat_object),
      infercnv_output_path = output_dir,
      top_n        = 10,
      bp_tolerance = 2e6
    ),
    error = function(e)
      warning("scSidekick: add_to_seurat() failed: ", e$message,
              "\n  Post-processing will use map_metadata_from_infercnv.txt directly.")
  )

  # ── Post-process map_metadata ─────────────────────────────────────────────
  seurat_object <- .nk_infercnv_postprocess(seurat_object, output_dir,
                                              cnv_score_threshold, cnv_k)

  # ── Standard figures ──────────────────────────────────────────────────────
  .nk_infercnv_plots(seurat_object, cell_type_col, output_dir)

  n_mal <- if ("CNV_Prediction" %in% colnames(seurat_object@meta.data))
    sum(seurat_object$CNV_Prediction == "Malignant", na.rm = TRUE) else NA
  message("scSidekick: inferCNV complete. Malignant cells predicted: ", n_mal,
          if (identical(cnv_score_threshold, "auto"))
            paste0("\n  Threshold was auto-set. Adjust cnv_k (currently ", cnv_k,
                   ") to tune sensitivity.")
          else
            paste0("\n  Manual threshold: ", cnv_score_threshold,
                   ". Inspect CNV_score distribution to validate."))

  seurat_object
}

# ── RunCopyKAT ────────────────────────────────────────────────────────────────

#' Run CopyKAT copy-number inference
#'
#' Fast, reference-optional CNV calling using the CopyKAT Bayesian model.
#' Useful as a quick first-pass or to validate [RunInferCNV()] results.
#' Adds `copykat_prediction` (`"aneuploid"` / `"diploid"` / `"not.defined"`)
#' and `copykat_score` to Seurat metadata.
#'
#' @param seurat_object A Seurat object.
#' @param cell_type_col Character. Metadata column used to color output plots.
#'   Default `"Assignment"`.
#' @param ref_cells Character vector of cell-type labels to use as the diploid
#'   reference. `NULL` runs in unsupervised mode (CopyKAT infers normals
#'   automatically). Default `NULL`.
#' @param output_dir Character or `NULL`. Directory for CopyKAT output files
#'   and plots. Walks up from [PrepObject()] when `NULL`.
#' @param n_cores Integer. Parallel threads. Default `4L`.
#' @param LOW_DR,UP_DR Numeric. Lower and upper bounds on the dropout rate
#'   filter passed to CopyKAT. Defaults `0.05` / `0.1`.
#' @param win.size Integer. Genomic window size (number of genes) for CNV
#'   smoothing. Default `25`.
#' @param caffeinate Logical. Prevent the machine from sleeping. Default `FALSE`.
#'
#' @return The input Seurat object with `copykat_prediction` and `copykat_score`
#'   added to metadata.
#' @seealso [RunInferCNV()]
#' @export
RunCopyKAT <- function(seurat_object,
                        cell_type_col = "Assignment",
                        ref_cells     = NULL,
                        output_dir    = NULL,
                        n_cores       = 4L,
                        LOW_DR        = 0.05,
                        UP_DR         = 0.1,
                        win.size      = 25L,
                        caffeinate    = FALSE) {

  if (caffeinate) { .caff <- .nk_caffeinate(); on.exit(.nk_decaffeinate(.caff), add = TRUE) }

  if (!requireNamespace("copykat", quietly = TRUE))
    stop("Package 'copykat' is required.\n",
         "  devtools::install_github('navinlabcode/copykat')")

  # ── Resolve output directory ───────────────────────────────────────────────
  output_dir <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL
  if (is.null(output_dir))
    stop("'output_dir' must be supplied or stored via PrepObject(output_dir = ...).")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ── Reference cell barcodes ────────────────────────────────────────────────
  ref_barcodes <- if (!is.null(ref_cells)) {
    rownames(seurat_object@meta.data)[
      seurat_object@meta.data[[cell_type_col]] %in% ref_cells
    ]
  } else ""   # empty string = unsupervised

  counts_matrix <- Seurat::GetAssayData(seurat_object, layer = "counts",
                                         assay = Seurat::DefaultAssay(seurat_object))

  # ── Run CopyKAT ───────────────────────────────────────────────────────────
  message("scSidekick: Running CopyKAT (",
          if (identical(ref_barcodes, "")) "unsupervised" else
            paste(length(ref_barcodes), "reference cells"), ") ...")

  ck <- copykat::copykat(
    rawmat     = as.matrix(counts_matrix),
    id.type    = "S",
    cell.line  = "no",
    ngene.chr  = 5,
    win.size   = as.integer(win.size),
    KS.cut     = 0.1,
    sam.name   = file.path(output_dir, "copykat"),
    distance   = "euclidean",
    norm.cell.names = ref_barcodes,
    n.cores    = as.integer(n_cores),
    LOW.DR     = LOW_DR,
    UP.DR      = UP_DR
  )

  # ── Add predictions to Seurat metadata ────────────────────────────────────
  if (!is.null(ck$prediction) && nrow(ck$prediction) > 0) {
    pred_df <- data.frame(
      copykat_prediction = ck$prediction$copykat.pred,
      row.names          = ck$prediction$cell.names,
      stringsAsFactors   = FALSE
    )
    seurat_object <- Seurat::AddMetaData(seurat_object, metadata = pred_df)
  }

  # CNA matrix row sums as a continuous score
  if (!is.null(ck$CNAmat) && ncol(ck$CNAmat) > 1) {
    cna_cells <- colnames(ck$CNAmat)[-1]   # first col is genomic position info
    cna_score <- colSums(abs(t(ck$CNAmat[, -1, drop = FALSE])))
    score_df  <- data.frame(
      copykat_score = as.numeric(cna_score),
      row.names     = cna_cells,
      stringsAsFactors = FALSE
    )
    seurat_object <- Seurat::AddMetaData(seurat_object, metadata = score_df)
  }

  # ── Plots ──────────────────────────────────────────────────────────────────
  if ("copykat_prediction" %in% colnames(seurat_object@meta.data)) {
    cell_colors <- .nk_colors(seurat_object, cell_type_col)
    pred_colors <- c(aneuploid = "#ab3000", diploid = "#2f4b7c",
                     not.defined = "gray70")

    p <- Seurat::DimPlot(seurat_object, group.by = "copykat_prediction",
                          cols = pred_colors, label = FALSE) +
         Seurat::DimPlot(seurat_object, group.by = cell_type_col,
                          cols = cell_colors, label = TRUE, label.size = 3) +
         patchwork::plot_layout(ncol = 2) &
         ggplot2::theme_classic(base_size = 11) &
         ggplot2::theme(legend.position = "bottom")

    grDevices::pdf(file.path(output_dir, "CopyKAT_dimplot.pdf"), width = 12, height = 5)
    print(p)
    grDevices::dev.off()
  }

  n_aneu <- if ("copykat_prediction" %in% colnames(seurat_object@meta.data))
    sum(seurat_object$copykat_prediction == "aneuploid", na.rm = TRUE) else NA
  message("scSidekick: CopyKAT complete. Aneuploid cells: ", n_aneu)

  seurat_object
}
