# =============================================================================
# scSidekick sample loading helpers
#
# LoadSamplesRNA()  - per-sample 10X loading, QC (doublets + miQC), and merge
#                     for scRNA-seq / snRNA-seq datasets.
# PlotQCMetrics()   - QC visualisation for the merged (pre-filtration) object:
#                     mito/ribo ratios, density plots, violin plots, barplots.
#
# Supports two metadata styles:
#   1. data.frame (recommended) - one row per sample; every column is added to
#      Seurat @meta.data automatically.
#   2. named lists / NULL - legacy parallel-vector style.
#
# Path construction uses {placeholder} templates so the function works with
# any Cellranger output layout.
# =============================================================================

# Internal: substitute {name} placeholders in a path string
.fill_path <- function(template, ...) {
  vals   <- list(...)
  result <- template
  for (nm in names(vals))
    result <- gsub(paste0("{", nm, "}"), vals[[nm]], result, fixed = TRUE)
  result
}

# Internal: auto-convert "MYSC220" → "MYSC_220" (insert _ before trailing digits)
.auto_run_id_x <- function(run_id) {
  gsub("([A-Za-z]+)([0-9]+)$", "\\1_\\2", run_id)
}


#' Load, QC, and merge 10X scRNA-seq samples
#'
#' Iterates over a list of samples, reads the 10X count data, optionally
#' adds per-cell annotation CSVs, attaches per-sample metadata, runs
#' DoubletFinder and miQC quality control, saves per-sample QC plots,
#' and finally merges all objects into a single Seurat object.
#'
#' @section Input methods:
#' Three ways to point the function at your data:
#' \enumerate{
#'   \item \strong{h5 files via templates} (default, \code{input_type = "h5"}):
#'     paths built from \code{base_dir} + \code{h5_template} using
#'     \code{Read10X_h5()}. Matches Cellranger-multi output.
#'   \item \strong{10X matrix folders via templates}
#'     (\code{input_type = "10x_dir"}): paths built from \code{base_dir} +
#'     \code{matrix_dir_template} and read with \code{Read10X()}. Each folder
#'     must contain \code{barcodes.tsv(.gz)}, \code{features.tsv(.gz)}, and
#'     \code{matrix.mtx(.gz)} (the Cellranger \code{filtered_feature_bc_matrix/}
#'     directory).
#'   \item \strong{Explicit paths} (\code{data_paths = ...}): a character vector
#'     of exact h5 files or matrix folders, one per sample, used verbatim.
#'     Use this when your samples are not organized in a consistent layout.
#'     \code{base_dir} and the templates are then ignored for path building.
#' }
#'
#' @section Metadata (recommended - `metadata` data.frame):
#' Pass a `data.frame` with one row per sample.  The column named by
#' `sample_id_col` must match `sample_ids` exactly.  Every other column is
#' added to `@@meta.data` automatically.
#'
#' ```r
#' sample_meta <- data.frame(
#'   sample_id = c("HC1","HC2","AD1","FTD1"),
#'   run_id    = c("MYSC220","MYSC221","MYSC222","MYSC223"),
#'   patient   = c("Firstable","Cododge","Transferism","Promiseable"),
#'   group     = c("HC","HC","AD","FTD"),
#'   sex       = c("M","F","M","F"),
#'   batch     = c(1L, 1L, 1L, 1L)
#' )
#' ```
#'
#' @section Path templates:
#' Paths are built from `base_dir`, `run_id`, and `run_id_x` using
#' `{placeholder}` substitution.  The defaults match Cellranger-multi output:
#' ```
#' {base_dir}/{run_id}/per_sample_outs/{run_id_x}/count/
#'   sample_filtered_feature_bc_matrix.h5
#' ```
#' Override `h5_template` or `cellanno_template` for other layouts.
#'
#' @param sample_ids Character vector of human-readable sample IDs
#'   (e.g. `c("HC1","HC2","AD1")`). Used as Seurat project names and as the
#'   key column when matching `metadata`.
#' @param run_ids Character vector of Cellranger run IDs corresponding to
#'   `sample_ids` (e.g. `paste0("MYSC", 220:223)`). Used to build h5 / matrix
#'   paths and as cell-id prefixes on merge. `NULL` (allowed only when
#'   `data_paths` is supplied) falls back to `sample_ids`.
#' @param input_type Character. How to read each sample's counts:
#'   `"h5"` (default) uses `Seurat::Read10X_h5()`; `"10x_dir"` uses
#'   `Seurat::Read10X()` on a matrix folder. Ignored if `data_paths` already
#'   points at the right kind of source (the type is still used to pick the
#'   reader, so set it to match your `data_paths`).
#' @param data_paths Character vector or `NULL`. Exact path to each sample's
#'   data source (an h5 file when `input_type = "h5"`, or a matrix folder when
#'   `input_type = "10x_dir"`), one entry per `sample_ids`. When supplied,
#'   `base_dir` and the path templates are bypassed entirely. Default `NULL`.
#' @param run_ids_x Character vector or `NULL`. The "underscored" form of
#'   `run_ids` used in `per_sample_outs/` subdirectory names
#'   (e.g. `paste0("MYSC_", 220:223)`). `NULL` auto-generates by inserting
#'   `_` before the trailing digit block: `"MYSC220"` → `"MYSC_220"`.
#' @param metadata Data frame or `NULL`. One row per sample; must contain a
#'   column matching `sample_id_col`. All other columns are attached to
#'   `@@meta.data`. `NULL` skips this step (metadata must be added manually
#'   afterwards).
#' @param sample_id_col Character. Name of the column in `metadata` that
#'   matches `sample_ids`. Default `"sample_id"`.
#' @param base_dir Character or `NULL`. Root directory containing the
#'   Cellranger run folders (one subfolder per `run_id`). Required unless
#'   `data_paths` is supplied.
#' @param h5_template Character. Path template for the h5 file
#'   (`input_type = "h5"`). Available placeholders: `{base_dir}`, `{run_id}`,
#'   `{run_id_x}`. Default matches Cellranger-multi per-sample output layout.
#' @param matrix_dir_template Character. Path template for the 10X matrix
#'   folder (`input_type = "10x_dir"`). Same placeholders as `h5_template`.
#'   Default `"{base_dir}/{run_id}/filtered_feature_bc_matrix"`.
#' @param cellanno_template Character or `NULL`. Path template for the
#'   per-cell annotation CSV (e.g. from Cellranger cell-type calling).
#'   `NULL` disables cell annotation loading. Files that do not exist are
#'   silently skipped per sample - it is safe to provide a template even when
#'   only some samples have the file.
#' @param output_dir Character. Directory for per-sample QC PDF plots.
#' @param robj_dir Character or `NULL`. Directory for saving the final merged
#'   Seurat object. `NULL` skips saving.
#' @param merged_filename Character. Filename for the saved merged object.
#'   Default `"merged_object.rds"`.
#' @param save_individual Logical. Save each fully-processed sample as its own
#'   `<sample_id>.rds` in `individual_dir` as it finishes. Combined with
#'   `resume = TRUE`, this lets a long load survive interruptions. Default
#'   `FALSE`.
#' @param resume Logical. Before processing a sample, check `individual_dir`
#'   for an existing `<sample_id>.rds`; if found, load it and skip all
#'   re-processing (read, QC, doublets, miQC) for that sample. Default `FALSE`.
#' @param individual_dir Character or `NULL`. Directory for the per-sample
#'   `.rds` files used by `save_individual` / `resume`. `NULL` defaults to
#'   `<robj_dir or output_dir>/Individual_Samples`.
#' @param min.features Integer. Minimum features per cell for
#'   `CreateSeuratObject()`. Default `100`.
#' @param mt_pattern Character. Regex pattern for mitochondrial genes.
#'   Default `"^MT-"` (human). Use `"^mt-"` for mouse.
#' @param run_doublet Logical. Run `SCP::db_scDblFinder()` doublet detection?
#'   Default `TRUE`.
#' @param doublet_rate Numeric or `NULL`. Expected doublet rate per 1 000
#'   cells loaded. `NULL` auto-computes as `ncol(seurat_object) / 1000 * 0.01`.
#' @param run_miqc Logical. Run `SeuratWrappers::RunMiQC()`? Default `TRUE`.
#' @param miqc_posterior Numeric. miQC posterior probability cutoff. Cells
#'   above this are flagged as low quality. Default `0.95`.
#' @param force_assay_v3 Logical. Convert the RNA assay to Seurat v3 format
#'   (`as(seurat_object[["RNA"]], "Assay")`) after creation? Needed for compatibility
#'   with older packages. Default `TRUE`.
#' @param merge_samples Logical. Merge all sample objects into one Seurat
#'   object at the end? Default `TRUE`.
#' @param add_cell_ids Character vector or `NULL`. Prefix added to each cell
#'   barcode on merge. `NULL` uses `run_ids`.
#' @param modality Character. Sequencing modality label stored in the merged
#'   object's metadata and `@misc$scSidekick_params`. One of `"scRNA-seq"`
#'   (default) or `"snRNA-seq"`.
#'
#' @return If `merge_samples = TRUE`, the merged Seurat object (invisibly).
#'   If `merge_samples = FALSE`, a named list of individual Seurat objects.
#' @export
LoadSamplesRNA <- function(
    sample_ids,
    run_ids            = NULL,
    run_ids_x          = NULL,
    input_type         = c("h5", "10x_dir"),
    data_paths         = NULL,
    metadata           = NULL,
    sample_id_col      = "sample_id",
    base_dir           = NULL,
    h5_template        = paste0("{base_dir}/{run_id}/per_sample_outs/",
                                "{run_id_x}/count/",
                                "sample_filtered_feature_bc_matrix.h5"),
    matrix_dir_template = "{base_dir}/{run_id}/filtered_feature_bc_matrix",
    cellanno_template  = paste0("{base_dir}/{run_id}/per_sample_outs/",
                                "{run_id_x}/count/",
                                "cell_types/cell_types.csv"),
    output_dir,
    robj_dir           = NULL,
    merged_filename    = "merged_object.rds",
    save_individual    = FALSE,
    resume             = FALSE,
    individual_dir     = NULL,
    min.features       = 100L,
    mt_pattern         = "^MT-",
    run_doublet        = TRUE,
    doublet_rate       = NULL,
    run_miqc           = TRUE,
    miqc_posterior     = 0.95,
    force_assay_v3     = TRUE,
    merge_samples      = TRUE,
    add_cell_ids       = NULL,
    modality           = c("scRNA-seq", "snRNA-seq")
) {
  modality   <- match.arg(modality)
  input_type <- match.arg(input_type)

  # ── Validate / prepare ──────────────────────────────────────────────────────
  # run_ids are optional when explicit data_paths are supplied (fall back to
  # the human-readable sample_ids for cell-id prefixes etc.).
  if (is.null(run_ids)) {
    if (is.null(data_paths))
      stop("`run_ids` is required unless you supply explicit `data_paths`.")
    run_ids <- sample_ids
  }
  stopifnot(length(sample_ids) == length(run_ids))

  # Explicit paths bypass templates; otherwise base_dir is required.
  use_explicit_paths <- !is.null(data_paths)
  if (use_explicit_paths) {
    if (length(data_paths) != length(sample_ids))
      stop("`data_paths` must have one entry per sample (",
           length(sample_ids), " expected, got ", length(data_paths), ").")
  } else if (is.null(base_dir)) {
    stop("`base_dir` is required unless you supply explicit `data_paths`.")
  }

  if (is.null(run_ids_x))
    run_ids_x <- .auto_run_id_x(run_ids)
  stopifnot(length(run_ids_x) == length(run_ids))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(robj_dir))
    dir.create(robj_dir, recursive = TRUE, showWarnings = FALSE)

  # Per-sample RDS directory for save_individual / resume
  if (save_individual || resume) {
    if (is.null(individual_dir))
      individual_dir <- file.path(robj_dir %||% output_dir, "Individual_Samples")
    dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Pre-process metadata data.frame (keyed by sample_id_col)
  meta_df <- if (!is.null(metadata) && is.data.frame(metadata)) {
    if (!sample_id_col %in% colnames(metadata))
      stop("'", sample_id_col, "' column not found in `metadata`.")
    rownames(metadata) <- metadata[[sample_id_col]]
    metadata
  } else NULL

  # Columns to skip when attaching metadata (sample_id_col + run_id already handled)
  skip_cols <- c(sample_id_col, "run_id", "run_id_x")

  # Check optional packages
  if (run_doublet && !requireNamespace("SCP", quietly = TRUE)) {
    warning("Package 'SCP' not installed - skipping doublet detection.")
    run_doublet <- FALSE
  }
  if (run_miqc &&
      (!requireNamespace("SeuratWrappers", quietly = TRUE) ||
       !requireNamespace("flexmix",        quietly = TRUE))) {
    warning("Package 'SeuratWrappers' / 'flexmix' not installed - skipping miQC.")
    run_miqc <- FALSE
  }

  # ── Per-sample loop ─────────────────────────────────────────────────────────
  objects_list <- vector("list", length(sample_ids))
  names(objects_list) <- sample_ids

  for (i in seq_along(sample_ids)) {
    sid   <- sample_ids[i]
    rid   <- run_ids[i]
    ridx  <- run_ids_x[i]

    message("\n========== ", sid, " (", rid, ") ==========")

    # ── 0. Resume: load cached per-sample object and skip re-processing ───────
    ind_path <- if (save_individual || resume)
      file.path(individual_dir, paste0(sid, ".rds")) else NULL
    if (resume && !is.null(ind_path) && file.exists(ind_path)) {
      message("  Resuming - loading cached object: ", basename(ind_path))
      objects_list[[sid]] <- tryCatch(
        readRDS(ind_path),
        error = function(e) {
          warning("Cached object for ", sid, " is unreadable, reprocessing: ",
                  conditionMessage(e)); NULL
        }
      )
      if (!is.null(objects_list[[sid]])) next   # only skip if the load succeeded
    }

    # ── 1. Read 10X counts (h5 file or matrix folder) ─────────────────────────
    # Resolve the path: explicit data_paths win; otherwise build from template.
    data_path <- if (use_explicit_paths) {
      data_paths[i]
    } else if (input_type == "h5") {
      .fill_path(h5_template, base_dir = base_dir,
                 run_id = rid, run_id_x = ridx)
    } else {
      .fill_path(matrix_dir_template, base_dir = base_dir,
                 run_id = rid, run_id_x = ridx)
    }

    # file.exists() is TRUE for both files and directories
    if (!file.exists(data_path)) {
      warning(if (input_type == "h5") "h5 file" else "matrix folder",
              " not found - skipping ", sid, ":\n  ", data_path)
      next
    }

    mat <- tryCatch(
      if (input_type == "h5")
        Seurat::Read10X_h5(data_path)
      else
        Seurat::Read10X(data.dir = data_path),
      error = function(e) {
        warning("Failed to read counts for ", sid, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(mat)) next

    # Read10X_h5 may return a list (e.g. multiome); take the Gene Expression slot
    if (is.list(mat)) {
      if ("Gene Expression" %in% names(mat)) {
        message("  h5 contains multiple modalities - using 'Gene Expression'.")
        mat <- mat[["Gene Expression"]]
      } else {
        mat <- mat[[1]]
      }
    }

    seurat_object <- Seurat::CreateSeuratObject(counts = mat,
                                      min.features = min.features,
                                      project      = sid)
    message("  Cells loaded: ", ncol(seurat_object))

    # ── 2. Cell annotation CSV (optional) ─────────────────────────────────────
    # Skipped in explicit-paths mode (no base_dir to build the template from).
    if (!is.null(cellanno_template) && !is.null(base_dir)) {
      ca_path <- .fill_path(cellanno_template, base_dir = base_dir,
                            run_id = rid, run_id_x = ridx)
      if (file.exists(ca_path)) {
        message("  Adding cell annotations from: ", basename(ca_path))
        ca <- tryCatch(
          utils::read.csv(ca_path, row.names = 1),
          error = function(e) {
            warning("Could not read cell annotation for ", sid, ": ",
                    conditionMessage(e)); NULL
          }
        )
        if (!is.null(ca))
          seurat_object <- Seurat::AddMetaData(seurat_object, metadata = ca)
      }
      # silently skip if file absent - per-sample annotations are optional
    }

    # ── 3. Per-sample metadata ─────────────────────────────────────────────────
    if (!is.null(meta_df) && sid %in% rownames(meta_df)) {
      row      <- meta_df[sid, , drop = FALSE]
      add_cols <- setdiff(colnames(row), skip_cols)
      for (col in add_cols)
        seurat_object[[col]] <- row[[col]]
    }
    seurat_object$Sample <- sid   # always add Sample column

    # ── 4. Normalize → PCA ────────────────────────────────────────────────────
    message("  Normalizing and running PCA...")
    seurat_object <- Seurat::NormalizeData(seurat_object, normalization.method = "LogNormalize",
                                 scale.factor = 10000, verbose = FALSE)
    seurat_object <- Seurat::FindVariableFeatures(seurat_object, nfeatures = 3000,
                                        selection.method = "vst", verbose = FALSE)
    seurat_object <- Seurat::ScaleData(seurat_object, verbose = FALSE)
    seurat_object <- Seurat::RunPCA(seurat_object, npcs = 50, verbose = FALSE)

    # ── 4b. Convert to Seurat v3 assay (compatibility) ────────────────────────
    if (force_assay_v3) {
      tryCatch(
        seurat_object[["RNA"]] <- as(object = seurat_object[["RNA"]], Class = "Assay"),
        error = function(e)
          message("  Note: assay v3 conversion skipped (", conditionMessage(e), ")")
      )
    }

    # ── 5. Doublet detection ──────────────────────────────────────────────────
    if (run_doublet) {
      message("  Running doublet detection...")
      db_rate <- if (is.null(doublet_rate))
        ncol(seurat_object) / 1000 * 0.01 else doublet_rate
      seurat_object <- tryCatch(
        SCP::db_scDblFinder(srt = seurat_object, assay = "RNA", db_rate = db_rate),
        error = function(e) {
          warning("db_scDblFinder failed for ", sid, ": ", conditionMessage(e))
          seurat_object
        }
      )
      # Standardise column name
      if ("db.scDblFinder_class" %in% colnames(seurat_object@meta.data))
        colnames(seurat_object@meta.data)[
          colnames(seurat_object@meta.data) == "db.scDblFinder_class"] <- "DoubletStatus"

      # Doublet DimPlot
      f_dbl <- file.path(output_dir, paste0(sid, " DoubletStatus plot.pdf"))
      if ("DoubletStatus" %in% colnames(seurat_object@meta.data)) {
        grDevices::pdf(f_dbl, width = 12, height = 4.5)
        try(print(Seurat::DimPlot(seurat_object, split.by = "DoubletStatus",
                                  order = TRUE, shuffle = TRUE)),
            silent = TRUE)
        grDevices::dev.off()
        .write_legend_sidecar(f_dbl, paste0(
          "UMAP of sample ", sid, " split by doublet classification from ",
          "SCP::db_scDblFinder. Each panel shows the cells assigned to one ",
          "DoubletStatus class (singlet vs. doublet); points are color-",
          "matched to the active cluster identity and plotted in shuffled ",
          "order so neither class is occluded."
        ))
        message("  Saved: ", f_dbl)
      }
    }

    # ── 6. Mitochondrial % + miQC ─────────────────────────────────────────────
    n_mt_genes <- sum(grepl(mt_pattern, rownames(seurat_object)))
    seurat_object[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = mt_pattern)

    # Guard: if mt_pattern matches no genes, percent.mt is all zero. miQC's
    # flexmix model then cannot fit (Log-likelihood: NaN) and silently falls
    # back to a percentile, producing a meaningless miQC.keep column. The usual
    # cause is a species mismatch: human mito genes are "^MT-" (uppercase),
    # mouse are "^mt-" (lowercase). Warn clearly and skip miQC for this sample.
    run_miqc_sample <- run_miqc
    if (n_mt_genes == 0) {
      # Detect whether a case-flipped pattern would have matched, so we can name
      # the correct one (human "^MT-" vs mouse "^mt-").
      alt <- if (grepl("MT", mt_pattern)) sub("MT", "mt", mt_pattern)
             else sub("mt", "MT", mt_pattern)
      alt_n  <- sum(grepl(alt, rownames(seurat_object)))
      suggest <- if (alt_n > 0)
        paste0(" Use mt_pattern = '", alt, "' (matches ", alt_n,
               " genes in this object).")
      else " No mitochondrial genes found under either case."
      warning("mt_pattern '", mt_pattern, "' matched 0 genes in ", sid,
              " - percent.mt is all zero, so miQC cannot fit and is skipped.",
              suggest, call. = FALSE)
      message("  [skip] miQC for ", sid, ": mt_pattern '", mt_pattern,
              "' matched 0 genes.",
              if (alt_n > 0) paste0(" Try mt_pattern = '", alt, "'.") else "")
      run_miqc_sample <- FALSE
    }

    if (run_miqc_sample) {
      message("  Running miQC (", n_mt_genes, " mito genes)...")
      seurat_object <- tryCatch(
        SeuratWrappers::RunMiQC(seurat_object,
                                percent.mt    = "percent.mt",
                                nFeature_RNA  = "nFeature_RNA",
                                posterior.cutoff = miqc_posterior,
                                model.slot    = "flexmix_model"),
        error = function(e) {
          warning("RunMiQC failed for ", sid, ": ", conditionMessage(e))
          seurat_object
        }
      )
    }

    # ── 7. Save per-sample object (so resume can skip it next time) ───────────
    if (save_individual && !is.null(ind_path)) {
      saveRDS(seurat_object, ind_path)
      message("  Saved per-sample object: ", ind_path)
    }

    objects_list[[sid]] <- seurat_object
    message("  Done: ", sid)
  }

  # Remove NULLs (samples that failed to load)
  objects_list <- Filter(Negate(is.null), objects_list)
  n_ok <- length(objects_list)
  message("\n", n_ok, "/", length(sample_ids),
          " samples loaded successfully.")

  if (n_ok == 0) {
    warning("No samples loaded - returning NULL.")
    return(invisible(NULL))
  }

  if (!merge_samples)
    return(invisible(objects_list))

  # ── Merge ──────────────────────────────────────────────────────────────────
  message("Merging ", n_ok, " samples...")
  ids_for_merge <- if (is.null(add_cell_ids))
    run_ids[sample_ids %in% names(objects_list)]
  else
    add_cell_ids

  merged <- if (n_ok == 1) {
    objects_list[[1]]
  } else {
    # merge() is a base S3 generic; the Seurat method (merge.Seurat) lives in
    # SeuratObject and is dispatched via the generic - there is no exported
    # Seurat::merge. Call the generic so dispatch picks merge.Seurat.
    merge(
      objects_list[[1]],
      y            = objects_list[-1],
      add.cell.ids = ids_for_merge
    )
  }
  message("Merged object: ", ncol(merged), " cells, ", nrow(merged), " genes.")

  # ── Store modality ─────────────────────────────────────────────────────────
  merged$modality <- modality
  merged@misc$scSidekick_params <- c(
    merged@misc$scSidekick_params,
    list(modality = modality)
  )

  # ── Save ───────────────────────────────────────────────────────────────────
  if (!is.null(robj_dir)) {
    save_path <- file.path(robj_dir, merged_filename)
    saveRDS(merged, save_path)
    message("Saved merged object: ", save_path)
  }

  invisible(merged)
}


# =============================================================================
# PlotQCMetrics()
# =============================================================================

#' Compute QC metrics and generate pre-filtration QC plots
#'
#' Adds `mitoRatio`, `riboRatio`, and `log10GenesPerUMI` to the object, then
#' saves four PDFs:
#' \enumerate{
#'   \item **Doublet barplot** - cells per sample colored by `DoubletStatus`
#'     (only if that column exists from `SCP::db_scDblFinder`).
#'   \item **miQC keep barplot** - cells per sample colored by `miQC.keep`
#'     (only if that column exists from `SeuratWrappers::RunMiQC`).
#'   \item **QC density + violin plots** - four metrics side-by-side with
#'     user-defined cutoff guidelines.
#'   \item **Cell count barplot** - cells per sample before any filtering.
#' }
#'
#' @param seurat_object Seurat object (merged, unfiltered).
#' @param sample_col Character. Metadata column holding sample identity.
#'   Default `"Sample"`.
#' @param sample_levels Character vector or `NULL`. Desired factor level order
#'   for `sample_col`. `NULL` leaves the current order unchanged.
#' @param group_col Character or `NULL`. Column used to color bars in the
#'   doublet and miQC barplots (e.g. `"orig.ident"`). `NULL` uses
#'   `sample_col`.
#' @param output_dir Character. Directory for PDF output.
#' @param object_name,subset_name Character. Prefix components for PDF filenames.
#' @param species Character. Controls the default mitochondrial and ribosomal
#'   gene patterns. One of:
#'   \itemize{
#'     \item `"human"` - `mt_pattern = "^MT-"`, `ribo_pattern = "^RP[LS]"`
#'     \item `"mouse"` - `mt_pattern = "^mt-"`, `ribo_pattern = "^Rp[ls]"`
#'   }
#' @param mt_pattern Character or `NULL`. Override the mitochondrial gene
#'   regex. `NULL` derives from `species`.
#' @param ribo_pattern Character or `NULL`. Override the ribosomal gene regex.
#'   `NULL` derives from `species`.
#' @param assay Character. Seurat assay for `PercentageFeatureSet()`.
#'   Default `"RNA"`.
#' @param sample_colors Named character vector mapping `sample_col` levels to
#'   colors. `NULL` auto-generates from `Nour_pal()`.
#' @param count_vlines Numeric vector. Vertical guide lines on the nCount_RNA
#'   density plot. Default `c(500, 50000)`.
#' @param feature_vlines Numeric vector. Guide lines on the nFeature_RNA
#'   density plot. Default `c(300, 80000)`.
#' @param mito_vline Numeric. Guide line on the mitoRatio density plot.
#'   Default `0.1`.
#' @param complexity_vline Numeric. Guide line on the log10GenesPerUMI
#'   density plot. Default `0.8`.
#'
#' @return The modified Seurat object with `mitoRatio`, `riboRatio`, and
#'   `log10GenesPerUMI` added to `@@meta.data` (invisibly).
#' @export
PlotQCMetrics <- function(
    seurat_object,
    sample_col        = "Sample",
    sample_levels     = NULL,
    group_col         = NULL,
    output_dir,
    object_name          = "",
    subset_name       = "",
    species           = "human",
    mt_pattern        = NULL,
    ribo_pattern      = NULL,
    assay             = "RNA",
    sample_colors     = NULL,
    count_vlines      = c(500, 50000),
    feature_vlines    = c(300, 80000),
    mito_vline        = 0.1,
    complexity_vline  = 0.8
) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file_prefix <- trimws(paste(object_name, subset_name))
  pref        <- function(fname) file.path(output_dir, paste0(file_prefix, " ", fname))

  # ── Species-specific gene patterns ─────────────────────────────────────────
  .mt_ribo <- list(
    human = list(mt = "^MT-",  ribo = "^RP[LS]"),
    mouse = list(mt = "^mt-",  ribo = "^Rp[ls]")
  )
  sp_key <- tolower(species)
  if (!sp_key %in% names(.mt_ribo))
    stop("species must be 'human' or 'mouse'. Got: '", species, "'.")
  if (is.null(mt_pattern))   mt_pattern   <- .mt_ribo[[sp_key]]$mt
  if (is.null(ribo_pattern)) ribo_pattern <- .mt_ribo[[sp_key]]$ribo
  if (is.null(group_col))    group_col    <- sample_col

  # ── Factor ordering ─────────────────────────────────────────────────────────
  if (!is.null(sample_levels))
    seurat_object@meta.data[[sample_col]] <- factor(seurat_object@meta.data[[sample_col]],
                                          levels = sample_levels)

  # ── Sample colors ──────────────────────────────────────────────────────────
  samp_lvls <- levels(factor(seurat_object@meta.data[[sample_col]]))
  if (is.null(sample_colors)) {
    sample_colors <- stats::setNames(Nour_pal("all")(length(samp_lvls)),
                                     samp_lvls)
  }

  # ── Compute QC metrics ──────────────────────────────────────────────────────
  message("Computing QC metrics (species: ", species, ")...")
  seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)

  seurat_object$mitoRatio <- Seurat::PercentageFeatureSet(seurat_object, pattern = mt_pattern,
                                                assay = assay) / 100
  seurat_object$riboRatio <- Seurat::PercentageFeatureSet(seurat_object, pattern = ribo_pattern,
                                                assay = assay) / 100

  CellInfo <- seurat_object@meta.data

  # ── Helper: clean barplot theme ─────────────────────────────────────────────
  .bar_theme <- function() {
    list(
      theme_NourMin(),
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 15),
        plot.margin  = ggplot2::unit(c(0.5, 0.5, 0.2, 0.5), "cm")
      ),
      ggplot2::labs(y = "Composition (Number of cells)", x = NULL)
    )
  }

  # ── 1. Doublet barplot ──────────────────────────────────────────────────────
  if ("DoubletStatus" %in% colnames(CellInfo)) {
    Maxp    <- max(as.matrix(table(CellInfo$DoubletStatus)))
    Heightp <- roundUpNice(Maxp)

    PDS <- CellInfo |>
      ggplot2::ggplot(ggplot2::aes(x = DoubletStatus,
                                   fill = .data[[group_col]])) +
      ggplot2::geom_bar(color = "black", width = 0.7) +
      ggplot2::scale_fill_manual(values = sample_colors) +
      .bar_theme() +
      ggplot2::scale_y_continuous(expand = c(0, 0),
                                  limits = c(0, Heightp),
                                  breaks = seq(0, Heightp, by = roundUpNice(Heightp / 5)))

    f_dbl <- pref("Barplot number of Doublets per sample.pdf")
    grDevices::pdf(f_dbl, width = 8, height = 5.5)
    print(PDS)
    grDevices::dev.off()
    .write_legend_sidecar(f_dbl, paste0(
      "Bar plot showing the number of cells classified as singlets or doublets ",
      "by SCP::db_scDblFinder in ", object_name, subset_name,
      ". Each bar is colored by sample of origin."
    ))
    message("  Saved: ", f_dbl)
  }

  # ── 2. miQC keep barplot ─────────────────────────────────────────────────────
  if ("miQC.keep" %in% colnames(CellInfo)) {
    Maxp    <- max(as.matrix(table(CellInfo[[sample_col]])))
    Heightp <- roundUpNice(Maxp)

    KeepP <- CellInfo |>
      ggplot2::ggplot(ggplot2::aes(x = miQC.keep,
                                   fill = .data[[group_col]])) +
      ggplot2::geom_bar(color = "black", width = 0.7) +
      ggplot2::scale_fill_manual(values = sample_colors) +
      .bar_theme() +
      ggplot2::scale_y_continuous(expand = c(0, 0),
                                  limits = c(0, Heightp),
                                  breaks = seq(0, Heightp, by = roundUpNice(Heightp / 5)))

    f_miqc <- pref("Barplot number of cells to keep per sample calc by miQC.pdf")
    grDevices::pdf(f_miqc, width = 8, height = 5.5)
    print(KeepP)
    grDevices::dev.off()
    .write_legend_sidecar(f_miqc, paste0(
      "Bar plot showing the number of cells flagged for retention or removal ",
      "by SeuratWrappers::RunMiQC in ", object_name, subset_name,
      ". Cells with a posterior probability of being low quality above the ",
      "cutoff threshold are marked 'discard'."
    ))
    message("  Saved: ", f_miqc)
  }

  # ── 3. Density + violin QC plots ─────────────────────────────────────────────
  message("  Building QC density and violin plots...")
  grp_aes <- ggplot2::aes(color = .data[[sample_col]],
                          fill  = .data[[sample_col]])

  x1 <- CellInfo |>
    ggplot2::ggplot(grp_aes) +
    ggplot2::aes(x = nCount_RNA) +
    ggplot2::geom_density(alpha = 0.2) + ggplot2::theme_classic() +
    ggplot2::ylab("Cell density") +
    ggplot2::geom_vline(xintercept = count_vlines) +
    ggplot2::scale_fill_manual(values = sample_colors) +
    ggplot2::scale_color_manual(values = sample_colors)

  x2 <- CellInfo |>
    ggplot2::ggplot(grp_aes) +
    ggplot2::aes(x = nFeature_RNA) +
    ggplot2::geom_density(alpha = 0.2) + ggplot2::theme_classic() +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept = feature_vlines) +
    ggplot2::scale_fill_manual(values = sample_colors) +
    ggplot2::scale_color_manual(values = sample_colors)

  x3 <- CellInfo |>
    ggplot2::ggplot(grp_aes) +
    ggplot2::aes(x = mitoRatio) +
    ggplot2::geom_density(alpha = 0.1) + ggplot2::theme_classic() +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept = mito_vline) +
    ggplot2::scale_fill_manual(values = sample_colors) +
    ggplot2::scale_color_manual(values = sample_colors)

  x4 <- CellInfo |>
    ggplot2::ggplot(grp_aes) +
    ggplot2::aes(x = log10GenesPerUMI) +
    ggplot2::geom_density(alpha = 0.2) + ggplot2::theme_classic() +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept = complexity_vline) +
    ggplot2::scale_fill_manual(values = sample_colors) +
    ggplot2::scale_color_manual(values = sample_colors)

  density_grid <- patchwork::wrap_plots(x1, x2, x3, x4,
                                        guides = "collect",
                                        ncol = 2, nrow = 2) &
    theme_NourMin()

  # Violin plots
  vln <- Seurat::VlnPlot(
    seurat_object,
    features = c("nFeature_RNA", "nCount_RNA", "mitoRatio", "log10GenesPerUMI"),
    ncol     = 2, pt.size = 0,
    cols     = sample_colors,
    group.by = sample_col
  ) & Seurat::NoLegend() & Seurat::RotatedAxis()

  f_qc <- pref("qualitycontrol.pdf")
  grDevices::pdf(f_qc, width = 10, height = 10)
  print(density_grid / patchwork::wrap_plots(vln))
  grDevices::dev.off()
  .write_legend_sidecar(f_qc, paste0(
    "Quality control metrics for all recovered ",
    if (species == "mouse") "cells" else "cells",
    " in ", object_name, subset_name, " prior to filtering. ",
    "Density distributions (top row) and violin plots (bottom row) show total UMI counts ",
    "(nCount_RNA), number of detected genes (nFeature_RNA), mitochondrial read fraction ",
    "(mitoRatio, pattern: ", mt_pattern, "), and transcriptional complexity ",
    "(log10GenesPerUMI), stratified by sample. ",
    "Vertical lines indicate the suggested filtering thresholds."
  ))
  message("  Saved: ", f_qc)

  # ── 4. Cell count barplot before filtration ──────────────────────────────────
  Max    <- max(as.matrix(table(CellInfo[[sample_col]])))
  Height <- roundUpNice(Max)

  P1 <- CellInfo |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[sample_col]],
                                 fill = .data[[sample_col]])) +
    ggplot2::scale_fill_manual(values = sample_colors) +
    ggplot2::geom_bar(color = "black", width = 0.7) +
    .bar_theme() +
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                limits = c(0, Height),
                                breaks = seq(0, Height, by = roundUpNice(Height / 5)))

  f_count <- pref("Barplot number of cells per Sample before filtration.pdf")
  grDevices::pdf(f_count, width = 10, height = 5.5)
  print(P1)
  grDevices::dev.off()
  .write_legend_sidecar(f_count, paste0(
    "Bar plot of total recovered cells per sample in ", object_name, subset_name,
    " prior to any quality-based filtering. ",
    "Each bar represents one sample (", sample_col, "), colored by sample identity. ",
    "Cell counts shown are pre-doublet-removal and pre-miQC filtration."
  ))
  message("  Saved: ", f_count)

  message("PlotQCMetrics complete. Four PDFs written to ", output_dir)
  invisible(seurat_object)
}
