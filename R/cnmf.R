# =============================================================================
# scSidekick - Consensus NMF (cNMF) gene expression programs
#
# RunCNMF()         - prepare + factorize + combine + k-selection diagnostic.
#                     Wraps dylkot/cNMF (Kotliar et al. 2019, eLife) through
#                     reticulate. This is the heavy step; it honors resume and
#                     caffeinate. Inspect the k-selection plot it writes, pick a
#                     K, then call GetCNMFPrograms().
# GetCNMFPrograms() - run consensus for a chosen K and density threshold, then
#                     pull usages (cells x K) and gene spectra (K x genes) back
#                     onto the Seurat object as a DimReduc + @misc$cnmf$results.
#
# cNMF is a Python package and is NOT an R dependency. Install it once into a
# Python / conda environment:
#     pip install cnmf
# and point RunCNMF(conda_env = "...") or RunCNMF(python = "...") at that
# environment. reticulate drives it from R.
#
# NOTE: the cNMF Python API names used below (prepare / factorize / combine /
# k_selection_plot / consensus / load_results) follow the documented dylkot/cNMF
# interface and should be validated against the installed cnmf version.
# =============================================================================

# ── Internal: configure reticulate and import the Python modules ──────────────
.nk_cnmf_setup <- function(conda_env = NULL, python = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required to run cNMF.\n",
         "  install.packages('reticulate')")

  # Was Python already bound this session? If so, use_condaenv/use_python below
  # are silently ignored (reticulate can only initialize one interpreter per
  # session), which is the usual cause of 'cnmf not found'.
  already_bound <- reticulate::py_available(initialize = FALSE)

  if (!is.null(conda_env)) {
    reticulate::use_condaenv(conda_env, required = TRUE)
  } else if (!is.null(python)) {
    reticulate::use_python(python, required = TRUE)
  }

  if (!reticulate::py_module_available("cnmf")) {
    active <- tryCatch(reticulate::py_config()$python, error = function(e) "unknown")
    hint <- if (already_bound && (!is.null(conda_env) || !is.null(python)))
      paste0(
        "  Python was ALREADY initialized in this R session (", active, "), so ",
        "the requested\n",
        "  environment could not be bound. RESTART R and, before any other ",
        "Python use, run:\n",
        "    Sys.setenv(RETICULATE_PYTHON = '",
        python %||% "/path/to/cnmf_env/bin/python", "')\n",
        "  (or reticulate::use_condaenv('", conda_env %||% "<name>",
        "', required = TRUE)) first, then call RunCNMF() again.")
    else
      paste0(
        "  Install it (e.g. `pip install cnmf`) and point RunCNMF() at the ",
        "environment with\n",
        "  conda_env = '<name>' or python = '/path/to/python'.\n",
        "  Active Python: ", active)
    stop("The Python package 'cnmf' is not available in the active ",
         "environment.\n", hint)
  }

  list(
    cnmf = reticulate::import("cnmf",    delay_load = FALSE),
    ad   = reticulate::import("anndata", delay_load = FALSE),
    pd   = reticulate::import("pandas",  delay_load = FALSE),
    np   = reticulate::import("numpy",   delay_load = FALSE)
  )
}

# ── Internal: render a cNMF diagnostic PNG into the active graphics device ─────
# (RStudio Plots pane). Silent no-op when not interactive or when no PNG reader
# is installed.
.nk_show_png <- function(path) {
  if (is.null(path) || !file.exists(path) || !interactive()) return(invisible(NULL))
  if (requireNamespace("png", quietly = TRUE) &&
      requireNamespace("grid", quietly = TRUE)) {
    img <- tryCatch(png::readPNG(path), error = function(e) NULL)
    if (!is.null(img)) {
      grid::grid.newpage()
      grid::grid.raster(img)
      return(invisible(NULL))
    }
  }
  if (requireNamespace("magick", quietly = TRUE)) {
    img <- tryCatch(magick::image_read(path), error = function(e) NULL)
    if (!is.null(img)) { print(img); return(invisible(NULL)) }
  }
  message("scSidekick: install 'png' or 'magick' to preview plots in RStudio. ",
          "Plot saved at: ", path)
  invisible(NULL)
}

# ── Internal: a space-free directory for cNMF to run in ───────────────────────
# cNMF shells out unquoted `cp` commands internally, which break on spaces in
# paths (common with Box / Dropbox / iCloud / "Cloud Storage" folders). When
# output_dir contains a space, run cNMF through a space-free symlink so all of
# cNMF's file operations see a clean path; the files physically live in the real
# output_dir. Returns the path cNMF should use.
.nk_cnmf_run_dir <- function(output_dir) {
  od <- normalizePath(output_dir, mustWork = FALSE)
  if (!grepl("[[:space:]]", od)) return(od)

  base <- tempdir()
  if (grepl("[[:space:]]", base)) base <- "/tmp"
  if (grepl("[[:space:]]", base))
    stop("cNMF cannot handle spaces in paths and no space-free scratch ",
         "directory was found. Please set output_dir to a path without spaces.")

  tag  <- gsub("[^A-Za-z0-9]+", "_", od)
  tag  <- substr(tag, max(1L, nchar(tag) - 80L), nchar(tag))
  link <- file.path(base, paste0("scSidekick_cnmf_", tag))

  # Sys.readlink returns "" for a non-symlink, NA for a non-existent path, and
  # the target for an existing symlink. The link is good only when it already
  # points at `od`; otherwise (re)create it.
  cur    <- suppressWarnings(Sys.readlink(link))
  cur_ok <- !is.na(cur) && nzchar(cur) &&
            normalizePath(cur, mustWork = FALSE) == od
  if (!isTRUE(cur_ok)) {
    if (file.exists(link) || (!is.na(cur) && nzchar(cur)))
      unlink(link, force = TRUE)
    if (!isTRUE(file.symlink(od, link)))
      stop("Could not create a space-free symlink for cNMF at '", link,
           "'. Please use an output_dir without spaces.")
  }
  message("scSidekick: output_dir contains a space; cNMF will run via a ",
          "space-free symlink:\n  ", link, "\n  -> ", od)
  link
}

# ── Internal: build a raw-counts AnnData (cells x genes) for cNMF ─────────────
# Optionally include meta.data columns in obs (needed as Harmony batch vars).
.nk_cnmf_build_adata <- function(seurat_object, mods, assay, genes_use = NULL,
                                 metadata_cols = NULL) {
  counts <- .get_layer_data(seurat_object, assay = assay, layer = "counts")
  if (is.null(counts) || nrow(counts) == 0)
    stop("Could not extract a 'counts' layer from assay '", assay, "'. ",
         "cNMF requires raw UMI counts, not normalized data.")

  if (!is.null(genes_use)) {
    keep <- intersect(genes_use, rownames(counts))
    if (length(keep) < 2)
      stop("Fewer than two of 'genes_use' are present in the object.")
    counts <- counts[keep, , drop = FALSE]
  }

  # cNMF expects cells x genes, raw counts. Transpose and hand a CSC sparse
  # matrix to reticulate (converted to scipy.sparse on the Python side).
  counts_t <- Matrix::t(methods::as(counts, "CsparseMatrix"))
  X_py     <- reticulate::r_to_py(counts_t)

  if (is.null(metadata_cols) || length(metadata_cols) == 0L) {
    obs <- mods$pd$DataFrame(index = reticulate::r_to_py(colnames(counts)))
  } else {
    miss <- setdiff(metadata_cols, colnames(seurat_object@meta.data))
    if (length(miss) > 0L)
      stop("metadata/batch column(s) not found: ", paste(miss, collapse = ", "))
    md <- seurat_object@meta.data[colnames(counts), metadata_cols, drop = FALSE]
    md[] <- lapply(md, as.character)   # categorical batch vars for Harmony
    obs <- reticulate::r_to_py(md)     # pandas index taken from row names (cells)
  }
  var <- mods$pd$DataFrame(index = reticulate::r_to_py(rownames(counts)))
  mods$ad$AnnData(X = X_py, obs = obs, var = var)
}


#' Run consensus NMF (cNMF) to discover gene expression programs
#'
#' @description
#' Wraps the dylkot/cNMF pipeline (Kotliar et al. 2019) through reticulate:
#' \code{prepare} (high-variance gene selection + normalization),
#' \code{factorize} (the heavy NMF replication step), \code{combine}, and the
#' \code{k_selection_plot} diagnostic. cNMF is inherently a two-stage,
#' human-in-the-loop method: run this function, inspect the k-selection plot it
#' writes, choose the number of programs \code{k}, then call
#' \code{\link{GetCNMFPrograms}} to extract the consensus result for that K.
#'
#' cNMF is a Python package, not an R dependency. Install it once
#' (\code{pip install cnmf}) and point this function at the environment with
#' \code{conda_env} or \code{python}.
#'
#' @param seurat_object A Seurat object containing raw counts.
#' @param assay Character or \code{NULL}. Assay to pull counts from. \code{NULL}
#'   uses the default assay.
#' @param name Character. Run name; cNMF writes all output under
#'   \code{output_dir/name/}. Default \code{"cnmf"}.
#' @param output_dir Character or \code{NULL}. Output directory. Walks up to the
#'   \code{output_dir} stored by \code{\link{PrepObject}} when \code{NULL}.
#' @param k_range Integer vector. Candidate numbers of programs (K) to factorize
#'   across. Default \code{5:15}.
#' @param n_iter Integer. NMF replicates per K. Default \code{100L}. Higher is
#'   more robust and much slower.
#' @param num_highvar_genes Integer. Number of overdispersed genes cNMF selects
#'   for factorization. Default \code{2000L}. Ignored when \code{genes_use} is
#'   supplied.
#' @param genes_use Character vector or \code{NULL}. Optional precomputed gene
#'   set to restrict the input to (bypasses cNMF's own high-variance selection).
#'   Ignored when \code{batch.by} is set.
#' @param batch.by Character vector or \code{NULL}. One or more \code{meta.data}
#'   columns to correct for batch effects before factorization. When supplied,
#'   cNMF's \code{Preprocess} Harmony-corrects the counts on these variables,
#'   learns program usages on the corrected data, and re-fits spectra on the
#'   uncorrected TP10K data (Kotliar batch-correction workflow). Requires the
#'   Python package \code{harmonypy} in the environment. \code{NULL} (default)
#'   runs without batch correction.
#' @param seed Integer. Random seed for reproducibility. Default \code{14L}.
#' @param n_workers Integer. Number of cNMF worker processes to run the
#'   \code{factorize} step in parallel. Default \code{1L} (single, in-process).
#'   When \code{> 1}, the env's \code{cnmf} command is launched once per worker
#'   and the (K x iteration) jobs are split across them; set this to the number
#'   of cores you can spare. Each worker writes a
#'   \code{<name>.factorize.worker_<i>.log} in \code{output_dir}.
#' @param conda_env Character or \code{NULL}. Name of the conda environment that
#'   has \code{cnmf} installed.
#' @param python Character or \code{NULL}. Path to a Python interpreter with
#'   \code{cnmf} installed (alternative to \code{conda_env}).
#' @param resume Logical. If \code{TRUE}, reuse work already on disk: a completed
#'   factorization (skip prepare/factorize/combine), and on the batch-correction
#'   path the existing Harmony-corrected files (skip the expensive correction
#'   step, e.g. to recover a run that crashed after Harmony). Default
#'   \code{FALSE}.
#' @param caffeinate Logical. Prevent the machine from sleeping during the run
#'   (macOS). Default \code{FALSE}.
#' @param show_plot Logical. Display the k-selection plot in the active graphics
#'   device (e.g. the RStudio Plots pane) when finished. Default \code{TRUE}.
#'   Requires the \code{png} (preferred) or \code{magick} package.
#'
#' @return The input Seurat object with run metadata stored in
#'   \code{seurat_object@misc$cnmf}. The k-selection plot is written to
#'   \code{output_dir/name/} with a \code{.legend} sidecar, and an
#'   \code{analysis_params.json} (parameters plus a draft methods paragraph) is
#'   written to \code{output_dir} for \code{\link{create_analysis_pptx}}.
#'
#' @details
#' \strong{Parallelism.} Factorization runtime scales with
#' \code{n_iter * length(k_range)} (e.g. the defaults are 100 x 11 = 1100 NMF
#' fits), not with object size. Set \code{n_workers > 1} to split those fits
#' across parallel cNMF worker processes on a multi-core machine. For
#' cluster-scale runs, cNMF's command-line workflow can distribute the same
#' workers across jobs; see the cNMF documentation.
#'
#' @seealso \code{\link{GetCNMFPrograms}}
#' @references Kotliar D, et al. (2019) Identifying gene expression programs of
#'   cell-type identity and cellular activity with single-cell RNA-Seq. eLife
#'   8:e43803.
#' @export
RunCNMF <- function(seurat_object,
                    assay             = NULL,
                    name              = "cnmf",
                    output_dir        = NULL,
                    k_range           = 5:15,
                    n_iter            = 100L,
                    num_highvar_genes = 2000L,
                    genes_use         = NULL,
                    batch.by          = NULL,
                    seed              = 14L,
                    n_workers         = 1L,
                    conda_env         = NULL,
                    python            = NULL,
                    resume            = FALSE,
                    caffeinate        = FALSE,
                    show_plot         = TRUE) {

  if (caffeinate) { .caff <- .nk_caffeinate(); on.exit(.nk_decaffeinate(.caff), add = TRUE) }

  assay <- assay %||% SeuratObject::DefaultAssay(seurat_object)

  output_dir <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL
  if (is.null(output_dir))
    stop("'output_dir' must be supplied or stored via PrepObject(output_dir = ...).")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  mods   <- .nk_cnmf_setup(conda_env, python)
  k_plot <- file.path(output_dir, name, paste0(name, ".k_selection.png"))

  if (resume && file.exists(k_plot)) {
    message("scSidekick: Resume - existing cNMF run found for '", name,
            "', skipping prepare/factorize.")
  } else {
    # cNMF runs against run_dir (= output_dir, or a space-free symlink to it).
    run_dir  <- .nk_cnmf_run_dir(output_dir)
    cnmf_obj <- mods$cnmf$cNMF(output_dir = run_dir, name = name)

    prep_args <- list(
      components = reticulate::r_to_py(as.integer(k_range)),
      n_iter     = as.integer(n_iter),
      seed       = as.integer(seed)
    )

    if (is.null(batch.by)) {
      # ── Standard path: write the counts AnnData and prepare from it ──────────
      message("scSidekick: Building cNMF input (", ncol(seurat_object), " cells)")
      adata      <- .nk_cnmf_build_adata(seurat_object, mods, assay, genes_use)
      input_h5ad <- file.path(run_dir, paste0(name, ".counts_input.h5ad"))
      adata$write(input_h5ad)
      prep_args$counts_fn <- input_h5ad
      if (is.null(genes_use)) {
        prep_args$num_highvar_genes <- as.integer(num_highvar_genes)
      } else {
        genes_file <- file.path(run_dir, paste0(name, ".genes_use.txt"))
        writeLines(intersect(genes_use, rownames(seurat_object)), genes_file)
        prep_args$genes_file <- genes_file
      }
      message("scSidekick: cNMF prepare - K = ", min(k_range), "-", max(k_range),
              ", ", n_iter, " iterations, ",
              if (is.null(genes_use)) paste0(num_highvar_genes, " high-var genes")
              else paste0(length(genes_use), " supplied genes"))
    } else {
      # ── Batch-corrected path: Harmony-correct counts via cNMF Preprocess, ────
      #    then prepare from the corrected-HVG + TP10K + HVG-list outputs. The
      #    Harmony step is expensive; with resume = TRUE we reuse the corrected
      #    files if they already exist (e.g. a run that crashed after Harmony).
      base      <- file.path(run_dir, paste0(name, ".batchcorrect"))
      corrected <- paste0(base, ".Corrected.HVG.Varnorm.h5ad")
      tp10k     <- paste0(base, ".TP10K.h5ad")
      hvgs_file <- paste0(base, ".Corrected.HVGs.txt")

      have_corrected <- file.exists(corrected) && file.exists(tp10k) &&
                        file.exists(hvgs_file)

      if (resume && have_corrected) {
        message("scSidekick: Resume - reusing existing Harmony-corrected files ",
                "(skipping batch correction).")
      } else {
        if (!is.null(genes_use))
          message("scSidekick: 'genes_use' is ignored when batch.by is set; ",
                  "Preprocess selects high-variance genes on the corrected data.")
        message("scSidekick: Building cNMF input (", ncol(seurat_object),
                " cells) with batch var(s): ", paste(batch.by, collapse = ", "))
        adata <- .nk_cnmf_build_adata(seurat_object, mods, assay, genes_use,
                                      metadata_cols = batch.by)
        message("scSidekick: Harmony batch correction on '",
                paste(batch.by, collapse = ", "), "' via cNMF Preprocess ",
                "(needs 'harmonypy' in the env) ...")
        pp   <- mods$cnmf$Preprocess(random_seed = as.integer(seed))
        hvar <- if (length(batch.by) == 1L) batch.by
                else reticulate::r_to_py(as.list(batch.by))
        pp$preprocess_for_cnmf(adata,
          harmony_vars     = hvar,
          n_top_rna_genes  = as.integer(num_highvar_genes),
          quantile_thresh  = 0.9999,
          makeplots        = TRUE,
          save_output_base = base)
      }
      prep_args$counts_fn  <- corrected
      prep_args$tpm_fn     <- tp10k
      prep_args$genes_file <- hvgs_file
      message("scSidekick: cNMF prepare on batch-corrected data - K = ",
              min(k_range), "-", max(k_range), ", ", n_iter, " iterations")
    }
    do.call(cnmf_obj$prepare, prep_args)

    n_workers <- max(1L, as.integer(n_workers))
    if (n_workers > 1L) {
      # Parallel factorize: launch n_workers cNMF worker processes via the
      # env's `cnmf` console script. cNMF splits the (K x iteration) jobs across
      # workers by --worker-index; `wait` blocks until all finish. We use the
      # CLI here (not reticulate) because reticulate is single-interpreter per
      # process and cannot run workers concurrently.
      cnmf_bin <- file.path(dirname(reticulate::py_exe()), "cnmf")
      if (!file.exists(cnmf_bin))
        stop("n_workers > 1 needs the 'cnmf' command at ", cnmf_bin,
             ".\n  Use n_workers = 1 for the in-process path.")
      message("scSidekick: cNMF factorize across ", n_workers,
              " parallel workers (logs: ", file.path(output_dir,
              paste0(name, ".factorize.worker_*.log")), ") ...")
      # Write the worker loop to a temp script and run that, so R's shell layer
      # does not re-split the inline `for ... do ... done` command.
      script <- tempfile(fileext = ".sh")
      writeLines(c(
        "#!/bin/bash",
        sprintf("for i in $(seq 0 %d); do", n_workers - 1L),
        sprintf(paste0("  %s factorize --output-dir %s --name %s ",
                       "--worker-index \"$i\" --total-workers %d ",
                       "> %s/%s.factorize.worker_\"$i\".log 2>&1 &"),
                shQuote(cnmf_bin), shQuote(run_dir), shQuote(name),
                n_workers, shQuote(run_dir), name),
        "done",
        "wait"), script)
      rc <- system2("bash", shQuote(script), wait = TRUE)
      unlink(script)
      if (!identical(as.integer(rc), 0L))
        warning("scSidekick: parallel factorize exit code ", rc,
                "; inspect the per-worker .log files in ", output_dir)
    } else {
      message("scSidekick: cNMF factorize (single worker; set n_workers > 1 ",
              "to parallelize) ...")
      cnmf_obj$factorize(worker_i = 0L, total_workers = 1L)
    }
    cnmf_obj$combine()
    cnmf_obj$k_selection_plot(close_fig = TRUE)
  }

  .write_legend_sidecar(k_plot, paste0(
    "cNMF k-selection diagnostic for run '", name, "': solution stability ",
    "(silhouette) and reconstruction error across K = ", min(k_range), "-",
    max(k_range), " (", n_iter, " iterations, ",
    if (is.null(genes_use)) paste0(num_highvar_genes, " high-variance genes")
    else paste0(length(genes_use), " supplied genes"),
    ", seed ", seed, "). Choose the K at the stability peak / error elbow and ",
    "pass it to GetCNMFPrograms(seurat_object, k = ...)."))

  seurat_object@misc$cnmf <- list(
    name              = name,
    output_dir        = output_dir,
    assay             = assay,
    k_range           = k_range,
    n_iter            = n_iter,
    num_highvar_genes = num_highvar_genes,
    seed              = seed,
    batch.by          = batch.by,
    k_selection_plot  = k_plot
  )

  # ── Auto-write method params so create_analysis_pptx() finds them ───────────
  .write_subdir_params(output_dir, list(
    date                = format(Sys.Date()),
    cnmf_name           = name,
    cnmf_assay          = assay,
    cnmf_k_range        = paste0(min(k_range), "-", max(k_range)),
    cnmf_n_iter         = n_iter,
    cnmf_highvar_genes  = if (is.null(genes_use)) num_highvar_genes
                          else paste0(length(genes_use), " supplied"),
    cnmf_seed           = seed,
    cnmf_n_workers      = n_workers,
    cnmf_batch_correct  = if (is.null(batch.by)) "none"
                          else paste(batch.by, collapse = ", "),
    methods_text        = paste0(
      "Gene expression programs were identified by consensus non-negative ",
      "matrix factorization (cNMF; Kotliar et al., 2019). Raw counts were ",
      "factorized over K = ", min(k_range), "-", max(k_range), " with ", n_iter,
      " replicates per K using ",
      if (is.null(genes_use)) paste0(num_highvar_genes, " high-variance genes")
      else paste0(length(genes_use), " supplied genes"),
      " (random seed ", seed, ")",
      if (!is.null(batch.by))
        paste0(", after Harmony batch correction on ",
               paste(batch.by, collapse = " and "),
               " (counts corrected, program spectra re-fit on the uncorrected ",
               "TP10K data)")
      else "",
      ". The number of programs was selected from the cNMF stability / error ",
      "diagnostic.")
  ))

  if (show_plot) .nk_show_png(k_plot)

  message("scSidekick: cNMF complete. Inspect '", k_plot,
          "', then run GetCNMFPrograms(obj, k = <chosen K>).")
  seurat_object
}


#' Extract cNMF consensus programs for a chosen K
#'
#' @description
#' Runs cNMF \code{consensus} for a selected number of programs \code{k} and a
#' \code{density_threshold} (chosen after inspecting the k-selection plot from
#' \code{\link{RunCNMF}}), then attaches the results to the Seurat object:
#' per-cell program usages become a dimensional reduction (and optionally
#' metadata columns), and gene spectra are stored in
#' \code{seurat_object@misc$cnmf$results}.
#'
#' @param seurat_object A Seurat object previously processed by
#'   \code{\link{RunCNMF}} (run metadata is read from \code{@misc$cnmf}).
#' @param k Integer. Number of programs to extract (the K you chose from the
#'   k-selection plot).
#' @param density_threshold Numeric. Local-density outlier filter for the
#'   consensus step. Lower is stricter; \code{2.0} keeps everything. Default
#'   \code{0.1}.
#' @param name Character or \code{NULL}. cNMF run name. Defaults to the name
#'   stored by \code{RunCNMF}.
#' @param output_dir Character or \code{NULL}. cNMF output directory. Defaults to
#'   the directory stored by \code{RunCNMF}.
#' @param reduction_name Character. Name of the \code{DimReduc} to create from
#'   the usage matrix. Default \code{"cnmf"}.
#' @param add_metadata Logical. Also add \code{cNMF_1 ... cNMF_k} usage columns
#'   to \code{meta.data}. Default \code{TRUE}.
#' @param top_n_genes Integer. Number of top genes per program to retain in the
#'   stored \code{top_genes} table. Default \code{50L}.
#' @param conda_env,python Character or \code{NULL}. Python environment with
#'   \code{cnmf} installed (as in \code{\link{RunCNMF}}).
#' @param caffeinate Logical. Prevent the machine from sleeping. Default
#'   \code{FALSE}.
#' @param show_plot Logical. Display the consensus clustergram in the active
#'   graphics device (e.g. the RStudio Plots pane). Default \code{TRUE}.
#'   Requires the \code{png} (preferred) or \code{magick} package.
#'
#' @return The Seurat object with a new \code{DimReduc} (\code{reduction_name}),
#'   optional usage metadata columns, and
#'   \code{seurat_object@misc$cnmf$results} holding \code{usage},
#'   \code{spectra_scores}, \code{spectra_tpm}, and \code{top_genes}. The
#'   consensus clustergram is written with a \code{.legend} sidecar, and the
#'   \code{analysis_params.json} in \code{output_dir} is updated with the
#'   consensus parameters and a finalized methods paragraph.
#'
#' @seealso \code{\link{RunCNMF}}
#' @export
GetCNMFPrograms <- function(seurat_object,
                            k,
                            density_threshold = 0.1,
                            name              = NULL,
                            output_dir        = NULL,
                            reduction_name    = "cnmf",
                            add_metadata      = TRUE,
                            top_n_genes       = 50L,
                            conda_env         = NULL,
                            python            = NULL,
                            caffeinate        = FALSE,
                            show_plot         = TRUE) {

  if (caffeinate) { .caff <- .nk_caffeinate(); on.exit(.nk_decaffeinate(.caff), add = TRUE) }

  meta <- seurat_object@misc$cnmf
  if (is.null(meta))
    stop("No cNMF run found on this object. Run RunCNMF() first.")
  name       <- name       %||% meta$name
  output_dir <- output_dir %||% meta$output_dir
  assay      <- meta$assay  %||% SeuratObject::DefaultAssay(seurat_object)

  mods     <- .nk_cnmf_setup(conda_env, python)
  run_dir  <- .nk_cnmf_run_dir(output_dir)   # space-free symlink if needed
  cnmf_obj <- mods$cnmf$cNMF(output_dir = run_dir, name = name)

  message("scSidekick: cNMF consensus - k = ", k,
          ", density_threshold = ", density_threshold)
  cnmf_obj$consensus(k                     = as.integer(k),
                     density_threshold     = density_threshold,
                     show_clustering       = TRUE,
                     close_clustergram_fig = TRUE)

  res <- cnmf_obj$load_results(K = as.integer(k),
                               density_threshold = density_threshold)

  # reticulate auto-converts the returned pandas DataFrames to R data.frames,
  # carrying the pandas index across as row names. Coerce to a matrix, falling
  # back to py_to_r only if a raw Python object comes through.
  .as_rmat <- function(x) {
    if (inherits(x, "python.builtin.object")) x <- reticulate::py_to_r(x)
    as.matrix(x)
  }
  usage_mat <- .as_rmat(res[[1]])   # cells x K, rownames = cell barcodes
  if (is.null(rownames(usage_mat)))
    stop("cNMF usage matrix has no cell names; cannot map back to the object.")
  colnames(usage_mat) <- paste0("cNMF_", seq_len(ncol(usage_mat)))

  # Build a full cells x K matrix (cells dropped as outliers get 0 usage).
  full <- matrix(0, nrow = ncol(seurat_object), ncol = ncol(usage_mat),
                 dimnames = list(colnames(seurat_object), colnames(usage_mat)))
  common <- intersect(rownames(full), rownames(usage_mat))
  full[common, ] <- usage_mat[common, , drop = FALSE]

  dr <- SeuratObject::CreateDimReducObject(
    embeddings = full, key = "cNMF_", assay = assay)
  seurat_object[[reduction_name]] <- dr

  if (add_metadata)
    for (j in seq_len(ncol(full)))
      seurat_object[[colnames(full)[j]]] <- full[, j]

  # Gene spectra (programs x genes) and top genes per program.
  spectra_scores <- .as_rmat(res[[2]])
  spectra_tpm    <- .as_rmat(res[[3]])
  top_genes      <- tryCatch({
    tg <- res[[4]]
    if (inherits(tg, "python.builtin.object")) tg <- reticulate::py_to_r(tg)
    utils::head(as.data.frame(tg), as.integer(top_n_genes))
  }, error = function(e) NULL)

  seurat_object@misc$cnmf$results <- list(
    k                 = k,
    density_threshold = density_threshold,
    usage             = full,
    spectra_scores    = spectra_scores,
    spectra_tpm       = spectra_tpm,
    top_genes         = top_genes
  )

  dt_str <- gsub("[.]", "_", format(density_threshold, trim = TRUE))
  clustergram <- file.path(output_dir, name,
    paste0(name, ".clustering.k_", k, ".dt_", dt_str, ".png"))
  .write_legend_sidecar(clustergram, paste0(
    "cNMF consensus clustergram for run '", name, "', k = ", k,
    " programs, density_threshold = ", density_threshold,
    ". Usages stored as reduction '", reduction_name, "' (cNMF_1..cNMF_", k,
    "); gene spectra and top genes per program in @misc$cnmf$results."))

  if (show_plot) .nk_show_png(clustergram)

  # ── Update method params: preserve the factorization fields written by
  #    RunCNMF, add the consensus step, and finalize the methods paragraph ─────
  meta_cnmf  <- seurat_object@misc$cnmf
  local_json <- file.path(output_dir, "analysis_params.json")
  prev <- if (requireNamespace("jsonlite", quietly = TRUE) &&
              file.exists(local_json))
    tryCatch(jsonlite::read_json(local_json, simplifyVector = TRUE),
             error = function(e) list()) else list()
  kr <- meta_cnmf$k_range
  full_methods <- paste0(
    "Gene expression programs were identified by consensus non-negative matrix ",
    "factorization (cNMF; Kotliar et al., 2019)",
    if (!is.null(meta_cnmf$batch.by))
      paste0(", with Harmony batch correction on ",
             paste(meta_cnmf$batch.by, collapse = " and "),
             " applied to the counts prior to factorization") else "",
    ". Raw counts were factorized",
    if (!is.null(kr)) paste0(" over K = ", min(kr), "-", max(kr)) else "",
    if (!is.null(meta_cnmf$n_iter))
      paste0(" with ", meta_cnmf$n_iter, " replicates per K") else "",
    if (!is.null(meta_cnmf$num_highvar_genes))
      paste0(" using ", meta_cnmf$num_highvar_genes, " high-variance genes") else "",
    ". A consensus solution with ", ncol(full), " programs (K = ", k,
    ") was selected, filtering outlier components at a local density threshold ",
    "of ", density_threshold, "; per-cell program usages were retained as a ",
    "dimensional reduction for downstream analysis.")
  .write_subdir_params(output_dir, utils::modifyList(prev, list(
    date                 = format(Sys.Date()),
    cnmf_k               = k,
    cnmf_n_programs      = ncol(full),
    cnmf_density_thresh  = density_threshold,
    cnmf_reduction       = reduction_name,
    methods_text         = full_methods
  )))

  message("scSidekick: Extracted ", ncol(full), " cNMF programs -> reduction '",
          reduction_name, "'", if (add_metadata) " and cNMF_* metadata" else "",
          ". Spectra in @misc$cnmf$results.")
  seurat_object
}


#' Top genes per cNMF program
#'
#' @description
#' Pulls the top-scoring genes for each gene expression program from the
#' spectra stored by \code{\link{GetCNMFPrograms}} in
#' \code{seurat_object@misc$cnmf$results}. Returns a tidy data frame (one row
#' per program-gene) ready to feed into \code{\link{GroupHeatmap}},
#' \code{\link{SplitDotPlot}}, or \code{\link{RunGSEA}} as a marker table.
#'
#' @param seurat_object A Seurat object processed by \code{\link{GetCNMFPrograms}}.
#' @param n Integer. Top genes per program to return. Default \code{20L}
#'   (capped at the number stored, set by \code{top_n_genes} in
#'   \code{GetCNMFPrograms}).
#' @param programs Integer/character vector or \code{NULL}. Programs to include,
#'   given as numbers (\code{c(1, 3)}) or labels (\code{c("cNMF_1", "cNMF_3")}).
#'   \code{NULL} (default) returns all.
#'
#' @return A data frame with columns \code{program} (\code{"cNMF_<k>"}),
#'   \code{gene}, and \code{rank}.
#'
#' @examples
#' \dontrun{
#' tg <- GetCNMFTopGenes(obj, n = 15)
#' # heatmap of program marker genes, blocked by program:
#' GroupHeatmap(obj, features = tg$gene, group.by = "CellType",
#'              feature_split = stats::setNames(tg$program, tg$gene))
#' }
#' @seealso \code{\link{GetCNMFPrograms}}, \code{\link{GroupHeatmap}}
#' @export
GetCNMFTopGenes <- function(seurat_object, n = 20L, programs = NULL) {
  res <- seurat_object@misc$cnmf$results
  if (is.null(res) || is.null(res$top_genes))
    stop("No cNMF top genes found. Run GetCNMFPrograms() first.")
  tg <- as.data.frame(res$top_genes, stringsAsFactors = FALSE)
  prog_cols <- colnames(tg)

  if (!is.null(programs)) {
    want <- gsub("^cNMF_", "", as.character(programs))
    keep <- prog_cols[prog_cols %in% want]
    if (length(keep) == 0L)
      stop("None of the requested programs found. Available: ",
           paste0("cNMF_", prog_cols, collapse = ", "))
    prog_cols <- keep
  }

  out <- do.call(rbind, lapply(prog_cols, function(p) {
    genes <- as.character(tg[[p]])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    genes <- utils::head(genes, n)
    if (length(genes) == 0L) return(NULL)
    data.frame(program = paste0("cNMF_", p),
               gene    = genes,
               rank    = seq_along(genes),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}
