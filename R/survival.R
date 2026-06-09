# =============================================================================
# survival.R  --  SurvPlot
#
# Kaplan-Meier survival analysis for bulk RNA-seq (or microarray) expression.
# Accepts either a pre-loaded expression matrix + metadata, or downloads and
# caches a TCGA/TARGET project via one of two backends:
#
#   recount3      (default) – pre-processed STAR counts served via Bioconductor;
#                             single HTTP download per project, very reliable.
#   TCGAbiolinks  (legacy)  – chunked downloads from the GDC API; prone to
#                             partial-download retry crashes on some networks.
#
# Two input modes
#   genes      – single-gene expression; one KM panel per gene
#   gene_lists – named list of gene vectors scored as mean expression (using
#                the same bin-control scorer as PlotPathwayButterfly); one KM
#                panel per signature
#
# Splitting strategies
#   "median"   – High/Low split at the within-group median (default)
#   "quartile" – four groups (Q1–Q4) via quantile()
#   "value"    – user-supplied fixed threshold (e.g. 0 for enrichment scores)
#   "optimal"  – auto-detect via survminer::surv_cutpoint()
#
# Faceting
#   facet_row / facet_col control how panels are arranged.  When both are
#   supplied the output is one PDF page per gene/signature, with panels
#   arranged in a (n_row_levels × n_col_levels) grid.  When neither is
#   supplied, all genes share a single row.
#
# Sidecar CSV
#   Written alongside every PDF: records the cutoff, group n, split method,
#   and log-rank p-value for every panel.  Serves as the audit trail for the
#   stratification strategy used in each figure.
# =============================================================================


# ── Internal: TCGA download / cache ──────────────────────────────────────────
#
# recount3 (default)
#   Downloads one pre-built RSE object per project — a single consolidated
#   genes × samples count matrix + colData clinical table.  No per-patient
#   files, no chunked downloads.  Much faster and more reliable than GDC.
#   Counts are converted to log2(TPM + 1) using gene lengths from rowRanges.
#   Clinical columns are prefixed "tcga.gdc_cases.*"; SurvPlot standardises
#   OS (days) and status (0/1) automatically and prints the column names so
#   the user knows what is available for faceting.
#
# TCGAbiolinks (fallback / microarray)
#   Downloads one file per patient from the GDC API — 391 files for TCGA-GBM.
#   Prone to partial-download crashes.  Used automatically for microarray
#   data_type since recount3 only holds RNA-seq.

#' @noRd
.surv_recount3 <- function(project, cache_dir) {
  for (pkg in c("recount3", "SummarizedExperiment"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(pkg, " is required.\n",
           "Install with: BiocManager::install('", pkg, "')")

  # recount3 uses bare cancer codes: "TCGA-GBM" → "GBM"
  is_target  <- grepl("^TARGET-", project, ignore.case = TRUE)
  proj_code  <- sub("^TCGA-|^TARGET-", "", project, ignore.case = TRUE)
  proj_home  <- if (is_target) "data_sources/target" else "data_sources/tcga"

  message("[SurvPlot] Fetching consolidated matrix for ", project,
          " via recount3 (one download, no per-patient files)...")

  rse <- recount3::create_rse_manual(
    project      = proj_code,
    project_home = proj_home,
    organism     = "human",
    annotation   = "gencode_v26",
    type         = "gene"
  )

  # ── counts → log2(TPM + 1) ───────────────────────────────────────────────
  counts <- SummarizedExperiment::assay(rse, "raw_counts")
  rlens  <- SummarizedExperiment::rowRanges(rse)$bp_length

  if (!is.null(rlens) && !all(is.na(rlens))) {
    rpk  <- counts / pmax(as.numeric(rlens) / 1000, 1)
    tpm  <- t(t(rpk) / (colSums(rpk, na.rm = TRUE) / 1e6))
    expr <- log2(tpm + 1)
  } else {
    # fall back to RPM when gene lengths are unavailable
    message("[SurvPlot] Gene lengths unavailable — using log2(RPM + 1) instead of TPM.")
    lib  <- colSums(counts)
    expr <- log2(t(t(counts) / (lib / 1e6)) + 1)
  }

  # gene symbols as rownames
  rdat <- SummarizedExperiment::rowData(rse)
  sym_col <- intersect(c("gene_name", "Symbol", "SYMBOL"), colnames(rdat))[1]
  if (!is.na(sym_col))
    rownames(expr) <- make.unique(as.character(rdat[[sym_col]]))

  # ── clinical metadata ─────────────────────────────────────────────────────
  clin <- as.data.frame(SummarizedExperiment::colData(rse))

  # recount3 TCGA columns are prefixed "tcga.gdc_cases."
  # Standardise OS and status from the most common column patterns
  death_col   <- grep("days_to_death$",          colnames(clin), value = TRUE)[1]
  followup_col <- grep("days_to_last_follow",     colnames(clin), value = TRUE)[1]
  vital_col   <- grep("vital_status$",            colnames(clin), value = TRUE)[1]
  gender_col  <- grep("demographic\\.gender$|\\bgender\\b", colnames(clin), value = TRUE)[1]

  if (!is.na(death_col) && !is.na(followup_col)) {
    d <- suppressWarnings(as.numeric(clin[[death_col]]))
    f <- suppressWarnings(as.numeric(clin[[followup_col]]))
    clin$OS <- ifelse(is.na(d), f, d)
  }
  if (!is.na(vital_col)) {
    clin$status <- as.integer(
      tolower(as.character(clin[[vital_col]])) %in% c("dead", "deceased")
    )
  }
  if (!is.na(gender_col) && !"gender" %in% colnames(clin))
    clin$gender <- clin[[gender_col]]

  # ── Enrich with richer clinical data (GDC comprehensive + paper subtypes) ──
  clin <- .surv_enrich_clinical(clin, project)

  # print available columns so users know what to pass to facet_row / facet_col
  message("[SurvPlot] Clinical columns available for facet_row / facet_col:\n  ",
          paste(sort(colnames(clin)), collapse = "\n  "))

  common <- intersect(colnames(expr), rownames(clin))
  list(expr_mat  = expr[, common, drop = FALSE],
       meta      = clin[common, , drop = FALSE],
       project   = project,
       data_type = "rnaseq",
       backend   = "recount3")
}

# ── Internal: enrich recount3 clinical data with GDC + paper subtypes ────────
#
# TCGAbiolinks::GDCquery_clinic() returns a comprehensive single-table clinical
# data frame (age, stage, race, treatment, etc.) from GDC via REST API — fast,
# no file downloads. TCGAquery_subtype() returns marker-paper molecular subtypes
# (IDH status, PAM50, CIMP, Pan-Glioma clusters, etc.).
# Both are merged into the recount3 colData by patient barcode.
#' @noRd
.surv_enrich_clinical <- function(clin, project) {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    message("[LoadTCGA] TCGAbiolinks not installed — metadata will use recount3 ",
            "basic clinical columns only. For richer metadata (IDH status, ",
            "molecular subtypes, etc.) install with: ",
            "BiocManager::install('TCGAbiolinks')")
    return(clin)
  }

  # Resolve patient barcodes (12-char TCGA-XX-XXXX format).
  # recount3 colData usually carries tcga.gdc_cases.submitter_id directly.
  pat_col <- grep("cases\\.submitter_id$", colnames(clin), value = TRUE)[1]
  if (!is.na(pat_col)) {
    patient_barcodes <- toupper(as.character(clin[[pat_col]]))
  } else {
    # Parse from rownames: dots → hyphens, take first 3 segments
    ids <- gsub("\\.", "-", rownames(clin))
    patient_barcodes <- toupper(
      sapply(strsplit(ids, "-"),
             function(x) paste(x[seq_len(min(3L, length(x)))], collapse = "-"))
    )
  }

  # ── GDC comprehensive clinical ────────────────────────────────────────────
  gdc_clin <- tryCatch(
    TCGAbiolinks::GDCquery_clinic(project = project, type = "clinical"),
    error = function(e) {
      message("[LoadTCGA] GDCquery_clinic failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(gdc_clin) && "submitter_id" %in% colnames(gdc_clin)) {
    gdc_clin$submitter_id <- toupper(gdc_clin$submitter_id)
    rownames(gdc_clin)    <- gdc_clin$submitter_id
    new_cols <- setdiff(colnames(gdc_clin), c(colnames(clin), "submitter_id"))
    if (length(new_cols) > 0L) {
      extra <- gdc_clin[patient_barcodes, new_cols, drop = FALSE]
      rownames(extra) <- rownames(clin)
      clin <- cbind(clin, extra)
      message("[LoadTCGA] Added ", length(new_cols),
              " GDC clinical columns (age, stage, race, treatment, ...).")
    }
  }

  # ── Molecular subtypes (paper_IDH.status, paper_Subtype, CIMP, etc.) ─────
  tumor_code <- sub("^TCGA-|^TARGET-", "", project, ignore.case = TRUE)
  subtypes <- tryCatch(
    suppressWarnings(TCGAbiolinks::TCGAquery_subtype(tumor = tumor_code)),
    error = function(e) NULL
  )
  if (!is.null(subtypes) && "patient" %in% colnames(subtypes)) {
    subtypes$patient <- toupper(subtypes$patient)
    rownames(subtypes) <- subtypes$patient
    new_cols <- setdiff(colnames(subtypes), c(colnames(clin), "patient"))
    if (length(new_cols) > 0L) {
      extra <- subtypes[patient_barcodes, new_cols, drop = FALSE]
      rownames(extra) <- rownames(clin)
      clin <- cbind(clin, extra)
      message("[LoadTCGA] Added ", length(new_cols),
              " molecular subtype columns (paper_IDH.status, paper_Subtype, etc.).")
    }
  }

  clin
}


#' @noRd
.surv_tcgabiolinks <- function(project, data_type, cache_dir) {
  for (pkg in c("TCGAbiolinks", "SummarizedExperiment"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(pkg, " is required.\n  Install: BiocManager::install('", pkg, "')")

  message("[SurvPlot] Downloading via TCGAbiolinks (one file per patient — ",
          "use tcga_backend = 'recount3' for a faster single-file download)...")

  if (data_type == "rnaseq") {
    query <- TCGAbiolinks::GDCquery(
      project       = project,
      data.category = "Transcriptome Profiling",
      data.type     = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
  } else {
    query <- TCGAbiolinks::GDCquery(
      project       = project,
      data.category = "Gene Expression",
      legacy        = TRUE
    )
  }

  TCGAbiolinks::GDCdownload(query, files.per.chunk = 5L)
  se <- TCGAbiolinks::GDCprepare(query)

  a_names <- SummarizedExperiment::assayNames(se)
  pick    <- intersect(c("tpm_unstrand", "fpkm_unstrand"), a_names)[1]
  pick    <- if (is.na(pick)) a_names[1] else pick
  expr    <- log2(SummarizedExperiment::assay(se, pick) + 1)

  rdat <- SummarizedExperiment::rowData(se)
  if ("gene_name" %in% colnames(rdat))
    rownames(expr) <- make.unique(as.character(rdat$gene_name))

  clin <- as.data.frame(SummarizedExperiment::colData(se))
  if ("days_to_death" %in% colnames(clin)) {
    clin$OS <- ifelse(
      is.na(clin$days_to_death),
      suppressWarnings(as.numeric(clin$days_to_last_follow_up)),
      suppressWarnings(as.numeric(clin$days_to_death))
    )
  }
  if ("vital_status" %in% colnames(clin))
    clin$status <- as.integer(
      tolower(clin$vital_status) %in% c("dead", "deceased")
    )

  common <- intersect(colnames(expr), rownames(clin))
  list(expr_mat  = expr[, common, drop = FALSE],
       meta      = clin[common, , drop = FALSE],
       project   = project,
       data_type = data_type,
       backend   = "TCGAbiolinks")
}

#' @noRd
.surv_tcga_get <- function(project, data_type, backend, cache_dir,
                           force_reload = FALSE) {
  if (is.null(cache_dir))
    cache_dir <- tools::R_user_dir("scSidekick", "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  cache_file <- file.path(
    cache_dir,
    paste0(gsub("[^A-Za-z0-9]", "_", project), "_", data_type,
           "_", backend, ".rds")
  )

  if (!force_reload && file.exists(cache_file)) {
    message("[SurvPlot] Loading cached data from:\n  ", cache_file,
            "\n  (pass force_reload = TRUE to re-download)")
    return(readRDS(cache_file))
  }

  if (force_reload && file.exists(cache_file)) {
    message("[LoadTCGA] force_reload = TRUE — ignoring cache, re-downloading...")
    file.remove(cache_file)
  }

  result <- if (backend == "recount3") {
    .surv_recount3(project, cache_dir)
  } else {
    .surv_tcgabiolinks(project, data_type, cache_dir)
  }

  saveRDS(result, cache_file)
  message("[SurvPlot] Cached to:\n  ", cache_file)
  result
}


# ── LoadTCGA ──────────────────────────────────────────────────────────────────

#' Download and cache a TCGA or TARGET expression dataset
#'
#' Downloads a consolidated expression matrix and clinical metadata for a
#' TCGA or TARGET project and caches the result locally as an \code{.rds}
#' file.  Subsequent calls with the same arguments load from the cache
#' instantly without re-downloading.
#'
#' The returned list can be passed directly to \code{\link{SurvMetaSummary}}
#' to explore available clinical columns before calling \code{\link{SurvPlot}}.
#'
#' @param project Character scalar.  Project ID, e.g. \code{"TCGA-GBM"},
#'   \code{"TCGA-LGG"}, \code{"TARGET-NBL"}.  See \code{\link{SurvPlot}} for
#'   the full list of available projects.
#' @param data_type \code{"rnaseq"} (default) or \code{"microarray"}.
#'   \code{"microarray"} forces \code{backend = "TCGAbiolinks"} automatically.
#' @param backend Download engine.  \code{"recount3"} (default) fetches a
#'   single pre-consolidated matrix — fast and reliable.
#'   \code{"TCGAbiolinks"} downloads one file per patient from the GDC API —
#'   use only when recount3 does not carry the project or for microarray data.
#' @param cache_dir Directory for the cached \code{.rds} file.  Defaults to
#'   the OS user cache directory for scSidekick.
#' @param force_reload Logical (default \code{FALSE}).  When \code{TRUE},
#'   deletes the existing \code{.rds} cache file and re-downloads the dataset
#'   from scratch.  Useful after adding TCGAbiolinks or when the upstream data
#'   have been updated.  The new download is cached as usual.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{expr_mat}}{Numeric matrix, genes × samples,
#'     \code{log2(TPM + 1)} values.}
#'   \item{\code{meta}}{Data frame of clinical metadata, one row per sample.}
#'   \item{\code{project}}{Project ID string.}
#'   \item{\code{data_type}}{Data modality used.}
#'   \item{\code{backend}}{Download engine used.}
#' }
#'
#' @examples
#' \dontrun{
#' gbm <- LoadTCGA("TCGA-GBM")
#' SurvMetaSummary(gbm$meta)
#' SurvPlot(expr_mat   = gbm$expr_mat,
#'          meta       = gbm$meta,
#'          genes      = c("YAP1", "S100A4"),
#'          time_col   = "OS",
#'          status_col = "status",
#'          facet_row  = "gender",
#'          facet_col  = "paper_IDH.status",
#'          output_file = "~/Desktop/GBM_survival/")
#' }
#'
#' @seealso \code{\link{SurvMetaSummary}}, \code{\link{SurvPlot}}
#' @export
LoadTCGA <- function(project,
                     data_type    = c("rnaseq", "microarray"),
                     backend      = c("recount3", "TCGAbiolinks"),
                     cache_dir    = NULL,
                     force_reload = FALSE) {
  data_type <- match.arg(data_type)
  backend   <- match.arg(backend)

  if (backend == "recount3" && data_type == "microarray") {
    message("[LoadTCGA] recount3 only supports RNA-seq. ",
            "Switching to backend = 'TCGAbiolinks' for microarray.")
    backend <- "TCGAbiolinks"
  }

  .surv_tcga_get(project, data_type, backend, cache_dir,
                 force_reload = force_reload)
}


# ── EnrichTCGAMeta ───────────────────────────────────────────────────────────

#' Enrich TCGA metadata with GDC clinical and molecular subtype columns
#'
#' Adds comprehensive clinical data and TCGA marker-paper molecular subtype
#' columns to an existing metadata data frame — without re-downloading the
#' expression matrix.  This is the fastest way to add \code{paper_IDH.status},
#' \code{paper_Subtype}, \code{paper_CIMP.status}, tumour stage, treatment
#' history, and other GDC variables to a dataset that was previously loaded
#' with \code{\link{LoadTCGA}}.
#'
#' Internally calls two fast GDC REST API endpoints via \pkg{TCGAbiolinks}
#' (no per-patient file downloads):
#' \itemize{
#'   \item \code{TCGAbiolinks::GDCquery_clinic()} — comprehensive clinical
#'     table (age at diagnosis, pathologic stage, race, vital status,
#'     treatment, etc.)
#'   \item \code{TCGAbiolinks::TCGAquery_subtype()} — TCGA marker-paper
#'     molecular subtypes (\code{paper_*} columns), available for most major
#'     cancer types.
#' }
#'
#' Requires \pkg{TCGAbiolinks}:
#' \code{BiocManager::install("TCGAbiolinks")}.
#'
#' @param tcga_data The list returned by \code{\link{LoadTCGA}}, or any named
#'   list with at least a \code{$meta} data frame and a \code{$project} string.
#' @param project Character scalar.  TCGA project ID (e.g. \code{"TCGA-GBM"}).
#'   Required only when \code{tcga_data} does not carry a \code{$project}
#'   element (e.g. when passing a raw metadata data frame).
#'
#' @return The input list with \code{$meta} replaced by the enriched data
#'   frame.  All original columns are preserved; new columns are appended.
#'   If only a metadata data frame was supplied, a data frame is returned.
#'
#' @examples
#' \dontrun{
#' # Typical workflow after first LoadTCGA call
#' gbm <- LoadTCGA("TCGA-GBM")
#' gbm <- EnrichTCGAMeta(gbm)
#' SurvMetaSummary(gbm$meta)   # now shows paper_IDH.status, stage, etc.
#'
#' # Or enrich a bare metadata data frame
#' meta_enriched <- EnrichTCGAMeta(gbm$meta, project = "TCGA-GBM")
#' }
#'
#' @seealso \code{\link{LoadTCGA}}, \code{\link{SurvMetaSummary}}
#' @export
EnrichTCGAMeta <- function(tcga_data, project = NULL) {
  # Accept either a LoadTCGA list or a bare meta data frame
  is_list <- is.list(tcga_data) && !is.data.frame(tcga_data)

  if (is_list) {
    meta    <- tcga_data$meta
    project <- project %||% tcga_data$project
  } else if (is.data.frame(tcga_data)) {
    meta <- tcga_data
  } else {
    stop("tcga_data must be the list returned by LoadTCGA(), or a metadata data frame.")
  }

  if (is.null(project) || !nzchar(project))
    stop("project must be supplied (e.g. 'TCGA-GBM') when tcga_data is a data frame.")

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
    stop("Package 'TCGAbiolinks' is required.\n",
         "Install with: BiocManager::install('TCGAbiolinks')")

  message("[EnrichTCGAMeta] Enriching metadata for ", project, "...")
  meta_enriched <- .surv_enrich_clinical(meta, project)

  n_new <- ncol(meta_enriched) - ncol(meta)
  message("[EnrichTCGAMeta] Added ", n_new, " new columns to metadata.")

  if (is_list) {
    tcga_data$meta <- meta_enriched
    return(invisible(tcga_data))
  } else {
    return(invisible(meta_enriched))
  }
}


# ── SurvMetaSummary ───────────────────────────────────────────────────────────

#' Summarise clinical metadata for survival analysis
#'
#' Prints a formatted overview of a clinical metadata data frame, identifying
#' which columns are suitable for \code{facet_row} / \code{facet_col} in
#' \code{\link{SurvPlot}}, which columns look like survival time or event
#' indicators, and flagging columns with high missingness.
#'
#' Designed to be called on the \code{$meta} element returned by
#' \code{\link{LoadTCGA}} before running \code{\link{SurvPlot}}.
#'
#' @param meta Data frame of sample-level clinical metadata.
#' @param max_levels Integer (default \code{10}).  Columns with more unique
#'   non-NA values than this are treated as continuous and shown separately.
#' @param na_threshold Numeric (default \code{0.5}).  Columns where more than
#'   this fraction of values are \code{NA} are flagged as high-missingness.
#' @param show_all Logical (default \code{FALSE}).  If \code{TRUE}, prints
#'   every column including those with > \code{max_levels} unique values.
#'   Useful for large data frames where you want to see numeric columns too.
#'
#' @return Invisibly returns a data frame with one row per column containing:
#'   \code{column}, \code{type}, \code{n_unique}, \code{pct_na},
#'   \code{category} (one of \code{"survival_time"}, \code{"survival_event"},
#'   \code{"good_facet"}, \code{"many_levels"}, \code{"high_na"},
#'   \code{"numeric"}), and \code{sample_values}.
#'
#' @examples
#' \dontrun{
#' gbm <- LoadTCGA("TCGA-GBM")
#' SurvMetaSummary(gbm$meta)
#'
#' # See all columns including numeric ones
#' SurvMetaSummary(gbm$meta, show_all = TRUE)
#' }
#'
#' @seealso \code{\link{LoadTCGA}}, \code{\link{SurvPlot}}
#' @export
SurvMetaSummary <- function(meta,
                             max_levels    = 10L,
                             na_threshold  = 0.5,
                             show_all      = FALSE) {
  n_samp <- nrow(meta)

  # ── classify every column ──────────────────────────────────────────────────
  col_info <- lapply(colnames(meta), function(cn) {
    x        <- meta[[cn]]
    n_na     <- sum(is.na(x))
    pct_na   <- n_na / n_samp
    is_num   <- is.numeric(x)
    n_unique <- length(unique(stats::na.omit(x)))

    # sample values: up to 5 unique non-NA values as a string
    vals <- unique(stats::na.omit(x))
    vals <- vals[seq_len(min(5L, length(vals)))]
    sample_vals <- paste(vals, collapse = ", ")
    if (length(unique(stats::na.omit(x))) > 5L)
      sample_vals <- paste0(sample_vals, ", ...")

    # category
    cn_l <- tolower(cn)
    cat <- dplyr::case_when(
      grepl("days_to|survival_time|os_time|^os$|\\btime\\b", cn_l) &
        is_num                                          ~ "survival_time",
      grepl("vital_status|^status$|event|deceased|dead", cn_l) &
        n_unique <= 3L                                  ~ "survival_event",
      pct_na > na_threshold                             ~ "high_na",
      !is_num && n_unique >= 2L && n_unique <= max_levels ~ "good_facet",
      is_num && n_unique >= 2L && n_unique <= max_levels  ~ "good_facet",
      is_num                                            ~ "numeric",
      n_unique > max_levels                             ~ "many_levels",
      TRUE                                              ~ "other"
    )

    data.frame(column       = cn,
               type         = class(x)[1],
               n_unique     = n_unique,
               pct_na       = round(pct_na * 100, 1),
               category     = cat,
               sample_values = sample_vals,
               stringsAsFactors = FALSE)
  })
  summary_df <- do.call(rbind, col_info)

  # ── print sections ─────────────────────────────────────────────────────────
  .section <- function(title, df, cols = c("column", "n_unique", "pct_na",
                                            "sample_values")) {
    if (nrow(df) == 0L) return(invisible(NULL))
    cat("\n", cli_rule(title), "\n", sep = "")
    print(df[, cols, drop = FALSE], row.names = FALSE, right = FALSE)
  }

  cli_rule <- function(x) {
    pad <- max(0L, 72L - nchar(x) - 4L)
    paste0("── ", x, " ", strrep("─", pad))
  }

  cat(sprintf("\n[SurvMetaSummary]  %d samples · %d columns\n",
              n_samp, ncol(meta)))

  .section(
    "Survival time columns  → use as time_col",
    summary_df[summary_df$category == "survival_time", ]
  )
  .section(
    "Survival event columns  → use as status_col",
    summary_df[summary_df$category == "survival_event", ]
  )
  .section(
    paste0("Good facet columns (2–", max_levels,
           " levels, <50% NA)  → use for facet_row / facet_col"),
    summary_df[summary_df$category == "good_facet", ]
  )
  .section(
    paste0("Many-level columns (>", max_levels,
           " unique values)  — consider binning before faceting"),
    summary_df[summary_df$category == "many_levels", ],
    cols = c("column", "type", "n_unique", "pct_na")
  )
  .section(
    "High-missingness columns (>50% NA)  — use with caution",
    summary_df[summary_df$category == "high_na", ],
    cols = c("column", "type", "n_unique", "pct_na")
  )
  if (show_all) {
    .section(
      "Numeric columns",
      summary_df[summary_df$category == "numeric", ],
      cols = c("column", "n_unique", "pct_na", "sample_values")
    )
  }

  cat("\n")
  invisible(summary_df)
}


# ── Internal: signature scorer ────────────────────────────────────────────────

#' Mean-expression score for each gene set (same logic as .bf_score)
#' Returns a features × samples matrix.
#' @noRd
.surv_score <- function(mat, gene_lists) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  s <- sapply(gene_lists, function(g) {
    g <- intersect(g, rownames(mat))
    if (length(g) == 0L) return(rep(NA_real_, ncol(mat)))
    colMeans(mat[g, , drop = FALSE])
  })
  if (is.null(dim(s)))
    s <- matrix(s, nrow = 1L, dimnames = list(names(gene_lists), colnames(mat)))
  t(s)  # features × samples
}


# ── Internal: classify samples into groups ────────────────────────────────────

#' @noRd
.surv_classify <- function(vals, method, split_value,
                            time_vec, status_vec, min_n) {
  n <- sum(!is.na(vals))
  empty <- list(labels = rep(NA_character_, length(vals)),
                cutoff = NA_real_, method = method, pval = NA_real_)
  if (n < min_n * 2L) return(empty)

  if (method == "median") {
    cutoff <- stats::median(vals, na.rm = TRUE)
    labels <- dplyr::if_else(is.na(vals), NA_character_,
                              dplyr::if_else(vals >= cutoff, "High", "Low"))
    return(list(labels = labels, cutoff = cutoff, method = "median",
                pval = NA_real_))
  }

  if (method == "value") {
    cutoff <- split_value
    labels <- dplyr::if_else(is.na(vals), NA_character_,
                              dplyr::if_else(vals >= cutoff, "High", "Low"))
    return(list(labels = labels, cutoff = cutoff, method = "value",
                pval = NA_real_))
  }

  if (method == "quartile") {
    qs     <- stats::quantile(vals, na.rm = TRUE)
    labels <- as.character(cut(vals, breaks = qs, include.lowest = TRUE,
                                labels = c("Q1", "Q2", "Q3", "Q4")))
    return(list(labels = labels, cutoff = qs, method = "quartile",
                pval = NA_real_))
  }

  if (method == "optimal") {
    if (!requireNamespace("survminer", quietly = TRUE))
      stop("survminer required for split_method = 'optimal'")
    tmp <- stats::na.omit(
      data.frame(val = vals, time = time_vec, status = status_vec)
    )
    if (nrow(tmp) < min_n * 2L) return(empty)
    result <- tryCatch({
      rc     <- survminer::surv_cutpoint(tmp, time = "time", event = "status",
                                          variables = "val", minprop = 0.2)
      cutoff <- rc$cutpoint$cutpoint
      list(
        labels = dplyr::if_else(is.na(vals), NA_character_,
                                dplyr::if_else(vals >= cutoff, "High", "Low")),
        cutoff = cutoff, method = "optimal", pval = NA_real_
      )
    }, error = function(e) {
      warning("[SurvPlot] surv_cutpoint failed — falling back to median: ",
              conditionMessage(e))
      .surv_classify(vals, "median", NULL, time_vec, status_vec, min_n)
    })
    return(result)
  }

  stop("Unknown split_method: '", method, "'")
}


# ── Internal: build one KM panel ─────────────────────────────────────────────

#' @noRd
.surv_one_panel <- function(dat, time_col, status_col, feat_col,
                             method, split_value, palette, show_pval,
                             show_n, title, min_n,
                             quartile_palette) {

  vals   <- dat[[feat_col]]
  t_vec  <- suppressWarnings(as.numeric(dat[[time_col]]))
  s_vec  <- suppressWarnings(as.numeric(dat[[status_col]]))

  cls <- .surv_classify(vals, method, split_value, t_vec, s_vec, min_n)
  if (all(is.na(cls$labels))) return(NULL)

  dat$SurvLevel <- factor(cls$labels)
  lv <- levels(dat$SurvLevel)
  if (length(lv) < 2L) return(NULL)

  # legend labels with n
  n_tab <- table(dat$SurvLevel)
  labs  <- if (show_n) {
    sapply(lv, function(x) paste0(x, " (n=", n_tab[x], ")"))
  } else lv

  pal <- if (method == "quartile") quartile_palette[seq_along(lv)] else palette[seq_along(lv)]

  # Use fixed column names .surv_time / .surv_status so the formula is always
  # a literal symbol — ggsurvplot reconstructs the formula from fit$call and
  # crashes when it finds a programmatically-built string instead of a symbol.
  dat$.surv_time   <- t_vec
  dat$.surv_status <- s_vec

  # compute log-rank p-value explicitly so we can store it in the sidecar
  sd   <- tryCatch(
    survival::survdiff(
      survival::Surv(.surv_time, .surv_status) ~ SurvLevel, data = dat
    ),
    error = function(e) NULL
  )
  pval <- if (!is.null(sd))
    stats::pchisq(sd$chisq, df = length(sd$n) - 1L, lower.tail = FALSE)
  else NA_real_

  fit <- survival::survfit(
    survival::Surv(.surv_time, .surv_status) ~ SurvLevel, data = dat
  )

  p <- survminer::ggsurvplot(
    fit,
    data          = dat,
    pval          = show_pval,
    palette       = pal,
    risk.table    = FALSE,
    legend.title  = feat_col,
    legend.labs   = unname(labs),
    font.main     = c(11, "bold"),
    font.x        = c(10, "bold"),
    font.y        = c(10, "bold"),
    font.tickslab = c(9, "plain"),
    ggtheme       = survminer::theme_survminer()
  )
  if (!is.null(title))
    p$plot <- p$plot + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold",
                                                         hjust = 0.5))

  list(ggsurvplot = p, cutoff = cls$cutoff, method = cls$method,
       n_total = nrow(dat), n_tab = n_tab, pval = pval)
}


# ── Internal: sidecar row ────────────────────────────────────────────────────

#' @noRd
.surv_sidecar_row <- function(feature, row_val, col_val, panel) {
  base <- data.frame(
    feature      = feature,
    facet_row    = if (is.na(row_val)) "all" else row_val,
    facet_col    = if (is.na(col_val)) "all" else col_val,
    split_method = if (is.null(panel)) NA_character_ else panel$method,
    n_total      = if (is.null(panel)) NA_integer_   else as.integer(panel$n_total),
    cutoff       = if (is.null(panel)) NA_character_ else {
      co <- panel$cutoff
      if (length(co) == 1L) as.character(round(co, 4L))
      else paste(round(co, 4L), collapse = "; ")
    },
    pval         = if (is.null(panel)) NA_real_ else panel$pval,
    stringsAsFactors = FALSE
  )
  if (!is.null(panel) && !is.null(panel$n_tab)) {
    for (lv in names(panel$n_tab))
      base[[paste0("n_", lv)]] <- as.integer(panel$n_tab[[lv]])
  }
  base
}


# ── Internal: multivariate Cox regression ────────────────────────────────────

#' @noRd
.surv_cox_fit <- function(dat, feat_col, time_col, status_col, covariates) {
  missing_cols <- setdiff(covariates, colnames(dat))
  if (length(missing_cols) > 0L)
    warning("[SurvPlot] Cox covariates not found in meta, skipped: ",
            paste(missing_cols, collapse = ", "), call. = FALSE)
  covariates <- intersect(covariates, colnames(dat))
  if (length(covariates) == 0L) return(NULL)

  t_vec <- suppressWarnings(as.numeric(dat[[time_col]]))
  s_vec <- suppressWarnings(as.numeric(dat[[status_col]]))

  all_vars <- c(feat_col, covariates)
  tmp <- dat[, all_vars, drop = FALSE]
  tmp$.surv_time   <- t_vec
  tmp$.surv_status <- s_vec
  tmp <- stats::na.omit(tmp)
  if (nrow(tmp) < 10L) return(NULL)

  # Ensure character covariates are factors so coxph can handle them
  for (v in covariates)
    if (is.character(tmp[[v]])) tmp[[v]] <- factor(tmp[[v]])

  rhs <- paste(
    sapply(all_vars, function(v) paste0("`", v, "`")),
    collapse = " + "
  )
  f <- stats::as.formula(
    paste0("survival::Surv(.surv_time, .surv_status) ~ ", rhs)
  )

  tryCatch({
    fit <- survival::coxph(f, data = tmp)
    list(fit = fit, data = tmp, n = nrow(tmp))
  }, error = function(e) {
    warning("[SurvPlot] coxph failed for '", feat_col, "': ",
            conditionMessage(e), call. = FALSE)
    NULL
  })
}


# ── Internal: build .legend text for one feature ─────────────────────────────

#' @noRd
.surv_legend_text <- function(feat, feature_type, feat_sidecar,
                               project, time_col, status_col,
                               split_method, split_per_group,
                               facet_row, facet_col,
                               has_cox, cox_covariates) {
  n_total <- if (nrow(feat_sidecar) > 0L && !is.na(feat_sidecar$n_total[1]))
               feat_sidecar$n_total[1] else NA_integer_

  feat_desc <- if (identical(feature_type, "gene"))
    paste0(feat, " gene expression (log2 TPM+1)")
  else
    paste0(feat, " gene-set signature (mean log2 expression)")

  split_desc <- switch(split_method,
    median   = paste0("median split",
                      if (split_per_group) " (per facet group)" else " (global)"),
    quartile = "quartile split (Q1–Q4)",
    value    = {
      cv <- feat_sidecar$cutoff[!is.na(feat_sidecar$cutoff)][1]
      paste0("fixed threshold at ",
             if (is.na(cv)) "user-specified value" else cv)
    },
    optimal  = "optimal log-rank cutpoint (survminer::surv_cutpoint)"
  )

  facet_desc <- if (!is.null(facet_row) && !is.null(facet_col))
    paste0("rows: ", facet_row, "; columns: ", facet_col)
  else if (!is.null(facet_row))
    paste0("rows: ", facet_row)
  else if (!is.null(facet_col))
    paste0("columns: ", facet_col)
  else
    "no faceting"

  valid <- feat_sidecar[!is.na(feat_sidecar$pval), ]
  pval_str <- if (nrow(valid) == 0L) {
    "log-rank p-values not available"
  } else {
    pv <- valid$pval
    if (length(pv) == 1L) sprintf("log-rank p = %.3g", pv)
    else sprintf("log-rank p range %.3g – %.3g", min(pv, na.rm = TRUE),
                 max(pv, na.rm = TRUE))
  }

  cohort <- if (!is.null(project)) project else "custom cohort"
  n_str  <- if (!is.na(n_total)) paste0(" (n = ", n_total, " samples)") else ""

  cox_str <- if (has_cox && length(cox_covariates) > 0L)
    paste0(" Multivariate Cox regression (covariates: ",
           paste(cox_covariates, collapse = ", "),
           ") shown on the following page.")
  else ""

  paste0(
    "Kaplan-Meier survival curves for ", feat_desc, ". ",
    "Cohort: ", cohort, n_str, ". ",
    "Survival endpoint: ", time_col, " (days); event indicator: ", status_col, ". ",
    "Stratification: ", split_desc, ". ",
    "Panel layout: ", facet_desc, ". ",
    pval_str, ".",
    cox_str,
    " Cutoffs and per-group sample sizes documented in the companion sidecar CSV."
  )
}


# ── Main exported function ────────────────────────────────────────────────────

#' Kaplan-Meier survival plots for bulk expression data
#'
#' @description
#' Generates publication-ready Kaplan-Meier survival plots for one or more
#' genes or gene-set signatures against bulk RNA-seq or microarray expression
#' data.  The function has two input modes:
#'
#' \enumerate{
#'   \item \strong{TCGA mode} — supply \code{tcga_project} and the function
#'     downloads, processes, and caches the expression matrix and clinical
#'     metadata automatically via \pkg{TCGAbiolinks}.  Subsequent calls with
#'     the same project load from the local cache instantly.
#'   \item \strong{Custom matrix mode} — supply your own \code{expr_mat}
#'     (genes × samples) and \code{meta} (samples × variables) data frames.
#'     This is the path for non-TCGA cohorts such as the Cavalli medulloblastoma
#'     dataset, PBTA, or any in-house bulk RNA-seq study.
#' }
#'
#' Samples are classified as High / Low (or Q1–Q4) using the chosen
#' \code{split_method}, and KM curves are drawn with log-rank p-values.
#' Results can be faceted across two clinical variables simultaneously (e.g.
#' sex as rows, tumour subtype as columns) producing one grid page per
#' gene/signature.  A sidecar CSV records every cutoff, group size, and
#' p-value for full reproducibility.
#'
#' @section TCGA projects:
#' Pass any of the following strings to \code{tcga_project}.  All IDs follow
#' the format \code{"TCGA-<CODE>"}.
#'
#' \strong{Brain / CNS}
#' \describe{
#'   \item{\code{"TCGA-GBM"}}{Glioblastoma multiforme (n ≈ 173)}
#'   \item{\code{"TCGA-LGG"}}{Brain Lower Grade Glioma (n ≈ 516)}
#' }
#'
#' \strong{Breast}
#' \describe{
#'   \item{\code{"TCGA-BRCA"}}{Breast invasive carcinoma (n ≈ 1 098)}
#' }
#'
#' \strong{Lung}
#' \describe{
#'   \item{\code{"TCGA-LUAD"}}{Lung adenocarcinoma (n ≈ 585)}
#'   \item{\code{"TCGA-LUSC"}}{Lung squamous cell carcinoma (n ≈ 504)}
#' }
#'
#' \strong{Colorectal}
#' \describe{
#'   \item{\code{"TCGA-COAD"}}{Colon adenocarcinoma (n ≈ 521)}
#'   \item{\code{"TCGA-READ"}}{Rectum adenocarcinoma (n ≈ 177)}
#' }
#'
#' \strong{Liver / Biliary}
#' \describe{
#'   \item{\code{"TCGA-LIHC"}}{Liver hepatocellular carcinoma (n ≈ 377)}
#'   \item{\code{"TCGA-CHOL"}}{Cholangiocarcinoma (n ≈ 51)}
#' }
#'
#' \strong{Kidney}
#' \describe{
#'   \item{\code{"TCGA-KIRC"}}{Kidney renal clear cell carcinoma (n ≈ 611)}
#'   \item{\code{"TCGA-KIRP"}}{Kidney renal papillary cell carcinoma (n ≈ 321)}
#'   \item{\code{"TCGA-KICH"}}{Kidney chromophobe (n ≈ 113)}
#' }
#'
#' \strong{Gynaecological}
#' \describe{
#'   \item{\code{"TCGA-OV"}}{Ovarian serous cystadenocarcinoma (n ≈ 608)}
#'   \item{\code{"TCGA-UCEC"}}{Uterine corpus endometrial carcinoma (n ≈ 587)}
#'   \item{\code{"TCGA-UCS"}}{Uterine carcinosarcoma (n ≈ 57)}
#'   \item{\code{"TCGA-CESC"}}{Cervical squamous cell carcinoma (n ≈ 309)}
#' }
#'
#' \strong{Skin}
#' \describe{
#'   \item{\code{"TCGA-SKCM"}}{Skin cutaneous melanoma (n ≈ 473)}
#'   \item{\code{"TCGA-UVM"}}{Uveal melanoma (n ≈ 80)}
#' }
#'
#' \strong{Gastrointestinal / Pancreatic}
#' \describe{
#'   \item{\code{"TCGA-STAD"}}{Stomach adenocarcinoma (n ≈ 443)}
#'   \item{\code{"TCGA-ESCA"}}{Esophageal carcinoma (n ≈ 185)}
#'   \item{\code{"TCGA-PAAD"}}{Pancreatic adenocarcinoma (n ≈ 185)}
#' }
#'
#' \strong{Thoracic}
#' \describe{
#'   \item{\code{"TCGA-MESO"}}{Mesothelioma (n ≈ 87)}
#'   \item{\code{"TCGA-THYM"}}{Thymoma (n ≈ 124)}
#' }
#'
#' \strong{Head and Neck / Thyroid}
#' \describe{
#'   \item{\code{"TCGA-HNSC"}}{Head and neck squamous cell carcinoma (n ≈ 528)}
#'   \item{\code{"TCGA-THCA"}}{Thyroid carcinoma (n ≈ 507)}
#' }
#'
#' \strong{Haematological}
#' \describe{
#'   \item{\code{"TCGA-LAML"}}{Acute myeloid leukaemia (n ≈ 200)}
#'   \item{\code{"TCGA-DLBC"}}{Diffuse large B-cell lymphoma (n ≈ 48)}
#' }
#'
#' \strong{Sarcoma / Mesenchymal}
#' \describe{
#'   \item{\code{"TCGA-SARC"}}{Sarcoma (n ≈ 265)}
#' }
#'
#' \strong{Urological}
#' \describe{
#'   \item{\code{"TCGA-BLCA"}}{Bladder urothelial carcinoma (n ≈ 433)}
#'   \item{\code{"TCGA-PRAD"}}{Prostate adenocarcinoma (n ≈ 551)}
#'   \item{\code{"TCGA-TGCT"}}{Testicular germ cell tumours (n ≈ 156)}
#' }
#'
#' \strong{Endocrine}
#' \describe{
#'   \item{\code{"TCGA-ACC"}}{Adrenocortical carcinoma (n ≈ 92)}
#'   \item{\code{"TCGA-PCPG"}}{Pheochromocytoma and paraganglioma (n ≈ 187)}
#' }
#'
#' \strong{Paediatric (TARGET programme)}
#' \describe{
#'   \item{\code{"TARGET-AML"}}{Acute myeloid leukaemia}
#'   \item{\code{"TARGET-NBL"}}{Neuroblastoma}
#'   \item{\code{"TARGET-OS"}}{Osteosarcoma}
#'   \item{\code{"TARGET-WT"}}{Wilms tumour}
#'   \item{\code{"TARGET-RT"}}{Rhabdoid tumour}
#'   \item{\code{"TARGET-CCSK"}}{Clear cell sarcoma of the kidney}
#' }
#'
#' @section Standard TCGA clinical columns:
#' After download, \code{SurvPlot} standardises \code{OS} (overall survival in
#' days, combining \code{days_to_death} and \code{days_to_last_follow_up}) and
#' \code{status} (1 = deceased, 0 = living/censored).  Common columns available
#' for \code{facet_row} / \code{facet_col} include:
#' \describe{
#'   \item{\code{gender}}{Reported sex (\code{"male"} / \code{"female"})}
#'   \item{\code{race}}{Self-reported race}
#'   \item{\code{age_at_index}}{Age at diagnosis (years)}
#'   \item{\code{tumor_stage} / \code{pathologic_stage}}{Pathological stage}
#'   \item{\code{paper_*}}{TCGA marker-paper columns, e.g.
#'     \code{paper_IDH.status}, \code{paper_Histology},
#'     \code{paper_CIMP.status}, \code{paper_Pan.Glioma.RNA.Expression.Cluster}}
#' }
#' Available columns vary by cancer type.  To inspect what is available after
#' the first download run: \code{colnames(result$sidecar)} or examine
#' \code{meta} from the returned list.
#'
#' @section Splitting strategies:
#' \describe{
#'   \item{\code{"median"} (default)}{Splits samples at the median expression /
#'     score value into \strong{High} and \strong{Low} groups.  When
#'     \code{split_per_group = TRUE} (recommended), the median is computed
#'     independently within each facet cell, so that a subgroup with
#'     uniformly low expression is still split 50/50 — matching the approach
#'     used in most published survival analyses.}
#'   \item{\code{"quartile"}}{Divides samples into four equally-sized groups
#'     (Q1 = lowest 25\%, Q4 = highest 25\%) using \code{quantile()}.  Produces
#'     four KM curves per panel coloured by \code{quartile_palette}.  Useful
#'     for detecting non-linear dose–response relationships.}
#'   \item{\code{"value"}}{Splits at a user-supplied fixed numeric threshold
#'     (\code{split_value}).  Applies the same cutoff across all groups and
#'     facets.  Set \code{split_value = 0} for gene-set enrichment scores
#'     computed by \code{gene_lists} mode (positive score = enriched).}
#'   \item{\code{"optimal"}}{Uses \code{survminer::surv_cutpoint()} to find
#'     the threshold that maximises the log-rank statistic.  Computationally
#'     intensive and prone to overfitting on small cohorts; not recommended
#'     for discovery.  Falls back to \code{"median"} if the optimisation
#'     fails.  Requires \pkg{survminer} \eqn{\geq} 0.4.9.}}
#'
#' @section Gene-set signature scoring:
#' When \code{gene_lists} is supplied instead of \code{genes}, the function
#' computes a per-sample score for each gene set as the mean log2 expression
#' across all genes in the set that are present in the matrix.  This is the
#' same mean-expression scorer used internally by \code{PlotPathwayButterfly};
#' no additional packages are required.  Missing genes are silently dropped
#' and a warning is issued if more than 50\% of a set is absent.
#'
#' For enrichment-style interpretation, use \code{split_method = "value"} with
#' \code{split_value = 0} (score \eqn{\geq} 0 = enriched) or
#' \code{split_method = "median"} (top half = enriched).
#'
#' @section Output layout:
#' \describe{
#'   \item{No faceting}{All genes / signatures appear on a single row in one
#'     page (or one PDF file per gene in directory mode).}
#'   \item{\code{facet_row} only}{Panels are stacked vertically, one row per
#'     level of the row variable.}
#'   \item{\code{facet_row} + \code{facet_col}}{One page per gene/signature.
#'     Panels are arranged in an
#'     (\emph{n\_row\_levels} × \emph{n\_col\_levels}) grid.  Cells with fewer
#'     than \code{2 × min_group_n} samples are left blank to preserve
#'     alignment.  PDF dimensions are auto-scaled to the grid size unless
#'     overridden with \code{pdf_width} / \code{pdf_height}.}
#' }
#'
#' @section Output file path:
#' \describe{
#'   \item{Single PDF}{Set \code{output_file} to a path ending in \code{.pdf}
#'     (e.g. \code{"~/results/survival.pdf"}).  All genes appear as separate
#'     pages in one file.}
#'   \item{Directory}{Set \code{output_file} to a folder path (existing or
#'     new; must \emph{not} end in \code{.pdf}).  One PDF per gene/signature
#'     is written, named \code{<gene>_survival.pdf}.  A single
#'     \code{survival_sidecar.csv} covering all genes is written to the same
#'     folder.}
#' }
#'
#' @section Sidecar CSV:
#' Written alongside every PDF output, the sidecar records one row per KM
#' panel with the following columns:
#' \describe{
#'   \item{\code{feature}}{Gene name or signature name}
#'   \item{\code{facet_row}, \code{facet_col}}{Level values (\code{"all"} when
#'     no faceting)}
#'   \item{\code{split_method}}{Method actually used (may differ from request
#'     if fallback occurred)}
#'   \item{\code{cutoff}}{Threshold applied; for \code{"quartile"} this is a
#'     semicolon-separated list of the five quantile boundaries}
#'   \item{\code{n_total}}{Total samples in the subset}
#'   \item{\code{n_High}, \code{n_Low} (or \code{n_Q1}–\code{n_Q4})}{Samples
#'     per group}
#'   \item{\code{pval}}{Log-rank p-value}
#' }
#'
#' @param expr_mat Numeric matrix, genes × samples.  Rownames must be gene
#'   symbols; colnames must be sample identifiers matching \code{rownames(meta)}.
#'   Values should be on a log scale (e.g. log2 TPM + 1 or log2 FPKM + 1).
#'   Supply either this + \code{meta}, or \code{tcga_project}.
#' @param meta Data frame of sample-level clinical metadata.  Rownames must
#'   match \code{colnames(expr_mat)}.  Must contain at minimum the columns
#'   named by \code{time_col} and \code{status_col}.
#' @param tcga_project Character scalar.  TCGA (or TARGET) project identifier.
#'   See the \strong{TCGA projects} section above for the full list of
#'   available IDs.  Requires \pkg{TCGAbiolinks} to be installed:
#'   \code{BiocManager::install("TCGAbiolinks")}.  Data are downloaded once
#'   and cached as an \code{.rds} file; subsequent calls load from the cache.
#' @param tcga_data_type Data modality to download from GDC.  One of:
#'   \describe{
#'     \item{\code{"rnaseq"} (default)}{STAR-aligned RNA-seq counts
#'       (\code{workflow.type = "STAR - Counts"}).  TPM values are extracted
#'       and log2-transformed (\code{log2(TPM + 1)}).  Falls back to FPKM if
#'       TPM is unavailable.}
#'     \item{\code{"microarray"}}{Legacy microarray data via the GDC legacy
#'       archive.  Use this for older TCGA accessions that pre-date RNA-seq,
#'       or for external cohorts processed on Affymetrix / Illumina arrays
#'       (e.g. the Cavalli medulloblastoma dataset).}
#'   }
#' @param tcga_backend Download engine for TCGA / TARGET data.  One of:
#'   \describe{
#'     \item{\code{"recount3"} (default)}{Downloads one pre-consolidated
#'       genes × samples count matrix per project — a single file, no
#'       per-patient downloads.  Counts are converted to
#'       \code{log2(TPM + 1)}.  Requires \pkg{recount3}:
#'       \code{BiocManager::install("recount3")}.  Automatically switches
#'       to \code{"TCGAbiolinks"} when \code{tcga_data_type = "microarray"}.}
#'     \item{\code{"TCGAbiolinks"}}{Downloads one file per patient from the
#'       GDC API (e.g. 391 files for TCGA-GBM).  Much slower and prone to
#'       partial-download crashes.  Use only when recount3 does not carry
#'       the project, or for microarray data.}
#'   }
#' @param tcga_cache_dir Path to a directory where downloaded data are cached
#'   as \code{.rds} files.  Defaults to the OS user cache directory for
#'   scSidekick (\code{tools::R_user_dir("scSidekick", "cache")}).  Set a
#'   shared network path to allow multiple users to reuse the same cache.
#'   Delete the \code{.rds} file to force a fresh download.
#' @param genes Character vector of gene symbols to analyse in single-gene
#'   mode.  Each gene produces one set of KM panels.  Genes absent from
#'   \code{expr_mat} trigger a warning and are skipped.
#' @param gene_lists Named list of character vectors for signature-score mode.
#'   Each element is a gene set; the name becomes the panel title.  Scores
#'   are the mean log2 expression of genes present in the matrix.  Example:
#'   \code{list(YAP_up = c("PRC1","HPCA"), YAP_dn = c("ZIC1","CHD7"))}.
#' @param time_col Name of the survival time column in \code{meta} (or in the
#'   standardised TCGA metadata).  For TCGA downloads this is \code{"OS"}
#'   (overall survival in days).  For disease-specific survival use the
#'   column produced by your preprocessing pipeline.
#' @param status_col Name of the event indicator column.  Must be numeric:
#'   \code{1} = event occurred (death / recurrence), \code{0} = censored.
#'   For TCGA downloads this is \code{"status"} (standardised from
#'   \code{vital_status} by \code{SurvPlot} automatically).
#' @param split_method Strategy for classifying samples into expression groups.
#'   One of \code{"median"} (default), \code{"quartile"}, \code{"value"}, or
#'   \code{"optimal"}.  See the \strong{Splitting strategies} section for
#'   full details of each option.
#' @param split_value Numeric cutoff used when \code{split_method = "value"}.
#'   Samples with expression / score \eqn{\geq} \code{split_value} are
#'   labelled \strong{High}; all others are \strong{Low}.  Commonly set to
#'   \code{0} for gene-set enrichment scores.
#' @param split_per_group Logical (default \code{TRUE}).  When \code{TRUE},
#'   the threshold (median / quantiles) is computed independently within
#'   each facet group.  This ensures a 50/50 split within every subtype or
#'   sex strata, preventing a globally-low subgroup from being assigned
#'   entirely to "Low".  Set to \code{FALSE} to use a single global threshold
#'   across all samples (equivalent to not splitting per group).
#' @param facet_row Name of a column in \code{meta} whose levels define the
#'   \emph{rows} of the output grid.  Common choices: \code{"gender"},
#'   \code{"race"}, \code{"age_group"}.  Set to \code{NULL} (default) for
#'   no row faceting.
#' @param facet_col Name of a column in \code{meta} whose levels define the
#'   \emph{columns} of the output grid.  Common choices: \code{"Subgroup"},
#'   \code{"paper_IDH.status"}, \code{"tumor_stage"}.  Set to \code{NULL}
#'   (default) for no column faceting.
#' @param output_file Output path.  Two modes depending on the path:
#'   \describe{
#'     \item{Ends in \code{.pdf}}{All genes are written as pages into a single
#'       multi-page PDF.}
#'     \item{Does not end in \code{.pdf}}{Treated as a directory.  One PDF per
#'       gene/signature is written (\code{<gene>_survival.pdf}), plus one
#'       shared \code{survival_sidecar.csv}.  The directory is created if it
#'       does not exist.}
#'   }
#'   Set to \code{NULL} (default) to suppress file output and return plots
#'   invisibly for further manipulation.
#' @param pdf_width,pdf_height Numeric.  PDF dimensions in inches.  When
#'   \code{NULL} (default) dimensions are auto-calculated from the number of
#'   genes and facet levels: approximately 4.5 in per column and 4.5 in per
#'   row.
#' @param palette Character vector of length \eqn{\geq 2}.  Colours for the
#'   \strong{High} and \strong{Low} groups respectively.  Default:
#'   \code{c("Red", "Blue")}.  Any R colour name or hex code is accepted.
#' @param quartile_palette Character vector of length \eqn{\geq 4}.  Colours
#'   for Q1 through Q4 when \code{split_method = "quartile"}.  Default is a
#'   blue-to-red diverging scale:
#'   \code{c("#2166AC", "#92C5DE", "#F4A582", "#D6604D")}.
#' @param show_pval Logical (default \code{TRUE}).  Display the log-rank
#'   p-value on each KM panel.
#' @param show_n Logical (default \code{TRUE}).  Append the group sample size
#'   to each legend label, e.g. \emph{"High (n = 47)"}.
#' @param write_sidecar Logical (default \code{TRUE}).  Write a companion CSV
#'   file recording the cutoff, group sizes, split method, and log-rank
#'   p-value for every panel.  Ignored when \code{output_file = NULL}.  See
#'   the \strong{Sidecar CSV} section for column definitions.
#' @param min_group_n Integer (default \code{5}).  Minimum number of samples
#'   required in each group after splitting.  Panels where any group falls
#'   below this threshold are skipped and recorded as \code{NA} in the sidecar.
#' @param multivariate_covariates Character vector of \code{meta} column names
#'   to include as covariates in a multivariate Cox proportional hazards model.
#'   When supplied, a forest plot (\code{survminer::ggforest}) is printed as an
#'   additional page after the KM panels for each gene / signature.  The Cox
#'   model uses continuous expression values (not High/Low labels) as the
#'   feature term, which is more powerful than treating it as a binary
#'   predictor.  All samples with complete data across the expression feature
#'   and covariates are used (no faceting for Cox).  Example:
#'   \code{c("gender", "age_at_index", "paper_IDH.status")}.  Set to
#'   \code{NULL} (default) to skip multivariate analysis.
#'
#' @return Invisibly returns a named list with two elements:
#' \describe{
#'   \item{\code{plots}}{Named list, one entry per gene/signature.  Each entry
#'     is itself a list of \code{ggsurvplot} objects (one per facet
#'     combination), or \code{NULL} for skipped panels.}
#'   \item{\code{sidecar}}{Data frame with one row per KM panel documenting
#'     the feature name, facet levels, cutoff, group sizes, and log-rank
#'     p-value.}
#' }
#'
#' @examples
#' \dontrun{
#' # ── Example 1: TCGA-GBM, two genes, sex × IDH status grid ────────────────
#' SurvPlot(
#'   tcga_project = "TCGA-GBM",
#'   genes        = c("YAP1", "S100A4"),
#'   time_col     = "OS",
#'   status_col   = "status",
#'   facet_row    = "gender",
#'   facet_col    = "paper_IDH.status",
#'   output_file  = "~/Desktop/GBM_survival/"   # directory → one PDF per gene
#' )
#'
#' # ── Example 2: TCGA-GBM, quartile split, all in one PDF ──────────────────
#' SurvPlot(
#'   tcga_project = "TCGA-GBM",
#'   genes        = c("YAP1", "CD276", "CD47"),
#'   split_method = "quartile",
#'   output_file  = "~/Desktop/GBM_quartiles.pdf"
#' )
#'
#' # ── Example 3: Gene-set signature scores, value split at 0 ───────────────
#' SurvPlot(
#'   tcga_project = "TCGA-LGG",
#'   gene_lists   = list(
#'     YAP_Activated = c("PRC1", "HPCA", "SRSF7", "TIMM17A", "PIM1"),
#'     YAP_Repressed = c("ZIC1", "ZIC2", "CHD7", "KIF5C", "DPYSL2")
#'   ),
#'   split_method  = "value",
#'   split_value   = 0,       # >= 0 = enriched
#'   facet_row     = "gender",
#'   facet_col     = "paper_Pan.Glioma.RNA.Expression.Cluster",
#'   output_file   = "~/Desktop/LGG_signatures/"
#' )
#'
#' # ── Example 4: Custom matrix (e.g. Cavalli medulloblastoma) ──────────────
#' Cav.EXP  <- t(read.table("Cavalli_expression.txt", header = TRUE,
#'                            row.names = 1))
#' Cav.meta <- read.table("Cavalli_pheno.txt", header = TRUE, row.names = 1)
#' # Cavalli metadata uses "survival" and "status" already
#' res <- SurvPlot(
#'   expr_mat     = Cav.EXP,
#'   meta         = Cav.meta,
#'   genes        = c("YAP1", "CD276", "CD47", "CCL5"),
#'   time_col     = "survival",
#'   status_col   = "status",
#'   facet_row    = "Gender",
#'   facet_col    = "Subgroup",
#'   split_method = "median",
#'   output_file  = "~/Desktop/Cavalli_survival/"
#' )
#' # Inspect the sidecar
#' head(res$sidecar)
#'
#' # ── Example 5: Paediatric cancer via TARGET ───────────────────────────────
#' SurvPlot(
#'   tcga_project  = "TARGET-NBL",
#'   genes         = c("MYCN", "ALK"),
#'   time_col      = "OS",
#'   status_col    = "status",
#'   split_method  = "optimal",
#'   output_file   = "~/Desktop/NBL_survival.pdf"
#' )
#'
#' # ── Example 6: No faceting, return plots for manual arrangement ───────────
#' res <- SurvPlot(
#'   tcga_project = "TCGA-SKCM",
#'   genes        = "S100A4",
#'   split_method = "median"
#' )
#' # Access the ggsurvplot object for the first panel
#' print(res$plots[["S100A4"]][[1]])
#' }
#'
#' @seealso
#' \code{\link{PlotPathwayButterfly}} for single-cell pathway state scoring
#' using the same mean-expression scorer.
#'
#' \href{https://bioconductor.org/packages/TCGAbiolinks/}{TCGAbiolinks} for
#' the GDC download engine used in TCGA mode.
#'
#' \href{https://rpkgs.datanovia.com/survminer/}{survminer} for the underlying
#' \code{ggsurvplot} plotting engine.
#'
#' @export
SurvPlot <- function(
  # ── Input ────────────────────────────────────────────────────────────────
  expr_mat         = NULL,
  meta             = NULL,
  tcga_project     = NULL,
  tcga_data_type   = c("rnaseq", "microarray"),
  tcga_backend     = c("recount3", "TCGAbiolinks"),
  tcga_cache_dir   = NULL,

  # ── Features ─────────────────────────────────────────────────────────────
  genes            = NULL,
  gene_lists       = NULL,

  # ── Survival columns ─────────────────────────────────────────────────────
  time_col         = "OS",
  status_col       = "status",

  # ── Splitting ────────────────────────────────────────────────────────────
  split_method     = c("median", "quartile", "value", "optimal"),
  split_value      = NULL,
  split_per_group  = TRUE,

  # ── Faceting ─────────────────────────────────────────────────────────────
  facet_row        = NULL,
  facet_col        = NULL,

  # ── Output ───────────────────────────────────────────────────────────────
  output_file      = NULL,
  pdf_width        = NULL,
  pdf_height       = NULL,
  palette          = c("Red", "Blue"),
  quartile_palette = c("#2166AC", "#92C5DE", "#F4A582", "#D6604D"),
  show_pval        = TRUE,
  show_n           = TRUE,
  write_sidecar            = TRUE,
  min_group_n              = 5L,
  multivariate_covariates  = NULL
) {

  # ── argument validation ───────────────────────────────────────────────────
  tcga_data_type <- match.arg(tcga_data_type)
  tcga_backend   <- match.arg(tcga_backend)
  split_method   <- match.arg(split_method)

  if (tcga_backend == "recount3" && tcga_data_type == "microarray") {
    message("[SurvPlot] recount3 only supports RNA-seq. ",
            "Switching to tcga_backend = 'TCGAbiolinks' for microarray.")
    tcga_backend <- "TCGAbiolinks"
  }

  if (is.null(tcga_project) && (is.null(expr_mat) || is.null(meta)))
    stop("Supply either tcga_project, or both expr_mat and meta.")
  if (!is.null(tcga_project) && (!is.null(expr_mat) || !is.null(meta)))
    warning("tcga_project supplied — ignoring expr_mat / meta.")
  if (is.null(genes) && is.null(gene_lists))
    stop("Supply genes (character vector) or gene_lists (named list).")
  if (!is.null(genes) && !is.null(gene_lists))
    stop("Supply genes or gene_lists, not both.")
  if (split_method == "value" && is.null(split_value))
    stop("split_value must be provided when split_method = 'value'.")

  for (pkg in c("survival", "survminer"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. Install with install.packages('", pkg, "').")

  # ── load / validate data ──────────────────────────────────────────────────
  if (!is.null(tcga_project)) {
    tcga <- .surv_tcga_get(tcga_project, tcga_data_type,
                           tcga_backend, tcga_cache_dir)
    expr_mat <- tcga$expr_mat
    meta     <- tcga$meta
  }

  common <- intersect(colnames(expr_mat), rownames(meta))
  if (length(common) == 0L)
    stop("No overlapping sample IDs between expr_mat columns and meta rownames.")
  expr_mat <- expr_mat[, common, drop = FALSE]
  meta     <- meta[common, , drop = FALSE]

  for (col in c(time_col, status_col))
    if (!col %in% colnames(meta))
      stop("Column '", col, "' not found in meta.")
  for (col in c(facet_row, facet_col))
    if (!is.null(col) && !col %in% colnames(meta))
      stop("Facet column '", col, "' not found in meta.")

  # ── build feature matrix ──────────────────────────────────────────────────
  if (!is.null(genes)) {
    missing <- setdiff(genes, rownames(expr_mat))
    if (length(missing))
      warning("[SurvPlot] Genes not found and will be skipped: ",
              paste(missing, collapse = ", "))
    features    <- intersect(genes, rownames(expr_mat))
    feature_mat <- expr_mat[features, , drop = FALSE]
    feature_type <- "gene"
  } else {
    feature_mat  <- .surv_score(expr_mat, gene_lists)
    features     <- rownames(feature_mat)
    feature_type <- "signature"
  }

  if (length(features) == 0L)
    stop("No valid features found in the expression matrix.")

  # safe column names (gene symbols may contain hyphens etc.)
  safe_feat  <- make.names(features)
  feat_map   <- stats::setNames(safe_feat, features)   # original → safe

  # merge into one data frame
  plot_dat <- cbind(
    meta,
    as.data.frame(t(feature_mat),
                  col.names = safe_feat,
                  check.names = FALSE)
  )
  colnames(plot_dat)[(ncol(meta) + 1L):ncol(plot_dat)] <- safe_feat

  # ── build facet combination table ────────────────────────────────────────
  # Row-major order: within each row level, iterate over all column levels.
  # This matches arrange_ggsurvplots(ncol = n_col_levels).
  if (!is.null(facet_row) && !is.null(facet_col)) {
    row_lvls <- sort(unique(stats::na.omit(as.character(plot_dat[[facet_row]]))))
    col_lvls <- sort(unique(stats::na.omit(as.character(plot_dat[[facet_col]]))))
    combos   <- expand.grid(r = row_lvls, c = col_lvls,
                            stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    # expand.grid: first arg varies fastest → c varies fastest → row-major ✓
    combos   <- combos[order(combos$r, combos$c), ]
  } else if (!is.null(facet_row)) {
    row_lvls <- sort(unique(stats::na.omit(as.character(plot_dat[[facet_row]]))))
    col_lvls <- NA_character_
    combos   <- data.frame(r = row_lvls, c = NA_character_,
                           stringsAsFactors = FALSE)
  } else if (!is.null(facet_col)) {
    col_lvls <- sort(unique(stats::na.omit(as.character(plot_dat[[facet_col]]))))
    row_lvls <- NA_character_
    combos   <- data.frame(r = NA_character_, c = col_lvls,
                           stringsAsFactors = FALSE)
  } else {
    row_lvls <- NA_character_
    col_lvls <- NA_character_
    combos   <- data.frame(r = NA_character_, c = NA_character_,
                           stringsAsFactors = FALSE)
  }

  n_row_lvls <- if (!is.null(facet_row)) length(row_lvls) else 1L
  n_col_lvls <- if (!is.null(facet_col)) length(col_lvls) else 1L
  has_facet  <- !is.null(facet_row) || !is.null(facet_col)

  # ── auto PDF dimensions ───────────────────────────────────────────────────
  if (is.null(pdf_width))
    pdf_width <- if (has_facet) max(5, 4.5 * n_col_lvls)
                 else           max(5, 4.5 * length(features))
  if (is.null(pdf_height))
    pdf_height <- if (has_facet) max(4, 4.5 * n_row_lvls)
                  else           5

  # ── main loop ─────────────────────────────────────────────────────────────
  all_plots     <- list()
  cox_results   <- list()
  sidecar_rows  <- list()

  for (feat in features) {
    sf         <- feat_map[[feat]]
    page_plots <- vector("list", nrow(combos))  # preserves NULL for blanks

    for (k in seq_len(nrow(combos))) {
      rv <- combos$r[k]
      cv <- combos$c[k]

      # subset rows
      sub <- plot_dat
      if (!is.na(rv))
        sub <- sub[!is.na(sub[[facet_row]]) & sub[[facet_row]] == rv, ]
      if (!is.na(cv))
        sub <- sub[!is.na(sub[[facet_col]]) & sub[[facet_col]] == cv, ]

      if (nrow(sub) < min_group_n * 2L) {
        page_plots[[k]] <- NULL
        sidecar_rows[[length(sidecar_rows) + 1L]] <-
          .surv_sidecar_row(feat, rv, cv, NULL)
        next
      }

      # panel title: "GeneName | RowLevel | ColLevel"
      title_parts <- feat
      if (!is.na(rv)) title_parts <- c(title_parts, rv)
      if (!is.na(cv)) title_parts <- c(title_parts, cv)
      panel_title <- paste(title_parts, collapse = " | ")

      result <- .surv_one_panel(
        dat              = sub,
        time_col         = time_col,
        status_col       = status_col,
        feat_col         = sf,
        method           = split_method,
        split_value      = split_value,
        palette          = palette,
        show_pval        = show_pval,
        show_n           = show_n,
        title            = panel_title,
        min_n            = min_group_n,
        quartile_palette = quartile_palette
      )

      page_plots[[k]] <- if (!is.null(result)) result$ggsurvplot else NULL
      sidecar_rows[[length(sidecar_rows) + 1L]] <-
        .surv_sidecar_row(feat, rv, cv, result)
    }

    all_plots[[feat]] <- page_plots

    # ── Multivariate Cox (whole cohort, continuous expression) ─────────────
    if (!is.null(multivariate_covariates)) {
      cox_results[[feat]] <- .surv_cox_fit(
        dat        = plot_dat,
        feat_col   = sf,
        time_col   = time_col,
        status_col = status_col,
        covariates = multivariate_covariates
      )
    }
  }

  # compile sidecar data frame
  sidecar_df <- do.call(
    dplyr::bind_rows,
    lapply(sidecar_rows, function(r) as.data.frame(r, stringsAsFactors = FALSE))
  )

  # ── write output ──────────────────────────────────────────────────────────
  if (!is.null(output_file)) {

    # Determine output mode:
    #   directory → output_file is an existing dir OR does not end in .pdf
    #   file      → output_file ends in .pdf
    is_dir_mode <- dir.exists(output_file) ||
                   !grepl("\\.pdf$", output_file, ignore.case = TRUE)

    if (is_dir_mode) {
      dir.create(output_file, showWarnings = FALSE, recursive = TRUE)
    }

    # arrange plots onto the current open device (caller manages pdf/dev.off)
    # arrange_ggsurvplots requires every element to be a ggsurvplot object and
    # cannot handle NULL / blank placeholders for empty facet cells.  We
    # therefore extract the $plot ggplot2 layer from each ggsurvplot and use
    # ggpubr::ggarrange, which accepts plain ggplot objects including blanks.
    .arrange_feat <- function(feat) {
      page  <- all_plots[[feat]]
      blank <- ggplot2::ggplot() + ggplot2::theme_void()

      if (has_facet) {
        valid <- Filter(Negate(is.null), page)
        if (length(valid) == 0L) return(invisible(NULL))
        plot_list <- lapply(page, function(p)
          if (is.null(p)) blank else p$plot
        )
        grid <- ggpubr::ggarrange(
          plotlist = plot_list,
          ncol     = n_col_lvls,
          nrow     = n_row_lvls
        )
        print(ggpubr::annotate_figure(
          grid,
          top = ggpubr::text_grob(feat, face = "bold", size = 13)
        ))
      } else {
        valid <- Filter(Negate(is.null), page)
        if (length(valid) == 0L) return(invisible(NULL))
        plot_list <- lapply(valid, function(p) p$plot)
        grid <- ggpubr::ggarrange(
          plotlist = plot_list,
          ncol     = length(plot_list),
          nrow     = 1L
        )
        print(grid)
      }
    }

    # ── shared metadata for JSON and .legend ─────────────────────────────────
    pkg_versions <- vapply(
      c("survival", "survminer", "recount3", "TCGAbiolinks"),
      function(p) tryCatch(as.character(utils::packageVersion(p)),
                           error = function(e) NA_character_),
      character(1L)
    )
    n_total_all  <- ncol(expr_mat)
    analysis_date <- format(Sys.Date())

    # helper: write forest plot to a separate PDF, return path
    .write_cox_pdf <- function(feat, cox_result, pdf_path_base, pw, ph) {
      cr <- cox_result
      if (is.null(cr)) return(NULL)
      cox_path <- sub("\\.pdf$", "_cox.pdf", pdf_path_base, ignore.case = TRUE)
      fp <- tryCatch(
        survminer::ggforest(
          cr$fit, data = cr$data,
          main = paste0(feat, "  —  Multivariate Cox  (n = ", cr$n, ")")
        ),
        error = function(e) {
          warning("[SurvPlot] ggforest failed for '", feat, "': ",
                  conditionMessage(e), call. = FALSE)
          NULL
        }
      )
      if (is.null(fp)) return(NULL)
      grDevices::pdf(cox_path, width = max(pw, 7), height = max(ph, 6),
                     onefile = FALSE)
      print(fp)
      grDevices::dev.off()
      message("[SurvPlot] Cox forest plot: ", cox_path)
      cox_path
    }

    # helper: build JSON list for one feature
    .surv_json_payload <- function(feat, feat_sidecar, cox_result,
                                   pdf_path, cox_path) {
      valid_panels <- feat_sidecar[!is.na(feat_sidecar$pval), ]
      pvals        <- if (nrow(valid_panels) > 0L) valid_panels$pval else NULL

      cox_summary <- if (!is.null(cox_result)) {
        sm <- tryCatch(
          as.data.frame(summary(cox_result$fit)$coefficients),
          error = function(e) NULL
        )
        if (!is.null(sm)) {
          list(
            n_samples   = cox_result$n,
            covariates  = as.list(multivariate_covariates),
            coefficients = lapply(seq_len(nrow(sm)), function(i)
              list(term      = rownames(sm)[i],
                   coef      = sm[i, "coef"],
                   exp_coef  = sm[i, "exp(coef)"],
                   se_coef   = sm[i, "se(coef)"],
                   pval      = sm[i, "Pr(>|z|)"])
            )
          )
        } else NULL
      } else NULL

      list(
        date             = analysis_date,
        function_name    = "SurvPlot",
        feature          = feat,
        feature_type     = feature_type,
        cohort           = if (!is.null(tcga_project)) tcga_project else "custom",
        backend          = if (!is.null(tcga_project)) tcga_backend else "custom",
        n_samples_total  = n_total_all,
        n_features       = length(features),
        time_col         = time_col,
        status_col       = status_col,
        split_method     = split_method,
        split_per_group  = split_per_group,
        facet_row        = if (is.null(facet_row)) NA_character_ else facet_row,
        facet_col        = if (is.null(facet_col)) NA_character_ else facet_col,
        n_panels         = nrow(feat_sidecar),
        pval_min         = if (!is.null(pvals)) min(pvals, na.rm = TRUE) else NA_real_,
        pval_max         = if (!is.null(pvals)) max(pvals, na.rm = TRUE) else NA_real_,
        cox              = cox_summary,
        pdf_path         = pdf_path,
        cox_pdf_path     = if (!is.null(cox_path)) cox_path else NA_character_,
        package_versions = as.list(pkg_versions),
        methods_text     = .surv_legend_text(
          feat            = feat,
          feature_type    = feature_type,
          feat_sidecar    = feat_sidecar,
          project         = tcga_project,
          time_col        = time_col,
          status_col      = status_col,
          split_method    = split_method,
          split_per_group = split_per_group,
          facet_row       = facet_row,
          facet_col       = facet_col,
          has_cox         = !is.null(cox_result),
          cox_covariates  = multivariate_covariates
        )
      )
    }

    if (is_dir_mode) {
      # Build a reusable facet tag for file names, e.g. "by_gender_x_Subtype"
      facet_tag <- if (!is.null(facet_row) && !is.null(facet_col))
        paste0("by_",
               gsub("[^A-Za-z0-9._-]", "_", facet_row), "_x_",
               gsub("[^A-Za-z0-9._-]", "_", facet_col))
      else if (!is.null(facet_row))
        paste0("by_", gsub("[^A-Za-z0-9._-]", "_", facet_row))
      else if (!is.null(facet_col))
        paste0("by_", gsub("[^A-Za-z0-9._-]", "_", facet_col))
      else
        NULL

      # ── directory mode: one PDF (+ optional Cox PDF) per feature ───────────
      for (feat in features) {
        safe_name    <- gsub("[^A-Za-z0-9._-]", "_", feat)
        fname_stem   <- if (!is.null(facet_tag))
          paste0(safe_name, "_", facet_tag, "_survival")
        else
          paste0(safe_name, "_survival")
        pdf_path     <- file.path(output_file, paste0(fname_stem, ".pdf"))
        feat_sidecar <- sidecar_df[sidecar_df$feature == feat, , drop = FALSE]

        grDevices::pdf(pdf_path, width = pdf_width, height = pdf_height,
                       onefile = FALSE)
        .arrange_feat(feat)
        grDevices::dev.off()
        message("[SurvPlot] Written: ", pdf_path)

        # .legend sidecar
        .write_legend_sidecar(
          pdf_path,
          .surv_legend_text(
            feat            = feat,
            feature_type    = feature_type,
            feat_sidecar    = feat_sidecar,
            project         = tcga_project,
            time_col        = time_col,
            status_col      = status_col,
            split_method    = split_method,
            split_per_group = split_per_group,
            facet_row       = facet_row,
            facet_col       = facet_col,
            has_cox         = !is.null(cox_results[[feat]]),
            cox_covariates  = multivariate_covariates
          )
        )

        # separate Cox PDF
        cox_path <- .write_cox_pdf(feat, cox_results[[feat]],
                                   pdf_path, pdf_width, pdf_height)

        # .json methods file
        if (requireNamespace("jsonlite", quietly = TRUE)) {
          json_path <- sub("\\.pdf$", ".json", pdf_path, ignore.case = TRUE)
          jsonlite::write_json(
            .surv_json_payload(feat, feat_sidecar, cox_results[[feat]],
                               pdf_path, cox_path),
            path = json_path, auto_unbox = TRUE, pretty = TRUE
          )
          message("[SurvPlot] Methods JSON: ", json_path)
        }
      }

      if (write_sidecar) {
        sc_file <- file.path(output_file, "survival_sidecar.csv")
        utils::write.csv(sidecar_df, sc_file, row.names = FALSE)
        message("[SurvPlot] Sidecar CSV: ", sc_file)
      }

    } else {
      # ── single multi-page PDF mode ─────────────────────────────────────────
      grDevices::pdf(output_file, width = pdf_width, height = pdf_height,
                     onefile = TRUE)
      for (feat in features) .arrange_feat(feat)
      grDevices::dev.off()
      message("[SurvPlot] PDF written to:\n  ", output_file)

      # .legend sidecar (one for the whole PDF — summarises all features)
      all_valid <- sidecar_df[!is.na(sidecar_df$pval), ]
      pval_str  <- if (nrow(all_valid) > 0L)
        sprintf("log-rank p range %.3g – %.3g",
                min(all_valid$pval, na.rm = TRUE),
                max(all_valid$pval, na.rm = TRUE))
      else "p-values not available"

      combined_legend <- paste0(
        "Kaplan-Meier survival analysis of ", length(features), " ",
        feature_type, "(s): ", paste(features, collapse = ", "), ". ",
        "Cohort: ", if (!is.null(tcga_project)) tcga_project else "custom",
        " (n = ", n_total_all, " samples). ",
        "Survival endpoint: ", time_col, " (days); event indicator: ", status_col, ". ",
        "Stratification: ", split_method,
        if (split_per_group) " (per facet group)" else " (global)", ". ",
        "Panel layout: ",
        if (!is.null(facet_row) && !is.null(facet_col))
          paste0("rows: ", facet_row, "; columns: ", facet_col)
        else if (!is.null(facet_row)) paste0("rows: ", facet_row)
        else if (!is.null(facet_col)) paste0("columns: ", facet_col)
        else "no faceting", ". ",
        pval_str, ".",
        if (!is.null(multivariate_covariates))
          paste0(" Multivariate Cox results in companion _cox.pdf file.")
        else "",
        " Full details in companion sidecar CSV."
      )
      .write_legend_sidecar(output_file, combined_legend)

      # separate Cox PDF (one PDF for all features)
      if (!is.null(multivariate_covariates) &&
          any(!vapply(cox_results, is.null, logical(1L)))) {
        cox_path <- sub("\\.pdf$", "_cox.pdf", output_file, ignore.case = TRUE)
        grDevices::pdf(cox_path, width = max(pdf_width, 7),
                       height = max(pdf_height, 6), onefile = TRUE)
        for (feat in features) {
          cr <- cox_results[[feat]]
          if (is.null(cr)) next
          sf <- feat_map[[feat]]
          fp <- tryCatch(
            survminer::ggforest(
              cr$fit, data = cr$data,
              main = paste0(feat, "  —  Multivariate Cox  (n = ",
                            cr$n, ")")
            ),
            error = function(e) NULL
          )
          if (!is.null(fp)) print(fp)
        }
        grDevices::dev.off()
        message("[SurvPlot] Cox forest plots: ", cox_path)
        .write_legend_sidecar(
          cox_path,
          paste0(
            "Multivariate Cox proportional hazards regression for ",
            length(features), " ", feature_type, "(s). ",
            "Covariates: ", paste(multivariate_covariates, collapse = ", "), ". ",
            "Cohort: ", if (!is.null(tcga_project)) tcga_project else "custom",
            " (n = ", n_total_all, " samples). ",
            "Each panel shows hazard ratios and 95% CI from survival::coxph, ",
            "visualised with survminer::ggforest."
          )
        )
      } else {
        cox_path <- NULL
      }

      # .json (one for the whole multi-feature PDF)
      if (requireNamespace("jsonlite", quietly = TRUE)) {
        json_path <- sub("\\.pdf$", ".json", output_file, ignore.case = TRUE)
        jsonlite::write_json(
          list(
            date             = analysis_date,
            function_name    = "SurvPlot",
            features         = as.list(features),
            feature_type     = feature_type,
            cohort           = if (!is.null(tcga_project)) tcga_project else "custom",
            backend          = if (!is.null(tcga_project)) tcga_backend else "custom",
            n_samples_total  = n_total_all,
            time_col         = time_col,
            status_col       = status_col,
            split_method     = split_method,
            split_per_group  = split_per_group,
            facet_row        = if (is.null(facet_row)) NA_character_ else facet_row,
            facet_col        = if (is.null(facet_col)) NA_character_ else facet_col,
            multivariate_covariates = if (is.null(multivariate_covariates))
                                        NA_character_ else as.list(multivariate_covariates),
            pdf_path         = output_file,
            cox_pdf_path     = if (!is.null(cox_path)) cox_path else NA_character_,
            package_versions = as.list(pkg_versions),
            methods_text     = combined_legend
          ),
          path = json_path, auto_unbox = TRUE, pretty = TRUE
        )
        message("[SurvPlot] Methods JSON: ", json_path)
      }

      if (write_sidecar) {
        sc_file <- sub("\\.pdf$", "_sidecar.csv", output_file, ignore.case = TRUE)
        utils::write.csv(sidecar_df, sc_file, row.names = FALSE)
        message("[SurvPlot] Sidecar CSV: ", sc_file)
      }
    }
  }

  invisible(list(plots = all_plots, sidecar = sidecar_df))
}
