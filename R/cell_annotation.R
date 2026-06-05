# =============================================================================
# scSidekick cell-type annotation helpers
#
# .sctype_gene_sets_prepare()  — bundled from scType (MIT)
# .sctype_score()              — bundled from scType (MIT)
# .get_sctype_db()             — locates / downloads ScTypeDB_full.xlsx
# CellTypeAssignmentHelper()   — runs rcellmarker + scType + SingleR,
#                                produces dotplots and UMAP summary figures
#
# scType source:  https://github.com/IanevskiAleksandr/sc-type
# scType citation: Ianevski A, Giri AK, Aittokallio T.
#   "Fully-automated and ultra-fast cell-type identification using specific
#    marker combinations from single-cell transcriptomic data."
#   Nature Communications, 2022. https://doi.org/10.1038/s41467-022-28803-w
# =============================================================================

# ---------------------------------------------------------------------------
# .sctype_gene_sets_prepare()
# Reads a ScTypeDB Excel file and returns named gene-set lists.
# Bundled verbatim from https://github.com/IanevskiAleksandr/sc-type (MIT).
# ---------------------------------------------------------------------------
.sctype_gene_sets_prepare <- function(path_to_db_file, cell_type) {
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Package 'openxlsx' is required for scType. ",
         "Install with: install.packages('openxlsx')")

  cm <- openxlsx::read.xlsx(path_to_db_file)
  cm <- cm[cm$tissueType %in% cell_type & !is.na(cm$tissueType), ]

  if (nrow(cm) == 0) {
    warning("No markers found for tissue type(s): ",
            paste(cell_type, collapse = ", "))
    return(list(gs_positive = list(), gs_negative = list()))
  }

  .parse <- function(x) {
    if (is.na(x) || trimws(x) == "") return(character(0))
    g <- trimws(unlist(strsplit(as.character(x), ",")))
    g[nchar(g) > 0]
  }

  ct   <- unique(cm$cellName)
  make_gs <- function(col) {
    stats::setNames(lapply(ct, function(n) {
      unique(unlist(lapply(cm[[col]][cm$cellName == n], .parse)))
    }), ct)
  }

  gs_pos <- Filter(function(x) length(x) > 0, make_gs("geneSymbolmore1"))
  gs_neg <- Filter(function(x) length(x) > 0, make_gs("geneSymbolmore2"))

  list(gs_positive = gs_pos, gs_negative = gs_neg)
}


# ---------------------------------------------------------------------------
# .sctype_score()
# Scores cells against gene sets.
# Bundled verbatim from https://github.com/IanevskiAleksandr/sc-type (MIT).
# ---------------------------------------------------------------------------
.sctype_score <- function(scRNAseqData, scaled = TRUE, gs, gs2 = NULL,
                           gene_names_to_uppercase = TRUE) {
  if (gene_names_to_uppercase) {
    rownames(scRNAseqData) <- toupper(rownames(scRNAseqData))
    gs  <- lapply(gs,  function(x) unique(toupper(x)))
    if (!is.null(gs2))
      gs2 <- lapply(gs2, function(x) unique(toupper(x)))
  }

  n_cells <- ncol(scRNAseqData)

  es <- vapply(names(gs), function(gset) {
    pos <- intersect(gs[[gset]], rownames(scRNAseqData))
    neg <- if (!is.null(gs2) && gset %in% names(gs2))
      intersect(gs2[[gset]], rownames(scRNAseqData)) else character(0)

    if (length(pos) == 0) return(rep(0, n_cells))

    pos_sc <- colSums(scRNAseqData[pos, , drop = FALSE]) / sqrt(length(pos))
    neg_sc <- if (length(neg) > 0)
      colSums(scRNAseqData[neg, , drop = FALSE]) / sqrt(length(neg))
    else
      rep(0, n_cells)

    pos_sc - neg_sc
  }, numeric(n_cells))

  # vapply produces n_cells × n_sets → transpose to n_sets × n_cells
  if (is.vector(es)) {
    matrix(es, nrow = 1,
           dimnames = list(names(gs), colnames(scRNAseqData)))
  } else {
    t(es)
  }
}


# ---------------------------------------------------------------------------
# .get_sctype_db()
# Returns a path to ScTypeDB_full.xlsx:
#   1. User-supplied path (sctype_db arg)
#   2. inst/extdata/ (bundled with scSidekick)
#   3. R user cache (download once on first use)
# ---------------------------------------------------------------------------
.get_sctype_db <- function(sctype_db = NULL) {
  # Helper: reject files that are 0 bytes (e.g. un-synced Dropbox placeholders)
  .valid_xlsx <- function(p) nzchar(p) && file.exists(p) && file.size(p) > 0

  if (!is.null(sctype_db) && .valid_xlsx(sctype_db)) return(sctype_db)

  bundled <- system.file("extdata", "ScTypeDB_full.xlsx", package = "scSidekick")
  if (.valid_xlsx(bundled)) return(bundled)

  cache_dir  <- tools::R_user_dir("scSidekick", "cache")
  cache_file <- file.path(cache_dir, "ScTypeDB_full.xlsx")

  if (!file.exists(cache_file)) {
    message("ScTypeDB_full.xlsx not found — downloading to user cache ",
            "(one-time setup)...")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    url <- paste0("https://raw.githubusercontent.com/",
                  "IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx")
    tryCatch(
      utils::download.file(url, cache_file, mode = "wb", quiet = FALSE),
      error = function(e) stop(
        "Could not download ScTypeDB_full.xlsx.\n",
        "Download it manually from:\n  ", url,
        "\nand place it at:\n  ", cache_file,
        "\nor pass sctype_db = '<path>' directly."
      )
    )
    message("Cached to: ", cache_file)
  }

  cache_file
}


# ---------------------------------------------------------------------------
# CellTypeAssignmentHelper
# ---------------------------------------------------------------------------

#' Automated cell-type annotation with rcellmarker, scType, and SingleR
#'
#' Runs up to three automated cell-type annotation methods on a Seurat object,
#' adds their results as new metadata columns, and produces two summary PDFs:
#' \enumerate{
#'   \item An `AutoAssignment umaps.pdf` with scType, rcellmarker, and cluster
#'         UMAPs side-by-side.
#'   \item An `AutoAssignment umaps with singleR.pdf` + `Assignment Helper.pdf`
#'         combining all annotations and canonical-marker dotplots.
#' }
#'
#' @param seurat_object A Seurat object.
#' @param markers Data frame of DE results from [presto::wilcoxauc()]. Required
#'   for rcellmarker annotation and for the top-marker dotplot. Expected columns:
#'   `feature`, `group`, `avgExpr`, `avg_log2FC` (or `logFC`), `auc`,
#'   `pval`, `padj`, `pct_in`, `pct_out`.
#' @param cluster_column Character. Metadata column holding cluster labels.
#'   Default `"Cluster"`.
#' @param output_dir Character. Directory for PDFs.
#' @param object_name Character. Prefix for PDF filenames.
#' @param subset_name Character. Second prefix component.
#' @param reduction Character. UMAP reduction name. Default `"umap"`.
#' @param cluster_colors Named character vector of colours for `cluster_column`
#'   levels. `NULL` auto-generates from `Nour_pal()`.
#' @param species Character. `"human"` (default) or `"mouse"`. Controls:
#'   (1) which gene column is read from the markers CSV (`Human` vs `Mouse`),
#'   (2) the rcellmarker `species` argument, and
#'   (3) which CellChatDB is loaded in downstream CellChat workflows.
#'
#' @section rcellmarker:
#' @param run_rcellmarker Logical. Run rcellmarker automated annotation?
#'   Requires the `rcellmarker` package. Default `TRUE`.
#'
#' @section scType:
#' @param run_sctype Logical. Run scType automated annotation? Requires
#'   `openxlsx`. Default `TRUE`.
#' @param sctype_tissues Character vector of tissue type(s) to query from the
#'   bundled `ScTypeDB_full.xlsx`. Combine multiple types for broader coverage.
#'   All valid values (from the bundled database):
#'   \itemize{
#'     \item `"Adrenal"` — adrenal gland
#'     \item `"Brain"` — cerebral cortex, broad neuronal/glial types
#'     \item `"Eye"` — retinal and eye-specific cell types
#'     \item `"Heart"` — cardiac cell types
#'     \item `"Hippocampus"` — hippocampal-specific subtypes (use with `"Brain"`
#'       for combined brain data)
#'     \item `"Immune system"` — broad immune cells: T, B, NK, myeloid, etc.
#'     \item `"Intestine"` — intestinal epithelial and stromal types
#'     \item `"Kidney"` — renal cell types
#'     \item `"Liver"` — hepatocytes, Kupffer cells, stellate cells
#'     \item `"Lung"` — airway epithelium, alveolar cells, immune
#'     \item `"Muscle"` — skeletal and smooth muscle
#'     \item `"Pancreas"` — islet and exocrine types
#'     \item `"Placenta"` — trophoblast and decidual types
#'     \item `"Spleen"` — splenic immune populations
#'     \item `"Stomach"` — gastric epithelial types
#'     \item `"Thymus"` — thymocyte and thymic stromal types
#'   }
#'   Recommended combinations: `c("Brain", "Hippocampus")` for brain snRNAseq;
#'   `c("Immune system", "Spleen")` for immune-rich tissues;
#'   `c("Immune system", "Brain")` for neuroinflammation datasets.
#' @param sctype_assay Character. Seurat assay to use for scoring.
#'   Must contain a `scale.data` (scaled) layer. Common choices:
#'   \itemize{
#'     \item `"RNA"` — standard full assay (default); use `sctype_layer = "scale.data"`
#'     \item `"sketch"` — BPCells sketch assay for large datasets; use
#'       `sctype_layer = "scale.data"` after `ScaleData()` on the sketch
#'   }
#' @param sctype_layer Character. Layer within `sctype_assay` that holds
#'   scaled expression values. Default `"scale.data"`. Must be a scaled
#'   (mean-centred) matrix — scType scores are not meaningful on raw counts.
#' @param sctype_db Character or `NULL`. Path to a `ScTypeDB_full.xlsx`
#'   database file. `NULL` (default) uses the copy bundled in
#'   `inst/extdata/ScTypeDB_full.xlsx`; if that is missing it is downloaded
#'   once from the scType GitHub and cached at
#'   `tools::R_user_dir("scSidekick", "cache")`.
#'
#' @section SingleR:
#' @param run_singler Logical. Run SingleR cluster-level annotation? Requires
#'   `SingleR`, `scuttle`, `celldex`, `SingleCellExperiment`. Default `TRUE`.
#' @param singler_ref A reference `SummarizedExperiment` for
#'   [SingleR::SingleR()], or `NULL` to use [celldex::BlueprintEncodeData()]
#'   (default). Choose based on tissue type and species:
#'   \describe{
#'     \item{\strong{Human references}:}{}
#'     \item{`celldex::BlueprintEncodeData()`}{Blueprint/ENCODE bulk RNA-seq.
#'       Broad cell types; good general-purpose starting point. **Default.**}
#'     \item{`celldex::HumanPrimaryCellAtlasData()`}{Human Primary Cell Atlas.
#'       Very broad coverage across many tissues; best for highly heterogeneous
#'       datasets with diverse cell populations.}
#'     \item{`celldex::MonacoImmuneData()`}{Monaco et al. Finest immune subtype
#'       resolution (29 populations). Best for PBMCs, spleen, or any
#'       immune-focused study.}
#'     \item{`celldex::DatabaseImmuneCellExpressionData()`}{DICE database.
#'       Good alternative immune reference with activated/resting states.}
#'     \item{`celldex::NovershternHematopoieticData()`}{Hematopoietic lineages.
#'       Best for bone marrow or stem-cell datasets.}
#'     \item{\strong{Mouse references}:}{}
#'     \item{`celldex::ImmGenData()`}{ImmGen consortium. The most comprehensive
#'       mouse immune cell reference (over 250 subtypes).}
#'     \item{`celldex::MouseRNAseqData()`}{Broad mouse tissue types from bulk
#'       RNA-seq. Good for non-immune mouse tissues.}
#'   }
#'   All celldex references are downloaded once and cached by
#'   `BiocFileCache`; subsequent calls are instant.
#' @param singler_col Character. Name for the new metadata column storing
#'   numbered SingleR labels (e.g. `"01.CD8+ T cells"`).
#'   Default `"SingleRAssignment"`.
#'
#' @section Dotplots:
#' @param markers_csv Character or `NULL`. Path to a canonical-marker CSV.
#'   `NULL` uses the CSV bundled with scSidekick (`inst/extdata/Markers.csv`).
#'   Required columns: `Type`, `Cell types`, plus `Human` and/or `Mouse`
#'   gene symbol columns. The bundled CSV covers brain and immune markers
#'   (see `marker_groups` for the full list of `Type` values).
#' @param marker_groups Named list controlling how many dotplots are produced
#'   and which `Type` values from the markers CSV appear in each. One
#'   [SplitDotPlot()] is generated per list entry. `NULL` produces a single
#'   dotplot with all markers combined.
#'
#'   **`Type` values in the bundled `Markers.csv`:**
#'   \itemize{
#'     \item `"Brain Stroma"` — Astrocytes, Oligodendrocytes, OPC, Microglia,
#'       Endothelial, Pericytes, Neurons (excitatory, inhibitory, dopaminergic,
#'       glutamatergic, serotonergic), VLMC, and stem/progenitor types
#'     \item `"Lymphocytes"` — T cells (naive, cytotoxic, Th1/2/17, Treg,
#'       exhausted, MAIT, GDT), B cells, NK cells
#'     \item `"Monocytes"` — CD14/CD16 monocytes, MDSCs
#'     \item `"Macrophages"` — tissue macrophages
#'     \item `"Dendritic cells"` — cDC1, cDC2, pDCs, LAMP3+ DCs
#'     \item `"Granulocytes"` — neutrophils, basophils, mast cells
#'     \item `"Bone Marrow Derived"` — BMD immune cells
#'     \item `"All"` — proliferating cells (Prol)
#'     \item `"Tumor"` — generic tumor marker
#'   }
#'
#'   Example for brain snRNAseq with separate immune panel:
#'   ```r
#'   marker_groups = list(
#'     "Brain"  = "Brain Stroma",
#'     "Immune" = c("Lymphocytes", "Monocytes", "Macrophages",
#'                  "Dendritic cells", "Granulocytes",
#'                  "Bone Marrow Derived", "All")
#'   )
#'   ```
#' @param top_n_per_cluster Integer. Number of top DE genes per cluster shown
#'   in the cluster-faceted dotplot (uses the `markers` data frame).
#'   Default `3`.
#'
#' @section UMAP panels:
#' @param extra_dimplot_cols Character vector of additional metadata columns
#'   to include as panels in the multi-UMAP summary figure (e.g. `"Subclass"`
#'   for SEAAD reference labels or `"Donor.ID"` for QC). Colours are
#'   auto-generated from `Nour_pal()` unless supplied via `extra_col_colors`.
#' @param extra_col_colors Named list of colour vectors for `extra_dimplot_cols`.
#'   Names must match `extra_dimplot_cols`. Example:
#'   `list(Subclass = subclass_colors)`.
#'
#' @return The modified Seurat object (returned invisibly), with:
#' \describe{
#'   \item{New metadata columns}{\code{scType_CellType},
#'     \code{rcellmarker_CellType}, and/or the column named by
#'     \code{singler_col}, depending on which methods ran.}
#'   \item{\code{seurat_object@misc$annotation_colors}}{Named list of colour vectors for
#'     each annotation column: \code{rcellmarker}, \code{sctype},
#'     \code{singler}, plus any \code{extra_dimplot_cols} colours. Retrieve
#'     with \code{seurat_object@misc$annotation_colors$sctype} etc.}
#' }
#' @export
CellTypeAssignmentHelper <- function(
    seurat_object,
    markers,
    cluster_column        = "Cluster",
    output_dir            = NULL,
    object_name           = "",
    subset_name           = "",
    reduction          = "umap",
    cluster_colors     = NULL,
    species            = "human",

    # rcellmarker
    run_rcellmarker    = TRUE,

    # scType
    run_sctype         = TRUE,
    sctype_tissues     = c("Immune system", "Brain"),
    sctype_assay       = "RNA",
    sctype_layer       = "scale.data",
    sctype_db          = NULL,

    # SingleR
    run_singler        = TRUE,
    singler_ref        = NULL,
    singler_col        = "SingleRAssignment",

    # Dotplots
    markers_csv        = NULL,
    marker_groups      = NULL,
    top_n_per_cluster  = 3L,

    # UMAP panels
    extra_dimplot_cols  = NULL,
    extra_col_colors    = NULL,

    # dotplot scale params
    col.min = 0, col.max = 2,
    scale.min = 0, scale.max = 100
) {

  # Walk up to PrepObject-stored defaults
  output_dir  <- output_dir %||% .nk_setting(seurat_object, "output_dir")
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""
  if (is.null(output_dir))
    stop("'output_dir' must be supplied, or stored via PrepObject(output_dir = ...).")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Build a clean filename prefix and use file.path() for all save paths so
  # the function works regardless of whether output_dir ends with "/" or not.
  file_prefix <- trimws(paste(object_name, subset_name))
  cols_out <- list()

  # Cluster colours: PrepObject stored → auto Nour_pal fallback
  if (is.null(cluster_colors)) {
    cluster_colors <- .nk_colors(seurat_object, cluster_column) %||% {
      cl_lvls <- if (is.factor(seurat_object@meta.data[[cluster_column]]))
        levels(seurat_object@meta.data[[cluster_column]])
      else sort(unique(as.character(seurat_object@meta.data[[cluster_column]])))
      n <- length(cl_lvls)
      stats::setNames(
        Nour_pal(if (n <= 8) "all" else "spectrum")(n),
        cl_lvls
      )
    }
  }

  # Gene column in markers CSV
  gene_col_csv <- if (species == "mouse") "Mouse" else "Human"

  # ── 1. rcellmarker ─────────────────────────────────────────────────────────
  rcm_cols <- NULL
  if (run_rcellmarker) {
    if (!requireNamespace("rcellmarker", quietly = TRUE)) {
      warning("Package 'rcellmarker' not installed — skipping rcellmarker.")
      run_rcellmarker <- FALSE
    }
  }

  if (run_rcellmarker) {
    message("Running rcellmarker...")
    tryCatch({
      # Map presto/wilcoxauc columns to the names rcellmarker expects
      m_rcm <- markers
      col_map <- c(feature = "gene", group = "cluster", mean = "avgExpr",
                   auc = "avg_log2FC", pvalue = "pval", padj = "p_val_adj",
                   pct_in = "pct_in", pct_out = "pct_out")
      for (from in names(col_map)) {
        to <- col_map[[from]]
        if (from %in% colnames(m_rcm) && !to %in% colnames(m_rcm))
          colnames(m_rcm)[colnames(m_rcm) == from] <- to
      }
      res_rcm <- rcellmarker::cellMarker(m_rcm, type = "seurat", species = species,
                                         topn = 1, keytype = "SYMBOL", weight = 0,
                                         padj = 0.05, tissue = NULL, minSize = 2)
      if (!is.null(res_rcm)) {
        xxy  <- res_rcm$cellType
        names(xxy) <- res_rcm$cluster
        xxy2 <- xxy[as.character(seurat_object@meta.data[[cluster_column]])]
        names(xxy2) <- rownames(seurat_object@meta.data)
        seurat_object$rcellmarker_CellType <- factor(xxy2)
        rcm_lvls <- levels(seurat_object$rcellmarker_CellType)
        rcm_cols <- stats::setNames(
          Nour_pal(if (length(rcm_lvls) <= 18L) "all" else "spectrum")(length(rcm_lvls)),
          rcm_lvls)
        cols_out$rcellmarker <- rcm_cols
      } else {
        run_rcellmarker <- FALSE
      }
    }, error = function(e) {
      warning("rcellmarker failed: ", conditionMessage(e), " — skipping.")
      run_rcellmarker <<- FALSE
    })
  }

  # ── 2. scType ──────────────────────────────────────────────────────────────
  sct_cols <- NULL
  if (run_sctype) {
    message("Running scType...")
    db_path <- tryCatch(.get_sctype_db(sctype_db), error = function(e) {
      warning("scType DB unavailable: ", conditionMessage(e)); NULL })

    if (!is.null(db_path)) {
      gs_list <- tryCatch(
        .sctype_gene_sets_prepare(db_path, sctype_tissues),
        error = function(e) {
          warning("scType gene_sets_prepare failed: ", conditionMessage(e)); NULL
        }
      )

      if (!is.null(gs_list) && length(gs_list$gs_positive) > 0) {
        expr_mat <- tryCatch(
          .get_layer_data(seurat_object, assay = sctype_assay, layer = sctype_layer),
          error = function(e) {
            warning("Could not get assay data for scType (assay='", sctype_assay,
                    "', layer='", sctype_layer, "'): ", conditionMessage(e)); NULL
          }
        )

        if (!is.null(expr_mat)) {
          es_max <- .sctype_score(
            scRNAseqData = expr_mat, scaled = TRUE,
            gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
          )

          # es_max columns are the cells in the assay used for scoring
          # (e.g. the 50k sketched cells). seurat_object@meta.data contains ALL cells.
          # Intersect to avoid "subscript out of bounds" when the assay is
          # a subset of the full object.
          scored_cells <- colnames(es_max)

          cl_levels <- unique(seurat_object@meta.data[[cluster_column]])
          cl_res <- do.call(rbind, lapply(cl_levels, function(cl) {
            meta      <- seurat_object@meta.data
            cells_cl  <- rownames(meta[meta[[cluster_column]] == cl, ])
            cells_cl  <- intersect(cells_cl, scored_cells)
            if (length(cells_cl) == 0)
              return(data.frame(cluster = cl, type = "Unknown",
                                scores = 0, ncells = 0, row.names = NULL))
            sc <- sort(rowSums(es_max[, cells_cl, drop = FALSE]),
                       decreasing = TRUE)
            head(data.frame(cluster = cl, type = names(sc),
                            scores  = sc,
                            ncells  = length(cells_cl),
                            row.names = NULL), 10)
          }))

          sctype_scores <- cl_res |>
            dplyr::group_by(cluster) |>
            dplyr::top_n(n = 1, wt = scores) |>
            dplyr::ungroup() |>
            dplyr::arrange(cluster)
          colnames(sctype_scores) <- c("cluster", "sctypeType",
                                       "sctypeScore", "sctypeNcells")

          xxy  <- sctype_scores$sctypeType
          names(xxy) <- sctype_scores$cluster
          xxy2 <- xxy[as.character(seurat_object@meta.data[[cluster_column]])]
          names(xxy2) <- rownames(seurat_object@meta.data)
          seurat_object$scType_CellType <- factor(xxy2)

          sct_lvls <- levels(seurat_object$scType_CellType)
          sct_cols <- stats::setNames(
            Nour_pal(if (length(sct_lvls) <= 18L) "all" else "spectrum")(length(sct_lvls)),
            sct_lvls)
          cols_out$sctype <- sct_cols
        } else { run_sctype <- FALSE }
      } else { run_sctype <- FALSE }
    } else { run_sctype <- FALSE }
  }

  # ── 3. AutoAssignment UMAP (scType + rcellmarker + Cluster) ───────────────
  {
    panels_auto <- list()
    .mk_dimplot <- function(grp, cols, ttl)
      Seurat::DimPlot(seurat_object, reduction = reduction, group.by = grp,
                      label = TRUE, cols = cols, pt.size = 0.3) +
        ggplot2::ggtitle(ttl) +
        Seurat::NoLegend()
    if (run_sctype)
      panels_auto[["scType"]] <- .mk_dimplot("scType_CellType", sct_cols, "scType")
    if (run_rcellmarker)
      panels_auto[["rcellmarker"]] <- .mk_dimplot("rcellmarker_CellType", rcm_cols, "rcellmarker")
    panels_auto[["Clusters"]] <- .mk_dimplot(cluster_column, cluster_colors, "Clusters")

    if (length(panels_auto) > 0) {
      f_auto <- file.path(output_dir, paste0(file_prefix, " AutoAssignment umaps.pdf"))
      grDevices::pdf(f_auto, width = 5 * length(panels_auto), height = 6)
      print(patchwork::wrap_plots(panels_auto, nrow = 1))
      grDevices::dev.off()
      message("  Saved: ", f_auto)
      .write_legend_sidecar(f_auto, paste0(
        "UMAP projections of ", object_name, subset_name,
        " coloured by automated cell-type annotation results. ",
        if (run_sctype)
          "scType scores cells against curated tissue-specific gene sets. " else "",
        if (run_rcellmarker)
          "rcellmarker assigns cluster-level labels from published marker databases. " else "",
        "Clusters are shown for reference."
      ))
    }
  }

  # ── 4. SingleR ─────────────────────────────────────────────────────────────
  singler_cols_out <- NULL
  p_singler        <- NULL

  if (run_singler) {
    for (pkg in c("SingleR", "scuttle", "celldex", "SingleCellExperiment")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        warning("Package '", pkg,
                "' not installed — skipping SingleR. Install from Bioconductor.")
        run_singler <- FALSE
        break
      }
    }
  }

  if (run_singler) {
    message("Running SingleR...")
    ref <- if (is.null(singler_ref)) {
      tryCatch(celldex::BlueprintEncodeData(),
               error = function(e) {
                 warning("Could not load BlueprintEncodeData: ", conditionMessage(e))
                 NULL
               })
    } else if (is.character(singler_ref)) {
      # Accept a celldex function name as a string, e.g. "HumanPrimaryCellAtlasData"
      fn_name <- singler_ref
      if (!grepl("::", fn_name, fixed = TRUE)) fn_name <- paste0("celldex::", fn_name)
      tryCatch(
        eval(parse(text = paste0(fn_name, "()"))),
        error = function(e) {
          warning("Could not load SingleR reference '", singler_ref,
                  "': ", conditionMessage(e), " — falling back to BlueprintEncodeData.")
          tryCatch(celldex::BlueprintEncodeData(), error = function(e2) NULL)
        }
      )
    } else singler_ref

    if (!is.null(ref)) {
      # Build a minimal SCE from the RNA assay only.
      # Seurat::as.SingleCellExperiment() attaches ALL assays as altExps,
      # which fails when different assays have different numbers of cells
      # (e.g. a full RNA assay + a 50k sketch assay in the same object).
      # We need only the normalised RNA counts for SingleR, so we create
      # the SCE manually and avoid the altExp dimension mismatch entirely.
      expr_singler <- tryCatch(
        .get_layer_data(seurat_object, assay = "RNA", layer = "data"),
        error = function(e) NULL
      )
      if (is.null(expr_singler)) {
        warning("Could not extract RNA data layer for SingleR — skipping.")
        run_singler <- FALSE
      }

      sce <- if (run_singler) {
        SingleCellExperiment::SingleCellExperiment(
          assays  = list(logcounts = expr_singler),
          colData = S4Vectors::DataFrame(
            row.names = colnames(expr_singler),
            cluster   = seurat_object@meta.data[colnames(expr_singler),
                                      cluster_column, drop = TRUE]
          )
        )
      } else NULL

      se_agg <- if (run_singler) tryCatch(
        scuttle::summarizeAssayByGroup(
          sce, ids = sce$cluster, assay.type = "logcounts"),
        error = function(e) { warning(conditionMessage(e)); NULL }
      ) else NULL

      if (!is.null(se_agg)) {
        anno <- tryCatch(
          SingleR::SingleR(se_agg, ref = ref, labels = ref$label.fine,
                           assay.type.test = "mean"),
          error = function(e) { warning(conditionMessage(e)); NULL })

        if (!is.null(anno)) {
          sr_assign        <- anno$labels
          names(sr_assign) <- rownames(anno)
          seurat_object[[singler_col]] <- as.character(sr_assign[seurat_object@meta.data[[cluster_column]]])

          # Number the labels: 01.CellType, 02.CellType …
          lvls   <- levels(as.factor(seurat_object[[singler_col]]))
          labels <- paste0(sprintf("%02d", seq_along(lvls)), ".", lvls)
          names(labels) <- lvls
          seurat_object[[singler_col]] <- factor(
            as.character(labels[seurat_object[[singler_col]]]))

          mycols_sr <- c(
            "#0D47A1","#3F81BD","#0F6FC6","#009DD9","#1B587C","#4BACC6",
            "#0BD0D9","#10CF9B","#7CCA62","#A5C249","#FFD000","#F79646",
            "#E65100","#FD817E","#FD625E","#C00000","#8872C4","#9B57D3",
            "#9030A0","#300890")
          singler_cols_out <- grDevices::colorRampPalette(mycols_sr)(
            length(levels(seurat_object[[singler_col]])))
          names(singler_cols_out) <- levels(seurat_object[[singler_col]])
          cols_out$singler <- singler_cols_out

          # Abbreviated numeric DimPlot
          xx <- Seurat::DimPlot(seurat_object, reduction = reduction,
                                group.by = singler_col, label = FALSE,
                                cols = singler_cols_out, label.size = 7) +
            ggplot2::guides(colour = ggplot2::guide_legend(
              override.aes = list(size = 5), ncol = 2,
              title.theme = ggplot2::element_text(size = 15, face = "bold"),
              title.position = "top",
              label.theme  = ggplot2::element_text(size = 10)))

          num_labels        <- sprintf("%02d", seq_along(lvls))
          names(num_labels) <- labels   # labels are already "01.X", "02.X"…
          xx[[1]][["data"]][["Label"]] <- factor(
            num_labels[xx[[1]][["data"]][[singler_col]]])

          p_singler <- Seurat::LabelClusters(
            plot = xx, id = "Label", repel = FALSE)
        } else { run_singler <- FALSE }
      } else { run_singler <- FALSE }
    } else { run_singler <- FALSE }
  }

  # ── 5. Multi-panel UMAP (Cluster + extra + SingleR + scType + rcellmarker) ─
  {
    umap_panels <- list()

    # Cluster (always first)
    umap_panels[["Cluster"]] <- Seurat::DimPlot(
      seurat_object, reduction = reduction, group.by = cluster_column,
      label = TRUE, cols = cluster_colors) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 5), ncol = 5))

    # Dataset-specific extra columns
    for (ecol in extra_dimplot_cols) {
      ecols_use <- if (!is.null(extra_col_colors) && ecol %in% names(extra_col_colors))
        extra_col_colors[[ecol]]
      else {
        el <- levels(as.factor(seurat_object@meta.data[[ecol]]))
        stats::setNames(Nour_pal("all")(length(el)), el)
      }
      cols_out[[ecol]] <- ecols_use
      umap_panels[[ecol]] <- Seurat::DimPlot(
        seurat_object, reduction = reduction, group.by = ecol,
        label = TRUE, cols = ecols_use) +
        ggplot2::guides(colour = ggplot2::guide_legend(
          override.aes = list(size = 5), ncol = 3))
    }

    # SingleR
    if (run_singler && !is.null(p_singler))
      umap_panels[[singler_col]] <- p_singler

    # scType
    if (run_sctype)
      umap_panels[["scType"]] <- Seurat::DimPlot(
        seurat_object, reduction = reduction, group.by = "scType_CellType",
        label = TRUE, cols = sct_cols) +
        ggplot2::guides(colour = ggplot2::guide_legend(
          override.aes = list(size = 5), ncol = 2))

    # rcellmarker
    if (run_rcellmarker)
      umap_panels[["rcellmarker"]] <- Seurat::DimPlot(
        seurat_object, reduction = reduction, group.by = "rcellmarker_CellType",
        label = TRUE, cols = rcm_cols) +
        ggplot2::guides(colour = ggplot2::guide_legend(
          override.aes = list(size = 5), ncol = 3))

    p1 <- patchwork::wrap_plots(umap_panels, ncol = length(umap_panels)) &
      theme_NourMin() &
      ggplot2::theme(legend.position = "bottom") &
      Seurat::FontSize(x.title = 20, y.title = 20, main = 20) &
      ggplot2::scale_y_continuous(breaks = NULL) &
      ggplot2::scale_x_continuous(breaks = NULL) &
      ggplot2::xlab("UMAP1") & ggplot2::ylab("UMAP2")

    f_sr <- file.path(output_dir, paste0(file_prefix, " AutoAssignment umaps with singleR.pdf"))
    grDevices::pdf(f_sr, width = 5 * length(umap_panels), height = 7)
    print(p1)
    grDevices::dev.off()
    message("  Saved: ", f_sr)
    .write_legend_sidecar(f_sr, paste0(
      "Multi-panel UMAP summary of cell-type annotation results for ",
      object_name, subset_name, ". Panels show cluster labels",
      if (length(extra_dimplot_cols) > 0)
        paste0(", ", paste(extra_dimplot_cols, collapse = ", ")) else "",
      if (run_singler) paste0(", ", singler_col,
                              " (numbered for readability)") else "",
      if (run_sctype)  ", scType automated annotation" else "",
      if (run_rcellmarker) ", rcellmarker automated annotation" else "",
      "."
    ))
  }

  # ── 6. Canonical-marker dotplots ───────────────────────────────────────────
  dotplot_list <- list()
  {
    csv_path <- if (is.null(markers_csv)) {
      system.file("extdata", "Markers.csv", package = "scSidekick")
    } else markers_csv

    if (!nzchar(csv_path) || !file.exists(csv_path)) {
      warning("Markers CSV not found — skipping canonical-marker dotplots.")
    } else {
      mdf_all <- utils::read.csv(csv_path, stringsAsFactors = FALSE)

      # Rename "Cell types" → "Cell.types" (read.csv already does this, but be safe)
      if ("Cell.types" %in% colnames(mdf_all)) {
        mdf_all$CellType <- mdf_all$Cell.types
      } else if ("Cell types" %in% colnames(mdf_all)) {
        mdf_all$CellType <- mdf_all[["Cell types"]]
      }

      if (!gene_col_csv %in% colnames(mdf_all)) {
        warning("Column '", gene_col_csv,
                "' not found in markers CSV — skipping canonical dotplots.")
      } else {
        # Determine groups
        groups_to_plot <- if (!is.null(marker_groups)) {
          marker_groups
        } else {
          # NULL → one group with everything
          list("All canonical markers" = unique(mdf_all$Type))
        }

        for (grp_name in names(groups_to_plot)) {
          type_filter <- groups_to_plot[[grp_name]]
          sub_df <- mdf_all[mdf_all$Type %in% type_filter &
                            !is.na(mdf_all[[gene_col_csv]]) &
                            nzchar(mdf_all[[gene_col_csv]]), ]
          if (nrow(sub_df) == 0) next

          markers_df <- data.frame(
            Genes    = sub_df[[gene_col_csv]],
            CellType = factor(sub_df$CellType,
                              levels = unique(sub_df$CellType)),
            stringsAsFactors = FALSE
          )
          dp <- tryCatch(
            SplitDotPlot(seurat_object, markers_df = markers_df,
                         gene_col = "Genes", gene_group_col = "CellType",
                         group.by = cluster_column, split.by = NULL,
                         col.min = col.min, col.max = col.max,
                         scale.min = scale.min, scale.max = scale.max),
            error = function(e) { warning(conditionMessage(e)); NULL }
          )
          if (!is.null(dp))
            dotplot_list[[grp_name]] <- dp
        }
      }
    }
  }

  # ── 7. Top DE markers dotplot (cluster-faceted) ────────────────────────────
  p_top_markers <- NULL
  {
    # Detect column names from the markers data frame
    fc_col <- if ("logFC" %in% colnames(markers)) "logFC" else
              if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else NULL
    grp_col <- if ("group" %in% colnames(markers)) "group" else
               if ("cluster" %in% colnames(markers)) "cluster" else NULL
    feat_col <- if ("feature" %in% colnames(markers)) "feature" else
                if ("gene" %in% colnames(markers)) "gene" else NULL
    padj_col <- if ("padj" %in% colnames(markers)) "padj" else
                if ("p_val_adj" %in% colnames(markers)) "p_val_adj" else NULL

    if (!is.null(fc_col) && !is.null(grp_col) && !is.null(feat_col)) {
      m2 <- markers
      colnames(m2)[colnames(m2) == fc_col]   <- "logFC"
      colnames(m2)[colnames(m2) == grp_col]  <- "group"
      colnames(m2)[colnames(m2) == feat_col] <- "feature"
      if (!is.null(padj_col))
        colnames(m2)[colnames(m2) == padj_col] <- "padj"

      top_dp <- m2
      if ("padj" %in% colnames(top_dp))
        top_dp <- dplyr::filter(top_dp, padj < 0.05)
      if ("avgExpr" %in% colnames(top_dp))
        top_dp <- top_dp |>
          dplyr::group_by(group) |>
          dplyr::top_n(n = 100, wt = avgExpr) |>
          dplyr::top_n(n = -50, wt = if ("padj" %in% colnames(top_dp)) padj else logFC) |>
          dplyr::ungroup()
      top_dp <- top_dp |>
        dplyr::group_by(group) |>
        dplyr::top_n(n = top_n_per_cluster, wt = logFC) |>
        dplyr::ungroup() |>
        dplyr::arrange(group)
      top_dp <- top_dp[!duplicated(top_dp$feature), ]

      if (nrow(top_dp) > 0) {
        top_df <- data.frame(Genes = top_dp$feature,
                             ClusterGroup = top_dp$group,
                             stringsAsFactors = FALSE)
        p_top_markers <- tryCatch(
          SplitDotPlot(seurat_object, markers_df = top_df,
                       gene_col = "Genes", gene_group_col = "ClusterGroup",
                       group.by = cluster_column, split.by = NULL,
                       col.min = col.min, col.max = col.max,
                       scale.min = scale.min, scale.max = scale.max),
          error = function(e) { warning(conditionMessage(e)); NULL }
        )
        if (!is.null(p_top_markers))
          dotplot_list[["Top DE markers"]] <- p_top_markers
      }
    }
  }

  # ── 8. Combined Assignment Helper PDF ──────────────────────────────────────
  # Stack everything vertically on ONE page using patchwork's `/` operator,
  # matching the original `p1 / p23/p22 / p2 / p3` pattern.
  # Each entry in dotplot_list is already a patchwork object from SplitDotPlot.
  {
    all_plot_list <- c(list(p1), unname(dotplot_list))
    all_plot_list <- Filter(Negate(is.null), all_plot_list)

    if (length(all_plot_list) > 0) {
      f_helper <- file.path(output_dir,
                            paste0("Assignment Helper ", file_prefix,
                                   " UMAPs and Dotplots.pdf"))

      # Combine into a single stacked figure
      combined <- Reduce(`/`, all_plot_list)

      # Height: ~7 in for the UMAP panel + ~8 in per dotplot panel
      total_h <- max(7, min(80, 7 + length(dotplot_list) * 8))

      grDevices::pdf(f_helper, width = 25, height = total_h)
      print(combined)
      grDevices::dev.off()
      message("  Saved: ", f_helper)

      .write_legend_sidecar(f_helper, paste0(
        "Combined annotation helper for ", object_name, subset_name, ". ",
        "Top panel: multi-method UMAP comparison (cluster labels, ",
        if (run_singler) paste0(singler_col, ", ") else "",
        if (run_sctype)  "scType, "      else "",
        if (run_rcellmarker) "rcellmarker, " else "",
        "and any additional metadata). ",
        "Subsequent panels: canonical marker dot plots",
        if ("Top DE markers" %in% names(dotplot_list))
          paste0(" and top-", top_n_per_cluster, " DE genes per cluster")
        else "",
        ". Dot size = % expressing cells; colour = mean scaled expression."
      ))
    }
  }

  seurat_object@misc$annotation_colors <- cols_out
  message("CellTypeAssignmentHelper complete.")
  invisible(seurat_object)
}
