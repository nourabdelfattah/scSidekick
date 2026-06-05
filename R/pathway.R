# =============================================================================
# scSidekick GSEA pathway analysis wrapper

# ---------------------------------------------------------------------------
# .gsea_ht()
# Internal helper: build one ComplexHeatmap NES heatmap with auto-sized PDF.
#
# Design choices:
#   - Diverging colour (blue-white-red) symmetric around 0 — NES sign matters
#   - Cell height and font shrink gracefully for large matrices
#   - PDF dimensions computed from content so pathway names are never clipped
#   - heatmap_params (named list) are merged via modifyList so the caller can
#     override any default without repeating all arguments
# ---------------------------------------------------------------------------
.gsea_ht <- function(mat, title, filepath, heatmap_params = list()) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("Package 'ComplexHeatmap' is required for GSEA heatmaps.")
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Package 'circlize' is required for GSEA heatmaps.")

  n_rows <- nrow(mat)
  n_cols <- ncol(mat)

  # Diverging colour centred on 0 — matches the biological meaning of NES
  nes_lim  <- max(1, min(4, max(abs(mat), na.rm = TRUE)))
  col_fun  <- circlize::colorRamp2(
    c(-nes_lim, 0, nes_lim),
    c("#2166ac", "white", "#b2182b")
  )

  # Adaptive cell size (pt): shrink for large matrices, cap for small ones
  cell_h_pt <- max(7,  min(14, 400 / max(n_rows, 1)))
  cell_w_pt <- max(12, min(25, 200 / max(n_cols, 1)))
  rn_fs     <- max(5,  min(9,  cell_h_pt * 0.75))   # row-name font size

  # PDF dimensions (inches): body + row names margin + legend + title padding
  rn_max_in <- max(nchar(rownames(mat)), na.rm = TRUE) * rn_fs * 0.50 / 72
  pdf_h <- min(40, max(3.5, n_rows * cell_h_pt / 72 + 1.5))
  pdf_w <- min(40, max(4.5, n_cols * cell_w_pt / 72 + rn_max_in + 2.5))

  default_args <- list(
    mat,
    name              = "NES",
    col               = col_fun,
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = grid::gpar(fontsize = rn_fs),
    column_names_gp   = grid::gpar(fontsize = 9),
    column_names_rot  = 45,
    column_title      = title,
    column_title_gp   = grid::gpar(fontsize = 10, fontface = "bold"),
    border            = FALSE,
    use_raster        = n_rows > 100 || n_cols > 50,
    width             = grid::unit(n_cols * cell_w_pt, "pt"),
    height            = grid::unit(n_rows * cell_h_pt, "pt"),
    heatmap_legend_param = list(
      title         = "NES",
      at            = c(-nes_lim, 0, nes_lim),
      labels        = c(sprintf("%.1f", -nes_lim), "0", sprintf("%.1f", nes_lim)),
      legend_height = grid::unit(30, "mm"),
      title_gp      = grid::gpar(fontsize = 9),
      labels_gp     = grid::gpar(fontsize = 8)
    )
  )

  ht_args <- utils::modifyList(default_args, heatmap_params)
  ht      <- do.call(ComplexHeatmap::Heatmap, ht_args)

  grDevices::pdf(filepath, width = pdf_w, height = pdf_h)
  ComplexHeatmap::draw(ht, padding = grid::unit(c(5, 5, 5, 5), "mm"))
  grDevices::dev.off()

  invisible(ht)
}
#
# RunGSEA — runs Wilcoxon rank-sum DE (via presto::wilcoxauc), then fgsea
#   against one or more MSigDB collections, and produces per-cluster lollipop
#   plots and NES heatmaps. Returns a nested results list and optionally writes
#   CSVs + PDFs to disk.
#
# Replaces the ad-hoc setwd()-based scripts with a clean function interface.
# =============================================================================

#' Run GSEA pathway analysis on a Seurat object
#'
#' Performs the following steps for each combination of `label.by` (outer
#' loop) and `split.by` (inner loop) levels:
#' \enumerate{
#'   \item Subset the object to the current level.
#'   \item Run Wilcoxon rank-sum DE with [presto::wilcoxauc()].
#'   \item For each pathway database in `pathway_sets`, run
#'         [fgsea::fgsea()] using AUC scores as ranking statistics.
#'   \item Generate a lollipop NES plot (top/bottom 10 pathways per cluster).
#'   \item Build NES heatmaps showing top `top_n` pathways per cluster.
#' }
#'
#' **No `setwd()` calls are made.** All output paths are constructed from
#' `output_dir`.
#'
#' @param seurat_object A Seurat object.
#' @param group.by Character. Metadata column whose levels are compared within
#'   each subset (used as DE groups). Default `"Assignment"`.
#' @param split.by Character or `NULL`. Metadata column whose levels each get
#'   their own independent GSEA (e.g., run GSEA separately per broad cell type
#'   such as `"GlobalAssignment"`). `NULL` runs a single GSEA on the whole
#'   object (or on each `label.by` subset), comparing all `group.by` levels
#'   without any further cell-type subsetting. Default `"GlobalAssignment"`.
#' @param label.by Character or `NULL`. Outer metadata column for further
#'   subsetting (e.g., condition). `NULL` treats the whole object as one group.
#' @param pathway_sets Named list of MSigDB gene-set databases to test. Each
#'   element is a list with a `category` field and an optional `subcategory`
#'   field; the element name becomes the short label used in output filenames
#'   and the return list. The default covers four curated collections:
#'   ```r
#'   list(
#'     Hallmark = list(category = "H"),
#'     KEGG     = list(category = "C2", subcategory = "CP:KEGG"),
#'     Reactome = list(category = "C2", subcategory = "CP:REACTOME"),
#'     WP       = list(category = "C2", subcategory = "CP:WIKIPATHWAYS")
#'   )
#'   ```
#'   Other useful collections (pass as additional list elements):
#'   \itemize{
#'     \item `list(category = "C2", subcategory = "CP:BIOCARTA")` — BioCarta
#'     \item `list(category = "C2", subcategory = "CP:PID")` — NCI Pathway
#'       Interaction Database
#'     \item `list(category = "C5", subcategory = "GO:BP")` — Gene Ontology
#'       Biological Process
#'     \item `list(category = "C5", subcategory = "GO:MF")` — GO Molecular
#'       Function
#'     \item `list(category = "C5", subcategory = "GO:CC")` — GO Cellular
#'       Component
#'     \item `list(category = "C7", subcategory = "IMMUNESIGDB")` — ImmuneSigDB
#'       (immune gene sets; recommended for immune datasets)
#'     \item `list(category = "C8")` — Cell-type signature gene sets
#'   }
#'   See `msigdbr::msigdbr_collections()` for the full catalogue.
#' @param species Character. Species passed to [msigdbr::msigdbr()]. Must
#'   match the species names recognised by that function. Common values:
#'   \itemize{
#'     \item `"Homo sapiens"` — human (default)
#'     \item `"Mus musculus"` — mouse
#'     \item `"Rattus norvegicus"` — rat
#'     \item `"Danio rerio"` — zebrafish
#'   }
#'   Run `msigdbr::msigdbr_species()` for the complete list.
#' @param assay Character. Seurat assay for DE via `presto::wilcoxauc()`.
#'   Default `"RNA"`. For BPCells sketch workflows use `"sketch"`.
#' @param min_cells Integer. Skip subsets with fewer cells. Default `20`.
#'   A separate guard also skips any `(label, split)` combination that has
#'   fewer than 2 levels of `group.by` after subsetting — this prevents
#'   `presto::wilcoxauc` from crashing when e.g. one `label.by` level contains
#'   only one `group.by` category (e.g., a Cognitive.Status group with a single
#'   Sex). A message is emitted for each skipped combination.
#' @param padj_thresh Numeric. Adjusted p-value threshold for DE and fgsea
#'   filtering. Default `0.05`.
#' @param logfc_thresh Numeric. Minimum log-fold-change for DE filtering.
#'   Default `0.1`.
#' @param top_n Integer. Number of top and bottom pathways per cluster shown
#'   in summary heatmaps. Default `5`.
#' @param output_dir Character or `NULL`. If provided, CSVs and PDFs are
#'   written here in a structured subdirectory tree. If `NULL`, results are
#'   only returned as an R list.
#' @param heatmap_params Named list of additional arguments forwarded to
#'   [ComplexHeatmap::Heatmap()] for both the full and summary NES heatmaps.
#'   Any default set by scSidekick can be overridden here — e.g.
#'   ```r
#'   heatmap_params = list(
#'     clustering_distance_rows = "pearson",
#'     row_names_gp = grid::gpar(fontsize = 6, fontface = "italic"),
#'     col = circlize::colorRamp2(c(-3, 0, 3),
#'                                c("navy", "white", "darkred"))
#'   )
#'   ```
#' @param resume Logical. If `TRUE` and `output_dir` is set, skip work that
#'   has already been completed:
#'   \itemize{
#'     \item **DE cache**: Wilcoxon DE results for each `(label, split)` pair
#'       are saved to `output_dir/label/de_cache_<split>_<label>.rds` and
#'       reloaded on subsequent runs instead of re-running `presto::wilcoxauc`.
#'     \item **DB-level skip**: If the full NES heatmap PDF for a
#'       `(label, split, db_name)` triple already exists, that database is
#'       skipped entirely (its slot in the return list is set to `NULL`).
#'   }
#'   Default `FALSE`.
#'
#' @return A nested named list structured as
#'   `results[[label]][[split]][[db_name]]` with elements:
#'   \describe{
#'     \item{`de_table`}{Full DE results from `presto::wilcoxauc`.}
#'     \item{`gsea_table`}{Full fgsea NES table (all pathways, all clusters).}
#'     \item{`nes_matrix`}{Wide NES matrix (pathways × clusters).}
#'     \item{`lollipop_plots`}{Named list of ggplot2 lollipop plots per cluster.}
#'   }
#'
#' @export
RunGSEA <- function(seurat_object,
                     group.by       = "Assignment",
                     split.by       = "GlobalAssignment",
                     label.by       = NULL,
                     pathway_sets   = list(
                       Hallmark  = list(category = "H"),
                       KEGG      = list(category = "C2", subcategory = "CP:KEGG"),
                       Reactome  = list(category = "C2", subcategory = "CP:REACTOME"),
                       WP        = list(category = "C2", subcategory = "CP:WIKIPATHWAYS")
                     ),
                     species        = "Homo sapiens",
                     assay          = "RNA",
                     min_cells      = 20L,
                     padj_thresh    = 0.05,
                     logfc_thresh   = 0.1,
                     top_n          = 5L,
                     output_dir     = NULL,
                     resume         = FALSE,
                     heatmap_params = list(row_names_side      = "left",
                                           show_row_dend       = FALSE,
                                           row_names_max_width = grid::unit(15, "cm"))) {

  # Validate required packages
  for (pkg in c("presto", "fgsea", "msigdbr", "ComplexHeatmap", "circlize")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. Install it with install.packages('", pkg, "').")
  }

  # Pre-load all fgsea gene sets so we only query msigdbr once per database
  message("Loading pathway gene sets from MSigDB...")
  fgsea_dbs <- lapply(names(pathway_sets), function(db_name) {
    ps  <- pathway_sets[[db_name]]
    mdf <- if (!is.null(ps$subcategory)) {
      msigdbr::msigdbr(species = species, category = ps$category,
                       subcategory = ps$subcategory)
    } else {
      msigdbr::msigdbr(species = species, category = ps$category)
    }
    split(mdf$gene_symbol, mdf$gs_name)
  })
  names(fgsea_dbs) <- names(pathway_sets)

  # Determine outer loop levels
  if (!is.null(label.by)) {
    label_levels <- levels(factor(seurat_object@meta.data[[label.by]]))
  } else {
    label_levels <- "All"
  }

  # ── Write method params to output_dir so create_analysis_pptx() finds them ─
  # Inherits shared fields (n_cells, resolution, dataset …) from the parent
  # analysis_params.json via .write_subdir_params(); adds GSEA-specific fields.
  if (!is.null(output_dir)) {
    db_list_str <- paste(names(pathway_sets), collapse = ", ")
    .write_subdir_params(output_dir, list(
      date              = format(Sys.Date()),
      gsea_group_by     = group.by,
      gsea_split_by     = split.by,
      gsea_label_by     = label.by,
      gsea_databases    = names(pathway_sets),
      gsea_padj_thresh  = padj_thresh,
      gsea_logfc_thresh = logfc_thresh,
      gsea_top_n        = top_n,
      gsea_species      = species,
      methods_text      = paste0(
        "Gene set enrichment analysis (GSEA) was performed using fgsea ",
        "(Korotkevich et al., 2021) with Wilcoxon rank-sum AUC pre-ranking ",
        "via presto. Differential expression was computed comparing ",
        group.by, " groups",
        if (!is.null(split.by))
          paste0(" within each ", split.by, " cell-type subset")
        else "",
        if (!is.null(label.by))
          paste0(", independently for each level of ", label.by)
        else "",
        ". The following MSigDB gene-set collections were tested: ",
        db_list_str, ". ",
        "Gene sets with fewer than 10 members were excluded. ",
        "Significance threshold: adjusted p < ", padj_thresh,
        "; minimum log-fold-change for DE pre-ranking: ", logfc_thresh, ". ",
        "Lollipop plots show the top/bottom 10 significant pathways per ",
        "comparison group ranked by Normalised Enrichment Score (NES); ",
        "summary heatmaps show the top & bottom ", top_n,
        " pathways per group."
      )
    ))
  }

  results <- list()

  for (lab in label_levels) {
    message("\n===== Label: ", lab, " =====")

    if (!is.null(label.by)) {
      Seurat::Idents(seurat_object) <- seurat_object@meta.data[[label.by]]
      obj_sub <- subset(seurat_object, idents = lab)
    } else {
      obj_sub <- seurat_object
    }

    # split.by = NULL means "no subsetting — run GSEA on all cells at once".
    # A single dummy level "All" is used so the loop body executes exactly once.
    use_split  <- !is.null(split.by) && split.by %in% colnames(obj_sub@meta.data)
    split_levels <- if (use_split)
      levels(factor(obj_sub@meta.data[[split.by]]))
    else
      "All"
    results[[lab]] <- list()

    for (sp in split_levels) {
      message("  Split: ", sp)

      if (use_split) {
        Seurat::Idents(obj_sub) <- obj_sub@meta.data[[split.by]]
        cells_sp <- Seurat::WhichCells(obj_sub, idents = sp)
      } else {
        cells_sp <- colnames(obj_sub)   # all cells — no subsetting
      }

      if (length(cells_sp) < min_cells) {
        message("    Skipping — fewer than ", min_cells, " cells.")
        next
      }

      obj_sp   <- subset(obj_sub, cells = cells_sp)
      Seurat::Idents(obj_sp) <- obj_sp@meta.data[[group.by]]
      clusters <- levels(factor(obj_sp@meta.data[[group.by]]))
      sp_safe  <- gsub("/", ".", sp)          # safe for filenames

      # ── Guard: need ≥ 2 group.by levels for DE + fgsea to run ─────────────
      # This fires when label.by subsetting leaves only one level of group.by
      # in the current split — e.g., one Cognitive.Status group has only one
      # Sex, so wilcoxauc cannot compare anything.
      if (length(clusters) < 2) {
        message("    Skipping ", sp, " [", lab, "] — only ",
                length(clusters), " level(s) of '", group.by, "' (",
                paste(clusters, collapse = ", "),
                "); need ≥ 2 to compare.")
        next
      }

      # ── DE via presto (cached to disk when output_dir is set) ─────────────
      # On resume = TRUE, reload the cached RDS instead of re-running wilcoxauc.
      de_cache_path <- if (!is.null(output_dir))
        file.path(output_dir, lab,
                  paste0("de_cache_", sp_safe, "_", lab, ".rds"))
      else NULL

      if (resume && !is.null(de_cache_path) && file.exists(de_cache_path)) {
        message("    Loading cached DE: ", basename(de_cache_path))
        de_all <- readRDS(de_cache_path)
      } else {
        de_all <- presto::wilcoxauc(obj_sp, group_by = group.by,
                                    seurat_assay = assay)
        if (!is.null(de_cache_path)) {
          dir.create(dirname(de_cache_path), recursive = TRUE,
                     showWarnings = FALSE)
          saveRDS(de_all, de_cache_path)
        }
      }

      results[[lab]][[sp]] <- list()

      for (db_name in names(fgsea_dbs)) {
        message("    DB: ", db_name)

        # ── Resume: skip this (lab, sp, db_name) if full heatmap exists ─────
        if (resume && !is.null(output_dir)) {
          full_hm_check <- file.path(
            output_dir, lab, db_name,
            paste0("GSEA ", db_name, " ", sp_safe, " ", lab,
                   " (NES)-Heatmap.pdf")
          )
          if (file.exists(full_hm_check)) {
            message("      Resuming: heatmap exists — skipping ", db_name)
            results[[lab]][[sp]][[db_name]] <- list(
              de_table       = de_all,
              nes_matrix     = NULL,   # not reconstructed from disk
              lollipop_plots = list()
            )
            next
          }
        }
        fgsea_sets <- fgsea_dbs[[db_name]]

        # db_dir is deterministic — hoist it so the cluster-level resume check
        # can reference it before the directory is necessarily created.
        db_dir <- if (!is.null(output_dir))
          file.path(output_dir, lab, db_name) else NULL

        # NES accumulator (pathways × clusters)
        all_pathways <- names(fgsea_sets)
        nes_mat      <- data.frame(pathway = all_pathways,
                                   stringsAsFactors = FALSE)
        lollipop_plots <- list()

        for (cl in clusters) {
          cl_safe <- gsub("/", ".", cl)

          # ── Resume: skip fgsea if per-cluster CSV already exists ───────────
          # The CSV is the finest-grained checkpoint — it's written immediately
          # after fgsea() completes.  Loading NES from it is fast and avoids
          # re-running the expensive fgsea call.
          # Note: lollipop plot objects are NOT restored from disk; the PDFs
          # already exist and the in-memory list is left empty for that cluster.
          if (resume && !is.null(db_dir)) {
            cl_csv_check <- file.path(db_dir,
                                      paste0(db_name, " ", cl_safe, " ",
                                             sp_safe, " ", lab, ".csv"))
            if (file.exists(cl_csv_check)) {
              message("        Cluster ", cl, ": CSV cached — skipping fgsea")
              cached <- tryCatch(
                utils::read.csv(cl_csv_check, row.names = 1,
                                stringsAsFactors = FALSE),
                error = function(e) NULL
              )
              if (!is.null(cached) &&
                  all(c("pathway", "NES") %in% colnames(cached))) {
                nes_col_df <- data.frame(
                  pathway = cached$pathway,
                  NES     = suppressWarnings(as.numeric(cached$NES)),
                  stringsAsFactors = FALSE
                )
                colnames(nes_col_df)[2] <- cl
                nes_mat <- dplyr::left_join(nes_mat, nes_col_df,
                                            by = "pathway")
              }
              next
            }
          }

          ranks <- de_all |>
            dplyr::filter(group == cl) |>
            dplyr::arrange(dplyr::desc(auc)) |>
            dplyr::select(feature, auc) |>
            tibble::deframe()

          gsea_res <- fgsea::fgsea(fgsea_sets, stats = ranks,
                                   minSize = 10)

          gsea_tidy <- gsea_res |>
            tibble::as_tibble() |>
            dplyr::arrange(dplyr::desc(NES)) |>
            dplyr::mutate(Enrichment = ifelse(NES > 0, "Up-regulated",
                                              "Down-regulated"))

          # Write per-cluster CSV (creates db_dir if needed)
          if (!is.null(db_dir)) {
            dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)
            cl_csv <- file.path(db_dir,
                                paste0(db_name, " ", cl_safe, " ", sp_safe,
                                       " ", lab, ".csv"))
            utils::write.csv(apply(gsea_tidy, 2, as.character), cl_csv)
          }

          # Lollipop plot (top/bottom 10 significant)
          sig        <- dplyr::filter(gsea_tidy, padj < padj_thresh)
          filt_res   <- dplyr::bind_rows(
            dplyr::slice_head(dplyr::filter(sig, Enrichment == "Up-regulated"),   n = 10),
            dplyr::slice_head(dplyr::filter(sig, Enrichment == "Down-regulated"), n = 10)
          )

          if (nrow(filt_res) > 0) {
            lp <- ggplot2::ggplot(filt_res,
                                  ggplot2::aes(reorder(pathway, NES), NES)) +
              ggplot2::geom_segment(
                ggplot2::aes(reorder(pathway, NES), xend = pathway,
                             y = 0, yend = NES)
              ) +
              ggplot2::geom_point(
                size = 3,
                ggplot2::aes(fill = Enrichment),
                shape = 21, stroke = 1
              ) +
              ggplot2::scale_fill_manual(
                values = c("Down-regulated" = "dodgerblue",
                           "Up-regulated"   = "firebrick")
              ) +
              ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = 0.2)
              ) +
              ggplot2::coord_flip() +
              ggplot2::labs(x = "Pathway", y = "NES",
                            title = paste(db_name, "-", cl, sp)) +
              theme_NourMin()

            lollipop_plots[[cl]] <- lp

            if (!is.null(db_dir)) {
              lp_path <- file.path(db_dir,
                                   paste0(db_name, " - ", lab, " ",
                                          cl_safe, " ", sp_safe, ".pdf"))
              ggplot2::ggsave(plot = lp, filename = lp_path,
                              width = 10, height = 10)
              .write_legend_sidecar(lp_path, paste0(
                "Lollipop plot showing the top up- and down-regulated gene sets ",
                "from the ", db_name, " pathway database in ", cl, " cells",
                if (nchar(sp) > 0) paste0(" (", sp, ")") else "",
                if (nchar(lab) > 0 && lab != "all") paste0(" from ", lab) else "",
                ", ranked by Normalised Enrichment Score (NES). ",
                "Gene set enrichment analysis was performed using the fgsea ",
                "algorithm with Wilcoxon rank-sum pre-ranking (presto). ",
                "Positive NES (red) indicates enrichment relative to all other ",
                "clusters; negative NES (blue) indicates depletion. ",
                "Only gene sets with adjusted p-value < ", padj_thresh, " are shown."
              ))
            }
          }

          # Accumulate NES column
          nes_col <- gsea_tidy |>
            dplyr::select(pathway, NES) |>
            dplyr::rename(!!cl := NES)
          nes_mat <- dplyr::left_join(nes_mat, nes_col, by = "pathway")
        }

        # Build NES matrix and heatmaps.
        # Convert to a base R *matrix* (not tibble/data.frame) so that
        # nes_mat[, ci] always returns a *named* vector — dplyr::left_join
        # can silently return a tibble that drops row names, which makes
        # names(sort(col_vec)) return NULL and breaks the top_n selection.
        pathway_names     <- nes_mat$pathway
        nes_mat           <- as.matrix(
          nes_mat[, setdiff(colnames(nes_mat), "pathway"), drop = FALSE]
        )
        rownames(nes_mat) <- pathway_names
        # Reorder columns to match cluster factor ordering
        col_order         <- clusters[clusters %in% colnames(nes_mat)]
        nes_mat           <- nes_mat[, col_order, drop = FALSE]
        nes_mat[is.na(nes_mat)] <- 0

        if (!is.null(db_dir) && ncol(nes_mat) > 0) {

          # Full heatmap — ComplexHeatmap with auto-sized PDF so pathway
          # names are never clipped and the matrix is never squished.
          full_hm_path <- file.path(db_dir,
                                    paste0("GSEA ", db_name, " ", sp_safe,
                                           " ", lab, " (NES)-Heatmap.pdf"))
          .gsea_ht(nes_mat,
                   title       = paste(db_name, sp_safe, lab, "(NES)"),
                   filepath    = full_hm_path,
                   heatmap_params = heatmap_params)
          .write_legend_sidecar(full_hm_path, paste0(
            "Heatmap of Normalised Enrichment Scores (NES) from GSEA using the ",
            db_name, " pathway database, showing all tested pathways across ",
            sp_safe, " cells",
            if (lab != "All") paste0(" (", lab, ")") else "", ". ",
            "Pathways (rows) are hierarchically clustered; columns represent ",
            "comparison groups (", group.by, "). Blue = depleted (negative NES); ",
            "red = enriched (positive NES). Colour scale is symmetric around 0 ",
            "and calibrated to the maximum |NES| in this matrix."
          ))

          # Summary heatmap — top/bottom top_n per cluster.
          summary_rows <- unique(unlist(lapply(seq_len(ncol(nes_mat)), function(ci) {
            col_vec <- nes_mat[, ci]
            c(names(utils::head(sort(col_vec, decreasing = TRUE), top_n)),
              names(utils::head(sort(col_vec),                     top_n)))
          })))
          if (length(summary_rows) > 0) {
            summary_rows    <- summary_rows[summary_rows %in% rownames(nes_mat)]
            sub_mat         <- nes_mat[summary_rows, , drop = FALSE]
            sub_mat[is.na(sub_mat)] <- 0
            summary_hm_path <- file.path(db_dir,
                                         paste0("GSEA ", db_name, " ", sp_safe,
                                                " ", lab, " top", top_n,
                                                "-Heatmap.pdf"))
            .gsea_ht(sub_mat,
                     title       = paste(db_name, sp_safe, lab,
                                         paste0("top & bottom ", top_n), "(NES)"),
                     filepath    = summary_hm_path,
                     heatmap_params = heatmap_params)
            .write_legend_sidecar(summary_hm_path, paste0(
              "Summary NES heatmap for the ", db_name, " pathway database in ",
              sp_safe, " cells",
              if (lab != "All") paste0(" (", lab, ")") else "", ". ",
              "Rows show the top ", top_n, " up-regulated and bottom ", top_n,
              " down-regulated pathways per comparison group, selected by NES rank. ",
              "Pathways are hierarchically clustered. Blue = depleted; red = enriched."
            ))
          }
        }

        results[[lab]][[sp]][[db_name]] <- list(
          de_table       = de_all,
          nes_matrix     = nes_mat,
          lollipop_plots = lollipop_plots
        )
      }
    }
  }

  message("\nRunGSEA complete.")
  results
}
