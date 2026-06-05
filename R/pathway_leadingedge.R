# =============================================================================
# scSidekick — Leading-Edge Visualization  (pathway_leadingedge.R)
#
# Exported:
#   StackedVlnPlot()        — stacked per-feature violin plots (patchwork)
#   VisualizeLeadingEdge()  — loads GSEA CSVs, extracts leading-edge genes,
#                             heatmap + optional DEG bonus plots
#
# Internal helpers:
#   .parse_leading_edge()   — parse serialized leadingEdge column from fgsea CSV
# =============================================================================


# -----------------------------------------------------------------------------
# .parse_leading_edge()
# Parses a single entry from a leadingEdge column written by fgsea to CSV.
# Handles: c("GENE1", "GENE2") format, comma-sep, or whitespace-sep strings.
# -----------------------------------------------------------------------------
.parse_leading_edge <- function(x) {
  if (is.na(x) || nchar(trimws(x)) == 0) return(character(0))
  x <- trimws(x)
  # c("A", "B", ...) format produced when fgsea list-column is saved with write.csv
  if (grepl("^c\\(", x)) {
    tryCatch(eval(parse(text = x)), error = function(e) character(0))
  } else if (grepl(",", x)) {
    trimws(strsplit(x, ",")[[1]])
  } else {
    trimws(strsplit(x, "\\s+")[[1]])
  }
}


# =============================================================================
# StackedVlnPlot
# =============================================================================

#' Stacked Violin Plot for Multiple Features
#'
#' @description
#' Builds a stacked vertical layout of per-feature violin plots using
#' \code{Seurat::VlnPlot}, assembled with \code{patchwork}. X-axis labels and
#' ticks are hidden on all panels except the bottom one; the gene name is shown
#' as the y-axis title (rotated 0\eqn{^\circ}) on each panel.
#'
#' @param seurat_object A Seurat object.
#' @param features Character vector of features (genes) to plot. One panel per
#'   feature.
#' @param pt.size Jitter point size (default \code{0} — no points).
#' @param plot.margin Plot margin around each panel as a
#'   \code{grid::unit} object (default \code{unit(c(-0.75, 0, -0.75, 0), "cm")}).
#'   Negative top/bottom values collapse space between stacked panels.
#' @param ... Additional arguments forwarded to \code{Seurat::VlnPlot} (e.g.
#'   \code{group.by}, \code{split.by}, \code{cols}).
#'
#' @return A \code{patchwork} object (one panel per feature, stacked in a
#'   single column).
#'
#' @export
StackedVlnPlot <- function(
    seurat_object,
    features,
    pt.size     = 0,
    plot.margin = grid::unit(c(-0.75, 0, -0.75, 0), "cm"),
    ...
) {

  .vln_theme <- ggplot2::theme(
    legend.position  = "none",
    axis.text.x      = ggplot2::element_blank(),
    axis.ticks.x     = ggplot2::element_blank(),
    axis.title.x     = ggplot2::element_blank(),
    axis.title.y     = ggplot2::element_text(
      size   = ggplot2::rel(1), angle = 0,
      face   = "bold", colour = "black",
      vjust  = 0.5
    ),
    axis.text.y      = ggplot2::element_text(size = ggplot2::rel(0.85),
                                              colour = "black"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line        = ggplot2::element_line(colour = "black", linewidth = 0.4),
    strip.background = ggplot2::element_rect(fill = NA, colour = NA),
    strip.text       = ggplot2::element_text(colour = "black",
                                              size = ggplot2::rel(1.1)),
    plot.margin      = plot.margin
  )

  .modify_vln <- function(feature) {
    Seurat::VlnPlot(seurat_object, features = feature, pt.size = pt.size, ...) +
      ggplot2::xlab("") +
      ggplot2::ylab(feature) +
      ggplot2::ggtitle("") +
      .vln_theme
  }

  .extract_max <- function(p) {
    ymax <- max(ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    ceiling(ymax)
  }

  plot_list <- lapply(features, .modify_vln)

  # Restore x-axis on the bottom panel only
  n <- length(plot_list)
  plot_list[[n]] <- plot_list[[n]] +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(
        colour = "black", angle = 25, hjust = 1,
        size   = ggplot2::rel(0.9), face = "bold"
      ),
      axis.ticks.x = ggplot2::element_line()
    )

  # Expand y-limit to max observed value across all panels
  ymaxs     <- vapply(plot_list, .extract_max, numeric(1))
  plot_list <- mapply(
    function(p, y) p + ggplot2::expand_limits(y = y),
    plot_list, ymaxs, SIMPLIFY = FALSE
  )

  patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
}


# =============================================================================
# VisualizeLeadingEdge
# =============================================================================

#' Visualize Leading-Edge Genes from GSEA Results
#'
#' @description
#' Loads one or more GSEA result CSV files (from \code{\link{RunGSEA}} or
#' \code{\link{RunGSEA_pseudobulk}}), searches for pathways by keyword or exact
#' name, extracts the union of leading-edge genes, and produces a per-cell
#' expression heatmap. Optionally, if a DEG data frame is supplied, the
#' leading-edge genes are intersected with significant DEGs and additional
#' diagnostic plots are produced (a second heatmap, dotplot, stacked violins,
#' and boxplots).
#'
#' @section Pathway search:
#' \describe{
#'   \item{\code{search_terms} as a character vector}{OR logic — any term
#'     matches, case-insensitive.}
#'   \item{\code{search_terms} as a list of character vectors}{Each element is
#'     an AND-group (all terms must appear in the same pathway name); results are
#'     OR'd across elements.}
#'   \item{\code{pathways}}{Exact pathway names (alternative to
#'     \code{search_terms}).}
#' }
#'
#' @param seurat_object A Seurat object.
#' @param gsea_files Character vector of CSV file paths (output of
#'   \code{RunGSEA} or \code{RunGSEA_pseudobulk}). Must each contain
#'   \code{pathway} and \code{leadingEdge} columns.
#' @param search_terms Keyword(s) to filter pathways. See \emph{Pathway search}
#'   above.
#' @param pathways Exact pathway names to select (alternative to
#'   \code{search_terms}).
#' @param group.by Metadata column for the primary grouping (required). Used for
#'   column annotation, dotplot, and violin grouping.
#' @param split.by Optional second metadata column for column annotation and
#'   violin splitting.
#' @param assay Assay to pull expression from (default \code{"RNA"}).
#' @param layer Layer to pull expression from (default \code{"data"}).
#' @param min_expr Genes whose row-mean expression is below this threshold are
#'   excluded from the heatmap (default \code{0.05}).
#' @param subset_cells Logical; if \code{TRUE}, subset \code{seurat_object} to
#'   \code{subset_values} in \code{subset_by} before analysis.
#' @param subset_by Metadata column used for subsetting (requires
#'   \code{subset_cells = TRUE}).
#' @param subset_values Character vector of values to retain in \code{subset_by}.
#' @param deg_df Optional DEG data frame (e.g. from presto or
#'   \code{FindMarkers}). When provided, leading-edge genes are intersected with
#'   significant DEGs and additional plots are generated.
#' @param deg_gene_column Column in \code{deg_df} with gene names (default
#'   \code{"feature"}).
#' @param deg_group_column Column in \code{deg_df} with group labels (default
#'   \code{"group"}).
#' @param deg_fc_column Column in \code{deg_df} with log-fold change (default
#'   \code{"logFC"}).
#' @param deg_padj_column Column in \code{deg_df} with adjusted p-value (default
#'   \code{"padj"}).
#' @param deg_padj_cutoff Adjusted p-value cut-off (default \code{0.05}).
#' @param deg_top_n Top N DEGs per group (by logFC) for bonus plots (default
#'   \code{20}).
#' @param add_dotplot Logical; produce a \code{Seurat::DotPlot} of top DEGs
#'   (default \code{TRUE}; requires \code{deg_df}).
#' @param add_violin Logical; produce stacked violin plots of top DEGs (default
#'   \code{TRUE}; requires \code{deg_df}).
#' @param add_boxplots Logical; produce boxplots of top DEGs (default
#'   \code{TRUE}; requires \code{deg_df}).
#' @param add_pvalues Logical; overlay Kruskal-Wallis p-value on boxplots via
#'   \code{ggpubr::stat_compare_means} (default \code{TRUE}).
#' @param add_feature_maps Logical; produce per-gene \code{FeaturePlot} PDFs
#'   for the top DEG genes (default \code{FALSE}; can be slow for large gene
#'   sets).
#' @param group_colors Named character vector of colors for \code{group.by}
#'   levels. Auto-assigned from \code{Nour_pal("all")} if \code{NULL}.
#' @param split_colors Named character vector of colors for \code{split.by}
#'   levels. Auto-assigned from \code{Nour_pal("spectrum")} if \code{NULL}.
#' @param heatmap_params Named list of additional arguments forwarded to
#'   \code{ComplexHeatmap::Heatmap()}. Any key supplied here overrides the
#'   default.
#' @param output_dir Directory where all output files are saved.
#' @param object_name Label used in output file names and subfolder (default
#'   \code{"Analysis"}).
#' @param subset_name Optional label appended to output file names (default
#'   \code{""}).
#'
#' @return Invisibly returns a named list:
#' \describe{
#'   \item{\code{genes_all}}{All leading-edge genes present in the Seurat object
#'     (after \code{min_expr} filtering).}
#'   \item{\code{genes_top}}{Top DEG genes within the leading edge, or
#'     \code{NULL} if no \code{deg_df} was provided.}
#'   \item{\code{top_degs}}{Filtered DEG data frame, or \code{NULL}.}
#' }
#'
#' @seealso \code{\link{RunGSEA}}, \code{\link{RunGSEA_pseudobulk}},
#'   \code{\link{StackedVlnPlot}}
#'
#' @export
VisualizeLeadingEdge <- function(
    seurat_object,
    gsea_files,
    search_terms    = NULL,
    pathways        = NULL,
    group.by,
    split.by        = NULL,
    assay           = "RNA",
    layer           = "data",
    min_expr        = 0.05,
    subset_cells      = FALSE,
    subset_by       = NULL,
    subset_values   = NULL,
    deg_df          = NULL,
    deg_gene_column    = "feature",
    deg_group_column   = "group",
    deg_fc_column      = "logFC",
    deg_padj_column    = "padj",
    deg_padj_cutoff = 0.05,
    deg_top_n       = 20,
    add_dotplot     = TRUE,
    add_violin      = TRUE,
    add_boxplots    = TRUE,
    add_pvalues     = TRUE,
    add_feature_maps = FALSE,
    group_colors    = NULL,
    split_colors    = NULL,
    heatmap_params  = list(
      row_names_side      = "left",
      show_row_dend       = FALSE,
      row_names_max_width = grid::unit(15, "cm")
    ),
    output_dir,
    object_name    = "Analysis",
    subset_name = ""
) {

  # ── 0. Validate inputs ─────────────────────────────────────────────────────
  if (missing(gsea_files) || length(gsea_files) == 0)
    stop("Provide at least one CSV path in 'gsea_files'.")
  missing_files <- gsea_files[!file.exists(gsea_files)]
  if (length(missing_files) > 0)
    stop("File(s) not found:\n  ", paste(missing_files, collapse = "\n  "))
  if (is.null(search_terms) && is.null(pathways))
    stop("Provide 'search_terms' (keyword search) or 'pathways' (exact names).")
  if (missing(output_dir))
    stop("'output_dir' is required.")

  # ── 1. Subset object ───────────────────────────────────────────────────────
  if (isTRUE(subset_cells)) {
    if (is.null(subset_by) || is.null(subset_values))
      stop("'subset_by' and 'subset_values' are required when subset_cells = TRUE.")
    keep <- seurat_object@meta.data[[subset_by]] %in% subset_values
    seurat_object  <- seurat_object[, keep]
    message("Subset to ", sum(keep), " cells where ", subset_by,
            " in {", paste(subset_values, collapse = ", "), "}.")
  }

  # ── 2. Load CSVs ──────────────────────────────────────────────────────────
  message("\nStep 1 | Loading ", length(gsea_files), " GSEA CSV(s)...")
  all_rows <- lapply(seq_along(gsea_files), function(i) {
    df <- utils::read.csv(gsea_files[i], stringsAsFactors = FALSE)
    required_cols <- c("pathway", "leadingEdge")
    missing_cols  <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0)
      stop("File '", basename(gsea_files[i]), "' is missing columns: ",
           paste(missing_cols, collapse = ", "))
    df$._src_file <- basename(gsea_files[i])
    df
  })
  all_df <- do.call(rbind, all_rows)

  # ── 3. Filter pathways ────────────────────────────────────────────────────
  if (!is.null(search_terms)) {
    hits   <- .apply_search_terms(all_df$pathway, search_terms)
    all_df <- all_df[hits, , drop = FALSE]
    term_label <- .format_search_terms(search_terms)
  } else {
    all_df <- all_df[all_df$pathway %in% pathways, , drop = FALSE]
    raw_lbl    <- paste(pathways, collapse = "_")
    term_label <- gsub("[^A-Za-z0-9_-]", "_", raw_lbl)
    if (nchar(term_label) > 60)
      term_label <- paste0(substr(term_label, 1, 57), "...")
  }

  if (nrow(all_df) == 0)
    stop("No pathways matched your search. Check 'search_terms' or 'pathways'.")

  message("  Found ", nrow(all_df), " matching row(s) across ",
          length(unique(all_df$._src_file)), " file(s):")
  for (pw in unique(all_df$pathway))
    message("    • ", pw)

  # ── 4. Extract leading-edge genes ─────────────────────────────────────────
  message("\nStep 2 | Extracting leading-edge genes...")
  le_genes <- unique(unlist(lapply(all_df$leadingEdge, .parse_leading_edge)))
  le_genes <- le_genes[!is.na(le_genes) & nchar(le_genes) > 0]
  message("  Union of leading-edge genes across all rows: ", length(le_genes))

  # ── 5. Create output subfolder ────────────────────────────────────────────
  lbl     <- gsub("[^A-Za-z0-9_-]", "_", term_label)
  lbl     <- gsub("_+", "_", lbl)
  pfx_parts <- c(object_name, subset_name)
  pfx     <- paste(pfx_parts[nchar(pfx_parts) > 0], collapse = "_")
  out_sub <- file.path(output_dir, paste0(pfx, "_LeadingEdge_", lbl))
  dir.create(out_sub, showWarnings = FALSE, recursive = TRUE)

  # ── 6. Pull expression matrix ─────────────────────────────────────────────
  message("\nStep 3 | Extracting expression matrix...")
  # .get_layer_data() handles Seurat v3/v5/BPCells transparently
  mat_full <- .get_layer_data(seurat_object, assay = assay, layer = layer)

  genes_present <- intersect(le_genes, rownames(mat_full))
  n_missing     <- length(le_genes) - length(genes_present)
  if (length(genes_present) == 0)
    stop("None of the leading-edge genes are present in the Seurat object (",
         assay, " assay, ", layer, " layer).")
  message("  ", length(genes_present), " / ", length(le_genes),
          " leading-edge genes found in the object",
          if (n_missing > 0) paste0(" (", n_missing, " not found — check species/assay)") else ".")

  # Save full gene list
  utils::write.csv(
    data.frame(gene       = le_genes,
               in_object  = le_genes %in% rownames(mat_full),
               stringsAsFactors = FALSE),
    file.path(out_sub, "leading_edge_genes.csv"),
    row.names = FALSE
  )

  mat_sub <- as.matrix(mat_full[genes_present, , drop = FALSE])

  # Filter low-mean genes
  keep_genes <- rowMeans(mat_sub) >= min_expr
  mat_sub    <- mat_sub[keep_genes, , drop = FALSE]
  message("  After min_expr (>= ", min_expr, ") filter: ", nrow(mat_sub), " genes retained.")

  if (nrow(mat_sub) == 0)
    stop("All genes removed by min_expr filter. Try lowering 'min_expr'.")

  genes_all <- rownames(mat_sub)

  # ── 7. Metadata — factor levels & column ordering ─────────────────────────
  meta <- seurat_object@meta.data

  .get_levels <- function(col_name) {
    col <- meta[[col_name]]
    if (is.factor(col)) levels(col) else unique(as.character(col))
  }

  group_levels <- .get_levels(group.by)
  meta[[group.by]] <- factor(as.character(meta[[group.by]]), levels = group_levels)

  if (!is.null(split.by)) {
    split_levels <- .get_levels(split.by)
    meta[[split.by]] <- factor(as.character(meta[[split.by]]), levels = split_levels)
    ord <- order(meta[[group.by]], meta[[split.by]])
  } else {
    split_levels <- NULL
    ord <- order(meta[[group.by]])
  }

  cell_ord <- intersect(rownames(meta)[ord], colnames(mat_sub))
  mat_ord  <- mat_sub[, cell_ord, drop = FALSE]

  # Row-center (subtract row mean)
  mat_centered <- mat_ord - rowMeans(mat_ord)
  mat_centered  <- stats::na.omit(mat_centered)

  # ── 8. Colors — priority: explicit arg > PrepObject stored > Nour_pal auto ──
  n_grp <- length(group_levels)
  n_spl <- length(split_levels)
  if (is.null(group_colors))
    group_colors <- .nk_colors(seurat_object, group.by) %||%
      stats::setNames(
        scSidekick::Nour_pal(if (n_grp <= 8) "all" else "spectrum")(n_grp),
        group_levels)
  if (!is.null(split.by) && is.null(split_colors))
    split_colors <- .nk_colors(seurat_object, split.by) %||%
      stats::setNames(scSidekick::Nour_pal("spectrum")(n_spl), split_levels)

  # ── 9. Column annotation & legend objects ────────────────────────────────
  annot_cols  <- c(group.by, if (!is.null(split.by)) split.by)
  annot_df    <- meta[cell_ord, annot_cols, drop = FALSE]
  col_colors  <- list()
  col_colors[[group.by]] <- group_colors
  if (!is.null(split.by)) col_colors[[split.by]] <- split_colors

  colanno <- ComplexHeatmap::columnAnnotation(
    df                   = annot_df,
    col                  = col_colors,
    show_annotation_name = TRUE,
    show_legend          = FALSE
  )

  lgd_group <- ComplexHeatmap::Legend(
    labels    = group_levels,
    title     = group.by,
    legend_gp = grid::gpar(fill = group_colors)
  )
  lgd_list <- list(lgd_group)
  if (!is.null(split.by)) {
    lgd_list <- c(lgd_list, list(
      ComplexHeatmap::Legend(
        labels    = split_levels,
        title     = split.by,
        legend_gp = grid::gpar(fill = split_colors)
      )
    ))
  }

  # ── 10. Shared heatmap arguments ──────────────────────────────────────────
  col_fun <- circlize::colorRamp2(c(-2, 0, 2),
                                   c("#007dd1", "white", "#ab3000"))

  .make_heatmap <- function(mat, title_str) {
    n_g <- nrow(mat)
    default_args <- list(
      matrix            = as.matrix(mat),
      name              = "logcounts\n(centered)",
      col               = col_fun,
      top_annotation    = colanno,
      cluster_rows      = FALSE,
      cluster_columns   = FALSE,
      show_column_names = FALSE,
      show_row_names    = TRUE,
      row_names_gp      = grid::gpar(fontsize = max(4, 30 / sqrt(n_g))),
      border            = FALSE,
      use_raster        = TRUE,
      column_title      = title_str,
      row_names_side      = "left",
      show_row_dend       = FALSE,
      row_names_max_width = grid::unit(15, "cm")
    )
    hm_args <- modifyList(default_args, heatmap_params)
    do.call(ComplexHeatmap::Heatmap, hm_args)
  }

  # ── 11. Heatmap A: all leading-edge genes ─────────────────────────────────
  message("\nStep 4 | Drawing heatmap — all leading-edge genes (",
          nrow(mat_centered), " genes × ", ncol(mat_centered), " cells)...")

  HM_all   <- .make_heatmap(mat_centered,
                             paste0("Leading Edge: ", term_label))
  pdf_all  <- file.path(out_sub, "heatmap_all_leadingedge.pdf")
  pdf(pdf_all,
      width  = 9,
      height = max(6, nrow(mat_centered) * 0.18 + 3))
  ComplexHeatmap::draw(HM_all,
                       heatmap_legend_list = lgd_list,
                       heatmap_legend_side = "right")
  grDevices::dev.off()
  message("  Saved: ", basename(pdf_all))

  .write_legend_sidecar(pdf_all, paste0(
    "Heatmap: Mean-centered log-normalized expression of leading-edge genes ",
    "from GSEA pathways matching [", term_label, "]. ",
    "Genes: union across all provided CSVs. ",
    "Color: blue = -2, white = 0, red = +2 (centered by row mean). ",
    "Genes with row mean < ", min_expr, " excluded. ",
    "Cells ordered by ", group.by,
    if (!is.null(split.by)) paste0(" then ", split.by) else "", "."
  ))

  # ── 12. DEG-filtered bonus plots ──────────────────────────────────────────
  genes_top <- NULL
  top_degs  <- NULL

  if (!is.null(deg_df)) {
    message("\nStep 5 | Filtering DEGs within leading edge...")

    # Subset to significant DEGs that are in the leading edge
    sig_mask <- (deg_df[[deg_gene_column]] %in% genes_all) &
      (!is.na(deg_df[[deg_padj_column]])) &
      (deg_df[[deg_padj_column]] < deg_padj_cutoff)
    deg_sig <- deg_df[sig_mask, , drop = FALSE]
    message("  Significant DEGs (padj < ", deg_padj_cutoff,
            ") within leading edge: ", nrow(deg_sig))

    if (nrow(deg_sig) > 0) {

      # Top N by logFC per group
      top_degs <- do.call(rbind, lapply(
        split(deg_sig, deg_sig[[deg_group_column]]),
        function(g) {
          g <- g[order(g[[deg_fc_column]], decreasing = TRUE), , drop = FALSE]
          g[seq_len(min(deg_top_n, nrow(g))), , drop = FALSE]
        }
      ))
      genes_top <- unique(top_degs[[deg_gene_column]])
      message("  Top DEGs (up to ", deg_top_n, " per group): ",
              length(genes_top), " unique genes.")

      # Save CSV
      top_csv <- file.path(out_sub, "top_degs_in_leadingedge.csv")
      utils::write.csv(top_degs, top_csv, row.names = FALSE)
      message("  Saved: ", basename(top_csv))

      # ── 12a. Heatmap B: top DEGs only ─────────────────────────────────
      genes_top_hm <- intersect(genes_top, rownames(mat_centered))
      if (length(genes_top_hm) > 0) {
        mat_top  <- mat_centered[genes_top_hm, , drop = FALSE]
        HM_top   <- .make_heatmap(mat_top,
                                   paste0("Top DEGs in Leading Edge: ", term_label))
        pdf_top  <- file.path(out_sub, "heatmap_top_degs.pdf")
        pdf(pdf_top,
            width  = 9,
            height = max(5, nrow(mat_top) * 0.22 + 3))
        ComplexHeatmap::draw(HM_top,
                             heatmap_legend_list = lgd_list,
                             heatmap_legend_side = "right")
        grDevices::dev.off()
        message("  Saved: ", basename(pdf_top))
        .write_legend_sidecar(pdf_top, paste0(
          "Heatmap: Top ", deg_top_n, " DEGs per group (padj < ",
          deg_padj_cutoff, ") within leading-edge genes of [",
          term_label, "]. Mean-centered log-normalized expression."
        ))
      }

      # ── 12b. DotPlot ────────────────────────────────────────────────────
      if (isTRUE(add_dotplot)) {
        genes_dp <- intersect(genes_top, rownames(seurat_object))
        if (length(genes_dp) > 0) {
          message("  Drawing DotPlot...")
          p_dot <- Seurat::DotPlot(
            seurat_object, group.by = group.by,
            features  = genes_dp,
            dot.scale = 5, scale = TRUE
          ) +
            ggplot2::scale_colour_viridis_c(option = "plasma") +
            ggplot2::theme_minimal() +
            Seurat::RotatedAxis() +
            ggplot2::theme(axis.title = ggplot2::element_blank())

          pdf_dot <- file.path(out_sub, "dotplot_top_degs.pdf")
          pdf(pdf_dot,
              width  = max(6, length(genes_dp) * 0.35 + 2),
              height = max(3, length(group_levels) * 0.5 + 1.5))
          print(p_dot)
          grDevices::dev.off()
          message("  Saved: ", basename(pdf_dot))
        }
      }

      # ── 12c. Stacked Violin ──────────────────────────────────────────────
      if (isTRUE(add_violin)) {
        genes_vln <- intersect(genes_top, rownames(seurat_object))
        if (length(genes_vln) > 0) {
          message("  Drawing stacked violin...")
          pdf_vln <- file.path(out_sub, "stacked_violin.pdf")
          vln_h   <- max(3, length(genes_vln) * 0.65 + 1)
          vln_w   <- max(4, length(group_levels) * 0.8 + 1)
          pdf(pdf_vln, width = vln_w, height = vln_h)

          # Panel 1: grouped by group.by
          p_v1 <- StackedVlnPlot(
            seurat_object, features = genes_vln,
            group.by = group.by, cols = group_colors, pt.size = 0
          )
          print(p_v1)

          # Panel 2: grouped by group.by, split by split.by (if present)
          if (!is.null(split.by)) {
            p_v2 <- StackedVlnPlot(
              seurat_object, features = genes_vln,
              group.by = group.by, split.by = split.by,
              cols = split_colors, pt.size = 0.01
            ) + Seurat::RotatedAxis()
            print(p_v2)
          }

          grDevices::dev.off()
          message("  Saved: ", basename(pdf_vln))
        }
      }

      # ── 12d. Boxplots per gene ───────────────────────────────────────────
      if (isTRUE(add_boxplots)) {
        genes_bx <- intersect(genes_top, rownames(seurat_object))
        if (length(genes_bx) > 0) {
          message("  Drawing boxplots (", length(genes_bx), " genes)...")

          # Build per-cell expression data frame (reuse already-materialised mat_full)
          expr_mat <- as.matrix(mat_full[genes_bx, , drop = FALSE])
          # Long format via base R (no tidyr needed)
          expr_lng <- do.call(rbind, lapply(genes_bx, function(g) {
            df_g <- data.frame(
              gene       = g,
              expression = as.numeric(expr_mat[g, ]),
              stringsAsFactors = FALSE
            )
            df_g[[group.by]] <- factor(
              as.character(meta[colnames(expr_mat), group.by]),
              levels = group_levels
            )
            if (!is.null(split.by))
              df_g[[split.by]] <- factor(
                as.character(meta[colnames(expr_mat), split.by]),
                levels = split_levels
              )
            df_g
          }))

          # X-axis variable and fill colors
          x_var  <- if (!is.null(split.by)) split.by else group.by
          x_cols <- if (!is.null(split.by)) split_colors else group_colors

          .bx_theme <- ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.line        = ggplot2::element_line(colour = "black", linewidth = 0.4),
            strip.background = ggplot2::element_rect(fill = "white", colour = "black",
                                                      linewidth = 0.4),
            strip.text       = ggplot2::element_text(face = "bold", colour = "black",
                                                      size = 8),
            plot.title       = ggplot2::element_text(face = "bold", colour = "black",
                                                      hjust = 0.5, size = 9),
            axis.text.y      = ggplot2::element_text(colour = "black", size = 6),
            axis.text.x      = ggplot2::element_text(colour = "black", angle = 25,
                                                      hjust = 1, size = 6.5,
                                                      face = "bold"),
            axis.title       = ggplot2::element_text(size = 8, face = "bold"),
            plot.margin      = ggplot2::unit(c(0.3, 0.5, 0.2, 0.5), "cm")
          )

          pdf_bx <- file.path(out_sub, "boxplots_top_degs.pdf")
          pdf(pdf_bx, width = 5, height = 5.5)
          for (g in genes_bx) {
            df_g <- expr_lng[expr_lng$gene == g, , drop = FALSE]
            p_bx <- ggplot2::ggplot(
              df_g,
              ggplot2::aes(
                x    = .data[[x_var]],
                y    = .data[["expression"]],
                fill = .data[[x_var]]
              )
            ) +
              ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
              ggplot2::geom_dotplot(
                binaxis  = "y", stackdir = "center",
                dotsize  = 0.4, fill = "black", alpha = 0.3
              ) +
              ggplot2::scale_fill_manual(values = x_cols) +
              ggplot2::labs(title = g, x = NULL, y = "Expression") +
              ggplot2::guides(fill = "none") +
              .bx_theme

            # Facet by group.by when split.by is the x-axis
            if (!is.null(split.by))
              p_bx <- p_bx +
                ggplot2::facet_wrap(
                  stats::as.formula(paste("~", group.by)),
                  scales = "free"
                )

            # p-values via ggpubr
            if (isTRUE(add_pvalues) &&
                requireNamespace("ggpubr", quietly = TRUE)) {
              p_bx <- p_bx +
                ggpubr::stat_compare_means(
                  method = if (length(unique(df_g[[x_var]])) == 2)
                    "wilcox.test" else "kruskal.test",
                  label  = "p.signif",
                  size   = 3
                )
            }
            print(p_bx)
          }
          grDevices::dev.off()
          message("  Saved: ", basename(pdf_bx))
        }
      }

    } else {
      message("  No significant DEGs found within leading-edge genes at padj < ",
              deg_padj_cutoff, ". Skipping bonus plots.")
    }
  }

  # ── 13. Feature maps (optional) ───────────────────────────────────────────
  if (isTRUE(add_feature_maps)) {
    fm_genes <- if (!is.null(genes_top)) genes_top else genes_all
    fm_genes <- intersect(fm_genes, rownames(seurat_object))
    if (length(fm_genes) > 0) {
      message("\nStep 6 | Drawing feature maps for ", length(fm_genes), " genes...")
      fm_dir <- file.path(out_sub, "feature_maps")
      dir.create(fm_dir, showWarnings = FALSE, recursive = TRUE)
      for (Y in fm_genes) {
        X <- Seurat::FeaturePlot(
          seurat_object, features = Y,
          split.by = split.by,
          order    = TRUE, label = FALSE,
          pt.size  = 0.1
        )
        n_panels <- if (!is.null(split.by)) length(split_levels) else 1
        pdf_fm <- file.path(fm_dir, paste0(Y, "_featuremap.pdf"))
        pdf(pdf_fm, width = n_panels * 2.5 + 0.5, height = 2.5)
        if (is.list(X)) {
          print(patchwork::wrap_plots(X, nrow = 1))
        } else {
          print(X)
        }
        grDevices::dev.off()
      }
      message("  Feature maps saved to: feature_maps/")
    }
  }

  # ── 14. Methods JSON ──────────────────────────────────────────────────────
  .write_subdir_params(out_sub, list(
    function_name        = "VisualizeLeadingEdge",
    gsea_files           = basename(gsea_files),
    search_terms         = if (!is.null(search_terms)) search_terms else pathways,
    pathways_matched     = unique(all_df$pathway),
    n_leading_edge_genes = length(le_genes),
    n_genes_in_object    = length(genes_all),
    group.by             = group.by,
    split.by             = split.by,
    assay                = assay,
    layer                = layer,
    min_expr             = min_expr,
    deg_padj_cutoff      = if (!is.null(deg_df)) deg_padj_cutoff else NA,
    deg_top_n            = if (!is.null(deg_df)) deg_top_n else NA,
    date                 = as.character(Sys.time())
  ))

  message("\nVisualizeLeadingEdge complete.\nOutput: ", out_sub)
  invisible(list(genes_all = genes_all,
                 genes_top = genes_top,
                 top_degs  = top_degs))
}
