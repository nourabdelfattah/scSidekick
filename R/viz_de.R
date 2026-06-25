# =============================================================================
# scSidekick - Differential Expression & Pathway Visualizations  (viz_de.R)
#
# Exported:
#   PlotVolcano()          - volcano plot from presto / FindMarkers output
#   PlotGSEAEnrichment()   - 3-panel running-score mountain plot from fgsea
#   PlotEnrichment()       - ORA bar / dot / lollipop from enrichR / clusterProfiler
# =============================================================================


# ── Shared internal helpers ──────────────────────────────────────────────────

# Auto-detect a column from a list of candidates (first match wins)
.detect_col <- function(df, candidates) {
  found <- intersect(candidates, colnames(df))
  if (length(found)) found[1] else NULL
}

# scSidekick DE theme (white panel, black axes, bold text)
.de_theme <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(color = "black", linewidth = 0.5),
      axis.text         = ggplot2::element_text(color = "black"),
      axis.title        = ggplot2::element_text(face = "bold", color = "black"),
      strip.background  = ggplot2::element_rect(fill = "white", color = "black",
                                                 linewidth = 0.5),
      strip.text        = ggplot2::element_text(face = "bold", color = "black"),
      plot.title        = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.key        = ggplot2::element_rect(fill = NA),
      plot.margin       = ggplot2::unit(c(0.4, 0.5, 0.3, 0.4), "cm")
    )
}

# Compute the running enrichment score vector from a ranked statistic and gene set
# full_set_size: the full gene-set size from fgsea output (the 'size' column).
# When provided, the miss penalty is computed against the full set rather than
# just the leading-edge subset, giving a more accurate curve shape.
.calc_running_es <- function(stats, leading_edge, full_set_size = NULL) {
  le <- intersect(leading_edge, names(stats))
  n  <- length(stats)
  k_le    <- length(le)
  k_total <- full_set_size %||% k_le
  if (k_le == 0L || n == 0L)
    return(data.frame(rank = seq_len(n), gene = names(stats),
                      score = 0, is_hit = FALSE))
  nr       <- sum(abs(stats[le]))
  if (nr == 0) nr <- 1  # guard against all-zero stats
  is_hit   <- names(stats) %in% le
  hit_step <- ifelse(is_hit,  abs(stats) / nr,             0)
  mis_step <- ifelse(!is_hit, -1 / max(n - k_total, 1L),  0)
  data.frame(
    rank   = seq_len(n),
    gene   = names(stats),
    score  = cumsum(hit_step + mis_step),
    is_hit = is_hit
  )
}

# ── RunGSEA directory helpers ────────────────────────────────────────────────

# Portable regex escape: uses fixed=TRUE per character, avoids character-class
# syntax that BSD regex (macOS) rejects (e.g. [][...]).
.re_escape <- function(x) {
  # Backslash must be handled first to avoid double-escaping
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  for (ch in c(".", "+", "*", "?", "^", "$", "{", "}", "(", ")", "|", "[", "]"))
    x <- gsub(ch, paste0("\\", ch), x, fixed = TRUE)
  x
}

# Discover all (label, sp, cluster) → CSV path mappings inside a RunGSEA output dir.
# File naming convention from RunGSEA:
#   CSV:      {db_name} {cl_safe} {sp_safe} {lab}.csv
#   de_cache: de_cache_{sp_safe}_{lab}.rds
.discover_gsea_files <- function(gsea_dir, database) {
  lab_dirs <- list.dirs(gsea_dir, full.names = TRUE, recursive = FALSE)
  out <- lapply(lab_dirs, function(lab_dir) {
    lab    <- basename(lab_dir)
    lab_re <- .re_escape(lab)

    # Which databases are present?
    db_dirs <- if (!is.null(database)) {
      file.path(lab_dir, database)
    } else {
      list.dirs(lab_dir, full.names = TRUE, recursive = FALSE)
    }
    db_dirs <- db_dirs[dir.exists(db_dirs)]

    lapply(db_dirs, function(db_dir) {
      db_nm <- basename(db_dir)
      db_re <- .re_escape(db_nm)

      # Discover known split levels from de_cache filenames
      # Pattern: de_cache_{sp_safe}_{lab}.rds
      de_files   <- list.files(lab_dir, pattern = "^de_cache_.*\\.rds$",
                                full.names = FALSE)
      sp_pattern <- paste0("^de_cache_(.+)_", lab_re, "\\.rds$")
      sp_safes   <- sub(sp_pattern, "\\1",
                        de_files[grepl(sp_pattern, de_files, perl = TRUE)])
      if (length(sp_safes) == 0) sp_safes <- "All"

      csvs <- list.files(db_dir, pattern = "\\.csv$", full.names = TRUE)

      lapply(csvs, function(csv_path) {
        base <- basename(csv_path)
        # Strip "{db_name} " prefix and " {lab}.csv" suffix to get
        # the middle part "{cl_safe} {sp_safe}"
        mid <- sub(paste0("^", db_re, " "), "", base, perl = TRUE)
        mid <- sub(paste0(" ", lab_re, "\\.csv$"), "", mid,
                   ignore.case = TRUE, perl = TRUE)

        # Match sp_safe from the right of mid
        sp_match <- NA_character_
        cl_match <- NA_character_
        for (sp in sp_safes) {
          sp_re <- .re_escape(sp)
          suffix <- paste0(" ", sp_re, "$")
          if (grepl(suffix, mid, perl = TRUE)) {
            sp_match <- sp
            cl_match <- sub(suffix, "", mid, perl = TRUE)
            break
          }
          if (identical(mid, sp)) {      # no split: sp == "All", cl is empty
            sp_match <- sp
            cl_match <- ""
            break
          }
        }

        list(
          lab      = lab,
          db       = db_nm,
          sp       = sp_match,
          cl       = cl_match,
          csv_path = csv_path,
          de_cache = file.path(lab_dir,
                               paste0("de_cache_", sp_match, "_", lab, ".rds"))
        )
      })
    })
  })

  # Flatten and drop unmatched entries (NA sp)
  flat <- unlist(unlist(out, recursive = FALSE), recursive = FALSE)
  flat[!sapply(flat, function(x) is.na(x$sp))]
}


# =============================================================================
# PlotVolcano
# =============================================================================

#' Volcano Plot of Differential Expression Results
#'
#' @description
#' Produces a publication-ready volcano plot from a DEG data frame (presto,
#' \code{FindMarkers}, or any table with logFC and adjusted p-value columns).
#' Automatically detects column names for common output formats.  Points are
#' colored by direction × significance; top genes per direction are labelled
#' with \code{ggrepel}.
#'
#' @param deg_df Data frame with differential expression results. Compatible
#'   with presto (\code{"feature"/"logFC"/"padj"}) and Seurat
#'   (\code{"gene"/"avg_log2FC"/"p_val_adj"}).
#' @param gene_column Column with gene names. Auto-detected if \code{NULL}.
#' @param fc_column Column with log fold-change. Auto-detected if \code{NULL}.
#' @param padj_column Column with adjusted p-value. Auto-detected if \code{NULL}.
#' @param group_column Optional metadata column to facet by (one panel per
#'   group level). \code{NULL} = no faceting.
#' @param fc_cutoff Absolute log fold-change threshold (vertical dashed lines).
#'   Default \code{0.25}.
#' @param padj_cutoff Adjusted p-value significance cutoff (horizontal dashed
#'   line). Default \code{0.05}.
#' @param top_n_labels Number of top genes to label per direction (up and down
#'   separately). Genes ranked by \eqn{-\log_{10}(padj) \times |logFC|}.
#'   Default \code{10}.
#' @param highlight_genes Character vector of specific genes to always label
#'   regardless of their ranking.
#' @param colors Named character vector with elements \code{"up"}, \code{"down"},
#'   and \code{"ns"} (not significant). Defaults to red/blue/gray.
#' @param point_size Point size. Default \code{1}.
#' @param point_alpha Point transparency. Default \code{0.6}.
#' @param max_overlaps Maximum label overlaps passed to \code{ggrepel}.
#'   Default \code{20}.
#' @param output_dir Directory to save a PDF. \code{NULL} = plot only, no save.
#' @param object_name Label prefix for the output file name.
#' @param subset_name Optional subset label appended to the file name.
#'
#' @return A \code{ggplot2} object (invisibly if saved). Printed to the active
#'   device.
#'
#' @export
PlotVolcano <- function(
    deg_df,
    gene_column     = NULL,
    fc_column       = NULL,
    padj_column     = NULL,
    group_column    = NULL,
    fc_cutoff       = 0.25,
    padj_cutoff     = 0.05,
    top_n_labels    = 10,
    highlight_genes = NULL,
    colors          = c(up = "#B40426", down = "#3B4CC0", ns = "gray70"),
    point_size      = 1,
    point_alpha     = 0.6,
    max_overlaps    = 20,
    output_dir      = NULL,
    object_name     = "",
    subset_name     = ""
) {
  if (!is.data.frame(deg_df) || nrow(deg_df) == 0)
    stop("'deg_df' must be a non-empty data frame.")

  # ── Auto-detect columns ──────────────────────────────────────────────────
  gene_col  <- gene_column  %||% .detect_col(deg_df, c("feature","gene","Gene","Symbol","name"))
  fc_col    <- fc_column    %||% .detect_col(deg_df, c("logFC","log2FC","avg_log2FC","log2FoldChange","FC"))
  padj_col  <- padj_column  %||% .detect_col(deg_df, c("padj","p_val_adj","p.adjust","FDR","adj.P.Val"))

  missing <- c(
    if (is.null(gene_col)) "gene column",
    if (is.null(fc_col))   "logFC column",
    if (is.null(padj_col)) "padj column"
  )
  if (length(missing))
    stop("Could not auto-detect: ", paste(missing, collapse = ", "),
         ". Please specify the relevant *_column parameter(s).")

  df <- deg_df
  df$._gene  <- as.character(df[[gene_col]])
  df$._fc    <- as.numeric(df[[fc_col]])
  df$._padj  <- as.numeric(df[[padj_col]])
  df$._padj  <- pmax(df$._padj, 1e-300)   # floor to avoid Inf in -log10
  df$._neglog10padj <- -log10(df$._padj)

  # ── Classify direction ───────────────────────────────────────────────────
  df$._dir <- dplyr::case_when(
    df$._fc >  fc_cutoff & df$._padj < padj_cutoff ~ "up",
    df$._fc < -fc_cutoff & df$._padj < padj_cutoff ~ "down",
    TRUE ~ "ns"
  )

  # ── Select genes to label ────────────────────────────────────────────────
  df$._score <- df$._neglog10padj * abs(df$._fc)
  label_df <- df[df$._dir != "ns", ]

  top_up   <- utils::head(
    label_df[order(-label_df$._score[label_df$._dir == "up"]), ][label_df[order(-label_df$._score), ]$._dir == "up", ],
    top_n_labels
  )
  top_dn   <- utils::head(
    label_df[order(-label_df$._score[label_df$._dir == "down"]), ][label_df[order(-label_df$._score), ]$._dir == "down", ],
    top_n_labels
  )

  # Simpler approach - top by score per direction
  if (!is.null(group_column) && group_column %in% colnames(df)) {
    label_genes <- do.call(rbind, lapply(
      split(df, df[[group_column]]),
      function(g) {
        up_g  <- utils::head(g[g$._dir == "up"  & order(g$._score[g$._dir == "up"],   decreasing = TRUE), ], top_n_labels)
        dn_g  <- utils::head(g[g$._dir == "down" & order(g$._score[g$._dir == "down"], decreasing = TRUE), ], top_n_labels)
        rbind(up_g, dn_g)
      }
    ))
  } else {
    up_idx <- which(df$._dir == "up")
    dn_idx <- which(df$._dir == "down")
    top_up_idx <- up_idx[order(df$._score[up_idx], decreasing = TRUE)[seq_len(min(top_n_labels, length(up_idx)))]]
    top_dn_idx <- dn_idx[order(df$._score[dn_idx], decreasing = TRUE)[seq_len(min(top_n_labels, length(dn_idx)))]]
    label_genes <- df[c(top_up_idx, top_dn_idx), ]
  }

  if (!is.null(highlight_genes)) {
    extra <- df[df$._gene %in% highlight_genes, ]
    label_genes <- unique(rbind(label_genes, extra))
  }

  # ── Symmetric x-axis ─────────────────────────────────────────────────────
  max_fc  <- max(abs(df$._fc), na.rm = TRUE) * 1.05
  max_nlp <- max(df$._neglog10padj, na.rm = TRUE) * 1.05

  # ── Colors ───────────────────────────────────────────────────────────────
  cols <- c(up = "#B40426", down = "#3B4CC0", ns = "gray70")
  cols[names(colors)] <- colors
  df$._dir <- factor(df$._dir, levels = c("up", "down", "ns"))

  # ── Count annotations ────────────────────────────────────────────────────
  n_up <- sum(df$._dir == "up",   na.rm = TRUE)
  n_dn <- sum(df$._dir == "down", na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x     = .data[["._fc"]],
    y     = .data[["._neglog10padj"]],
    color = .data[["._dir"]]
  )) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_manual(
      values = cols,
      labels = c(up   = paste0("Up (", n_up, ")"),
                 down = paste0("Down (", n_dn, ")"),
                 ns   = "Not significant"),
      name   = NULL
    ) +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
                        linetype = "dashed", color = "gray40", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff),
                        linetype = "dashed", color = "gray40", linewidth = 0.4) +
    ggplot2::xlim(-max_fc, max_fc) +
    ggplot2::ylim(0, max_nlp) +
    ggplot2::labs(x = "log₂ Fold Change",
                  y = expression(-log[10](p[adj]))) +
    .de_theme()

  # Labels
  if (nrow(label_genes) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data        = label_genes,
      ggplot2::aes(label = .data[["._gene"]]),
      size        = 2.8,
      fontface    = "bold",
      color       = "black",
      box.padding = 0.35,
      max.overlaps = max_overlaps,
      segment.size = 0.2,
      show.legend = FALSE
    )
  }

  # Facet
  if (!is.null(group_column) && group_column %in% colnames(df)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", group_column)))
  }

  # ── Save ─────────────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    pfx   <- paste(c(object_name, subset_name)[nchar(c(object_name, subset_name)) > 0],
                   collapse = "_")
    grp_n <- if (!is.null(group_column)) length(unique(df[[group_column]])) else 1
    pdf(file.path(output_dir, paste0(pfx, " volcano.pdf")),
        width  = max(5, grp_n * 4),
        height = 5)
    print(p)
    grDevices::dev.off()
  }

  print(p)
  invisible(p)
}


# =============================================================================
# PlotGSEAEnrichment - internal mountain-plot builder
# =============================================================================

# Build one 3-panel mountain-plot stack for a single pathway.
# ranked_stats: named numeric vector sorted descending (gene → AUC/logFC)
# fgsea_row:   one-row data frame from fgsea CSV (must have pathway, NES, leadingEdge, size)
# col_pos/neg: colors for positive/negative NES
# title_extra: appended to the pathway title (e.g. " | NRF1_mRNA")
.one_mountain <- function(ranked_stats, fgsea_row, col_pos, col_neg,
                           title_extra = "") {
  nes   <- suppressWarnings(as.numeric(fgsea_row$NES))
  pval  <- suppressWarnings(as.numeric(
    if ("padj" %in% colnames(fgsea_row)) fgsea_row$padj else fgsea_row$pval
  ))
  fset_size <- suppressWarnings(as.integer(fgsea_row$size))
  plab      <- if (!is.na(pval) && pval < 0.001)
    formatC(pval, format = "e", digits = 2)
  else
    formatC(as.numeric(pval), digits = 3, format = "fg")
  curve_col <- if (!is.na(nes) && nes >= 0) col_pos else col_neg

  # Parse leadingEdge
  le_raw <- fgsea_row$leadingEdge
  le <- tryCatch(eval(parse(text = le_raw)),
                 error = function(e) trimws(strsplit(le_raw, "[,\\s]+")[[1]]))

  es_df  <- .calc_running_es(ranked_stats, le,
                              full_set_size = if (!is.na(fset_size)) fset_size else NULL)
  peak_i <- which.max(abs(es_df$score))

  pw_title <- paste0(fgsea_row$pathway, title_extra)

  # Panel 1: running ES
  p1 <- ggplot2::ggplot(es_df, ggplot2::aes(x = rank, y = score)) +
    ggplot2::geom_area(fill = curve_col, alpha = 0.12) +
    ggplot2::geom_line(color = curve_col, linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "gray40") +
    ggplot2::geom_vline(xintercept = peak_i, linetype = "dotted",
                        color = "gray40", linewidth = 0.4) +
    ggplot2::annotate("text",
                      x     = peak_i + nrow(es_df) * 0.02,
                      y     = es_df$score[peak_i],
                      label = paste0("NES = ", round(nes, 3),
                                     "\npadj = ", plab),
                      hjust = 0, vjust = 0.5,
                      size  = 2.6, fontface = "bold", color = "gray20") +
    ggplot2::labs(title = pw_title, x = NULL, y = "ES") +
    .de_theme(base_size = 9) +
    ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title   = ggplot2::element_text(size = 8, face = "bold"))

  # Panel 2: leading-edge rug
  rug_df <- es_df[es_df$is_hit, ]
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_vline(data = rug_df,
                        ggplot2::aes(xintercept = rank),
                        color = curve_col, linewidth = 0.25, alpha = 0.8) +
    ggplot2::xlim(1, nrow(es_df)) +
    ggplot2::labs(x = NULL, y = NULL) +
    .de_theme(base_size = 9) +
    ggplot2::theme(axis.text  = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())

  # Panel 3: AUC gradient bar
  stat_df <- data.frame(rank = seq_along(ranked_stats), stat = ranked_stats)
  p3 <- ggplot2::ggplot(stat_df,
                         ggplot2::aes(x = rank, y = 1, fill = stat)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = c(col_neg, "white", col_pos),
      name    = "AUC",
      guide   = ggplot2::guide_colorbar(barwidth = 3, barheight = 0.35,
                                         title.position = "top",
                                         title.hjust    = 0.5)
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      name   = paste0("Gene rank (n=", length(ranked_stats), ")")
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(y = NULL) +
    .de_theme(base_size = 8) +
    ggplot2::theme(axis.text.y     = ggplot2::element_blank(),
                   axis.ticks.y    = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.text     = ggplot2::element_text(size = 6))

  patchwork::wrap_plots(p1, p2, p3, ncol = 1, heights = c(4, 0.5, 1.1))
}


# =============================================================================
# PlotGSEAEnrichment
# =============================================================================

#' GSEA Running Score Mountain Plots - Direct or from RunGSEA Output Directory
#'
#' @description
#' Two operating modes:
#'
#' \strong{Direct mode} - pass \code{fgsea_result} and \code{ranked_stats}
#' directly (e.g. from an ad-hoc fgsea run or a single CSV).
#'
#' \strong{RunGSEA directory mode} - point to the \code{output_dir} of a
#' completed \code{\link{RunGSEA}} run via \code{gsea_dir}.  The function
#' auto-discovers all CSV files, searches pathway names with keyword search
#' (same AND/OR logic as \code{\link{RunSCssGSEA}}), loads the ranked AUC
#' vectors from the cached DE RDS files, and produces a grid layout:
#' \itemize{
#'   \item \strong{One PDF per split.by level} (cell type) - each PDF may have
#'     multiple pages when there are more groups than \code{max_cols}
#'   \item \strong{Columns} = group.by levels (comparison groups, e.g. Male / Female)
#'     - up to \code{max_cols} per page; additional groups wrap to new pages
#'   \item \strong{Rows} = matched pathways (one 3-panel mountain stack per pathway)
#' }
#'
#' @param fgsea_result \emph{Direct mode only.} Data frame with fgsea output or
#'   a path to a single RunGSEA CSV.  Must have \code{pathway}, \code{NES},
#'   \code{leadingEdge}, and \code{size} columns.
#' @param ranked_stats \emph{Direct mode only.} Named numeric vector (gene →
#'   AUC or logFC), sorted descending.
#' @param gsea_dir \emph{RunGSEA directory mode.} Path to the \code{output_dir}
#'   used in \code{\link{RunGSEA}}.
#' @param database Database subfolder name (e.g. \code{"Hallmark"},
#'   \code{"C2"}). \code{NULL} uses the first database found.
#' @param search_terms Keyword search for pathways - same AND/OR logic as
#'   \code{\link{RunSCssGSEA}}: a character vector applies OR; a list of
#'   character vectors applies AND within each element, OR across elements.
#'   Case-insensitive.
#' @param pathways Exact pathway names to plot (alternative to
#'   \code{search_terms}).
#' @param top_n \emph{Direct mode only.} When neither \code{pathways} nor
#'   \code{search_terms} is given, plot the top N by \code{|NES|}.  Default
#'   \code{5}.
#' @param label_levels \emph{RunGSEA directory mode.} Which label.by levels to
#'   include. \code{NULL} = all.
#' @param group_levels \emph{RunGSEA directory mode.} Which group.by (cluster)
#'   levels to include as pages. \code{NULL} = all.
#' @param split_levels \emph{RunGSEA directory mode.} Which split.by (cell type)
#'   levels to include as pages. \code{NULL} = all.
#' @param max_cols \emph{RunGSEA directory mode.} Maximum number of group.by
#'   levels (columns) per page. When the total number of groups exceeds this
#'   value the groups are chunked across multiple pages within the same PDF.
#'   Default \code{6}.
#' @param colors Named vector with \code{"pos"} (positive NES) and \code{"neg"}
#'   (negative NES) color entries. Default red / blue.
#' @param ncol \emph{Direct mode only.} Number of columns in the output layout.
#'   Default \code{1}.
#' @param output_dir Directory to save PDFs. \code{NULL} = no save.
#' @param object_name Label prefix for output file names.
#' @param subset_name Optional subset label.
#'
#' @return \emph{Direct mode:} a single \code{patchwork} object.
#'   \emph{RunGSEA directory mode:} invisibly, a nested list
#'   \code{[[lab]][[cl]][[pw]][[sp]]} of mountain-plot patchwork objects;
#'   PDFs are written to \code{output_dir} (one per lab×cl combination).
#'
#' @seealso \code{\link{RunGSEA}}, \code{\link{RunEnrichment}}
#'
#' @export
PlotGSEAEnrichment <- function(
    # Direct mode
    fgsea_result  = NULL,
    ranked_stats  = NULL,
    # RunGSEA directory mode
    gsea_dir      = NULL,
    database      = NULL,
    # Pathway selection (both modes)
    search_terms  = NULL,
    pathways      = NULL,
    top_n         = 5,
    # RunGSEA mode filters
    label_levels  = NULL,
    group_levels  = NULL,
    split_levels  = NULL,
    # RunGSEA directory mode layout
    max_cols      = 6,
    # Aesthetics
    colors        = c(pos = "#B40426", neg = "#3B4CC0"),
    ncol          = 1,
    output_dir    = NULL,
    object_name   = "",
    subset_name   = ""
) {
  col_pos <- colors[["pos"]] %||% "#B40426"
  col_neg <- colors[["neg"]] %||% "#3B4CC0"

  # ════════════════════════════════════════════════════════════════════════════
  # MODE 1: RunGSEA directory
  # ════════════════════════════════════════════════════════════════════════════
  if (!is.null(gsea_dir)) {
    if (!dir.exists(gsea_dir))
      stop("'gsea_dir' does not exist: ", gsea_dir)
    if (is.null(search_terms) && is.null(pathways))
      stop("Provide 'search_terms' or 'pathways' to select pathways.")

    message("Discovering RunGSEA output files in:\n  ", gsea_dir)
    file_map <- .discover_gsea_files(gsea_dir, database)
    if (length(file_map) == 0)
      stop("No GSEA CSV files found. Check 'gsea_dir' and 'database'.")

    # Filter by label / group / split
    if (!is.null(label_levels))
      file_map <- Filter(function(x) x$lab %in% label_levels, file_map)
    if (!is.null(group_levels))
      file_map <- Filter(function(x) x$cl  %in% group_levels, file_map)
    if (!is.null(split_levels))
      file_map <- Filter(function(x) x$sp  %in% split_levels, file_map)

    # Discover all pathway names to apply search
    all_pws <- unique(unlist(lapply(file_map, function(f) {
      tryCatch({
        df <- utils::read.csv(f$csv_path, stringsAsFactors = FALSE)
        df$pathway
      }, error = function(e) character(0))
    })))

    # Filter pathways
    if (!is.null(search_terms)) {
      matched_pws <- all_pws[.apply_search_terms(all_pws, search_terms)]
      term_label  <- .format_search_terms(search_terms)
    } else {
      matched_pws <- intersect(pathways, all_pws)
      term_label  <- paste(pathways, collapse = " | ")
    }
    if (length(matched_pws) == 0)
      stop("No pathways matched the search. Found ", length(all_pws),
           " total pathways across all CSVs.")

    message("  Matched ", length(matched_pws), " pathway(s): ",
            paste(utils::head(matched_pws, 5), collapse = ", "),
            if (length(matched_pws) > 5) " ..." else "")

    # Layout:
    #   Page    = split.by level (cell type)     - one PDF per cell type
    #   Columns = group.by levels (Male/Female)  - up to max_cols per page
    #   Rows    = matched pathways               - all on the same page
    #
    # When n_groups > max_cols the groups are chunked: each chunk becomes a
    # separate page in the PDF titled "CellType [Lab] - groups 1-6 of 32".

    labs_found <- unique(sapply(file_map, `[[`, "lab"))

    if (!is.null(output_dir))
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    pfx <- paste(c(object_name, subset_name)[nchar(c(object_name, subset_name)) > 0],
                 collapse = "_")

    .safe_name <- function(x) gsub("[/\\\\:*?\"<>|]", "_", x)

    # ── Pre-load all CSVs and de_caches once (avoid re-reading per pathway) ──
    # Each de_cache_{sp}_{lab}.rds contains all groups for that (sp, lab) pair.
    # Reading it once and slicing per group is much faster than N_group re-reads.
    message("  Pre-loading ranked stats...")
    de_store    <- list()   # de_cache path  → full DE data.frame
    ranks_store <- list()   # "<de_cache>\x1f<cl>" -> named numeric vector
    csv_store   <- list()   # csv_path → data.frame

    for (entry in file_map) {
      if (!entry$de_cache %in% names(de_store)) {
        de_store[[entry$de_cache]] <- tryCatch(
          readRDS(entry$de_cache),
          error = function(e) {
            message("  Cannot read de_cache: ", basename(entry$de_cache)); NULL
          }
        )
      }
      rkey <- paste0(entry$de_cache, "\x1f", entry$cl)
      if (!rkey %in% names(ranks_store)) {
        de <- de_store[[entry$de_cache]]
        ranks_store[[rkey]] <- if (!is.null(de))
          tryCatch(
            de |> dplyr::filter(group == entry$cl) |>
                  dplyr::arrange(dplyr::desc(auc)) |>
                  dplyr::select(feature, auc) |>
                  tibble::deframe(),
            error = function(e) NULL
          )
        else NULL
      }
      if (!entry$csv_path %in% names(csv_store)) {
        csv_store[[entry$csv_path]] <- tryCatch(
          utils::read.csv(entry$csv_path, stringsAsFactors = FALSE),
          error = function(e) {
            message("  Cannot read: ", basename(entry$csv_path)); NULL
          }
        )
      }
    }
    rm(de_store)   # free full DE tables now that ranks are extracted

    # ── One PDF per pathway; rows = cell types, columns = comparison groups ──
    all_plots <- list()

    for (pw in matched_pws) {
      pw_safe <- .safe_name(pw)

      for (lab in labs_found) {
        lab_map  <- Filter(function(x) x$lab == lab, file_map)
        sps_lab  <- unique(sapply(lab_map, `[[`, "sp"))   # cell types → rows
        cls_lab  <- unique(sapply(lab_map, `[[`, "cl"))   # groups     → columns
        lab_safe <- .safe_name(lab)

        cl_batches <- split(cls_lab, ceiling(seq_along(cls_lab) / max_cols))
        n_batches  <- length(cl_batches)

        pdf_name <- paste0(
          pfx, if (nchar(pfx)) " " else "",
          pw_safe,
          if (length(labs_found) > 1) paste0(" [", lab_safe, "]") else "",
          " GSEA mountain.pdf"
        )
        pdf_path <- if (!is.null(output_dir))
          file.path(output_dir, pdf_name) else NULL

        if (!is.null(pdf_path)) {
          pdf(pdf_path,
              width  = min(length(cls_lab), max_cols) * 4.5,
              height = length(sps_lab) * 5.5 + 0.8)
          on.exit(grDevices::dev.off(), add = TRUE)
        }

        for (b_idx in seq_along(cl_batches)) {
          cls_batch <- cl_batches[[b_idx]]

          # One column per comparison group; rows = cell types
          page_cols <- lapply(cls_batch, function(cl) {
            sp_plots <- lapply(sps_lab, function(sp) {
              entry <- Filter(function(x) x$sp == sp && x$cl == cl, lab_map)
              if (length(entry) == 0) return(NULL)
              entry <- entry[[1]]

              csv_df <- csv_store[[entry$csv_path]]
              if (is.null(csv_df)) return(NULL)
              rows <- csv_df[csv_df$pathway == pw, , drop = FALSE]
              if (nrow(rows) == 0) return(NULL)

              rkey   <- paste0(entry$de_cache, "\x1f", cl)
              ranks  <- ranks_store[[rkey]]
              if (is.null(ranks) || length(ranks) == 0) return(NULL)

              .one_mountain(ranks, rows[1, ], col_pos, col_neg,
                            title_extra = paste0("\n", sp))
            })
            sp_plots <- Filter(Negate(is.null), sp_plots)
            if (length(sp_plots) == 0) return(NULL)
            patchwork::wrap_plots(sp_plots, ncol = 1)
          })

          page_cols <- Filter(Negate(is.null), page_cols)
          if (length(page_cols) == 0) next

          batch_lbl <- if (n_batches > 1)
            paste0(" - groups ", min(which(cls_lab %in% cls_batch)),
                   "-", max(which(cls_lab %in% cls_batch)),
                   " of ", length(cls_lab))
          else ""

          page <- patchwork::wrap_plots(page_cols, nrow = 1) +
            patchwork::plot_annotation(
              title    = paste0(pw, if (nchar(lab) && lab != "All")
                                      paste0("  [", lab, "]") else "",
                                batch_lbl),
              subtitle = paste0("Rows: ", paste(sps_lab, collapse = " | ")),
              theme    = ggplot2::theme(
                plot.title    = ggplot2::element_text(face = "bold", size = 11),
                plot.subtitle = ggplot2::element_text(size = 8, color = "gray40")
              )
            )

          print(page)
          all_plots[[pw]][[lab]][[b_idx]] <- page
        }  # end batch loop

        if (!is.null(pdf_path)) {
          on.exit(NULL, add = FALSE)
          grDevices::dev.off()
          message("  Saved: ", basename(pdf_path))
        }
      }  # end lab loop
    }  # end pathway loop

    return(invisible(all_plots))
  }

  # ════════════════════════════════════════════════════════════════════════════
  # MODE 2: Direct - fgsea_result + ranked_stats
  # ════════════════════════════════════════════════════════════════════════════
  if (is.null(fgsea_result) || is.null(ranked_stats))
    stop("Provide either 'gsea_dir' (RunGSEA directory mode) or ",
         "both 'fgsea_result' and 'ranked_stats' (direct mode).")

  if (is.character(fgsea_result) && length(fgsea_result) == 1 &&
      file.exists(fgsea_result))
    fgsea_result <- utils::read.csv(fgsea_result, stringsAsFactors = FALSE)
  if (!is.data.frame(fgsea_result))
    stop("'fgsea_result' must be a data frame or a CSV file path.")

  ranked_stats <- sort(ranked_stats, decreasing = TRUE)

  # Select pathways
  if (!is.null(search_terms)) {
    idx <- .apply_search_terms(fgsea_result$pathway, search_terms)
    pws <- fgsea_result$pathway[idx]
  } else if (!is.null(pathways)) {
    pws <- intersect(pathways, fgsea_result$pathway)
  } else {
    fgsea_result <- fgsea_result[order(abs(as.numeric(fgsea_result$NES)),
                                       decreasing = TRUE,
                                       na.last = TRUE), ]
    pws <- utils::head(fgsea_result$pathway, top_n)
  }
  if (length(pws) == 0)
    stop("No pathways selected. Check 'search_terms' or 'pathways'.")

  pw_list <- lapply(pws, function(pw) {
    row <- fgsea_result[fgsea_result$pathway == pw, ][1, ]
    .one_mountain(ranked_stats, row, col_pos, col_neg)
  })
  pw_list <- Filter(Negate(is.null), pw_list)

  combined <- patchwork::wrap_plots(pw_list, ncol = ncol)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    pfx    <- paste(c(object_name, subset_name)[nchar(c(object_name, subset_name)) > 0],
                    collapse = "_")
    pdf(file.path(output_dir, paste0(pfx, if (nchar(pfx)) " " else "",
                                      "GSEA mountain.pdf")),
        width  = ncol * 4.5,
        height = ceiling(length(pw_list) / ncol) * 5.5)
    print(combined)
    grDevices::dev.off()
  }

  print(combined)
  invisible(combined)
}


# =============================================================================
# PlotEnrichment
# =============================================================================

#' Over-Representation Analysis Visualization
#'
#' @description
#' Produces bar, dot, or lollipop plots from ORA enrichment results (enrichR,
#' clusterProfiler, gprofiler2, or any table with term names and p-values).
#' Column names are auto-detected for the most common output formats; override
#' with the \code{*_column} parameters when needed.
#'
#' @param enrichment_df Data frame with ORA results, or a character path to a
#'   CSV file.
#' @param term_column Column with term/pathway names. Auto-detected from common
#'   names (\code{"term"}, \code{"Term"}, \code{"Description"}, etc.).
#' @param padj_column Adjusted p-value column. Auto-detected.
#' @param gene_count_column Column with overlap gene count. Auto-detected.
#' @param gene_ratio_column Column with gene ratio or fold enrichment.
#'   Auto-detected; computed from \code{gene_count_column / background} if
#'   absent.
#' @param group_column Optional column to facet by (e.g. cluster or condition).
#'   Produces one panel per group.
#' @param top_n Maximum number of terms to show per group. Default \code{20}.
#' @param plot_type One of \code{"dot"} (default), \code{"bar"}, or
#'   \code{"lollipop"}.
#' @param x_by What drives the x-axis (or bar length). One of
#'   \code{"gene_ratio"} (default), \code{"gene_count"}, or \code{"-log10padj"}.
#' @param color_by What drives fill/color. One of \code{"padj"} (default) or
#'   \code{"gene_count"}.
#' @param colors Character vector of 2-5 colors for the gradient (low padj →
#'   high padj, i.e. significant → not). Default diverging red-blue.
#' @param max_term_length Maximum characters for term label before wrapping.
#'   Default \code{50}.
#' @param output_dir Directory to save a PDF. \code{NULL} = no save.
#' @param object_name Label prefix for output file names.
#' @param subset_name Optional subset label.
#'
#' @return A \code{ggplot2} or \code{patchwork} object.
#'
#' @export
PlotEnrichment <- function(
    enrichment_df,
    term_column       = NULL,
    padj_column       = NULL,
    gene_count_column = NULL,
    gene_ratio_column = NULL,
    group_column      = NULL,
    top_n             = 20,
    plot_type         = "dot",
    x_by              = "gene_ratio",
    color_by          = "padj",
    colors            = c("#B40426", "#FDDBC7", "#D1E5F0", "#3B4CC0"),
    max_term_length   = 50,
    output_dir        = NULL,
    object_name       = "",
    subset_name       = ""
) {
  plot_type <- match.arg(plot_type, c("dot", "bar", "lollipop"))
  x_by      <- match.arg(x_by,      c("gene_ratio", "gene_count", "-log10padj"))
  color_by  <- match.arg(color_by,  c("padj", "gene_count"))

  # ── Load if path ──────────────────────────────────────────────────────────
  if (is.character(enrichment_df) && length(enrichment_df) == 1 &&
      file.exists(enrichment_df))
    enrichment_df <- utils::read.csv(enrichment_df, stringsAsFactors = FALSE)
  if (!is.data.frame(enrichment_df))
    stop("'enrichment_df' must be a data frame or a CSV file path.")

  # ── Auto-detect columns ───────────────────────────────────────────────────
  term_col  <- term_column %||%
    .detect_col(enrichment_df, c("term","Term","Description","pathway","category",
                                  "Term.name","name","ID"))
  padj_col  <- padj_column %||%
    .detect_col(enrichment_df, c("padj","p.adjust","qvalue","FDR","p_val_adj",
                                   "Adjusted.P.value","adjusted_p_value","P.adj"))
  count_col <- gene_count_column %||%
    .detect_col(enrichment_df, c("gene_count","Count","count","overlap","Overlap",
                                   "Genes","ngenes","n_genes"))
  ratio_col <- gene_ratio_column %||%
    .detect_col(enrichment_df, c("gene_ratio","GeneRatio","RichFactor",
                                   "FoldEnrichment","fold_enrichment","enrichment_ratio"))

  if (is.null(term_col))
    stop("Could not auto-detect a term/pathway name column. ",
         "Set 'term_column' explicitly.")
  if (is.null(padj_col))
    stop("Could not auto-detect an adjusted p-value column. ",
         "Set 'padj_column' explicitly.")

  df <- enrichment_df
  df$._term  <- as.character(df[[term_col]])
  df$._padj  <- pmax(as.numeric(df[[padj_col]]), 1e-300)
  df$._nlp   <- -log10(df$._padj)

  if (!is.null(count_col)) {
    raw <- df[[count_col]]
    # enrichR returns "5/200" format; extract the numerator
    if (is.character(raw) && any(grepl("/", raw, fixed = TRUE)))
      raw <- sub("/.*", "", raw)
    df$._count <- suppressWarnings(as.numeric(raw))
  }
  if (!is.null(ratio_col)) {
    raw <- df[[ratio_col]]
    if (is.character(raw) && any(grepl("/", raw, fixed = TRUE))) {
      parts <- strsplit(raw, "/", fixed = TRUE)
      raw   <- sapply(parts, function(x) as.numeric(x[1]) / as.numeric(x[2]))
    }
    df$._ratio <- suppressWarnings(as.numeric(raw))
  }

  # Compute ratio from count if missing
  if (is.null(ratio_col) && !is.null(count_col)) {
    # Use count directly when ratio unavailable
    df$._ratio <- df$._count
    if (x_by == "gene_ratio") x_by <- "gene_count"
  }

  # X variable
  df$._x <- switch(x_by,
    gene_ratio  = if (!is.null(df$._ratio)) df$._ratio else df$._count,
    gene_count  = if (!is.null(df$._count)) df$._count else df$._nlp,
    `-log10padj`= df$._nlp
  )
  x_label <- switch(x_by,
    gene_ratio  = "Gene Ratio",
    gene_count  = "Gene Count",
    `-log10padj`= expression(-log[10](p[adj]))
  )

  # Color variable
  df$._color <- switch(color_by,
    padj       = df$._padj,
    gene_count = if (!is.null(df$._count)) df$._count else df$._nlp
  )
  color_label <- switch(color_by, padj = "p-adj", gene_count = "Gene Count")

  # Reverse gradient direction for padj (low padj = significant = warmer)
  col_scale <- ggplot2::scale_color_gradientn(
    colors  = if (color_by == "padj") colors else rev(colors),
    name     = color_label,
    trans    = if (color_by == "padj") "log10" else "identity"
  )
  fill_scale <- ggplot2::scale_fill_gradientn(
    colors  = if (color_by == "padj") colors else rev(colors),
    name     = color_label,
    trans    = if (color_by == "padj") "log10" else "identity"
  )

  # ── Build per-group plots ─────────────────────────────────────────────────
  .build_one <- function(g_df, group_label = NULL) {
    # Top N by padj
    g_df <- g_df[order(g_df$._padj), ]
    g_df <- utils::head(g_df, top_n)
    # Wrap long term names
    g_df$._term <- stringr::str_wrap(g_df$._term, width = max_term_length)
    g_df$._term <- factor(g_df$._term, levels = rev(unique(g_df$._term)))

    p <- ggplot2::ggplot(g_df, ggplot2::aes(
      x    = .data[["._x"]],
      y    = .data[["._term"]],
      color = .data[["._color"]],
      fill  = .data[["._color"]]
    ))

    if (plot_type == "dot") {
      size_aes <- if (!is.null(count_col))
        ggplot2::aes(size = .data[["._count"]]) else ggplot2::aes(size = 3)
      p <- p +
        ggplot2::geom_point(size_aes, shape = 21, stroke = 0.4, color = "gray30") +
        ggplot2::scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
        fill_scale

    } else if (plot_type == "bar") {
      p <- p +
        ggplot2::geom_bar(stat = "identity", width = 0.7,
                           color = "gray30", linewidth = 0.3) +
        fill_scale

    } else {  # lollipop
      p <- p +
        ggplot2::geom_segment(
          ggplot2::aes(x = 0, xend = .data[["._x"]],
                       y = .data[["._term"]], yend = .data[["._term"]]),
          color = "gray60", linewidth = 0.5
        ) +
        ggplot2::geom_point(size = 4, shape = 21, stroke = 0.4, color = "gray30") +
        fill_scale
    }

    p <- p +
      ggplot2::labs(
        title = group_label,
        x     = x_label,
        y     = NULL
      ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
      .de_theme() +
      ggplot2::theme(
        legend.position  = "right",
        legend.key.height = ggplot2::unit(0.4, "cm"),
        axis.text.y      = ggplot2::element_text(size = 8)
      )

    p
  }

  if (!is.null(group_column) && group_column %in% colnames(df)) {
    groups   <- split(df, df[[group_column]])
    plot_lst <- mapply(.build_one, groups, names(groups), SIMPLIFY = FALSE)
    out_plot <- patchwork::wrap_plots(plot_lst, ncol = 1)
  } else {
    out_plot <- .build_one(df)
  }

  # ── Save ─────────────────────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    pfx  <- paste(c(object_name, subset_name)[nchar(c(object_name, subset_name)) > 0],
                  collapse = "_")
    n_grp <- if (!is.null(group_column) && group_column %in% colnames(df))
      length(unique(df[[group_column]])) else 1
    pdf(file.path(output_dir, paste0(pfx, " enrichment ", plot_type, ".pdf")),
        width  = 7,
        height = max(4, min(top_n, nrow(df)) * 0.28 + 2) * n_grp)
    print(out_plot)
    grDevices::dev.off()
  }

  print(out_plot)
  invisible(out_plot)
}


# =============================================================================
# RunEnrichment
# =============================================================================

#' Run Over-Representation Analysis and Plot Results
#'
#' @description
#' Runs over-representation analysis (ORA) using \pkg{enrichR} on a gene list
#' or on significant DEGs from a data frame, then automatically calls
#' \code{\link{PlotEnrichment}} to visualize the results.  When
#' \code{group_column} is supplied, analysis is run separately per group.
#'
#' @param features Character vector of gene symbols to test.  Provide either
#'   this or \code{deg_df}, not both.
#' @param deg_df Data frame of DE results (presto or FindMarkers).  Genes are
#'   filtered by \code{fc_cutoff}, \code{padj_cutoff}, and \code{direction}
#'   before testing.
#' @param gene_column Gene name column in \code{deg_df}. Auto-detected if
#'   \code{NULL}.
#' @param fc_column logFC column in \code{deg_df}. Auto-detected if \code{NULL}.
#' @param padj_column Adjusted p-value column. Auto-detected if \code{NULL}.
#' @param group_column Column in \code{deg_df} to loop over (e.g.
#'   \code{"group"}). \code{NULL} = test all significant genes together.
#' @param fc_cutoff Minimum absolute logFC to include a gene. Default
#'   \code{0.25}.
#' @param padj_cutoff Maximum adjusted p-value cutoff. Default \code{0.05}.
#' @param direction Which DEGs to test: \code{"up"} (logFC > 0), \code{"down"}
#'   (logFC < 0), or \code{"both"}. Default \code{"up"}.
#' @param databases Character vector of enrichR database names to query.
#'   See \code{enrichR::listEnrichrDbs()} for all options.  Default:
#'   \code{c("MSigDB_Hallmark_2020", "GO_Biological_Process_2023", "KEGG_2021_Human")}.
#' @param top_n_plot Maximum number of terms to show per database. Default
#'   \code{20}.
#' @param plot_type Type of plot passed to \code{\link{PlotEnrichment}}.
#'   Default \code{"dot"}.
#' @param padj_cutoff_plot Significance cutoff for filtering enrichment results
#'   before plotting. Default \code{0.05}.
#' @param output_dir Directory to save PDFs and result CSVs.  \code{NULL} =
#'   no save, print only.
#' @param object_name Label prefix for output file names.
#' @param subset_name Optional subset label.
#'
#' @return A nested list: \code{[[group]][[database]]} containing enrichR
#'   result data frames.  Invisibly when \code{output_dir} is set.
#'
#' @seealso \code{\link{PlotEnrichment}}, \code{\link{PlotVolcano}}
#'
#' @export
RunEnrichment <- function(
    features          = NULL,
    deg_df            = NULL,
    gene_column       = NULL,
    fc_column         = NULL,
    padj_column       = NULL,
    group_column      = NULL,
    fc_cutoff         = 0.25,
    padj_cutoff       = 0.05,
    direction         = "up",
    databases         = c("MSigDB_Hallmark_2020",
                          "GO_Biological_Process_2023",
                          "KEGG_2021_Human"),
    top_n_plot        = 20,
    plot_type         = "dot",
    padj_cutoff_plot  = 0.05,
    output_dir        = NULL,
    object_name       = "",
    subset_name       = ""
) {
  direction <- match.arg(direction, c("up", "down", "both"))

  if (!requireNamespace("enrichR", quietly = TRUE))
    stop("Package 'enrichR' is required.\n",
         "Install with: install.packages('enrichR')")
  # enrichR's connection setup runs in .onAttach(), which only fires with
  # library(), not requireNamespace(). Warn early if the session is not ready.
  if (is.null(getOption("enrichR.base.address")))
    stop("enrichR is not initialized. Run library(enrichR) once at the top of ",
         "your script before calling RunEnrichment().")

  if (is.null(features) && is.null(deg_df))
    stop("Provide either 'features' (gene list) or 'deg_df' (DEG data frame).")

  # ── Build gene lists per group ────────────────────────────────────────────
  if (!is.null(features)) {
    gene_lists <- list(All = as.character(features))
  } else {
    gene_col  <- gene_column  %||% .detect_col(deg_df, c("feature","gene","Gene","Symbol"))
    fc_col    <- fc_column    %||% .detect_col(deg_df, c("logFC","log2FC","avg_log2FC"))
    padj_col  <- padj_column  %||% .detect_col(deg_df, c("padj","p_val_adj","p.adjust","FDR"))
    if (is.null(gene_col)) stop("Cannot auto-detect gene column. Set 'gene_column'.")

    # Filter for significant genes
    sig <- deg_df
    if (!is.null(padj_col) && padj_col %in% colnames(sig))
      sig <- sig[!is.na(sig[[padj_col]]) & sig[[padj_col]] < padj_cutoff, ]
    if (!is.null(fc_col) && fc_col %in% colnames(sig)) {
      if (direction == "up")
        sig <- sig[sig[[fc_col]] >  fc_cutoff, ]
      else if (direction == "down")
        sig <- sig[sig[[fc_col]] < -fc_cutoff, ]
      else
        sig <- sig[abs(sig[[fc_col]]) > fc_cutoff, ]
    }

    if (!is.null(group_column) && group_column %in% colnames(sig)) {
      gene_lists <- lapply(
        split(sig, sig[[group_column]]),
        function(g) unique(as.character(g[[gene_col]]))
      )
    } else {
      gene_lists <- list(All = unique(as.character(sig[[gene_col]])))
    }
  }

  # ── Run enrichR per group ─────────────────────────────────────────────────
  enrichR::setEnrichrSite("Enrichr")

  all_results <- list()
  pfx <- paste(c(object_name, subset_name)[nchar(c(object_name, subset_name)) > 0],
               collapse = "_")

  for (grp in names(gene_lists)) {
    genes <- gene_lists[[grp]]
    message("\nGroup: ", grp, " - ", length(genes), " gene(s)")

    if (length(genes) == 0) {
      message("  No genes after filtering - skipping.")
      next
    }

    res <- tryCatch(
      enrichR::enrichr(genes, databases),
      error = function(e) {
        message("  enrichR failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(res)) next

    all_results[[grp]] <- res

    for (db in names(res)) {
      df_res <- res[[db]]
      if (is.null(df_res) || nrow(df_res) == 0) next

      # Filter by significance and take top_n
      padj_col_enr <- .detect_col(df_res, c("Adjusted.P.value","p.adjust","padj","FDR"))
      if (!is.null(padj_col_enr)) {
        df_res <- df_res[!is.na(df_res[[padj_col_enr]]) &
                           df_res[[padj_col_enr]] < padj_cutoff_plot, ]
      }
      df_res <- utils::head(df_res, top_n_plot)
      if (nrow(df_res) == 0) {
        message("  ", db, ": no significant terms at padj < ", padj_cutoff_plot)
        next
      }

      message("  ", db, ": ", nrow(df_res), " significant term(s)")

      # Save CSV
      if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        grp_safe <- gsub("[/\\\\:*?\"<>|\\s]", "_", grp)
        db_safe  <- gsub("[/\\\\:*?\"<>|\\s]", "_", db)
        csv_path <- file.path(output_dir,
                               paste0(pfx, if (nchar(pfx)) " " else "",
                                      grp_safe, " ", db_safe, " enrichment.csv"))
        utils::write.csv(df_res, csv_path, row.names = FALSE)
      }

      # Plot
      PlotEnrichment(
        enrichment_df     = df_res,
        term_column       = .detect_col(df_res, c("Term","term","Description","name")),
        padj_column       = padj_col_enr,
        gene_count_column = .detect_col(df_res, c("Overlap","Count","gene_count","overlap")),
        top_n             = top_n_plot,
        plot_type         = plot_type,
        output_dir        = output_dir,
        object_name       = paste0(pfx, if (nchar(pfx)) " " else "", grp, " ", db),
        subset_name       = ""
      )
    }
  }

  message("\nRunEnrichment complete.")
  invisible(all_results)
}
