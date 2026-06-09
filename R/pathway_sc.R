# =============================================================================
# scSidekick - single-cell ssGSEA  (pathway_sc.R)
#
# Exported:
#   RunSCssGSEA()          - per-cell ssGSEA via GSVA, significance testing,
#                            ComplexHeatmap visualization
#
# Internal helpers:
#   .build_gene_sets_sc()  - gene-set preparation (4 input modes)
#   .run_gsva_ssgsea()     - GSVA v1 / v2 version-aware wrapper
#   .sc_significance()     - ANOVA or Kruskal-Wallis per pathway
#   .sc_ssgsea_heatmap()   - ComplexHeatmap mean z-score heatmap
#
# Gene-set input modes
#   1. Named list          gene_sets = list(SetA = c("G1","G2"), ...)
#   2. DEG data frame      deg_df + deg_gene_column + deg_group_column (1 or 2 cols)
#   3. MSigDB library      gene_set_library = "H" / "C2" / "C8" ...
#   4. MSigDB term search  search_terms = c("APOPTOSIS", "CELL_DEATH")   # OR
#                          search_terms = list(c("T_CELL","ACTIVATION"))  # AND
#                          search_terms = list(c("T_CELL","ACTIVATION"), "APOPTOSIS") # (A AND B) OR C
#   Modes 3 and 4 can be combined.
#
# When deg_group_column has TWO columns the two values are joined internally
# with "__NKSEP__" so the first column can be recovered as a heatmap row-split
# category.  This is transparent to the user - just pass
#   deg_group_column = c("celltype", "cytokine")
# and the heatmap rows will automatically be grouped by celltype.
# =============================================================================


# -----------------------------------------------------------------------------
# .apply_search_terms()
# Applies AND / OR search logic to a character vector of gene-set names.
#
# search_terms can be:
#   character vector  → OR across terms  (any term matches)
#   list of elements  → union (OR) across elements; each element is:
#       single string     → simple match
#       character vector  → AND: ALL terms must appear in the same name
#
# All matching is case-insensitive.
#
# Examples:
#   c("APOPTOSIS","CELL_DEATH")               → APOPTOSIS OR CELL_DEATH
#   list(c("T_CELL","ACTIVATION"))            → T_CELL AND ACTIVATION
#   list(c("T_CELL","ACTIVATION"),"APOPTOSIS")→ (T_CELL AND ACTIVATION) OR APOPTOSIS
# -----------------------------------------------------------------------------
.apply_search_terms <- function(gs_names, search_terms) {
  .match_one <- function(terms, names) {
    Reduce(`&`, lapply(terms, function(t)
      grepl(t, names, ignore.case = TRUE, perl = TRUE)
    ))
  }

  if (is.character(search_terms)) {
    # Simple OR - each term is an independent OR element
    Reduce(`|`, lapply(search_terms, function(t) .match_one(t, gs_names)))
  } else if (is.list(search_terms)) {
    # Union of AND-groups
    Reduce(`|`, lapply(search_terms, function(el)
      .match_one(as.character(el), gs_names)
    ))
  } else {
    stop("`search_terms` must be a character vector or a list of character vectors.")
  }
}


# -----------------------------------------------------------------------------
# .format_search_terms()
# Human-readable representation of search logic (for messages and file labels).
# -----------------------------------------------------------------------------
.format_search_terms <- function(search_terms, sep_and = " AND ",
                                  sep_or = " OR ", wrap_and = TRUE) {
  fmt_group <- function(el) {
    el <- toupper(as.character(el))
    if (length(el) == 1) el
    else if (wrap_and) paste0("(", paste(el, collapse = sep_and), ")")
    else paste(el, collapse = sep_and)
  }
  if (is.character(search_terms)) {
    paste(toupper(search_terms), collapse = sep_or)
  } else {
    paste(vapply(search_terms, fmt_group, character(1)), collapse = sep_or)
  }
}


# -----------------------------------------------------------------------------
# .build_gene_sets_sc()
# Returns list(gene_sets = <named list>, feature_cats = <data.frame or NULL>)
# feature_cats has columns "pathway" and "category" - used for heatmap row split
# -----------------------------------------------------------------------------
.build_gene_sets_sc <- function(
    gene_sets            = NULL,
    deg_df               = NULL,
    deg_gene_column         = "feature",
    deg_group_column        = "group",
    deg_fc_column           = "logFC",
    deg_padj_column         = "padj",
    deg_padj_cutoff      = 0.05,
    deg_top_n            = 20,
    gene_set_library     = NULL,
    gene_set_subcategory = NULL,
    search_terms         = NULL,
    species              = "Homo sapiens"
) {

  # ── Mode 1: direct named list ───────────────────────────────────────────────
  if (!is.null(gene_sets)) {
    if (!is.list(gene_sets) || is.null(names(gene_sets)))
      stop("`gene_sets` must be a *named* list of character vectors.")
    return(list(gene_sets = gene_sets, feature_cats = NULL))
  }

  # ── Mode 2: DEG data frame ──────────────────────────────────────────────────
  if (!is.null(deg_df)) {
    missing_cols <- setdiff(c(deg_gene_column, deg_group_column, deg_fc_column),
                            colnames(deg_df))
    if (length(missing_cols) > 0)
      stop("Columns not found in deg_df: ",
           paste(missing_cols, collapse = ", "))

    # Filter by adjusted p-value
    if (!is.null(deg_padj_column) && deg_padj_column %in% colnames(deg_df)) {
      deg_df <- deg_df[!is.na(deg_df[[deg_padj_column]]) &
                         deg_df[[deg_padj_column]] < deg_padj_cutoff, ]
      if (nrow(deg_df) == 0)
        stop("No DEGs remain after filtering at padj < ", deg_padj_cutoff,
             ". Check deg_padj_column and deg_padj_cutoff.")
    }

    # Build group key
    if (length(deg_group_column) == 1) {
      deg_df$.nk_group <- as.character(deg_df[[deg_group_column]])
      two_cols <- FALSE
    } else if (length(deg_group_column) == 2) {
      deg_df$.nk_group <- paste(
        as.character(deg_df[[deg_group_column[1]]]),
        as.character(deg_df[[deg_group_column[2]]]),
        sep = "__NKSEP__"
      )
      two_cols <- TRUE
    } else {
      stop("`deg_group_column` must be a character vector of length 1 or 2.")
    }

    # Top n by FC per group (descending)
    deg_df <- deg_df[order(deg_df[[deg_fc_column]], decreasing = TRUE), ]
    deg_df <- do.call(rbind, lapply(
      split(deg_df, deg_df$.nk_group),
      function(d) utils::head(d, deg_top_n)
    ))

    gs_out <- split(as.character(deg_df[[deg_gene_column]]),
                    f = deg_df$.nk_group)

    # Build feature_cats when two group columns were used
    feature_cats <- if (two_cols) {
      parts <- strsplit(names(gs_out), "__NKSEP__", fixed = TRUE)
      data.frame(
        pathway  = names(gs_out),
        category = vapply(parts, `[`, character(1), 1),
        stringsAsFactors = FALSE
      )
    } else NULL

    message("  ", length(gs_out), " gene sets built from DEG data frame.")
    return(list(gene_sets = gs_out, feature_cats = feature_cats))
  }

  # ── Mode 3 / 4: MSigDB via msigdbr ─────────────────────────────────────────
  if (!is.null(gene_set_library) || !is.null(search_terms)) {
    if (!requireNamespace("msigdbr", quietly = TRUE))
      stop("Package 'msigdbr' is required for MSigDB gene sets.\n",
           "Install with: BiocManager::install('msigdbr')")

    message("  Fetching gene sets from MSigDB",
            if (!is.null(gene_set_library))
              paste0(" (", gene_set_library,
                     if (!is.null(gene_set_subcategory))
                       paste0("/", gene_set_subcategory) else "",
                     ")") else " (all collections)",
            "...")
    m_df <- .msigdbr_get(
      species     = species,
      category    = gene_set_library,
      subcategory = gene_set_subcategory
    )

    if (!is.null(search_terms)) {
      logic_str <- .format_search_terms(search_terms)
      message("  Search logic: ", logic_str)
      keep  <- .apply_search_terms(m_df$gs_name, search_terms)
      m_df  <- m_df[keep, ]
      if (nrow(m_df) == 0)
        stop("No gene sets matched: ", logic_str,
             "\nTip: check spelling, try broader terms, or inspect ",
             "names(msigdbr::msigdbr(...)).")
      message("  ", length(unique(m_df$gs_name)), " gene sets matched.")
    } else {
      message("  ", length(unique(m_df$gs_name)), " gene sets fetched.")
    }

    gs_out <- split(m_df$gene_symbol, m_df$gs_name)
    return(list(gene_sets = gs_out, feature_cats = NULL))
  }

  stop("Must provide one of: `gene_sets`, `deg_df`, `gene_set_library`, ",
       "or `search_terms`.")
}


# -----------------------------------------------------------------------------
# .run_gsva_ssgsea()
# Wraps GSVA::gsva() / GSVA::ssgseaParam() handling both v1 and v2 APIs.
# Returns a pathways × cells matrix.
# -----------------------------------------------------------------------------
.run_gsva_ssgsea <- function(expr_mat, gene_sets, normalize = TRUE,
                              min_size = 5, cores = 1) {
  if (!requireNamespace("GSVA", quietly = TRUE))
    stop("Package 'GSVA' is required.\n",
         "Install with: BiocManager::install('GSVA')")

  # BiocParallel backend - progressbar = TRUE gives per-gene-set progress in GSVA v2
  BPPARAM <- tryCatch({
    if (!requireNamespace("BiocParallel", quietly = TRUE))
      return(NULL)
    if (cores > 1) {
      if (.Platform$OS.type == "windows")
        BiocParallel::SnowParam(cores, progressbar = TRUE)
      else
        BiocParallel::MulticoreParam(cores, progressbar = TRUE)
    } else {
      BiocParallel::SerialParam(progressbar = TRUE)
    }
  }, error = function(e) NULL)

  gsva_v <- utils::packageVersion("GSVA")

  if (gsva_v >= "1.46") {
    # GSVA v2 API (Bioconductor >= 3.15)
    param <- GSVA::ssgseaParam(
      exprData  = expr_mat,
      geneSets  = gene_sets,
      normalize = normalize,
      minSize   = min_size
    )
    if (!is.null(BPPARAM)) {
      GSVA::gsva(param, BPPARAM = BPPARAM, verbose = FALSE)
    } else {
      GSVA::gsva(param, verbose = FALSE)
    }
  } else {
    # GSVA v1 legacy API
    GSVA::gsva(
      as.matrix(expr_mat),
      gene_sets,
      method      = "ssgsea",
      ssgsea.norm = normalize,
      min.sz      = min_size,
      parallel.sz = cores,
      verbose     = FALSE
    )
  }
}


# -----------------------------------------------------------------------------
# .run_ucell()
# UCell backend: U statistic (AUC of top-maxRank gene recovery curve).
# ~15-30x faster than GSVA on large scRNA-seq objects because:
#   - Only partially sorts each cell (top maxRank genes, default 1500)
#   - Chunks cells internally so all gene sets are scored together per chunk
# Returns pathways × cells matrix (same convention as .run_gsva_ssgsea).
# -----------------------------------------------------------------------------
.run_ucell <- function(expr_mat, gene_sets, max_rank = 1500L,
                        chunk_size = 1000L, cores = 1L) {
  if (!requireNamespace("UCell", quietly = TRUE))
    stop("Package 'UCell' is required for method = 'ucell'.\n",
         "Install with: BiocManager::install('UCell')")

  BPPARAM <- tryCatch({
    if (!requireNamespace("BiocParallel", quietly = TRUE)) return(NULL)
    if (cores > 1L) {
      if (.Platform$OS.type == "windows")
        BiocParallel::SnowParam(cores, progressbar = TRUE)
      else
        BiocParallel::MulticoreParam(cores, progressbar = TRUE)
    } else {
      BiocParallel::SerialParam(progressbar = TRUE)
    }
  }, error = function(e) NULL)

  args <- list(
    matrix     = expr_mat,
    features   = gene_sets,
    maxRank    = as.integer(max_rank),
    chunk.size = as.integer(chunk_size),
    ncores     = as.integer(cores)
  )
  if (!is.null(BPPARAM)) args$BPPARAM <- BPPARAM

  # UCell returns cells × pathways - transpose to pathways × cells (GSVA convention)
  t(do.call(UCell::ScoreSignatures_UCell, args))
}


# -----------------------------------------------------------------------------
# .sc_significance()
# One-way ANOVA or Kruskal-Wallis per pathway.
# Returns data.frame ordered by p_adj: pathway, statistic, p_value, p_adj
# -----------------------------------------------------------------------------
.sc_significance <- function(scores_mat, group_vec, fit = "ANOVA") {
  fit <- match.arg(fit, c("ANOVA", "Kruskal"))

  res <- do.call(rbind, lapply(rownames(scores_mat), function(pw) {
    sc <- as.numeric(scores_mat[pw, ])
    df <- data.frame(score = sc, group = factor(group_vec))

    out <- if (fit == "ANOVA") {
      aov_s <- tryCatch(
        summary(stats::aov(score ~ group, data = df))[[1]],
        error = function(e) NULL
      )
      if (is.null(aov_s))
        return(data.frame(pathway = pw, statistic = NA_real_,
                          p_value = NA_real_, stringsAsFactors = FALSE))
      data.frame(pathway   = pw,
                 statistic = aov_s["group", "F value"],
                 p_value   = aov_s["group", "Pr(>F)"],
                 stringsAsFactors = FALSE)
    } else {
      kw <- tryCatch(
        stats::kruskal.test(score ~ group, data = df),
        error = function(e) NULL
      )
      if (is.null(kw))
        return(data.frame(pathway = pw, statistic = NA_real_,
                          p_value = NA_real_, stringsAsFactors = FALSE))
      data.frame(pathway   = pw,
                 statistic = as.numeric(kw$statistic),
                 p_value   = kw$p.value,
                 stringsAsFactors = FALSE)
    }
    out
  }))

  res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
  res[order(res$p_adj, na.last = TRUE), ]
}


# -----------------------------------------------------------------------------
# .sc_ssgsea_heatmap()
# ComplexHeatmap of row-z-scored mean ssGSEA enrichment scores per group.
# Columns = group.by levels; optional column split by split.by.
# Row split by feature_cats$category if provided.
# heatmap_params merged via modifyList - any Heatmap() arg can be overridden.
# split_levels preserves user-defined factor order for split.by columns.
# -----------------------------------------------------------------------------
.sc_ssgsea_heatmap <- function(scores_mat, meta_df, group_by,
                                split_by       = NULL,
                                split_levels   = NULL,
                                group_levels   = NULL,
                                feature_cats   = NULL,
                                show_pws,
                                group_colors   = NULL,
                                split_colors   = NULL,
                                heatmap_params = list(),
                                heatmap_colors = NULL,
                                pdf_path) {

  pws_use <- show_pws[show_pws %in% rownames(scores_mat)]
  if (length(pws_use) == 0) {
    message("  No pathways to display in heatmap - skipping.")
    return(invisible(NULL))
  }

  cells_use <- intersect(colnames(scores_mat), rownames(meta_df))
  if (length(cells_use) == 0) {
    message("  No overlapping cells between scores and metadata - skipping heatmap.")
    return(invisible(NULL))
  }

  sc_sub <- scores_mat[pws_use, cells_use, drop = FALSE]
  meta   <- meta_df[cells_use, , drop = FALSE]

  # ── Build mean matrix ───────────────────────────────────────────────────────
  has_split <- !is.null(split_by) && split_by %in% colnames(meta)

  if (has_split) {
    grp <- as.character(meta[[group_by]])
    spl <- as.character(meta[[split_by]])
    combos <- unique(data.frame(grp = grp, spl = spl, stringsAsFactors = FALSE))

    # Order split columns by user-defined levels; fall back to first-appearance
    spl_lvl_ord <- if (!is.null(split_levels) && length(split_levels) > 0)
      split_levels else unique(spl)
    # Order groups within each split by user-defined levels; fall back to first-appearance
    grp_lvl_ord <- if (!is.null(group_levels) && length(group_levels) > 0)
      group_levels else unique(grp)

    combos$spl_f <- factor(combos$spl, levels = spl_lvl_ord)
    combos$grp_f <- factor(combos$grp, levels = grp_lvl_ord)
    combos <- combos[order(combos$spl_f, combos$grp_f), ]

    mean_mat <- do.call(cbind, lapply(seq_len(nrow(combos)), function(i) {
      idx <- which(grp == combos$grp[i] & spl == combos$spl[i])
      if (length(idx) == 0) return(rep(NA_real_, length(pws_use)))
      rowMeans(sc_sub[, idx, drop = FALSE], na.rm = TRUE)
    }))
    colnames(mean_mat) <- paste(combos$grp, combos$spl, sep = "\n")
    rownames(mean_mat) <- pws_use
    col_split <- combos$spl_f   # factor with user-defined level order

    ha_colors <- list()
    if (!is.null(group_colors)) ha_colors[[group_by]] <- group_colors
    top_annot <- ComplexHeatmap::HeatmapAnnotation(
      df  = stats::setNames(as.data.frame(combos$grp), group_by),
      col = ha_colors,
      show_annotation_name = TRUE,
      annotation_name_side = "left",
      annotation_name_gp   = grid::gpar(fontsize = 8)
    )
  } else {
    # Respect user-defined group levels; keep only those present in data
    present_grps <- unique(as.character(meta[[group_by]]))
    grps <- if (!is.null(group_levels) && length(group_levels) > 0)
      group_levels[group_levels %in% present_grps]
    else
      present_grps   # first-appearance order, not sorted
    mean_mat <- do.call(cbind, lapply(grps, function(g) {
      idx <- which(as.character(meta[[group_by]]) == g)
      rowMeans(sc_sub[, idx, drop = FALSE], na.rm = TRUE)
    }))
    colnames(mean_mat) <- grps
    rownames(mean_mat) <- pws_use
    col_split <- NULL

    ha_colors <- list()
    if (!is.null(group_colors)) ha_colors[[group_by]] <- group_colors[grps]
    top_annot <- if (length(ha_colors) > 0) {
      ComplexHeatmap::HeatmapAnnotation(
        df  = stats::setNames(as.data.frame(grps), group_by),
        col = ha_colors,
        show_annotation_name = TRUE,
        annotation_name_side = "left",
        annotation_name_gp   = grid::gpar(fontsize = 8)
      )
    } else NULL
  }

  # ── Row z-score ─────────────────────────────────────────────────────────────
  z_mat <- t(scale(t(mean_mat)))
  z_mat[is.nan(z_mat) | is.na(z_mat)] <- 0

  # ── Row split ───────────────────────────────────────────────────────────────
  row_split <- if (!is.null(feature_cats) && nrow(feature_cats) > 0) {
    cats <- feature_cats$category[match(pws_use, feature_cats$pathway)]
    factor(cats, levels = unique(feature_cats$category))
  } else NULL

  # ── Auto-size PDF ───────────────────────────────────────────────────────────
  n_rows  <- nrow(z_mat)
  n_cols  <- ncol(z_mat)
  cell_h  <- 12   # points per row
  cell_w  <- max(25, min(55, 400 / max(n_cols, 1)))
  pdf_h   <- max(4,  min(80, n_rows * cell_h / 72 + 3.5))
  pdf_w   <- max(6,  min(30, n_cols * cell_w / 72 + 5.5))

  col_fun <- if (is.null(heatmap_colors)) {
    circlize::colorRamp2(c(-2, -1, 0, 1, 2),
                         c("#007dd1", "#b3d9f5", "white", "#f5c08a", "#ab3000"))
  } else {
    heatmap_colors
  }

  # ── Build heatmap via modifyList so user params override gracefully ──────────
  default_args <- list(
    z_mat,
    name              = "Z-score",
    col               = col_fun,
    top_annotation    = top_annot,
    row_split         = row_split,
    column_split      = col_split,
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_side      = "left",              # default: left avoids legend overlap
    show_row_dend       = FALSE,               # default: cleaner without dendrogram
    row_names_max_width = grid::unit(15, "cm"),# default: cap long names at 15 cm
    row_names_gp        = grid::gpar(fontsize = max(6, min(10, 120 / n_rows))),
    column_names_gp   = grid::gpar(fontsize = 9),
    row_title_gp      = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title_gp   = grid::gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(
      title          = "Row Z-score",
      title_gp       = grid::gpar(fontsize = 9, fontface = "bold"),
      labels_gp      = grid::gpar(fontsize = 8),
      legend_height  = grid::unit(3, "cm")
    )
  )
  ht_args <- utils::modifyList(default_args, heatmap_params)
  ht      <- do.call(ComplexHeatmap::Heatmap, ht_args)

  grDevices::pdf(pdf_path, width = pdf_w, height = pdf_h)
  ComplexHeatmap::draw(ht,
                       merge_legends           = TRUE,
                       heatmap_legend_side     = "right",
                       annotation_legend_side  = "right")
  grDevices::dev.off()
  message("  Saved: ", pdf_path)
  invisible(ht)
}


# -----------------------------------------------------------------------------
# .sc_ssgsea_boxplots()
# One boxplot per pathway, all printed to a single multi-page PDF.
# x-axis = split.by (if set); facets = group.by; else x = group.by (no facets).
# -----------------------------------------------------------------------------
.sc_ssgsea_boxplots <- function(scores_df, group_by, split_by = NULL,
                                 group_levels = NULL, split_levels = NULL,
                                 show_pws, group_colors = NULL,
                                 split_colors = NULL, add_pvalues = TRUE,
                                 pdf_path) {
  pws_use <- show_pws[show_pws %in% colnames(scores_df)]
  if (length(pws_use) == 0) {
    message("  No pathways to plot - skipping boxplots.")
    return(invisible(NULL))
  }

  has_split   <- !is.null(split_by) && split_by %in% colnames(scores_df)
  has_ggpubr  <- requireNamespace("ggpubr", quietly = TRUE)

  if (add_pvalues && !has_ggpubr)
    message("  Note: ggpubr not installed - p-values skipped. ",
            "Install with: install.packages('ggpubr')")

  # Auto PDF dimensions
  n_facets <- if (!is.null(group_levels)) length(group_levels)
              else length(unique(as.character(scores_df[[group_by]])))
  pdf_w    <- if (has_split) min(24, max(7, n_facets * 3.0 + 3)) else 7
  pdf_h    <- 5.5

  # User-preferred theme: clean, black axes, white strip background
  .bx_theme <-
    ggplot2::theme(
      panel.grid.major    = ggplot2::element_blank(),
      panel.grid.minor    = ggplot2::element_blank(),
      panel.background    = ggplot2::element_blank(),
      axis.line           = ggplot2::element_line(color = "black", linewidth = 0.4),
      strip.background    = ggplot2::element_rect(fill = "white", color = "black",
                                                   linewidth = 0.4),
      strip.text          = ggplot2::element_text(face = "bold", color = "black",
                                                   size = 8),
      plot.title          = ggplot2::element_text(face = "bold", color = "black",
                                                   hjust = 0.5, size = 8),
      axis.text.y         = ggplot2::element_text(color = "black", size = 6),
      axis.text.x         = ggplot2::element_text(color = "black", angle = 25,
                                                   hjust = 1, size = 6.5,
                                                   face = "bold"),
      axis.title          = ggplot2::element_text(size = 8, face = "bold"),
      legend.title        = ggplot2::element_text(size = 8, face = "bold"),
      legend.text         = ggplot2::element_text(size = 7),
      plot.margin         = ggplot2::unit(c(0.3, 0.5, 0.2, 0.5), "cm")
    )

  grDevices::pdf(pdf_path, width = pdf_w, height = pdf_h)
  on.exit(grDevices::dev.off())

  for (pw in pws_use) {
    plot_df <- data.frame(
      score = as.numeric(scores_df[[pw]]),
      group = as.character(scores_df[[group_by]]),
      stringsAsFactors = FALSE
    )
    if (!is.null(group_levels))
      plot_df$group <- factor(plot_df$group, levels = group_levels)
    else
      plot_df$group <- factor(plot_df$group)

    # Determine x-axis variable and its levels (for pairwise comparisons)
    if (has_split) {
      plot_df$split_var <- as.character(scores_df[[split_by]])
      if (!is.null(split_levels))
        plot_df$split_var <- factor(plot_df$split_var, levels = split_levels)
      else
        plot_df$split_var <- factor(plot_df$split_var)

      x_lvls   <- levels(plot_df$split_var)
      fill_col <- "split_var"
      x_lab    <- split_by
      fill_vals <- split_colors

      p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(x = split_var, y = score, fill = split_var)
          ) +
          ggplot2::scale_y_continuous(
            limits = NULL,
            expand = ggplot2::expansion(mult = 0.2, add = 0)
          ) +
          ggplot2::geom_boxplot(color = "black", width = 0.5, lwd = 0.3,
                                 outlier.size = 0.3, outlier.alpha = 0.3) +
          ggplot2::geom_dotplot(binaxis = "y", stackdir = "center",
                                 position = ggplot2::position_dodge(1),
                                 dotsize = 0.3) +
          ggplot2::facet_wrap(~ group, nrow = 1, scales = "free") +
          ggplot2::labs(title = pw, x = NULL,
                        y = "ssGSEA enrichment score", fill = split_by) +
          .bx_theme
    } else {
      x_lvls    <- levels(plot_df$group)
      fill_vals <- group_colors

      p <- ggplot2::ggplot(
            plot_df,
            ggplot2::aes(x = group, y = score, fill = group)
          ) +
          ggplot2::scale_y_continuous(
            limits = NULL,
            expand = ggplot2::expansion(mult = 0.2, add = 0)
          ) +
          ggplot2::geom_boxplot(color = "black", width = 0.5, lwd = 0.3,
                                 outlier.size = 0.3, outlier.alpha = 0.3) +
          ggplot2::geom_dotplot(binaxis = "y", stackdir = "center",
                                 position = ggplot2::position_dodge(1),
                                 dotsize = 0.3) +
          ggplot2::labs(title = pw, x = NULL,
                        y = "ssGSEA enrichment score", fill = group_by) +
          .bx_theme
    }

    # Apply fill colors
    if (!is.null(fill_vals))
      p <- p + ggplot2::scale_fill_manual(values = fill_vals)

    # p-values via ggpubr (graceful skip if not installed or only 1 level)
    if (add_pvalues && has_ggpubr && length(x_lvls) >= 2) {
      my_comparisons <- utils::combn(x_lvls, 2, simplify = FALSE)
      p <- p + ggpubr::stat_compare_means(
        comparisons = my_comparisons,
        method      = "t.test",
        label       = "p.format",
        size        = 2,
        vjust       = 0.5
      )
    }

    print(p)
  }

  message("  Saved: ", pdf_path, " (", length(pws_use), " pathway(s))")
  invisible(NULL)
}


# =============================================================================
# RunSCssGSEA - main exported function
# =============================================================================

#' Single-cell ssGSEA with significance testing and heatmap visualization
#'
#' Runs per-cell single-sample Gene Set Enrichment Analysis (ssGSEA) on a
#' Seurat object using \pkg{GSVA} (compatible with both GSVA v1 and v2 APIs),
#' computes group-level significance (ANOVA or Kruskal-Wallis), saves enrichment
#' scores and significance tables as CSVs, and produces a ComplexHeatmap PDF of
#' mean z-scored enrichment per group.
#'
#' Unlike pseudobulk ssGSEA (see [RunGSEA_pseudobulk()]), this function scores
#' every cell individually rather than aggregating by sample first. It is
#' appropriate when you have no meaningful biological replicates (sample IDs),
#' or when you want to use the single-cell resolution for visualization.
#' For datasets with proper replicates, pseudobulk ssGSEA is statistically
#' more rigorous.
#'
#' @section Gene-set input (choose exactly one mode):
#'
#' **Mode 1 - pre-built named list**
#' ```r
#' gene_sets = list(SetA = c("GENE1", "GENE2"), SetB = c("GENE3"))
#' ```
#'
#' **Mode 2 - DEG data frame**\cr
#' Converts a marker/DEG table into one gene set per group. The `deg_group_column`
#' argument can name **one column** (simple group, e.g. `"celltype"`) or **two
#' columns** (e.g. `c("celltype", "cytokine")`). When two columns are supplied,
#' gene set names are built as `"celltype__NKSEP__cytokine"` internally and the
#' heatmap rows are automatically split by the first column (celltype).
#'
#' **Mode 3 - MSigDB collection**
#' ```r
#' gene_set_library = "H"                        # Hallmark
#' gene_set_library = "C2", gene_set_subcategory = "CP:KEGG"
#' ```
#'
#' **Mode 4 - MSigDB term search** (fetch all gene sets matching a pattern)
#' ```r
#' search_terms = c("CELL_DEATH", "APOPTOSIS")
#' # optionally restrict to a collection:
#' gene_set_library = "C2", search_terms = "CELL_DEATH"
#' ```
#'
#' @param seurat_object A Seurat object.
#' @param group.by Character. Metadata column used as the primary grouping
#'   variable for the heatmap columns and significance testing
#'   (e.g. `"Assignment"`, `"celltype"`).
#' @param split.by Character or `NULL`. Optional secondary metadata column for
#'   splitting heatmap columns (e.g. `"Health"`, `"Condition"`). When set,
#'   significance is also tested **within each split level** in addition to
#'   the pooled test.
#' @param subset_cells Logical. If `TRUE`, subset the Seurat object to
#'   `subset_values` before scoring. Default `FALSE`.
#' @param subset_by Character or `NULL`. Metadata column used for subsetting
#'   (required when `subset_cells = TRUE`). `NULL` treats `subset_values` as
#'   cell barcodes directly.
#' @param subset_values Character vector. Values to retain in `subset_by`, or
#'   cell barcodes when `subset_by = NULL`.
#'
#' @param gene_sets Named list of gene vectors (Mode 1). `NULL` if using another
#'   mode.
#' @param deg_df Data frame of DEG results (Mode 2). `NULL` otherwise.
#'   Compatible with standard outputs from [presto::wilcoxauc()], edgeR, or any
#'   table with gene, group, and log-fold-change columns.
#' @param deg_gene_column Character. Column in `deg_df` holding gene symbols.
#'   Default `"feature"`.
#' @param deg_group_column Character scalar **or** length-2 character vector.
#'   Column(s) in `deg_df` defining the gene-set groups. One column → one gene
#'   set per unique value. Two columns → gene sets named
#'   `"col1__NKSEP__col2"`; heatmap rows split by `col1` automatically.
#'   Default `"group"`.
#' @param deg_fc_column Character. Log-fold-change column for selecting top genes.
#'   Default `"logFC"`.
#' @param deg_padj_column Character or `NULL`. Adjusted p-value column for pre-
#'   filtering. `NULL` skips filtering. Default `"padj"`.
#' @param deg_padj_cutoff Numeric. Adjusted p-value threshold for DEG
#'   pre-filtering. Default `0.05`.
#' @param deg_top_n Integer. Top N genes per group (by `deg_fc_column`) to include
#'   in each gene set. Default `20`.
#'
#' @param gene_set_library Character or `NULL`. MSigDB collection code (Mode 3
#'   or 4). Common options:
#'   \itemize{
#'     \item `"H"` - Hallmark (50 coherent biological processes)
#'     \item `"C2"` - Curated: canonical pathways (combine with
#'       `gene_set_subcategory = "CP:KEGG"`, `"CP:REACTOME"`,
#'       `"CP:WIKIPATHWAYS"`, `"CP:BIOCARTA"`, `"CP:PID"`)
#'     \item `"C5"` - Gene Ontology (`subcategory = "GO:BP"`,
#'       `"GO:MF"`, `"GO:CC"`)
#'     \item `"C6"` - Oncogenic signatures
#'     \item `"C7"` - Immunologic signatures (IMMUNESIGDB)
#'     \item `"C8"` - Cell-type signatures
#'   }
#'   `NULL` fetches all collections (slow; use with `search_terms`).
#' @param gene_set_subcategory Character or `NULL`. MSigDB sub-collection code.
#'   Default `NULL` (entire collection).
#' @param search_terms Search filter applied to MSigDB gene-set names (case-
#'   insensitive, regex). Accepts two forms:
#'
#'   **Character vector - OR search.** Any gene set whose name matches at least
#'   one term is kept: \code{search_terms = c("APOPTOSIS", "CELL_DEATH")}.
#'
#'   **List of character vectors - AND/OR search.** Each list element is an
#'   AND-group: all terms in that element must appear in the same pathway name
#'   (order and separation do not matter). Results are the union (OR) across
#'   AND-groups. Examples:
#'   \itemize{
#'     \item \code{list(c("T_CELL", "ACTIVATION"))} - both words in same name
#'     \item \code{list(c("T_CELL", "ACTIVATION"), "APOPTOSIS")} - previous OR APOPTOSIS
#'     \item \code{list(c("T_CELL","ACTIVATION"), c("B_CELL","DIFFERENTIATION"))}
#'       - two independent AND-groups, union of results
#'   }
#'   \code{NULL} (default) keeps all gene sets in the fetched collection.
#' @param species Character. Species for MSigDB fetching. Passed to
#'   [msigdbr::msigdbr()]. Default `"Homo sapiens"`. Use `"Mus musculus"` for
#'   mouse datasets. See `msigdbr::msigdbr_species()` for all supported species.
#'
#' @param method Character. Scoring backend. `"gsva"` (default) uses the GSVA
#'   package's ssGSEA algorithm; `"ucell"` uses the UCell package's rank-based
#'   U statistic, which is faster and memory-efficient for large datasets but
#'   does not produce normalized enrichment scores.
#' @param downsample Integer or `NULL`. Maximum cells per `group.by` identity
#'   to use for scoring. ssGSEA is \eqn{O(n_{\text{cells}})} - downsampling to
#'   500-2000 cells per group dramatically reduces runtime with minimal loss of
#'   group-level signal. `NULL` uses all cells. Default `NULL`.
#' @param min.size Integer. Minimum gene-set size after intersecting with the
#'   expression matrix. Gene sets below this threshold are silently dropped.
#'   Default `5`.
#' @param ssgsea.norm Logical. Normalize ssGSEA enrichment scores to the
#'   `[0, 1]` range? Default `TRUE` (matches `escape::enrichIt` default).
#'   Ignored when `method = "ucell"`.
#' @param ucell_max_rank Integer. UCell only: only the top `ucell_max_rank`
#'   genes per cell (by expression rank) are considered when computing the
#'   U statistic. Lower values reduce memory use. Default `1500L`.
#' @param cores Integer. Number of parallel cores. For `method = "gsva"`,
#'   requires \pkg{BiocParallel}; falls back to serial if unavailable.
#'   For `method = "ucell"`, sets the number of UCell threads. Default `1`.
#' @param chunk_size Integer. UCell only: number of cells processed per
#'   parallel chunk. Tune to balance memory and speed. Default `1000L`.
#' @param resume Logical. If `TRUE` and a combined-scores CSV checkpoint from a
#'   previous run exists in `output_dir`, load it and skip scoring. Useful for
#'   re-running plots after adjusting significance thresholds. Default `FALSE`.
#'
#' @param fit Character. Significance test applied per pathway:
#'   `"ANOVA"` (default, parametric F-test) or `"Kruskal"` (non-parametric
#'   Kruskal-Wallis). Both use BH adjustment for multiple testing.
#' @param sig_group_by Character or `NULL`. Metadata column to use as the
#'   grouping factor for significance testing. `NULL` (default) uses `group.by`.
#'   Override when your display grouping (`group.by`) differs from your
#'   statistical grouping (e.g. display by fine cell type, test by broad
#'   disease category).
#' @param p_cutoff Numeric. BH-adjusted p-value threshold for selecting
#'   pathways to display in the heatmap. Default `0.05`.
#' @param show_only_significant Logical. If `TRUE`, heatmap shows only
#'   pathways with `p_adj < p_cutoff`. Falls back to top 30 by statistic if no
#'   pathways pass the threshold. If `FALSE`, all pathways are shown. Default
#'   `TRUE`.
#' @param group_colors Named character vector of colors for `group.by` levels.
#'   `NULL` (default) auto-generates colors from the Nour palette.
#' @param split_colors Named character vector of colors for `split.by` levels.
#'   `NULL` (default) auto-generates colors from the Nour palette. Ignored
#'   when `split.by = NULL`.
#' @param heatmap_params Named list of additional arguments forwarded to
#'   [ComplexHeatmap::Heatmap()]. Overrides internal defaults for row label
#'   placement, dendrogram display, and row-name width.
#' @param heatmap_colors A `circlize::colorRamp2` color function for the
#'   heatmap fill scale. `NULL` (default) uses a 5-stop blue-white-red
#'   diverging palette at breaks `c(-2, -1, 0, 1, 2)` (row z-scores):
#'   \preformatted{
#'   circlize::colorRamp2(c(-2, -1, 0, 1, 2),
#'                        c("#007dd1", "#b3d9f5", "white", "#f5c08a", "#ab3000"))
#'   }
#'   Supply any `colorRamp2` object to override. The breaks you set are used
#'   as-is - for row z-score heatmaps, symmetric breaks around 0 are
#'   recommended. Examples:
#'   \itemize{
#'     \item \strong{Classic RdBu (3-stop):}
#'       `circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))`
#'     \item \strong{Purple-green (PRGn):}
#'       `circlize::colorRamp2(c(-2, 0, 2), RColorBrewer::brewer.pal(3, "PRGn"))`
#'     \item \strong{Viridis plasma (for raw non-z-scored scores):}
#'       `circlize::colorRamp2(seq(0, 1, length.out = 9), viridis::plasma(9))`
#'     \item \strong{Tighter scale for subtle signals:}
#'       `circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#007dd1", "#b3d9f5", "white", "#f5c08a", "#ab3000"))`
#'   }
#' @param add_boxplots Logical. If `TRUE` (default), generate per-pathway
#'   boxplot pages showing per-cell score distributions across `group.by`
#'   groups (and split by `split.by` if set).
#' @param add_pvalues Logical. If `TRUE` (default), add BH-adjusted p-value
#'   brackets to boxplots using `ggpubr::stat_compare_means()`. Requires
#'   `add_boxplots = TRUE`.
#'
#' @param add_to_object Logical. If `FALSE` (default), returns the enrichment
#'   scores data frame invisibly and does not modify the Seurat object (scores
#'   are saved to CSV only - recommended for large datasets to avoid metadata
#'   bloat). If `TRUE`, adds one metadata column per gene set to `seurat_object@meta.data`
#'   and returns the modified Seurat object. When `downsample` is also set, cells
#'   not included in the downsampled scoring set receive `NA`.
#' @param output_dir Character. Directory for output files (created if absent).
#' @param object_name Character. First component of output file name prefix.
#'   Default `"Analysis"`.
#' @param subset_name Character. Second component of output file name prefix.
#'   Default `""`.
#' @param caffeinate Logical. If `TRUE`, prevents the Mac from sleeping during
#'   the run using the `caffeinate` system command. Default `FALSE`.
#'
#' @return
#' If `add_to_object = FALSE` (default): invisibly returns a data frame of
#' enrichment scores (cells × gene sets) with additional columns for
#' `group.by` and `split.by` identity.
#'
#' If `add_to_object = TRUE`: invisibly returns the Seurat object with one
#' metadata column per gene set.
#'
#' **Files saved to `output_dir`:**
#' \itemize{
#'   \item `<prefix> ssGSEA_scores.csv` - per-cell enrichment scores + group
#'     labels
#'   \item `<prefix> ssGSEA_significance.csv` - per-pathway significance
#'     (pooled across all cells)
#'   \item `<prefix> ssGSEA_significance_by_<split.by>.csv` - per-pathway
#'     significance within each `split.by` level (only when `split.by` is set)
#'   \item `<prefix> ssGSEA heatmap.pdf` - ComplexHeatmap of mean z-scored
#'     enrichment per group
#' }
#'
#' @export
RunSCssGSEA <- function(
    seurat_object,
    group.by,
    split.by               = NULL,

    # Optional cell subsetting (applied before any analysis)
    subset_cells             = FALSE,
    subset_by              = NULL,    # metadata column; NULL → subset_values are barcodes
    subset_values          = NULL,    # factor levels to keep, OR cell barcodes

    # Gene set input (pick one mode)
    gene_sets              = NULL,

    deg_df                 = NULL,
    deg_gene_column           = "feature",
    deg_group_column          = "group",
    deg_fc_column             = "logFC",
    deg_padj_column           = "padj",
    deg_padj_cutoff        = 0.05,
    deg_top_n              = 20,

    gene_set_library       = NULL,
    gene_set_subcategory   = NULL,
    search_terms           = NULL,
    species                = "Homo sapiens",

    # Scoring parameters
    method                 = "gsva",  # "gsva" (ssGSEA, default) or "ucell" (U statistic, faster)
    downsample             = NULL,
    min.size               = 5,
    ssgsea.norm            = TRUE,    # GSVA only: normalize scores to [0,1]
    ucell_max_rank         = 1500L,   # UCell only: consider only top N ranked genes per cell
    cores                  = 1,
    chunk_size             = 1000L,   # UCell only: cells per parallel chunk
    resume                 = FALSE,   # load existing combined-scores checkpoint if found

    # Significance testing
    fit                    = "ANOVA",
    sig_group_by           = NULL,
    p_cutoff               = 0.05,
    show_only_significant  = TRUE,

    # Colors (NULL = auto-generated from Nour palette)
    group_colors           = NULL,
    split_colors           = NULL,

    # Heatmap
    heatmap_params         = list(row_names_side      = "left",
                                  show_row_dend       = FALSE,
                                  row_names_max_width = grid::unit(15, "cm")),
    heatmap_colors         = NULL,

    # Boxplots
    add_boxplots           = TRUE,
    add_pvalues            = TRUE,    # add t-test p-values via ggpubr

    # Output
    add_to_object             = FALSE,
    output_dir,
    object_name               = "Analysis",
    subset_name               = "",
    caffeinate                = FALSE
) {

  if (caffeinate) { .caff <- .nk_caffeinate(); on.exit(.nk_decaffeinate(.caff), add = TRUE) }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  file_prefix <- trimws(paste(object_name, subset_name))
  fit         <- match.arg(fit, c("ANOVA", "Kruskal"))

  # ── 0. Optional cell subsetting ─────────────────────────────────────────────
  if (isTRUE(subset_cells)) {
    if (is.null(subset_values))
      stop("`subset_values` must be supplied when `subset_cells = TRUE`.")

    if (!is.null(subset_by)) {
      # Subset by metadata column values
      if (!subset_by %in% colnames(seurat_object@meta.data))
        stop("`subset_by` column '", subset_by, "' not found in @meta.data.")
      keep_cells <- rownames(seurat_object@meta.data)[
        as.character(seurat_object@meta.data[[subset_by]]) %in% as.character(subset_values)
      ]
      message("Subsetting by '", subset_by, "' (",
              paste(subset_values, collapse = ", "), ") - keeping ",
              length(keep_cells), " / ", ncol(seurat_object), " cells.")
    } else {
      # subset_values are treated as cell barcodes directly
      keep_cells <- intersect(as.character(subset_values), colnames(seurat_object))
      message("Subsetting by cell barcodes - keeping ",
              length(keep_cells), " / ", ncol(seurat_object), " cells.")
    }

    if (length(keep_cells) == 0)
      stop("No cells remain after subsetting. Check `subset_by` and `subset_values`.")

    seurat_object <- seurat_object[, keep_cells]
  }

  # ── 1. Validate metadata columns ────────────────────────────────────────────
  if (!group.by %in% colnames(seurat_object@meta.data))
    stop("`group.by` column '", group.by, "' not found in @meta.data.")
  if (!is.null(split.by) && !split.by %in% colnames(seurat_object@meta.data))
    stop("`split.by` column '", split.by, "' not found in @meta.data.")

  sig_col <- if (!is.null(sig_group_by)) sig_group_by else group.by
  if (!sig_col %in% colnames(seurat_object@meta.data))
    stop("`sig_group_by` column '", sig_col, "' not found in @meta.data.")

  # ── 2. Build gene sets ───────────────────────────────────────────────────────
  message("Preparing gene sets...")
  gs_res <- .build_gene_sets_sc(
    gene_sets            = gene_sets,
    deg_df               = deg_df,
    deg_gene_column         = deg_gene_column,
    deg_group_column        = deg_group_column,
    deg_fc_column           = deg_fc_column,
    deg_padj_column         = deg_padj_column,
    deg_padj_cutoff      = deg_padj_cutoff,
    deg_top_n            = deg_top_n,
    gene_set_library     = gene_set_library,
    gene_set_subcategory = gene_set_subcategory,
    search_terms         = search_terms,
    species              = species
  )
  gene_sets_use <- gs_res$gene_sets
  feature_cats  <- gs_res$feature_cats

  if (length(gene_sets_use) == 0)
    stop("No gene sets available after preparation.")

  # ── 2b. Preview gene sets + interactive confirmation ─────────────────────────
  {
    n_gs    <- length(gene_sets_use)
    preview <- utils::head(names(gene_sets_use), 10L)

    message("\n  Found ", n_gs, " gene set(s).")
    message("  First ", length(preview), ":")
    for (i in seq_along(preview))
      message("    ", formatC(i, width = 2), ". ", preview[i])
    if (n_gs > 10L)
      message("  ... and ", n_gs - 10L, " more.")
    message("")

    # In an interactive session ask before committing to (potentially) a long
    # computation.  In non-interactive / Rmd knit / batch mode, proceed silently.
    if (interactive()) {
      ans <- trimws(readline(
        "  Proceed? [y = all  |  <number> = keep first N  |  n = abort]: "
      ))

      if (grepl("^[Nn]", ans)) {
        stop(
          "Aborted by user.\n",
          "Refine gene_set_library, search_terms, or gene_sets and rerun.",
          call. = FALSE
        )
      } else if (grepl("^[0-9]+$", ans)) {
        n_keep <- as.integer(ans)
        if (n_keep < 1L)
          stop("Must keep at least 1 gene set.", call. = FALSE)
        if (n_keep < n_gs) {
          message("  Keeping first ", n_keep, " of ", n_gs, " gene set(s).")
          gene_sets_use <- gene_sets_use[seq_len(n_keep)]
          if (!is.null(feature_cats))
            feature_cats <- feature_cats[
              feature_cats$pathway %in% names(gene_sets_use), , drop = FALSE]
        } else {
          message("  ", n_keep, " >= total (", n_gs, ") - keeping all.")
        }
      } else if (!grepl("^[Yy]", ans)) {
        stop("Unrecognized input - enter y, n, or a number.", call. = FALSE)
      }
      # y → fall through, use all
    }
  }

  # ── 2c. Build descriptive label for output file names ────────────────────────
  # Ensures rerunning with different gene sets / terms never silently overwrites
  # previous results.
  {
    .clean <- function(x) gsub("[^A-Za-z0-9]", "_",
                                gsub("_+", "_",
                                     gsub("^_|_$", "", x)))
    gs_label <- if (!is.null(search_terms)) {
      # Mode 4: reflect AND/OR logic in the filename
      .clean(.format_search_terms(search_terms, sep_and = "_AND_",
                                   sep_or = "_OR_", wrap_and = FALSE))
    } else if (!is.null(gene_set_library)) {
      # Mode 3: library code + optional subcategory
      base <- .clean(gene_set_library)
      if (!is.null(gene_set_subcategory))
        paste0(base, "_", .clean(gene_set_subcategory))
      else
        base
    } else if (!is.null(deg_df)) {
      # Mode 2: group column name(s) from the DEG data frame
      .clean(paste(deg_group_column, collapse = "_"))
    } else {
      # Mode 1: user-supplied named list - label as "custom"
      "custom"
    }
    # Create a dedicated subfolder - each gene-set run gets its own directory
    # so reruns with different terms never overwrite previous results.
    sc_ssgsea_dir <- file.path(output_dir, paste0("SC ssGSEA ", gs_label))
    dir.create(sc_ssgsea_dir, recursive = TRUE, showWarnings = FALSE)
    message("  Output folder: ", sc_ssgsea_dir)

    # File prefix carries the gene-set label
    file_prefix <- trimws(paste(file_prefix, gs_label))
  }

  # ── 3. Downsample ────────────────────────────────────────────────────────────
  if (!is.null(downsample)) {
    message("Downsampling to ", downsample,
            " cells per '", group.by, "' group...")
    Seurat::Idents(seurat_object) <- seurat_object@meta.data[[group.by]]
    obj_use <- subset(seurat_object, downsample = downsample)
    message("  ", ncol(obj_use), " cells after downsampling (from ",
            ncol(seurat_object), ").")
  } else {
    obj_use <- seurat_object
  }

  # ── 4. Extract normalized expression ────────────────────────────────────────
  # .get_layer_data() handles Seurat v3 (slot=) / v5 (layer=) / BPCells lazy
  # matrices (subsets rows first, then materialises to dgCMatrix).
  message("Extracting RNA normalized counts (data layer)...")
  expr_mat <- .get_layer_data(obj_use, assay = "RNA", layer = "data")

  # ── 5. Score cells ────────────────────────────────────────────────────────────
  #
  # Two backends:
  #   "gsva"  - full ssGSEA via GSVA.  Ranks ALL genes per cell once per call,
  #             then scores all gene sets against that ranking.  Slow on large
  #             objects but gives true ssGSEA scores.  BiocParallel progress bar
  #             shows per-gene-set progress.
  #   "ucell" - U statistic (AUC of top-maxRank recovery curve).  Only partially
  #             sorts each cell (top ucell_max_rank genes), chunks by CELLS
  #             (not gene sets), and scores all gene sets per cell-chunk.
  #             ~15-30x faster on typical scRNA-seq data; scores in [0,1].
  #
  # Both save a combined-scores RDS checkpoint after completion.
  # resume = TRUE loads the checkpoint and skips scoring entirely.
  {
    method    <- match.arg(method, c("gsva", "ucell"))
    cache_dir <- file.path(sc_ssgsea_dir, ".cache")
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    ckpt_all  <- file.path(cache_dir, "scores_combined.rds")

    n_gs_score   <- length(gene_sets_use)
    n_cells_score <- ncol(expr_mat)

    # ── Smart suggestion: warn when GSVA will be painfully slow ─────────────
    if (method == "gsva" && n_cells_score > 20000 && n_gs_score > 30) {
      message(
        "\n  Heads-up: ", n_cells_score, " cells × ", n_gs_score,
        " gene sets with GSVA (ssGSEA) can take 30-90+ minutes.\n",
        "  Tip: try method = 'ucell' for ~15-30x faster scoring with\n",
        "  similar biological results (U statistic, scores in [0,1]).\n",
        "  Install: BiocManager::install('UCell')\n"
      )
    }

    if (isTRUE(resume) && file.exists(ckpt_all)) {
      message("Resuming: loading scores checkpoint (", ckpt_all, ")...")
      scores <- readRDS(ckpt_all)
      message("  Loaded: ", nrow(scores), " pathways × ",
              ncol(scores), " cells.")
    } else {
      t0 <- proc.time()[["elapsed"]]

      if (method == "ucell") {
        message("Scoring with UCell: ", n_cells_score, " cells × ",
                n_gs_score, " gene sets",
                " (top ", ucell_max_rank, " genes, ",
                chunk_size, " cells/chunk",
                if (cores > 1) paste0(", ", cores, " cores") else "",
                ")...")
        scores <- tryCatch(
          .run_ucell(expr_mat, gene_sets_use,
                     max_rank   = ucell_max_rank,
                     chunk_size = chunk_size,
                     cores      = cores),
          error = function(e) stop("UCell failed: ", conditionMessage(e))
        )
      } else {
        message("Scoring with GSVA (ssGSEA): ", n_cells_score, " cells × ",
                n_gs_score, " gene sets",
                if (cores > 1) paste0(" (", cores, " cores)") else "",
                "...")
        message("  Progress bar below tracks gene sets (if BiocParallel available).")
        scores <- tryCatch(
          .run_gsva_ssgsea(expr_mat, gene_sets_use,
                           normalize = ssgsea.norm,
                           min_size  = min.size,
                           cores     = cores),
          error = function(e) stop("GSVA ssGSEA failed: ", conditionMessage(e))
        )
      }

      t_elapsed <- proc.time()[["elapsed"]] - t0
      message(sprintf("  Scoring complete in %.1f s: %d pathways × %d cells.",
                      t_elapsed, nrow(scores), ncol(scores)))

      saveRDS(scores, ckpt_all)
      message("  Checkpoint saved (resume = TRUE will skip scoring next run).")
    }
  }
  # scores is pathways × cells

  # ── 6. Save scores CSV ───────────────────────────────────────────────────────
  cells_scored <- colnames(scores)
  meta_scored  <- obj_use@meta.data[cells_scored, , drop = FALSE]

  scores_df <- as.data.frame(t(scores))           # cells × pathways
  scores_df$.cell_barcode  <- rownames(scores_df)
  scores_df[[group.by]]    <- meta_scored[[group.by]]
  if (!is.null(split.by) && split.by %in% colnames(meta_scored))
    scores_df[[split.by]]  <- meta_scored[[split.by]]

  csv_scores <- file.path(sc_ssgsea_dir,
                           paste0(file_prefix, " ssGSEA_scores.csv"))
  utils::write.csv(scores_df, csv_scores, row.names = FALSE)
  message("  Saved: ", csv_scores)

  # ── 7. Significance testing ──────────────────────────────────────────────────
  message("Running significance tests (", fit, ") per pathway...")

  sig_vec   <- as.character(meta_scored[[sig_col]])
  grp_tally <- table(sig_vec)
  valid_grps <- names(grp_tally[grp_tally >= 2])

  sig_df <- if (length(valid_grps) < 2) {
    warning("Fewer than 2 groups with ≥2 cells in '", sig_col,
            "' - significance testing skipped.")
    data.frame(pathway   = rownames(scores),
               statistic = NA_real_,
               p_value   = NA_real_,
               p_adj     = NA_real_,
               stringsAsFactors = FALSE)
  } else {
    idx <- which(sig_vec %in% valid_grps)
    .sc_significance(scores[, idx, drop = FALSE], sig_vec[idx], fit = fit)
  }

  csv_sig <- file.path(sc_ssgsea_dir,
                        paste0(file_prefix, " ssGSEA_significance.csv"))
  utils::write.csv(sig_df, csv_sig, row.names = FALSE)
  message("  Saved: ", csv_sig)

  # Within-split significance
  if (!is.null(split.by) && split.by %in% colnames(meta_scored)) {
    spl_vec      <- as.character(meta_scored[[split.by]])
    sig_by_split <- do.call(rbind, lapply(unique(spl_vec), function(lv) {
      idx  <- which(spl_vec == lv)
      sv   <- sig_vec[idx]
      vg   <- names(table(sv)[table(sv) >= 2])
      if (length(vg) < 2) return(NULL)
      res  <- .sc_significance(
        scores[, idx[sv %in% vg], drop = FALSE],
        sv[sv %in% vg],
        fit = fit
      )
      res[[split.by]] <- lv
      res
    }))
    if (!is.null(sig_by_split) && nrow(sig_by_split) > 0) {
      csv_sig2 <- file.path(sc_ssgsea_dir,
                             paste0(file_prefix,
                                    " ssGSEA_significance_by_", split.by,
                                    ".csv"))
      utils::write.csv(sig_by_split, csv_sig2, row.names = FALSE)
      message("  Saved: ", csv_sig2)
    }
  }

  # ── 8. Select pathways for display ──────────────────────────────────────────
  if (show_only_significant && !all(is.na(sig_df$p_adj))) {
    show_pws <- sig_df$pathway[!is.na(sig_df$p_adj) &
                                  sig_df$p_adj < p_cutoff]
    if (length(show_pws) == 0) {
      n_fb <- min(30L, nrow(sig_df))
      message("  No pathways significant at p_adj < ", p_cutoff,
              "; falling back to top ", n_fb, " by test statistic.")
      show_pws <- utils::head(
        sig_df$pathway[order(sig_df$p_adj, na.last = TRUE)], n_fb)
    } else {
      message("  ", length(show_pws), " significant pathway(s) at p_adj < ",
              p_cutoff, ".")
    }
  } else {
    show_pws <- rownames(scores)
  }

  # ── 9. Resolve group / split colors ────────────────────────────────────────
  # Factor levels preserved; user-supplied colors take priority - only
  # auto-generate what the user did not provide.
  .lvls <- function(col) if (is.factor(col)) levels(col) else unique(as.character(col))
  grp_lvls <- .lvls(seurat_object@meta.data[[group.by]])
  spl_lvls <- if (!is.null(split.by) && split.by %in% colnames(seurat_object@meta.data))
    .lvls(seurat_object@meta.data[[split.by]]) else character(0)

  n_grp <- length(grp_lvls)
  n_spl <- length(spl_lvls)

  # Color resolution priority:
  #   1. Explicit argument passed by user
  #   2. PrepObject stored colors (@misc$nk_settings$colors)
  #   3. Auto Nour_pal: group.by → "all" (≤8) or "spectrum" (>8)
  #                     split.by → "spectrum" (always distinct)
  if (is.null(group_colors))
    group_colors <- .nk_colors(seurat_object, group.by) %||%
      stats::setNames(Nour_pal(if (n_grp <= 8) "all" else "spectrum")(max(n_grp, 1L)), grp_lvls)
  if (is.null(split_colors) && n_spl > 0)
    split_colors <- .nk_colors(seurat_object, split.by) %||%
      stats::setNames(Nour_pal("spectrum")(n_spl), spl_lvls)

  # ── 10. Heatmap ──────────────────────────────────────────────────────────────
  message("Generating ssGSEA heatmap...")
  pdf_ht <- file.path(sc_ssgsea_dir,
                       paste0(file_prefix, " ssGSEA heatmap.pdf"))
  .sc_ssgsea_heatmap(
    scores_mat     = scores,
    meta_df        = meta_scored,
    group_by       = group.by,
    split_by       = split.by,
    split_levels   = spl_lvls,
    group_levels   = grp_lvls,
    feature_cats   = feature_cats,
    show_pws       = show_pws,
    group_colors   = group_colors,
    split_colors   = split_colors,
    heatmap_params = heatmap_params,
    heatmap_colors = heatmap_colors,
    pdf_path       = pdf_ht
  )
  .write_legend_sidecar(pdf_ht, paste0(
    "Heatmap of mean row-z-scored ",
    if (method == "ucell") "UCell" else "ssGSEA",
    " enrichment scores per ", group.by,
    if (!is.null(split.by)) paste0(", split by ", split.by) else "",
    ". Scores computed with ",
    if (method == "ucell")
      paste0("UCell (rank-based U statistic; max rank: ", ucell_max_rank, ")")
    else
      paste0("GSVA ssGSEA", if (ssgsea.norm) " (normalized)" else ""),
    ". ", length(show_pws), " pathway(s) shown",
    if (show_only_significant)
      paste0(" (BH-adjusted p < ", p_cutoff, " by one-way ", fit, ")")
    else "",
    ". Rows: gene sets; columns: group means (one column per ", group.by, " level). ",
    "Color scale: row z-score, symmetric around 0 (blue = low, white = 0, red = high)."
  ))

  # ── 10b. Boxplots ─────────────────────────────────────────────────────────────
  if (add_boxplots && length(show_pws) > 0) {
    message("Generating boxplots for ", length(show_pws), " pathway(s)...")
    pdf_bx <- file.path(sc_ssgsea_dir,
                         paste0(file_prefix, " ssGSEA boxplots.pdf"))
    .sc_ssgsea_boxplots(
      scores_df    = scores_df,
      group_by     = group.by,
      split_by     = split.by,
      group_levels = grp_lvls,
      split_levels = spl_lvls,
      show_pws     = show_pws,
      group_colors = group_colors,
      split_colors = split_colors,
      add_pvalues  = add_pvalues,
      pdf_path     = pdf_bx
    )
    .write_legend_sidecar(pdf_bx, paste0(
      "Boxplots of per-cell ",
      if (method == "ucell") "UCell" else "ssGSEA",
      " enrichment scores for ",
      length(show_pws), " pathway(s)",
      if (show_only_significant)
        paste0(" (selected at BH-adjusted p < ", p_cutoff, " by one-way ", fit, ")")
      else "",
      ". Each page = one pathway. x-axis: ",
      if (!is.null(split.by))
        paste0(split.by, "; facets: ", group.by)
      else
        group.by,
      "; y-axis: ",
      if (method == "ucell") "UCell score (per cell)." else "ssGSEA enrichment score (per cell).",
      if (add_pvalues) " BH-adjusted p-value brackets shown (ggpubr)." else ""
    ))
  }

  # ── 10c. Write methods JSON ───────────────────────────────────────────────────
  {
    gs_desc <- if (!is.null(search_terms)) {
      paste0("MSigDB gene sets matching: ", .format_search_terms(search_terms))
    } else if (!is.null(gene_set_library)) {
      paste0("MSigDB ", gene_set_library,
             if (!is.null(gene_set_subcategory))
               paste0("/", gene_set_subcategory) else "")
    } else if (!is.null(deg_df)) {
      paste0("custom gene sets from DEG data frame (group column(s): ",
             paste(deg_group_column, collapse = " / "), ")")
    } else {
      paste0("user-supplied named list (", length(gene_sets_use), " sets)")
    }

    n_cells_scored <- ncol(obj_use)
    n_cells_total  <- ncol(seurat_object)

    methods_text <- paste0(
      "Single-cell gene set scoring was performed on ",
      n_cells_scored,
      if (n_cells_scored < n_cells_total)
        paste0(" cells (downsampled from ", n_cells_total, ")")
      else " cells",
      " from ", trimws(paste(object_name, subset_name)),
      " using ", length(gene_sets_use), " gene set(s) (", gs_desc, "). ",
      if (method == "ucell")
        paste0(
          "Per-cell scores were computed using the UCell package ",
          "(rank-based U statistic; maximum rank: ", ucell_max_rank,
          "; minimum gene-set size: ", min.size, "). "
        )
      else
        paste0(
          "Per-cell enrichment scores were computed using the GSVA package ",
          "(ssGSEA method",
          if (ssgsea.norm) ", normalized to [0,1]" else "",
          "; minimum gene-set size: ", min.size, "). "
        ),
      "Differential enrichment across ", group.by, " groups",
      if (!is.null(split.by))
        paste0(" (and within each ", split.by, " level)")
      else "",
      " was assessed using one-way ", fit,
      " with Benjamini-Hochberg correction. ",
      length(show_pws),
      if (show_only_significant)
        paste0(" pathway(s) with adjusted p < ", p_cutoff, " are highlighted.")
      else
        " pathway(s) shown (all)."
    )

    tryCatch(
      .write_subdir_params(
        output_dir   = sc_ssgsea_dir,
        extra_params = list(
          date                  = format(Sys.Date()),
          sc_ssgsea_group_by    = group.by,
          sc_ssgsea_split_by    = if (is.null(split.by)) "none" else split.by,
          sc_ssgsea_method      = method,
          sc_ssgsea_assay       = "RNA",
          sc_ssgsea_gene_sets   = gs_label,
          sc_ssgsea_n_gene_sets = length(gene_sets_use),
          sc_ssgsea_min_size    = min.size,
          sc_ssgsea_norm        = if (method == "gsva") ssgsea.norm else NA,
          sc_ssgsea_ucell_max_rank = if (method == "ucell") ucell_max_rank else NA,
          sc_ssgsea_fit         = fit,
          sc_ssgsea_p_cutoff    = p_cutoff,
          sc_ssgsea_n_cells     = n_cells_scored,
          sc_ssgsea_downsample  = if (is.null(downsample)) "none"
                                  else as.character(downsample),
          sc_ssgsea_subset_cells  = isTRUE(subset_cells),
          sc_ssgsea_subset_by     = if (isTRUE(subset_cells)) subset_by else NA,
          sc_ssgsea_subset_values = if (isTRUE(subset_cells)) subset_values else NA,
          methods_text          = methods_text
        )
      ),
      error = function(e)
        warning("Could not write analysis_params.json: ", conditionMessage(e))
    )
  }

  # ── 11. Return ───────────────────────────────────────────────────────────────
  if (add_to_object) {
    message("Adding scores to Seurat object @meta.data...")
    scores_t   <- as.data.frame(t(scores))   # cells × pathways
    n_all      <- ncol(seurat_object)
    n_scored   <- nrow(scores_t)
    all_cells  <- colnames(seurat_object)

    for (pw in colnames(scores_t)) {
      vals           <- rep(NA_real_, n_all)
      names(vals)    <- all_cells
      matched        <- intersect(rownames(scores_t), all_cells)
      vals[matched]  <- scores_t[matched, pw]
      seurat_object@meta.data[[pw]] <- vals
    }

    if (n_scored < n_all)
      message("  Note: ", n_scored, "/", n_all, " cells scored; ",
              n_all - n_scored,
              " cells have NA (not in downsampled set).")

    message("RunSCssGSEA complete.")
    return(invisible(seurat_object))
  }

  message("RunSCssGSEA complete.")
  invisible(scores_df)
}
