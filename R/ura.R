# =============================================================================
# RunURA — Upstream Regulator Analysis via decoupleR + CollecTRI
# =============================================================================
#
# Public function : RunURA()
# Private helpers : .ura_scores_to_matrix(), .ura_contrasts(),
#                   .ura_differential(), .ura_plot_heatmap(),
#                   .ura_plot_bubble()
# =============================================================================


# ── Private helpers ───────────────────────────────────────────────────────────

# Normalize any CollecTRI-style data frame to exactly 3 columns: source, target, mor.
# Handles all OmnipathR column naming variants across versions.
# Returns NULL only if source or target columns cannot be identified.
# min_tfs: if > 0, also returns NULL when unique sources < min_tfs (for cache validation).
.ura_normalize_net <- function(net, min_tfs = 0L) {
  if (!is.data.frame(net) || nrow(net) == 0L) return(NULL)
  cn <- colnames(net)

  # source
  if (!"source" %in% cn) {
    alt <- intersect(c("source_genesymbol", "genesymbol_of_source",
                       "tf", "TF", "regulator"), cn)
    if (!length(alt)) return(NULL)
    net[["source"]] <- net[[alt[1L]]]
  }

  # target
  if (!"target" %in% cn) {
    alt <- intersect(c("target_genesymbol", "genesymbol_of_target",
                       "gene", "target_gene"), cn)
    if (!length(alt)) return(NULL)
    net[["target"]] <- net[[alt[1L]]]
  }

  # mor (mode of regulation)
  # OmniPath TSV downloads use Python-style "True"/"False" strings, not R
  # logical TRUE/FALSE or integer 1/0 — handle all variants.
  .parse_omni_bool <- function(x) {
    as.integer(trimws(as.character(x)) %in% c("1", "TRUE", "True", "true", "T"))
  }

  if (!"mor" %in% cn) {
    if ("weight" %in% cn) {
      net[["mor"]] <- as.numeric(net[["weight"]])
    } else if (all(c("is_stimulation", "is_inhibition") %in% cn)) {
      # stim=1,inhib=0 → mor= 1 (activation)
      # stim=0,inhib=1 → mor=-1 (inhibition)
      # stim=0,inhib=0 → mor= 0 (undirected — filtered below)
      net[["mor"]] <- .parse_omni_bool(net[["is_stimulation"]]) -
                      .parse_omni_bool(net[["is_inhibition"]])
    } else {
      net[["mor"]] <- 1L  # assume activation if no sign info
    }
  }

  out <- data.frame(
    source = as.character(net[["source"]]),
    target = as.character(net[["target"]]),
    mor    = as.numeric(net[["mor"]]),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$source) & nzchar(out$target) & is.finite(out$mor) & out$mor != 0, ]

  # Deduplicate (source, target) pairs — CollecTRI often has repeated edges
  # from multiple evidence sources.  Average mor; drop contradictory pairs (→ 0).
  edge_key <- paste(out$source, out$target, sep = "\t")
  if (anyDuplicated(edge_key)) {
    out <- aggregate(mor ~ source + target, data = out, FUN = mean)
    out <- out[out$mor != 0, ]
  }

  if (min_tfs > 0L && length(unique(out$source)) < min_tfs) return(NULL)
  out
}


# Robust CollecTRI fetch — attempts in order:
#  0. ~/.scSidekick/collectri_{organism}.rds  (user cache, populated on first success)
#  1. inst/extdata bundle (if shipped with the package)
#  2. Standard get_collectri()
#  3. OmnipathR evidence options forced FALSE (fixes logical(0) static-fallback bug)
#  4. OmnipathR disk cache (raw file may survive even if post-processing crashed)
# Every result passes through .ura_normalize_net() — columns are always standardised
# and partial/malformed fetches (< 50 TFs) are rejected.
# On success, the normalized network is saved to the user cache.
.ura_get_collectri <- function(organism) {

  user_cache <- path.expand(
    file.path("~", ".scSidekick", paste0("collectri_", organism, ".rds"))
  )

  .try_net <- function(fn) {
    net <- tryCatch(suppressMessages(fn()), error = function(e) NULL)
    .ura_normalize_net(net, min_tfs = 50L)
  }

  .save_cache <- function(net) {
    tryCatch({
      dir.create(dirname(user_cache), recursive = TRUE, showWarnings = FALSE)
      saveRDS(net, user_cache)
      message("scSidekick RunURA: CollecTRI cached → ", user_cache)
    }, error = function(e) NULL)
  }

  # 0 — user cache (~/.scSidekick/)
  if (file.exists(user_cache)) {
    net <- tryCatch(.ura_normalize_net(readRDS(user_cache), min_tfs = 50L), error = function(e) NULL)
    if (!is.null(net)) {
      message("scSidekick RunURA: Using cached CollecTRI (",
              length(unique(net$source)), " TFs, ", user_cache, ")")
      return(net)
    }
    # Cache is bad — delete it so we don't keep loading it
    message("scSidekick RunURA: Cached network invalid/incomplete — deleting and re-fetching.")
    unlink(user_cache)
  }

  # 1 — inst/extdata bundle
  pkg_file <- system.file("extdata",
                           paste0("collectri_", organism, ".rds"),
                           package = "scSidekick")
  if (nzchar(pkg_file) && file.exists(pkg_file)) {
    net <- tryCatch(.ura_normalize_net(readRDS(pkg_file), min_tfs = 50L), error = function(e) NULL)
    if (!is.null(net)) {
      message("scSidekick RunURA: Using bundled CollecTRI.")
      .save_cache(net); return(net)
    }
  }

  # 2 — plain get_collectri()
  net <- .try_net(function()
    decoupleR::get_collectri(organism = organism, split_complexes = FALSE))
  if (!is.null(net)) { .save_cache(net); return(net) }

  # 3 — evidence options disabled (fixes logical(0) bug in OmnipathR ≤1.12)
  for (sc in c(FALSE, TRUE)) {
    net <- .try_net(function() {
      old <- options(OmnipathR.evidences        = FALSE,
                     OmnipathR.keep_evidences   = FALSE,
                     OmnipathR.strict_evidences = FALSE)
      on.exit(options(old), add = TRUE)
      decoupleR::get_collectri(organism = organism, split_complexes = sc)
    })
    if (!is.null(net)) { .save_cache(net); return(net) }
  }

  # 4 — OmnipathR disk cache (raw file may survive even if post-processing crashed)
  net <- tryCatch({
    cache_dir <- tools::R_user_dir("OmnipathR", "cache")
    cached    <- list.files(cache_dir, pattern = "collectri",
                            recursive = TRUE, full.names = TRUE)
    if (!length(cached)) stop("no cached files")
    f   <- cached[order(file.mtime(cached), decreasing = TRUE)[1L]]
    raw <- tryCatch(
      utils::read.table(f, header = TRUE, sep = "\t", quote = "",
                        stringsAsFactors = FALSE),
      error = function(e) as.data.frame(readRDS(f))
    )
    .ura_normalize_net(raw, min_tfs = 50L)
  }, error = function(e) NULL)
  if (!is.null(net)) { .save_cache(net); return(net) }

  # 5 — Direct OmniPath REST download (bypasses OmnipathR entirely)
  taxon <- if (organism == "human") "9606" else "10090"
  net <- tryCatch({
    url <- paste0(
      "https://omnipathdb.org/interactions",
      "?datasets=collectri&organisms=", taxon,
      "&genesymbols=1",
      "&fields=source_genesymbol,target_genesymbol,is_stimulation,is_inhibition"
    )
    tmp <- tempfile(fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    utils::download.file(url, tmp, quiet = TRUE, method = "auto",
                         extra = "--max-time 30")
    raw <- utils::read.table(tmp, header = TRUE, sep = "\t",
                             quote = "", stringsAsFactors = FALSE)
    .ura_normalize_net(raw, min_tfs = 50L)
  }, error = function(e) NULL)
  if (!is.null(net)) { .save_cache(net); return(net) }

  stop(
    "Could not fetch CollecTRI — OmniPath server unreachable and all fallbacks failed.\n\n",
    "PERMANENT FIX — bundle CollecTRI inside scSidekick (run once from any machine\n",
    "where OmniPath is accessible, e.g. outside your institution network):\n",
    "  source('data-raw/fetch_collectri.R')   # creates inst/extdata/ bundles\n",
    "  devtools::load_all()                    # pick up the bundle\n\n",
    "QUICK FIX — cache when reachable, reuse forever:\n",
    "  net <- decoupleR::get_collectri(organism = '", organism, "')\n",
    "  saveRDS(net, '", user_cache, "')\n\n",
    "DIRECT DOWNLOAD (if curl works in your terminal):\n",
    "  curl 'https://omnipathdb.org/interactions?datasets=collectri&organisms=",
    taxon, "&genesymbols=1' -o /tmp/collectri.tsv\n",
    "  # Then in R: net <- read.table('/tmp/collectri.tsv', header=TRUE, sep='\\t')\n",
    "  #            saveRDS(net, '", user_cache, "')"
  )
}


# Pivot decoupleR long-format scores → TF × pseudobulk-sample matrix
.ura_scores_to_matrix <- function(scores) {
  tfs      <- sort(unique(as.character(scores$source)))
  conds    <- unique(as.character(scores$condition))
  tf_idx   <- match(as.character(scores$source),    tfs)
  cond_idx <- match(as.character(scores$condition), conds)
  mat      <- matrix(NA_real_, nrow = length(tfs), ncol = length(conds),
                     dimnames = list(tfs, conds))
  mat[cbind(tf_idx, cond_idx)] <- as.numeric(scores$score)
  mat
}


# Normalize contrast argument to a list of c("test","ref") pairs.
# Accepts: c("A","B")  or  list(c("A","B"), c("C","B"))
.ura_contrasts <- function(contrast) {
  if (is.null(contrast)) return(NULL)
  if (is.character(contrast)) {
    if (length(contrast) != 2L)
      stop("A single contrast must be c('test', 'reference').")
    return(list(contrast))
  }
  if (is.list(contrast)) {
    bad <- !vapply(contrast, function(x) is.character(x) && length(x) == 2L, logical(1L))
    if (any(bad))
      stop("Each element of contrast must be c('test', 'reference').")
    return(contrast)
  }
  stop("'contrast' must be c('A','B') or list(c('A','B'), c('C','B')).")
}


# Differential TF activity for ONE contrast pair, stratified by group.by ×
# split.by combinations.  Uses limma when ≥3 donors each side, t-test otherwise.
# Returns data.frame: TF, cell_type, split_level, mean_a, mean_b,
#                     activity_diff, t, p.value, p.adj, contrast
.ura_differential <- function(act_mat, col_meta,
                               group.by, split.by,
                               contrast.by, contrast,
                               has_replicates) {

  # Build stratification keys: one entry per (group × split) combination
  if (!is.null(group.by) && !is.null(split.by)) {
    strat_keys <- paste(col_meta[[group.by]], col_meta[[split.by]], sep = "\t")
  } else if (!is.null(group.by)) {
    strat_keys <- col_meta[[group.by]]
  } else {
    strat_keys <- rep("All", nrow(col_meta))
  }

  strats <- split(col_meta, strat_keys)

  de_list <- lapply(names(strats), function(key) {
    sm <- strats[[key]]

    ct_cols <- intersect(sm$pb_col, colnames(act_mat))
    if (length(ct_cols) < 2L) return(NULL)

    ct_act  <- act_mat[, ct_cols, drop = FALSE]
    sm      <- sm[sm$pb_col %in% ct_cols, , drop = FALSE]

    grp_a <- intersect(sm$pb_col[sm[[contrast.by]] == contrast[1]], colnames(ct_act))
    grp_b <- intersect(sm$pb_col[sm[[contrast.by]] == contrast[2]], colnames(ct_act))

    if (length(grp_a) == 0L || length(grp_b) == 0L) return(NULL)

    mean_a <- rowMeans(ct_act[, grp_a, drop = FALSE], na.rm = TRUE)
    mean_b <- rowMeans(ct_act[, grp_b, drop = FALSE], na.rm = TRUE)
    diff   <- mean_a - mean_b

    # Parse group / split labels from strat key
    key_parts  <- strsplit(key, "\t", fixed = TRUE)[[1L]]
    cell_label <- key_parts[1L]
    split_label <- if (length(key_parts) > 1L) key_parts[2L] else NA_character_

    if (has_replicates &&
        length(grp_a) >= 3L && length(grp_b) >= 3L &&
        requireNamespace("limma", quietly = TRUE)) {

      sub    <- ct_act[, c(grp_a, grp_b), drop = FALSE]
      gvec   <- factor(c(rep(contrast[1], length(grp_a)),
                         rep(contrast[2], length(grp_b))),
                       levels = contrast)
      design <- model.matrix(~ gvec)
      fit    <- limma::lmFit(sub, design)
      fit    <- limma::eBayes(fit)
      top    <- limma::topTable(fit, coef = 2L, n = Inf, sort.by = "none")

      data.frame(
        TF            = rownames(top),
        cell_type     = cell_label,
        split_level   = split_label,
        mean_a        = mean_a[rownames(top)],
        mean_b        = mean_b[rownames(top)],
        activity_diff = top$logFC,
        t             = top$t,
        p.value       = top$P.Value,
        p.adj         = stats::p.adjust(top$P.Value, method = "BH"),
        stringsAsFactors = FALSE
      )

    } else {
      p_vals <- vapply(rownames(ct_act), function(tf) {
        a <- ct_act[tf, grp_a, drop = TRUE]
        b <- ct_act[tf, grp_b, drop = TRUE]
        a <- a[is.finite(a)]; b <- b[is.finite(b)]
        if (length(a) < 2L || length(b) < 2L) return(NA_real_)
        tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_)
      }, numeric(1L))

      data.frame(
        TF            = rownames(ct_act),
        cell_type     = cell_label,
        split_level   = split_label,
        mean_a        = mean_a,
        mean_b        = mean_b,
        activity_diff = diff,
        t             = NA_real_,
        p.value       = p_vals,
        p.adj         = stats::p.adjust(p_vals, method = "BH"),
        stringsAsFactors = FALSE
      )
    }
  })

  out <- do.call(rbind, Filter(Negate(is.null), de_list))
  if (!is.null(out)) {
    out$contrast  <- paste0(contrast[1], "_vs_", contrast[2])
    rownames(out) <- NULL
  }
  out
}


# Heatmap for URA results — one call per contrast when DE is available.
#
# DE mode  : Rows = top TFs, Columns = CellType (× split.by if set)
#            Values = activity_diff (test − ref).  Black dot = p.adj < 0.05.
#            Column groups titled by split.by level (e.g. "Female" / "Male").
# No-DE    : z-scored mean activity across groups.
.ura_plot_heatmap <- function(mean_act, n_top, de = NULL,
                               split.by, group.by,
                               contrast_name = NULL,
                               output_dir, object_name, subset_name) {

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) return(invisible(NULL))
  if (!requireNamespace("circlize",       quietly = TRUE)) return(invisible(NULL))

  use_de <- !is.null(de) && nrow(de) > 0L

  # ── TF selection ─────────────────────────────────────────────────────────
  if (use_de) {
    by_tf   <- tapply(abs(de$activity_diff), de$TF, max, na.rm = TRUE)
    top_tfs <- names(sort(by_tf, decreasing = TRUE))[seq_len(min(n_top, length(by_tf)))]
  } else {
    tf_var  <- apply(mean_act, 1L, stats::var, na.rm = TRUE)
    top_tfs <- names(sort(tf_var, decreasing = TRUE))[seq_len(min(n_top, nrow(mean_act)))]
  }

  if (use_de) {
    # ── DE mode: pivot activity_diff → TF × (cell_type × split) ───────────
    de_top    <- de[de$TF %in% top_tfs, , drop = FALSE]
    has_split <- !all(is.na(de_top$split_level))

    col_key <- if (has_split)
      paste(de_top$cell_type, de_top$split_level, sep = "\t")
    else
      de_top$cell_type

    cols <- unique(col_key)
    mat  <- matrix(NA_real_, nrow = length(top_tfs), ncol = length(cols),
                   dimnames = list(top_tfs, cols))
    pmat <- matrix(NA_real_, nrow = length(top_tfs), ncol = length(cols),
                   dimnames = list(top_tfs, cols))
    for (i in seq_len(nrow(de_top))) {
      tf <- de_top$TF[i]; cl <- col_key[i]
      if (tf %in% top_tfs && cl %in% cols) {
        mat[tf,  cl] <- de_top$activity_diff[i]
        pmat[tf, cl] <- de_top$p.adj[i]
      }
    }

    # Column split: split.by levels become group titles automatically
    col_split  <- if (has_split) sub(".*\t", "", cols) else NULL
    col_labels <- sub("\t.*", "", cols)

    # Contrast label for legend title  e.g. "Dementia.AD − NoDementia.Control"
    ctr_str <- if (!is.null(contrast_name))
      gsub("_vs_", " −\n", contrast_name) else "test − ref"

    lim     <- max(max(abs(mat), na.rm = TRUE), 0.5)
    col_fun <- circlize::colorRamp2(c(-lim, 0, lim),
                                    c("#003F5C", "white", "#F37388"))

    ht <- ComplexHeatmap::Heatmap(
      mat,
      name              = paste0("Activity\ndiff\n(", ctr_str, ")"),
      col               = col_fun,
      na_col            = "gray90",
      column_labels     = col_labels,
      column_split      = col_split,
      # No column_title here — let ComplexHeatmap use the split values
      # (e.g. "Female" / "Male") as each group's title automatically
      cluster_rows      = TRUE,
      cluster_columns   = FALSE,
      cluster_row_slices = FALSE,
      show_row_dend     = FALSE,
      show_row_names    = TRUE,
      show_column_names = TRUE,
      row_names_gp      = grid::gpar(fontsize = 8),
      column_names_gp   = grid::gpar(fontsize = 9),
      column_title_gp   = grid::gpar(fontsize = 10, fontface = "bold"),
      column_names_rot  = 45,
      rect_gp           = grid::gpar(col = "white", lwd = 0.5),
      cell_fun = function(j, i, x, y, width, height, fill) {
        p <- pmat[i, j]
        if (!is.na(p) && p < 0.05)
          grid::grid.points(x, y, pch = 16,
                            size  = grid::unit(1.5, "mm"),
                            gp    = grid::gpar(col = "black"))
      }
    )

    fig_title <- paste0("TF Activity — ", gsub("_vs_", " vs ", contrast_name %||% ""))

  } else {
    # ── No-contrast mode: z-scored mean activity ────────────────────────────
    mat_plot <- mean_act[intersect(top_tfs, rownames(mean_act)), , drop = FALSE]
    mat      <- t(scale(t(mat_plot)))
    mat[is.nan(mat)] <- 0

    col_split  <- if (!is.null(split.by)) sub(".*\t", "", colnames(mat)) else NULL
    colnames(mat) <- sub("\t.*", "", colnames(mat))

    col_fun <- circlize::colorRamp2(c(-2, 0, 2),
                                    c("#003F5C", "white", "#F37388"))

    ht <- ComplexHeatmap::Heatmap(
      mat,
      name              = "Activity\n(z-score)",
      col               = col_fun,
      column_split      = col_split,
      column_title_gp   = grid::gpar(fontsize = 10, fontface = "bold"),
      cluster_rows      = TRUE,
      cluster_columns   = FALSE,
      cluster_row_slices = FALSE,
      show_row_dend     = FALSE,
      show_row_names    = TRUE,
      show_column_names = TRUE,
      row_names_gp      = grid::gpar(fontsize = 8),
      column_names_gp   = grid::gpar(fontsize = 9),
      column_names_rot  = 45,
      rect_gp           = grid::gpar(col = "white", lwd = 0.5)
    )

    fig_title <- "TF Activity (CollecTRI)"
  }

  # ── File name: include group.by, split.by, contrast ──────────────────────
  tag_parts <- c(
    object_name,
    subset_name,
    "URA",
    if (!is.null(group.by))   paste0("by_", group.by),
    if (!is.null(split.by))   paste0("split_", split.by),
    if (!is.null(contrast_name)) gsub("_vs_", "-vs-", contrast_name),
    "heatmap"
  )
  fname    <- gsub("[^A-Za-z0-9._-]", "_", paste(Filter(nzchar, tag_parts), collapse = "_"))
  pdf_path <- file.path(output_dir, paste0(fname, ".pdf"))

  pdf_w <- max(4, ncol(mat) * 0.55 + 3.5)
  pdf_h <- max(4, nrow(mat) * 0.25 + 2.5)

  grDevices::pdf(pdf_path, width = pdf_w, height = pdf_h)
  ComplexHeatmap::draw(ht,
    # Overall figure title above everything (separate from per-group column titles)
    column_title        = fig_title,
    column_title_gp     = grid::gpar(fontsize = 11, fontface = "bold"),
    annotation_legend_side = "right",
    heatmap_legend_side    = "right"
  )
  grDevices::dev.off()
  message("scSidekick RunURA: Heatmap → ", pdf_path)

  # ── Legend sidecar ────────────────────────────────────────────────────────
  legend_text <- paste0(
    "Heatmap of transcription factor (TF) activity scores inferred by decoupleR ULM ",
    "using the CollecTRI network. ",
    if (use_de) paste0(
      "Color represents the difference in mean TF activity between ",
      gsub("_vs_", " and ", contrast_name %||% "groups"),
      " (positive = more active in test group, negative = less active). ",
      "Black dots indicate p.adj < 0.05 (Benjamini-Hochberg). "
    ) else
      "Color represents z-scored mean TF activity across groups. ",
    if (!is.null(split.by)) paste0("Columns are split by ", split.by, ". "),
    paste0("Top ", nrow(mat), " TFs selected by maximum |activity difference| across cell types.")
  )
  tryCatch(
    log_figure_legend(output_dir, basename(pdf_path), legend_text),
    error = function(e) NULL
  )

  invisible(pdf_path)
}


# Matrix bubble plot for differential TF activity.
# Layout: x = cell type, y = TF, color = activity_diff, size = -log10(p.adj).
# Faded points = non-significant (p.adj ≥ 0.05).
# Facets: contrast (columns) × split.by level (rows, when present).
# Multi-page PDF: max 40 TFs per page; overflow prints on subsequent pages.
.ura_plot_bubble <- function(de_tbl, n_top,
                              group.by = NULL, split.by = NULL,
                              output_dir, object_name, subset_name) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))

  de_use <- de_tbl[!is.na(de_tbl$p.value) & is.finite(de_tbl$activity_diff), ]
  if (nrow(de_use) == 0L) return(invisible(NULL))

  has_split      <- !all(is.na(de_use$split_level))
  multi_contrast <- length(unique(de_use$contrast)) > 1L

  # ── Global top N TFs by max |activity_diff| across all groups ─────────────
  by_tf   <- tapply(abs(de_use$activity_diff), de_use$TF, max, na.rm = TRUE)
  top_tfs <- names(sort(by_tf, decreasing = TRUE))[seq_len(min(n_top, length(by_tf)))]
  plot_df <- de_use[de_use$TF %in% top_tfs, , drop = FALSE]

  # TF order: most-inhibited at top, most-activated at bottom (matches heatmap)
  tf_mean    <- tapply(plot_df$activity_diff, plot_df$TF, mean, na.rm = TRUE)
  tf_ordered <- names(sort(tf_mean, decreasing = TRUE))  # top = most activated
  plot_df$TF <- factor(plot_df$TF, levels = rev(tf_ordered))  # rev → activated at top of y

  # Cell type order: sort alphabetically for a stable axis
  ct_levels      <- sort(unique(as.character(plot_df$cell_type)))
  plot_df$cell_type <- factor(plot_df$cell_type, levels = ct_levels)

  # Contrast label (strip underscores for display)
  plot_df$contrast_label <- gsub("_vs_", " vs ", plot_df$contrast)
  if (has_split) plot_df$split_label <- as.character(plot_df$split_level)

  plot_df$log10p <- pmin(-log10(plot_df$p.adj + 1e-300), 10)
  plot_df$signif <- !is.na(plot_df$p.adj) & plot_df$p.adj < 0.05

  lim <- max(abs(plot_df$activity_diff), na.rm = TRUE)
  lim <- max(lim, 0.2)

  # ── Build base plot (data swapped per page via %+%) ───────────────────────
  base_p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x     = cell_type,
                 y     = TF,
                 color = activity_diff,
                 size  = log10p,
                 alpha = signif)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient2(
      low      = "#003F5C",
      mid      = "white",
      high     = "#F37388",
      midpoint = 0,
      limits   = c(-lim, lim),
      name     = "Activity\ndifference"
    ) +
    ggplot2::scale_size_continuous(
      name  = expression(-log[10](p[adj])),
      range = c(0.5, 5)
    ) +
    ggplot2::scale_alpha_manual(
      values = c("FALSE" = 0.15, "TRUE" = 0.9),
      guide  = "none"
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = "Upstream Regulator Analysis") +
    theme_NourMin() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y      = ggplot2::element_text(size = 7),
      strip.text       = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.major = ggplot2::element_line(color = "gray92", linewidth = 0.3),
      panel.border     = ggplot2::element_rect(color = "gray80", fill = NA,
                                               linewidth = 0.4)
    )

  # ── Faceting ──────────────────────────────────────────────────────────────
  if (multi_contrast && has_split) {
    base_p <- base_p +
      ggplot2::facet_grid(split_label ~ contrast_label,
                          scales = "free_x", space = "free_x")
  } else if (multi_contrast) {
    base_p <- base_p +
      ggplot2::facet_wrap(~ contrast_label, nrow = 1, scales = "free_x")
  } else if (has_split) {
    base_p <- base_p +
      ggplot2::facet_wrap(~ split_label, nrow = 1, scales = "free_x")
  }

  # ── PDF dimensions ────────────────────────────────────────────────────────
  n_ct     <- length(ct_levels)
  n_ctr    <- length(unique(plot_df$contrast))
  n_split  <- if (has_split) length(unique(plot_df$split_level)) else 1L
  n_facet_cols <- if (multi_contrast) n_ctr else if (has_split) n_split else 1L
  n_facet_rows <- if (multi_contrast && has_split) n_split else 1L

  pdf_w <- max(6, n_ct * 0.5 * n_facet_cols + 4.5)

  # ── Multi-page: 40 TFs per page ───────────────────────────────────────────
  max_per_page <- 40L
  all_tfs  <- levels(plot_df$TF)          # already ordered top → bottom
  n_pages  <- ceiling(length(all_tfs) / max_per_page)
  tf_pages <- split(all_tfs,
                    ceiling(seq_along(all_tfs) / max_per_page))

  row_h    <- 0.22 * n_facet_rows         # inches per TF row
  header_h <- 2.0  + 1.2 * n_facet_rows  # title + strip + axis labels
  pdf_h    <- max(4, max_per_page * row_h + header_h)

  # ── File name ─────────────────────────────────────────────────────────────
  tag_parts <- c(
    object_name, subset_name, "URA",
    if (!is.null(group.by)) paste0("by_", group.by),
    if (!is.null(split.by)) paste0("split_", split.by),
    "bubble"
  )
  fname    <- gsub("[^A-Za-z0-9._-]", "_",
                   paste(Filter(nzchar, tag_parts), collapse = "_"))
  pdf_path <- file.path(output_dir, paste0(fname, ".pdf"))

  # ── Print one page per TF chunk ───────────────────────────────────────────
  grDevices::pdf(pdf_path, width = pdf_w, height = pdf_h, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page_tfs in tf_pages) {
    chunk_df      <- droplevels(plot_df[as.character(plot_df$TF) %in% page_tfs, ])
    chunk_df$TF   <- factor(as.character(chunk_df$TF), levels = page_tfs)
    p_page        <- base_p
    p_page$data   <- chunk_df
    print(p_page)
  }
  on.exit(NULL, add = FALSE)   # clear the guard — dev.off() ran normally
  grDevices::dev.off()

  message("scSidekick RunURA: Bubble plot → ", pdf_path,
          if (n_pages > 1L) paste0(" (", n_pages, " pages)") else "")

  # ── Legend sidecar ────────────────────────────────────────────────────────
  contrasts_str <- paste(unique(de_use$contrast), collapse = "; ")
  legend_text <- paste0(
    "Matrix bubble plot of differentially active transcription factors (TFs) ",
    "inferred by decoupleR ULM using the CollecTRI network. ",
    "X-axis: cell types. Y-axis: top ", length(all_tfs), " TFs ranked by maximum ",
    "activity difference across groups. ",
    "Point color: activity difference (red = more active in test group, ",
    "blue = less active). ",
    "Point size: -log10(adjusted p-value, BH). ",
    "Faded points: p.adj ≥ 0.05 (non-significant). ",
    if (has_split) paste0("Panels split by ", split.by, ". ") else "",
    "Contrasts: ", contrasts_str, ". ",
    if (n_pages > 1L) paste0("Printed across ", n_pages, " pages (", max_per_page,
                              " TFs per page).") else ""
  )
  tryCatch(
    log_figure_legend(output_dir, basename(pdf_path), legend_text),
    error = function(e) NULL
  )

  invisible(pdf_path)
}


# =============================================================================
# RunURA (exported)
# =============================================================================

#' Upstream Regulator Analysis via decoupleR + CollecTRI
#'
#' Infers transcription factor activity from pseudobulk gene expression using
#' \pkg{decoupleR} and the CollecTRI TF-target database.  Optionally tests for
#' differential TF activity between one or more condition pairs (limma when
#' \eqn{\geq 3} donors per side; t-test otherwise).
#'
#' @param seurat_object A Seurat object.
#' @param group.by Character or \code{NULL}.  Metadata column that defines the
#'   primary grouping (e.g. \code{"CellType"}, \code{"Cluster"}).
#'   \code{NULL} treats all cells as one group.
#' @param donor.by Character or \code{NULL}.  Replicates column
#'   (e.g. \code{"Donor.ID"}).  Enables limma-based differential testing when
#'   \eqn{\geq 3} donors per contrast arm.  \code{NULL} falls back to t-test.
#' @param split.by Character or \code{NULL}.  Additional stratification column
#'   (e.g. \code{"Sex"}).  Included in the pseudobulk grouping, generates
#'   \code{group × split} columns in the activity matrix, and facets the
#'   bubble plot.
#' @param contrast.by Character or \code{NULL}.  Condition column for DE
#'   (e.g. \code{"Diagnosis"}).
#' @param contrast A single \code{c("test", "ref")} pair or a list of such
#'   pairs, e.g. \code{list(c("AD","Control"), c("MCI","Control"))}.
#'   Required when \code{contrast.by} is set.
#' @param group.levels Character vector or \code{NULL}.  Restrict to specific
#'   levels of \code{group.by}.  \code{NULL} = all levels.
#' @param assay Character.  Seurat assay.  Default \code{"RNA"}.
#' @param layer Character.  Assay layer.  Default \code{"data"}.
#' @param use_pseudobulk Logical.  When \code{TRUE} (default) TF activity is
#'   inferred from pseudobulk expression averages — statistically rigorous for
#'   DE testing.  When \code{FALSE}, decoupleR runs directly on the full
#'   single-cell matrix (all cells × genes) and activity scores are then
#'   averaged per group; simpler, works without \code{donor.by}, but
#'   pseudoreplication means DE results should be interpreted with caution.
#' @param network A \code{data.frame} with columns \code{source}, \code{target},
#'   \code{mor}.  Fetched from CollecTRI when \code{NULL}.
#' @param method One of \code{"ulm"} (default), \code{"viper"}, \code{"mlm"}.
#' @param organism \code{"human"} (default) or \code{"mouse"}.
#' @param n_top Top N TFs to include in plots.  Default \code{25L}.
#' @param min_cells Minimum cells per pseudobulk group.  Default \code{10L}.
#' @param output_dir Directory for PDFs.  \code{NULL} = no files saved.
#' @param object_name Character.  Label prefix for output file names.
#' @param subset_name Character.  Optional subset label appended to output
#'   file names.
#' @param plot Logical.  Save plots.  Default \code{TRUE}.
#' @param caffeinate Logical.  Prevent the machine from sleeping during the run
#'   (macOS only; uses \code{caffeinate -i}).  Default \code{FALSE}.
#' @param force Logical.  When \code{TRUE}, recomputes pseudobulk even if a
#'   cached result is already stored in \code{seurat_object@misc$pseudobulk}.
#'   Default \code{FALSE}.
#' @param max_cells Integer or \code{NULL}.  When \code{use_pseudobulk = FALSE}
#'   and the object has more cells than this value, a stratified random sample
#'   of \code{max_cells} cells is drawn before running decoupleR.  Preserves
#'   group proportions.  \code{NULL} (default) uses all cells.
#'
#' @return Invisibly, a list:
#' \describe{
#'   \item{\code{activity}}{TF × group mean-activity matrix.}
#'   \item{\code{scores}}{Long-format decoupleR output (all pseudobulk samples).}
#'   \item{\code{de}}{Differential table, or \code{NULL}.  Columns:
#'     \code{TF}, \code{cell_type}, \code{split_level}, \code{mean_a},
#'     \code{mean_b}, \code{activity_diff}, \code{t}, \code{p.value},
#'     \code{p.adj}, \code{contrast}.}
#'   \item{\code{network}}{TF-target network used.}
#' }
#'
#' @export
RunURA <- function(
    seurat_object,
    group.by        = NULL,
    donor.by        = NULL,
    split.by        = NULL,
    contrast.by     = NULL,
    contrast        = NULL,
    group.levels    = NULL,
    assay           = "RNA",
    layer           = "data",
    use_pseudobulk  = TRUE,
    force           = FALSE,
    network         = NULL,
    method          = "ulm",
    organism        = "human",
    n_top           = 25L,
    min_cells       = 10L,
    output_dir      = NULL,
    object_name     = "",
    subset_name     = "",
    plot            = TRUE,
    caffeinate      = FALSE,
    max_cells       = NULL
) {

  # ── 0. Package check ──────────────────────────────────────────────────────
  if (!requireNamespace("decoupleR", quietly = TRUE))
    stop("decoupleR required.  Install: BiocManager::install('decoupleR')")

  # ── 0b. Caffeinate ────────────────────────────────────────────────────────
  caff_handle <- if (isTRUE(caffeinate)) .nk_caffeinate() else NULL
  on.exit(.nk_decaffeinate(caff_handle), add = TRUE)

  method   <- match.arg(method,   c("ulm", "viper", "mlm"))
  organism <- match.arg(organism, c("human", "mouse"))

  # Normalize contrast to a list of pairs (or NULL)
  contrast_list <- .ura_contrasts(contrast)
  if (!is.null(contrast.by) && is.null(contrast_list))
    stop("'contrast' is required when 'contrast.by' is set.")

  # ── 1. PrepObject walk-up ─────────────────────────────────────────────────
  output_dir  <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""
  subset_name <- if (nchar(subset_name) > 0) subset_name else
    .nk_setting(seurat_object, "subset_name") %||% ""

  # ── 2. Validate metadata columns ──────────────────────────────────────────
  meta     <- seurat_object@meta.data
  req_cols <- c(group.by, donor.by, split.by, contrast.by)
  bad_cols <- setdiff(req_cols, colnames(meta))
  if (length(bad_cols))
    stop("Column(s) not found in metadata: ", paste(bad_cols, collapse = ", "))

  # ── 3. Filter to requested group levels ───────────────────────────────────
  if (!is.null(group.levels) && !is.null(group.by)) {
    keep <- meta[[group.by]] %in% group.levels
    if (!any(keep))
      stop("No cells match group.levels in column '", group.by, "'.")
    seurat_object <- seurat_object[, keep]
    meta          <- seurat_object@meta.data
    message("scSidekick RunURA: Retained ", sum(keep), " cell(s) – ",
            paste(group.levels, collapse = ", "), ".")
  }

  # ── 4. Handle NULL group.by ───────────────────────────────────────────────
  using_dummy <- is.null(group.by)
  if (using_dummy) {
    seurat_object$.__ura_grp__. <- "All"
    group.by <- ".__ura_grp__."
    meta     <- seurat_object@meta.data
  }

  pb_group <- unique(c(group.by, donor.by, split.by, contrast.by))

  if (use_pseudobulk) {
    cache_key <- paste(sort(pb_group), collapse = "+")
    cached_pb <- seurat_object@misc$pseudobulk[[cache_key]]
    n_cached_genes <- if (!is.null(cached_pb))
      length(setdiff(colnames(cached_pb), c(pb_group, "n_cells"))) else 0L

    if (!force && n_cached_genes >= 1000L) {
      # ── 5. Use existing pseudobulk cache ──────────────────────────────────
      message("scSidekick RunURA: Using cached pseudobulk (",
              nrow(cached_pb), " groups, ", n_cached_genes, " genes). ",
              "Pass force = TRUE to recompute.")
    } else {
      # ── 5. Compute pseudobulk ─────────────────────────────────────────────
      message("scSidekick RunURA: Computing pseudobulk (",
              paste(pb_group, collapse = " × "), ")...")
      seurat_object <- ComputePseudobulk(
        seurat_object,
        group.by  = pb_group,
        genes     = NULL,
        assay     = assay,
        layer     = layer,
        min_cells = min_cells,
        force     = force,
        verbose   = FALSE
      )
    }

    pb_df <- seurat_object@misc$pseudobulk[[cache_key]]
    if (is.null(pb_df) || nrow(pb_df) == 0L)
      stop("Pseudobulk returned no groups – check min_cells and metadata.")

    # ── 6. Pseudobulk df → genes × samples matrix ───────────────────────────
    gene_cols <- setdiff(colnames(pb_df), c(pb_group, "n_cells"))
    if (length(gene_cols) < 50L)
      warning("Only ", length(gene_cols), " gene(s) in pseudobulk. ",
              "CollecTRI needs several hundred for robust TF inference.")

    pb_mat           <- t(as.matrix(pb_df[, gene_cols, drop = FALSE]))
    colnames(pb_mat) <- .pb_key(pb_df, pb_group)

    # ── 7. Column metadata (one row per pseudobulk sample) ──────────────────
    col_meta           <- pb_df[, pb_group, drop = FALSE]
    col_meta$pb_col    <- .pb_key(col_meta, pb_group)
    rownames(col_meta) <- NULL

  } else {
    # ── 5-7 (single-cell mode) — run decoupleR on all cells directly ────────
    message("scSidekick RunURA: Using single-cell expression matrix (use_pseudobulk = FALSE)...")

    # Stratified downsampling — preserves group proportions, cuts dense matrix size.
    # ULM scores are per-cell so sampling doesn't bias the activity estimates.
    n_cells_total <- ncol(seurat_object)
    if (!is.null(max_cells) && n_cells_total > max_cells) {
      strat_key <- apply(
        seurat_object@meta.data[, pb_group, drop = FALSE],
        1L, paste, collapse = "\t"
      )
      set.seed(42L)
      keep_cells <- unlist(lapply(split(colnames(seurat_object), strat_key), function(bc) {
        n_keep <- max(1L, round(max_cells * length(bc) / n_cells_total))
        if (length(bc) <= n_keep) bc else sample(bc, n_keep)
      }), use.names = FALSE)
      seurat_object <- seurat_object[, keep_cells]
      message("scSidekick RunURA: Downsampled ", n_cells_total, " → ",
              length(keep_cells), " cells (stratified, max_cells = ", max_cells, ").")
    }

    pb_mat <- tryCatch(
      as.matrix(SeuratObject::LayerData(seurat_object, assay = assay, layer = layer)),
      error = function(e)
        as.matrix(Seurat::GetAssayData(seurat_object, assay = assay, slot = layer))
    )  # genes × cells — dense

    n_genes_sc <- nrow(pb_mat)
    if (n_genes_sc < 50L)
      warning("Only ", n_genes_sc, " gene(s) in matrix. ",
              "CollecTRI needs several hundred for robust TF inference.")

    # col_meta: one row per cell (cell barcodes as pb_col)
    col_meta           <- seurat_object@meta.data[, pb_group, drop = FALSE]
    col_meta$pb_col    <- rownames(seurat_object@meta.data)
    rownames(col_meta) <- NULL
  }

  # ── 8. CollecTRI network ──────────────────────────────────────────────────
  if (is.null(network)) {
    message("scSidekick RunURA: Fetching CollecTRI (", organism, ")...")
    network <- .ura_get_collectri(organism)
  } else {
    # Normalize user-supplied network — handles weight/mor/column name variants
    raw_cols <- if (is.data.frame(network)) colnames(network) else class(network)
    network  <- .ura_normalize_net(network)  # no min_tfs — user controls their own net
    if (is.null(network))
      stop(
        "Supplied 'network' could not be parsed.\n",
        "  Columns present: ", paste(raw_cols, collapse = ", "), "\n",
        "  Need at minimum: source (or source_genesymbol), target (or target_genesymbol), ",
        "and mor (or weight / is_stimulation+is_inhibition).\n\n",
        "To fetch a fresh network when OmniPath is reachable:\n",
        "  net <- decoupleR::get_collectri(organism = 'human')\n",
        "  saveRDS(net, '~/.scSidekick/collectri_human.rds')"
      )
    n_tfs_supplied <- length(unique(network$source))
    if (n_tfs_supplied < 50L)
      warning("Supplied network has only ", n_tfs_supplied, " unique TF(s). ",
              "CollecTRI normally has ~500 — results may be unreliable.")
  }
  net_use <- network[network$target %in% rownames(pb_mat), ]
  n_tfs   <- length(unique(net_use$source))
  message("scSidekick RunURA: ", n_tfs, " TF(s) with targets in this panel.")
  if (n_tfs < 5L)
    stop("Only ", n_tfs, " TF(s) overlap your gene panel — ",
         "check that gene names in the network match rownames(pb_mat). ",
         "Network target sample: ",
         paste(head(network$target, 5), collapse = ", "))

  # ── 9. TF activity inference ──────────────────────────────────────────────
  # decoupleR::check_nas_infs() requires a plain numeric matrix.
  pb_mat <- base::matrix(
    as.numeric(pb_mat),
    nrow = nrow(pb_mat), ncol = ncol(pb_mat),
    dimnames = dimnames(pb_mat)
  )
  pb_mat[!is.finite(pb_mat)] <- 0

  message("scSidekick RunURA: Running ", toupper(method), " on ",
          ncol(pb_mat), " samples × ", nrow(pb_mat), " genes...")

  run_fn <- switch(method,
    ulm   = decoupleR::run_ulm,
    viper = decoupleR::run_viper,
    mlm   = decoupleR::run_mlm
  )
  scores_raw <- run_fn(mat = pb_mat, net = net_use,
                       .source = "source", .target = "target", .mor = "mor",
                       minsize = 5L)
  scores <- scores_raw[scores_raw$statistic == method, ]
  if (nrow(scores) == 0L)
    stop("decoupleR returned no scores for method '", method, "'.")

  # ── 10. TF × sample activity matrix ───────────────────────────────────────
  act_mat <- .ura_scores_to_matrix(scores)
  act_mat <- act_mat[, intersect(colnames(act_mat), col_meta$pb_col), drop = FALSE]

  # ── 10b. Single-cell mode: aggregate per-donor for DE (proper replicates) ─
  # In single-cell mode act_mat is TF × cells.  For differential testing we
  # need biological replicates, so average activity per (donor × group × ...)
  # to get a pseudobulk-of-activities matrix.  mean_act (heatmap) still uses
  # per-cell averages per group — no change there.
  if (!use_pseudobulk && !is.null(donor.by)) {
    agg_key <- apply(
      col_meta[match(colnames(act_mat), col_meta$pb_col), pb_group, drop = FALSE],
      1L, paste, collapse = "\t"
    )
    unique_agg <- unique(agg_key)
    agg_cols   <- lapply(unique_agg, function(k) {
      cells <- colnames(act_mat)[agg_key == k]
      if (length(cells) == 1L) act_mat[, cells, drop = TRUE]
      else rowMeans(act_mat[, cells, drop = FALSE], na.rm = TRUE)
    })
    act_mat_de           <- do.call(cbind, agg_cols)
    colnames(act_mat_de) <- unique_agg

    # Rebuild col_meta to match the aggregated columns
    col_meta_de        <- do.call(rbind, lapply(unique_agg, function(k) {
      row_idx <- which(agg_key == k)[1L]
      col_meta[match(colnames(act_mat)[row_idx], col_meta$pb_col), pb_group, drop = FALSE]
    }))
    col_meta_de$pb_col <- unique_agg
    rownames(col_meta_de) <- NULL
  } else {
    act_mat_de  <- act_mat
    col_meta_de <- col_meta
  }

  # ── 11. Mean activity per (group × split × contrast) combination ──────────
  # Column labels include contrast.by so AD and Control columns are distinct.
  grp_parts <- list(as.character(col_meta[[group.by]]))
  if (!is.null(contrast.by)) grp_parts[[length(grp_parts) + 1L]] <- as.character(col_meta[[contrast.by]])
  if (!is.null(split.by))    grp_parts[[length(grp_parts) + 1L]] <- as.character(col_meta[[split.by]])
  col_meta$grp_label <- do.call(paste, c(grp_parts, list(sep = "\t")))

  grp_levels <- unique(col_meta$grp_label)
  mean_cols  <- lapply(grp_levels, function(g) {
    cols <- intersect(col_meta$pb_col[col_meta$grp_label == g], colnames(act_mat))
    if (length(cols) == 0L) return(NULL)
    if (length(cols) == 1L) act_mat[, cols, drop = TRUE]
    else rowMeans(act_mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  non_null  <- !vapply(mean_cols, is.null, logical(1L))
  mean_act  <- do.call(cbind, mean_cols[non_null])
  colnames(mean_act) <- grp_levels[non_null]

  # ── 12. Differential TF activity ─────────────────────────────────────────
  de_results <- NULL
  if (!is.null(contrast.by) && !is.null(contrast_list)) {
    message("scSidekick RunURA: Differential testing (",
            length(contrast_list), " contrast pair(s))...")
    if (!use_pseudobulk && is.null(donor.by))
      message("  Note: use_pseudobulk = FALSE with no donor.by — cells treated as ",
              "replicates; DE p-values reflect pseudoreplication.")
    de_list <- lapply(contrast_list, function(ctr) {
      .ura_differential(
        act_mat        = act_mat_de,
        col_meta       = col_meta_de,
        group.by       = if (using_dummy) NULL else group.by,
        split.by       = split.by,
        contrast.by    = contrast.by,
        contrast       = ctr,
        has_replicates = !is.null(donor.by)
      )
    })
    de_results <- do.call(rbind, Filter(Negate(is.null), de_list))
    if (!is.null(de_results)) {
      rownames(de_results) <- NULL
      n_sig <- sum(de_results$p.adj < 0.05, na.rm = TRUE)
      message("scSidekick RunURA: ", n_sig, " TF × group row(s) with p.adj < 0.05.")
    }
  }

  # ── 13. Plots ──────────────────────────────────────────────────────────────
  if (isTRUE(plot) && !is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    if (!is.null(de_results) && nrow(de_results) > 0L) {
      # One heatmap per contrast
      for (ctr in unique(de_results$contrast)) {
        .ura_plot_heatmap(
          mean_act      = mean_act,
          n_top         = n_top,
          de            = de_results[de_results$contrast == ctr, ],
          split.by      = split.by,
          group.by      = if (using_dummy) NULL else group.by,
          contrast_name = ctr,
          output_dir    = output_dir,
          object_name   = object_name,
          subset_name   = subset_name
        )
      }
      # Bubble plot (all contrasts together, faceted by contrast)
      .ura_plot_bubble(
        de_tbl      = de_results,
        n_top       = n_top,
        group.by    = if (using_dummy) NULL else group.by,
        split.by    = split.by,
        output_dir  = output_dir,
        object_name = object_name,
        subset_name = subset_name
      )
    } else {
      # No contrast — single activity heatmap
      .ura_plot_heatmap(
        mean_act      = mean_act,
        n_top         = n_top,
        de            = NULL,
        split.by      = split.by,
        group.by      = if (using_dummy) NULL else group.by,
        contrast_name = NULL,
        output_dir    = output_dir,
        object_name   = object_name,
        subset_name   = subset_name
      )
    }
  }

  # ── 14. Return ─────────────────────────────────────────────────────────────
  invisible(list(
    activity = mean_act,
    scores   = scores,
    de       = de_results,
    network  = net_use
  ))
}
