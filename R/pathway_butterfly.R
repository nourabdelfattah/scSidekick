# =============================================================================
# pathway_butterfly.R  --  PlotPathwayButterfly
#
# Cell state hierarchy "butterfly" scatter plot.
#
# Cells are scored against exactly four user-supplied pathways and mapped onto
# a 2-D coordinate system (hierarchy transform) where:
#   Y > 0  → cells prefer the two "top-quadrant" pathways
#   Y < 0  → cells prefer the two "bottom-quadrant" pathways
#   X      → left-right balance within the dominant pair
#
# Output: one large overview panel (all cells colored by group.by) flanked by
# one small panel per group level (group highlighted; rest in gray).  When
# split.by is supplied, a separate figure is produced for each split level.
#
# Scoring logic is adapted from the scrabble R package (Suva Lab) with bug
# fixes applied.  The helpers are prefixed .bf_ and kept internal.
# =============================================================================


# ── scrabble-derived scoring helpers (internal) ───────────────────────────────

.bf_prepare_mat <- function(mat, FUN = rowMeans) {
  x <- if (!is.null(dim(mat))) FUN(mat) else mat
  if (is.numeric(x)) sort(x) else x
}

.bf_name_binids <- function(binids, x) {
  names(binids) <- if (!is.null(names(x))) names(x) else as.character(x)
  binids
}

.bf_bin <- function(mat = NULL, x = NULL, breaks = 30L, FUN = rowMeans) {
  if (!is.null(mat)) x <- .bf_prepare_mat(mat, FUN = FUN)
  if (is.null(x)) stop("Supply mat or x to .bf_bin().")
  binids <- cut(seq_along(x), breaks = breaks,
                labels = FALSE, include.lowest = TRUE)
  .bf_name_binids(binids, x)
}

.bf_match_bins <- function(group, bins) {
  binid <- bins[group]
  sapply(binid, function(id) names(bins)[bins == id], simplify = FALSE)
}

.bf_check_sample_size <- function(bins, n, replace) {
  if (any(lengths(bins) < n) && !replace)
    stop("A bin has fewer than n=", n, " genes and replace=FALSE. ",
         "Reduce n or set replace=TRUE.")
}

.bf_sample_bins <- function(bins, n = 100L, replace = FALSE) {
  .bf_check_sample_size(bins, n, replace)
  unlist(sapply(bins, function(b) sample(b, size = n, replace = replace),
                simplify = FALSE), use.names = FALSE)
}

.bf_col_center <- function(m) scale(as.matrix(m), center = TRUE, scale = FALSE)

.bf_score <- function(mat, groups, center = FALSE) {
  # Coerce to dense matrix to ensure colMeans dispatches correctly regardless
  # of sparse/lazy matrix subclass (mat is already gene-subset, so this is cheap)
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  s_mat <- sapply(groups, function(p) {
    p <- intersect(p, rownames(mat))
    if (length(p) == 0L) return(rep(NA_real_, ncol(mat)))
    colMeans(mat[p, , drop = FALSE])
  })
  rownames(s_mat) <- colnames(mat)
  if (center) s_mat <- .bf_col_center(s_mat)
  s_mat
}

.bf_hierarchy <- function(m, quadrant_order, log_scale = TRUE) {
  dat <- as.data.frame(m[, quadrant_order, drop = FALSE])
  colnames(dat) <- c("bl", "br", "tl", "tr")
  rows <- rownames(m)
  dat <- dplyr::mutate(dat,
    bottom   = pmax(bl, br),
    top      = pmax(tl, tr),
    b.center = br - bl,
    t.center = tr - tl,
    x        = ifelse(bottom > top, b.center, t.center),
    x.scaled = sign(x) * log2(abs(x) + 1),
    y        = top - bottom,
    y.scaled = sign(y) * log2(abs(y) + 1)
  )
  dat <- if (log_scale) dplyr::transmute(dat, X = x.scaled, Y = y.scaled)
         else           dplyr::transmute(dat, X = x,        Y = y)
  rownames(dat) <- rows
  dat
}

.bf_shorten_label <- function(s) {
  s <- sub(
    "^(HALLMARK|KEGG|REACTOME|GOBP|GOMF|GOCC|WP|BIOCARTA|PID|IMMUNESIGDB|VAX)_",
    "", s, ignore.case = TRUE
  )
  tools::toTitleCase(tolower(gsub("_", " ", s)))
}


# ── Build a single butterfly panel ───────────────────────────────────────────

# plot_df    : data.frame with columns X, Y, Group
# highlight  : character(1) group to foreground; NULL = overview (all colored)
# col_map    : named color vector for all groups
# qlabels    : character(4) corner labels in order bl, br, tl, tr
# is_large   : logical — TRUE for the overview panel (bigger text / border)
# title_size : plot title font size (pt)
.bf_panel <- function(plot_df,
                      highlight,
                      col_map,
                      x_lim, y_lim,
                      qlabels,
                      title,
                      point_size,
                      point_alpha,
                      lbl_size,
                      title_size,
                      is_large) {

  x_range <- diff(x_lim)
  y_range <- diff(y_lim)
  box_w   <- x_range * 0.32
  box_h   <- y_range * 0.20

  if (is.null(highlight)) {
    # Overview: all cells colored by group, randomize draw order
    plot_df$disp_col <- col_map[as.character(plot_df$Group)]
    plot_df          <- plot_df[sample(nrow(plot_df)), ]
    pts <- ggplot2::geom_point(
      data  = plot_df,
      ggplot2::aes(x = X, y = Y, color = disp_col),
      size  = point_size,
      alpha = point_alpha
    )
  } else {
    bg_df           <- plot_df[plot_df$Group != highlight, ]
    fg_df           <- plot_df[plot_df$Group == highlight, ]
    fg_df$disp_col  <- col_map[as.character(fg_df$Group)]
    pts <- list(
      ggplot2::geom_point(
        data  = bg_df,
        ggplot2::aes(x = X, y = Y),
        color = "gray85",
        size  = point_size * 0.55,
        alpha = 0.18
      ),
      ggplot2::geom_point(
        data  = fg_df,
        ggplot2::aes(x = X, y = Y, color = disp_col),
        size  = point_size,
        alpha = point_alpha
      )
    )
  }

  ggplot2::ggplot() +
    pts +
    ggplot2::scale_color_identity() +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = x_lim) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = y_lim) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
      panel.border     = ggplot2::element_rect(
                           color     = "black",
                           fill      = NA,
                           linewidth = if (is_large) 1.5 else 1.0),
      plot.title       = ggplot2::element_text(
                           size  = title_size,
                           face  = "bold"),
      axis.ticks       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_blank(),
      legend.position  = "none",
      plot.margin      = ggplot2::margin(4, 4, 4, 4)
    ) +
    # Bottom-left corner
    ggplot2::annotate("rect",
      xmin = x_lim[1],          xmax = x_lim[1] + box_w,
      ymin = y_lim[1],          ymax = y_lim[1] + box_h,
      fill = "black") +
    ggplot2::annotate("text",
      x = x_lim[1] + box_w / 2, y = y_lim[1] + box_h / 2,
      label = qlabels[1], color = "white", fontface = "bold", size = lbl_size) +
    # Bottom-right corner
    ggplot2::annotate("rect",
      xmin = x_lim[2] - box_w,  xmax = x_lim[2],
      ymin = y_lim[1],          ymax = y_lim[1] + box_h,
      fill = "black") +
    ggplot2::annotate("text",
      x = x_lim[2] - box_w / 2, y = y_lim[1] + box_h / 2,
      label = qlabels[2], color = "white", fontface = "bold", size = lbl_size) +
    # Top-left corner
    ggplot2::annotate("rect",
      xmin = x_lim[1],          xmax = x_lim[1] + box_w,
      ymin = y_lim[2] - box_h,  ymax = y_lim[2],
      fill = "black") +
    ggplot2::annotate("text",
      x = x_lim[1] + box_w / 2, y = y_lim[2] - box_h / 2,
      label = qlabels[3], color = "white", fontface = "bold", size = lbl_size) +
    # Top-right corner
    ggplot2::annotate("rect",
      xmin = x_lim[2] - box_w,  xmax = x_lim[2],
      ymin = y_lim[2] - box_h,  ymax = y_lim[2],
      fill = "black") +
    ggplot2::annotate("text",
      x = x_lim[2] - box_w / 2, y = y_lim[2] - box_h / 2,
      label = qlabels[4], color = "white", fontface = "bold", size = lbl_size)
}


# =============================================================================
# PlotPathwayButterfly  – main exported function
# =============================================================================

#' Pathway state hierarchy butterfly plot
#'
#' Scores cells against four pathways, maps them to a 2-D hierarchy coordinate
#' system, and produces a butterfly layout: one large overview panel (all cells
#' colored by \code{group.by}) plus one small panel per group level (group
#' highlighted in color; remaining cells in gray).  When \code{split.by} is
#' supplied, a separate figure is produced for each split level.
#'
#' @section Quadrant layout:
#' \code{quadrant_order} specifies the four pathway names in the order
#' **bottom-left, bottom-right, top-left, top-right**.  The Y-axis separates
#' top from bottom; the X-axis gives the left-right balance within whichever
#' axis dominates for each cell.
#'
#' @section Scoring:
#' Pathway scores are computed as the mean expression across each pathway's
#' genes per cell (adapted from the scrabble R package).  Scores are
#' column-centered when \code{center = TRUE} (default) so that the mean score
#' of each pathway across all cells is zero before the hierarchy transform.
#'
#' @param seurat_object A Seurat object.
#' @param group.by Character. Metadata column for cell grouping and coloring.
#' @param quadrant_order Character vector of exactly four pathway names, in
#'   order: bottom-left, bottom-right, top-left, top-right.  Names must match
#'   either \code{names(gene_sets)} (custom mode) or pathway names in the
#'   requested MSigDB collection.
#' @param split.by Character or \code{NULL}. Metadata column whose levels each
#'   generate a separate butterfly figure.  Default: \code{NULL}.
#' @param gene_sets Named list of character vectors — one element per pathway.
#'   When \code{NULL} (default), gene sets are fetched from MSigDB via
#'   \code{msigdbr} using \code{pathway_sets} and filtered to
#'   \code{quadrant_order}.
#' @param pathway_sets Named list specifying MSigDB collection(s) to search.
#'   Each element is a list with \code{category} and optional
#'   \code{subcategory}.  Only used when \code{gene_sets = NULL}.
#'   Default: \code{list(Hallmark = list(category = "H"))}.
#' @param species Species string for \code{msigdbr}.
#'   Default: \code{"Homo sapiens"}.
#' @param quadrant_labels Character vector of four display labels for the
#'   corner boxes (order matches \code{quadrant_order}: bl, br, tl, tr).
#'   Defaults to auto-shortened pathway names.
#' @param group_colors Named character vector mapping group levels to colors.
#'   Defaults to PrepObject colors or auto-generated palette.
#' @param assay Assay to extract counts from.  Default: \code{"RNA"}.
#' @param layer Layer / slot name.  Default: \code{"counts"}.
#' @param subsample Logical.  When \code{TRUE}, randomly sample \code{n_cells}
#'   cells before scoring.  Default: \code{FALSE}.
#' @param n_cells Integer.  Number of cells to sample when
#'   \code{subsample = TRUE}.  Default: \code{20000L}.
#' @param subsample_warn_at Integer or \code{NULL}.  If the cell count exceeds
#'   this value and \code{subsample = FALSE}, a message is emitted.
#'   Set \code{NULL} to silence.  Default: \code{50000L}.
#' @param center Logical.  Column-center scores before hierarchy transform.
#'   Default: \code{TRUE}.
#' @param log_scale Logical.  Apply log2 scaling to hierarchy coordinates.
#'   Default: \code{TRUE}.
#' @param overview_title Character.  Title for the overview panel.  Defaults
#'   to the object name stored by \code{PrepObject}.
#' @param overview_size Numeric.  Width of the overview panel relative to one
#'   small panel column.  Default: \code{2}.
#' @param panel_ncol Integer or \code{NULL}.  Number of columns in the small-
#'   panel grid.  Defaults to \code{ceiling(n_groups / 2)}.
#' @param point_size Numeric.  Point size.  Default: \code{0.3}.
#' @param point_alpha Numeric.  Point alpha.  Default: \code{0.6}.
#' @param label_size_large Numeric.  Corner label size for the overview panel.
#'   Default: \code{5}.
#' @param label_size_small Numeric.  Corner label size for small panels.
#'   Default: \code{2.5}.
#' @param title_size_large Numeric.  Panel title font size for the overview.
#'   Default: \code{12}.
#' @param title_size_small Numeric.  Panel title font size for small panels.
#'   Default: \code{8}.
#' @param width Numeric or \code{NULL}.  PDF width in inches.  When \code{NULL}
#'   (default), width is auto-sized as \code{(overview_size + panel_ncol) * 1.5}.
#' @param height Numeric or \code{NULL}.  PDF height in inches.  When \code{NULL}
#'   (default), height is auto-sized as \code{panel_nrow * 2.8 + 0.4}.
#' @param file_name Character or \code{NULL}.  Output file name (without
#'   directory).  Defaults to \code{"<object_name> pathway butterfly.pdf"}.
#' @param output_dir Character or \code{NULL}.  Directory for saving the PDF.
#'   Falls back to the PrepObject output directory.  When \code{NULL}, the
#'   plot is printed to the active device.
#' @param caffeinate Logical.  Keep the machine awake during computation.
#'   Default: \code{FALSE}.
#'
#' @return Invisibly, a named list (one element per split level) each
#'   containing \code{hierarchy} (the X/Y coordinate data frame) and
#'   \code{scores} (the raw pathway score matrix, cells × pathways).
#'
#' @export
PlotPathwayButterfly <- function(
    seurat_object,
    group.by,
    quadrant_order,
    split.by          = NULL,
    gene_sets         = NULL,
    pathway_sets      = list(Hallmark = list(category = "H")),
    species           = "Homo sapiens",
    quadrant_labels   = NULL,
    group_colors      = NULL,
    assay             = "RNA",
    layer             = "counts",
    subsample         = FALSE,
    n_cells           = 20000L,
    subsample_warn_at = 50000L,
    center            = TRUE,
    log_scale         = TRUE,
    overview_title    = NULL,
    overview_size     = 2,
    panel_ncol        = NULL,
    point_size        = 0.3,
    point_alpha       = 0.6,
    label_size_large  = 5,
    label_size_small  = 2.5,
    title_size_large  = 12,
    title_size_small  = 8,
    width             = NULL,
    height            = NULL,
    file_name         = NULL,
    output_dir        = NULL,
    caffeinate        = FALSE) {

  if (caffeinate) {
    .caff <- .nk_caffeinate()
    on.exit(.nk_decaffeinate(.caff), add = TRUE)
  }

  output_dir <- output_dir %||%
    if (.nk_autosave(seurat_object)) .nk_setting(seurat_object, "output_dir") else NULL

  # ── Package checks ───────────────────────────────────────────────────────────
  use_msigdb <- is.null(gene_sets)
  if (use_msigdb && !requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required for MSigDB mode. ",
         "Install with install.packages('msigdbr'), or supply gene_sets directly.")

  # ── Validate quadrant_order ──────────────────────────────────────────────────
  if (!is.character(quadrant_order) || length(quadrant_order) != 4L)
    stop("quadrant_order must be a character vector of exactly 4 pathway names ",
         "(order: bottom-left, bottom-right, top-left, top-right).")
  if (anyDuplicated(quadrant_order))
    stop("quadrant_order must contain four distinct pathway names.")

  # ── Load / validate gene sets ────────────────────────────────────────────────
  if (!use_msigdb) {
    if (!is.list(gene_sets) || is.null(names(gene_sets)))
      stop("gene_sets must be a named list of gene vectors.")
    missing_q <- quadrant_order[!quadrant_order %in% names(gene_sets)]
    if (length(missing_q))
      stop("Names in quadrant_order not found in gene_sets: ",
           paste(missing_q, collapse = ", "))
    gene_sets_4 <- gene_sets[quadrant_order]
  } else {
    message("Loading pathway gene sets from MSigDB...")
    all_sets <- do.call(c, unname(lapply(pathway_sets, function(ps) {
      mdf <- .msigdbr_get(
        species     = species,
        category    = ps$category,
        subcategory = ps[["subcategory"]]
      )
      split(mdf$gene_symbol, mdf$gs_name)
    })))
    missing_q <- quadrant_order[!quadrant_order %in% names(all_sets)]
    if (length(missing_q))
      stop("Pathway(s) not found in the requested MSigDB collection: ",
           paste(missing_q, collapse = ", "), "\n",
           "Check quadrant_order names against msigdbr output.")
    gene_sets_4 <- all_sets[quadrant_order]
  }

  # ── Quadrant display labels ──────────────────────────────────────────────────
  if (is.null(quadrant_labels)) {
    quadrant_labels <- stats::setNames(
      sapply(quadrant_order, .bf_shorten_label),
      quadrant_order
    )
  } else {
    if (length(quadrant_labels) != 4L)
      stop("quadrant_labels must have exactly 4 elements.")
    quadrant_labels <- stats::setNames(quadrant_labels, quadrant_order)
  }
  # qlabels vector in quadrant_order slot order: bl, br, tl, tr
  qlabels <- unname(quadrant_labels[quadrant_order])

  # ── Metadata ─────────────────────────────────────────────────────────────────
  meta <- seurat_object@meta.data
  if (!group.by %in% colnames(meta))
    stop("group.by column '", group.by, "' not found in metadata.")
  if (!is.null(split.by) && !split.by %in% colnames(meta))
    stop("split.by column '", split.by, "' not found in metadata.")

  grp_vec  <- as.character(meta[[group.by]])
  grp_lvls <- if (is.factor(meta[[group.by]])) levels(meta[[group.by]])
               else sort(unique(grp_vec))

  # ── Colors ───────────────────────────────────────────────────────────────────
  col_map <- group_colors %||%
    .nk_colors(seurat_object, group.by) %||%
    stats::setNames(
      Nour_pal(if (length(grp_lvls) <= 8) "all" else "spectrum")(length(grp_lvls)),
      grp_lvls
    )

  # ── Extract expression matrix (pathway genes only) ───────────────────────────
  all_genes <- unique(unlist(gene_sets_4, use.names = FALSE))
  message("Extracting counts for ", length(all_genes), " pathway genes...")
  mat_full <- .get_layer_data(seurat_object,
                               assay    = assay,
                               layer    = layer,
                               features = all_genes)

  # ── Split levels ──────────────────────────────────────────────────────────────
  if (!is.null(split.by)) {
    sp_vec  <- as.character(meta[[split.by]])
    sp_lvls <- if (is.factor(meta[[split.by]])) levels(meta[[split.by]])
                else sort(unique(sp_vec))
  } else {
    sp_lvls <- "All"
  }

  # ── Auto panel_ncol and PDF dimensions ───────────────────────────────────────
  n_groups       <- length(grp_lvls)
  panel_ncol_eff <- panel_ncol %||% max(1L, ceiling(n_groups / 2L))
  panel_nrow_eff <- ceiling(n_groups / panel_ncol_eff)
  # 1.5 in per small-panel column (overview counts as overview_size columns);
  # 2.8 in per row + 0.4 in title margin
  eff_width  <- width  %||% (overview_size + panel_ncol_eff) * 1.5
  eff_height <- height %||% panel_nrow_eff * 2.8 + 0.4

  # ── Overview title and file name defaults ─────────────────────────────────────
  obj_name <- .nk_setting(seurat_object, "object_name") %||%
              .nk_setting(seurat_object, "subset_name") %||%
              "Cells"
  overview_title_eff <- overview_title %||% obj_name

  # Quadrant tag for file names: first word of each display label
  q_tag      <- paste(sapply(qlabels, function(s) sub(" .*", "", s)), collapse = "-")
  fname_base <- file_name %||%
                paste0(obj_name, " butterfly by ", group.by, " [", q_tag, "]")
  # strip .pdf suffix if user added it — we append it below
  fname_base <- sub("\\.pdf$", "", fname_base, ignore.case = TRUE)

  # ── Main loop over split levels ───────────────────────────────────────────────
  results <- list()

  for (sp in sp_lvls) {

    # Subset cells
    if (sp == "All") {
      cell_ids <- colnames(mat_full)
    } else {
      cell_ids <- rownames(meta)[as.character(meta[[split.by]]) == sp]
      cell_ids <- intersect(cell_ids, colnames(mat_full))
    }
    if (length(cell_ids) == 0L) {
      message("Skipping split level '", sp, "': no cells found.")
      next
    }

    mat     <- mat_full[, cell_ids, drop = FALSE]
    sp_meta <- meta[cell_ids, , drop = FALSE]

    # Subsample warning / downsampling
    n_sp <- ncol(mat)
    if (!subsample && !is.null(subsample_warn_at) && n_sp > subsample_warn_at) {
      message(
        "PlotPathwayButterfly: ", n_sp, " cells",
        if (sp != "All") paste0(" in '", sp, "'") else "",
        ". Scoring may be slow — consider subsample = TRUE, n_cells = ",
        n_cells, "."
      )
    }
    if (subsample && n_sp > n_cells) {
      keep    <- sample(cell_ids, n_cells, replace = FALSE)
      mat     <- mat[, keep, drop = FALSE]
      sp_meta <- sp_meta[keep, , drop = FALSE]
      message("Subsampled to ", n_cells, " cells",
              if (sp != "All") paste0(" for '", sp, "'") else "", ".")
    }

    # Filter gene sets to genes present in this matrix
    gene_sets_filt <- lapply(gene_sets_4, function(g) intersect(g, rownames(mat)))
    n_overlap <- lengths(gene_sets_filt)
    if (any(n_overlap < 5L))
      warning(
        "Low gene overlap for: ",
        paste(paste0(names(n_overlap)[n_overlap < 5L], " (n=",
                     n_overlap[n_overlap < 5L], ")"), collapse = ", "),
        ". Scores may be unreliable.",
        call. = FALSE
      )

    # Score
    message("Scoring ", ncol(mat), " cells against ",
            length(gene_sets_filt), " pathways...")
    scores <- .bf_score(mat = mat, groups = gene_sets_filt, center = center)

    # Hierarchy coordinates
    h <- .bf_hierarchy(m = scores, quadrant_order = quadrant_order,
                       log_scale = log_scale)

    # Merge Group into plot data frame
    plot_df        <- as.data.frame(h)
    plot_df$Group  <- as.character(sp_meta[rownames(plot_df), group.by])
    plot_df$Group[is.na(plot_df$Group)] <- "Unknown"

    # Axis limits with asymmetric padding (more vertical headroom for corner boxes)
    x_pad <- diff(range(h$X)) * 0.12
    y_pad <- diff(range(h$Y)) * 0.20
    x_lim <- c(min(h$X) - x_pad, max(h$X) + x_pad)
    y_lim <- c(min(h$Y) - y_pad, max(h$Y) + y_pad)

    # ── Build overview panel ─────────────────────────────────────────────────
    ov_title <- if (sp == "All") overview_title_eff
                else paste0(overview_title_eff, "\n", sp)

    overview_plot <- .bf_panel(
      plot_df     = plot_df,
      highlight   = NULL,
      col_map     = col_map,
      x_lim       = x_lim,
      y_lim       = y_lim,
      qlabels     = qlabels,
      title       = ov_title,
      point_size  = point_size,
      point_alpha = point_alpha,
      lbl_size    = label_size_large,
      title_size  = title_size_large,
      is_large    = TRUE
    )

    # ── Build small per-group panels ─────────────────────────────────────────
    small_panels <- lapply(grp_lvls, function(grp) {
      .bf_panel(
        plot_df     = plot_df,
        highlight   = grp,
        col_map     = col_map,
        x_lim       = x_lim,
        y_lim       = y_lim,
        qlabels     = qlabels,
        title       = grp,
        point_size  = point_size,
        point_alpha = point_alpha,
        lbl_size    = label_size_small,
        title_size  = title_size_small,
        is_large    = FALSE
      )
    })

    # Pad small panel list to fill the grid
    n_total <- panel_ncol_eff * panel_nrow_eff
    if (n_total > n_groups) {
      blanks       <- rep(list(.blank_panel()), n_total - n_groups)
      small_panels <- c(small_panels, blanks)
    }

    # ── Arrange layout ────────────────────────────────────────────────────────
    right_grid <- ggpubr::ggarrange(
      plotlist = small_panels,
      ncol     = panel_ncol_eff,
      nrow     = panel_nrow_eff
    )
    final_plot <- ggpubr::ggarrange(
      overview_plot, right_grid,
      ncol   = 2L,
      widths = c(overview_size, panel_ncol_eff)
    )

    # ── Cells per group (computed once, used in legend + JSON) ───────────────
    cells_per_group <- sort(table(sp_meta[[group.by]]), decreasing = TRUE)
    cpg_str <- paste(
      paste0(names(cells_per_group), " (n=", as.integer(cells_per_group), ")"),
      collapse = ", "
    )

    # ── Save or print ─────────────────────────────────────────────────────────
    results[[sp]] <- list(hierarchy = h, scores = scores,
                          cells_per_group = as.list(cells_per_group))

    if (!is.null(output_dir)) {
      sp_suffix <- if (sp == "All") "" else paste0(" - ", sp)
      fpath     <- file.path(output_dir, paste0(fname_base, sp_suffix, ".pdf"))
      grDevices::pdf(fpath, width = eff_width, height = eff_height)
      print(final_plot)
      grDevices::dev.off()

      # ── .legend sidecar ───────────────────────────────────────────────────
      .write_legend_sidecar(fpath, paste0(
        "Butterfly plot of single-cell pathway state hierarchy",
        if (sp != "All") paste0(" — ", split.by, ": ", sp) else "",
        ". Cells scored against four pathways: ",
        paste(paste0(quadrant_order, " (", lengths(gene_sets_filt), " genes)"),
              collapse = "; "), ". ",
        "Scores are mean ", layer, " expression across each pathway's genes",
        if (center) ", column-centered across cells" else "",
        ". Hierarchy coordinates: Y-axis separates ",
        qlabels[3], " / ", qlabels[4], " (top, Y > 0) from ",
        qlabels[1], " / ", qlabels[2], " (bottom, Y < 0); ",
        "X-axis gives left-right balance within the dominant axis",
        if (log_scale) " (log2-scaled)" else "",
        ". Total cells: ", ncol(mat), ". ",
        "Cells per ", group.by, ": ", cpg_str, ". ",
        "Left panel: overview of all cells colored by ", group.by, ". ",
        "Right panels: one panel per ", group.by,
        " level (group highlighted in color; remaining cells in gray)."
      ))

      # ── .json methods sidecar ─────────────────────────────────────────────
      if (requireNamespace("jsonlite", quietly = TRUE)) {
        json_path <- sub("\\.pdf$", ".json", fpath, ignore.case = TRUE)
        methods_text <- paste0(
          ncol(mat), " cells",
          if (sp != "All") paste0(" (", split.by, ": ", sp, ")") else "",
          " were scored against four gene sets — ",
          paste(paste0(quadrant_order, " [", lengths(gene_sets_filt), " genes]"),
                collapse = "; "),
          " — using mean ", layer, " expression across each pathway's genes",
          if (center) " (column-centered across cells)" else "",
          ". Hierarchy coordinates were computed: Y = pmax(tl,tr) - pmax(bl,br); ",
          "X = left-right balance within the dominant axis",
          if (log_scale) " (log2-scaled)" else "",
          ". Quadrant assignment (bl, br, tl, tr): ",
          paste(paste0(qlabels, " = ", quadrant_order), collapse = "; "),
          ". Cells per ", group.by, ": ", cpg_str, "."
        )
        jsonlite::write_json(
          list(
            date              = format(Sys.Date()),
            function_name     = "PlotPathwayButterfly",
            object            = obj_name,
            split_by          = if (is.null(split.by)) NA_character_ else split.by,
            split_level       = if (sp == "All") NA_character_ else sp,
            group_by          = group.by,
            assay             = assay,
            layer             = layer,
            quadrant_order    = as.list(quadrant_order),
            quadrant_labels   = as.list(stats::setNames(qlabels, c("bl","br","tl","tr"))),
            gene_overlap      = as.list(lengths(gene_sets_filt)),
            n_cells_total     = ncol(mat),
            n_cells_per_group = as.list(as.integer(cells_per_group[grp_lvls[grp_lvls %in% names(cells_per_group)]])),
            subsampled        = isTRUE(subsample),
            center            = center,
            log_scale         = log_scale,
            methods_text      = methods_text
          ),
          path        = json_path,
          auto_unbox  = TRUE,
          pretty      = TRUE
        )
        message("Methods JSON: ", json_path)
      }

      message("Saved: ", fpath)
    } else {
      print(final_plot)
      message("Cells per ", group.by, ":\n",
              paste(paste0("  ", names(cells_per_group), ": ",
                           as.integer(cells_per_group)), collapse = "\n"))
    }
  }

  invisible(results)
}
