# =============================================================================
# scSidekick spatial transcriptomics visualization
#
# GenerateSpatialFeatureMaps — batch spatial feature plots (genes OR metadata)
# GenerateSpatialDimMaps     — batch spatial cell-type/cluster dim plots
# GenerateMasterGeneMaps     — per-gene 4-column PDF: UMAP + spatial side-by-side
#
# All functions share the same auto-sizing logic for spatial point sizes
# (median nearest-neighbor distance scaled to plot area) and support both
# Visium V1 (VisiumV1) and Visium HD / Xenium V2 (FOV) architectures.
# =============================================================================

# Internal: compute auto point size from spot coordinates
.auto_pt_size <- function(coords, scale_constant = 100, max_size = 5,
                           min_size = 0.3) {
  robust_w   <- max(coords[, 1]) - min(coords[, 1])
  robust_h   <- max(coords[, 2]) - min(coords[, 2])
  equiv_dim  <- sqrt(robust_w * robust_h)

  n_spots  <- nrow(coords)
  n_sample <- min(n_spots, 500L)
  idx      <- sample.int(n_spots, n_sample)

  min_dists <- vapply(idx, function(j) {
    d_sq    <- (coords[, 1] - coords[j, 1])^2 + (coords[, 2] - coords[j, 2])^2
    d_sq[j] <- Inf
    sqrt(min(d_sq))
  }, numeric(1))

  spots_across <- equiv_dim / stats::median(min_dists)
  max(min(scale_constant / spots_across, max_size), min_size)
}

# Internal: extract coords and cell names from a Seurat image slot
.img_coords <- function(img_obj) {
  if (inherits(img_obj, "FOV")) {
    cen   <- img_obj@boundaries[["centroids"]]
    list(coords = cen@coords, cells = cen@cells)
  } else if (inherits(img_obj, "VisiumV1")) {
    cd    <- img_obj@coordinates
    list(coords = as.matrix(cd[, c("imagerow", "imagecol")]),
         cells  = rownames(cd))
  } else {
    stop("Unrecognised image class: ", class(img_obj))
  }
}

# --------------------------------------------------------------------------- #
# GenerateSpatialFeatureMaps                                                   #
# --------------------------------------------------------------------------- #

#' Batch spatial feature plots for genes and/or metadata columns
#'
#' Loops over `features` (gene names or numeric metadata columns) and
#' produces [Seurat::SpatialFeaturePlot()] panels for each image in the
#' object. A shared color scale (0 → max) is applied across all images for
#' the same feature. Panels are arranged in a grid and a shared legend is
#' attached.
#'
#' Supports two layout modes:
#' \describe{
#'   \item{`"auto"`}{Automatic 1- or 2-row grid.}
#'   \item{`"metadata"`}{Groups images into rows based on a slide-level
#'     metadata variable (`row_var`). The variable must have a single
#'     unique value per image.}
#' }
#'
#' @param seurat_object A Seurat object with spatial images.
#' @param features Character vector of gene names or numeric metadata
#'   column names.
#' @param layout_method Character. `"auto"` (default) or `"metadata"`.
#' @param row_var Character or `NULL`. Slide-level metadata column for
#'   grouping images into rows (required when `layout_method = "metadata"`).
#' @param images_to_plot Character vector or `NULL`. Subset of image names
#'   to include. `NULL` uses all images.
#' @param remove_outliers Logical. Crop outlier cells before plotting.
#'   Default `FALSE`.
#' @param outlier_prob Numeric. Quantile to trim when `remove_outliers = TRUE`.
#'   Default `0.01`.
#' @param size_override Scalar or named numeric list. A plain number (e.g.
#'   \code{size_override = 1.5}) applies to all images; a named list applies
#'   per image. Overrides auto-sizing.
#' @param colors Character vector. Gradient colors from low to high.
#' @param calcptsizesc Numeric. Scaling constant for auto point-size
#'   calculation. Default \code{100}.
#' @param output_dir Character or `NULL`. Output directory for PDFs. Returns
#'   last plot when `NULL`.
#' @param object_name Character. Label appended to PDF filenames. Default `""`.
#'
#' @return Invisibly returns the last plot when `output_dir = NULL`.
#' @export
GenerateSpatialFeatureMaps <- function(seurat_object,
                                        features,
                                        layout_method  = "auto",
                                        row_var        = NULL,
                                        images_to_plot = NULL,
                                        remove_outliers = FALSE,
                                        outlier_prob   = 0.01,
                                        size_override  = NULL,
                                        colors         = c("#053061", "#2166AC",
                                                           "#D1E5F0", "#FDDBC7",
                                                           "#F4A582", "#D6604D",
                                                           "#B2182B", "#67001F"),
                                        calcptsizesc   = 3800,
                                        output_dir     = NULL,
                                        object_name    = "",
                                        subset_name    = "") {

  if (!layout_method %in% c("auto", "metadata"))
    stop("layout_method must be 'auto' or 'metadata'.")
  if (layout_method == "metadata" && is.null(row_var))
    stop("'row_var' must be specified for layout_method = 'metadata'.")

  # Walk up to PrepObject-stored defaults when not explicitly supplied
  output_dir  <- output_dir  %||% .nk_setting(seurat_object, "output_dir")
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  all_imgs  <- Seurat::Images(seurat_object)
  img_names <- if (!is.null(images_to_plot)) {

    intersect(images_to_plot, all_imgs)
  } else all_imgs
  if (length(img_names) == 0) stop("No matching images found.")

  last_plot <- NULL

  for (i in seq_along(features)) {
    feat <- features[i]

    # Detect if feature is a gene or a metadata column
    if (feat %in% rownames(seurat_object)) {
      feat_vals        <- .get_layer_data(seurat_object, layer = "data", features = feat)[feat, ]
      is_numeric_feat  <- TRUE
    } else if (feat %in% colnames(seurat_object@meta.data)) {
      feat_vals        <- seurat_object@meta.data[[feat]]
      is_numeric_feat  <- is.numeric(feat_vals)
    } else {
      warning("Feature '", feat, "' not found — skipping.")
      next
    }

    shared_scale <- if (is_numeric_feat) {
      max_exp <- max(feat_vals, na.rm = TRUE)
      col_lim <- round(max(max_exp, 0.01), digits = 2)
      ggplot2::scale_fill_gradientn(colors = colors,
                                    limits = c(0, col_lim),
                                    na.value = "transparent")
    } else NULL

    plot_list <- list()

    for (img in img_names) {
      img_obj  <- seurat_object@images[[img]]
      ic       <- .img_coords(img_obj)
      coords   <- ic$coords
      cells    <- ic$cells

      if (layout_method == "metadata") {
        vals <- unique(as.character(seurat_object@meta.data[cells, row_var]))
        if (length(vals) > 1)
          stop("Image '", img, "': '", row_var, "' has multiple values. ",
               "row_var must be a slide-level variable.")
      }

      if (remove_outliers) {
        qX      <- stats::quantile(coords[, 1], c(outlier_prob, 1 - outlier_prob))
        qY      <- stats::quantile(coords[, 2], c(outlier_prob, 1 - outlier_prob))
        keep    <- coords[, 1] >= qX[1] & coords[, 1] <= qX[2] &
          coords[, 2] >= qY[1] & coords[, 2] <= qY[2]
        cells   <- cells[keep]
        coords  <- coords[keep, , drop = FALSE]
        plot_obj <- subset(seurat_object, cells = cells)
      } else {
        plot_obj <- seurat_object
      }

      pt_size <- if (!is.null(size_override)) {
        if (length(size_override) == 1L && is.null(names(size_override))) as.numeric(size_override)
        else if (img %in% names(size_override)) size_override[[img]]
        else .auto_pt_size(coords, scale_constant = calcptsizesc)
      } else {
        .auto_pt_size(coords, scale_constant = calcptsizesc)
      }

      p <- Seurat::SpatialFeaturePlot(
        plot_obj, features = feat, images = img,
        image.alpha = 0.5, pt.size.factor = pt_size,
        crop = TRUE, alpha = c(1, 1)
      ) +
        ggplot2::ggtitle(img) +
        ggplot2::theme(plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 14))

      if (!is.null(shared_scale)) p <- p + shared_scale
      plot_list[[img]] <- p
    }

    # Extract legend
    legend_guide <- if (is_numeric_feat) {
      ggplot2::guides(fill = ggplot2::guide_colourbar(
        title = feat, title.position = "top",
        title.theme = ggplot2::element_text(size = 12, face = "bold")
      ))
    } else {
      ggplot2::guides(fill = ggplot2::guide_legend(
        title = feat, title.position = "top",
        title.theme = ggplot2::element_text(size = 12, face = "bold"),
        override.aes = list(size = 5)
      ))
    }
    PlotLegend <- ggpubr::get_legend(plot_list[[1]] + legend_guide)

    plot_list <- lapply(plot_list, function(pp)
      pp + ggplot2::theme_void() + Seurat::NoLegend() +
        ggplot2::theme(plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 14)))

    # Assemble grid
    if (layout_method == "metadata") {
      row_plot_lists <- list()
      for (img in names(plot_list)) {
        first_cell <- .img_coords(seurat_object@images[[img]])$cells[1]
        r_val      <- as.character(seurat_object@meta.data[first_cell, row_var])
        row_plot_lists[[r_val]] <- c(row_plot_lists[[r_val]],
                                     list(plot_list[[img]]))
      }
      max_row_len <- max(lengths(row_plot_lists))
      list_of_rows <- list()
      valid_rows   <- 0
      for (r in names(row_plot_lists)) {
        rp <- row_plot_lists[[r]]
        if (length(rp) == 0) next
        valid_rows <- valid_rows + 1
        while (length(rp) < max_row_len)
          rp[[length(rp) + 1]] <- ggplot2::ggplot() + ggplot2::theme_void()
        single_row <- ggpubr::ggarrange(plotlist = rp, nrow = 1, ncol = max_row_len)
        list_of_rows[[length(list_of_rows) + 1]] <- ggpubr::annotate_figure(
          single_row,
          left = ggpubr::text_grob(r, rot = 90, face = "bold", size = 16)
        )
      }
      n_rows    <- valid_rows
      n_cols    <- max_row_len
      inner_grid <- ggpubr::ggarrange(plotlist = list_of_rows, nrow = n_rows, ncol = 1)
    } else {
      n_plots   <- length(plot_list)
      n_rows    <- if (n_plots <= 4) 1L else 2L
      n_cols    <- ceiling(n_plots / n_rows)
      inner_grid <- ggpubr::ggarrange(plotlist = plot_list, nrow = n_rows, ncol = n_cols)
    }

    CustomPlot <- ggpubr::ggarrange(inner_grid, PlotLegend,
                                     nrow = 1, ncol = 2,
                                     widths = c(n_cols, 0.5))

    if (!is.null(output_dir)) {
      pdf_h <- n_rows * 5
      pdf_w <- (n_cols * 4) + 2
      fname <- file.path(output_dir,
                         paste0(feat, " SpatialFeatureMap ",
                                object_name,
                                if (nchar(subset_name) > 0) paste0(" ", subset_name) else "",
                                ".pdf"))
      grDevices::pdf(fname, width = pdf_w, height = pdf_h)
      print(CustomPlot)
      grDevices::dev.off()
      .write_legend_sidecar(fname, paste0(
        "Spatial feature plot showing the distribution of ", feat,
        " across ", length(img_names), " tissue section(s)",
        if (nchar(object_name) > 0) paste0(" in ", object_name) else "", ". ",
        "Each spot on the tissue section represents a spatial barcode, coloured ",
        "on a continuous scale from low (dark blue) to high (dark red) according ",
        "to the ", if (feat %in% rownames(seurat_object)) "log-normalised expression"
                   else "value", " of ", feat, ".",
        if (remove_outliers) paste0(
          " Outlier spots above the ", round((1 - outlier_prob) * 100),
          "th percentile are excluded from the colour scale to improve ",
          "visualisation of the main expression range.") else ""
      ))
    } else {
      last_plot <- CustomPlot
    }
    message(sprintf("Spatial feature map: %s (%d of %d) — %dx%d",
                    feat, i, length(features), n_rows, n_cols))
  }
  invisible(last_plot)
}

# --------------------------------------------------------------------------- #
# GenerateSpatialDimMaps                                                       #
# --------------------------------------------------------------------------- #

#' Batch spatial dimensionality-reduction (cell-type) plots
#'
#' Loops over `group_by_vars` and produces [Seurat::SpatialDimPlot()] panels
#' for each image. Supports the same `"auto"` / `"metadata"` layout modes as
#' [GenerateSpatialFeatureMaps()].
#'
#' @param seurat_object A Seurat object with spatial images.
#' @param group_by_vars Character vector. Metadata columns to plot (one PDF
#'   per variable).
#' @param layout_method Character. `"auto"` or `"metadata"`. Default `"auto"`.
#' @param row_var Character or `NULL`. Slide-level metadata column for row
#'   grouping (required for `"metadata"` layout).
#' @param images_to_plot Character vector or `NULL`. Images to include.
#' @param remove_outliers Logical. Trim outlier cells. Default `FALSE`.
#' @param outlier_prob Numeric. Quantile trim level. Default `0.01`.
#' @param size_override Scalar or named numeric list. A plain number applies
#'   to all images; a named list applies per image.
#' @param cluster_colors Named character vector. Colors for each group level.
#'   `NULL` uses Seurat defaults.
#' @param output_dir Character or `NULL`. PDF output directory. Returns last
#'   plot when `NULL`.
#' @param rowannsize Numeric. Row annotation text size. Default `16`.
#' @param imgalpha Numeric. Image transparency (0–1). Default `1`.
#' @param alpha Numeric. Point transparency (0–1). Default `1`.
#' @param uniform_size Logical. When `TRUE` all images receive the same
#'   auto-computed point size (the median across images), preventing tiny
#'   tissue fragments from getting huge dots. Default `FALSE`.
#' @param object_name Character. Label for PDF filenames. Default `""`.
#' @param group.by Alias for `group_by_vars` (preferred name).
#'
#' @return Invisibly returns the last plot when `output_dir = NULL`.
#' @export
GenerateSpatialDimMaps <- function(seurat_object,
                                    group_by_vars  = NULL,
                                    layout_method  = "auto",
                                    row_var        = NULL,
                                    images_to_plot = NULL,
                                    remove_outliers = FALSE,
                                    outlier_prob   = 0.01,
                                    size_override  = NULL,
                                    uniform_size   = FALSE,
                                    colors         = NULL,
                                    output_dir     = NULL,
                                    rowannsize     = 16,
                                    imgalpha       = 1,
                                    alpha          = 1,
                                    object_name    = "",
                                    # preferred aliases
                                    group.by       = NULL,
                                    cluster_colors = NULL) {   # deprecated → use colors

  # Accept group.by as the primary alias
  if (is.null(group_by_vars) && !is.null(group.by)) group_by_vars <- group.by
  # Deprecated alias
  if (!is.null(cluster_colors) && is.null(colors)) colors <- cluster_colors
  if (is.null(group_by_vars)) group_by_vars <- "seurat_clusters"

  if (!layout_method %in% c("auto", "metadata"))
    stop("layout_method must be 'auto' or 'metadata'.")
  if (layout_method == "metadata" && is.null(row_var))
    stop("'row_var' must be specified for layout_method = 'metadata'.")

  # Walk up to PrepObject-stored defaults when not explicitly supplied
  output_dir  <- output_dir  %||% .nk_setting(seurat_object, "output_dir")
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  all_imgs  <- Seurat::Images(seurat_object)
  img_names <- if (!is.null(images_to_plot)) intersect(images_to_plot, all_imgs) else all_imgs
  if (length(img_names) == 0) stop("No matching images found.")

  last_plot <- NULL

  for (i in seq_along(group_by_vars)) {
    grp <- group_by_vars[i]

    # ── Resolve colors once per variable ──────────────────────────────────
    grp_colors <- colors %||% .nk_colors(seurat_object, grp)

    # Determine all levels (for factor-setting and legend)
    all_lvls <- if (!is.null(grp_colors)) {
      names(grp_colors)
    } else {
      sort(unique(as.character(seurat_object@meta.data[[grp]])))
    }

    # ── Pre-compute point sizes (optionally uniform across images) ─────────
    raw_sizes <- vapply(img_names, function(img) {
      if (!is.null(size_override)) {
        if (length(size_override) == 1L && is.null(names(size_override))) {
          as.numeric(size_override)                    # scalar: apply to all
        } else if (img %in% names(size_override)) {
          size_override[[img]]                         # named: per-image
        } else {
          .auto_pt_size(.img_coords(seurat_object@images[[img]])$coords)
        }
      } else {
        .auto_pt_size(.img_coords(seurat_object@images[[img]])$coords)
      }
    }, numeric(1))
    if (uniform_size) raw_sizes[] <- stats::median(raw_sizes)
    names(raw_sizes) <- img_names

    # ── Build legend from the full color vector (not from a subset plot) ───
    # This guarantees every level appears in the legend regardless of which
    # cell types happen to be present in the first image.
    lgd_dat   <- data.frame(
      x     = 1L,
      y     = seq_along(all_lvls),
      label = factor(all_lvls, levels = all_lvls)
    )
    lgd_colors <- if (!is.null(grp_colors)) grp_colors else
      stats::setNames(
        Nour_pal(if (length(all_lvls) <= 8) "all" else "spectrum")(length(all_lvls)),
        all_lvls
      )
    lgd_plot <- ggplot2::ggplot(lgd_dat,
        ggplot2::aes(x = x, y = y, fill = label)) +
      ggplot2::geom_point(shape = 21, size = 5) +
      ggplot2::scale_fill_manual(values = lgd_colors) +
      ggplot2::guides(fill = ggplot2::guide_legend(
        title          = grp,
        override.aes   = list(size = 5),
        title.theme    = ggplot2::element_text(size = 15, face = "bold"),
        title.position = "top",
        label.theme    = ggplot2::element_text(size = 15),
        ncol           = 9
      )) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "bottom")
    PlotLegend <- ggpubr::get_legend(lgd_plot)

    plot_list <- list()

    for (img in img_names) {
      img_obj <- seurat_object@images[[img]]
      ic      <- .img_coords(img_obj)
      coords  <- ic$coords
      cells   <- ic$cells

      if (layout_method == "metadata") {
        vals <- unique(as.character(seurat_object@meta.data[cells, row_var]))
        if (length(vals) > 1)
          stop("Image '", img, "': row_var has multiple values.")
      }

      if (remove_outliers) {
        qX     <- stats::quantile(coords[, 1], c(outlier_prob, 1 - outlier_prob))
        qY     <- stats::quantile(coords[, 2], c(outlier_prob, 1 - outlier_prob))
        keep   <- coords[, 1] >= qX[1] & coords[, 1] <= qX[2] &
          coords[, 2] >= qY[1] & coords[, 2] <= qY[2]
        cells  <- cells[keep]
        plot_obj <- subset(seurat_object, cells = cells)
      } else {
        plot_obj <- seurat_object
      }

      # Ensure all factor levels are present in the (possibly subsetted) object
      # so SpatialDimPlot's internal scale covers all categories
      plot_obj@meta.data[[grp]] <- factor(
        as.character(plot_obj@meta.data[[grp]]),
        levels = all_lvls
      )

      pt_size <- raw_sizes[[img]]

      # Pass colors via `cols` to avoid "scale already present" warnings
      p_args <- list(
        object         = plot_obj,
        group.by       = grp,
        images         = img,
        image.alpha    = imgalpha,
        pt.size.factor = pt_size,
        crop           = TRUE,
        alpha          = alpha
      )
      if (!is.null(grp_colors)) p_args[["cols"]] <- grp_colors
      p <- do.call(Seurat::SpatialDimPlot, p_args) +
        ggplot2::ggtitle(img) +
        ggplot2::theme(plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 14))

      plot_list[[img]] <- p
    }

    plot_list <- lapply(plot_list, function(pp)
      pp + ggplot2::theme_void() + Seurat::NoLegend() +
        ggplot2::theme(plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 14)))

    if (layout_method == "metadata") {
      row_plot_lists <- list()
      for (img in names(plot_list)) {
        first_cell <- .img_coords(seurat_object@images[[img]])$cells[1]
        r_val      <- as.character(seurat_object@meta.data[first_cell, row_var])
        row_plot_lists[[r_val]] <- c(row_plot_lists[[r_val]], list(plot_list[[img]]))
      }
      max_row_len  <- max(lengths(row_plot_lists))
      list_of_rows <- list()
      valid_rows   <- 0
      for (r in names(row_plot_lists)) {
        rp <- row_plot_lists[[r]]
        if (length(rp) == 0) next
        valid_rows <- valid_rows + 1
        while (length(rp) < max_row_len)
          rp[[length(rp) + 1]] <- ggplot2::ggplot() + ggplot2::theme_void()
        single_row <- ggpubr::ggarrange(plotlist = rp, nrow = 1, ncol = max_row_len)
        list_of_rows[[length(list_of_rows) + 1]] <- ggpubr::annotate_figure(
          single_row,
          left = ggpubr::text_grob(r, rot = 90, face = "bold", size = rowannsize)
        )
      }
      n_rows    <- valid_rows
      n_cols    <- max_row_len
      inner_grid <- ggpubr::ggarrange(plotlist = list_of_rows, nrow = n_rows, ncol = 1)
    } else {
      n_plots   <- length(plot_list)
      n_rows    <- if (n_plots <= 4) 1L else 2L
      n_cols    <- ceiling(n_plots / n_rows)
      inner_grid <- ggpubr::ggarrange(plotlist = plot_list, nrow = n_rows, ncol = n_cols)
    }

    CustomPlot <- ggpubr::ggarrange(inner_grid, PlotLegend,
                                     nrow = 2, ncol = 1,
                                     heights = c(1, 0.2))

    if (!is.null(output_dir)) {
      pdf_h <- (n_rows * 4) + 1.5
      pdf_w <- max(n_cols * 4, 12)
      fname <- file.path(output_dir,
                         paste0(grp, " SpatialDimMap ", object_name, ".pdf"))
      grDevices::pdf(fname, width = pdf_w, height = pdf_h)
      print(CustomPlot)
      grDevices::dev.off()
      .write_legend_sidecar(fname, paste0(
        "Spatial dimension map showing ", grp, " assignments projected onto ",
        "tissue section coordinates across ", length(img_names), " section(s)",
        if (nchar(object_name) > 0) paste0(" in ", object_name) else "", ". ",
        "Each spot is coloured according to its group identity, enabling direct ",
        "comparison of transcriptional cluster boundaries with tissue anatomy. ",
        "The shared colour legend is displayed beneath the section grid."
      ))
    } else {
      last_plot <- CustomPlot
    }
    message(sprintf("Spatial DimMap: %s (%d of %d) — %dx%d",
                    grp, i, length(group_by_vars), n_rows, n_cols))
  }
  invisible(last_plot)
}

# --------------------------------------------------------------------------- #
# GenerateMasterGeneMaps                                                       #
# --------------------------------------------------------------------------- #

#' 4-column master gene map: UMAP clusters + UMAP expression + Spatial clusters + Spatial expression
#'
#' For each gene, generates one PDF with one row per spatial image, and four
#' columns: UMAP colored by clusters, UMAP colored by gene expression, spatial
#' plot colored by clusters, and spatial plot colored by gene expression. A
#' single pair of legends (categorical + continuous) is attached to the right.
#'
#' Accepts a vector for \code{group.by}: one PDF set is produced per variable.
#'
#' @param seurat_object A Seurat object with both \code{"umap"} reduction and
#'   spatial images.
#' @param features Character vector of gene names.
#' @param imgalpha Numeric. Tissue image transparency. Default \code{1}.
#' @param alpha Numeric. Spot transparency. Default \code{1}.
#' @param group.by Character (or vector). Metadata column(s) for cluster
#'   coloring. Default \code{"seurat_clusters"}. When a vector is supplied,
#'   the function recurses and produces one PDF set per variable.
#' @param images_to_plot Character vector or \code{NULL}. Images to include.
#' @param remove_outliers Logical. Default \code{FALSE}.
#' @param outlier_prob Numeric. Default \code{0.01}.
#' @param size_override Scalar or named numeric list. A plain number applies
#'   to all images; a named list applies per image.
#' @param colors Named character vector. Colors for \code{group.by} levels.
#'   \code{NULL} auto-resolves from PrepObject or \code{Nour_pal}.
#' @param gene_colors Character vector. Gradient for gene expression.
#' @param output_dir Character. Output directory for PDFs.
#' @param object_name Character. Label for PDF filenames. Default \code{""}.
#' @param cluster_col Deprecated alias for \code{group.by}.
#' @param cluster_colors Deprecated alias for \code{colors}.
#'
#' @return \code{NULL} invisibly (always writes to disk).
#' @export
GenerateMasterGeneMaps <- function(seurat_object,
                                    features,
                                    imgalpha        = 1,
                                    alpha           = 1,
                                    group.by        = "seurat_clusters",
                                    images_to_plot  = NULL,
                                    remove_outliers = FALSE,
                                    outlier_prob    = 0.01,
                                    size_override   = NULL,
                                    colors          = NULL,
                                    gene_colors     = c("#053061", "#2166AC",
                                                        "#D1E5F0", "#FDDBC7",
                                                        "#F4A582", "#D6604D",
                                                        "#B2182B", "#67001F"),
                                    output_dir,
                                    object_name     = "",
                                    # deprecated aliases
                                    cluster_col     = NULL,
                                    cluster_colors  = NULL) {

  # ── Deprecated aliases ─────────────────────────────────────────────────────
  if (!is.null(cluster_col)    && identical(group.by, "seurat_clusters"))
    group.by <- cluster_col
  if (!is.null(cluster_colors) && is.null(colors))
    colors <- cluster_colors

  # ── Walk-up PrepObject defaults ────────────────────────────────────────────
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  # ── Multiple group.by: recurse once per variable ───────────────────────────
  if (length(group.by) > 1L) {
    for (gb in group.by) {
      GenerateMasterGeneMaps(
        seurat_object   = seurat_object,
        features        = features,
        imgalpha        = imgalpha,
        alpha           = alpha,
        group.by        = gb,
        images_to_plot  = images_to_plot,
        remove_outliers = remove_outliers,
        outlier_prob    = outlier_prob,
        size_override   = size_override,
        colors          = colors,
        gene_colors     = gene_colors,
        output_dir      = output_dir,
        object_name     = object_name
      )
    }
    return(invisible(NULL))
  }

  # ── Resolve group colors ───────────────────────────────────────────────────
  # Use local alias `grp_colors` to avoid shadowing ggplot2 `colors` argument
  grp_colors <- colors %||% .nk_colors(seurat_object, group.by)

  all_imgs  <- Seurat::Images(seurat_object)
  img_names <- if (!is.null(images_to_plot)) intersect(images_to_plot, all_imgs) else all_imgs
  if (length(img_names) == 0) stop("No matching images found.")

  for (gene_idx in seq_along(features)) {
    gene       <- features[gene_idx]
    max_exp    <- max(.get_layer_data(seurat_object, layer = "data", features = gene)[gene, ])
    color_lim  <- round(max(max_exp, 0.01), digits = 2)

    umap_gene_scale    <- ggplot2::scale_color_gradientn(
      colors = gene_colors, limits = c(0, color_lim), na.value = "lightgrey")
    spatial_gene_scale <- ggplot2::scale_fill_gradientn(
      colors = gene_colors, limits = c(0, color_lim), na.value = "transparent")

    list_of_rows  <- list()
    ClusterLegend <- NULL
    GeneLegend    <- NULL

    for (i in seq_along(img_names)) {
      img     <- img_names[i]
      ic      <- .img_coords(seurat_object@images[[img]])
      coords  <- ic$coords
      cells   <- ic$cells

      if (remove_outliers) {
        qX      <- stats::quantile(coords[, 1], c(outlier_prob, 1 - outlier_prob))
        qY      <- stats::quantile(coords[, 2], c(outlier_prob, 1 - outlier_prob))
        keep    <- coords[, 1] >= qX[1] & coords[, 1] <= qX[2] &
          coords[, 2] >= qY[1] & coords[, 2] <= qY[2]
        cells  <- cells[keep]
        coords <- coords[keep, , drop = FALSE]
      }
      plot_obj <- subset(seurat_object, cells = cells)

      pt_size <- if (!is.null(size_override)) {
        if (length(size_override) == 1L && is.null(names(size_override))) as.numeric(size_override)
        else if (img %in% names(size_override)) size_override[[img]]
        else .auto_pt_size(coords)
      } else {
        .auto_pt_size(coords)
      }

      p1 <- Seurat::DimPlot(plot_obj, reduction = "umap",
                            group.by = group.by, label = TRUE) +
        theme_NourMin() + Seurat::NoAxes() + Seurat::NoLegend()
      if (!is.null(grp_colors))
        p1 <- p1 + ggplot2::scale_color_manual(values = grp_colors, drop = FALSE)
      if (i == 1)
        p1 <- p1 + ggplot2::ggtitle("UMAP: Clusters") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

      p2 <- Seurat::FeaturePlot(plot_obj, features = gene,
                                reduction = "umap", order = TRUE) +
        umap_gene_scale + theme_NourMin() + Seurat::NoAxes() + Seurat::NoLegend()
      if (i == 1)
        p2 <- p2 + ggplot2::ggtitle("UMAP: Expression") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

      # Pass colors via cols= to avoid duplicate scale warnings
      p3_args <- list(object = plot_obj, group.by = group.by, images = img,
                      image.alpha = imgalpha, pt.size.factor = pt_size,
                      crop = TRUE, alpha = alpha)
      if (!is.null(grp_colors)) p3_args[["cols"]] <- grp_colors
      p3 <- do.call(Seurat::SpatialDimPlot, p3_args) +
        ggplot2::theme_void() + Seurat::NoLegend()
      if (i == 1)
        p3 <- p3 + ggplot2::ggtitle("Spatial: Clusters") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

      p4 <- Seurat::SpatialFeaturePlot(plot_obj, features = gene,
                                       images = img, image.alpha = imgalpha,
                                       pt.size.factor = pt_size,
                                       crop = TRUE, alpha = alpha) +
        spatial_gene_scale + ggplot2::theme_void() + Seurat::NoLegend()
      if (i == 1)
        p4 <- p4 + ggplot2::ggtitle("Spatial: Expression") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

      if (i == 1) {
        dummy_cat <- Seurat::DimPlot(plot_obj, group.by = group.by) +
          ggplot2::guides(color = ggplot2::guide_legend(
            title          = group.by,
            override.aes   = list(size = 5),
            title.theme    = ggplot2::element_text(size = 14, face = "bold"),
            label.theme    = ggplot2::element_text(size = 12)))
        if (!is.null(grp_colors))
          dummy_cat <- dummy_cat +
          ggplot2::scale_color_manual(values = grp_colors, drop = FALSE)
        ClusterLegend <- ggpubr::get_legend(dummy_cat)

        dummy_cont <- Seurat::FeaturePlot(plot_obj, features = gene) +
          umap_gene_scale +
          ggplot2::guides(color = ggplot2::guide_colourbar(
            title          = gene,
            title.position = "top",
            title.theme    = ggplot2::element_text(size = 14, face = "bold"),
            barheight      = ggplot2::unit(3, "in")))
        GeneLegend <- ggpubr::get_legend(dummy_cont)
      }

      row_grid <- ggpubr::ggarrange(p1, p2, p3, p4, nrow = 1, ncol = 4)
      list_of_rows[[img]] <- ggpubr::annotate_figure(
        row_grid,
        left = ggpubr::text_grob(img, rot = 90, face = "bold", size = 16)
      )
    }

    n_rows       <- length(list_of_rows)
    main_grid    <- ggpubr::ggarrange(plotlist = list_of_rows, nrow = n_rows, ncol = 1)
    legends_comb <- ggpubr::ggarrange(ClusterLegend, GeneLegend,
                                       nrow = 2, ncol = 1, heights = c(1, 1))
    CustomPlot   <- ggpubr::ggarrange(main_grid, legends_comb,
                                       nrow = 1, ncol = 2, widths = c(4, 0.5))
    final_plot   <- ggpubr::annotate_figure(
      CustomPlot,
      top = ggpubr::text_grob(paste0("Gene Master Map: ", gene),
                              face = "bold", size = 22, color = "#B2182B")
    )

    fname <- file.path(output_dir,
                       paste0(gene, "_", group.by,
                              " Master Gene Map ", object_name, ".pdf"))
    grDevices::pdf(fname, width = 18, height = (n_rows * 3.5) + 1)
    print(final_plot)
    grDevices::dev.off()
    .write_legend_sidecar(fname, paste0(
      "Master spatial gene map showing the co-localisation of ", gene,
      " expression with ", group.by, " identity across ", length(img_names),
      " tissue section(s)",
      if (nchar(object_name) > 0) paste0(" in ", object_name) else "", ". ",
      "The leftmost panel of each row shows spot-level group assignments (",
      group.by, "); subsequent panels display ", gene,
      " expression on a continuous colour scale, facilitating direct spatial ",
      "comparison of gene expression patterns with tissue architecture."
    ))
    message(sprintf("Master gene map: %s / %s (%d of %d)",
                    gene, group.by, gene_idx, length(features)))
  }
  invisible(NULL)
}

# =============================================================================
# Short-name aliases  (preferred calling convention)
# =============================================================================

#' @describeIn GenerateSpatialFeatureMaps Preferred short alias.
#' @export
PlotSpatialFeaturePlots <- GenerateSpatialFeatureMaps

#' @describeIn GenerateSpatialDimMaps Preferred short alias.
#' @export
PlotSpatialDimPlots <- GenerateSpatialDimMaps

#' @describeIn GenerateMasterGeneMaps Preferred short alias.
#' @export
PlotMasterMaps <- GenerateMasterGeneMaps
