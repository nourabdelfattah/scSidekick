# =============================================================================
# scSidekick spatial transcriptomics visualization
#
# PlotSpatialFeaturePlots - batch spatial feature plots (genes OR metadata)
# PlotSpatialDimPlots     - batch spatial cell-type/cluster dim plots
# PlotMasterMaps          - per-gene 4-column PDF: UMAP + spatial side-by-side
#
# All functions share the same auto-sizing logic for spatial point sizes
# (median nearest-neighbor distance scaled to plot area) and support both
# Visium V1 (VisiumV1) and Visium HD / Xenium V2 (FOV) architectures.
# =============================================================================

# Internal: compute auto point size from spot coordinates.
# spot_radius: the @radius value from the image's centroids slot (in full-res px).
# Calibrated so that VisiumV2 (radius ~90px, scale_constant=100) is unchanged.
# For VisiumHD (radius ~1.7px) the constant auto-scales to ~3500, matching
# the empirically validated range for 8-micron bin data. VisiumV1 (no radius)
# passes NULL and uses scale_constant as-is.
.auto_pt_size <- function(coords, scale_constant = 100, max_size = 50,
                           min_size = 0.3, spot_radius = NULL) {
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

  median_nn    <- stats::median(min_dists)
  spots_across <- equiv_dim / median_nn

  eff_const <- if (!is.null(spot_radius) && spot_radius > 0 && spot_radius < 90)
    scale_constant * (90 / spot_radius)^0.9
  else
    as.numeric(scale_constant)

  raw_size <- eff_const / spots_across
  max(min(raw_size, max_size), min_size)
}

# Internal: extract coords, cell names, and spot radius from a Seurat image slot.
# VisiumV2 and FOV both store centroids in full-res pixel coordinates.
# VisiumV1 stores spots in full-res coordinates but has no centroids slot
# (radius returned as NULL so .auto_pt_size() uses scale_constant unchanged).
.img_coords <- function(img_obj) {
  if (inherits(img_obj, c("FOV", "VisiumV2"))) {
    cen <- img_obj@boundaries[["centroids"]]
    list(coords = cen@coords, cells = cen@cells, radius = cen@radius)
  } else if (inherits(img_obj, "VisiumV1")) {
    cd <- img_obj@coordinates
    list(coords = as.matrix(cd[, c("imagerow", "imagecol")]),
         cells  = rownames(cd),
         radius = NULL)
  } else {
    stop("Unrecognized image class: ", class(img_obj))
  }
}

# Internal: categorical legend panel size (width, height in inches)
.cat_legend_dims <- function(n_lvls, display_labels, legend.direction, legendnrow) {
  max_lbl <- max(nchar(as.character(display_labels)), 1L)
  if (identical(legend.direction, "horizontal")) {
    n_rows  <- min(legendnrow, n_lvls)
    per_col <- ceiling(n_lvls / n_rows)
    c(w = per_col * (0.45 + max_lbl * 0.075) + 1.0, h = 0.5 + n_rows * 0.4)
  } else {
    c(w = (0.45 + max_lbl * 0.075) + 1.0, h = 0.5 + n_lvls * 0.35)
  }
}

# Internal: ggarrange + legend panel; sets nk_pdf_dims attr in inches
.ggarrange_with_legend <- function(inner, lgd, legend.side,
                                    grid_h, grid_w, legend_h, legend_w) {
  if (legend.side %in% c("bottom", "top")) {
    parts <- if (legend.side == "bottom") list(inner, lgd) else list(lgd, inner)
    ht    <- if (legend.side == "bottom") c(grid_h, legend_h) else c(legend_h, grid_h)
    p <- ggpubr::ggarrange(plotlist = parts, nrow = 2L, ncol = 1L, heights = ht)
    attr(p, "nk_pdf_dims") <- c(max(grid_w, legend_w), grid_h + legend_h)
  } else {
    parts <- if (legend.side == "right") list(inner, lgd) else list(lgd, inner)
    ws    <- if (legend.side == "right") c(grid_w, legend_w) else c(legend_w, grid_w)
    p <- ggpubr::ggarrange(plotlist = parts, nrow = 1L, ncol = 2L, widths = ws)
    attr(p, "nk_pdf_dims") <- c(grid_w + legend_w, max(grid_h, legend_h))
  }
  p
}

# --------------------------------------------------------------------------- #
# PlotSpatialFeaturePlots                                                      #
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
#'     metadata variable (`row.by`). The variable must have a single
#'     unique value per image.}
#' }
#'
#' @param seurat_object A Seurat object with spatial images.
#' @param features Character vector of gene names or numeric metadata
#'   column names.
#' @param layout_method Character. `"auto"` (default) or `"metadata"`.
#' @param row.by Character or `NULL`. Slide-level metadata column for
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
#' @param legend.side Character. Where to place the legend panel relative to
#'   the image grid: \code{"right"} (default), \code{"left"}, \code{"bottom"},
#'   or \code{"top"}.
#' @param legend.direction Character. Internal layout of legend entries:
#'   \code{"vertical"} (default, entries stacked) or \code{"horizontal"}
#'   (entries arranged in rows). Also controls colorbar orientation.
#' @param join_plots Logical. When `features` has more than one gene, combine
#'   all per-feature figures into a single page/PDF instead of one each.
#'   Default `FALSE`.
#' @param join_nrow,join_ncol Integer or `NULL`. Grid layout for `join_plots`.
#'   `NULL` auto-arranges into a near-square grid.
#' @param output_dir Character or `NULL`. Output directory for PDFs. Returns
#'   last plot when `NULL`.
#' @param object_name Character. Label appended to PDF filenames. Default `""`.
#' @param subset_name Character. Optional second label appended to PDF
#'   filenames after \code{object_name}. Default `""`.
#'
#' @return Invisibly returns the last plot when `output_dir = NULL`. When
#'   `join_plots = TRUE`, invisibly returns the combined tiled figure.
#' @export
PlotSpatialFeaturePlots <- function(seurat_object,
                                        features,
                                        layout_method  = "auto",
                                        row.by        = NULL,
                                        images_to_plot = NULL,
                                        remove_outliers = FALSE,
                                        outlier_prob   = 0.01,
                                        size_override  = NULL,
                                        colors         = c("#053061", "#2166AC",
                                                           "#D1E5F0", "#FDDBC7",
                                                           "#F4A582", "#D6604D",
                                                           "#B2182B", "#67001F"),
                                        calcptsizesc   = 100,
                                        legend.side    = "right",
                                        legend.direction = "vertical",
                                        join_plots     = FALSE,
                                        join_nrow      = NULL,
                                        join_ncol      = NULL,
                                        output_dir     = NULL,
                                        object_name    = "",
                                        subset_name    = "") {

  if (!layout_method %in% c("auto", "metadata"))
    stop("layout_method must be 'auto' or 'metadata'.")
  if (layout_method == "metadata" && is.null(row.by))
    stop("'row.by' must be specified for layout_method = 'metadata'.")

  # Walk up to PrepObject-stored defaults only when the caller did not pass
  # output_dir at all. An explicit output_dir = NULL means "show, don't save."
  if (missing(output_dir))
    output_dir <- if (.nk_autosave(seurat_object))
      .nk_setting(seurat_object, "output_dir") else NULL
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  .nk_warn_donor(seurat_object)
  all_imgs  <- Seurat::Images(seurat_object)
  img_names <- if (!is.null(images_to_plot)) {

    intersect(images_to_plot, all_imgs)
  } else all_imgs
  if (length(img_names) == 0) stop("No matching images found.")

  last_plot    <- NULL
  do_join      <- isTRUE(join_plots) && length(features) > 1L
  join_collect <- if (do_join) list() else NULL

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
      warning("Feature '", feat, "' not found - skipping.")
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
        vals <- unique(as.character(seurat_object@meta.data[cells, row.by]))
        if (length(vals) > 1)
          stop("Image '", img, "': '", row.by, "' has multiple values. ",
               "row.by must be a slide-level variable.")
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
        else .auto_pt_size(coords, scale_constant = calcptsizesc, spot_radius = ic$radius)
      } else {
        .auto_pt_size(coords, scale_constant = calcptsizesc, spot_radius = ic$radius)
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
    lgd_pos <- if (legend.side %in% c("bottom", "top")) "bottom" else "right"
    legend_guide <- if (is_numeric_feat) {
      ggplot2::guides(fill = ggplot2::guide_colorbar(
        title          = feat, title.position = "top",
        title.theme    = ggplot2::element_text(size = 12, face = "bold"),
        direction      = legend.direction
      ))
    } else {
      ggplot2::guides(fill = ggplot2::guide_legend(
        title          = feat, title.position = "top",
        title.theme    = ggplot2::element_text(size = 12, face = "bold"),
        override.aes   = list(size = 5),
        nrow           = if (legend.direction == "horizontal") 3L else NULL,
        ncol           = if (legend.direction == "vertical")   1L else NULL
      ))
    }
    PlotLegend <- ggpubr::get_legend(
      plot_list[[1]] + legend_guide +
        ggplot2::theme(legend.position = lgd_pos,
                       legend.direction = legend.direction)
    )

    plot_list <- lapply(plot_list, function(pp)
      pp + ggplot2::theme_void() + Seurat::NoLegend() +
        ggplot2::theme(plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 14)))

    # Assemble grid
    if (layout_method == "metadata") {
      row_plot_lists <- list()
      for (img in names(plot_list)) {
        first_cell <- .img_coords(seurat_object@images[[img]])$cells[1]
        r_val      <- as.character(seurat_object@meta.data[first_cell, row.by])
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

    grid_h_in  <- n_rows * 5
    grid_w_in  <- (n_cols * 4) + 2
    if (is_numeric_feat) {
      if (legend.direction == "vertical") {
        legend_w <- 1.5; legend_h <- min(grid_h_in * 0.7, 5.0)
      } else {
        legend_w <- max(n_cols * 2.5, 5.0); legend_h <- 1.2
      }
    } else {
      n_lvls   <- length(unique(as.character(seurat_object@meta.data[[feat]])))
      ld       <- .cat_legend_dims(n_lvls, unique(as.character(seurat_object@meta.data[[feat]])),
                                    legend.direction, legendnrow = 3L)
      legend_w <- ld[["w"]]; legend_h <- ld[["h"]]
    }

    CustomPlot <- .ggarrange_with_legend(inner_grid, PlotLegend, legend.side,
                                          grid_h_in, grid_w_in, legend_h, legend_w)
    pdf_h <- attr(CustomPlot, "nk_pdf_dims")[2]
    pdf_w <- attr(CustomPlot, "nk_pdf_dims")[1]

    if (do_join) {
      # join_plots: defer saving/printing until every feature's panel is built
      attr(CustomPlot, "nk_pdf_dims") <- c(pdf_w, pdf_h)
      join_collect[[feat]] <- CustomPlot
    } else if (!is.null(output_dir)) {
      fname <- file.path(output_dir,
                         paste0(feat, " SpatialFeatureMap ",
                                object_name,
                                if (nchar(subset_name) > 0) paste0(" ", subset_name) else "",
                                ".pdf"))
      grDevices::pdf(fname, width = pdf_w, height = pdf_h)
      print(CustomPlot)
      grDevices::dev.off()
      sp_ctx <- .nk_legend_context(seurat_object)
      sp_unit <- sp_ctx$unit   # "spots" / "bins" per assay_type
      .write_legend_sidecar(fname, paste0(
        "Spatial transcriptomics map of ", sp_ctx$n_obs, " ", sp_unit,
        if (!is.null(sp_ctx$n_donors))
          paste0(" from ", format(sp_ctx$n_donors, big.mark = ","), " tissue sections")
        else paste0(" across ", length(img_names), " tissue section(s)"),
        if (nchar(sp_ctx$obj_name) > 0) paste0(" [", sp_ctx$obj_name, "]") else "",
        ", showing ",
        if (feat %in% rownames(seurat_object)) "log-normalized expression" else "values",
        " of ", feat, ". ",
        "Each ", sp_unit, " is colored on a continuous scale from low (dark blue) ",
        "to high (dark red).",
        if (remove_outliers) paste0(
          " Outlier ", sp_unit, " above the ", round((1 - outlier_prob) * 100),
          "th percentile are excluded from the color scale to highlight ",
          "the main expression range.") else ""
      ))
    } else {
      print(CustomPlot)
      last_plot <- CustomPlot
    }
    message(sprintf("Spatial feature map: %s (%d of %d) - %dx%d",
                    feat, i, length(features), n_rows, n_cols))
  }

  # ── join_plots: tile every per-feature plot into one figure ───────────────
  if (do_join) {
    n <- length(join_collect)
    if (is.null(join_ncol) && is.null(join_nrow)) {
      ncol_j <- ceiling(sqrt(n)); nrow_j <- ceiling(n / ncol_j)
    } else if (is.null(join_ncol)) {
      nrow_j <- join_nrow; ncol_j <- ceiling(n / nrow_j)
    } else {
      ncol_j <- join_ncol; nrow_j <- ceiling(n / ncol_j)
    }

    combined <- patchwork::wrap_plots(
      lapply(join_collect, function(p) patchwork::wrap_elements(full = p)),
      nrow = nrow_j, ncol = ncol_j
    )

    one_dim <- attr(join_collect[[1]], "nk_pdf_dims") %||% c(12, 6)
    pdf_w   <- ncol_j * one_dim[1]
    pdf_h   <- nrow_j * one_dim[2]

    if (!is.null(output_dir)) {
      fname <- file.path(output_dir,
                         paste0(paste(features, collapse = "-"),
                                " SpatialFeatureMaps ", object_name,
                                if (nchar(subset_name) > 0) paste0(" ", subset_name) else "",
                                ".pdf"))
      grDevices::pdf(fname, width = pdf_w, height = pdf_h)
      print(combined)
      grDevices::dev.off()
      spj_ctx <- .nk_legend_context(seurat_object)
      spj_unit <- spj_ctx$unit
      .write_legend_sidecar(fname, paste0(
        "Spatial transcriptomics maps of ", spj_ctx$n_obs, " ", spj_unit,
        if (!is.null(spj_ctx$n_donors))
          paste0(" from ", format(spj_ctx$n_donors, big.mark = ","), " tissue sections")
        else paste0(" across ", length(img_names), " tissue section(s)"),
        if (nchar(spj_ctx$obj_name) > 0) paste0(" [", spj_ctx$obj_name, "]") else "",
        ", showing log-normalized expression of ",
        paste(features, collapse = ", "),
        ", arranged in a ", nrow_j, " x ", ncol_j, " grid. ",
        "Each panel uses its own continuous color scale (dark blue = low, ",
        "dark red = high) capped at that feature's maximum value."
      ))
      return(invisible(combined))
    }

    print(combined)
    return(invisible(combined))
  }

  invisible(last_plot)
}

# --------------------------------------------------------------------------- #
# PlotSpatialDimPlots                                                          #
# --------------------------------------------------------------------------- #

#' Batch spatial dimensionality-reduction (cell-type) plots
#'
#' Produces [Seurat::SpatialDimPlot()] panels for each image in the object,
#' one figure per metadata variable in \code{group_by_vars}. Supports the same
#' \code{"auto"} / \code{"metadata"} layout modes as
#' [PlotSpatialFeaturePlots()]. When multiple variables are supplied,
#' \code{join_plots = TRUE} tiles them into a single page.
#'
#' @param seurat_object A Seurat object with spatial images.
#' @param group_by_vars Character vector. Metadata columns to plot.
#' @param layout_method Character. \code{"auto"} or \code{"metadata"}.
#'   Default \code{"auto"}.
#' @param row.by Character or \code{NULL}. Slide-level metadata column for row
#'   grouping (required for \code{"metadata"} layout).
#' @param images_to_plot Character vector or \code{NULL}. Images to include.
#' @param remove_outliers Logical. Trim outlier cells. Default \code{FALSE}.
#' @param outlier_prob Numeric. Quantile trim level. Default \code{0.01}.
#' @param size_override Scalar or named numeric list. A plain number applies
#'   to all images; a named list applies per image.
#' @param colors Named character vector. Colors for each group level.
#'   \code{NULL} auto-resolves from PrepObject or \code{Nour_pal}.
#'   Supersedes the deprecated \code{cluster_colors}.
#' @param number_labels Logical. Prefix each legend label with a zero-padded
#'   index (e.g., "01. CellType"). Default \code{FALSE}.
#' @param legendnrow Integer or vector. Number of rows the categorical legend
#'   wraps into; PDF dimensions scale automatically to fit. When
#'   \code{group_by_vars} has multiple variables, pass a vector aligned
#'   positionally (e.g. \code{c(3, 1)}) or named by variable (e.g.
#'   \code{c(seurat_clusters = 3, Diagnosis = 1)}) to set a different row
#'   count per variable; a single value applies to all. Default \code{3}.
#' @param legend.side Character. Where to place the legend panel: \code{"bottom"}
#'   (default), \code{"top"}, \code{"left"}, or \code{"right"}.
#' @param legend.direction Character. Internal layout of legend entries:
#'   \code{"horizontal"} (default, entries in rows) or \code{"vertical"}
#'   (entries stacked in one column).
#' @param join_plots Logical. When \code{group_by_vars} has more than one
#'   variable, combine all per-variable figures into a single page instead of
#'   one page each. Default \code{FALSE}.
#' @param join_nrow,join_ncol Integer or \code{NULL}. Grid layout for
#'   \code{join_plots}. \code{NULL} auto-arranges into a near-square grid.
#' @param output_dir Character or \code{NULL}. PDF output directory. When
#'   \code{NULL} the plot is printed to the active device and returned
#'   invisibly.
#' @param rowannsize Numeric. Font size for the row-group annotation label.
#'   Default \code{16}.
#' @param imgalpha Numeric. Tissue image transparency. Default \code{1}.
#' @param alpha Numeric. Spot transparency. Default \code{1}.
#' @param uniform_size Logical. Use the median point size across all images.
#'   Default \code{FALSE}.
#' @param calcptsizesc Numeric. Scaling constant for automatic point-size
#'   calculation; larger values yield larger spots. Radius-aware, so the same
#'   value works across Visium and Visium HD. Default \code{100}.
#' @param object_name Character. Label appended to PDF filenames. Default
#'   \code{""}.
#' @param group.by Character vector. Preferred alias for \code{group_by_vars}.
#' @param cluster_colors Named character vector. Deprecated alias for
#'   \code{colors}.
#'
#' @return Invisibly returns the plot, or a named list of plots when
#'   \code{join_plots = FALSE} and multiple variables are supplied.
#' @export
PlotSpatialDimPlots <- function(seurat_object,
                                    group_by_vars   = NULL,
                                    layout_method   = "auto",
                                    row.by          = NULL,
                                    images_to_plot  = NULL,
                                    remove_outliers = FALSE,
                                    outlier_prob    = 0.01,
                                    size_override   = NULL,
                                    uniform_size    = FALSE,
                                    calcptsizesc    = 100,
                                    colors          = NULL,
                                    number_labels   = FALSE,
                                    legendnrow      = 3,
                                    legend.side     = "bottom",
                                    legend.direction = "horizontal",
                                    join_plots      = FALSE,
                                    join_nrow       = NULL,
                                    join_ncol       = NULL,
                                    output_dir      = NULL,
                                    rowannsize      = 16,
                                    imgalpha        = 1,
                                    alpha           = 1,
                                    object_name     = "",
                                    group.by        = NULL,
                                    cluster_colors  = NULL) {

  if (is.null(group_by_vars) && !is.null(group.by)) group_by_vars <- group.by
  if (!is.null(cluster_colors) && is.null(colors)) colors <- cluster_colors
  if (is.null(group_by_vars)) group_by_vars <- "seurat_clusters"

  if (!layout_method %in% c("auto", "metadata"))
    stop("layout_method must be 'auto' or 'metadata'.")
  if (layout_method == "metadata" && is.null(row.by))
    stop("'row.by' must be specified for layout_method = 'metadata'.")

  # Walk up to PrepObject-stored defaults only when the caller omitted it
  # entirely. NA_character_ is a sentinel from join recursion meaning "no
  # save"; explicit NULL means "show, don't save."
  if (identical(output_dir, NA_character_)) {
    output_dir <- NULL
  } else if (missing(output_dir)) {
    output_dir <- if (.nk_autosave(seurat_object))
      .nk_setting(seurat_object, "output_dir") else NULL
  }
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""
  .nk_warn_donor(seurat_object)

  all_imgs  <- Seurat::Images(seurat_object)
  img_names <- if (!is.null(images_to_plot)) intersect(images_to_plot, all_imgs) else all_imgs
  if (length(img_names) == 0) stop("No matching images found.")

  # ── Multiple group.by: recurse once per variable ───────────────────────────
  if (length(group_by_vars) > 1L) {
    sub_dir <- if (isTRUE(join_plots)) NA_character_ else output_dir

    plots <- lapply(seq_along(group_by_vars), function(i) {
      gb <- group_by_vars[i]
      # legendnrow may be a single number (applies to all), a vector aligned
      # positionally with group_by_vars, or a vector named by variable
      gb_legendnrow <- if (!is.null(names(legendnrow)) && gb %in% names(legendnrow)) {
        legendnrow[[gb]]
      } else {
        legendnrow[[((i - 1) %% length(legendnrow)) + 1]]
      }
      PlotSpatialDimPlots(
        seurat_object   = seurat_object,
        group_by_vars   = gb,
        layout_method   = layout_method,
        row.by          = row.by,
        images_to_plot  = images_to_plot,
        remove_outliers = remove_outliers,
        outlier_prob    = outlier_prob,
        size_override   = size_override,
        uniform_size    = uniform_size,
        colors          = colors,
        number_labels    = number_labels,
        legendnrow       = gb_legendnrow,
        legend.side      = legend.side,
        legend.direction = legend.direction,
        output_dir       = sub_dir,
        rowannsize       = rowannsize,
        imgalpha        = imgalpha,
        alpha           = alpha,
        object_name     = object_name
      )
    })
    names(plots) <- group_by_vars

    if (!isTRUE(join_plots)) return(invisible(plots))

    # ── Join: tile every per-variable plot into one figure ───────────────────
    n <- length(plots)
    if (is.null(join_ncol) && is.null(join_nrow)) {
      ncol_j <- ceiling(sqrt(n)); nrow_j <- ceiling(n / ncol_j)
    } else if (is.null(join_ncol)) {
      nrow_j <- join_nrow; ncol_j <- ceiling(n / nrow_j)
    } else {
      ncol_j <- join_ncol; nrow_j <- ceiling(n / ncol_j)
    }

    combined <- patchwork::wrap_plots(
      lapply(plots, function(p) patchwork::wrap_elements(full = p)),
      nrow = nrow_j, ncol = ncol_j
    )

    if (!is.null(output_dir)) {
      one_dim <- attr(plots[[1]], "nk_pdf_dims") %||% c(12, 8)
      pdf_w   <- ncol_j * one_dim[1]
      pdf_h   <- nrow_j * one_dim[2]
      fname   <- file.path(output_dir,
                           paste0(paste(group_by_vars, collapse = "-"),
                                  " SpatialDimMaps",
                                  if (nchar(object_name) > 0) paste0(" ", object_name) else "",
                                  ".pdf"))
      grDevices::pdf(fname, width = pdf_w, height = pdf_h)
      print(combined)
      grDevices::dev.off()
      sdj_ctx <- .nk_legend_context(seurat_object)
      sdj_unit <- sdj_ctx$unit
      .write_legend_sidecar(fname, paste0(
        "Spatial transcriptomics maps of ", sdj_ctx$n_obs, " ", sdj_unit,
        if (!is.null(sdj_ctx$n_donors))
          paste0(" from ", format(sdj_ctx$n_donors, big.mark = ","), " tissue sections")
        else paste0(" across ", length(img_names), " tissue section(s)"),
        if (nchar(sdj_ctx$obj_name) > 0) paste0(" [", sdj_ctx$obj_name, "]") else "",
        ", color-coded by ", paste(group_by_vars, collapse = " and "),
        ", arranged in a ", nrow_j, " x ", ncol_j, " grid. ",
        "Each column corresponds to one metadata variable; each row to one tissue section."
      ))
      return(invisible(combined))
    }

    print(combined)
    return(invisible(combined))
  }

  # ── Single variable ────────────────────────────────────────────────────────
  grp <- group_by_vars

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
          ic2 <- .img_coords(seurat_object@images[[img]])
          .auto_pt_size(ic2$coords, scale_constant = calcptsizesc, spot_radius = ic2$radius)
        }
      } else {
        ic2 <- .img_coords(seurat_object@images[[img]])
        .auto_pt_size(ic2$coords, scale_constant = calcptsizesc, spot_radius = ic2$radius)
      }
    }, numeric(1))
    if (uniform_size) raw_sizes[] <- stats::median(raw_sizes)
    names(raw_sizes) <- img_names

    # ── Build legend from the full color vector (not from a subset plot) ───
    # This guarantees every level appears in the legend regardless of which
    # cell types happen to be present in the first image.
    lgd_colors <- if (!is.null(grp_colors)) grp_colors else
      stats::setNames(
        Nour_pal(if (length(all_lvls) <= 8) "all" else "spectrum")(length(all_lvls)),
        all_lvls
      )
    # number_labels: show "01. Level", "02. Level" in the legend
    display_lvls <- if (number_labels) {
      stats::setNames(
        paste0(sprintf("%02d", seq_along(all_lvls)), ". ", all_lvls),
        all_lvls
      )
    } else {
      stats::setNames(all_lvls, all_lvls)
    }
    lgd_dat <- data.frame(
      x     = 1L,
      y     = seq_along(all_lvls),
      label = factor(display_lvls[all_lvls], levels = display_lvls)
    )
    # Map display labels to original colors (colors are keyed by original names)
    lgd_colors_display <- stats::setNames(lgd_colors, display_lvls[names(lgd_colors)])

    lgd_extract_pos <- if (legend.side %in% c("bottom", "top")) "bottom" else "right"
    lgd_plot <- ggplot2::ggplot(lgd_dat,
        ggplot2::aes(x = x, y = y, fill = label)) +
      ggplot2::geom_point(shape = 21, size = 5) +
      ggplot2::scale_fill_manual(values = lgd_colors_display) +
      ggplot2::guides(fill = ggplot2::guide_legend(
        title          = grp,
        override.aes   = list(size = 5),
        title.theme    = ggplot2::element_text(size = 15, face = "bold"),
        title.position = "top",
        label.theme    = ggplot2::element_text(size = 15),
        nrow           = if (legend.direction == "horizontal") legendnrow else NULL,
        ncol           = if (legend.direction == "vertical")   1L         else NULL
      )) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position  = lgd_extract_pos,
                     legend.direction = legend.direction)
    PlotLegend <- ggpubr::get_legend(lgd_plot)

    plot_list <- list()

    for (img in img_names) {
      img_obj <- seurat_object@images[[img]]
      ic      <- .img_coords(img_obj)
      coords  <- ic$coords
      cells   <- ic$cells

      if (layout_method == "metadata") {
        vals <- unique(as.character(seurat_object@meta.data[cells, row.by]))
        if (length(vals) > 1)
          stop("Image '", img, "': row.by has multiple values.")
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
        r_val      <- as.character(seurat_object@meta.data[first_cell, row.by])
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

    # ── Legend and grid dimensions ─────────────────────────────────────────────
    ld        <- .cat_legend_dims(length(all_lvls), display_lvls, legend.direction, legendnrow)
    legend_w_in <- ld[["w"]]
    legend_h_in <- ld[["h"]]
    grid_h_in   <- n_rows * 4
    grid_w_in   <- max(n_cols * 4, 12)

    CustomPlot <- .ggarrange_with_legend(inner_grid, PlotLegend, legend.side,
                                          grid_h_in, grid_w_in, legend_h_in, legend_w_in)

  pdf_w <- attr(CustomPlot, "nk_pdf_dims")[1]
  pdf_h <- attr(CustomPlot, "nk_pdf_dims")[2]

  if (!is.null(output_dir)) {
    fname <- file.path(output_dir,
                       paste0(grp, " SpatialDimMap ", object_name, ".pdf"))
    grDevices::pdf(fname, width = pdf_w, height = pdf_h)
    print(CustomPlot)
    grDevices::dev.off()
    sd_ctx <- .nk_legend_context(seurat_object)
    sd_unit <- sd_ctx$unit
    .write_legend_sidecar(fname, paste0(
      "Spatial transcriptomics map of ", sd_ctx$n_obs, " ", sd_unit,
      if (!is.null(sd_ctx$n_donors))
        paste0(" from ", format(sd_ctx$n_donors, big.mark = ","), " tissue sections")
      else paste0(" across ", length(img_names), " tissue section(s)"),
      if (nchar(sd_ctx$obj_name) > 0) paste0(" [", sd_ctx$obj_name, "]") else "",
      ", color-coded by ", grp, ". ",
      "Each ", sd_unit, " is colored by its assigned ", grp, " identity, enabling ",
      "direct comparison of transcriptional boundaries with tissue anatomy. ",
      "The shared color legend is displayed ",
      if (legend.side %in% c("bottom", "top")) legend.side else "to the right",
      " of the section grid."
    ))
    message(sprintf("Spatial DimMap: %s - %dx%d", grp, n_rows, n_cols))
    return(invisible(CustomPlot))
  }

  message(sprintf("Spatial DimMap: %s - %dx%d", grp, n_rows, n_cols))
  print(CustomPlot)
  invisible(CustomPlot)
}

# --------------------------------------------------------------------------- #
# PlotMasterMaps                                                               #
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
#' @param calcptsizesc Numeric. Scaling constant for automatic point-size
#'   calculation; larger values yield larger spots. Radius-aware, so the same
#'   value works across Visium and Visium HD. Default \code{100}.
#' @param colors Named character vector. Colors for \code{group.by} levels.
#'   \code{NULL} auto-resolves from PrepObject or \code{Nour_pal}.
#' @param gene_colors Character vector. Gradient for gene expression.
#' @param legend.side Character. Where to place the paired legends: \code{"right"}
#'   (default, stacked vertically) or \code{"bottom"} (side by side).
#' @param legend.direction Character. Internal layout of the categorical
#'   cluster legend: \code{"vertical"} (default) or \code{"horizontal"}.
#' @param output_dir Character or \code{NULL}. Output directory for PDFs.
#'   When \code{NULL} the plot is printed to the active device and returned
#'   invisibly. Walk-up from PrepObject applies when omitted entirely.
#' @param object_name Character. Label for PDF filenames. Default \code{""}.
#' @param cluster_col Deprecated alias for \code{group.by}.
#' @param cluster_colors Deprecated alias for \code{colors}.
#'
#' @return Invisibly returns the last assembled plot when \code{output_dir =
#'   NULL}; otherwise writes PDFs and returns \code{NULL} invisibly.
#' @export
PlotMasterMaps <- function(seurat_object,
                                    features,
                                    imgalpha        = 1,
                                    alpha           = 1,
                                    group.by        = "seurat_clusters",
                                    images_to_plot  = NULL,
                                    remove_outliers = FALSE,
                                    outlier_prob    = 0.01,
                                    size_override   = NULL,
                                    calcptsizesc    = 100,
                                    colors          = NULL,
                                    gene_colors     = c("#053061", "#2166AC",
                                                        "#D1E5F0", "#FDDBC7",
                                                        "#F4A582", "#D6604D",
                                                        "#B2182B", "#67001F"),
                                    legend.side      = "right",
                                    legend.direction = "vertical",
                                    output_dir      = NULL,
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
  if (missing(output_dir))
    output_dir <- if (.nk_autosave(seurat_object))
      .nk_setting(seurat_object, "output_dir") else NULL
  object_name <- if (nchar(object_name) > 0) object_name else
    .nk_setting(seurat_object, "object_name") %||% ""

  # ── Multiple group.by: recurse once per variable ───────────────────────────
  if (length(group.by) > 1L) {
    for (gb in group.by) {
      PlotMasterMaps(
        seurat_object   = seurat_object,
        features        = features,
        imgalpha        = imgalpha,
        alpha           = alpha,
        group.by        = gb,
        images_to_plot  = images_to_plot,
        remove_outliers = remove_outliers,
        outlier_prob    = outlier_prob,
        size_override   = size_override,
        colors           = colors,
        gene_colors      = gene_colors,
        legend.side      = legend.side,
        legend.direction = legend.direction,
        output_dir       = output_dir,
        object_name      = object_name
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

  last_plot <- NULL

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
        else .auto_pt_size(coords, scale_constant = calcptsizesc, spot_radius = ic$radius)
      } else {
        .auto_pt_size(coords, scale_constant = calcptsizesc, spot_radius = ic$radius)
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
        lgd_pos   <- if (legend.side %in% c("bottom", "top")) "bottom" else "right"
        dummy_cat <- Seurat::DimPlot(plot_obj, group.by = group.by) +
          ggplot2::guides(color = ggplot2::guide_legend(
            title          = group.by,
            override.aes   = list(size = 5),
            title.theme    = ggplot2::element_text(size = 14, face = "bold"),
            label.theme    = ggplot2::element_text(size = 12),
            ncol           = if (legend.direction == "vertical")   1L else NULL,
            nrow           = if (legend.direction == "horizontal") 3L else NULL)) +
          ggplot2::theme(legend.position  = lgd_pos,
                         legend.direction = legend.direction)
        if (!is.null(grp_colors))
          dummy_cat <- dummy_cat +
          ggplot2::scale_color_manual(values = grp_colors, drop = FALSE)
        ClusterLegend <- ggpubr::get_legend(dummy_cat)

        bar_dir   <- if (legend.side %in% c("bottom", "top")) "horizontal" else "vertical"
        bar_h     <- if (bar_dir == "vertical")   ggplot2::unit(3, "in") else ggplot2::unit(0.5, "in")
        bar_w     <- if (bar_dir == "horizontal") ggplot2::unit(4, "in") else ggplot2::unit(0.5, "in")
        dummy_cont <- Seurat::FeaturePlot(plot_obj, features = gene) +
          umap_gene_scale +
          ggplot2::guides(color = ggplot2::guide_colorbar(
            title          = gene,
            title.position = "top",
            title.theme    = ggplot2::element_text(size = 14, face = "bold"),
            direction      = bar_dir,
            barheight      = bar_h,
            barwidth       = bar_w)) +
          ggplot2::theme(legend.position = lgd_pos)
        GeneLegend <- ggpubr::get_legend(dummy_cont)
      }

      row_grid <- ggpubr::ggarrange(p1, p2, p3, p4, nrow = 1, ncol = 4)
      list_of_rows[[img]] <- ggpubr::annotate_figure(
        row_grid,
        left = ggpubr::text_grob(img, rot = 90, face = "bold", size = 16)
      )
    }

    n_rows    <- length(list_of_rows)
    main_grid <- ggpubr::ggarrange(plotlist = list_of_rows, nrow = n_rows, ncol = 1)

    # Combine legends and assemble with main grid based on legend.side
    if (legend.side %in% c("right", "left")) {
      legends_comb <- ggpubr::ggarrange(ClusterLegend, GeneLegend,
                                         nrow = 2, ncol = 1, heights = c(1, 1))
      parts  <- if (legend.side == "right") list(main_grid, legends_comb)
                 else list(legends_comb, main_grid)
      ws     <- if (legend.side == "right") c(4, 0.5) else c(0.5, 4)
      CustomPlot <- ggpubr::ggarrange(plotlist = parts, nrow = 1, ncol = 2, widths = ws)
      pdf_w_val  <- 18; pdf_h_val <- (n_rows * 3.5) + 1
    } else {
      legends_comb <- ggpubr::ggarrange(ClusterLegend, GeneLegend,
                                         nrow = 1, ncol = 2, widths = c(1, 0.4))
      parts  <- if (legend.side == "bottom") list(main_grid, legends_comb)
                 else list(legends_comb, main_grid)
      ht     <- if (legend.side == "bottom") c(n_rows * 3.5, 2) else c(2, n_rows * 3.5)
      CustomPlot <- ggpubr::ggarrange(plotlist = parts, nrow = 2, ncol = 1, heights = ht)
      pdf_w_val  <- 18; pdf_h_val <- (n_rows * 3.5) + 3
    }

    final_plot   <- ggpubr::annotate_figure(
      CustomPlot,
      top = ggpubr::text_grob(paste0("Gene Master Map: ", gene),
                              face = "bold", size = 22, color = "#B2182B")
    )

    if (!is.null(output_dir)) {
      fname <- file.path(output_dir,
                         paste0(gene, "_", group.by,
                                " Master Gene Map ", object_name, ".pdf"))
      grDevices::pdf(fname, width = pdf_w_val, height = pdf_h_val)
      print(final_plot)
      grDevices::dev.off()
      mg_ctx  <- .nk_legend_context(seurat_object)
      mg_unit <- mg_ctx$unit
      .write_legend_sidecar(fname, paste0(
        "Master gene map for ", gene, " showing co-localization of expression ",
        "with ", group.by, " identity. ",
        "Data: ", mg_ctx$n_obs, " ", mg_unit,
        if (!is.null(mg_ctx$n_donors))
          paste0(" from ", format(mg_ctx$n_donors, big.mark = ","), " tissue sections")
        else paste0(" across ", length(img_names), " tissue section(s)"),
        if (nchar(mg_ctx$obj_name) > 0) paste0(" [", mg_ctx$obj_name, "]") else "",
        ". Each row corresponds to one tissue section and contains four panels: ",
        "(1) UMAP colored by ", group.by, " clusters, ",
        "(2) UMAP colored by ", gene, " log-normalized expression, ",
        "(3) spatial map of ", mg_unit, " colored by ", group.by, " clusters, and ",
        "(4) spatial map of ", mg_unit, " colored by ", gene, " expression ",
        "(dark blue = low, dark red = high). ",
        "Paired categorical and continuous legends are shown ",
        if (legend.side %in% c("bottom", "top")) legend.side else "to the right", "."
      ))
    } else {
      print(final_plot)
      last_plot <- final_plot
    }
    message(sprintf("Master gene map: %s / %s (%d of %d)",
                    gene, group.by, gene_idx, length(features)))
  }
  invisible(last_plot)
}

# =============================================================================
# Long-form aliases  (kept for backward compatibility)
# =============================================================================

#' @describeIn PlotSpatialFeaturePlots Long-form alias.
#' @export
GenerateSpatialFeatureMaps <- PlotSpatialFeaturePlots

#' @describeIn PlotSpatialDimPlots Long-form alias.
#' @export
GenerateSpatialDimMaps <- PlotSpatialDimPlots

#' @describeIn PlotMasterMaps Long-form alias.
#' @export
GenerateMasterGeneMaps <- PlotMasterMaps
