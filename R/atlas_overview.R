# =============================================================================
# scSidekick - Atlas Overview Wheel  (atlas_overview.R)
#
# Exported:
#   PlotAtlasWheel()  - circular "dataset map" showing groups as polar wedges,
#                       annotated with study names and summary statistics.
#                       Works with Seurat objects or plain data.frames.
#                       Supports facet.by to build one wheel per group level.
# =============================================================================


#' Dataset Overview Wheel
#'
#' @description
#' Creates a circular "atlas wheel" that summarises a large multi-cohort
#' dataset.  Each wedge represents one level of \code{group.by}
#' (e.g. cancer type, tissue, condition) and is annotated with the names of
#' contributing studies.  Summary statistics (total cells, donors, groups,
#' and optionally normal-tissue types) are placed in the center.
#'
#' Works with both **Seurat objects** and plain \code{data.frame}s.
#'
#' When \code{facet.by} is supplied, one wheel is built per level, all
#' assembled side-by-side with \pkg{patchwork}.  The filename includes
#' \code{"by_{facet.by}"}.  Outer label margins are auto-expanded based on the
#' longest group-level name.
#'
#' @param meta A Seurat object **or** a plain \code{data.frame}/\code{tibble}
#'   of per-cell (or per-sample) metadata.
#' @param group.by Character.  Primary grouping column - one wedge per level
#'   (e.g. \code{"CancerType"}).
#' @param study_by Character or \code{NULL}.  Column whose values are pasted
#'   inside each wedge as a multi-line study label (e.g. \code{"Study"}).
#'   When \code{NULL}, no inner text is drawn.
#' @param patient_by Character or \code{NULL}.  Column identifying individual
#'   donors / patients - used to count unique donors in the summary statistics.
#' @param facet.by Character or \code{NULL}.  Metadata column to facet by -
#'   one wheel per level (e.g. \code{"Tissue"}).  Default \code{NULL}.
#' @param ncol_facet Integer or \code{NULL}.  Number of wheel columns when
#'   \code{facet.by} is used.  Auto-computed (max 2) if \code{NULL}.
#' @param exclude_groups Character vector or \code{NULL}.  Levels of
#'   \code{group.by} to exclude from the wheel itself, but whose cells are
#'   still counted in the center statistics.  Typical use:
#'   \code{exclude_groups = "Healthy"}.
#' @param normal_col,normal_val Character or \code{NULL}.  Column and value
#'   used to identify normal (non-disease) samples
#'   (e.g. \code{normal_col = "Health"}, \code{normal_val = "Healthy"}).
#' @param normal_count_by Character or \code{NULL}.  Column to count distinct
#'   normal tissue types (e.g. \code{"Organ"}).
#' @param title Character.  Main title text in the center of the wheel.
#'   Supports \code{"\\n"} for line breaks.
#' @param cell_label Character.  Descriptor appended to the cell count
#'   (e.g. \code{"Myeloid Cells"}).  Default \code{"Cells"}.
#' @param show_stats Logical.  Show cell / donor / group counts in the center.
#'   Default \code{TRUE}.
#' @param colors Named character vector mapping \code{group.by} levels to
#'   colors.  Auto-assigned if \code{NULL}.
#' @param show_circles Logical.  Draw a white placeholder circle at the inner
#'   base of each wedge.  Default \code{TRUE}.
#' @param circle_size Numeric.  \code{size} of placeholder circles.
#'   Default \code{18}.
#' @param wedge_alpha Numeric.  Transparency of the wedge fill.  Default
#'   \code{0.5}.
#' @param y_inner Numeric.  Inner y edge of the wedge ring (top of the center
#'   hole).  Default \code{-1}.
#' @param y_outer Numeric.  Outer y edge of the wedge ring.  Default \code{10}.
#' @param y_text Numeric.  y position of the study-name text inside each
#'   wedge.  Default \code{7}.
#' @param y_label Numeric.  y position of the outer group-name labels.
#'   Default \code{11}.
#' @param y_circle Numeric.  y position of the inner placeholder circles.
#'   Default \code{1.5}.
#' @param y_min Numeric.  Bottom of the y-axis range.  Default \code{-15}.
#' @param y_center_title Numeric.  y position of the center title annotation.
#'   Default \code{-5.5}.
#' @param y_center_stats Numeric.  y position of the center statistics.
#'   Default \code{-13}.
#' @param label_size Numeric.  Font size of outer group-name labels.
#'   Default \code{4}.
#' @param study_text_size Numeric.  Font size of study-name text inside
#'   wedges.  Default \code{2.2}.
#' @param study_lineheight Numeric.  Line height for multi-line study labels.
#'   Default \code{0.9}.
#' @param title_size Numeric.  Font size of the center title.  Default
#'   \code{7}.
#' @param stats_size Numeric.  Font size of the center statistics text.
#'   Default \code{5}.
#' @param title_color Character.  Color for the center title and circle
#'   borders.  Default \code{"#1e4e5f"}.
#' @param output_dir Directory to save a PDF.  \code{NULL} = no save.
#' @param object_name Label prefix for the output filename.
#' @param pdf_width,pdf_height PDF dimensions per wheel in inches.  Default
#'   \code{9 x 9}.
#'
#' @return A \code{ggplot2} or \code{patchwork} object (invisibly).
#' @export
PlotAtlasWheel <- function(
    meta,
    group.by,
    study_by          = NULL,
    patient_by        = NULL,
    facet.by          = NULL,
    ncol_facet        = NULL,
    exclude_groups    = NULL,
    normal_col        = NULL,
    normal_val        = NULL,
    normal_count_by   = NULL,
    title             = "Dataset Overview",
    cell_label        = "Cells",
    show_stats        = TRUE,
    colors            = NULL,
    show_circles      = TRUE,
    circle_size       = 18,
    wedge_alpha       = 0.5,
    y_inner           = -1,
    y_outer           = 10,
    y_text            = 7,
    y_label           = 11,
    y_circle          = 1.5,
    y_min             = -15,
    y_center_title    = -5.5,
    y_center_stats    = -13,
    label_size        = 4,
    study_text_size   = 2.2,
    study_lineheight  = 0.9,
    title_size        = 7,
    stats_size        = 5,
    title_color       = "#1e4e5f",
    output_dir        = NULL,
    object_name       = "",
    pdf_width         = 9,
    pdf_height        = 9
) {

  # ── Extract metadata ──────────────────────────────────────────────────────
  md <- if (inherits(meta, "Seurat")) meta@meta.data else as.data.frame(meta)

  for (col in c(group.by, study_by, patient_by, facet.by,
                normal_col, normal_count_by)) {
    if (!is.null(col) && !col %in% colnames(md))
      stop("'", col, "' not found in metadata.")
  }

  md <- md[!is.na(md[[group.by]]), , drop = FALSE]
  if (!is.null(facet.by))
    md <- md[!is.na(md[[facet.by]]), , drop = FALSE]

  # ── All group levels (for color assignment - consistent across facets) ────
  all_group_levels <- if (is.factor(md[[group.by]])) levels(md[[group.by]])
                      else sort(unique(as.character(md[[group.by]])))

  # ── Colors ────────────────────────────────────────────────────────────────
  if (is.null(colors)) {
    colors <- if (inherits(meta, "Seurat")) .nk_colors(meta, group.by) else NULL
    if (is.null(colors)) {
      n_col  <- length(all_group_levels)
      colors <- stats::setNames(
        Nour_pal(if (n_col <= 8) "all" else "spectrum")(n_col),
        all_group_levels)
    }
  }

  # ── Auto-pad: outer labels can extend beyond polar boundary ───────────────
  # ~0.35 mm per character per size unit; add buffer for the y_label offset
  outer_pad <- max(5, max(nchar(all_group_levels)) * label_size * 0.35)

  # ── Inner wheel builder (one metadata slice → one ggplot) ─────────────────
  .build_one_wheel <- function(md_sub, wheel_title = title) {

    # Global stats from this slice (before exclude_groups filter)
    total_cells  <- nrow(md_sub)
    total_donors <- if (!is.null(patient_by))
      length(unique(md_sub[[patient_by]])) else NA_integer_

    n_normal_types <- NA_integer_
    if (!is.null(normal_col) && !is.null(normal_val)) {
      md_norm <- md_sub[as.character(md_sub[[normal_col]]) %in% normal_val, , drop = FALSE]
      n_normal_types <- if (!is.null(normal_count_by))
        length(unique(md_norm[[normal_count_by]])) else nrow(md_norm)
    }

    # Filter wheel data
    md_wheel <- if (!is.null(exclude_groups))
      md_sub[!as.character(md_sub[[group.by]]) %in% exclude_groups, , drop = FALSE]
    else md_sub

    wheel_grp_char   <- as.character(md_wheel[[group.by]])
    wheel_grp_levels <- if (is.factor(md_wheel[[group.by]])) {
      intersect(levels(md_wheel[[group.by]]), unique(wheel_grp_char))
    } else sort(unique(wheel_grp_char))

    n_groups <- length(wheel_grp_levels)
    if (n_groups == 0L) return(NULL)

    # Aggregate to per-wedge data frame
    if (!is.null(study_by)) {
      agg <- md_wheel |>
        dplyr::group_by(dplyr::across(dplyr::all_of(c(study_by, group.by)))) |>
        dplyr::summarize(
          n_cells    = dplyr::n(),
          n_patients = if (!is.null(patient_by))
            dplyr::n_distinct(.data[[patient_by]]) else NA_integer_,
          .groups    = "drop"
        )
      df <- agg |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group.by))) |>
        dplyr::summarize(
          n_cells    = sum(n_cells),
          n_patients = if (!is.null(patient_by)) sum(n_patients) else NA_integer_,
          Studies    = paste(unique(.data[[study_by]]), collapse = "\n"),
          .groups    = "drop"
        )
    } else {
      df <- md_wheel |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group.by))) |>
        dplyr::summarize(
          n_cells    = dplyr::n(),
          n_patients = if (!is.null(patient_by))
            dplyr::n_distinct(.data[[patient_by]]) else NA_integer_,
          Studies    = "",
          .groups    = "drop"
        )
    }

    df[[group.by]] <- factor(as.character(df[[group.by]]), levels = wheel_grp_levels)
    df <- df[order(df[[group.by]]), ]
    df$id <- seq_len(nrow(df))

    # Angle calculations for rotated outer labels
    n_bars    <- nrow(df)
    raw_angle <- 90 - 360 * (df$id - 0.5) / n_bars
    df$angle  <- ifelse(raw_angle < -90, raw_angle + 180, raw_angle)
    df$hjust  <- ifelse(raw_angle < -90, 1, 0)

    # Center statistics text
    stats_lines <- c(
      paste0(format(total_cells, big.mark = ","), " ", cell_label),
      if (!is.na(total_donors))
        paste0(format(total_donors, big.mark = ","), " Donors"),
      paste0(n_groups, " Groups"),
      if (!is.na(n_normal_types) && n_normal_types > 0L)
        paste0(n_normal_types, " Normal tissue types")
    )
    stats_text <- paste(stats_lines, collapse = "\n\n")

    # ── Build ggplot ──────────────────────────────────────────────────────
    p <- ggplot2::ggplot(df) +

      # Colored wedges
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = id - 0.5, xmax = id + 0.5,
          ymin = y_inner,  ymax = y_outer,
          fill = .data[[group.by]]
        ),
        color = "white", alpha = wedge_alpha
      ) +
      ggplot2::scale_fill_manual(values = colors) +

      # Study names inside wedges
      ggplot2::geom_text(
        ggplot2::aes(
          x     = id,
          y     = y_text,
          label = Studies,
          angle = angle
        ),
        fontface   = "bold",
        size       = study_text_size,
        lineheight = study_lineheight,
        color     = "black"
      ) +

      # Outer group-name labels
      ggplot2::geom_text(
        ggplot2::aes(
          x     = id,
          y     = y_label,
          label = .data[[group.by]],
          angle = angle,
          hjust = hjust
        ),
        size     = label_size,
        fontface = "bold"
      ) +

      ggplot2::scale_x_continuous(
        limits = c(0.5, n_bars + 0.5),
        expand = c(0, 0)
      ) +
      ggplot2::ylim(y_min, y_label + 5) +
      ggplot2::coord_polar(theta = "x", clip = "off") +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "none",
        plot.margin     = ggplot2::unit(rep(outer_pad, 4), "mm")
      )

    # Inner placeholder circles
    if (isTRUE(show_circles))
      p <- p + ggplot2::geom_point(
        ggplot2::aes(x = id, y = y_circle),
        shape  = 21, fill = "white",
        color  = title_color, size = circle_size, stroke = 1
      )

    # Center title
    p <- p + ggplot2::annotate(
      "text",
      x          = 0.5, y = y_center_title,
      label      = wheel_title,
      fontface   = "bold",
      size       = title_size,
      color      = title_color,
      lineheight = 1
    )

    # Center statistics
    if (isTRUE(show_stats))
      p <- p + ggplot2::annotate(
        "text",
        x          = 0.5, y = y_center_stats,
        label      = stats_text,
        size       = stats_size,
        color      = "black",
        lineheight = 1
      )

    p
  }

  # ── Build plot(s) ─────────────────────────────────────────────────────────
  if (!is.null(facet.by)) {
    facet_lvls <- if (is.factor(md[[facet.by]])) levels(md[[facet.by]])
                  else sort(unique(as.character(md[[facet.by]])))

    plot_list <- Filter(Negate(is.null), lapply(facet_lvls, function(fv)
      .build_one_wheel(
        md_sub      = md[as.character(md[[facet.by]]) == fv, , drop = FALSE],
        wheel_title = paste0(title, "\n(", fv, ")")
      )
    ))

    ncol_f   <- ncol_facet %||% min(length(plot_list), 2L)
    nrow_f   <- ceiling(length(plot_list) / ncol_f)
    combined <- patchwork::wrap_plots(plot_list, ncol = ncol_f)
    pdf_w    <- pdf_width  * ncol_f
    pdf_h    <- pdf_height * nrow_f
  } else {
    combined <- .build_one_wheel(md)
    pdf_w    <- pdf_width
    pdf_h    <- pdf_height
  }

  # ── Filename ──────────────────────────────────────────────────────────────
  parts <- c(
    if (nchar(object_name) > 0) object_name,
    group.by,
    if (!is.null(facet.by)) paste0("by_", facet.by),
    "AtlasWheel"
  )
  fname <- paste(parts, collapse = "_")

  # ── Save + .legend sidecar ────────────────────────────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))
    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    print(combined)
    grDevices::dev.off()
    message("scSidekick: Saved ", basename(fpath),
            " (", round(pdf_w, 1), " x ", round(pdf_h, 1), " in)")

    .write_legend_sidecar(fpath, paste0(
      "Circular dataset overview wheel showing ",
      group.by, " composition",
      if (!is.null(study_by)) paste0(", with contributing studies (", study_by, ") labeled inside each wedge") else "",
      ". Each wedge represents one ", group.by, " level; wedge width is uniform (does not encode cell count).",
      if (show_stats) paste0(" Center shows total cell count, donor count, and group count.") else "",
      if (!is.null(patient_by)) paste0(" Donor counts are based on unique ", patient_by, " values.") else "",
      if (!is.null(facet.by))
        paste0(" Separate wheels are shown for each level of ", facet.by, ".")
      else "",
      if (nchar(object_name) > 0) paste0(" Dataset: ", object_name, ".") else ""
    ))
  }

  print(combined)
  invisible(combined)
}
