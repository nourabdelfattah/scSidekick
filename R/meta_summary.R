# =============================================================================
# meta_summary.R
#
# SummarizeMetadata() - collapse cell-level metadata to donor / patient level
# PlotMetaSummary()   - faceted bar chart of metadata distributions
# =============================================================================


# =============================================================================
# SummarizeMetadata
# =============================================================================

#' Summarize cell-level metadata to donor / patient level
#'
#' Groups cells by one or more ID columns (e.g. donor ID, sample ID) and
#' returns a wide data frame with one row per unique ID combination.
#'
#' \strong{Column handling:}
#' \describe{
#'   \item{Numeric columns}{Summarized with \code{numeric_func} (default:
#'     \code{mean}).}
#'   \item{Categorical columns - simple}{When every cell in an ID group shares
#'     the same value (e.g. every cell from donor D1 has \code{Sex = "Male"}),
#'     that single value is returned as-is.}
#'   \item{Categorical columns - complex}{When cells in an ID group have
#'     multiple values (e.g. a donor has several cell types), the result is a
#'     formatted string listing all values and their cell counts in descending
#'     order: \code{"CD4 T (n=1 420), NK (n=380), B (n=200)"}.}
#' }
#'
#' @param data A Seurat object or a plain data frame / tibble.
#' @param id_columns Character vector of one or more column names that
#'   uniquely identify a donor or sample (e.g.
#'   \code{c("Donor.ID")} or \code{c("Donor.ID", "Visit")}).
#'   One output row is produced per unique combination.
#' @param numeric_func Function to summarize numeric columns.  Must accept a
#'   numeric vector and return a single value.  Default \code{mean}.  Common
#'   alternatives: \code{median}, \code{max}, \code{function(x) sd(x, na.rm = TRUE)}.
#' @param numeric_na_rm Logical.  Pass \code{na.rm = TRUE} to
#'   \code{numeric_func}.  Default \code{TRUE}.
#' @param pivot Logical.  If \code{TRUE}, the result is pivoted to long format
#'   (one row per ID × variable combination).  Default \code{FALSE} (wide).
#'
#' @return A data frame (wide by default).  Always includes \code{n_cells}
#'   (number of cells per ID combination).
#' @export
#'
#' @examples
#' \dontrun{
#' # One ID column
#' donor_df <- SummarizeMetadata(SeuratObj, id_columns = "Donor.ID")
#'
#' # Multiple ID columns (one row per donor × visit combination)
#' donor_df <- SummarizeMetadata(SeuratObj,
#'   id_columns   = c("Donor.ID", "Visit"),
#'   numeric_func = median)
#'
#' # Long format for easy downstream plotting
#' long_df <- SummarizeMetadata(SeuratObj, id_columns = "Donor.ID", pivot = TRUE)
#' }
SummarizeMetadata <- function(data,
                               id_columns,
                               numeric_func  = mean,
                               numeric_na_rm = TRUE,
                               pivot         = FALSE) {

  # ── 1. Extract metadata ───────────────────────────────────────────────────
  meta <- if (inherits(data, "Seurat")) data@meta.data else as.data.frame(data)

  # ── 2. Validate id_columns ────────────────────────────────────────────────
  missing_ids <- setdiff(id_columns, colnames(meta))
  if (length(missing_ids) > 0L)
    stop("The following id_columns were not found: ",
         paste(missing_ids, collapse = ", "))

  # ── 3. Identify column types ──────────────────────────────────────────────
  other_cols <- setdiff(colnames(meta), id_columns)
  num_cols   <- other_cols[vapply(meta[other_cols], is.numeric, logical(1))]
  cat_cols   <- setdiff(other_cols, num_cols)

  func_name  <- tryCatch(deparse(substitute(numeric_func)),
                         error = function(e) "custom function")
  message("scSidekick SummarizeMetadata: ",
          format(nrow(meta), big.mark = ","), " cells × ",
          length(id_columns), " ID column(s). ",
          length(num_cols), " numeric [", func_name, "], ",
          length(cat_cols), " categorical.")

  # ── 4. Base: n_cells per ID combination ──────────────────────────────────
  result <- meta |>
    dplyr::group_by(dplyr::across(dplyr::all_of(id_columns))) |>
    dplyr::summarize(n_cells = dplyr::n(), .groups = "drop")

  # ── 5. Numeric columns ────────────────────────────────────────────────────
  if (length(num_cols) > 0L) {
    num_summary <- meta |>
      dplyr::group_by(dplyr::across(dplyr::all_of(id_columns))) |>
      dplyr::summarize(
        dplyr::across(
          dplyr::all_of(num_cols),
          \(v) numeric_func(v, na.rm = numeric_na_rm)
        ),
        .groups = "drop"
      )
    result <- dplyr::left_join(result, num_summary, by = id_columns)
  }

  # ── 6. Categorical columns ────────────────────────────────────────────────
  for (col in cat_cols) {
    # Count cells per (id_combination × column_value)
    val_counts <- meta |>
      dplyr::count(dplyr::across(dplyr::all_of(c(id_columns, col))),
                   name = ".n_cells_inner")

    col_summary <- val_counts |>
      dplyr::group_by(dplyr::across(dplyr::all_of(id_columns))) |>
      dplyr::summarize(
        !!col := {
          vals   <- as.character(.data[[col]])
          counts <- .data$.n_cells_inner
          if (dplyr::n() == 1L) {
            # Simple: every cell in this ID group has the same value
            vals[1L]
          } else {
            # Complex: multiple values - list all with counts, sorted descending
            ord <- order(counts, decreasing = TRUE)
            paste(paste0(vals[ord], " (n=", counts[ord], ")"), collapse = ", ")
          }
        },
        .groups = "drop"
      )

    result <- dplyr::left_join(result, col_summary, by = id_columns)
  }

  # ── 7. Optional pivot to long ─────────────────────────────────────────────
  if (pivot) {
    result <- tidyr::pivot_longer(
      result,
      cols      = -dplyr::all_of(c(id_columns, "n_cells")),
      names_to  = "variable",
      values_to = "value"
    )
  }

  result
}


# =============================================================================
# PlotMetaSummary
# =============================================================================

#' Faceted bar chart summarizing metadata distributions across patients
#'
#' Deduplicates cell-level metadata to donor / patient level and plots a
#' multi-panel bar chart: one column per metadata variable, one row per
#' level of an optional row-splitting variable (e.g. Sex).  Bars are stacked
#' by \code{fill_variable} and each segment is optionally annotated with its
#' count.
#'
#' @param data A Seurat object or a plain data frame.
#' @param id_column Character.  The column that uniquely identifies a donor or
#'   patient (e.g. \code{"Donor.ID"}).  Set to \code{NULL} to count cells
#'   directly without deduplication.
#' @param variables Character vector.  Metadata columns to display as
#'   \strong{column facets} on the x-axis (e.g.
#'   \code{c("Braak", "APOE.Genotype", "CERAD.score")}).
#' @param fill_variable Character.  Metadata column used to \strong{color-fill}
#'   the stacked bars (e.g. \code{"Dementia.AD"}).  When \code{NULL}, bars are
#'   shown in a single neutral color with no legend.
#' @param row_variable Character.  Metadata column used to create
#'   \strong{row facets} (e.g. \code{"Sex"}).  \code{NULL} = no row splitting.
#' @param exclude Named list of values to \strong{exclude} before plotting.
#'   Each name is a column name and each value is a character vector of levels
#'   to drop.  Example: \code{list(Cognitive.Status = "Reference")}.
#' @param percent Logical.  If \code{TRUE}, convert counts to percentages
#'   before plotting.  When \code{fill_variable} is set, each bar sums to
#'   100\% (fill segments show their share).  When \code{fill_variable} is
#'   \code{NULL}, bars show the percentage of the panel total.  Segment labels
#'   automatically switch to \code{"X\%"} format.  Default \code{FALSE}.
#' @param count_unit One of \code{"auto"} (default), \code{"donors"},
#'   or \code{"cells"}.  \code{"auto"} counts donors when \code{id_column} is
#'   provided and more than one unique donor exists per group; otherwise counts
#'   cells.
#' @param show_counts Logical.  Annotate each bar segment with its count.
#'   Default \code{TRUE}.
#' @param min_count_to_label Integer.  Suppress the count label on segments
#'   with fewer than this many observations.  Default \code{1}.
#' @param label_size Numeric.  Text size for count labels (passed to
#'   \code{geom_text}).  Default \code{3}.
#' @param colors Named character vector of colors for \code{fill_variable}
#'   levels.  \code{NULL} auto-resolves from \code{PrepObject} (if a Seurat
#'   object is passed) or \code{SelectColors()}.
#' @param y_label Character or \code{NULL}.  Override the y-axis label.
#'   \code{NULL} derives it from \code{count_unit}.
#' @param output_dir Character or \code{NULL}.  Directory to save a PDF.
#'   \code{NULL} returns the plot without saving; for a Seurat object it first
#'   walks up to the \code{output_dir} stored by \code{\link{PrepObject}}
#'   (unless \code{AutoSavePlots = FALSE}).
#' @param object_name Character.  Prefix for the auto-generated file name.
#'   Falls back to the \code{object_name} stored by \code{\link{PrepObject}}.
#' @param file_name Character or \code{NULL}.  Base name (no extension) for the
#'   saved PDF.  \code{NULL} (default) auto-deduces from \code{object_name},
#'   \code{fill_variable}, and \code{variables}.
#' @param return_data Logical.  Also return the deduplicated donor-level data
#'   frame.  Default \code{FALSE}.
#' @param column_variable Character or \code{NULL}.  Metadata column used to
#'   create \strong{column facets within each panel}.  When supplied, one plot
#'   is produced per entry in \code{variables} (multi-page PDF), and each plot
#'   is faceted by this variable.  \code{NULL} = no column faceting (default).
#' @param alpha Numeric.  Fill transparency for bars (0-1).  Default
#'   \code{0.85}.
#' @param return_flextable Logical.  Also return a list of formatted
#'   \code{flextable} objects: \code{$donor_table} (wide donor-level data) and
#'   \code{$crosstab} (cross-tabulation per variable, rows separated by
#'   variable).  Requires the \pkg{flextable} package.  Default \code{FALSE}.
#'
#' @return When only the plot is requested (default), a \code{ggplot2} object.
#'   When \code{return_data} or \code{return_flextable} are \code{TRUE}, a
#'   named list with elements \code{$plot}, optionally \code{$data}, and
#'   optionally \code{$flextable}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' PlotMetaSummary(SeuratObj,
#'   id_column     = "Donor.ID",
#'   variables     = c("Braak", "APOE.Genotype", "CERAD.score"),
#'   fill_variable = "Dementia.AD",
#'   row_variable  = "Sex",
#'   exclude       = list(Cognitive.Status = "Reference"))
#'
#' # Also return the donor table and a cross-tab flextable
#' out <- PlotMetaSummary(SeuratObj,
#'   id_column        = "Donor.ID",
#'   variables        = c("Braak", "CERAD.score"),
#'   fill_variable    = "Dementia.AD",
#'   return_data      = TRUE,
#'   return_flextable = TRUE)
#' out$plot
#' out$data
#' out$flextable$crosstab
#' }
PlotMetaSummary <- function(data,
                             id_column          = NULL,
                             variables,
                             fill_variable      = NULL,
                             row_variable       = NULL,
                             column_variable    = NULL,
                             exclude            = NULL,
                             percent            = FALSE,
                             count_unit         = c("auto", "donors", "cells"),
                             show_counts        = TRUE,
                             min_count_to_label = 1L,
                             label_size         = 3,
                             alpha              = 0.85,
                             colors             = NULL,
                             y_label            = NULL,
                             output_dir         = NULL,
                             object_name        = "",
                             file_name          = NULL,
                             return_data        = FALSE,
                             return_flextable   = FALSE) {

  count_unit <- match.arg(count_unit)

  # ── 0. Walk up PrepObject defaults (Seurat only) ──────────────────────────
  if (inherits(data, "Seurat")) {
    output_dir  <- output_dir %||%
      if (.nk_autosave(data)) .nk_setting(data, "output_dir") else NULL
    object_name <- if (nchar(object_name) > 0) object_name else
      .nk_setting(data, "object_name") %||% ""
  }

  # ── 1. Extract metadata ───────────────────────────────────────────────────
  meta <- if (inherits(data, "Seurat")) data@meta.data else as.data.frame(data)

  # ── 2. Validate columns ───────────────────────────────────────────────────
  all_needed  <- unique(c(id_column, variables, fill_variable, row_variable, column_variable))
  missing_col <- setdiff(all_needed, colnames(meta))
  if (length(missing_col) > 0L)
    stop("Column(s) not found in the data: ",
         paste(missing_col, collapse = ", "))

  # ── 2b. Split numeric vs categorical variables ───────────────────────────
  numeric_vars <- variables[vapply(variables,
    function(v) v %in% colnames(meta) && is.numeric(meta[[v]]), logical(1))]
  cat_vars <- setdiff(variables, numeric_vars)

  numeric_plot <- NULL
  if (length(numeric_vars) > 0L) {
    feat_group <- fill_variable %||% id_column
    if (is.null(feat_group))
      stop("Cannot redirect numeric variables to PlotFeature(): ",
           "supply fill_variable or id_column so PlotFeature knows how to group.")
    message("scSidekick: ", paste(numeric_vars, collapse = ", "),
            " is numeric - passing to PlotFeature(group_by = \"",
            feat_group, "\") automatically.")
    numeric_plot <- tryCatch(
      PlotFeature(
        data        = data,
        features    = numeric_vars,
        group_by    = feat_group,
        split_by    = row_variable,
        colors      = colors,
        output_dir  = output_dir,
        object_name = object_name
      ),
      error = function(e) {
        warning("PlotFeature() failed for numeric variables: ",
                conditionMessage(e))
        NULL
      }
    )
  }

  # If ALL variables were numeric, return early
  if (length(cat_vars) == 0L) {
    if (!is.null(numeric_plot)) return(numeric_plot)
    stop("All variables are numeric. Use PlotFeature() directly.")
  }
  variables <- cat_vars   # proceed with categorical variables only

  # ── 3. Apply exclusions ───────────────────────────────────────────────────
  if (!is.null(exclude)) {
    for (col in names(exclude)) {
      if (!col %in% colnames(meta)) {
        warning("Exclusion column '", col, "' not found - skipping.")
        next
      }
      meta <- meta[!as.character(meta[[col]]) %in% as.character(exclude[[col]]),
                   , drop = FALSE]
    }
    message("scSidekick PlotMetaSummary: ",
            format(nrow(meta), big.mark = ","),
            " rows remain after exclusions.")
  }

  # ── 4. Select only the needed columns ────────────────────────────────────
  keep_cols <- intersect(all_needed, colnames(meta))
  meta      <- meta[, keep_cols, drop = FALSE]

  # ── 5. Determine count_unit ───────────────────────────────────────────────
  if (count_unit == "auto") {
    if (!is.null(id_column)) {
      check_cols <- intersect(c(variables, fill_variable, row_variable),
                              colnames(meta))
      max_donors <- meta |>
        dplyr::group_by(dplyr::across(dplyr::all_of(check_cols))) |>
        dplyr::summarize(
          .n = dplyr::n_distinct(.data[[id_column]]),
          .groups = "drop"
        ) |>
        dplyr::pull(.n) |>
        max(na.rm = TRUE)

      if (max_donors <= 1L) {
        message("scSidekick: Each group has only one unique '", id_column,
                "' - counting cells instead of donors. ",
                "Pass count_unit = \"donors\" to override.")
        count_unit <- "cells"
      } else {
        count_unit <- "donors"
      }
    } else {
      count_unit <- "cells"
    }
  }

  # ── 6. Deduplicate to donor level or keep all cells ───────────────────────
  # When id_column + count_unit = "cells", keep all cells but aggregate to
  # per-donor means with SE error bars in step 8 (use_donor_avg path).
  use_donor_avg <- count_unit == "cells" && !is.null(id_column)

  if (count_unit == "donors" && !is.null(id_column)) {
    meta_dedup    <- dplyr::distinct(meta)
    auto_y_label  <- paste0("Number of ", id_column, "s")
  } else {
    meta_dedup    <- meta
    auto_y_label  <- if (use_donor_avg)
      paste0("Mean cells per ", id_column)
    else
      "Number of cells"
  }
  # User-supplied y_label overrides everything; otherwise respect percent mode
  y_label <- y_label %||%
    if (percent)
      if (!is.null(fill_variable)) "Percentage (%)" else "Percent of total (%)"
    else
      auto_y_label

  # ── 7. Pivot to long format on variables ─────────────────────────────────
  # Collect the per-variable level order BEFORE pivoting so the x-axis
  # respects existing factor levels (or alphabetical order for non-factors).
  val_levels <- unlist(lapply(variables, function(v) {
    col <- meta_dedup[[v]]
    if (is.factor(col)) levels(col) else sort(unique(as.character(col)))
  }), use.names = FALSE)
  val_levels <- unique(val_levels)   # deduplicate while preserving order

  long_df <- tidyr::pivot_longer(
    meta_dedup,
    cols      = dplyr::all_of(variables),
    names_to  = "Variable",
    values_to = "Value"
  )
  long_df$Variable <- factor(long_df$Variable, levels = variables)
  long_df$Value    <- factor(as.character(long_df$Value), levels = val_levels)

  # ── 8. Aggregate counts per segment ──────────────────────────────────────
  group_cols <- c("Variable", "Value")
  if (!is.null(fill_variable))   group_cols <- c(group_cols, fill_variable)
  if (!is.null(row_variable))    group_cols <- c(group_cols, row_variable)
  if (!is.null(column_variable)) group_cols <- c(group_cols, column_variable)

  if (use_donor_avg) {
    # Count cells per donor × group segment, then average across donors
    per_donor <- long_df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, id_column)))) |>
      dplyr::summarize(n_cells = dplyr::n(), .groups = "drop")

    agg_df <- per_donor |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarize(count    = mean(n_cells),
                       n_donors = dplyr::n(),
                       .groups  = "drop")

    # Total bar mean ± SE per x-position for error bars (ignoring fill splits)
    total_cols <- c("Variable", "Value")
    if (!is.null(row_variable))    total_cols <- c(total_cols, row_variable)
    if (!is.null(column_variable)) total_cols <- c(total_cols, column_variable)

    errorbar_df <- long_df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(total_cols, id_column)))) |>
      dplyr::summarize(n_total = dplyr::n(), .groups = "drop") |>
      dplyr::group_by(dplyr::across(dplyr::all_of(total_cols))) |>
      dplyr::summarize(
        total_mean = mean(n_total),
        total_se   = stats::sd(n_total) / sqrt(dplyr::n()),
        .groups    = "drop"
      )
  } else {
    agg_df <- long_df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarize(count = dplyr::n(), .groups = "drop")
  }

  # ── 8b. Convert to percentages (optional) ────────────────────────────────
  if (percent) {
    if (!is.null(fill_variable)) {
      # Each x-axis bar sums to 100%: divide within (Variable, Value, row, col)
      pct_group <- c("Variable", "Value")
      if (!is.null(row_variable))    pct_group <- c(pct_group, row_variable)
      if (!is.null(column_variable)) pct_group <- c(pct_group, column_variable)
      agg_df <- agg_df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(pct_group))) |>
        dplyr::mutate(pct = count / sum(count) * 100) |>
        dplyr::ungroup()
    } else {
      # No fill: show each bar as % of the panel total
      pct_group <- "Variable"
      if (!is.null(row_variable))    pct_group <- c(pct_group, row_variable)
      if (!is.null(column_variable)) pct_group <- c(pct_group, column_variable)
      agg_df <- agg_df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(pct_group))) |>
        dplyr::mutate(pct = count / sum(count) * 100) |>
        dplyr::ungroup()
    }
  }

  # ── 9. Resolve colors ─────────────────────────────────────────────────────
  if (!is.null(fill_variable)) {
    fill_colors <- colors
    if (is.null(fill_colors) && inherits(data, "Seurat"))
      fill_colors <- tryCatch(.nk_colors(data, fill_variable),
                              error = function(e) NULL)
    if (is.null(fill_colors))
      fill_colors <- SelectColors(agg_df[[fill_variable]], palette = "all")

    fill_lvls               <- names(fill_colors)
    agg_df[[fill_variable]] <- factor(as.character(agg_df[[fill_variable]]),
                                      levels = fill_lvls)
  }
  var_colors <- if (is.null(fill_variable))
    SelectColors(agg_df$Variable, palette = "spectrum") else NULL
  fmt_count  <- if (use_donor_avg) function(x) round(x) else function(x) x

  # ── 10-13. Plot builder (called once per page) ────────────────────────────
  .build_p <- function(df, err_df, facet_form, plot_title = NULL) {
    # Build base bar plot — use literal column names in aes() to avoid
    # deferred-evaluation issues when the plot is printed outside scope.
    if (is.null(fill_variable)) {
      pi <- if (percent)
        ggplot2::ggplot(df, ggplot2::aes(x = Value, y = pct, fill = Variable)) +
        ggplot2::geom_bar(stat = "identity", alpha = alpha) +
        ggplot2::scale_fill_manual(values = var_colors, name = "Variable")
      else
        ggplot2::ggplot(df, ggplot2::aes(x = Value, y = count, fill = Variable)) +
        ggplot2::geom_bar(stat = "identity", alpha = alpha) +
        ggplot2::scale_fill_manual(values = var_colors, name = "Variable")
    } else {
      pi <- if (percent)
        ggplot2::ggplot(df, ggplot2::aes(x = Value, y = pct,
                                         fill = .data[[fill_variable]])) +
        ggplot2::geom_bar(stat = "identity", alpha = alpha) +
        ggplot2::scale_fill_manual(values = fill_colors, name = fill_variable)
      else
        ggplot2::ggplot(df, ggplot2::aes(x = Value, y = count,
                                         fill = .data[[fill_variable]])) +
        ggplot2::geom_bar(stat = "identity", alpha = alpha) +
        ggplot2::scale_fill_manual(values = fill_colors, name = fill_variable)
    }

    # Error bars: only in raw-count mode (not meaningful after % rescaling)
    if (use_donor_avg && !percent && !is.null(err_df) && nrow(err_df) > 0L) {
      pi <- pi + ggplot2::geom_errorbar(
        data        = err_df,
        ggplot2::aes(x    = Value,
                     y    = total_mean,
                     ymin = pmax(total_mean - total_se, 0),
                     ymax = total_mean + total_se),
        width       = 0.2,
        linewidth   = 0.6,
        color      = "gray20",
        inherit.aes = FALSE
      )
    }

    pi <- pi + ggplot2::labs(x = "Value", y = y_label, title = plot_title)

    # Count / percent labels
    if (show_counts) {
      ldf <- df[df$count >= min_count_to_label, ]
      if (!is.null(fill_variable)) {
        pi <- pi + ggplot2::geom_text(
          data     = ldf,
          ggplot2::aes(label = if (percent) paste0(round(pct, 1), "%")
                               else         fmt_count(count)),
          position = ggplot2::position_stack(vjust = 0.5),
          size     = label_size, color = "white", fontface = "bold"
        )
      } else {
        pi <- pi + ggplot2::geom_text(
          data  = ldf,
          ggplot2::aes(label = if (percent) paste0(round(pct, 1), "%")
                               else         fmt_count(count)),
          vjust = -0.3,
          size  = label_size, color = "black"
        )
      }
    }

    pi +
      ggplot2::facet_grid(facet_form, scales = "free_x", space = "free_x") +
      theme_NourMin() +
      ggplot2::theme(
        axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
        strip.text    = ggplot2::element_text(face = "bold"),
        panel.spacing = ggplot2::unit(0.3, "lines"),
        plot.margin   = ggplot2::margin(t = 5, r = 10, b = 5, l = 15, unit = "mm"),
        axis.title.y  = ggplot2::element_text(angle = 90, vjust = 0.5, size = 11)
      ) +
      ggplot2::coord_cartesian(clip = "off")
  }

  # ── 12. Build plot list ───────────────────────────────────────────────────
  # column_variable → one plot per variable (multi-page PDF / list return)
  # no column_variable → single plot with Variable as column facet
  if (is.null(column_variable)) {
    ff  <- if (!is.null(row_variable))
      stats::as.formula(paste(row_variable, "~ Variable"))
    else
      stats::as.formula(". ~ Variable")
    err <- if (use_donor_avg && !percent) errorbar_df else NULL
    plots <- list(.build_p(agg_df, err, ff))
  } else {
    plots <- lapply(variables, function(var) {
      ff     <- if (!is.null(row_variable))
        stats::as.formula(paste(row_variable, "~", column_variable))
      else
        stats::as.formula(paste(". ~", column_variable))
      sub_df  <- agg_df[agg_df$Variable == var, , drop = FALSE]
      err_sub <- if (use_donor_avg && !percent)
        errorbar_df[errorbar_df$Variable == var, , drop = FALSE]
      else
        NULL
      .build_p(sub_df, err_sub, ff, plot_title = var)
    })
  }
  p <- plots[[1L]]

  # ── 14. Auto-save PDF if output_dir is available ─────────────────────────
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    n_row_facets <- if (!is.null(row_variable))
      length(unique(as.character(meta[[row_variable]]))) else 1L
    pdf_h <- n_row_facets * 3 + 1.5

    if (!is.null(column_variable)) {
      n_col_facets <- length(unique(as.character(meta[[column_variable]])))
      pdf_w        <- max(4.0, n_col_facets * 2.5 + 2.5)
    } else {
      n_x_ticks    <- vapply(variables, function(v)
        length(unique(as.character(meta[[v]]))), integer(1))
      pdf_w        <- sum(pmax(n_x_ticks * 0.35, 1.5)) + 2.5
    }

    if (!is.null(file_name) && nzchar(file_name)) {
      base <- file_name
    } else {
      parts <- c(
        if (nchar(object_name) > 0) object_name,
        if (!is.null(fill_variable)) fill_variable,
        paste(variables, collapse = "_"),
        if (!is.null(column_variable)) column_variable,
        if (percent) "Pct",
        "MetaSummary"
      )
      base <- paste(parts, collapse = "_")
    }
    fname <- gsub("[^A-Za-z0-9._-]", "_", base)
    fpath <- file.path(output_dir, paste0(fname, ".pdf"))

    grDevices::pdf(fpath, width = pdf_w, height = pdf_h)
    for (pl in plots) print(pl)
    grDevices::dev.off()
    message("scSidekick: Saved to ", fpath,
            " (", round(pdf_w, 1), " × ", round(pdf_h, 1), " in, ",
            length(plots), " page(s))")

    .write_legend_sidecar(fpath, paste0(
      "Stacked bar chart summarizing the distribution of categorical metadata variables",
      if (!is.null(id_column)) paste0(" deduplicated to one row per ", id_column) else
        " at the cell level",
      ". ",
      "Y-axis shows ", y_label, ". ",
      "Each column panel corresponds to a separate metadata variable (",
      paste(variables, collapse = ", "), "). ",
      if (!is.null(fill_variable))
        paste0("Bars are stacked and colored by ", fill_variable,
               ", with each segment labeled with its count. ")
      else
        "Each bar represents the total count per category with no fill grouping. ",
      if (!is.null(row_variable))
        paste0("Rows are split by ", row_variable, ". ")
      else "",
      if (!is.null(column_variable))
        paste0("When column_variable is set, one page is produced per variable in 'variables',",
               " with bars faceted by ", column_variable, ". ")
      else "",
      if (percent) "Values are shown as percentages rather than raw counts. " else "",
      if (!is.null(exclude) && length(exclude) > 0)
        paste0("The following groups were excluded before plotting: ",
               paste(mapply(function(col, vals)
                 paste0(col, " = ", paste(vals, collapse = ", ")),
                 names(exclude), exclude), collapse = "; "), ". ")
      else "",
      if (!is.null(object_name) && nchar(object_name) > 0)
        paste0("Dataset: ", object_name, ".")
      else ""
    ))
  }

  # ── 15. Assemble output ────────────────────────────────────────────────────
  has_numeric  <- !is.null(numeric_plot)
  multi_plots  <- !is.null(column_variable) && length(plots) > 1L

  if (!return_data && !return_flextable && !has_numeric) {
    return(if (multi_plots) plots else p)
  }

  result <- list(
    categorical_plot = if (multi_plots) plots else p,
    numeric_plot     = numeric_plot
  )

  if (return_data)
    result$data <- meta_dedup

  if (return_flextable) {
    if (!requireNamespace("flextable", quietly = TRUE)) {
      warning("Package 'flextable' is not installed - skipping flextable output. ",
              "Install with: install.packages(\"flextable\")")
    } else {
      # Donor-level wide table - only sensible when meta_dedup is small
      # (i.e. deduplicated to donors). With count_unit = "cells" meta_dedup
      # can be millions of rows; skip and warn rather than hang.
      if (nrow(meta_dedup) > 5000L) {
        warning("scSidekick: 'donor_table' flextable skipped - ",
                format(nrow(meta_dedup), big.mark = ","),
                " rows is too large. Set count_unit = \"donors\" and supply ",
                "id_column to get a meaningful per-donor table.")
        donor_ft <- NULL
      } else {
        donor_ft <- flextable::flextable(as.data.frame(meta_dedup)) |>
          flextable::autofit() |>
          flextable::theme_vanilla()
      }

      # Cross-tab per variable, all stacked into one flextable
      if (!is.null(fill_variable)) {
        ct_list <- lapply(variables, function(var) {
          tab <- meta_dedup |>
            dplyr::count(.data[[var]], .data[[fill_variable]]) |>
            tidyr::pivot_wider(
              names_from  = dplyr::all_of(fill_variable),
              values_from = n,
              values_fill = 0L
            )
          colnames(tab)[1] <- "Value"
          tab$Variable     <- var
          tab[, c("Variable", "Value",
                  setdiff(colnames(tab), c("Variable", "Value")))]
        })

        ct_combined <- do.call(rbind, ct_list)

        # Row positions where a separator line should appear (after each variable's block)
        sep_rows <- cumsum(vapply(ct_list, nrow, integer(1)))
        sep_rows <- sep_rows[-length(sep_rows)]   # no line after the last block

        crosstab_ft <- flextable::flextable(as.data.frame(ct_combined)) |>
          flextable::merge_v(j = "Variable") |>
          flextable::autofit() |>
          flextable::theme_vanilla()

        if (length(sep_rows) > 0L)
          crosstab_ft <- flextable::hline(crosstab_ft, i = sep_rows)

        result$flextable <- list(donor_table = donor_ft,
                                 crosstab    = crosstab_ft)
      } else {
        result$flextable <- list(donor_table = donor_ft)
      }
    }
  }

  result
}
