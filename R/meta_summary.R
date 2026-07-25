# =============================================================================
# meta_summary.R
#
# SummarizeMetadata() - collapse cell-level metadata to donor / patient level
# PlotMetaSummary()   - faceted bar chart of metadata distributions
# =============================================================================


# =============================================================================
# Private stat helpers (used by PlotMetaSummary)
# =============================================================================

.fmt_pval <- function(p) {
  if (is.na(p) || !is.finite(p)) return(NA_character_)
  if (p < 0.001) "< 0.001" else sprintf("%.3f", p)
}

# Categorical: Fisher's exact or chi-squared
.cat_pval <- function(tab, method) {
  tryCatch({
    if (method == "parametric") {
      stats::chisq.test(tab, correct = FALSE)$p.value
    } else if (method == "nonparametric") {
      use_sim <- nrow(tab) > 2L || ncol(tab) > 2L
      stats::fisher.test(tab, simulate.p.value = use_sim, B = 2000L)$p.value
    } else {
      expected <- stats::chisq.test(tab, correct = FALSE)$expected
      if (any(expected < 5, na.rm = TRUE)) {
        use_sim <- nrow(tab) > 2L || ncol(tab) > 2L
        stats::fisher.test(tab, simulate.p.value = use_sim, B = 2000L)$p.value
      } else {
        stats::chisq.test(tab, correct = FALSE)$p.value
      }
    }
  }, error = function(e) NA_real_)
}

# Numeric: Wilcoxon/Kruskal-Wallis or t-test/ANOVA
.num_pval <- function(vals, groups, method) {
  groups    <- factor(as.character(groups))
  n_groups  <- nlevels(groups)
  keep      <- !is.na(vals) & !is.na(groups)
  vals      <- vals[keep]; groups <- groups[keep]
  if (length(vals) < 2L || n_groups < 2L) return(NA_real_)

  use_np <- switch(method,
    nonparametric = TRUE,
    parametric    = FALSE,
    {
      g_split <- split(vals, groups)
      sw_ok   <- vapply(g_split, function(x)
        length(x) >= 3L && stats::shapiro.test(x)$p.value > 0.05,
        logical(1L))
      !all(sw_ok)
    }
  )

  tryCatch({
    if (n_groups == 2L) {
      g <- split(vals, groups)
      if (use_np) stats::wilcox.test(g[[1L]], g[[2L]], exact = FALSE)$p.value
      else        stats::t.test(g[[1L]], g[[2L]])$p.value
    } else {
      df_tmp <- data.frame(v = vals, g = groups)
      if (use_np) stats::kruskal.test(v ~ g, data = df_tmp)$p.value
      else        summary(stats::aov(v ~ g, data = df_tmp))[[1L]][["Pr(>F)"]][1L]
    }
  }, error = function(e) NA_real_)
}

# Column header label reflecting the test used
.pval_col_label <- function(method, n_groups, type = c("cat", "num")) {
  type <- match.arg(type)
  if (method == "auto") return("p.value")
  if (type == "cat") {
    if (method == "nonparametric") "p (Fisher)" else "p (Chi-sq)"
  } else {
    if (method == "nonparametric")
      if (n_groups == 2L) "p (Wilcoxon)" else "p (Kruskal-Wallis)"
    else
      if (n_groups == 2L) "p (t-test)" else "p (ANOVA)"
  }
}

# One-paragraph, table-specific explanation of exactly which tests were used,
# built from what is actually present in THIS table (categorical vars, numeric
# vars, group count, stratification) - not a generic boilerplate sentence.
# Used both as the in-table footer and as the .legend sidecar text.
.stats_footer_text <- function(method, has_cat, has_num, n_groups, has_strata,
                               cat_display, numeric_stat, run_stats) {
  parts <- character(0)

  if (run_stats) {
    test_bits <- character(0)
    if (has_cat) {
      test_bits <- c(test_bits, switch(method,
        nonparametric = "Fisher's exact test",
        parametric    = "Pearson's chi-squared test",
        auto          = "Pearson's chi-squared test (Fisher's exact test when any expected cell count < 5)",
        none          = NULL
      ))
    }
    if (has_num) {
      test_bits <- c(test_bits, switch(method,
        nonparametric = if (n_groups == 2L) "Wilcoxon rank-sum test" else "Kruskal-Wallis test",
        parametric    = if (n_groups == 2L) "Welch's t-test" else "one-way ANOVA",
        auto          = if (n_groups == 2L)
          "Welch's t-test (or Wilcoxon rank-sum test when a group fails a Shapiro-Wilk normality check)"
        else
          "one-way ANOVA (or Kruskal-Wallis test when a group fails a Shapiro-Wilk normality check)",
        none          = NULL
      ))
    }
    if (length(test_bits)) {
      lbl <- character(0)
      if (has_cat) lbl <- c(lbl, "categorical variables")
      if (has_num) lbl <- c(lbl, "numeric variables")
      stat_sentence <- paste0(
        "Statistical tests: ",
        paste(mapply(function(t, l) paste0(t, " (", l, ")"), test_bits, lbl),
              collapse = "; "), ".")
      parts <- c(parts, stat_sentence)
    }
    if (has_strata)
      parts <- c(parts,
        "Each stratum's p-value was computed independently within that stratum (not pooled across strata).")
    parts <- c(parts, "Bold p-values indicate p < 0.05.")
  }

  fmt_bit <- if (has_cat) switch(cat_display,
    n_pct    = "n (% of the column's group total)",
    pct_n    = "% (n) of the column's group total",
    fraction = "n/N (N = the column's group total)",
    n        = "n (raw count)"
  ) else NULL
  num_bit <- if (has_num) switch(numeric_stat,
    median_iqr = "median [IQR]",
    "mean ± SD"
  ) else NULL
  fmt_sentence <- paste0(
    "Values shown as: ",
    paste(c(if (!is.null(fmt_bit)) paste0(fmt_bit, " for categorical variables"),
            if (!is.null(num_bit)) paste0(num_bit, " for numeric variables")),
          collapse = "; "), ".")
  parts <- c(fmt_sentence, parts)

  paste(parts, collapse = " ")
}

# =============================================================================
# .build_meta_table1()
# Assemble a stratified "Table 1": every variable (categorical + numeric) in a
# single flextable, grouped by `fill_var`, optionally split into side-by-side
# blocks by `strata_var`, each block carrying its own p-value and group N.
# Categorical cells show `cat_display` (n / % column-wise); numeric rows show
# `numeric_stat` (mean +/- sd or median [IQR]) at whatever level `dd` is (the
# caller passes a donor-level frame when id_column is set). Returns a flextable
# or NULL when there is nothing to show.
# =============================================================================
.build_meta_table1 <- function(dd, all_vars, num_vars, fill_var, strata_var,
                               method, cat_display, numeric_stat, run_stats) {
  if (is.null(fill_var) || length(all_vars) == 0L) return(NULL)
  is_num <- function(v) v %in% num_vars

  fill_lvls <- if (is.factor(dd[[fill_var]])) levels(dd[[fill_var]])
               else sort(unique(as.character(dd[[fill_var]])))
  fill_lvls <- intersect(fill_lvls, unique(as.character(dd[[fill_var]][!is.na(dd[[fill_var]])])))
  if (length(fill_lvls) == 0L) return(NULL)

  strata <- if (is.null(strata_var)) list(All = dd) else {
    sl <- if (is.factor(dd[[strata_var]])) levels(dd[[strata_var]])
          else sort(unique(as.character(dd[[strata_var]])))
    sl <- intersect(sl, unique(as.character(dd[[strata_var]][!is.na(dd[[strata_var]])])))
    stats::setNames(lapply(sl, function(s)
      dd[!is.na(dd[[strata_var]]) & as.character(dd[[strata_var]]) == s, , drop = FALSE]), sl)
  }

  # skeleton rows (Variable, Value) in display order
  skel <- do.call(rbind, lapply(all_vars, function(v) {
    if (is_num(v)) data.frame(Variable = v, Value = "", stringsAsFactors = FALSE)
    else {
      lv <- if (is.factor(dd[[v]])) levels(dd[[v]])
            else sort(unique(as.character(dd[[v]][!is.na(dd[[v]])])))
      lv <- intersect(lv, unique(as.character(dd[[v]][!is.na(dd[[v]])])))
      data.frame(Variable = v, Value = lv, stringsAsFactors = FALSE)
    }
  }))

  fmt_cell <- function(n, N) {
    if (is.na(N) || N == 0L) return("0")
    pct <- round(100 * n / N)
    switch(cat_display,
      n_pct    = sprintf("%d (%d%%)", n, pct),
      pct_n    = sprintf("%d%% (%d)", pct, n),
      fraction = sprintf("%d/%d", n, N),
      n        = sprintf("%d", n))
  }
  fmt_num <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return("")
    if (numeric_stat == "median_iqr") {
      q <- stats::quantile(x, c(0.25, 0.5, 0.75))
      sprintf("%.2f [%.2f, %.2f]", q[2L], q[1L], q[3L])
    } else sprintf("%.2f ± %.2f", mean(x), stats::sd(x))
  }

  all_cat <- !any(vapply(all_vars, is_num, logical(1)))
  all_num <-  all(vapply(all_vars, is_num, logical(1)))
  p_lab   <- if (!run_stats) NULL
             else if (all_cat) .pval_col_label(method, length(fill_lvls), "cat")
             else if (all_num) .pval_col_label(method, length(fill_lvls), "num")
             else "p-value"

  df     <- skel
  blocks <- list()
  for (s in names(strata)) {
    sub <- strata[[s]]
    Ns  <- vapply(fill_lvls, function(f)
      sum(!is.na(sub[[fill_var]]) & as.character(sub[[fill_var]]) == f), integer(1))
    fcols <- character(0); flabels <- character(0)
    for (fi in seq_along(fill_lvls)) {
      f  <- fill_lvls[fi]
      sf <- sub[!is.na(sub[[fill_var]]) & as.character(sub[[fill_var]]) == f, , drop = FALSE]
      col <- vapply(seq_len(nrow(skel)), function(i) {
        v <- skel$Variable[i]; L <- skel$Value[i]
        if (is_num(v)) fmt_num(sf[[v]])
        else fmt_cell(sum(!is.na(sf[[v]]) & as.character(sf[[v]]) == L), sum(!is.na(sf[[v]])))
      }, character(1))
      cn <- paste(s, f, sep = "___")
      df[[cn]] <- col
      fcols   <- c(fcols, cn)
      flabels <- c(flabels, paste0(f, "\n(n=", Ns[fi], ")"))
    }
    pcn <- NULL
    if (run_stats) {
      pcol <- rep("", nrow(skel))
      for (v in all_vars) {
        idx <- which(skel$Variable == v)[1L]
        p <- if (is_num(v)) .num_pval(sub[[v]], sub[[fill_var]], method)
             else .cat_pval(table(as.character(sub[[v]]), as.character(sub[[fill_var]])), method)
        ps <- .fmt_pval(p)
        pcol[idx] <- if (is.na(ps)) "" else ps
      }
      pcn <- paste(s, "p", sep = "___")
      df[[pcn]] <- pcol
    }
    blocks[[s]] <- list(fcols = fcols, flabels = flabels, pcol = pcn,
                        span = length(fcols) + as.integer(run_stats))
  }

  ft <- flextable::flextable(as.data.frame(df, check.names = FALSE))
  lab_map <- list()
  for (s in names(blocks)) {
    b <- blocks[[s]]
    for (i in seq_along(b$fcols)) lab_map[[ b$fcols[i] ]] <- b$flabels[i]
    if (!is.null(b$pcol)) lab_map[[ b$pcol ]] <- p_lab
  }
  if (length(lab_map)) ft <- do.call(flextable::set_header_labels, c(list(ft), lab_map))

  if (!is.null(strata_var)) {
    top_vals   <- c("", "", names(blocks))
    top_widths <- c(1L, 1L, vapply(names(blocks), function(s) blocks[[s]]$span, integer(1)))
    ft <- flextable::add_header_row(ft, values = top_vals, colwidths = top_widths, top = TRUE)
    ft <- flextable::align(ft, i = 1, part = "header", align = "center")
  }

  ft <- flextable::merge_v(ft, j = "Variable")
  var_ends <- cumsum(rle(as.character(df$Variable))$lengths)
  var_ends <- var_ends[-length(var_ends)]
  if (length(var_ends)) ft <- flextable::hline(ft, i = var_ends)

  if (run_stats) {
    for (s in names(blocks)) {
      pcn <- blocks[[s]]$pcol
      if (is.null(pcn)) next
      pv  <- suppressWarnings(as.numeric(gsub("< ", "", df[[pcn]])))
      sig <- which(!is.na(pv) & pv < 0.05)
      if (length(sig)) ft <- flextable::bold(ft, i = sig, j = pcn)
    }
  }
  ft <- flextable::valign(ft, j = "Variable", valign = "top", part = "body")
  ft <- flextable::autofit(flextable::theme_vanilla(ft))

  # Table-specific legend: which test(s) were used (this is what was missing -
  # a mixed cat+numeric table previously showed a bare, unexplained "p-value"
  # header with no indication of which test each row actually used). Added
  # both as an in-table footer (travels with the .docx) and as an attribute
  # (reused verbatim for the .legend sidecar by .meta_emit_tables()).
  legend_txt <- .stats_footer_text(
    method        = method,
    has_cat       = any(!vapply(all_vars, is_num, logical(1))),
    has_num       = any(vapply(all_vars, is_num, logical(1))),
    n_groups      = length(fill_lvls),
    has_strata    = !is.null(strata_var),
    cat_display   = cat_display,
    numeric_stat  = numeric_stat,
    run_stats     = run_stats
  )
  ft <- flextable::add_footer_lines(ft, values = legend_txt)
  ft <- flextable::fontsize(ft, size = 7, part = "footer")
  ft <- flextable::italic(ft, part = "footer")
  attr(ft, "legend_text") <- legend_txt
  ft
}

# =============================================================================
# .meta_emit_tables()
# Orchestrates the flextable outputs for PlotMetaSummary: builds the unified
# Table 1, optionally the separate categorical / numeric tables, and autosaves
# each to .docx when output_dir is set. Returns the $flextable list.
# =============================================================================
.meta_emit_tables <- function(dd, all_vars, num_vars, fill_var, strata_var,
                              method, cat_display, numeric_stat, separate_tables,
                              output_dir, object_name, file_name, id_column = NULL) {
  run_stats <- method != "none" && !is.null(fill_var)
  cat_vars  <- setdiff(all_vars, num_vars)

  out <- list(summary_table = .build_meta_table1(
    dd, all_vars, num_vars, fill_var, strata_var, method, cat_display,
    numeric_stat, run_stats))

  if (isTRUE(separate_tables)) {
    if (length(cat_vars) > 0L)
      out$crosstab <- .build_meta_table1(dd, cat_vars, character(0), fill_var,
        strata_var, method, cat_display, numeric_stat, run_stats)
    if (length(num_vars) > 0L)
      out$numeric_table <- .build_meta_table1(dd, num_vars, num_vars, fill_var,
        strata_var, method, cat_display, numeric_stat, run_stats)
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    # Build a descriptive, unique name that mirrors the companion plot's
    # (variables_fill_strata), so the table and its plot sit side by side and
    # different variable/grouping combinations never overwrite each other.
    base <- if (!is.null(file_name) && nzchar(file_name)) file_name else
      paste(c(if (nchar(object_name) > 0) object_name,
              paste(all_vars, collapse = "_"),
              fill_var, strata_var, "MetaSummary"),
            collapse = "_")
    base <- gsub("[^A-Za-z0-9._-]", "_", base)
    for (nm in names(out)) {
      if (is.null(out[[nm]])) next
      suffix <- if (nm == "summary_table") "" else paste0("_", nm)
      fp <- file.path(output_dir, paste0(base, suffix, ".docx"))
      ok <- tryCatch({ flextable::save_as_docx(out[[nm]], path = fp); TRUE },
                     error = function(e) {
                       warning("Could not save ", nm, " to .docx: ",
                               conditionMessage(e)); FALSE })
      if (ok) {
        message("scSidekick: Saved table to ", fp)
        # .legend sidecar, matching every other autosaved output in the
        # package. Leads with what the table shows, then the exact test/format
        # legend already computed for (and printed as a footer inside) the
        # table itself, so the two never drift out of sync.
        tbl_vars <- switch(nm,
          crosstab      = cat_vars,
          numeric_table = num_vars,
          all_vars)
        n_unit <- if (!is.null(id_column))
          paste0("unique '", id_column, "' values") else "rows"
        opening <- paste0(
          "Table summarizing ", paste(tbl_vars, collapse = ", "),
          " by ", fill_var,
          if (!is.null(strata_var)) paste0(", stratified by ", strata_var) else "",
          ". N = ", format(nrow(dd), big.mark = ","), " ", n_unit,
          if (nchar(object_name) > 0) paste0(" [", object_name, "]") else "", "."
        )
        .write_legend_sidecar(fp, paste(opening, attr(out[[nm]], "legend_text")))
      }
    }
  }
  out
}


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
#' @param spread Character or \code{NULL}.  When \code{"sd"} or \code{"sem"},
#'   adds a companion column for each numeric variable (e.g. \code{Age_sd} or
#'   \code{Age_sem}).  At the patient level these are computed across cells;
#'   when \code{group.by} is also supplied they reflect the spread of patient
#'   means across patients in each group.  Default \code{NULL} (no spread).
#' @param columns Character vector of column names to include in the summary.
#'   \code{NULL} (default) includes all columns except \code{id_columns}.
#'   \code{id_columns} are always retained regardless of this setting.
#' @param group.by Optional column name to aggregate the patient-level result
#'   up to a group level.  When supplied, a two-stage aggregation is performed:
#'   (1) cells are collapsed to one row per \code{id_columns} combination
#'   (current behavior), then (2) those patient-level rows are averaged by
#'   \code{group.by}, giving each patient equal weight regardless of cell
#'   count.  The returned data frame has one row per \code{group.by} level and
#'   gains an \code{n_patients} column.  Default \code{NULL} (single-stage).
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
                               columns       = NULL,
                               group.by      = NULL,
                               numeric_func  = mean,
                               numeric_na_rm = TRUE,
                               spread        = NULL,
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
  if (!is.null(columns)) {
    bad <- setdiff(columns, colnames(meta))
    if (length(bad) > 0L)
      stop("The following columns were not found: ", paste(bad, collapse = ", "))
    other_cols <- intersect(other_cols, columns)
  }
  num_cols   <- other_cols[vapply(meta[other_cols], is.numeric, logical(1))]
  cat_cols   <- setdiff(other_cols, num_cols)

  if (!is.null(spread) && !spread %in% c("sd", "sem"))
    stop("'spread' must be NULL, \"sd\", or \"sem\".")

  .spread_fn <- if (identical(spread, "sem"))
    \(v) { v <- v[!is.na(v)]; if (length(v) < 2L) NA_real_ else sd(v) / sqrt(length(v)) }
  else
    \(v) sd(v, na.rm = TRUE)

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

    if (!is.null(spread)) {
      spread_summary <- meta |>
        dplyr::group_by(dplyr::across(dplyr::all_of(id_columns))) |>
        dplyr::summarize(
          dplyr::across(dplyr::all_of(num_cols), .spread_fn),
          .groups = "drop"
        ) |>
        dplyr::rename_with(\(x) paste0(x, "_", spread), dplyr::all_of(num_cols))
      result <- dplyr::left_join(result, spread_summary, by = id_columns)
    }
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

  # ── 7. Optional second-stage: average patient rows by group.by ───────────
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(meta))
      stop("'group.by' column '", group.by, "' not found in metadata.")
    if (group.by %in% id_columns)
      stop("'group.by' column '", group.by, "' cannot also be an id_column.")

    # Attach group.by to patient-level result if not already present
    if (!group.by %in% colnames(result)) {
      grp_map <- unique(meta[, c(id_columns, group.by), drop = FALSE])
      result  <- dplyr::left_join(result, grp_map, by = id_columns)
    }

    # Identify numeric vs categorical among patient-level summary columns.
    # Exclude any _sd/_sem columns added in stage 1 - they will be recomputed
    # at the group level rather than averaged (mean-of-SDs is not meaningful).
    pat_cols      <- setdiff(colnames(result), c(id_columns, group.by, "n_cells"))
    spread_suffix <- if (!is.null(spread)) paste0("_", spread) else ""
    spread_pat    <- if (!is.null(spread))
      grep(paste0(spread_suffix, "$"), pat_cols, value = TRUE) else character(0L)
    pat_num_cols  <- setdiff(pat_cols[vapply(result[pat_cols], is.numeric, logical(1L))],
                             spread_pat)
    pat_cat_cols  <- setdiff(pat_cols, c(pat_num_cols, spread_pat))

    grp_result <- result |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group.by))) |>
      dplyr::summarize(
        n_patients = dplyr::n(),
        n_cells    = sum(n_cells, na.rm = TRUE),
        dplyr::across(dplyr::all_of(pat_num_cols), \(v) mean(v, na.rm = TRUE)),
        .groups    = "drop"
      )

    if (!is.null(spread) && length(pat_num_cols) > 0L) {
      grp_spread <- result |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group.by))) |>
        dplyr::summarize(
          dplyr::across(dplyr::all_of(pat_num_cols), .spread_fn),
          .groups = "drop"
        ) |>
        dplyr::rename_with(\(x) paste0(x, "_", spread), dplyr::all_of(pat_num_cols))
      grp_result <- dplyr::left_join(grp_result, grp_spread, by = group.by)
    }

    for (col in pat_cat_cols) {
      val_counts <- result |>
        dplyr::count(dplyr::across(dplyr::all_of(c(group.by, col))),
                     name = ".n_pat_inner")

      col_summary <- val_counts |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group.by))) |>
        dplyr::summarize(
          !!col := {
            vals   <- as.character(.data[[col]])
            counts <- .data$.n_pat_inner
            if (dplyr::n() == 1L) {
              vals[1L]
            } else {
              ord <- order(counts, decreasing = TRUE)
              paste(paste0(vals[ord], " (n=", counts[ord], ")"), collapse = ", ")
            }
          },
          .groups = "drop"
        )
      grp_result <- dplyr::left_join(grp_result, col_summary, by = group.by)
    }

    message("scSidekick SummarizeMetadata: aggregated to ",
            nrow(grp_result), " group(s) via '", group.by, "'.")
    result <- grp_result
  }

  # ── 8. Optional pivot to long ─────────────────────────────────────────────
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
#'   \code{flextable} objects: \code{$donor_table} (wide donor-level data),
#'   \code{$crosstab} (n per category level per \code{fill_variable} group),
#'   and \code{$numeric_table} (mean ± SD per numeric variable per
#'   \code{fill_variable} group, with each donor weighted equally).
#'   Requires the \pkg{flextable} package.  Default \code{FALSE}.
#' @param stats.method Character.  Statistical test to add as a \code{p.value}
#'   column in the flextable outputs (requires \code{return_flextable = TRUE}
#'   and \code{fill_variable}).
#'   \describe{
#'     \item{\code{"nonparametric"} (default)}{Fisher's exact for categorical,
#'       Wilcoxon / Kruskal-Wallis for numeric.}
#'     \item{\code{"parametric"}}{Chi-squared for categorical,
#'       t-test / ANOVA for numeric.}
#'     \item{\code{"auto"}}{Parametric when Shapiro-Wilk passes (numeric) or
#'       expected counts >= 5 (categorical); otherwise non-parametric.}
#'     \item{\code{"none"}}{No statistics added.}
#'   }
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
                             return_flextable   = FALSE,
                             numeric_stat       = c("mean_sd", "median_iqr"),
                             cat_display        = c("n_pct", "pct_n", "fraction", "n"),
                             separate_tables    = FALSE,
                             stats.method       = c("nonparametric", "parametric", "auto", "none"),
                             pdf.width          = NULL,
                             pdf.height         = NULL) {

  count_unit   <- match.arg(count_unit)
  stats.method <- match.arg(stats.method)
  numeric_stat <- match.arg(numeric_stat)
  cat_display  <- match.arg(cat_display)

  # ── 0. Walk up PrepObject defaults (Seurat only) ──────────────────────────
  if (inherits(data, "Seurat")) {
    if (missing(output_dir))
      output_dir <- if (.nk_autosave(data)) .nk_setting(data, "output_dir") else NULL
    object_name <- if (nchar(object_name) > 0) object_name else
      .nk_setting(data, "object_name") %||% ""
    if (is.null(id_column)) .nk_warn_donor(data)
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
  all_variables <- variables   # original order (mixed) - used for the table
  numeric_vars <- variables[vapply(variables,
    function(v) v %in% colnames(meta) && is.numeric(meta[[v]]), logical(1))]
  cat_vars <- setdiff(variables, numeric_vars)

  # Numeric variables are summarized IN the table (mean +/- sd or median [IQR]).
  # A companion PLOT routes to a donor-level function (PlotPseudoBulk) when an
  # id_column defines donors, or a cell-level one (PlotFeature) when it does not,
  # so the language and statistics match the requested unit.
  numeric_plot <- NULL
  if (length(numeric_vars) > 0L && !is.null(fill_variable)) {
    numeric_plot <- tryCatch(
      if (!is.null(id_column)) {
        message("scSidekick: ", paste(numeric_vars, collapse = ", "),
                " numeric - donor-level summary; companion plot via PlotPseudoBulk().")
        # add_stats = FALSE: the p-values belong in the Table 1 (donor-level,
        # correct test); the companion plot is illustrative only, so no brackets.
        PlotPseudoBulk(seurat_object = data, features = numeric_vars,
                       group.by = fill_variable, donor.by = id_column,
                       split.by = row_variable, colors = colors, add_stats = FALSE,
                       output_dir = output_dir, object_name = object_name)
      } else {
        message("scSidekick: ", paste(numeric_vars, collapse = ", "),
                " numeric, no id_column - cell-level plot via PlotFeature().")
        # add_stats = FALSE: cell-level tests are pseudoreplicated; keep the plot
        # illustrative and leave inference to the (donor-level) table.
        PlotFeature(data = data, features = numeric_vars, group.by = fill_variable,
                    split.by = row_variable, colors = colors, add_stats = FALSE,
                    output_dir = output_dir, object_name = object_name)
      },
      error = function(e) {
        warning("Numeric companion plot failed: ", conditionMessage(e)); NULL })
  }

  # ── All-numeric: no categorical bars. Build the summary table and return. ──
  if (length(cat_vars) == 0L) {
    result <- list(categorical_plot = NULL, numeric_plot = numeric_plot)
    if (return_flextable && requireNamespace("flextable", quietly = TRUE) &&
        !is.null(fill_variable)) {
      dd <- meta
      if (!is.null(exclude))
        for (col in names(exclude))
          if (col %in% colnames(dd))
            dd <- dd[!as.character(dd[[col]]) %in% as.character(exclude[[col]]),
                     , drop = FALSE]
      need <- intersect(c(id_column, all_variables, fill_variable, row_variable),
                        colnames(dd))
      dd   <- dd[, need, drop = FALSE]
      if (!is.null(id_column)) dd <- dplyr::distinct(dd)
      result$flextable <- .meta_emit_tables(dd, all_variables, numeric_vars,
        fill_variable, row_variable, stats.method, cat_display, numeric_stat,
        separate_tables, output_dir, object_name, file_name, id_column)
    }
    if (return_data) result$data <- meta
    if (!return_flextable && !return_data) return(numeric_plot)
    return(result)
  }

  variables <- cat_vars   # categorical bar-chart path uses only these

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
    pdf_h <- pdf.height %||% (n_row_facets * 3 + 1.5)

    if (!is.null(column_variable)) {
      n_col_facets <- length(unique(as.character(meta[[column_variable]])))
      pdf_w        <- pdf.width %||% max(4.0, n_col_facets * 2.5 + 2.5)
    } else {
      n_x_ticks    <- vapply(variables, function(v)
        length(unique(as.character(meta[[v]]))), integer(1))
      pdf_w        <- pdf.width %||% (sum(pmax(n_x_ticks * 0.35, 1.5)) + 2.5)
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

    ms_n_obs    <- format(nrow(meta), big.mark = ",")
    ms_unit     <- if (inherits(data, "Seurat")) .nk_unit_label(data) else "cells"
    ms_n_donors <- if (!is.null(id_column) && id_column %in% colnames(meta))
      length(unique(meta[[id_column]]))
    else if (inherits(data, "Seurat")) {
      don_col <- .nk_setting(data, "donor.by")
      if (!is.null(don_col) && don_col %in% colnames(meta))
        length(unique(meta[[don_col]])) else NULL
    } else NULL
    .write_legend_sidecar(fpath, paste0(
      "Stacked bar chart of ", ms_n_obs, " ", ms_unit,
      if (!is.null(ms_n_donors))
        paste0(" from ", format(ms_n_donors, big.mark = ","),
               " ", if (!is.null(id_column)) id_column else "donors")
      else "",
      if (nchar(object_name) > 0) paste0(" [", object_name, "]") else "",
      " summarizing the distribution of categorical metadata variables",
      if (!is.null(id_column)) paste0(" (deduplicated to one row per ", id_column, ")") else "",
      ". Y-axis shows ", y_label, ". ",
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
        paste0("One page is produced per variable, with bars faceted by ",
               column_variable, ". ")
      else "",
      if (percent) "Values are shown as percentages rather than raw counts. " else "",
      if (!is.null(exclude) && length(exclude) > 0)
        paste0("The following groups were excluded: ",
               paste(mapply(function(col, vals)
                 paste0(col, " = ", paste(vals, collapse = ", ")),
                 names(exclude), exclude), collapse = "; "), ". ")
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
    } else if (is.null(fill_variable)) {
      warning("scSidekick: the flextable summary needs a fill_variable to ",
              "group by - skipping table output.")
    } else {
      # Donor-level frame (dedup on id_column) so the Table 1 is per donor, then
      # emit the unified table (+ optional separate tables) and autosave to .docx.
      need <- intersect(c(id_column, all_variables, fill_variable, row_variable),
                        colnames(meta))
      dd   <- meta[, need, drop = FALSE]
      if (!is.null(id_column)) dd <- dplyr::distinct(dd)
      result$flextable <- .meta_emit_tables(dd, all_variables, numeric_vars,
        fill_variable, row_variable, stats.method, cat_display, numeric_stat,
        separate_tables, output_dir, object_name, file_name, id_column)
    }
  }

  result
}
