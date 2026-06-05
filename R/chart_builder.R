# =============================================================================
# scSidekick — Chart Builder launcher
#
# ChartBuilder — open the standalone HTML/JS Chart Builder app in a browser
#   with an R data frame (or a Seurat object's metadata) pre-loaded, so the
#   user lands directly on the data-setup screen instead of dragging a file in.
#
# Mechanism (launcher pattern, one-way R -> browser):
#   1. Coerce the input to a clean data frame.
#   2. Snapshot the Chart Builder app into a private temp session directory.
#   3. Serialise the data frame to `nourkit_payload.js` (a global
#      NOURKIT_PAYLOAD = { filename, rows: [...] }).
#   4. Inject a <script> tag for the payload into the copied HTML.
#   5. Open the copied HTML; an auto-load hook in the app routes the payload
#      through the same code path as a manual CSV drop.
#
# The finished figure is exported from within the app (SVG/PNG/PDF); this is
# an interactive, human-in-the-loop launcher, not a headless renderer.
# =============================================================================

#' Open the Chart Builder app pre-loaded with an R data frame
#'
#' Launches the standalone Chart Builder (an HTML/CSS/JS Plotly app) in the
#' default browser with `data` already loaded, dropping the user straight onto
#' the data-setup screen. This is a one-way, interactive launcher: data flows
#' from R into the browser, and the finished chart is exported from within the
#' app. There is no automated round-trip back to R.
#'
#' The app is copied into a fresh temporary directory on each call, so the
#' source app is never modified and concurrent sessions do not collide. The app
#' loads its charting libraries (Plotly, PapaParse, SheetJS) from CDNs, so an
#' internet connection is required for it to render.
#'
#' @param data A `data.frame` (or object coercible with [as.data.frame()]),
#'   or a Seurat object — in which case its `meta.data` slot is used.
#' @param app_dir Character. Path to the Chart Builder app directory (the folder
#'   containing `chart_builder.html`). Defaults to
#'   `getOption("nourkit.chart_builder_dir")`, falling back to
#'   `~/Desktop/Claude Tools/Chart_Builder`.
#' @param filename Character. Name shown in the app for the loaded dataset.
#'   Defaults to the deparsed name of `data` with a `.csv` suffix.
#' @param max_rows Integer or `NULL`. If set and `data` has more rows, a random
#'   sample of `max_rows` rows is sent (keeps the payload small). `NULL`
#'   (default) sends everything; the app itself samples large data for
#'   plotting. A warning is issued for very large frames.
#' @param launch Logical. Open the browser? Default `TRUE`. `FALSE` builds the
#'   session directory and returns its path without opening (useful for tests).
#'
#' @return (Invisibly) the path to the temporary session directory containing
#'   the launched copy of the app.
#' @keywords internal
ChartBuilder <- function(data,
                         app_dir  = getOption(
                           "nourkit.chart_builder_dir",
                           "~/Desktop/Claude Tools/Chart_Builder"
                         ),
                         filename = NULL,
                         max_rows = NULL,
                         launch   = TRUE) {

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("ChartBuilder() requires the 'jsonlite' package. ",
         "Install it with install.packages('jsonlite').", call. = FALSE)
  }

  # ── Default filename from the symbol the user passed ───────────────────────
  if (is.null(filename)) {
    nm <- deparse(substitute(data))
    if (length(nm) != 1L || !grepl("^[A-Za-z][A-Za-z0-9._]*$", nm)) nm <- "data"
    filename <- paste0(nm, ".csv")
  }

  # ── Coerce input to a clean data frame ─────────────────────────────────────
  df <- .nk_cb_as_dataframe(data)
  if (nrow(df) == 0L) stop("`data` has no rows to plot.", call. = FALSE)
  if (ncol(df) == 0L) stop("`data` has no columns to plot.", call. = FALSE)

  # Optional R-side downsampling
  if (!is.null(max_rows) && nrow(df) > max_rows) {
    df <- df[sort(sample.int(nrow(df), max_rows)), , drop = FALSE]
  } else if (nrow(df) > 2e5) {
    warning("`data` has ", nrow(df), " rows; the payload file may be large. ",
            "Consider setting `max_rows` to downsample before launching.",
            call. = FALSE)
  }

  # ── Resolve the app directory ──────────────────────────────────────────────
  app_dir   <- path.expand(app_dir)
  html_src  <- file.path(app_dir, "chart_builder.html")
  if (!dir.exists(app_dir) || !file.exists(html_src)) {
    stop("Could not find 'chart_builder.html' in:\n  ", app_dir, "\n",
         "Set the location with options(nourkit.chart_builder_dir = '/path/to/Chart_Builder').",
         call. = FALSE)
  }

  # ── Snapshot the app into a private temp session dir ───────────────────────
  session_dir <- tempfile("nourkit_chartbuilder_")
  dir.create(session_dir)
  assets <- list.files(app_dir, pattern = "\\.(html|js|css)$", full.names = TRUE)
  file.copy(assets, session_dir, overwrite = TRUE)

  # ── Write the data payload ────────────────────────────────────────────────
  rows_json <- jsonlite::toJSON(
    df, dataframe = "rows", na = "null", null = "null", auto_unbox = TRUE
  )
  payload <- paste0(
    "// Auto-generated by scSidekick::ChartBuilder() — do not edit.\n",
    "var NOURKIT_PAYLOAD = {\n",
    "  filename: ", jsonlite::toJSON(filename, auto_unbox = TRUE), ",\n",
    "  rows: ", rows_json, "\n",
    "};\n"
  )
  writeLines(payload, file.path(session_dir, "nourkit_payload.js"), useBytes = TRUE)

  # ── Inject the payload <script> into the copied HTML ───────────────────────
  html_dst  <- file.path(session_dir, "chart_builder.html")
  html      <- readLines(html_dst, warn = FALSE, encoding = "UTF-8")
  anchor    <- grep("<script[^>]+src=[\"']config\\.js[\"']", html)[1]
  tag       <- "<script src=\"nourkit_payload.js\"></script>"
  if (is.na(anchor)) {
    # Fall back to inserting before </body>
    anchor <- grep("</body>", html, ignore.case = TRUE)[1]
    if (is.na(anchor)) stop("Could not find an injection point in chart_builder.html.",
                            call. = FALSE)
  }
  html <- append(html, tag, after = anchor - 1L)
  writeLines(html, html_dst, useBytes = TRUE)

  # ── Launch ─────────────────────────────────────────────────────────────────
  if (launch) {
    utils::browseURL(html_dst)
    message("Chart Builder launched with '", filename, "' (",
            nrow(df), " rows x ", ncol(df), " cols).")
  }

  invisible(session_dir)
}

# ── Internal: coerce supported inputs to a plotting-friendly data frame ───────
.nk_cb_as_dataframe <- function(data) {
  # Seurat object -> cell metadata
  if (methods::is(data, "Seurat")) {
    md <- data@meta.data
    md <- cbind(cell = rownames(md), md)
    rownames(md) <- NULL
    data <- md
  }

  df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)

  # Make column-typed values JSON/Plotly-friendly:
  #   factors -> character, Dates/times -> ISO character, list-cols collapsed.
  df[] <- lapply(df, function(col) {
    if (is.factor(col))                       return(as.character(col))
    if (inherits(col, c("Date", "POSIXct")))  return(format(col))
    if (is.list(col))                         return(vapply(col, function(x)
      paste(as.character(x), collapse = "; "), character(1)))
    col
  })
  df
}
