# =============================================================================
# scSidekick session logger
#
# start_session_logger() - opt-in passive logger.  Call once per session
#   (or add to ~/.Rprofile) to silently record every top-level command,
#   auto-save ggplot outputs, and capture base-R plots.
#
# Logs are written to ~/.R/session_logs/.  Use summarize_r_session() to
# generate a markdown + clean-script summary afterwards.
# =============================================================================

# Package-level environment that holds the logger state.
# Persists for the lifetime of the R session once logger is started.
.nk_logger_env <- new.env(parent = emptyenv())
.nk_logger_env$active       <- FALSE
.nk_logger_env$log_file     <- NULL
.nk_logger_env$plot_dir     <- NULL
.nk_logger_env$session_id   <- NULL
.nk_logger_env$plot_n       <- 0L
.nk_logger_env$last_gg_id   <- ""
.nk_logger_env$base_plt_rdy <- FALSE


#' Start the scSidekick passive session logger
#'
#' Attaches an R task callback that silently logs every top-level expression
#' evaluated in the session, together with its status (OK / ERR), timestamp,
#' and working directory.  ggplot2 figures are auto-saved as PNGs whenever a
#' new plot is detected; base-R graphics are captured on `plot.new`.
#'
#' Logs are written to `<log_root>/<session_id>_session.log` and plots to
#' `<log_root>/<session_id>_plots/`.  Call [summarize_r_session()] afterwards
#' to generate a markdown report and a clean reconstructed R script.
#'
#' **Recommended usage:** add the following line to your `~/.Rprofile` so the
#' logger starts automatically in every new R session:
#' ```r
#' if (requireNamespace("scSidekick", quietly = TRUE)) scSidekick::start_session_logger()
#' ```
#'
#' @param log_root Character.  Parent directory for log files.
#'   Defaults to `~/.R/session_logs`.
#' @param quiet Logical.  If `FALSE` (default) prints a startup message showing
#'   the log file path.
#'
#' @return Invisibly returns the path to the log file.
#' @export
#' @seealso [summarize_r_session()]
start_session_logger <- function(log_root = NULL,
                                 quiet    = FALSE) {

  e <- .nk_logger_env

  if (isTRUE(e$active)) {
    if (!quiet) message("[scSidekick logger] Already running → ", e$log_file)
    return(invisible(e$log_file))
  }

  if (is.null(log_root)) log_root <- path.expand("~/.R/session_logs")
  dir.create(log_root, recursive = TRUE, showWarnings = FALSE)

  e$session_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  e$log_file   <- file.path(log_root,
                            paste0(e$session_id, "_session.log"))
  e$plot_dir   <- file.path(log_root,
                            paste0(e$session_id, "_plots"))
  dir.create(e$plot_dir, showWarnings = FALSE)
  e$plot_n       <- 0L
  e$last_gg_id   <- ""
  e$base_plt_rdy <- FALSE
  e$active       <- TRUE

  # ---- helpers (closures over e) ----

  .nk_write <- function(...) {
    tryCatch(cat(paste0(..., "\n"), file = e$log_file, append = TRUE),
             error = function(err) NULL)
  }

  .nk_ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  .nk_next_fname <- function() {
    e$plot_n <- e$plot_n + 1L
    file.path(e$plot_dir,
              sprintf("plot%03d_%s.png", e$plot_n,
                      format(Sys.time(), "%H%M%S")))
  }

  .nk_save_gg <- function(p) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible())
    f <- .nk_next_fname()
    tryCatch({
      suppressMessages(ggplot2::ggsave(f, plot = p, width = 10,
                                       height = 7, dpi = 150))
      .nk_write("PLOT\t", .nk_ts(), "\tggplot\t", f)
    }, error = function(err)
      .nk_write("PLOT_ERR\t", .nk_ts(), "\t", conditionMessage(err)))
  }

  .nk_save_base <- function() {
    if (grDevices::dev.cur() == 1L) return(invisible())
    f <- .nk_next_fname()
    tryCatch({
      p <- grDevices::recordPlot()
      if (is.null(p) || length(p[[1]]) == 0) return(invisible())
      grDevices::png(f, width = 1500, height = 1050, res = 150)
      tryCatch(grDevices::replayPlot(p),
               error = function(err) NULL)
      grDevices::dev.off()
      .nk_write("PLOT\t", .nk_ts(), "\tbase_r\t", f)
    }, error = function(err) {
      tryCatch(grDevices::dev.off(), error = function(e2) NULL)
    })
  }

  # ---- task callback ----

  .nk_callback <- function(expr, value, ok, visible) {
    status    <- if (isTRUE(ok)) "OK" else "ERR"
    wd        <- tryCatch(getwd(), error = function(e2) "?")
    expr_text <- tryCatch(
      gsub("[\t\n\r]", " ", paste(deparse(expr), collapse = " ")),
      error = function(e2) "<parse-err>"
    )
    .nk_write("CMD\t", .nk_ts(), "\t", status, "\t", wd, "\t", expr_text)

    if (isNamespaceLoaded("ggplot2")) {
      lp <- tryCatch(ggplot2::last_plot(), error = function(x) NULL)
      if (!is.null(lp)) {
        lp_id <- tryCatch(
          paste(class(lp), length(lp$layers),
                paste(deparse(lp$mapping), collapse = ""), collapse = "|"),
          error = function(e2) ""
        )
        if (nchar(lp_id) > 0 && !identical(lp_id, e$last_gg_id)) {
          e$last_gg_id <- lp_id
          .nk_save_gg(lp)
        }
      }
    }
    TRUE
  }

  # ---- write header ----
  .nk_write("SESSION_START\t", e$session_id, "\t", as.character(Sys.time()))
  .nk_write("R_VERSION\t",    R.version$version.string)
  .nk_write("WD\t",           .nk_ts(), "\t", getwd())

  addTaskCallback(.nk_callback, name = "scSidekickLogger")

  # Capture the previous base-R plot before each new one starts
  setHook("plot.new", function() {
    if (e$base_plt_rdy) .nk_save_base()
    e$base_plt_rdy <- TRUE
  }, action = "append")

  # On session exit: save last base-R plot + write end marker
  reg.finalizer(.GlobalEnv, function(env) {
    tryCatch({
      if (!is.null(e$log_file) && isTRUE(e$active)) {
        if (e$base_plt_rdy) .nk_save_base()
        cat(paste0("SESSION_END\t", as.character(Sys.time()), "\n"),
            file = e$log_file, append = TRUE)
      }
    }, error = function(err) NULL)
  }, onexit = TRUE)

  if (!quiet)
    message("[scSidekick logger] ", e$session_id, " → ", e$log_file)

  invisible(e$log_file)
}


#' Return the current session log file path
#'
#' @return Character path, or `NULL` if the logger has not been started.
#' @export
#' @seealso [start_session_logger()]
logger_log_file <- function() .nk_logger_env$log_file

#' Return the current session plot directory
#'
#' @return Character path, or `NULL` if the logger has not been started.
#' @export
#' @seealso [start_session_logger()]
logger_plot_dir <- function() .nk_logger_env$plot_dir
