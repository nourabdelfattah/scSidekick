# =============================================================================
# scSidekick session summary
#
# summarize_r_session() — rule-based summarizer for R session logs produced
#   by start_session_logger() or any ~/.Rhistory file.
#
# Outputs a markdown summary + a deduplicated clean R script.
# =============================================================================

#' Summarize an R session log
#'
#' Parses a structured session log created by [start_session_logger()] (or a
#' plain `.Rhistory` file) and writes two files into `output_dir`:
#'
#' * `session_YYYYMMDD_summary.md` — a markdown report with phase-by-phase
#'   summaries, error counts, and a plot inventory.
#' * `session_YYYYMMDD_clean.R`   — a deduplicated, phase-annotated R script
#'   containing only the successful commands.
#'
#' @param log_file Character or `NULL`.  Path to a `.log` file written by
#'   [start_session_logger()].  When `NULL` (default) the most recent log in
#'   `~/.R/session_logs/` is used.
#' @param rhistory Character or `NULL`.  Path to an `.Rhistory` file to
#'   summarize instead of a structured log.  Takes precedence over `log_file`.
#' @param output_dir Character or `NULL`.  Directory for the output files.
#'   Defaults to the directory containing `log_file` (or `rhistory`), or the
#'   current working directory.
#' @param open Logical.  If `TRUE` (default) and the session is interactive,
#'   opens the markdown file in the editor after writing.
#'
#' @return Invisibly returns a named list with `$md`, `$r`, and `$data`.
#' @export
summarize_r_session <- function(log_file  = NULL,
                                rhistory  = NULL,
                                output_dir = NULL,
                                open      = TRUE) {

  if (!is.null(rhistory)) {
    df <- .sr_parse_rhistory(path.expand(rhistory))
  } else {
    if (is.null(log_file)) {
      log_root <- path.expand("~/.R/session_logs")
      logs <- list.files(log_root, "_session\\.log$", full.names = TRUE)
      if (!length(logs))
        stop("No session logs in ", log_root,
             ".\nIs the session logger running? ",
             "Call scSidekick::start_session_logger() or add it to ~/.Rprofile.")
      log_file <- sort(logs, decreasing = TRUE)[1]
      message("Using: ", log_file)
    }
    df <- .sr_parse_log(log_file)
  }

  df$cat    <- vapply(df$command, .sr_categorize,   character(1))
  df$label  <- vapply(df$cat,    .sr_cat_label,     character(1))
  df$params <- mapply(.sr_extract_params, df$command, df$cat,
                      USE.NAMES = FALSE, SIMPLIFY = TRUE)

  phases <- .sr_detect_phases(df)

  if (is.null(output_dir)) {
    if (!is.null(log_file)) {
      output_dir <- dirname(log_file)
    } else if (!is.null(rhistory)) {
      output_dir <- dirname(path.expand(rhistory))
    } else {
      output_dir <- getwd()
    }
  }

  ts_vals  <- df$timestamp[!is.na(df$timestamp)]
  date_str <- if (length(ts_vals)) format(min(ts_vals), "%Y-%m-%d") else
                format(Sys.Date())
  base    <- file.path(output_dir, paste0("session_", gsub("-", "", date_str)))
  md_file <- paste0(base, "_summary.md")
  r_file  <- paste0(base, "_clean.R")

  writeLines(.sr_gen_markdown(df, phases), md_file)
  writeLines(.sr_gen_clean_r(df, phases),  r_file)

  message("Written:\n  ", md_file, "\n  ", r_file)
  if (open && interactive())
    tryCatch(utils::file.edit(md_file), error = function(e) NULL)

  invisible(list(md = md_file, r = r_file, data = df))
}

# ==============================================================
# Parsers
# ==============================================================

.sr_parse_log <- function(path) {
  lines     <- readLines(path, warn = FALSE)
  cmd_lines <- grep("^CMD\t", lines, value = TRUE)

  if (!length(cmd_lines)) {
    message("No CMD lines found; treating as raw history")
    return(.sr_parse_rhistory(path))
  }

  parts <- strsplit(cmd_lines, "\t")
  df <- data.frame(
    timestamp = as.POSIXct(vapply(parts, `[`, character(1), 2),
                           format = "%Y-%m-%d %H:%M:%S"),
    status    = vapply(parts, `[`, character(1), 3),
    wd        = vapply(parts, `[`, character(1), 4),
    command   = vapply(parts, function(x)
                  paste(x[5:length(x)], collapse = "\t"), character(1)),
    stringsAsFactors = FALSE
  )

  plot_lines <- grep("^PLOT\t", lines, value = TRUE)
  attr(df, "n_plots")    <- length(plot_lines)
  attr(df, "plot_files") <- if (length(plot_lines))
    vapply(strsplit(plot_lines, "\t"), `[`, character(1), 4) else character(0)

  start_l <- grep("^SESSION_START\t", lines, value = TRUE)
  end_l   <- grep("^SESSION_END\t",   lines, value = TRUE)
  attr(df, "session_id")  <- if (length(start_l))
    strsplit(start_l[1], "\t")[[1]][2] else "?"
  attr(df, "session_end") <- if (length(end_l))
    strsplit(end_l[1], "\t")[[1]][2] else ""
  df
}

.sr_parse_rhistory <- function(path) {
  lines  <- readLines(path, warn = FALSE)
  lines  <- lines[nchar(trimws(lines)) > 0]
  joined <- .sr_join_continuations(lines)
  df <- data.frame(
    timestamp = as.POSIXct(NA),
    status    = NA_character_,
    wd        = NA_character_,
    command   = joined,
    stringsAsFactors = FALSE
  )
  attr(df, "n_plots")    <- 0L
  attr(df, "plot_files") <- character(0)
  attr(df, "session_id") <- "rhistory"
  df
}

.sr_join_continuations <- function(lines) {
  result  <- character(0)
  current <- ""
  END_PAT   <- "[+|,%>({\\[]\\s*$"
  START_PAT <- "^\\s*[+|)\\],]"
  for (ln in lines) {
    if (nchar(current) == 0) {
      current <- ln
    } else if (grepl(END_PAT, current) || grepl(START_PAT, ln)) {
      current <- paste(current, trimws(ln))
    } else {
      result  <- c(result, current)
      current <- ln
    }
  }
  if (nchar(current) > 0) result <- c(result, current)
  result
}

# ==============================================================
# Categorization table — ordered, first match wins
# ==============================================================

.SR_CATS <- data.frame(
  name  = c(
    "packages","load_obj","create_seurat","open_bpcells","read_file",
    "preproc","pca","batch","umap","clustering",
    "qc","markers","dim_viz","feat_viz","dot_vln","heatmap",
    "ggplot","graph_viz","table_insp","view",
    "setwd","save_out","annot","subset_op","integrate",
    "knn_proj","func_def","dir_setup","trajectory","ccc",
    "recover","comment","other"
  ),
  pat = c(
    "^(library|require|libraries)\\(",
    "^(load|readRDS)\\(",
    "CreateSeuratObject|CreateAssay5Object",
    "open_matrix_dir",
    "read\\.csv\\(|read_csv\\(|fread\\(",
    "NormalizeData|ScaleData|FindVariableFeatures|SCTransform|SketchData|CellCycleScoring",
    "RunPCA|Determine_PCs|ElbowPlot|VizDimLoadings|print.*\\[\\[.pca",
    "RunHarmony|harmony::",
    "RunUMAP|RunTSNE",
    "FindNeighbors|FindClusters",
    "PercentageFeatureSet|mitoRatio|Doublet\\.score|table\\(|range\\(",
    "FindAllMarkers|FindMarkers|top_n\\(",
    "DimPlot|SpatialDimPlot",
    "FeaturePlot|SpatialFeaturePlot",
    "DotPlot|VlnPlot|RidgePlot",
    "DoHeatmap|pheatmap|Heatmap\\(",
    "^ggplot\\(",
    "ggraph|igraph|visNetwork",
    "^head\\(|^str\\(|^dim\\(|^summary\\(|^colnames\\(|^rownames\\(",
    "^View\\(",
    "^setwd\\(",
    "saveRDS\\(|^save\\(|write\\.csv|ggsave\\(|savehistory|pdf\\(|png\\(|svg\\(",
    "Idents\\s*(<-|\\(.*=)|RenameIdents|AddMetaData",
    "subset\\(|WhichCells\\(",
    "merge\\(.*Seurat|IntegrateData|FindIntegrationAnchors",
    "knnreg|CreateDimReducObject",
    "^function\\s*\\(|<-\\s*function\\s*\\(",
    "dir\\.create\\(",
    "monocle|pseudotime|slingshot|learn_graph",
    "scTensor|CellChat|CellPhoneDB|LIANA",
    "rsrecovr|recovr\\(",
    "^#",
    "."
  ),
  label = c(
    "Package loading","Object loading","Seurat object creation","BPCells matrix",
    "File reading","Preprocessing","PCA","Batch correction","UMAP/embedding",
    "Clustering","QC","Marker analysis","UMAP visualization","Feature visualization",
    "Dot/Violin plot","Heatmap","Custom ggplot","Graph visualization","Data inspection",
    "Data viewer","Directory change","Saving","Cell annotation","Cell subsetting",
    "Dataset integration","KNN projection","Function definition","Directory setup",
    "Trajectory analysis","Cell-cell communication","Session recovery","Comment","Other"
  ),
  major = c(
    TRUE,TRUE,TRUE,TRUE,FALSE,
    TRUE,TRUE,TRUE,TRUE,TRUE,
    TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,
    FALSE,FALSE,FALSE,FALSE,
    TRUE,TRUE,TRUE,TRUE,TRUE,
    TRUE,TRUE,FALSE,TRUE,TRUE,
    FALSE,FALSE,FALSE
  ),
  stringsAsFactors = FALSE
)

.sr_categorize <- function(cmd) {
  cmd <- trimws(cmd)
  for (i in seq_len(nrow(.SR_CATS))) {
    if (grepl(.SR_CATS$pat[i], cmd, perl = TRUE)) return(.SR_CATS$name[i])
  }
  "other"
}

.sr_cat_label <- function(cat) {
  idx <- match(cat, .SR_CATS$name)
  if (!is.na(idx)) .SR_CATS$label[idx] else "Other"
}

.sr_cat_major <- function(cat) {
  idx <- match(cat, .SR_CATS$name)
  ifelse(!is.na(idx), .SR_CATS$major[idx], FALSE)
}

# ==============================================================
# Parameter extraction
# ==============================================================

.sr_extract_params <- function(cmd, cat) {
  tryCatch(switch(cat,
    batch = {
      vars <- regmatches(cmd, gregexpr('"[A-Za-z_][^"]*"', cmd))[[1]]
      vars <- gsub('"', '', vars[!vars %in% c("sketch","RNA","SCT")])
      if (length(vars)) paste("vars:", paste(vars, collapse = " + ")) else ""
    },
    clustering = {
      res <- sub("resolution\\s*=\\s*", "res=",
                 .sr_re1(cmd, "resolution\\s*=\\s*[0-9.]+"))
      alg <- if (grepl("algorithm\\s*=\\s*2", cmd)) "Leiden"
             else if (grepl("algorithm\\s*=\\s*1", cmd)) "Louvain" else ""
      paste(c(res, alg)[nchar(c(res, alg)) > 0], collapse = ", ")
    },
    umap = {
      dims <- .sr_re1(cmd, "dims\\s*=\\s*1:\\d+")
      nb   <- .sr_re1(cmd, "n\\.neighbors\\s*=\\s*\\d+")
      paste(c(dims, nb)[nchar(c(dims, nb)) > 0], collapse = ", ")
    },
    preproc = {
      nf <- .sr_re1(cmd, "nfeatures\\s*=\\s*\\d+")
      nc <- .sr_re1(cmd, "ncells\\s*=\\s*\\d+")
      if (nchar(nf)) sub("nfeatures\\s*=\\s*", "nfeatures=", nf)
      else if (nchar(nc)) {
        n <- as.numeric(gsub("[^0-9]", "", nc))
        paste0("ncells=", prettyNum(n, big.mark = ","))
      } else ""
    },
    setwd = {
      p <- .sr_re1(cmd, '"[^"]*"')
      gsub('"', '', p)
    },
    markers = {
      lfc  <- sub("logfc\\.threshold\\s*=\\s*", "logFC>",
                  .sr_re1(cmd, "logfc\\.threshold\\s*=\\s*[0-9.]+"))
      test <- sub('test\\.use\\s*=\\s*"', '',
                  gsub('"', '', .sr_re1(cmd, 'test\\.use\\s*=\\s*"[^"]+"')))
      paste(c(lfc, test)[nchar(c(lfc, test)) > 0], collapse = ", ")
    },
    ""
  ), error = function(e) "")
}

.sr_re1 <- function(s, pat) {
  m <- regmatches(s, regexpr(pat, s, perl = TRUE))
  if (length(m)) m else ""
}

# ==============================================================
# Phase detection
# ==============================================================

.COMPAT_GROUPS <- list(
  c("packages","recover"),
  c("load_obj","create_seurat","open_bpcells","read_file"),
  c("pca","batch","umap"),
  c("dim_viz","feat_viz","dot_vln","heatmap","ggplot","graph_viz"),
  c("table_insp","view","qc"),
  c("dir_setup","save_out")
)

.sr_compatible <- function(a, b) {
  any(vapply(.COMPAT_GROUPS, function(g) a %in% g && b %in% g, logical(1)))
}

.sr_detect_phases <- function(df) {
  SKIP_CATS <- c("comment","other","view","recover","table_insp","dir_setup")
  phases <- list()
  cur    <- list(cats = character(0), rows = integer(0),
                 start_t = df$timestamp[1], wd = df$wd[1])

  for (i in seq_len(nrow(df))) {
    cat_i    <- df$cat[i]
    boundary <- FALSE

    if (cat_i == "setwd") {
      boundary <- TRUE
    } else if (.sr_cat_major(cat_i) && length(cur$cats) > 0) {
      last_major <- utils::tail(cur$cats[.sr_cat_major(cur$cats)], 1)
      if (length(last_major) > 0 && last_major != cat_i &&
          !.sr_compatible(last_major, cat_i)) {
        boundary <- TRUE
      }
    }

    if (boundary && length(cur$rows) > 0) {
      phases <- c(phases, list(cur))
      cur <- list(cats = character(0), rows = integer(0),
                  start_t = df$timestamp[i], wd = df$wd[i])
    }

    cur$cats <- c(cur$cats, cat_i)
    cur$rows <- c(cur$rows, i)
  }
  if (length(cur$rows) > 0) phases <- c(phases, list(cur))
  phases
}

# ==============================================================
# Phase titles
# ==============================================================

.SR_PHASE_TITLES <- c(
  packages    = "Library Setup",
  load_obj    = "Data Loading",
  create_seurat = "Seurat Object Creation",
  open_bpcells  = "BPCells Data Loading",
  read_file   = "File Reading",
  preproc     = "Preprocessing",
  pca         = "PCA",
  batch       = "Batch Correction",
  umap        = "Dimensionality Reduction",
  clustering  = "Clustering",
  qc          = "Quality Control",
  markers     = "Marker Analysis",
  dim_viz     = "UMAP Visualization",
  feat_viz    = "Feature Visualization",
  dot_vln     = "Dot/Violin Plots",
  heatmap     = "Heatmap",
  ggplot      = "Custom Visualization",
  graph_viz   = "Graph Visualization",
  setwd       = "Project / Directory Change",
  save_out    = "Saving Output",
  annot       = "Cell Annotation",
  subset_op   = "Cell Subsetting",
  integrate   = "Dataset Integration",
  knn_proj    = "KNN Projection",
  func_def    = "Function Definition",
  trajectory  = "Trajectory Analysis",
  ccc         = "Cell-Cell Communication",
  table_insp  = "Data Inspection"
)

.sr_phase_title <- function(cats) {
  useful <- cats[!cats %in% c("comment","other","view","recover",
                               "table_insp","dir_setup","read_file")]
  if (!length(useful)) return("Miscellaneous")
  dom   <- names(sort(table(useful), decreasing = TRUE))[1]
  title <- .SR_PHASE_TITLES[dom]
  if (!is.na(title)) title else "Analysis"
}

# ==============================================================
# Bullet generators per category
# ==============================================================

.sr_bullets_for_phase <- function(ph_rows) {
  bullets <- character(0)
  by_cat  <- split(ph_rows, ph_rows$cat)

  for (cat_name in names(by_cat)) {
    rows    <- by_cat[[cat_name]]
    ok_mask <- is.na(rows$status) | rows$status == "OK"
    cmds    <- rows$command[ok_mask]
    if (!length(cmds)) next

    b <- switch(cat_name,
      packages = {
        pkgs <- unique(gsub('"', '', unlist(
          regmatches(rows$command, gregexpr('"[A-Za-z][^"]*"', rows$command))
        )))
        pkgs <- pkgs[!pkgs %in% c("")]
        if (length(pkgs))
          paste0("Loaded packages: ",
                 paste(utils::head(pkgs, 15), collapse = ", "),
                 if (length(pkgs) > 15) ", ..." else "")
        else NULL
      },
      load_obj = {
        fnames <- gsub('"', '', .sr_re1(cmds, '"[^"]+"'))
        if (any(nchar(fnames) > 0))
          paste0("Loaded: ", paste(basename(fnames[nchar(fnames) > 0]),
                                   collapse = ", "))
        else paste0("Loaded R objects (", length(cmds), ")")
      },
      open_bpcells = {
        d <- gsub('"', '', .sr_re1(cmds[1], '"[^"]+"'))
        paste0("Opened on-disk BPCells matrix: `",
               if (nchar(d)) d else "?", "`")
      },
      create_seurat = {
        if (any(grepl("CreateSeuratObject", cmds)))
          "Created Seurat object"
        else if (any(grepl("CreateAssay5Object", cmds)))
          "Created Assay5 object (BPCells-backed)"
        else NULL
      },
      read_file = {
        f <- gsub('"', '', .sr_re1(cmds[1], '"[^"]*\\.csv"'))
        if (nchar(f)) paste0("Read CSV: `", basename(f), "`") else NULL
      },
      preproc = {
        steps <- character(0)
        if (any(grepl("NormalizeData", cmds)))
          steps <- c(steps, "normalized")
        if (any(grepl("FindVariableFeatures", cmds))) {
          nf <- gsub("[^0-9]", "",
                     .sr_re1(cmds[grep("FindVariableFeatures", cmds)[1]],
                             "nfeatures\\s*=\\s*\\d+"))
          steps <- c(steps, paste0("found variable features",
                     if (nchar(nf)) paste0(" (n=", nf, ")") else ""))
        }
        if (any(grepl("ScaleData", cmds))) steps <- c(steps, "scaled data")
        if (any(grepl("SketchData", cmds))) {
          nc <- gsub("[^0-9]", "",
                     .sr_re1(cmds[grep("SketchData", cmds)[1]],
                             "ncells\\s*=\\s*\\d+"))
          steps <- c(steps, paste0("sketched ",
                     if (nchar(nc)) prettyNum(as.numeric(nc), big.mark = ",")
                     else "?", " cells"))
        }
        if (length(steps))
          paste0("Preprocessed: ", paste(steps, collapse = " → "))
        else "Preprocessing steps"
      },
      pca = {
        steps <- character(0)
        if (any(grepl("RunPCA",        cmds))) steps <- c(steps, "ran PCA")
        if (any(grepl("Determine_PCs", cmds))) steps <- c(steps, "selected PCs")
        if (any(grepl("ElbowPlot",     cmds))) steps <- c(steps, "elbow plot")
        if (length(steps)) paste0("PCA: ", paste(steps, collapse = ", "))
        else "PCA"
      },
      batch = {
        var_sets <- lapply(cmds, function(cmd) {
          v <- gsub('"', '', unlist(
            regmatches(cmd, gregexpr('"[A-Za-z_][^"]*"', cmd))
          ))
          v[!v %in% c("sketch","RNA","SCT","data")]
        })
        var_sets <- unique(lapply(var_sets, function(v) v[nchar(v) > 0]))
        if (length(var_sets)) {
          paste(vapply(var_sets, function(v)
            paste0("Harmony correction: ", paste(v, collapse = " + ")),
            character(1)), collapse = "; ")
        } else "Harmony batch correction"
      },
      umap = {
        dims  <- .sr_re1(cmds[1], "dims\\s*=\\s*1:\\d+")
        nb    <- .sr_re1(cmds[1], "n\\.neighbors\\s*=\\s*\\d+")
        extra <- paste(c(dims, nb)[nchar(c(dims, nb)) > 0], collapse = ", ")
        paste0("UMAP", if (nchar(extra)) paste0(" (", extra, ")") else "")
      },
      clustering = {
        parts <- character(0)
        if (any(grepl("FindNeighbors", cmds))) {
          dims <- .sr_re1(cmds[grep("FindNeighbors", cmds)[1]],
                          "dims\\s*=\\s*1:\\d+")
          parts <- c(parts, paste0("built KNN graph",
                     if (nchar(dims)) paste0(" (", dims, ")") else ""))
        }
        if (any(grepl("FindClusters", cmds))) {
          cmd1 <- cmds[grep("FindClusters", cmds)[1]]
          res  <- sub("resolution\\s*=\\s*", "res=",
                      .sr_re1(cmd1, "resolution\\s*=\\s*[0-9.]+"))
          alg  <- if (grepl("algorithm\\s*=\\s*2", cmd1)) " [Leiden]"
                  else if (grepl("algorithm\\s*=\\s*1", cmd1)) " [Louvain]" else ""
          parts <- c(parts, paste0("clustered",
                     if (nchar(res)) paste0(" (", res, alg, ")") else ""))
        }
        paste(parts, collapse = "; ")
      },
      qc = {
        subs <- character(0)
        if (any(grepl("PercentageFeatureSet", cmds)))
          subs <- c(subs, "computed mito %")
        if (any(grepl("table\\(", cmds))) {
          tc <- cmds[grepl("table\\(", cmds)]
          if (any(grepl("mito|Used|Doublet", tc, ignore.case = TRUE)))
            subs <- c(subs, "cross-tabulated QC flags")
          else subs <- c(subs, "tabulated counts")
        }
        if (any(grepl("range\\(", cmds)))
          subs <- c(subs, "checked value ranges")
        if (any(grepl("Doublet", cmds, ignore.case = TRUE)))
          subs <- c(subs, "inspected doublet scores")
        if (!length(subs)) subs <- "QC checks"
        paste(subs, collapse = "; ")
      },
      markers = {
        p <- character(0)
        if (any(grepl("FindAllMarkers", cmds))) {
          lfc <- sub("logfc\\.threshold\\s*=\\s*", "logFC>",
                     .sr_re1(cmds[grep("FindAllMarkers", cmds)[1]],
                             "logfc\\.threshold\\s*=\\s*[0-9.]+"))
          p <- c(p, paste0("FindAllMarkers",
                 if (nchar(lfc)) paste0(" (", lfc, ")") else ""))
        }
        if (any(grepl("top_n\\(", cmds))) {
          n <- gsub("[^0-9]", "",
                    .sr_re1(cmds[grep("top_n", cmds)[1]], "n\\s*=\\s*\\d+"))
          p <- c(p, paste0("extracted top", n, " markers by log2FC"))
        }
        if (length(p)) paste(p, collapse = "; ") else "Marker analysis"
      },
      dim_viz = {
        gb <- unique(unlist(regmatches(cmds,
               gregexpr('group\\.by\\s*=\\s*"[^"]*"', cmds))))
        gb <- gsub('group\\.by\\s*=\\s*"|"', '', gb)
        n  <- length(cmds)
        if (length(gb))
          paste0("DimPlot by: ", paste(gb, collapse = ", "),
                 sprintf(" (%d calls)", n))
        else sprintf("Generated %d DimPlot(s)", n)
      },
      feat_viz = {
        feat <- unique(unlist(
          regmatches(cmds, gregexpr('"[A-Z][A-Z0-9a-z]+[0-9]?"', cmds))
        ))
        feat <- gsub('"', '', feat)
        feat <- feat[nchar(feat) > 1 &
                     !feat %in% c("RNA","sketch","TRUE","FALSE","order")]
        if (length(feat))
          paste0("FeaturePlot: ",
                 paste(utils::head(feat, 10), collapse = ", "),
                 if (length(feat) > 10) ", ..." else "")
        else paste0("FeaturePlot (", length(cmds), " calls)")
      },
      dot_vln = {
        if (any(grepl("DotPlot", cmds))) "DotPlot for marker genes"
        else paste0("VlnPlot (", length(cmds), " calls)")
      },
      ggplot = {
        geoms <- unique(unlist(regmatches(cmds, gregexpr("geom_[a-z_]+", cmds))))
        if (length(geoms))
          paste0("Custom ggplot (", paste(geoms, collapse = ", "), ")")
        else paste0("Custom ggplot (", length(unique(cmds)), " unique figures)")
      },
      graph_viz = "Graph visualization (ggraph/igraph)",
      setwd = {
        d <- gsub('"', '', .sr_re1(cmds[length(cmds)], '"[^"]*"'))
        if (nchar(d)) paste0("Switched to: `", d, "`") else NULL
      },
      annot     = paste0("Cell annotation / identity operations (", length(cmds), ")"),
      subset_op = {
        n <- length(unique(cmds))
        paste0("Subsetted cells (", n, " operation", if (n > 1) "s" else "", ")")
      },
      save_out = {
        fnames <- unlist(regmatches(cmds,
                    gregexpr('"[^"]+\\.[a-zA-Z]{1,5}"', cmds)))
        fnames <- basename(gsub('"', '', fnames))
        fnames <- fnames[nchar(fnames) > 0]
        if (length(fnames))
          paste0("Saved: ", paste(fnames, collapse = ", "))
        else "Saved output"
      },
      integrate  = paste0("Merged/integrated Seurat objects (", length(cmds), ")"),
      knn_proj   = "KNN projection onto reference UMAP",
      func_def   = "Defined custom function(s)",
      trajectory = paste0("Trajectory analysis (", length(cmds), " calls)"),
      ccc        = "Cell-cell communication analysis",
      NULL
    )

    if (!is.null(b) && nchar(b) > 0) bullets <- c(bullets, b)
  }

  if (!length(bullets)) bullets <- "Miscellaneous commands"
  bullets
}

# ==============================================================
# Markdown generator
# ==============================================================

.sr_gen_markdown <- function(df, phases) {
  n_cmds   <- nrow(df)
  n_errors <- sum(!is.na(df$status) & df$status == "ERR")
  n_plots  <- attr(df, "n_plots") %||% 0L

  ts_vals  <- df$timestamp[!is.na(df$timestamp)]
  date_str <- if (length(ts_vals)) format(min(ts_vals), "%Y-%m-%d") else
                format(Sys.Date())
  duration <- if (length(ts_vals) >= 2) {
    d <- as.numeric(difftime(max(ts_vals), min(ts_vals), units = "mins"))
    sprintf("~%d min", round(d))
  } else "unknown"

  wds      <- unique(df$wd[!is.na(df$wd) & df$wd != "?"])
  all_cats <- unique(df$cat)
  overview <- character(0)

  if (any(c("load_obj","create_seurat","open_bpcells") %in% all_cats))
    overview <- c(overview, "loaded data")
  if (any(c("preproc","pca","batch","umap") %in% all_cats))
    overview <- c(overview, "ran preprocessing and dimensionality reduction")
  if ("clustering" %in% all_cats)  overview <- c(overview, "clustered cells")
  if ("qc"         %in% all_cats)  overview <- c(overview, "performed QC")
  if ("markers"    %in% all_cats)  overview <- c(overview, "ran marker analysis")
  if (any(c("dim_viz","feat_viz","ggplot","dot_vln") %in% all_cats))
    overview <- c(overview, "generated visualizations")

  out <- c(
    paste0("# R Session Summary — ", date_str),
    "",
    sprintf("**Duration:** %s | **Commands:** %d | **Errors:** %d | **Plots captured:** %d",
            duration, n_cmds, n_errors, n_plots),
    ""
  )

  if (length(wds)) {
    out <- c(out, "**Working directories:**")
    for (w in wds) out <- c(out, paste0("- `", w, "`"))
    out <- c(out, "")
  }

  if (length(overview)) {
    out <- c(out, "## Overview", "",
             paste0("This session ", paste(overview, collapse = ", "), "."),
             "")
  }

  out <- c(out, "## Analysis Phases", "")

  SKIP_CATS <- c("comment","other","view","recover")
  phase_num <- 0L
  for (ph in phases) {
    ph_rows <- df[ph$rows, ]
    useful  <- ph_rows[!ph_rows$cat %in% SKIP_CATS, ]
    if (nrow(useful) == 0) next

    phase_num <- phase_num + 1L
    title     <- .sr_phase_title(useful$cat)

    ts_ph  <- ph_rows$timestamp[!is.na(ph_rows$timestamp)]
    t_hdr  <- if (length(ts_ph)) format(ts_ph[1], "(%H:%M)") else ""

    out <- c(out, sprintf("### Phase %d: %s %s", phase_num, title, t_hdr), "")
    for (b in .sr_bullets_for_phase(useful)) out <- c(out, paste0("- ", b))

    err_rows <- ph_rows[!is.na(ph_rows$status) & ph_rows$status == "ERR", ]
    if (nrow(err_rows) > 0) {
      out <- c(out, "", "  **Errors in this phase:**")
      for (j in seq_len(min(nrow(err_rows), 5))) {
        s <- substr(err_rows$command[j], 1, 90)
        if (nchar(err_rows$command[j]) > 90) s <- paste0(s, "...")
        out <- c(out, paste0("  - `", s, "`"))
      }
    }
    out <- c(out, "")
  }

  pf <- attr(df, "plot_files")
  if (length(pf) > 0) {
    out <- c(out, "## Plots Captured", "")
    for (i in seq_along(pf)) out <- c(out, sprintf("%d. `%s`", i, basename(pf[i])))
    out <- c(out, "")
  }

  if (n_errors > 0) {
    err_all <- df[!is.na(df$status) & df$status == "ERR", ]
    out <- c(out, "## All Errors", "")
    for (i in seq_len(min(nrow(err_all), 15))) {
      s <- substr(err_all$command[i], 1, 100)
      if (nchar(err_all$command[i]) > 100) s <- paste0(s, "...")
      out <- c(out, paste0("- `", s, "`"))
    }
    out <- c(out, "")
  }

  pkg_cmds <- df$command[df$cat == "packages"]
  if (length(pkg_cmds)) {
    pkgs <- unique(gsub('"', '', unlist(
      regmatches(pkg_cmds, gregexpr('"[A-Za-z][^"]*"', pkg_cmds))
    )))
    pkgs <- pkgs[nchar(pkgs) > 0]
    if (length(pkgs)) {
      out <- c(out, "## Packages Used", "",
               paste(pkgs, collapse = ", "), "")
    }
  }

  out <- c(out, "---",
           paste0("*Generated by scSidekick::summarize_r_session() — ",
                  format(Sys.time(), "%Y-%m-%d %H:%M"), "*"))
  out
}

# ==============================================================
# Clean R script generator
# ==============================================================

.sr_gen_clean_r <- function(df, phases) {
  SKIP_CATS <- c("comment","other","view","recover","table_insp","dir_setup")

  out <- c(
    "# ============================================================",
    paste0("# Reconstructed R script — session: ",
           attr(df, "session_id") %||% "?"),
    paste0("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "# ============================================================",
    ""
  )

  phase_num <- 0L
  for (ph in phases) {
    ph_rows <- df[ph$rows, ]
    useful  <- ph_rows[!ph_rows$cat %in% SKIP_CATS, ]
    useful  <- useful[is.na(useful$status) | useful$status == "OK", ]
    if (nrow(useful) == 0) next

    phase_num <- phase_num + 1L
    title <- .sr_phase_title(useful$cat)
    out <- c(out,
             paste0("# ---- Phase ", phase_num, ": ", title, " ----"),
             "")

    cmds <- useful$command
    keep <- c(TRUE, cmds[-1] != cmds[-length(cmds)])  # remove consecutive dups
    cmds <- cmds[keep]

    for (cmd in cmds) {
      cmd <- gsub("\\\\n", "\n", cmd)
      out <- c(out, cmd, "")
    }
  }
  out
}
