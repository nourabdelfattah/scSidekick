# =============================================================================
# scSidekick - project scaffolding  (generate_directories.R)
#
# GenerateDirectories() - create the standard scSidekick project folder tree
#   (R_Objects / R_Code / Output / Output/<Subset> / per-analysis subfolders)
#   in one call and return every path as a named list, so downstream functions
#   and PrepObject() can be pointed at them without re-pasting paths by hand.
# =============================================================================


#' Create a standard project directory tree
#'
#' @description
#' Builds the conventional scSidekick analysis folder layout under `base_dir`
#' and returns every created path as a named list.  Replaces the usual pile of
#' `paste0()` + `dir.create()` lines at the top of an analysis script.
#'
#' The tree created is:
#' ```
#' base_dir/
#'   R_Objects/                  # saved Seurat objects (.rds)
#'   R_Code/                     # analysis scripts
#'   Output/                     # all figures / tables
#'     <subset>/                 # one analysis branch (e.g. "All_Clusters")
#'       Featuremaps/            # (and any other `output_subdirs`)
#' ```
#'
#' @param base_dir Character. Project root directory. Created if absent.
#' @param subset Character. Name of the analysis branch under `Output/`
#'   (e.g. `"All_Clusters"`, `"Tcells"`). Default `"All_Clusters"`.
#' @param object_name Character or `NULL`. Project label. Stored on the
#'   returned list (and used by downstream functions for file prefixes); does
#'   not affect the folder structure. Default `NULL`.
#' @param resolution Numeric or `NULL`. Clustering resolution. Stored on the
#'   returned list for reference only. Default `NULL`.
#' @param output_subdirs Character vector. Subfolders to create *inside* the
#'   `Output/<subset>/` branch. Default `"Featuremaps"`. Pass `character(0)`
#'   for none, or e.g. `c("Featuremaps", "DotPlots", "Pathways")`.
#' @param extra_top Character vector or `NULL`. Additional top-level folders to
#'   create directly under `base_dir` (beyond `R_Objects`, `R_Code`,
#'   `Output`). Default `NULL`.
#' @param set_wd Logical. Set the working directory to `base_dir`? Default
#'   `FALSE`.
#' @param verbose Logical. Message each folder as it is created? Default
#'   `TRUE`.
#'
#' @return Invisibly, a named list of absolute paths:
#'   \describe{
#'     \item{`base_dir`}{the project root}
#'     \item{`robj_dir`}{`base_dir/R_Objects`}
#'     \item{`rcode_dir`}{`base_dir/R_Code`}
#'     \item{`output_root`}{`base_dir/Output`}
#'     \item{`output_dir`}{`base_dir/Output/<subset>` - the per-analysis branch}
#'     \item{one entry per `output_subdirs`}{keyed by a lower-case, snake-case
#'       version of the folder name, e.g. `featuremaps`}
#'     \item{one entry per `extra_top`}{keyed the same way}
#'     \item{`object_name`, `subset`, `resolution`}{echoed back for convenience}
#'   }
#'   Use `output_dir` directly with [PrepObject()] or any plotting function.
#'
#' @examples
#' \dontrun{
#' dirs <- GenerateDirectories(
#'   base_dir    = "~/Projects/YAP_scRNAseq",
#'   subset      = "All_Clusters",
#'   object_name = "YAP project",
#'   resolution  = 1,
#'   output_subdirs = c("Featuremaps", "DotPlots", "Pathways"))
#'
#' # Wire the paths straight into the rest of the workflow:
#' obj <- PrepObject(obj, variables = "celltype.l1",
#'                   output_dir = dirs$output_dir, object_name = dirs$object_name)
#' saveRDS(obj, file.path(dirs$robj_dir, "obj.rds"))
#' RunGSEA(obj, output_dir = dirs$pathways)
#' }
#'
#' @seealso [PrepObject()], [LoadSamplesRNA()]
#' @export
GenerateDirectories <- function(
    base_dir,
    subset         = "All_Clusters",
    object_name    = NULL,
    resolution     = NULL,
    output_subdirs = "Featuremaps",
    extra_top      = NULL,
    set_wd         = FALSE,
    verbose        = TRUE
) {
  if (missing(base_dir) || !is.character(base_dir) || length(base_dir) != 1L)
    stop("`base_dir` must be a single directory path.")

  base_dir <- path.expand(base_dir)

  # Internal: create one directory (recursive), optionally announce it.
  .mk <- function(path) {
    created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (verbose)
      message(if (created) "  created " else "  exists  ", path)
    path
  }

  # Internal: turn "Featuremaps" / "Dot Plots" into a valid list key "featuremaps"
  .key <- function(x) {
    k <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
    gsub("^_|_$", "", k)
  }

  if (verbose) message("Generating project directories under: ", base_dir)

  # ── Top-level tree ──────────────────────────────────────────────────────────
  paths <- list(
    base_dir    = .mk(base_dir),
    robj_dir    = .mk(file.path(base_dir, "R_Objects")),
    rcode_dir   = .mk(file.path(base_dir, "R_Code")),
    output_root = .mk(file.path(base_dir, "Output"))
  )

  # ── Per-analysis branch: Output/<subset> ────────────────────────────────────
  output_dir <- .mk(file.path(paths$output_root, subset))
  paths$output_dir <- output_dir

  # ── Subfolders inside the analysis branch ───────────────────────────────────
  for (sd in output_subdirs) {
    key <- .key(sd)
    if (nzchar(key)) paths[[key]] <- .mk(file.path(output_dir, sd))
  }

  # ── Extra top-level folders (optional) ──────────────────────────────────────
  for (td in (extra_top %||% character(0))) {
    key <- .key(td)
    if (nzchar(key)) paths[[key]] <- .mk(file.path(base_dir, td))
  }

  # ── Echo back the descriptive fields for convenience ────────────────────────
  paths$object_name <- object_name
  paths$subset      <- subset
  paths$resolution  <- resolution

  if (isTRUE(set_wd)) {
    setwd(base_dir)
    if (verbose) message("Working directory set to: ", base_dir)
  }

  if (verbose) message("Done. ", length(output_subdirs) + 4L,
                       " core folders ready.")

  invisible(paths)
}
