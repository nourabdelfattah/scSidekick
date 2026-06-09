# =============================================================================
# scSidekick feature annotation
#
# AnnotateFeatures  - adds genomic metadata (chr, start, end, biotype …) to
#                     feature-level metadata in a Seurat object.
# GetFeatures       - filter genes by any annotation column; returns a
#                     character vector usable anywhere downstream.
# GetGeneOrder      - returns an inferCNV-ready 3-column data.frame
#                     (chr | start | end) sorted genomically.
#
# Annotation source hierarchy (automatic, with clear messaging):
#   1. Local cache  (~/.cache/R/scSidekick/<species>_annotation.rds)
#   2. biomaRt      (Ensembl REST API — fast when available)
#   3. AnnotationHub / EnsDb  (downloads once, Bioconductor cache)
#   4. Precomputed table downloaded from the scSidekick GitHub release
#   5. User-supplied file  (annotation_file = "path/to/file.tsv")
# =============================================================================

# ── Internal species registry ─────────────────────────────────────────────────

.nk_species_info <- list(
  mouse = list(
    mart    = "mmusculus_gene_ensembl",
    ah_q    = c("EnsDb", "Mus musculus"),
    precomp = "mouse_annotation.rds"
  ),
  human = list(
    mart    = "hsapiens_gene_ensembl",
    ah_q    = c("EnsDb", "Homo sapiens"),
    precomp = "human_annotation.rds"
  ),
  zebrafish = list(
    mart    = "drerio_gene_ensembl",
    ah_q    = c("EnsDb", "Danio rerio"),
    precomp = "zebrafish_annotation.rds"
  )
)

.nk_precomp_url <- paste0(
  "https://github.com/nourabdelfattah/scSidekick/",
  "releases/download/annotation-v1/"
)

# ── Cache helpers ─────────────────────────────────────────────────────────────

.nk_cache_dir  <- function() tools::R_user_dir("scSidekick", "cache")

.nk_ann_cache_path <- function(species)
  file.path(.nk_cache_dir(), paste0(species, "_annotation.rds"))

.nk_ann_load_cache <- function(species) {
  path <- .nk_ann_cache_path(species)
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

.nk_ann_save_cache <- function(species, ann) {
  dir.create(.nk_cache_dir(), recursive = TRUE, showWarnings = FALSE)
  tryCatch(saveRDS(ann, .nk_ann_cache_path(species)), error = function(e) NULL)
  invisible(NULL)
}

# ── Annotation standardiser ───────────────────────────────────────────────────
# Converts any source's column names → our standard nk_* names.

.nk_standardise_annotation <- function(df, col_map) {
  for (src in names(col_map)) {
    tgt <- col_map[[src]]
    if (src %in% colnames(df)) {
      colnames(df)[colnames(df) == src] <- tgt
    }
  }
  keep <- c(grep("^nk_", colnames(df), value = TRUE))
  df[, keep, drop = FALSE]
}

# ── Source 2: biomaRt ─────────────────────────────────────────────────────────

.nk_ann_from_biomart <- function(species) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) return(NULL)
  info <- .nk_species_info[[species]]

  mart <- tryCatch(
    biomaRt::useMart("ensembl", dataset = info$mart),
    error = function(e) NULL
  )
  if (is.null(mart)) return(NULL)

  ann <- tryCatch(
    biomaRt::getBM(
      attributes = c(
        "external_gene_name", "ensembl_gene_id",
        "chromosome_name",    "start_position",
        "end_position",       "strand",
        "gene_biotype",       "description"
      ),
      mart = mart
    ),
    error = function(e) NULL
  )
  if (is.null(ann) || nrow(ann) == 0) return(NULL)

  ann <- ann[nzchar(ann$external_gene_name), ]
  # Prefer X over Y for PAR genes present on both chromosomes
  chr_rank <- ifelse(ann$chromosome_name == "X", 1L,
              ifelse(ann$chromosome_name == "Y", 2L, 0L))
  ann <- ann[order(chr_rank), ]
  ann <- ann[!duplicated(ann$external_gene_name), ]
  rownames(ann) <- ann$external_gene_name

  .nk_standardise_annotation(ann, c(
    ensembl_gene_id  = "nk_ensembl_id",
    chromosome_name  = "nk_chr",
    start_position   = "nk_start",
    end_position     = "nk_end",
    strand           = "nk_strand",
    gene_biotype     = "nk_biotype",
    description      = "nk_description"
  ))
}

# ── Source 3: AnnotationHub / EnsDb ──────────────────────────────────────────

.nk_ann_from_ensdb <- function(species) {
  if (!requireNamespace("AnnotationHub",  quietly = TRUE)) return(NULL)
  if (!requireNamespace("ensembldb",      quietly = TRUE)) return(NULL)
  info <- .nk_species_info[[species]]

  ah <- tryCatch(AnnotationHub::AnnotationHub(), error = function(e) NULL)
  if (is.null(ah)) return(NULL)

  q <- tryCatch(AnnotationHub::query(ah, info$ah_q), error = function(e) NULL)
  if (is.null(q) || length(q) == 0) return(NULL)

  ensdb <- tryCatch(ah[[tail(names(q), 1)]], error = function(e) NULL)
  if (is.null(ensdb)) return(NULL)

  ann <- tryCatch(
    ensembldb::genes(ensdb, return.type = "data.frame"),
    error = function(e) NULL
  )
  if (is.null(ann) || nrow(ann) == 0) return(NULL)

  ann <- ann[nzchar(ann$gene_name), ]
  # Prefer X over Y for PAR genes present on both chromosomes
  chr_rank <- ifelse(ann$seq_name == "X", 1L,
              ifelse(ann$seq_name == "Y", 2L, 0L))
  ann <- ann[order(chr_rank), ]
  ann <- ann[!duplicated(ann$gene_name), ]
  rownames(ann) <- ann$gene_name

  .nk_standardise_annotation(ann, c(
    gene_id       = "nk_ensembl_id",
    seq_name      = "nk_chr",
    gene_seq_start = "nk_start",
    gene_seq_end  = "nk_end",
    seq_strand    = "nk_strand",
    gene_biotype  = "nk_biotype",
    description   = "nk_description"
  ))
}

# ── Source 4: download precomputed table ──────────────────────────────────────

.nk_ann_download_precomputed <- function(species) {
  info <- .nk_species_info[[species]]
  url  <- paste0(.nk_precomp_url, info$precomp)
  tmp  <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  ok <- tryCatch({
    utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
    TRUE
  }, error = function(e) FALSE)

  if (!ok) return(NULL)
  tryCatch(readRDS(tmp), error = function(e) NULL)
}

# ── Source 5: user-supplied file ──────────────────────────────────────────────
# Accepts any TSV/CSV where rownames are gene symbols and columns include
# chromosome, start, end (flexible column name aliases accepted).

.nk_ann_from_file <- function(path) {
  ann <- tryCatch(
    read.delim(path, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Cannot read annotation file: ", path, "\n", e$message)
  )

  aliases <- list(
    nk_chr         = c("chr", "chromosome", "chromosome_name", "seqnames", "seq_name"),
    nk_start       = c("start", "start_position", "gene_seq_start"),
    nk_end         = c("end",   "end_position",   "gene_seq_end"),
    nk_strand      = c("strand", "seq_strand"),
    nk_biotype     = c("gene_biotype", "biotype", "gene_type"),
    nk_ensembl_id  = c("ensembl_gene_id", "gene_id", "ensembl_id"),
    nk_description = c("description")
  )

  result <- data.frame(row.names = rownames(ann), stringsAsFactors = FALSE)
  for (nk_col in names(aliases)) {
    hit <- intersect(aliases[[nk_col]], colnames(ann))
    if (length(hit) > 0) result[[nk_col]] <- ann[[hit[1]]]
  }
  result
}

# ── Master fetch (tries sources in order, caches success) ─────────────────────

.nk_fetch_annotation <- function(species) {
  # 1. Cache
  ann <- .nk_ann_load_cache(species)
  if (!is.null(ann)) {
    message("scSidekick: Loaded ", species, " annotation from cache.")
    return(ann)
  }

  # 2. biomaRt
  message("scSidekick: Querying biomaRt for ", species, " gene annotations",
          " (this runs once and is then cached) ...")
  ann <- .nk_ann_from_biomart(species)
  if (!is.null(ann)) {
    message("scSidekick: biomaRt succeeded. Caching to ", .nk_ann_cache_path(species))
    .nk_ann_save_cache(species, ann)
    return(ann)
  }
  message("scSidekick: biomaRt unavailable. Trying AnnotationHub ...")

  # 3. AnnotationHub / EnsDb
  ann <- .nk_ann_from_ensdb(species)
  if (!is.null(ann)) {
    message("scSidekick: AnnotationHub succeeded. Caching.")
    .nk_ann_save_cache(species, ann)
    return(ann)
  }
  message("scSidekick: AnnotationHub unavailable. Downloading precomputed table ...")

  # 4. Download precomputed
  ann <- .nk_ann_download_precomputed(species)
  if (!is.null(ann)) {
    message("scSidekick: Downloaded precomputed table. Caching.")
    .nk_ann_save_cache(species, ann)
    return(ann)
  }

  message(
    "scSidekick: All automatic annotation sources failed.\n",
    "  Supply your own file: AnnotateFeatures(..., annotation_file = 'path/to/file.tsv')\n",
    "  Expected columns: chr/chromosome, start, end  (gene symbol as rowname)."
  )
  NULL
}

# ── AnnotateFeatures ──────────────────────────────────────────────────────────

#' Annotate Seurat features with genomic metadata
#'
#' Adds chromosome location, gene biotype, Ensembl ID, and other genomic
#' metadata as columns to the feature-level metadata of a Seurat object.
#' Annotation is fetched from (in order): local cache → biomaRt →
#' AnnotationHub/EnsDb → precomputed download → user-supplied file.
#' Successful results are cached in `~/.cache/R/scSidekick/` so the network
#' call is only made once per species.
#'
#' @param seurat_object A Seurat object.
#' @param species Character. One of `"mouse"`, `"human"`, `"zebrafish"`.
#' @param assay Character. Assay to annotate. Defaults to `DefaultAssay()`.
#' @param annotation_file Character or `NULL`. Path to a user-supplied TSV/CSV
#'   (gene symbol as rowname; any recognized chromosome/start/end column names
#'   are accepted). Bypasses all automatic sources when supplied.
#' @param force Logical. Re-annotate even if `nk_chr` column is already
#'   present. Default `FALSE`.
#'
#' @return The Seurat object with new columns in its feature metadata:
#'   `nk_chr`, `nk_start`, `nk_end`, `nk_strand`, `nk_biotype`,
#'   `nk_ensembl_id`, `nk_description`.
#' @seealso [GetFeatures()], [GetGeneOrder()]
#' @export
AnnotateFeatures <- function(seurat_object,
                              species,
                              assay           = NULL,
                              annotation_file = NULL,
                              force           = FALSE) {

  species <- match.arg(species, names(.nk_species_info))
  assay   <- assay %||% Seurat::DefaultAssay(seurat_object)

  existing <- seurat_object@misc$nk_annotations[[assay]]

  if (!force && !is.null(existing) && "nk_chr" %in% colnames(existing)) {
    message("scSidekick: Features already annotated. Use force = TRUE to re-annotate.")
    return(seurat_object)
  }

  ann <- if (!is.null(annotation_file)) {
    message("scSidekick: Using user-supplied annotation file.")
    .nk_ann_from_file(annotation_file)
  } else {
    .nk_fetch_annotation(species)
  }

  if (is.null(ann)) return(seurat_object)

  genes     <- rownames(seurat_object@assays[[assay]])
  feat_meta <- data.frame(row.names = genes, stringsAsFactors = FALSE)
  common    <- intersect(genes, rownames(ann))
  nk_cols   <- grep("^nk_", colnames(ann), value = TRUE)

  for (col in nk_cols) {
    feat_meta[[col]] <- NA
    feat_meta[common, col] <- ann[common, col]
  }

  if (is.null(seurat_object@misc$nk_annotations))
    seurat_object@misc$nk_annotations <- list()
  seurat_object@misc$nk_annotations[[assay]] <- feat_meta
  seurat_object@misc$nk_annotation_species   <- species

  n_ann <- sum(!is.na(feat_meta$nk_chr))
  message("scSidekick: Annotated ", n_ann, " / ", nrow(feat_meta), " features with genomic coordinates.")
  seurat_object
}

# ── GetFeatures ───────────────────────────────────────────────────────────────

#' Filter features by genomic annotation
#'
#' Returns a character vector of gene names matching any combination of
#' chromosome, biotype, or name pattern filters. Requires [AnnotateFeatures()]
#' to have been run first.
#'
#' @param seurat_object A Seurat object with feature annotations.
#' @param chromosome Numeric or character vector of chromosomes to keep, e.g.
#'   `1:19`, `c("X","Y")`, `"MT"`. Both `"1"` and `"chr1"` are accepted.
#'   `NULL` keeps all chromosomes.
#' @param gene_biotype Character vector of Ensembl biotype(s) to keep, e.g.
#'   `"protein_coding"`, `c("protein_coding","lncRNA")`. `NULL` keeps all.
#' @param pattern Character. Regular expression matched against gene names.
#'   `NULL` skips pattern filtering.
#' @param assay Character. Which assay's feature metadata to use.
#' @param invert Logical. Return genes that do NOT match the filters.
#'   Default `FALSE`.
#'
#' @return A character vector of gene names passing all filters.
#' @seealso [AnnotateFeatures()], [GetGeneOrder()]
#' @export
GetFeatures <- function(seurat_object,
                         chromosome   = NULL,
                         gene_biotype = NULL,
                         pattern      = NULL,
                         assay        = NULL,
                         invert       = FALSE) {

  assay     <- assay %||% Seurat::DefaultAssay(seurat_object)
  feat_meta <- seurat_object@misc$nk_annotations[[assay]]

  if (is.null(feat_meta) || !"nk_chr" %in% colnames(feat_meta))
    stop("No feature annotations found. Run AnnotateFeatures() first.")

  keep <- rep(TRUE, nrow(feat_meta))

  if (!is.null(chromosome)) {
    chrs <- gsub("^chr", "", as.character(chromosome), ignore.case = TRUE)
    keep <- keep & !is.na(feat_meta$nk_chr) & feat_meta$nk_chr %in% chrs
  }

  if (!is.null(gene_biotype)) {
    keep <- keep & !is.na(feat_meta$nk_biotype) & feat_meta$nk_biotype %in% gene_biotype
  }

  if (!is.null(pattern)) {
    keep <- keep & grepl(pattern, rownames(feat_meta))
  }

  if (invert) keep <- !keep
  rownames(feat_meta)[keep]
}

# ── GetGeneOrder ──────────────────────────────────────────────────────────────

#' Get genomically-ordered gene position table
#'
#' Returns a data.frame of gene positions sorted by chromosome then start
#' position, suitable for direct use as the `gene_order_file` argument in
#' [infercnv::CreateInfercnvObject()].
#'
#' @param seurat_object A Seurat object with feature annotations from
#'   [AnnotateFeatures()].
#' @param chromosomes Numeric or character vector of chromosomes to include,
#'   e.g. `1:19` for autosomes. `NULL` includes all annotated chromosomes.
#' @param gene_biotype Character. Biotype filter applied before building the
#'   table. Default `"protein_coding"`.
#' @param assay Character. Assay to use.
#' @param chr_prefix Character. Prefix prepended to chromosome numbers in the
#'   output (e.g. `"chr"` → `"chr1"`). Default `"chr"`.
#'
#' @return A `data.frame` with columns `chr`, `start`, `end` and gene symbols
#'   as rownames, sorted genomically.
#' @seealso [AnnotateFeatures()], [GetFeatures()]
#' @export
GetGeneOrder <- function(seurat_object,
                          chromosomes  = NULL,
                          gene_biotype = "protein_coding",
                          assay        = NULL,
                          chr_prefix   = "chr") {

  assay     <- assay %||% Seurat::DefaultAssay(seurat_object)
  feat_meta <- seurat_object@misc$nk_annotations[[assay]]

  if (is.null(feat_meta) || !"nk_chr" %in% colnames(feat_meta))
    stop("No feature annotations found. Run AnnotateFeatures() first.")

  # Start with genes that have coordinates
  ann <- feat_meta[!is.na(feat_meta$nk_chr) & !is.na(feat_meta$nk_start), ]

  if (!is.null(chromosomes)) {
    chrs <- gsub("^chr", "", as.character(chromosomes), ignore.case = TRUE)
    ann  <- ann[ann$nk_chr %in% chrs, ]
  }

  if (!is.null(gene_biotype) && "nk_biotype" %in% colnames(ann)) {
    ann <- ann[!is.na(ann$nk_biotype) & ann$nk_biotype %in% gene_biotype, ]
  }

  if (nrow(ann) == 0)
    stop("No genes remain after filtering. Check chromosome and biotype arguments.")

  # Sort chromosomes numerically where possible, non-numeric (X, Y, MT) at end
  chr_levels <- unique(ann$nk_chr)
  chr_num    <- suppressWarnings(as.numeric(chr_levels))
  chr_levels <- chr_levels[order(ifelse(is.na(chr_num), Inf, chr_num))]
  ann$nk_chr <- factor(ann$nk_chr, levels = chr_levels)
  ann        <- ann[order(ann$nk_chr, ann$nk_start), ]

  data.frame(
    chr   = paste0(chr_prefix, as.character(ann$nk_chr)),
    start = as.integer(ann$nk_start),
    end   = as.integer(ann$nk_end),
    row.names        = rownames(ann),
    stringsAsFactors = FALSE
  )
}

# ── ClearAnnotationCache ──────────────────────────────────────────────────────

#' Clear the local feature annotation cache
#'
#' Removes cached annotation files from `~/.cache/R/scSidekick/` so that the
#' next call to [AnnotateFeatures()] re-fetches from source. Useful when you
#' want to update to a newer Ensembl version.
#'
#' @param species Character vector. Species to clear. `NULL` clears all cached
#'   annotations.
#'
#' @return Invisibly returns the paths deleted.
#' @export
ClearAnnotationCache <- function(species = NULL) {
  cache_dir <- .nk_cache_dir()
  if (!dir.exists(cache_dir)) { message("Cache directory does not exist."); return(invisible(character(0))) }

  if (is.null(species)) {
    paths <- list.files(cache_dir, pattern = "_annotation\\.rds$", full.names = TRUE)
  } else {
    paths <- file.path(cache_dir, paste0(species, "_annotation.rds"))
    paths <- paths[file.exists(paths)]
  }

  if (length(paths) == 0) { message("No cached annotations found."); return(invisible(character(0))) }

  unlink(paths)
  message("scSidekick: Cleared ", length(paths), " cached annotation file(s).")
  invisible(paths)
}
