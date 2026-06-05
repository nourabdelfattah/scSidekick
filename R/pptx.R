# =============================================================================
# scSidekick PPTX builder
#
# log_analysis_params()   — reads a Seurat object + writes analysis_params.json
# create_analysis_pptx()  — converts all PDFs in an output folder to a PPTX
#
# Required (Suggests): officer, magick, jsonlite
# =============================================================================

# ---------------------------------------------------------------------------
# Section-grouping rules — first match wins, order matters
# ---------------------------------------------------------------------------
.pptx_section_rules <- list(
  # ── Quality Control ────────────────────────────────────────────────────────
  list(section = "Quality Control",
       pattern = "qualitycontrol|doublet|miQC|cells per sample|filtration|barplot number"),

  # ── Dimensionality Reduction ───────────────────────────────────────────────
  list(section = "Dimensionality Reduction",
       pattern = "harmony|umap before and after|elbow|PCA|RNA UMAPs"),

  # ── Pathway Analysis (GSEA) — must come before Clustering so that GSEA
  #    heatmaps (which contain "heatmap") don't fall into the Clustering rule.
  list(section = "Pathway Analysis",
       pattern = "GSEA|\\bKEGG\\b|Reactome|WikiPathway|Hallmark"),

  # ── CellChat — before Clustering for the same reason ("heatmap" in names).
  #    Covers RunCellChat, CompareCellChat, and RankCellChatPathways outputs.
  list(section = "CellChat",
       pattern = paste0(
         "Number of Interactions|Interaction Weights|Individual Groups|",
         "cell-cell communication|signaling comparison|bubble plot|",
         "dominant sender receiver|communication patterns|CellChat_"
       )),

  # ── Cluster Markers (renamed from Clustering for clarity) ─────────────────
  list(section = "Cluster Markers",
       pattern = "heatmap|markers|cluster"),

  # ── Dotplots ───────────────────────────────────────────────────────────────
  list(section = "Dotplots",
       pattern = "dotplot"),

  # ── Cell Type Annotation ───────────────────────────────────────────────────
  list(section = "Cell Type Annotation",
       pattern = "AutoAssignment|singleR|Assignment Helper|Assignment umap|featureheatmap"),

  # ── All UMAPs ──────────────────────────────────────────────────────────────
  list(section = "All UMAPs",
       pattern = "All umaps|all umap"),

  # ── Composition ────────────────────────────────────────────────────────────
  list(section = "Composition",
       pattern = "percentage|iteration with bar|chord|barplot|rose|trend|alluvial"),

  # ── Feature Maps ───────────────────────────────────────────────────────────
  list(section = "Feature Maps",
       pattern = "featuremap"),

  # ── Other (catch-all) ──────────────────────────────────────────────────────
  list(section = "Other",
       pattern = ".")
)

.pptx_assign_section <- function(filename) {
  bn <- basename(filename)
  for (rule in .pptx_section_rules) {
    if (grepl(rule$pattern, bn, ignore.case = TRUE)) return(rule$section)
  }
  "Other"
}

.pptx_clean_title <- function(filepath, obj_name = "", subset_name = "") {
  bn <- tools::file_path_sans_ext(basename(filepath))
  if (nchar(obj_name) > 0)    bn <- gsub(obj_name,    "", bn, fixed = TRUE)
  if (nchar(subset_name) > 0) bn <- gsub(subset_name, "", bn, fixed = TRUE)
  bn <- gsub("\\bres[0-9\\.]+\\b", "", bn, ignore.case = TRUE)
  bn <- trimws(gsub("\\s+", " ", bn))
  paste0(toupper(substr(bn, 1, 1)), substr(bn, 2, nchar(bn)))
}

.pptx_make_footer <- function(params) {
  parts <- character(0)
  if (!is.null(params$resolution))
    parts <- c(parts, paste0("Resolution: ", params$resolution))
  if (!is.null(params$pcs_used))
    parts <- c(parts, paste0("PCs: ", params$pcs_used))
  if (!is.null(params$n_clusters))
    parts <- c(parts, paste0("Clusters: ", params$n_clusters))
  if (!is.null(params$harmony_vars)) {
    hv    <- paste(unlist(params$harmony_vars), collapse = ", ")
    parts <- c(parts, paste0("Harmony: ", hv))
  }
  if (!is.null(params$n_cells_final))
    parts <- c(parts, paste0("Cells: ", format(params$n_cells_final, big.mark = ",")))
  if (!is.null(params$date))
    parts <- c(parts, paste0("Date: ", params$date))
  paste(parts, collapse = "   |   ")
}

.pptx_pdf_to_png <- function(pdf_path, density = 150) {
  img_path <- tempfile(fileext = ".png")
  tryCatch({
    img <- magick::image_read_pdf(pdf_path, pages = 1, density = density)
    img <- magick::image_convert(img, "png")
    magick::image_write(img, img_path)
    img_path
  }, error = function(e) {
    message("  [skip] Could not convert: ", basename(pdf_path),
            " — ", conditionMessage(e))
    NULL
  })
}

# ---------------------------------------------------------------------------
# .write_legend_sidecar()
# Internal helper used by every scSidekick plotting function that saves PDFs.
# Writes a plain-text .legend file beside the PDF so create_analysis_pptx()
# can pick it up automatically.
# ---------------------------------------------------------------------------
.write_legend_sidecar <- function(pdf_path, text) {
  if (is.null(text) || !nzchar(trimws(text))) return(invisible(NULL))
  sidecar <- paste0(tools::file_path_sans_ext(pdf_path), ".legend")
  tryCatch(writeLines(trimws(text), sidecar),
           error = function(e) NULL)
  invisible(sidecar)
}

# ---------------------------------------------------------------------------
# .pptx_legend_templates
# Named-pattern legend templates for external (non-scSidekick) figures.
# Patterns are regex, tested against the PDF basename (case-insensitive).
# First match wins — put more specific patterns before broader ones.
# Placeholders: {n_cells}, {n_samples}, {n_clusters}, {resolution}, {pcs},
#               {harmony_vars}, {groups}, {dataset}, {subset}
# ---------------------------------------------------------------------------
.pptx_legend_templates <- list(

  # ---- Quality Control ----
  list(
    pattern = "doubletfinder|doublet",
    legend  = paste0(
      "Bar plot showing the number of predicted doublets and singlets per sample ",
      "({n_samples} samples total), identified using DoubletFinder at an expected ",
      "doublet frequency of 7.5%. Bars are colored by sample identity. Only singlets ",
      "were retained for downstream analysis."
    )
  ),
  list(
    pattern = "miQC|cells to keep",
    legend  = paste0(
      "Bar plot showing the number of {cell_or_nucleus} classified for retention ('keep') or ",
      "removal per sample based on miQC quality scoring. miQC fits a probabilistic ",
      "mixture model on mitochondrial content and gene count distributions to flag ",
      "low-quality or damaged cells. Only cells labelled 'keep' were carried forward."
    )
  ),
  list(
    pattern = "before filtration",
    legend  = paste0(
      "Bar plot of total recovered {cell_or_nucleus} per sample prior to quality-based filtering ",
      "({n_samples} samples). Each bar represents the total number of {cell_or_nucleus} captured ",
      "for each donor or sample before any exclusion criteria were applied."
    )
  ),
  list(
    pattern = "after filtration.*barplot|barplot.*after filtration",
    legend  = paste0(
      "Bar plot of retained {cell_or_nucleus} per sample following doublet removal (DoubletFinder) ",
      "and probabilistic quality filtering (miQC). A total of {n_cells} {cell_or_nucleus} across ",
      "{n_samples} samples were retained for downstream analysis."
    )
  ),
  list(
    pattern = "qualitycontrol after filtration|qualitycontrol after",
    legend  = paste0(
      "Quality control metrics for {n_cells} retained {cell_or_nucleus} following doublet removal ",
      "and miQC-based filtering. Density distributions (top row) and violin plots ",
      "(bottom row) show total UMI counts (nCount_RNA), number of detected genes ",
      "(nFeature_RNA), mitochondrial read fraction (mitoRatio), and transcriptional ",
      "complexity (log10 genes per UMI), stratified by sample. Distributions confirm ",
      "successful removal of low-quality cells."
    )
  ),
  list(
    pattern = "qualitycontrol",
    legend  = paste0(
      "Quality control metrics for all recovered {cell_or_nucleus} prior to filtering. Density ",
      "distributions (top row) and violin plots (bottom row) display total UMI counts ",
      "(nCount_RNA), number of detected genes (nFeature_RNA), mitochondrial read fraction ",
      "(mitoRatio), and transcriptional complexity (log10 genes per UMI), stratified by ",
      "sample ({n_samples} samples). Dashed vertical lines indicate the filtering ",
      "thresholds subsequently applied."
    )
  ),

  # ---- Dimensionality Reduction ----
  list(
    pattern = "umap before and after harmony|before and after harmony",
    legend  = paste0(
      "UMAP visualizations of {cell_or_nucleus} colored by sample identity ({n_samples} samples), ",
      "before (left) and after (right) Harmony batch correction. Harmony integration ",
      "was applied across {harmony_vars}. Convergence of sample distributions following ",
      "integration indicates successful removal of technical batch variation while ",
      "preserving biological signal."
    )
  ),

  # ---- Clustering: Heatmaps ----
  list(
    pattern = "featureheatmap.*wikipathway|featureheatmap.*wiki",
    legend  = paste0(
      "Heatmap of the top 20 differentially expressed genes per cluster ",
      "({n_clusters} clusters; resolution {resolution}), annotated with WikiPathway ",
      "pathway enrichment terms. Rows represent genes grouped by cluster of origin; ",
      "columns represent individual cells. Expression is displayed as scaled ",
      "log-normalized counts (blue: low, red: high). Additional gene-level annotations ",
      "indicate transcription factor (TF) status and cell surface protein atlas (CSPA) ",
      "membership."
    )
  ),
  list(
    pattern = "featureheatmap",
    legend  = paste0(
      "Heatmap of the top 20 differentially expressed genes per cluster ",
      "({n_clusters} clusters; resolution {resolution}), annotated with Gene Ontology ",
      "Biological Process (GO:BP) enrichment terms. Rows represent genes grouped by the ",
      "cluster for which they were identified as top markers; columns represent individual ",
      "cells. Expression is displayed as scaled log-normalized counts (blue: low, red: high). ",
      "TF: transcription factor annotation; CSPA: cell surface protein atlas annotation."
    )
  ),
  list(
    pattern = "heatmap.*top|top.*heatmap",
    legend  = paste0(
      "Heatmap of the top 20 differentially expressed marker genes per cluster ",
      "({n_clusters} clusters; resolution {resolution}), identified by Wilcoxon rank-sum ",
      "test (Presto). Gene expression is displayed as mean-centred log-normalized counts ",
      "(blue: below mean, red: above mean). Cells are ordered by cluster, then by sample ",
      "identity. Column annotations (top) indicate cluster and sample of origin. Row ",
      "annotations (right) label the highest-ranked marker gene for each cluster."
    )
  ),

  # ---- Clustering: Dotplots ----
  list(
    pattern = "dotplot.*top|top.*dotplot",
    legend  = paste0(
      "Dot plot showing expression of the top differentially expressed marker genes across ",
      "{n_clusters} clusters (resolution {resolution}). Dot size indicates the percentage ",
      "of {cell_or_nucleus} expressing each gene within the cluster; dot colour reflects the mean ",
      "scaled expression level (Viridis plasma colour scale). Genes are grouped and ",
      "labelled by the cluster for which they were identified as top markers."
    )
  ),
  list(
    pattern = "dotplot",
    legend  = paste0(
      "Dot plot of canonical marker gene expression across {n_clusters} clusters ",
      "(resolution {resolution}). Dot size indicates the percentage of expressing cells; ",
      "colour reflects mean scaled expression (Viridis plasma scale). Gene groups ",
      "(x-axis panels) correspond to established cell-type marker sets."
    )
  ),

  # ---- Cell Type Annotation ----
  list(
    pattern = "AutoAssignment.*singleR|singleR.*AutoAssignment",
    legend  = paste0(
      "UMAP projections of {n_cells} {cell_or_nucleus} with automated cell-type annotations from ",
      "four complementary methods: unsupervised cluster identity (Louvain algorithm, ",
      "resolution {resolution}); SingleR against the Human Primary Cell Atlas; scType; ",
      "and rcellmarker. Concordance across methods was used to guide subsequent manual ",
      "cell-type assignment. {n_clusters} clusters were identified following Harmony ",
      "integration across {harmony_vars}."
    )
  ),
  list(
    pattern = "AutoAssignment",
    legend  = paste0(
      "UMAP projections of {n_cells} {cell_or_nucleus} showing automated cell-type annotations by ",
      "scType (left) and rcellmarker (right) alongside unsupervised cluster labels ",
      "(resolution {resolution}). Annotations were used to guide manual review and ",
      "assignment of the {n_clusters} identified clusters."
    )
  ),
  list(
    pattern = "Assignment Helper",
    legend  = paste0(
      "Summary panel supporting manual cell-type assignment. Top row: UMAP projections ",
      "colored by automated annotations (scType, SingleR against the Human Primary Cell ",
      "Atlas, rcellmarker, and unsupervised Louvain clusters at resolution {resolution}). ",
      "Bottom rows: dot plots of canonical marker genes across the {n_clusters} clusters, ",
      "used to validate automated annotations and resolve ambiguous assignments."
    )
  ),

  # ---- All UMAPs ----
  list(
    pattern = "All umaps|all umap",
    legend  = paste0(
      "Multi-panel UMAP visualization of {n_cells} {cell_or_nucleus} from {n_samples} samples ",
      "following Harmony batch correction across {harmony_vars}. Panels display: cluster ",
      "identity (Louvain algorithm, resolution {resolution}; {n_clusters} clusters), ",
      "sample of origin, experimental group ({groups}), fine-grained cell-type assignment, ",
      "and global cell-type assignment. Cell counts are annotated per group."
    )
  ),

  # ---- Assignment UMAPs ----
  list(
    pattern = "Assignment umap.*sample|umap.*iteration.*sample",
    legend  = paste0(
      "UMAP projections of manually assigned cell types faceted by sample ({n_samples} ",
      "samples). Each panel represents one donor, enabling assessment of inter-sample ",
      "heterogeneity in cell-type composition and reproducibility of cluster assignments ",
      "across the cohort. Cell types are coloured consistently across panels. Cluster ",
      "labels are overlaid."
    )
  ),
  list(
    pattern = "Assignment umap",
    legend  = paste0(
      "UMAP projection of {n_cells} {cell_or_nucleus} coloured by manually assigned cell-type ",
      "identity ({n_clusters} clusters; resolution {resolution}). Cell-type assignments ",
      "were determined by integration of automated annotation tools (scType, SingleR, ",
      "rcellmarker) with canonical marker gene expression."
    )
  ),

  # ---- Composition ----
  list(
    pattern = "iteration.*bar.*group|bar.*group",
    legend  = paste0(
      "Combined compositional and spatial visualization of cell-type distributions across ",
      "disease groups ({groups}). Left: proportional bar charts showing the percentage of ",
      "each manually assigned cell type per group. Right: UMAP projections faceted by group ",
      "with cell-type labels overlaid. Enables direct comparison of cell-type composition ",
      "across {groups} in the {dataset} cohort ({n_cells} {cell_or_nucleus}, {n_samples} samples)."
    )
  ),
  list(
    pattern = "iteration.*bar|bar.*iteration",
    legend  = paste0(
      "Combined compositional and spatial visualization. Left: stacked bar charts showing ",
      "the proportion of each manually assigned cell type per sample ({n_samples} samples). ",
      "Right: UMAP projections faceted by sample with cell-type labels overlaid. Based on ",
      "{n_cells} {cell_or_nucleus} assigned to {n_clusters} clusters (resolution {resolution}) ",
      "following Harmony integration across {harmony_vars}."
    )
  ),
  list(
    pattern = "percentage|composition|chord|alluvial|rose|trend",
    legend  = paste0(
      "Compositional analysis of cell-type proportions across samples and experimental ",
      "groups ({groups}). Panels include proportional trend plots, rose (polar) plots, ",
      "and a chord diagram illustrating the relationship between cell-type assignment ",
      "and sample identity. Analysis based on {n_cells} {cell_or_nucleus} from {n_samples} samples."
    )
  ),

  # ---- Feature Maps ----
  list(
    pattern = "featuremap",
    legend  = paste0(
      "UMAP feature plots showing log-normalized expression of individual marker genes ",
      "across {n_cells} {cell_or_nucleus}, split by experimental group ({groups}). Within each panel, ",
      "cells are ordered by expression level to highlight positive populations. Colour ",
      "scale indicates log-normalized expression (blue: low/absent, red: high)."
    )
  )
)

# ---------------------------------------------------------------------------
# .pptx_assay_info()
# Converts the assay_type string stored in analysis_params.json into a
# human-readable long name and the correct word for individual observations.
# Recognised values: "scRNAseq", "snRNAseq", "scATACseq", "scMultiome",
# "Spatial".  Unknown values are passed through as-is.
# ---------------------------------------------------------------------------
.pptx_assay_info <- function(params) {
  at <- params$assay_type %||% "scRNAseq"
  switch(at,
    snRNAseq   = list(long = "single-nucleus RNA-seq",  cell = "nuclei"),
    scRNAseq   = list(long = "single-cell RNA-seq",     cell = "cells"),
    scATACseq  = list(long = "single-cell ATAC-seq",    cell = "cells"),
    scMultiome = list(long = "single-cell multiome",            cell = "cells"),
    Spatial    = list(long = "spatial transcriptomics",          cell = "spots"),
    Visium     = list(long = "Visium spatial transcriptomics",   cell = "spots"),
    VisiumHD   = list(long = "Visium HD spatial transcriptomics", cell = "bins"),
    Xenium     = list(long = "Xenium in situ transcriptomics",   cell = "cells"),
    list(long = at, cell = "cells")   # passthrough for custom / unknown types
  )
}

# ---------------------------------------------------------------------------
# .pptx_fill_template()
# Substitutes {placeholder} tokens in a legend template string using the
# values stored in the params list (from analysis_params.json).
# ---------------------------------------------------------------------------
.pptx_fill_template <- function(template, params) {
  n_cells  <- if (!is.null(params$n_cells_final))
    format(as.integer(params$n_cells_final), big.mark = ",") else "N"
  n_samp   <- params$n_samples  %||% "N"
  resol    <- params$resolution %||% "N"
  n_clust  <- params$n_clusters %||% "N"
  pcs      <- params$pcs_used   %||% "N"
  hv       <- if (!is.null(params$harmony_vars))
    paste(unlist(params$harmony_vars), collapse = ", ") else "N/A"
  grps     <- if (!is.null(params$groups))
    paste(unlist(params$groups), collapse = ", ") else "N/A"
  dataset  <- params$dataset %||% "the dataset"
  subset_n <- params$subset  %||% ""

  ai <- .pptx_assay_info(params)   # assay-type derived strings

  subs <- list(
    "{n_cells}"         = as.character(n_cells),
    "{n_samples}"       = as.character(n_samp),
    "{resolution}"      = as.character(resol),
    "{n_clusters}"      = as.character(n_clust),
    "{pcs}"             = as.character(pcs),
    "{harmony_vars}"    = hv,
    "{groups}"          = grps,
    "{dataset}"         = dataset,
    "{subset}"          = subset_n,
    "{assay_type}"      = ai$long,   # e.g. "single-nucleus RNA-seq"
    "{cell_or_nucleus}" = ai$cell    # e.g. "{cell_or_nucleus}" or "cells" or "spots"
  )
  for (nm in names(subs))
    template <- gsub(nm, subs[[nm]], template, fixed = TRUE)
  template
}

# ---------------------------------------------------------------------------
# .generate_fallback_legend()
# Two-tier legend generation for PDFs that have no .legend sidecar file:
#   Tier 1 — try specific filename-pattern templates (.pptx_legend_templates)
#             with {placeholder} substitution from the params JSON.
#   Tier 2 — fall back to section-level generic prose if no pattern matches.
# ---------------------------------------------------------------------------
.generate_fallback_legend <- function(filename, params) {

  bn <- basename(filename)

  # ---- Tier 1: specific pattern templates ----
  for (tmpl in .pptx_legend_templates) {
    if (grepl(tmpl$pattern, bn, ignore.case = TRUE))
      return(.pptx_fill_template(tmpl$legend, params))
  }

  # ---- Tier 2: section-level generic prose ----
  section <- .pptx_assign_section(filename)

  ds  <- trimws(paste(params$dataset %||% "", params$subset %||% ""))
  if (!nzchar(ds)) ds <- "the dataset"

  pcs  <- params$pcs_used
  res  <- params$resolution
  ncl  <- params$n_clusters
  nc   <- if (!is.null(params$n_cells_final))
            format(params$n_cells_final, big.mark = ",") else NULL
  harm <- if (!is.null(params$harmony_vars))
            paste(unlist(params$harmony_vars), collapse = " and ") else NULL

  ai        <- .pptx_assay_info(params)
  at_long   <- ai$long    # e.g. "single-nucleus RNA-seq"
  cell_word <- ai$cell    # e.g. "nuclei", "cells", "spots"

  pc_str  <- if (!is.null(pcs)) paste0("the top ", pcs, " principal components")
             else "principal components"
  res_str <- if (!is.null(res)) paste0("resolution ", res) else "the selected resolution"
  ncl_str <- if (!is.null(ncl)) paste0(ncl, " clusters") else "clusters"
  nc_str  <- if (!is.null(nc))  paste0(nc,  " ", cell_word) else cell_word

  switch(section,

    "Quality Control" = paste0(
      "Quality control metrics for ", ds, ". Cell-level features including ",
      "mitochondrial content, library complexity, and doublet scores were assessed ",
      "to identify and exclude low-quality cells prior to downstream analysis."
    ),

    "Dimensionality Reduction" = {
      harm_sent <- if (!is.null(harm))
        paste0("Harmony integration was applied across ", harm,
               " to correct for batch effects. ")
      else ""
      paste0(
        "Dimensionality reduction of ", ds, ". Principal component analysis (PCA) ",
        "was performed on highly variable genes. ", harm_sent,
        "UMAP embedding was computed on ", pc_str,
        " to enable two-dimensional visualisation of the transcriptional landscape."
      )
    },

    "Cluster Markers" = paste0(
      "Unsupervised clustering of ", ds, " at ", res_str, ", yielding ", ncl_str, ". ",
      "A shared nearest-neighbour graph was constructed on ", pc_str,
      " and cell communities were identified using the Louvain algorithm.",
      if (!is.null(nc)) paste0(" The dataset comprises ", nc_str, ".") else ""
    ),

    "Dotplots" = paste0(
      "Dot plot showing scaled gene expression and the proportion of expressing ",
      "cells across cell populations in ", ds, ". Dot size reflects the fraction ",
      "of cells with detectable expression; dot colour reflects mean scaled expression ",
      "within each group."
    ),

    "Cell Type Annotation" = paste0(
      "Cell type annotation of ", ds, ". Cluster identities were assigned based on ",
      "canonical marker gene expression profiles and, where applicable, automated ",
      "label transfer using reference datasets."
    ),

    "All UMAPs" = paste0(
      "UMAP embeddings of ", ds, " showing the distribution of ", nc_str,
      " across ", ncl_str, ". Single cells are projected onto a two-dimensional ",
      "representation computed from ", pc_str,
      ", coloured by the indicated metadata variable."
    ),

    "Composition" = paste0(
      "Compositional analysis of cell type proportions across conditions in ", ds, ". ",
      "Bars or arcs represent the relative abundance of each cell population ",
      "per sample or group, enabling detection of condition-associated shifts in ",
      "cellular composition."
    ),

    "Feature Maps" = paste0(
      "UMAP feature plots showing gene expression across cell populations in ", ds, ". ",
      "Log-normalised expression values are displayed on a continuous colour scale ",
      "from low (dark blue) to high (dark red). Higher-expressing cells are plotted ",
      "on top to highlight positive populations."
    ),

    "Pathway Analysis" = {
      db_str <- if (!is.null(params$gsea_databases))
        paste(unlist(params$gsea_databases), collapse = ", ")
      else "MSigDB gene-set collections"
      grp_str <- if (!is.null(params$gsea_group_by)) params$gsea_group_by else "groups"
      spl_str <- if (!is.null(params$gsea_split_by)) params$gsea_split_by else "cell types"
      paste0(
        "Gene set enrichment analysis (GSEA) result for ", ds, ". ",
        "Differential expression between ", grp_str, " groups was computed ",
        "within each ", spl_str, " cell-type subset using a Wilcoxon rank-sum ",
        "test (presto), and enrichment scores (NES) were calculated with fgsea ",
        "against the following gene-set collections: ", db_str, ". ",
        "Positive NES (red/warm) indicates up-regulation relative to the ",
        "comparison group; negative NES (blue/cool) indicates down-regulation."
      )
    },

    "CellChat" = {
      grp_str  <- if (!is.null(params$cellchat_group_col))
        params$cellchat_group_col else "experimental groups"
      id_str   <- if (!is.null(params$cellchat_ident_col))
        params$cellchat_ident_col else "cell type identities"
      sp_str   <- if (!is.null(params$cellchat_species))
        params$cellchat_species else "the"
      paste0(
        "Cell-cell communication figure for ", ds, ". ",
        "CellChat was applied to ", sp_str, " ", at_long, " data, ",
        "with cell-type identities drawn from '", id_str, "' and conditions ",
        "defined by '", grp_str, "'. ",
        "Ligand-receptor interactions were inferred from CellChatDB using a ",
        "mass-action communication probability model. Circle and chord diagrams ",
        "show the number and weight of predicted interactions between cell populations; ",
        "heatmaps show pairwise communication probability for individual signalling pathways."
      )
    },

    # default / "Other"
    paste0("Figure generated during the analysis of ", ds, ".")
  )
}

# ---------------------------------------------------------------------------

#' Log Seurat analysis parameters to a JSON file
#'
#' Reads the command history of a Seurat object to extract key analysis
#' parameters (resolution, PCs, Harmony variables, cell count, cluster count)
#' and writes them to `analysis_params.json` in `output_dir`.  When
#' `ExtractMethods()` is available in the session the full structured summary
#' it returns is used; otherwise parameters are read directly from
#' `obj@commands`.
#'
#' Call this once, after `RunHarmony` / `FindClusters` / `RunUMAP` are
#' complete.  Then call [create_analysis_pptx()] at the end of the Rmd.
#'
#' @param obj A Seurat object.
#' @param output_dir Character. Directory where figures are saved (and where
#'   `analysis_params.json` will be written).
#' @param dataset Character. Human-readable dataset label
#'   (e.g. `"SEAAD snRNAseq"`).
#' @param subset_name Character. Subset label (e.g. `"All_Clusters"`).
#' @param assay_type Character. Sequencing modality. Controls the word used for
#'   individual observations in auto-generated legends and methods text.
#'   Recognised values (case-sensitive):
#'   \itemize{
#'     \item `"scRNAseq"` → "single-cell RNA-seq" / "cells" (default)
#'     \item `"snRNAseq"` → "single-nucleus RNA-seq" / "nuclei"
#'     \item `"scATACseq"` → "single-cell ATAC-seq" / "cells"
#'     \item `"scMultiome"` → "single-cell multiome" / "cells"
#'     \item `"Spatial"` → "spatial transcriptomics" / "spots"
#'     \item `"Visium"` → "Visium spatial transcriptomics" / "spots"
#'     \item `"VisiumHD"` → "Visium HD spatial transcriptomics" / "bins"
#'     \item `"Xenium"` → "Xenium in situ transcriptomics" / "cells"
#'   }
#'   Any other string is passed through as-is (long name) with "cells" as the
#'   observation word.
#' @param params_json Character. Full path for the JSON file.  Defaults to
#'   `file.path(output_dir, "analysis_params.json")`.
#'
#' @return Invisibly returns the logged parameter list.
#' @export
log_analysis_params <- function(obj,
                                output_dir,
                                dataset     = "",
                                subset_name = "",
                                assay_type  = "scRNAseq",
                                params_json = file.path(output_dir,
                                                        "analysis_params.json")) {

  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required. Install with: install.packages('jsonlite')")

  params_log <- list()

  # ---- Try ExtractMethods (uses obj@commands) ----
  em <- NULL
  if (exists("ExtractMethods", mode = "function")) {
    em <- tryCatch(
      ExtractMethods(obj, cite_seurat = TRUE),
      error = function(e) {
        message("ExtractMethods() failed: ", conditionMessage(e),
                "\nFalling back to inline extraction.")
        NULL
      }
    )
  }

  if (!is.null(em)) {
    s <- em$summary
    params_log$resolution    <- s$clustering$resolution
    params_log$pcs_used      <- s$umap$dims
    params_log$harmony_vars  <- s$integration$params$group_by
    params_log$n_cells_final <- s$n_cells
    params_log$n_clusters    <- s$clustering$n_clusters
    params_log$methods_text  <- em$methods_text

  } else {
    cmds  <- names(obj@commands)
    get_p <- function(cmd, param)
      if (cmd %in% cmds) obj@commands[[cmd]]@params[[param]] else NULL

    umap_dims <- get_p("RunUMAP",       "dims")
    fn_dims   <- get_p("FindNeighbors", "dims")
    pcs_used  <- if (!is.null(umap_dims)) max(umap_dims) else
                 if (!is.null(fn_dims))   max(fn_dims)   else NULL

    params_log$resolution    <- get_p("FindClusters", "resolution")
    params_log$pcs_used      <- pcs_used
    params_log$harmony_vars  <- if ("RunHarmony" %in% cmds)
                                  get_p("RunHarmony", "group.by.vars") else NULL
    params_log$n_cells_final <- ncol(obj)
    params_log$n_clusters    <- if ("seurat_clusters" %in% colnames(obj@meta.data))
                                  length(levels(obj@meta.data$seurat_clusters)) else NULL
  }

  # ---- Always add these ----
  params_log$dataset     <- dataset
  params_log$subset      <- subset_name
  params_log$assay_type  <- assay_type
  params_log$date        <- format(Sys.Date(), "%Y-%m-%d")

  # Sample count
  for (col in c("Sample", "sample", "orig.ident", "Donor.ID")) {
    if (col %in% colnames(obj@meta.data)) {
      params_log$n_samples <- length(unique(obj@meta.data[[col]]))
      break
    }
  }

  # Groups
  for (col in c("Group", "group", "condition", "Condition", "Treatment")) {
    if (col %in% colnames(obj@meta.data)) {
      params_log$groups <- levels(factor(obj@meta.data[[col]]))
      break
    }
  }

  params_log <- Filter(Negate(is.null), params_log)
  jsonlite::write_json(params_log, path = params_json, pretty = TRUE,
                       auto_unbox = TRUE)
  message("Params logged to: ", params_json)
  invisible(params_log)
}


#' Manually log a figure legend for a saved PDF
#'
#' Writes a `.legend` sidecar file next to an existing PDF so that
#' [create_analysis_pptx()] can include the text automatically on the
#' corresponding slide.  Use this for figures generated outside scSidekick
#' (e.g. doublet score plots, external QC tools) or to override an
#' automatically generated legend.
#'
#' @param out_dir Character.  Directory containing the PDF.
#' @param filename Character.  PDF filename (with or without `.pdf` extension).
#' @param text Character.  The legend text to attach.  Should be one or two
#'   complete sentences suitable for a figure caption.
#'
#' @return Invisibly returns the path of the written `.legend` file.
#' @export
log_figure_legend <- function(out_dir, filename, text) {
  fname <- if (grepl("\\.pdf$", filename, ignore.case = TRUE))
    filename else paste0(filename, ".pdf")
  pdf_path <- file.path(out_dir, fname)
  sidecar  <- .write_legend_sidecar(pdf_path, text)
  if (!is.null(sidecar))
    message("Legend written: ", sidecar)
  invisible(sidecar)
}


# Walk up the directory tree from `start_dir` until analysis_params.json is
# found or the filesystem root is reached.  Returns the path or NULL.
.find_params_json <- function(start_dir) {
  dir <- normalizePath(start_dir, mustWork = FALSE)
  repeat {
    candidate <- file.path(dir, "analysis_params.json")
    if (file.exists(candidate)) return(candidate)
    parent <- dirname(dir)
    if (parent == dir) break   # reached root
    dir <- parent
  }
  NULL
}

#' Build a PowerPoint summary from an analysis output folder
#'
#' Collects all PDFs in `output_dir`, converts them to images with
#' [magick::image_read_pdf()], groups them into labelled sections based on
#' filename patterns, and assembles a polished `.pptx` via the
#' [officer][officer::officer-package] package.
#'
#' A parameters / methods overview slide is generated automatically from
#' `analysis_params.json` (written by [log_analysis_params()]).
#'
#' @param output_dir Character. Folder containing PDF figures.
#' @param params_json Character or `NULL`. Path to `analysis_params.json`.
#'   When `NULL` (default) the function walks up the directory tree from
#'   `output_dir` until it finds the file — this lets you call
#'   `create_analysis_pptx` on a sub-folder (e.g., the CellChat or GSEA
#'   output directory) and still get the shared analysis parameters from the
#'   parent folder.
#' @param out_pptx Character or `NULL`.  Where to save the `.pptx`.  When
#'   `NULL` (default) a filename is auto-generated from the dataset/subset
#'   labels and the current date, saved inside `output_dir`.
#' @param slide_format Character. Slide aspect ratio preset: `"standard"`
#'   (4:3, 10 × 7.5 in, default) or `"widescreen"` (16:9, 13.33 × 7.5 in).
#'   Override with explicit `slide_width` / `slide_height` if needed.
#' @param slide_width,slide_height Numeric. Explicit slide dimensions in
#'   inches. When provided these override `slide_format`.
#' @param pdf_density Integer. DPI used when rasterising PDFs.  Increase
#'   (e.g. `200`) for sharper figures at the cost of speed. Default `150`.
#' @param include Character vector or `NULL`. Section names to include; all
#'   other sections are dropped. When `NULL` (default) all sections are
#'   included. PDFs are assigned to sections by filename pattern matching —
#'   see `.pptx_section_rules` for the rules. Valid section names:
#'   \itemize{
#'     \item `"Quality Control"` — QC density plots, doublet scores, filtration
#'       barplots
#'     \item `"Dimensionality Reduction"` — PCA elbow, UMAP before/after
#'       Harmony, RNA UMAPs
#'     \item `"Cluster Markers"` — ComplexHeatmap heatmaps, presto marker
#'       tables
#'     \item `"Dotplots"` — any PDF with "dotplot" in the filename
#'     \item `"Cell Type Annotation"` — AutoAssignment UMAPs, SingleR, scType,
#'       featureheatmaps, assignment helper
#'     \item `"All UMAPs"` — the combined multi-panel UMAP PDF
#'     \item `"Composition"` — barplots, chord diagrams, rose/trend plots,
#'       iteration-with-bar figures
#'     \item `"Feature Maps"` — GenerateFeatureMaps per-gene UMAPs
#'     \item `"Pathway Analysis"` — GSEA heatmaps and lollipop PDFs
#'     \item `"CellChat"` — RunCellChat, CompareCellChat, and
#'       RankCellChatPathways output PDFs
#'     \item `"Other"` — anything not matched by the rules above
#'   }
#' @param exclude Character vector or `NULL`. Section names to drop (applied
#'   after `include`). `NULL` excludes nothing. Same valid values as `include`.
#'   Typical use: `exclude = c("Feature Maps", "CellChat")` to keep the main
#'   analysis compact while running separate sub-folder PPTXs for those
#'   heavy sections.
#'
#' @return Invisibly returns the path to the saved `.pptx` file.
#' @export
create_analysis_pptx <- function(
    output_dir,
    params_json      = NULL,
    out_pptx         = NULL,
    slide_format     = "standard",
    slide_width      = NULL,
    slide_height     = NULL,
    pdf_density      = 150,
    include_legends  = TRUE,
    legend_font_size = 8,
    include          = NULL,
    exclude          = NULL
) {
  for (pkg in c("officer", "magick", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. Install with: install.packages('",
           pkg, "')")
  }

  # ---- slide dimensions ----
  format_dims <- switch(
    match.arg(slide_format, c("standard", "widescreen")),
    standard   = c(10.00, 7.5),
    widescreen = c(13.33, 7.5)
  )
  if (is.null(slide_width))  slide_width  <- format_dims[1]
  if (is.null(slide_height)) slide_height <- format_dims[2]

  # ---- params: walk up tree if no explicit path given ----
  if (is.null(params_json))
    params_json <- .find_params_json(output_dir)

  if (!is.null(params_json) && file.exists(params_json)) {
    params <- jsonlite::read_json(params_json, simplifyVector = TRUE)
    if (normalizePath(dirname(params_json)) !=
        normalizePath(output_dir))
      message("Using params from: ", params_json)
  } else {
    params <- list(dataset = basename(output_dir), subset = "",
                   date = format(Sys.Date()))
    warning("No analysis_params.json found — call log_analysis_params() ",
            "before create_analysis_pptx().")
  }

  obj_name    <- params$dataset %||% ""
  subset_name <- params$subset  %||% ""
  footer_text <- .pptx_make_footer(params)

  # ---- collect & group PDFs ----
  pdfs <- list.files(output_dir, pattern = "\\.pdf$", full.names = TRUE,
                     recursive = TRUE)
  if (length(pdfs) == 0) stop("No PDF files found in: ", output_dir)

  df <- data.frame(
    path    = pdfs,
    section = sapply(pdfs, .pptx_assign_section),
    title   = sapply(pdfs, .pptx_clean_title,
                     obj_name = obj_name, subset_name = subset_name),
    stringsAsFactors = FALSE
  )

  section_order <- c(
    "Quality Control", "Dimensionality Reduction", "Cluster Markers",
    "Dotplots", "Cell Type Annotation", "All UMAPs", "Composition",
    "Feature Maps", "Pathway Analysis", "CellChat", "Other"
  )
  df$section <- factor(df$section, levels = section_order)
  df <- df[order(df$section, df$path), ]

  # ---- include / exclude filter ----
  if (!is.null(include)) {
    df <- df[as.character(df$section) %in% include, , drop = FALSE]
    if (nrow(df) == 0)
      stop("No figures match the 'include' sections: ",
           paste(include, collapse = ", "))
  }
  if (!is.null(exclude)) {
    df <- df[!as.character(df$section) %in% exclude, , drop = FALSE]
    if (nrow(df) == 0)
      stop("All figures were removed by the 'exclude' sections: ",
           paste(exclude, collapse = ", "))
  }
  df <- droplevels(df)

  n_figs <- nrow(df)
  n_sec  <- length(unique(df$section))
  message("Building PPTX: ", n_figs, " figures across ", n_sec, " sections",
          if (!is.null(include) || !is.null(exclude)) " (filtered)" else "",
          "...")

  # ---- typography ----
  head_fp    <- officer::fp_text(font.size = 22, bold = TRUE,  color = "#1F3864")
  subhead_fp <- officer::fp_text(font.size = 13, bold = TRUE,  color = "#2E4057")
  body_fp    <- officer::fp_text(font.size = 12, bold = FALSE, color = "#2C3E50")
  methods_fp <- officer::fp_text(font.size = 10, bold = FALSE, color = "#444444")
  title_fp   <- officer::fp_text(font.size = 16, bold = TRUE,  color = "#1F3864")
  footer_fp  <- officer::fp_text(font.size = 9,  italic = TRUE, color = "#888888")
  legend_fp  <- officer::fp_text(font.size = legend_font_size,
                                 bold = FALSE, italic = FALSE, color = "#222222")

  # Space reserved for legend box (inches). Zero when legends are off.
  LEG_H  <- if (include_legends) 0.78 else 0   # text box height
  LEG_GAP <- if (include_legends) 0.05 else 0   # gap between image and box

  # ---- build presentation ----
  prs <- officer::read_pptx()

  # -- Title slide --
  prs <- officer::add_slide(prs, layout = "Title Slide", master = "Office Theme")
  prs <- officer::ph_with(prs,
    value    = trimws(paste(obj_name, subset_name)),
    location = officer::ph_location_type(type = "ctrTitle"))
  prs <- officer::ph_with(prs,
    value    = paste0("Single-Cell Analysis Summary\n",
                      params$date %||% format(Sys.Date())),
    location = officer::ph_location_type(type = "subTitle"))

  # -- Parameters + Methods overview slide --
  prs <- officer::add_slide(prs, layout = "Blank", master = "Office Theme")

  param_fields <- list(
    list(label = "Dataset",      value = trimws(paste(obj_name, subset_name))),
    list(label = "Resolution",   value = params$resolution),
    list(label = "PCs used",     value = params$pcs_used),
    list(label = "Clusters",     value = params$n_clusters),
    list(label = "Harmony vars",
         value = if (!is.null(params$harmony_vars))
                   paste(unlist(params$harmony_vars), collapse = ", ") else NULL),
    list(label = "Cells",
         value = if (!is.null(params$n_cells_final))
                   format(params$n_cells_final, big.mark = ",") else NULL),
    list(label = "Samples",      value = params$n_samples),
    list(label = "Groups",
         value = if (!is.null(params$groups))
                   paste(unlist(params$groups), collapse = ", ") else NULL),
    list(label = "Date",         value = params$date)
  )

  overview_blocks <- list(
    officer::fpar(officer::ftext("Analysis Parameters", head_fp)),
    officer::fpar(officer::ftext("", body_fp))
  )
  for (pf in param_fields) {
    v <- pf$value
    if (!is.null(v) && length(v) > 0 && !all(is.na(v))) {
      overview_blocks[[length(overview_blocks) + 1]] <- officer::fpar(
        officer::ftext(paste0(pf$label, ":  "), subhead_fp),
        officer::ftext(as.character(v), body_fp)
      )
    }
  }

  if (!is.null(params$methods_text) && nchar(params$methods_text) > 0) {
    overview_blocks <- c(overview_blocks, list(
      officer::fpar(officer::ftext("", body_fp)),
      officer::fpar(officer::ftext("Methods paragraph (draft):", subhead_fp)),
      officer::fpar(officer::ftext(params$methods_text, methods_fp))
    ))
  }

  prs <- officer::ph_with(prs,
    value    = do.call(officer::block_list, overview_blocks),
    location = officer::ph_location(left = 0.6, top = 0.4,
                                    width  = slide_width  - 1.2,
                                    height = slide_height - 0.7))

  # -- Content slides --
  sections_done <- character(0)

  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    sec <- as.character(row$section)

    if (!sec %in% sections_done) {
      prs <- officer::add_slide(prs, layout = "Section Header",
                                master = "Office Theme")
      prs <- officer::ph_with(prs, value = sec,
                              location = officer::ph_location_type(type = "title"))
      sections_done <- c(sections_done, sec)
      message("  [", sec, "]")
    }

    img_path <- .pptx_pdf_to_png(row$path, density = pdf_density)
    if (is.null(img_path)) next

    message("    ", basename(row$path))

    # ---- Read or generate the legend sentence ----
    slide_legend <- ""
    if (include_legends) {
      sidecar <- paste0(tools::file_path_sans_ext(row$path), ".legend")
      if (file.exists(sidecar)) {
        slide_legend <- trimws(paste(readLines(sidecar, warn = FALSE),
                                     collapse = " "))
      } else {
        slide_legend <- tryCatch(
          .generate_fallback_legend(row$path, params),
          error = function(e) ""
        )
      }
      # Hard cap to avoid overflow (≈5 lines at 8pt in the box)
      if (nchar(slide_legend) > 750)
        slide_legend <- paste0(substr(slide_legend, 1, 747), "…")
    }

    # ---- Layout dimensions ----
    show_leg  <- include_legends && nzchar(slide_legend)
    avail_w   <- slide_width  - 0.5          # horizontal space for figure
    avail_h   <- slide_height - 1.1 - (if (show_leg) LEG_H + LEG_GAP else 0)

    # Preserve the image's natural aspect ratio — scale to fit within
    # avail_w × avail_h without stretching or squishing.
    img_info  <- tryCatch(
      magick::image_info(magick::image_read(img_path)),
      error = function(e) NULL
    )
    if (!is.null(img_info) && img_info$width > 0 && img_info$height > 0) {
      img_aspect <- img_info$width / img_info$height
      box_aspect <- avail_w / avail_h
      if (img_aspect >= box_aspect) {
        # Image is wider than the box → constrain by width
        render_w <- avail_w
        render_h <- avail_w / img_aspect
      } else {
        # Image is taller than the box → constrain by height
        render_h <- avail_h
        render_w <- avail_h * img_aspect
      }
    } else {
      render_w <- avail_w
      render_h <- avail_h
    }
    # Center the scaled image within the available space
    render_left <- 0.25 + (avail_w - render_w) / 2
    render_top  <- 0.58  + (avail_h - render_h) / 2
    leg_top     <- 0.58  + avail_h + LEG_GAP   # legend below full avail box

    prs <- officer::add_slide(prs, layout = "Blank", master = "Office Theme")

    # Title
    prs <- officer::ph_with(prs,
      value    = officer::fpar(officer::ftext(row$title, title_fp)),
      location = officer::ph_location(left = 0.25, top = 0.08,
                                      width = slide_width - 0.5, height = 0.48))

    # Figure — placed at its natural aspect-ratio-correct size, centred
    prs <- officer::ph_with(prs,
      value    = officer::external_img(img_path,
                                       width  = render_w,
                                       height = render_h),
      location = officer::ph_location(left   = render_left,
                                      top    = render_top,
                                      width  = render_w,
                                      height = render_h))

    # Legend sentence
    if (show_leg) {
      prs <- officer::ph_with(prs,
        value    = officer::fpar(officer::ftext(slide_legend, legend_fp)),
        location = officer::ph_location(left = 0.25, top = leg_top,
                                        width = slide_width - 0.5, height = LEG_H))
    }

    # Footer params
    if (nchar(footer_text) > 0) {
      prs <- officer::ph_with(prs,
        value    = officer::fpar(officer::ftext(footer_text, footer_fp)),
        location = officer::ph_location(left = 0.25,
                                        top  = slide_height - 0.4,
                                        width  = slide_width - 0.5,
                                        height = 0.35))
    }
  }

  # ---- save ----
  if (is.null(out_pptx)) {
    safe     <- gsub("[^A-Za-z0-9_-]", "_", trimws(paste(obj_name, subset_name)))
    safe     <- gsub("_+", "_", safe)
    out_pptx <- file.path(output_dir,
                          paste0(safe, "_summary_",
                                 format(Sys.Date(), "%Y%m%d"), ".pptx"))
  }

  print(prs, target = out_pptx)
  message("\nSaved: ", out_pptx)
  invisible(out_pptx)
}
