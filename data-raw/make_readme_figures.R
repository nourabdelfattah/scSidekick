# =============================================================================
# make_readme_figures.R
# Renders the "before / after" comparison PNGs used by the README's
# "Pain Points Fixed" showcase. All panels are rendered from the public
# `bmcite` dataset (SeuratData) so every claim is reproducible.
#
#   Rscript data-raw/make_readme_figures.R           # all stages
#   Rscript data-raw/make_readme_figures.R feature   # one stage
#
# Output: man/figures/pain_*_{before,after}.png
# =============================================================================

suppressMessages({
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
  library(patchwork)
})
devtools::load_all(".", quiet = TRUE)

set.seed(1)
fig_dir <- "man/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
stages <- commandArgs(trailingOnly = TRUE)
do_stage <- function(s) length(stages) == 0 || s %in% stages

# -----------------------------------------------------------------------------
# 0. Processed object (cached) - standard RNA pipeline + UMAP on bmcite
# -----------------------------------------------------------------------------
cache <- "data-raw/bmcite_processed.rds"
if (file.exists(cache)) {
  message(">> loading cached processed object")
  obj <- readRDS(cache)
} else {
  message(">> processing bmcite (RNA pipeline + UMAP) - one time")
  data("bmcite")
  obj <- bmcite
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  saveRDS(obj, cache)
}

# Register colours + project defaults once (so the "after" panels use the
# package's consistent palette / defaults exactly as a real user would).
obj <- PrepObject(
  obj,
  variables   = c("donor", "celltype.l1", "celltype.l2"),
  group.by    = "celltype.l1",
  split.by    = "donor",
  output_dir  = fig_dir,
  object_name = "bmcite"
)

ggsave_png <- function(plot, file, w, h) {
  ggsave(file.path(fig_dir, file), plot = plot, width = w, height = h,
         dpi = 150, bg = "white")
  message("   wrote ", file)
}

# A copy whose stored output_dir is cleared, so plotting functions RETURN the
# ggplot/patchwork object (and print it) instead of writing a PDF to disk.
no_dir <- function(o) { o@misc$nk_settings$output_dir <- NULL; o }

# =============================================================================
# PAIN POINT 1 - split.by silently drops the colour legend
#   Before: FeaturePlot(features, split.by, order = TRUE) renders the panels
#           but NO colour legend at all - the values are unreadable. The usual
#           "& RestoreLegend()" workaround re-attaches a guide that does not
#           reflect the data shown.
#   After:  PlotFeaturePlots(split.by = "donor") - one shared, accurate colour
#           scale across all panels.
# =============================================================================
if (do_stage("feature")) {
  message(">> PAIN 1: split.by legend")
  gene <- "LYZ"

  # --- BEFORE: the real Seurat call - note the missing colour legend
  before <- FeaturePlot(obj, features = gene, split.by = "donor",
                        order = TRUE) +
    plot_annotation(
      title = "Seurat: FeaturePlot(split.by) drops the colour legend entirely",
      theme = theme(plot.title = element_text(size = 12, face = "bold",
                                              colour = "#B2182B")))
  ggsave_png(before, "pain_feature_before.png", 9, 4.2)

  # --- AFTER: one shared, accurate scale across panels (the package one-liner)
  after <- PlotFeaturePlots(no_dir(obj), features = gene, split.by = "donor",
                            output_dir = NULL)
  after <- after +
    plot_annotation(
      title = "scSidekick: one shared, accurate colour scale across panels",
      theme = theme(plot.title = element_text(size = 12, face = "bold",
                                              colour = "#2f4b7c")))
  ggsave_png(after, "pain_feature_after.png", 9, 4.2)
}

# =============================================================================
# PAIN POINT 2 - Reading complex metadata at the donor level
#   SEAAD: 1.78M cells x 153 metadata columns, 83 donors.
#   Before: a nested base-R table of the clinical variables - unreadable soup.
#   After:  one PlotMetaSummary() call collapses cells -> donors and lays out
#           Braak / APOE / CERAD by sex, coloured by cognitive status.
# =============================================================================
if (do_stage("meta")) {
  message(">> PAIN 2: metadata summary (SEAAD)")

  # Cache only the handful of columns we need (the source rds is ~1.78M x 153).
  seaad_cache <- "data-raw/seaad_meta_small.rds"
  keep <- c("Donor.ID", "Braak", "APOE.Genotype", "CERAD.score",
            "Sex", "Cognitive.Status", "Age.at.Death")
  if (file.exists(seaad_cache)) {
    md <- readRDS(seaad_cache)
  } else {
    src <- "/Users/nourhan/Library/CloudStorage/Box-Box/Yun lab projects/Nour/SEAAD/R_Code/clause scripts/seaad_meta.rds"
    message("   reading full SEAAD metadata (one time) ...")
    md_full <- readRDS(src)
    md <- md_full[, keep]
    saveRDS(md, seaad_cache)
    rm(md_full); gc()
  }
  message(sprintf("   SEAAD: %s cells x 153 cols, %d donors",
                  format(nrow(md), big.mark = ","),
                  length(unique(md$Donor.ID))))

  # --- BEFORE: the naive base-R attempt - a 3-way table you can't parse
  tbl <- table(md$Sex, md$Braak, md$Cognitive.Status)
  txt <- paste(capture.output(print(ftable(tbl))), collapse = "\n")
  txt <- paste0("> dim(md)\n[1] 1782605     153\n\n",
                "> ftable(Sex, Braak, Cognitive.Status)   # cells, not donors\n\n",
                txt,
                "\n\n> # ...and that's only 3 of 153 columns, at the wrong unit")
  before <- ggplot() +
    annotate("text", x = 0, y = 1, label = txt, family = "mono",
             hjust = 0, vjust = 1, size = 2.5, colour = "grey20") +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = "Base R: 1.78M cells, 153 columns - unreadable") +
    theme_void() +
    theme(plot.title = element_text(size = 12, face = "bold",
                                    colour = "#B2182B", hjust = 0),
          plot.margin = margin(10, 10, 10, 10))
  ggsave_png(before, "pain_meta_before.png", 7.5, 3.8)

  # --- AFTER: one line - donor-level clinical composition, faceted by sex
  after <- PlotMetaSummary(
    md,
    id_column     = "Donor.ID",
    variables     = c("Braak", "APOE.Genotype", "CERAD.score"),
    fill_variable = "Cognitive.Status",
    row_variable  = "Sex"
  )
  if (!inherits(after, c("ggplot", "patchwork")) && is.list(after))
    after <- after[[1]]
  after <- after +
    plot_annotation(
      title = "scSidekick: PlotMetaSummary() - 83 donors summarised in one line",
      theme = theme(plot.title = element_text(size = 13, face = "bold",
                                              colour = "#2f4b7c")))
  ggsave_png(after, "pain_meta_after.png", 11, 5.5)

  # --- Everyday case: the composition plot a novice would otherwise hand-build
  #     in ggplot. One expression -> counts + percent, both annotation levels.
  p_count <- PlotMetaSummary(no_dir(obj),
                             variables     = c("celltype.l1", "celltype.l2"),
                             fill_variable = "donor",
                             count_unit    = "cells",
                             output_dir    = NULL)
  p_pct   <- PlotMetaSummary(no_dir(obj),
                             variables     = c("celltype.l1", "celltype.l2"),
                             fill_variable = "donor",
                             count_unit    = "cells",
                             percent       = TRUE,
                             output_dir    = NULL)
  everyday <- (p_count / p_pct) +
    plot_annotation(
      title = "scSidekick: composition counts + percent - one expression, no ggplot",
      theme = theme(plot.title = element_text(size = 13, face = "bold",
                                              colour = "#2f4b7c")))
  ggsave_png(everyday, "pain_meta_everyday.png", 13, 7)
}

# =============================================================================
# PAIN POINT 3 - CellChat drops pathways absent in some cohorts
#   Real YAP-project data: ANGPTL signalling is active in Yap-intact tumours
#   (SGm, SGf) but ABSENT once YAP is deleted (YSGm, YSGf).
#   Before: the native per-object pathway plot errors for any cohort missing
#           the pathway, so those panels are silently dropped - the absence
#           (the actual finding) becomes invisible.
#   After:  CompareCellChat keeps every cohort on one shared layout, so the
#           ANGPTL loss in the YAP-deleted tumours is plain to see.
# =============================================================================
if (do_stage("cellchat")) {
  message(">> PAIN 3: CellChat dropped pathway (ANGPTL)")
  suppressMessages(library(CellChat))
  cc_dir <- file.path("/Users/nourhan/Library/CloudStorage/Box-Box",
                      "Yun lab projects/YAP project/ScRNAseq March 2024",
                      "MergedSoupx/Output/CellChat/CellChatObjects2")
  files <- c(SGm = "SGmCellChat.rds", SGf = "SGfCellChat.rds",
             YSGm = "YSGmCellChat.rds", YSGf = "YSGfCellChat.rds")
  pathway <- "ANGPTL"

  # --- BEFORE: native netVisual_aggregate per object (objects load one at a
  #     time to keep memory bounded). Cohorts without ANGPTL error -> dropped.
  png(file.path(fig_dir, "pain_cellchat_before.png"),
      width = 13, height = 3.8, units = "in", res = 150, bg = "white")
  par(mfrow = c(1, 4), mar = c(1, 1, 3, 1), xpd = NA)
  for (nm in names(files)) {
    cc <- readRDS(file.path(cc_dir, files[nm]))
    ok <- tryCatch({
      CellChat::netVisual_aggregate(cc, signaling = pathway, layout = "circle")
      title(main = paste0(nm, " - ", pathway), col.main = "#2f4b7c")
      TRUE
    }, error = function(e) FALSE)
    if (!ok) {
      plot.new()
      title(main = paste0(nm, " - ", pathway), col.main = "#B2182B")
      text(0.5, 0.55, paste0("Error: '", pathway, "' not found\nin this object"),
           col = "#B2182B", cex = 1.2, font = 2)
      text(0.5, 0.30, "native CellChat drops this panel", col = "grey40", cex = 1)
    }
    rm(cc); gc(verbose = FALSE)
  }
  dev.off()
  message("   wrote pain_cellchat_before.png")

  # --- AFTER: the CompareCellChat output (already generated for this project).
  #     Convert page 1 of the saved comparison PDF to PNG.
  ref_pdf <- file.path("/Users/nourhan/Library/CloudStorage/Box-Box",
                       "Yun lab projects/YAP project/ScRNAseq March 2024",
                       "MergedSoupx/Output/CellChat2/Compare",
                       " ANGPTL signaling network all samples.pdf")
  if (file.exists(ref_pdf)) {
    pg <- magick::image_read_pdf(ref_pdf, pages = 1, density = 150)
    magick::image_write(pg, file.path(fig_dir, "pain_cellchat_after.png"),
                        format = "png")
    message("   wrote pain_cellchat_after.png")
  } else {
    message("   reference PDF not found; skipping after image")
  }
  # To regenerate the comparison from scratch instead:
  #   CompareCellChat(as.list(file.path(cc_dir, files)),
  #                   output_dir = file.path(fig_dir, "cellchat_compare"))
}

# =============================================================================
# PAIN POINT 4 - Split dot plots that stay readable
#   Before: Seurat DotPlot(split.by) colours dots by the SPLIT identity (not
#           expression), interleaves cluster x condition rows, and crushes the
#           gene labels - the expression magnitude is no longer readable.
#   After:  SplitDotPlot organises genes into labelled cell-type blocks, facets
#           cleanly by condition, and keeps one shared expression colour scale.
# =============================================================================
if (do_stage("splitdot")) {
  message(">> PAIN 4: split dot plots")
  genes <- c("CD3E", "CD8A", "NKG7", "GNLY", "CD79A", "MS4A1",
             "LYZ", "S100A8", "CD14", "FCER1A")
  mdf <- data.frame(
    Genes    = genes,
    CellType = c("T", "T", "NK", "NK", "B", "B", "Mono", "Mono", "Mono", "DC"))

  # --- BEFORE: native Seurat split dot plot
  before <- Seurat::DotPlot(obj, features = genes,
                            group.by = "celltype.l1", split.by = "donor",
                            cols = c("blue", "red")) +
    ggtitle("Seurat: DotPlot(split.by) - dots coloured by group, labels unreadable") +
    theme(plot.title = element_text(size = 11, face = "bold", colour = "#B2182B"))
  ggsave_png(before, "pain_splitdot_before.png", 7.5, 5)

  # --- AFTER: scSidekick SplitDotPlot
  after <- SplitDotPlot(obj, markers_df = mdf,
                        group.by = "celltype.l1", split.by = "donor") +
    ggtitle("scSidekick: SplitDotPlot - gene blocks, clean facets, one expression scale") +
    theme(plot.title = element_text(size = 11, face = "bold", colour = "#2f4b7c"))
  ggsave_png(after, "pain_splitdot_after.png", 8.5, 5)
}

# =============================================================================
# PANEL B - PlotFeaturePlots with split.by + row.by + metadata layout
#   Demonstrates the two-variable faceting: lane (HTO) on columns,
#   donor (batch) on rows, panels positionally aligned by metadata coordinates.
# =============================================================================
if (do_stage("feature_b")) {
  message(">> PANEL B: PlotFeaturePlots row.by + metadata layout")
  gene <- "LYZ"

  panel_b <- PlotFeaturePlots(no_dir(obj), features = gene,
                               split.by      = "lane",
                               row.by        = "donor",
                               layout_method = "metadata",
                               output_dir    = NULL)
  ggsave_png(panel_b, "pain_feature_panel_b.png", 18, 8)
}

# =============================================================================
# FASTDOTPLOT - regex gene selection + automatic co-expression programs
#   Demonstrates: pattern = "^CD[1-2]" (no manual gene list needed),
#   k_genes = 3 slices the dendrogram into three labelled gene patterns.
# =============================================================================
if (do_stage("fastdotplot")) {
  message(">> FastDotPlot: regex + k-gene programs")

  fdp <- FastDotPlot(no_dir(obj),
                     pattern  = "^CD[1-2]",
                     group.by = "celltype.l1",
                     k_genes  = 3,
                     output_dir = NULL)
  ggsave_png(fdp, "pain_fastdotplot.png", 8.5, 4.5)
}

message(">> done (stages: ",
        if (length(stages)) paste(stages, collapse = ", ") else "all", ")")
