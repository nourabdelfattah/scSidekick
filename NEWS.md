# scSidekick 1.1.0

## Pathway analysis

- `RunGSEA()` NES heatmaps now show significance stars (`*`/`**`/`***` by
  `padj`), previously computed but never drawn.
- `RunGSEA()`, `RunGSEA_pseudobulk()`, and `RunSCssGSEA()` now write every
  run into its own timestamped subfolder under `output_dir` by default
  (`timestamp = TRUE`), with a `resume` / `resume_folder` system that scans
  for a previous run with matching parameters and reuses it - `resume = TRUE`
  prompts when several match, `resume = "last"` picks the most recent
  automatically. The literal function call (as typed, with argument values)
  is also now recorded in each run's methods JSON.
- `RunGSEA_pseudobulk()` gains `use.padj` to choose raw p vs. BH-adjusted p
  as both the significance-star basis and the heatmap-selection cutoff for
  its main NES heatmaps (the ssGSEA sub-step already had this); legend text
  now correctly reflects whichever basis was actually used instead of always
  saying "BH-adjusted."
- `RunGSEA_pseudobulk()` adds a design-rank diagnostic
  (`stats::alias()`-based) that explains *why* a given `covariates` +
  `contrast.by`/`group.by` combination produces "coefficients not
  estimable" - typically a confound between a covariate and the contrast,
  not a bug.
- `RunSCssGSEA()`:
  - `search_terms` no longer gets silently bypassed when `gene_sets` is
    supplied manually - it now filters custom gene set names too, matching
    `RunGSEA()` / `RunGSEA_pseudobulk()` behavior.
  - New `pathway_sets` argument fetches and pools **multiple** MSigDB
    collections in one run (the same named-list format `RunGSEA()` uses),
    instead of being limited to a single `gene_set_library` /
    `gene_set_subcategory` pair.
  - New `assay` argument (default `"RNA"`) makes the function actually
    BPCells/sketch-aware as documented - previously hardcoded to `"RNA"`
    with no way to score a Seurat v5 sketch-subsampled assay.
  - New `use.padj` toggle and `pathways` argument (manual pathway override,
    with UCell `"_UCell"`-suffix-tolerant matching) for the significance
    heatmap.
  - Corrected boxplot legend text: pairwise p-value brackets are raw
    (unadjusted) t-test p-values via `ggpubr`, not BH-corrected as previously
    (incorrectly) stated.

## Cell-cell communication

- `RunCellChat()`'s per-pathway multi-panel PDF gains a gene-level
  ligand-receptor chord diagram (`CellChat::netVisual_chord_gene()`) showing
  individual ligand/receptor genes by sender/receiver cell type, alongside
  the existing pathway-aggregate chord diagram.

## Visualization

- `PlotMetaSummary()` Table 1 output now includes a table-specific legend
  (in-table footnote + `.legend` sidecar) stating exactly which statistical
  test was used per variable type, instead of an unexplained bare p-value
  column. New `numeric_stat`, `cat_display`, and `separate_tables` arguments.
- `GroupHeatmap()` gains a `timestamp` argument to version repeated PDF saves
  instead of overwriting; `legend_side` now moves both the color legend and
  annotation legends together.
- `PlotFeatureTrend()` gains `donor.by` (collapse to one point per donor
  before fitting - the honest way to relate a donor-level variable to a
  per-cell value), `add_cor` / `cor_method` (correlation coefficient and
  p-value annotation via `ggpubr::stat_cor()`), and `show_points` /
  `point_size` / `point_alpha` (overlay raw observations and train the
  y-axis on the full data range).

## Bug fixes

- Fixed a docs/export regression in `R/pseudobulk.R` where a private helper
  had absorbed `PlotPseudoBulk()`'s roxygen block, causing `PlotPseudoBulk()`
  to silently lose its `@export` and documentation.
- `PlotFeaturePlots()`, `PlotMasterMaps()`, `PlotSpatialDimPlots()`, and
  `PlotSpatialFeaturePlots()` are the primary documented names again (`?`
  now resolves directly for each), with `GenerateFeatureMaps()`,
  `GenerateMasterGeneMaps()`, `GenerateSpatialDimMaps()`, and
  `GenerateSpatialFeatureMaps()` kept as long-form aliases - matching how
  they're actually referenced everywhere else in the package (README,
  vignettes). Previously the doc pages were filed under the long-form
  names, and a stray `man/.gitignore` had excluded 3 of the 4 doc pages
  from git entirely.

## Site

- Added a Google Search Console verification tag to the pkgdown site config.

# scSidekick 1.0.0

Initial release.

## Project setup

- `PrepObject()` stores project-level settings (output directory, object name,
  assay type, donor column, colors) directly on the Seurat object so every
  downstream function can resolve them automatically.
- `GenerateDirectories()` creates a standard analysis folder tree.
- `LoadSamplesRNA()` + `PlotQCMetrics()` load and QC-plot one or more 10x
  libraries, with doublet detection (scDblFinder) and ambient RNA filtering
  (miQC).

## Dimensionality reduction and UMAP

- `Determine_nDims()` selects the number of PCs via elbow or JackStraw.
- `PlotDimPlots()` batch-generates labeled and split UMAP PDFs.
- `PlotGridUMAP()` arranges multiple feature UMAPs in a grid.

## Visualization

- `FastDotPlot()` / `FastDotPlot2()` — dot plots with flexible grouping and
  split panels.
- `SplitDotPlot()` / `SplitDotPlot2()` — side-by-side dot plots across
  conditions.
- `StackedVlnPlot()` — compact stacked violin plot.
- `GroupHeatmap()` — pseudobulk mean-expression heatmap by group.
- `GenerateFeatureMaps()` / `GenerateMasterGeneMaps()` / `PlotMultiFeature()`
  — batch feature expression maps.
- `PlotVolcano()` — styled volcano plot for differential expression results.

## Spatial transcriptomics

- `GenerateSpatialDimMaps()` and `GenerateSpatialFeatureMaps()` produce
  publication-ready spatial plots for Visium and VisiumHD data with correct
  spot/bin terminology.

## Composition and metadata

- `PlotComposition()` — stacked bar charts of cell type composition.
- `PlotPieUMAP()` — pie-chart overlays on UMAP embeddings.
- `PlotAlluvial()` — alluvial / Sankey diagrams across conditions.
- `PlotRose()` — rose (polar bar) charts.
- `PlotChord()` — chord diagrams via circlize.
- `PlotTrendAndUMAP()` / `PlotTrendLabeled()` — composition trend plots.
- `PlotAtlasWheel()` — circular dataset overview wheel.
- `SummarizeMetadata()` / `ChartBuilder()` — flexible metadata summary charts.

## Cell annotation

- `CellTypeAssignmentHelper()` — supervised annotation with SingleR + manual
  override, writes labeled UMAP PDFs.
- `AnnotateFeatures()` — marker-gene-based annotation overlay.
- `CheckSex()` — infers biological sex from XIST / Y-chromosome expression.

## Pathway analysis

- `RunGSEA()` — GSEA on cluster markers using fgsea; supports msigdbr and
  enrichR databases; writes lollipop and NES heatmap PDFs.
- `RunGSEA_pseudobulk()` — pseudobulk GSEA with edgeR DE; lollipop, NES
  heatmap, and ssGSEA scoring outputs.
- `RunSCssGSEA()` — single-cell ssGSEA scoring with heatmap and box plot.
- `VisualizeLeadingEdge()` — heatmap of leading-edge genes for selected
  pathways.
- `PlotPathwayButterfly()` — two-axis cell state hierarchy scatter plot
  (adapted from Suva Lab scrabble scorer).
- `RunEnrichment()` / `PlotEnrichment()` — over-representation analysis.
- `PlotGSEAEnrichment()` — enrichment curve plots for individual pathways.

## CellChat

- `RunCellChat()` — end-to-end CellChat pipeline with circle plots, bubble
  plots, and communication pattern outputs.
- `CompareCellChat()` — differential interaction comparison across conditions.
- `RankCellChatPathways()` — pathway-level signaling strength ranking.
- `RenameCellTypeInCC()` — rename a cell type across a CellChat list.
- `PlotPctHeatmap()` — percent-expressing heatmap for receptor/ligand genes.

## Copy number variation

- `RunCopyKAT()` — CopyKAT CNV inference with result integration onto the
  Seurat object.
- `RunInferCNV()` — inferCNV pipeline wrapper.

## Coexpression

- `CorrelateGene()` + `PlotCorrelation()` — gene-gene and gene-metadata
  correlation analysis.

## Consensus NMF

- `RunCNMF()` — cNMF factorization over a K range via reticulate; writes
  k-selection diagnostic PNG.
- `GetCNMFPrograms()` — consensus solution for a chosen K; stores usage matrix
  as a DimReduc and gene spectra in `@misc$cnmf`.
- `GetCNMFTopGenes()` — extract top genes per cNMF program.

## Upstream regulator analysis

- `RunURA()` — upstream regulator analysis on DE gene lists.

## Survival analysis

- `SurvPlot()` — Kaplan-Meier survival analysis for bulk RNA-seq or TCGA data;
  supports gene-level and gene-set-signature stratification, median/quartile/
  optimal splitting, row x column faceting by clinical metadata, and
  multivariate Cox regression.
- `LoadTCGA()` / `EnrichTCGAMeta()` — download and annotate TCGA cohorts via
  recount3 or TCGAbiolinks.
- `SurvMetaSummary()` — tabular summary of clinical metadata for a cohort.

## Reporting and reproducibility

- `ExtractMethods()` — auto-generates a Methods paragraph from settings stored
  on the Seurat object via `log_analysis_params()`.
- `log_analysis_params()` / `log_figure_legend()` — write analysis parameters
  and figure legends to JSON/text files alongside output PDFs.
- `create_analysis_pptx()` — assembles a PowerPoint summary deck from a folder
  of PDFs; supports section filtering and aspect-ratio-preserving layout.
- `start_session_logger()` / `summarize_r_session()` — passive session logging
  and post-hoc markdown + R-script summarizer.
- Every plotting function writes a `.txt` legend sidecar alongside each output
  PDF containing a full figure legend with observation counts (cells, nuclei,
  spots, or bins), donor/sample counts, and dataset context.

## Colors and theme

- `Nour_pal()`, `scale_color_Nour()`, `scale_fill_Nour()` — custom color
  palettes (`Nour18`, `Nour20`) with ggplot2 scale functions.
- `theme_NourMin()` — minimalist ggplot2 theme.
- `ShowColors()` / `GetColors()` — inspect and retrieve per-object color
  assignments.
