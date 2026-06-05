```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           ___   _      _         _  __  _        _   
 ___  __  / __| (_)  __| |  ___  | |/ / (_)  __  | |__
(_-< / _| \__ \ | | / _` | / -_) | ' <  | | / _| | / /
/__/ \__| |___/ |_| \__,_| \___| |_|\_\ |_| \__| |_\_\

 v0.1.0  ·  Your New Best Friend in Visualization

 ✨  It's a good day to make pretty figures!  ✨
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

**scSidekick** is a personal single-cell RNA-seq and spatial transcriptomics
toolkit built on top of [Seurat](https://satijalab.org/seurat/). It wraps
repetitive boilerplate into one-liners, enforces consistent colors and layouts
across an entire project, and automates the last mile from analysis to slide deck.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("nourabdelfattah/scSidekick")
```

---

## What's inside

| Category | Key functions |
|---|---|
| **Setup & Colors** | `PrepObject()`, `ShowColors()`, `GetColors()`, `Nour_pal()` |
| **Dimensionality Reduction** | `Determine_nDims()` |
| **UMAP & Composition** | `PlotDimPlots()`, `PlotTrendLabeled()`, `PlotComposition()` |
| **Feature Expression** | `PlotFeaturePlots()`, `PlotMultiFeature()` |
| **Heatmaps & Dot Plots** | `GroupHeatmap()`, `SplitDotPlot()`, `FastDotPlot()` |
| **Cell Annotation** | `CellTypeAssignmentHelper()` |
| **Spatial** | `PlotSpatialDimPlots()`, `PlotSpatialFeaturePlots()`, `PlotMasterMaps()` |
| **Pathway Analysis** | `RunGSEA()` |
| **Cell-Cell Communication** | `RunCellChat()`, `CompareCellChat()` |
| **QC & Loading** | `PlotQCMetrics()`, `LoadSamplesRNA()` |
| **Reporting** | `create_analysis_pptx()`, `ExtractMethods()`, `theme_NourMin()` |

---

## Quick start

```r
library(scSidekick)
library(Seurat)

# After clustering, register colors and defaults once:
SeuratObj <- PrepObject(SeuratObj,
  variables   = c("Sample", "Group", "seurat_clusters"),
  group.by    = "seurat_clusters",
  split.by    = "Sample",
  output_dir  = "./Figures",
  object_name = "MyProject")

# Every downstream function reads those defaults automatically:
PlotDimPlots(SeuratObj)               # UMAP split by Sample, colored by cluster
GroupHeatmap(SeuratObj, features = marker_genes)
PlotFeaturePlots(SeuratObj, features = c("CD3E", "CD79A", "LYZ"))
RunGSEA(SeuratObj, group.by = "seurat_clusters", output_dir = "./Figures/Pathways")
create_analysis_pptx(output_dir = "./Figures", object_name = "MyProject")
```

---

## Documentation

Full vignettes and function reference at
**[nourabdelfattah.github.io/scSidekick](https://nourabdelfattah.github.io/scSidekick/)**

| Vignette | Topic |
|---|---|
| 01 | PrepObject, ShowColors, GetColors |
| 02 | Determine_nDims |
| 03 | PlotDimPlots — UMAP visualization |
| 04 | PlotFeaturePlots, PlotMultiFeature |
| 05 | GroupHeatmap, SplitDotPlot, FastDotPlot |
| 06 | CellTypeAssignmentHelper |
| 07 | Spatial transcriptomics |
| 08 | RunGSEA — pathway analysis |
| 09 | RunCellChat, CompareCellChat |
| 10 | Reporting utilities |
