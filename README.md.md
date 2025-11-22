# Baron et al. 2016 pancreas scRNA-seq reproduction

This repository contains R scripts to reproduce the main analyses and figures
(Fig.1–4, human and mouse) from **Baron et al., 2016, Cell Metabolism**  
(GEO: **GSE84133**, human & mouse pancreas UMI counts).

The pipeline is implemented in base R + Seurat and is split into four scripts.

---

## 1. Data

Download the raw UMI count matrices (`*.csv.gz`) from GEO series **GSE84133**.

Place all files under a directory, for example:

```text
GSE84133_RAW/
  GSM2230757_human1_UMIcounts.csv.gz
  GSM2230758_human2_UMIcounts.csv.gz
  ...
  GSM2230767_mouse1_UMIcounts.csv.gz
  ...
```

File names should contain the strings **"human"** or **"mouse"**, so the scripts
can automatically separate human and mouse samples.

You can place this project anywhere on your machine, for example:

```text
C:/Users/<username>/Documents/Baron2016_pancreas_scRNA_repro
```

Inside that folder, create a subfolder and put all raw files under:

```text
C:/Users/<username>/Documents/Baron2016_pancreas_scRNA_repro/GSE84133_RAW
```

Then set `PROJECT_DIR` in each script (if needed) to this project root.

---

## 2. Dependencies

Main R packages:

- R (>= 4.3)
- Seurat (>= 5.0)
- data.table
- Matrix
- dplyr
- ggplot2
- pheatmap
- patchwork
- tidyr

Install in R:

```r
install.packages(c(
  "Seurat", "data.table", "Matrix", "dplyr",
  "ggplot2", "pheatmap", "patchwork", "tidyr"
))
```

If you use a personal library path, set it at the top of your R session, for example:

```r
dir.create("D:/Rlibs_new", showWarnings = FALSE)
.libPaths(c("D:/Rlibs_new", .libPaths()))
```

---

## 3. Pipeline overview

Assume your project root is:

```text
C:/Users/<username>/Documents/Baron2016_pancreas_scRNA_repro
```

and you have this structure:

```text
Baron2016_pancreas_scRNA_repro/
  ├── README.md
  ├── .gitignore
  ├── GSE84133_RAW/
  └── scripts/
        ├── 01_build_seurat.R
        ├── 02_fig1_markers.R
        ├── 03_fig2_corr.R
        └── 04_fig3_ductal_fig4_beta_PCA.R
```

All scripts should be executed from the **project root** (working directory):

```r
setwd("C:/Users/<username>/Documents/Baron2016_pancreas_scRNA_repro")
```

or using command line (if you run via `Rscript`):

```bash
cd C:/Users/<username>/Documents/Baron2016_pancreas_scRNA_repro
```

### Step 1 — Build Seurat objects & annotate cell types

Script: `scripts/01_build_seurat.R`

This script:

- reads all `GSE84133_RAW/*.csv.gz`,
- constructs human and mouse Seurat objects,
- performs QC (nFeature, nCount, percent.mt, per-sample MAD filters),
- normalizes, finds HVGs, PCA, neighbors, clustering, UMAP/tSNE,
- annotates major pancreas cell types (Alpha, Beta, Delta, PP, Ductal, Acinar,
  Endothelial, Stellate, Immune, Unknown) using canonical marker scores,

and writes:

- `results_human_paper/Human_seurat_celltyped.rds`
- `results_mouse_paper/Mouse_seurat_celltyped.rds` (if mouse is available)
- `results_human_paper/run_metadata_build_seurat.rds` (config + session info)

Run in R:

```r
setwd("C:/Users/<username>/Documents/Baron2016_pancreas_scRNA_repro")
source("scripts/01_build_seurat.R")
```

or via command line:

```bash
Rscript scripts/01_build_seurat.R
```

---

### Step 2 — Fig.1 UMAP/tSNE + canonical marker dotplots + top-marker heatmaps

Script: `scripts/02_fig1_markers.R`

This script generates:

**Human**

- `results_human_paper/Human_UMAP_celltype_paper.pdf`
- `results_human_paper/Human_tSNE_celltype_paper.pdf`
- `results_human_paper/Human_DotPlot_markers_paper.pdf`
- `results_human_paper/Human_TopMarkers_Heatmap_paper.pdf`
- `results_human_paper/Markers_Human_byCelltype.csv`

**Mouse** (if present)

- `results_mouse_paper/Mouse_UMAP_celltype_paper.pdf`
- `results_mouse_paper/Mouse_tSNE_celltype_paper.pdf`
- `results_mouse_paper/Mouse_DotPlot_markers_paper.pdf`
- `results_mouse_paper/Mouse_TopMarkers_Heatmap_paper.pdf`
- `results_mouse_paper/Markers_Mouse_byCelltype.csv`

Run:

```r
setwd("C:/Users/<username>/Documents/Baron2016_pancreas_scRNA_repro")
source("scripts/02_fig1_markers.R")
```

---

### Step 3 — Fig.2 donor × cell type correlation heatmaps

Script: `scripts/03_fig2_corr.R`

This script:

- computes average log-normalized expression for each **donor × cell type**
  (using `sample` and `celltype_final`),
- calculates Pearson correlations between donor–celltype profiles,
- plots clustered correlation heatmaps with donor/celltype annotations.

Outputs:

- `results_human_paper/Human_Fig2A_donor_celltype_corr_heatmap.pdf`
- `results_mouse_paper/Mouse_Fig2C_donor_celltype_corr_heatmap.pdf`
- `results_human_paper/Human_Mouse_Fig2_commonCT_corr_heatmap.pdf`
  (human + mouse common cell types)

Run:

```r
setwd("C:/Users/<username>/Documents/Baron2016_pancreas_scRNA_repro")
source("scripts/03_fig2_corr.R")
```

---

### Step 4 — Fig.3 ductal PCA & Fig.4 beta PCA / PC1 profiles

Script: `scripts/04_fig3_ductal_fig4_beta_PCA.R`

This script:

- subsets **Ductal** cells and runs PCA,
  - Fig.3A: PCA plots colored by canonical ductal markers
  - Fig.3B: moving average of marker expression along PC1
  - Fig.3C: heatmap of genes with highest PC1 loadings
- subsets **Beta** cells and runs PCA,
  - Fig.4A: PC1 vs PC2 colored by `nCount_RNA`
  - Fig.4B / 4C: moving average & heatmap along PC1 for ER-stress / β-identity
    genes

Outputs (human):

- `results_human_paper/Human_Fig3A_Ductal_PCA_PC1_PC2.pdf`
- `results_human_paper/Human_Fig3B_Ductal_PC1_Markers_MA.pdf`
- `results_human_paper/Human_Fig3C_Ductal_PC1_loading_heatmap.pdf`
- `results_human_paper/Human_Fig4A_Beta_PCA_PC1_PC2_nCount.pdf`
- `results_human_paper/Human_Fig4B_Beta_PC1_Markers_MA.pdf`
- `results_human_paper/Human_Fig4C_Beta_PC1_loading_heatmap.pdf`

Outputs (mouse, if data present) with analogous filenames under
`results_mouse_paper/`.

Run:

```r
setwd("C:/Users/<username>/Documents/Baron2016_pancreas_scRNA_repro")
source("scripts/04_fig3_ductal_fig4_beta_PCA.R")
```

---

## 4. Reproducibility notes

- All scripts set a global random seed (`set.seed(123)`).
- QC thresholds (min/max genes, percent.mt, MAD trimming) are recorded in
  `run_metadata_build_seurat.rds`.
- To adapt to different machines (laptop / server), change `PROJECT_DIR` at the
  top of each script to your own path.
- The code has been tested on Linux + R 4.3 + Seurat 5.0. On Windows, make sure
  you have a recent R version and sufficient memory for the full dataset.
