# Yixianyan (Pinyin for "pancreatitis")

Analysis code for: *Integrated metabolomics and single-cell transcriptomics
characterize etiology-specific immune-metabolic alterations in early acute
pancreatitis*

This repository contains exactly the code used to produce the manuscript's
reported numeric results for Figure 4 (single-marker and combined-model
ABP-vs-HTGP classification AUCs) and Figure 5C / Table S10 (per-cluster
differential-expression summaries). Each script's own header documents the
values it reproduces and how they were verified against the manuscript.

## Contents

- **`01_single_marker_AUC_reproduction.R`** -- reproduces the three
  single-marker AUCs reported in Figure 4B (Triglyceride, 2-(propylthio)
  nicotinic acid, and the third metabolite) directly from Table S1 and
  Table S2, to 6 decimal places.

- **`02_Figure5C_and_TableS10_regeneration.R`** -- regenerates Figure 5C /
  Figure S12A and Table S10 from the per-cluster Seurat `FindMarkers`
  differential-expression output (ABP vs Control and HTGP vs Control),
  using one consistent filter (p_val_adj < 0.01 & |avg_log2FC| > 0.25) for
  both the figure and the table.

- **`03_Figure4_script_and_output.R`** -- the authors' own
  original R script that was actually used to fit the Figure 4 combined
  and sensitivity classification models (LASSO variable selection +
  Firth-penalized logistic regression + bootstrap internal validation).
  Comments have been translated from Chinese to English and local file
  paths replaced with generic placeholders; the analysis logic is
  otherwise unchanged. Its real, previously-generated output workbook is
  included alongside it as `Etiology_Results_authors_original_output.xlsx`.

- **`04_Figure4_panelBC_regeneration.R`** -- regenerates Figure 4 panels B
  and C (the ROC curves) directly from Table S1 and Table S2, using the
  real, fixed model coefficients recovered from script 03's saved output.
  Reproduces the apparent AUCs shown in the figure: Combined model
  0.955026, Triglyceride 0.881, 2-(propylthio)nicotinic acid 0.724, third
  metabolite 0.804, sensitivity model 0.895503.

## Data

Raw metabolomic and single-cell RNA-seq data are deposited in the NCBI SRA
database under accession number PRJNA1203495. Processed data tables
(Table S1-S10) are provided as Supplementary Materials with the
manuscript. Each script's header states which supplementary table(s) it
expects as input and the exact column names it reads.

## Requirements

Each script's header comment documents the R packages it requires. In
general: `Seurat`, `glmnet`, `logistf`, `pROC`, `boot`, `readxl`, `dplyr`,
`tidyr`, `ggplot2`, `openxlsx`, `ResourceSelection`.

## License

MIT. See `LICENSE`.
