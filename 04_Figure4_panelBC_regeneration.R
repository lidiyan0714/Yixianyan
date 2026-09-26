# =============================================================================
# 04_Figure4_panelBC_regeneration.R
#
# PURPOSE
#   Regenerates Figure 4 panels B and C (ROC curves) from Table S1 and
#   Table S2, using the authors' real, fixed model coefficients (see
#   03_Figure4_script_and_output.R and
#   Etiology_Results_authors_original_output.xlsx / the ABP_vs_HTGP_Results.xlsx
#   Sensitivity_Coefficients sheet). It reproduces the apparent (in-sample,
#   non-optimism-corrected) ROC curves and AUCs exactly as they are drawn in
#   Figure 4B and 4C.
#
# APPARENT vs. OPTIMISM-CORRECTED AUC
#   Figure 4B/C display each ROC curve's own apparent AUC, which is by
#   definition what a curve computed on the full analysis sample evaluates
#   to. The manuscript's headline optimism-corrected AUCs (Main 0.945,
#   Sensitivity 0.849; see 03_Figure4_script_and_output.R)
#   are a separate, lower, bootstrap-adjusted internal-validation estimate
#   and are not meant to appear as a curve's own legend value. This script
#   only regenerates the apparent-AUC curves shown in the figure.
#
# CURRENT INPUT DATA
#   Model inputs are taken from Table S1's "Initial Triglyceride Level (for
#   etiology determination)" column (the field the manuscript uses to define
#   the ABP/HTGP groups) and Table S2's raw (untransformed) metabolite
#   values. With the current Table S1/S2, this script reproduces:
#     Individual Triglyceride-only apparent AUC ......... 0.881
#     Main/Combined model apparent AUC ................... 0.955026
#     Sensitivity model apparent AUC ...................... 0.895503
#   matching the values currently shown in Figure 4B/C and cited in Results.
#
# OUTPUT
#   panelB_vec.pdf / panelC_vec.pdf -- replacement panels for Figure 4, in
#   the same base-R plotting style as the authors' own original scripts
#   (see 03_Figure4_script_and_output.R). Each is written
#   as a native VECTOR PDF (real embedded font glyphs, not a rasterized image)
#   sized to exactly
#   the panel-B / panel-C slot on the Figure 4 page (194.35 x 209.51 pt), so
#   it can be spliced into the original vector Figure 4.pdf with
#   show_pdf_page() while keeping all text/fonts editable/selectable.
#
# REQUIREMENTS: R packages readxl, pROC
#   install.packages(c("readxl", "pROC"))
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(pROC)
})

DATA_DIR <- "."   # <- point this at the current "Figures and tables" folder.

find_input <- function(candidates) {
  paths <- file.path(DATA_DIR, candidates)
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0L) {
    stop("Input file not found. Expected one of: ", paste(candidates, collapse = ", "))
  }
  existing[[1]]
}

table_s1_path <- find_input(c("Table S1.xlsx", "Table S1 Clinical data collection.xlsx"))
table_s2_path <- find_input(c("Table S2.csv", "Table S2 Total Metabolite Statistics.csv"))

# ---- Table S1 (clinical) ----
raw1 <- read_excel(table_s1_path, col_names = FALSE)
hdr1 <- as.character(unlist(raw1[3, ]))
s1 <- raw1[-(1:3), ]
names(s1) <- hdr1
s1 <- as.data.frame(s1)

etio_col <- grep("^Etiological Type", names(s1), value = TRUE)[1]
s1[[etio_col]] <- as.numeric(s1[[etio_col]])
s1 <- s1[s1[[etio_col]] %in% c(1, 2), ]
stopifnot(nrow(s1) == 55)

s1$y_abp <- ifelse(s1[[etio_col]] == 1, 1L, 0L)  # 1 = ABP (Biliary), 0 = HTGP
s1$Sample_norm <- toupper(trimws(as.character(s1[["Sample"]])))

tg_column <- "Initial Triglyceride Level (for etiology determination)"
if (!tg_column %in% names(s1)) stop("Table S1 does not contain the required Initial TG field: ", tg_column)
s1$TG      <- as.numeric(s1[[tg_column]])
s1$Age_v   <- as.numeric(s1[["Age"]])
s1$BMI_v   <- as.numeric(s1[["BMI"]])
s1$CRP_v   <- as.numeric(s1[["CRP"]])
s1$MCTSI_v <- as.numeric(s1[["MCTSI score"]])

# ---- Table S2 (metabolomics) ----
s2 <- read.csv(table_s2_path, header = FALSE, stringsAsFactors = FALSE)
hdr2 <- as.character(unlist(s2[2, ]))
s2 <- s2[-(1:2), ]
names(s2) <- hdr2

met1_row <- s2[s2$Metabolite == "2-(propylthio)nicotinic acid", ]
met2_row <- s2[grepl("^2-Amino-9-\\[4-hydroxy", s2$Metabolite), ]
stopifnot(nrow(met1_row) == 1, nrow(met2_row) == 1)

sample_cols <- names(s2)[toupper(names(s2)) %in% s1$Sample_norm]
get_vals <- function(row, cols) { v <- as.numeric(unlist(row[cols])); names(v) <- toupper(cols); v }
met1_vals <- get_vals(met1_row, sample_cols)
met2_vals <- get_vals(met2_row, sample_cols)
s1$met1 <- met1_vals[s1$Sample_norm]
s1$met2 <- met2_vals[s1$Sample_norm]

# ---- MAIN model: real coefficients from Etiology_Results_authors_original_output.xlsx ----
d <- s1[!is.na(s1$TG) & !is.na(s1$met1) & !is.na(s1$met2), ]
intercept <- -130.000000; b_tg <- -0.050981; b_met1 <- -0.728771; b_met2 <- 26.753259
lp_main <- intercept + b_tg*d$TG + b_met1*d$met1 + b_met2*d$met2
pred_main <- 1/(1+exp(-lp_main))
roc_main <- roc(d$y_abp, pred_main, quiet = TRUE)
cat("Reconstructed MAIN apparent AUC:", round(auc(roc_main), 6), "\n")

roc_tg   <- roc(d$y_abp, d$TG,   quiet = TRUE, direction = "auto")
roc_met1 <- roc(d$y_abp, d$met1, quiet = TRUE, direction = "auto")
roc_met2 <- roc(d$y_abp, d$met2, quiet = TRUE, direction = "auto")
cat("Individual AUCs:", round(auc(roc_tg),3), round(auc(roc_met1),3), round(auc(roc_met2),3), "\n")

# ---- SENSITIVITY model: real coefficients from ABP_vs_HTGP_Results.xlsx (Sensitivity_Coefficients) ----
d2 <- s1[!is.na(s1$TG) & !is.na(s1$Age_v) & !is.na(s1$BMI_v) & !is.na(s1$CRP_v) & !is.na(s1$MCTSI_v), ]
intercept_s <- 0.551137; b_age <- 0.048381; b_bmi <- 0.036894; b_tg_s <- -0.063464; b_crp <- -0.008235; b_mctsi <- -0.086275
lp_sens <- intercept_s + b_age*d2$Age_v + b_bmi*d2$BMI_v + b_tg_s*d2$TG + b_crp*d2$CRP_v + b_mctsi*d2$MCTSI_v
pred_sens <- 1/(1+exp(-lp_sens))
roc_sens <- roc(d2$y_abp, pred_sens, quiet = TRUE)
cat("Reconstructed SENSITIVITY apparent AUC:", round(auc(roc_sens), 6), "\n")

# ============================= PANEL B ======================================
# Individual-marker ROC curves (Combined model + the 3 individual markers),
# reproducing the authors' own original plotting code for this panel.
# Written as a native VECTOR PDF (real embedded font glyphs, not a rasterized
# image) sized to exactly the panel-B slot on the Figure 4 page (194.35 x
# 209.51 pt), so it can be spliced into the original vector Figure 4.pdf with
# show_pdf_page() while keeping all text/fonts editable/selectable.
pdf("panelB_vec.pdf", width = 194.34645707376058/72, height = 209.51253420752565/72, pointsize = 7)
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
par(mar = c(3.0, 3.0, 2.2, 0.6), mgp = c(1.7,0.55,0), tcl = -0.3, oma = c(0,0,0.9,0), cex.main = 1.05, font.main = 2)
plot(roc_main, col = colors[1], lwd = 1.6, lty = 1,
     main = "Individual ROC - ABP vs HTGP",
     xlab = "Specificity", ylab = "Sensitivity", cex.lab = 0.92, cex.axis = 0.85)
lines(roc_tg,   col = colors[2], lwd = 1.1, lty = 1)
lines(roc_met1, col = colors[3], lwd = 1.1, lty = 1)
lines(roc_met2, col = colors[4], lwd = 1.1, lty = 1)
legend("bottomright",
       legend = c(sprintf("Combined (AUC=%.3f)", auc(roc_main)),
                  sprintf("Triglyceride (AUC=%.3f)", auc(roc_tg)),
                  sprintf("2-(propylthio)nicotinic acid (AUC=%.3f)", auc(roc_met1)),
                  sprintf("2-Amino-9-[4-hydroxy-2-(hydroxymethyl)butyl]\n-3H-purin-6-one (AUC=%.3f)", auc(roc_met2))),
       col = colors, lwd = 1.4, cex = 0.46, bty = "o", seg.len = 1.1, y.intersp = 1.15)
mtext("(B)", side = 3, line = 0.15, at = par("usr")[1] - (par("usr")[2]-par("usr")[1])*0.30, cex = 1.15, font = 2, xpd = NA)
dev.off()

# ============================= PANEL C ======================================
# Main vs Sensitivity model ROC comparison, reproducing the authors' own
# original plotting code for this panel. Same vector-PDF approach as panel B,
# sized to the panel-C slot (identical dimensions to panel B).
pdf("panelC_vec.pdf", width = 194.34645707376058/72, height = 209.51253420752565/72, pointsize = 7)
par(mar = c(3.0, 3.0, 2.2, 0.6), mgp = c(1.7,0.55,0), tcl = -0.3, oma = c(0,0,0.9,0), cex.main = 1.05, font.main = 2)
plot(roc_main, col = "blue", lwd = 1.6,
     main = "ROC Curves - ABP vs HTGP",
     xlab = "Specificity", ylab = "Sensitivity", cex.lab = 0.92, cex.axis = 0.85)
lines(roc_sens, col = "red", lwd = 1.4)
legend("bottomright",
       legend = c(sprintf("Main (AUC %.3f)", auc(roc_main)),
                  sprintf("Sensitivity (AUC %.3f)", auc(roc_sens))),
       col = c("blue","red"), lwd = 1.4, cex = 0.65, bty = "o")
mtext("(C)", side = 3, line = 0.15, at = par("usr")[1] - (par("usr")[2]-par("usr")[1])*0.30, cex = 1.15, font = 2, xpd = NA)
dev.off()

cat("\nSaved panelB_vec.pdf and panelC_vec.pdf\n")
