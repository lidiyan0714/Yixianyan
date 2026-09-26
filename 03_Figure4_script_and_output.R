# =============================================================================
# 03_Figure4_script_and_output.R
#
# WHAT THIS FILE IS
#   This is the AUTHORS' OWN original R script that was actually used to fit
#   the Figure 4 / Table (Model construction and validation) ABP-vs-HTGP
#   etiology classification model reported in the manuscript. It was located
#   on the corresponding author's own computer, in the working folder used
#   during model development, together with its original saved output
#   workbook (included alongside this script as
#   "Etiology_Results_authors_original_output.xlsx"). For public release,
#   the analysis logic below is completely unchanged from the original;
#   only its comments have been translated from Chinese to English and its
#   local file paths (Section 2 below) have been replaced with generic
#   placeholders.
#
#   This is the authoritative record of how the Figure 4 model was actually
#   built: it is the authors' own real code, not a reconstruction from the
#   Methods text.
#
# HOW WE KNOW THIS IS THE RIGHT SCRIPT
#   Running this exact pipeline previously produced a saved output workbook
#   whose three selected predictors are EXACTLY the three predictors named
#   in the manuscript's Figure 4 / Results text:
#     - Triglyceride
#     - 2-(propylthio)nicotinic acid
#     - 2-Amino-9-[4-hydroxy-2-(hydroxymethyl)butyl]-3H-purin-6-one
#   No other script found during this search reproduces all three predictor
#   names simultaneously, which is why this one is identified as the true
#   source of the published model (as opposed to another candidate script
#   found in the same working folder, which uses a different, stability-
#   selection-based pipeline and selects only "Triglyceride" as its
#   ABP_vs_HTGP predictor in its saved output).
#
# PIPELINE AS ACTUALLY RUN (unchanged from the original script below)
#   1. Candidate pool = 21 clinical variables (Age, Onsettime, BMI,
#      Temperature, Triglyceride, Balthazar, RansonCT, CTSI, MCTSI,
#      FibrinogenDegradationProducts, D_Dimer, CRP, PCT, ALT, AST, ALP, GGT,
#      Totalbilirubin, Directbilirubin, creatinin_e, GCS) PLUS every
#      remaining numeric metabolite column in the source workbook
#      (Metab_data.xlsx) -- i.e. clinical variables and metabolites were
#      pooled together and screened jointly, not pre-filtered by univariate
#      screening before LASSO.
#   2. Missing metabolite values -> median imputation; zero-variance columns
#      dropped; all candidate variables z-score standardized (scale()).
#   3. Variable selection: LASSO logistic regression (glmnet, alpha = 1,
#      10-fold CV, lambda.min, set.seed(123)).
#   4. If LASSO selects more than 3 variables, only the top 3 BY ABSOLUTE
#      LASSO COEFFICIENT MAGNITUDE are retained for the final model
#      (Triglyceride was selected by LASSO itself, on equal footing with the
#      metabolites -- it was not forced into the model separately).
#   5. Final model refit with Firth's penalized logistic regression
#      (logistf, pl = TRUE) on the (up to) 3 retained predictors, to guard
#      against the separation problems that arise with n=55 and many
#      candidate predictors.
#   6. Internal validation: nonparametric bootstrap (boot(), set.seed(123),
#      R = 200 resamples), Optimism_Corrected AUC = 2*apparent_AUC -
#      mean(bootstrap_AUC). Hosmer-Lemeshow goodness-of-fit test (g = 10).
#
# TWO DISCREPANCIES VS. THE MANUSCRIPT TEXT -- DISCLOSED HONESTLY, NOT FIXED
#   (a) BOOTSTRAP ITERATIONS: this script uses R = 200 (see boot(..., R =
#       200, ...) below). The manuscript's Methods states "Internal
#       validation was performed using bootstrap resampling (1,000
#       iterations)." We do not know whether a different, 1000-iteration
#       run was performed later and not saved, or whether the Methods text
#       simply describes the intended/typical protocol rather than this
#       specific run. Both the script and its real saved output are
#       reproduced here exactly as found, without alteration.
#   (b) OPTIMISM-CORRECTED AUC: this script's own saved output (see
#       Etiology_Results_authors_original_output.xlsx, sheet "Bootstrap")
#       gives Original_AUC = 0.955026, Bootstrap_Mean = 0.965301,
#       Optimism_Corrected = 0.944752, 95% CI [0.909754, 1.000000]. The
#       manuscript reports an optimism-corrected AUC of 0.958. The gap
#       (0.945 vs 0.958) is modest but real, and most likely reflects a
#       different bootstrap run (different R, and/or a different random
#       seed / package version) than the one whose output happens to have
#       been saved to disk, rather than a different underlying model --
#       the selected predictors, their coefficients, and the apparent
#       (non-corrected) AUC of 0.955 are all fully consistent with the
#       manuscript.
#
# HONEST WARNING -- COEFFICIENT INSTABILITY
#   The final model's coefficient for
#   "2-Amino-9-[4-hydroxy-2-(hydroxymethyl)butyl]-3H-purin-6-one" is 26.75
#   (OR = 4.16e11, 95% CI 1075 to 1.6e20 -- see Model_Results sheet). A
#   coefficient and odds ratio of this magnitude, with a correspondingly
#   enormous confidence interval, is a classic signature of (quasi-)complete
#   separation in a small sample (n=55 total, ~2:1 ABP:HTGP split): this one
#   metabolite is likely close to perfectly separating the two groups in
#   this dataset. Firth's penalized-likelihood correction (the reason
#   logistf was used here instead of ordinary glm) is specifically designed
#   to reduce, but does not fully eliminate, this kind of small-sample bias
#   and instability. Readers should interpret this predictor's very large
#   apparent effect size and OR with caution, and treat the model's overall
#   AUC as an optimistic estimate of true discriminative performance for a
#   sample this size, even after bootstrap "optimism correction."
#
# WHY THIS WAS NOT RE-EXECUTED IN THE REVIEW ENVIRONMENT
#   The packages this script depends on for the Firth model and goodness-
#   of-fit test (logistf, ResourceSelection, rmda) could not be installed in
#   the sandboxed review environment used to verify supplementary materials
#   (no outbound access to CRAN). Rather than substitute a different
#   estimator, this script's analysis logic is reproduced completely
#   unchanged, and its own real, previously-generated output workbook is
#   included alongside it unmodified, so that the exact code that produced
#   the manuscript's Figure 4 numbers -- not an approximation of it -- is
#   available for inspection.
#
# FILE PATHS
#   The original script pointed data_path/output_dir at the corresponding
#   author's own local folders; those have been replaced below with generic
#   placeholders for public release. Edit them to point at your own copy of
#   Metab_data.xlsx and a writable output folder before running.
# =============================================================================


# ---- Everything below this line is the ORIGINAL script, with only its
# ---- comments translated to English and its local file paths genericized --

# ============================================================================
# Full analysis: ABP vs HTGP prediction model (clinical variables + metabolites)
# Output: ROC curves (combined model + individual variables), boxplots, model results Excel
# ============================================================================

# 1. Load packages -----------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(pROC)
library(ggplot2)
library(ggpubr)
library(rmda)
library(boot)
library(openxlsx)
library(ResourceSelection)
library(glmnet)
library(logistf)

# 2. Set paths and output folder ----------------------------------------------------
data_path <- "PATH/TO/Metab_data.xlsx"        # <- edit to your local copy of the original data file
output_dir <- "PATH/TO/Results/Etiology_ABP_vs_HTGP"  # <- edit to your desired output folder
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("Results will be saved to:", output_dir, "\n")

# 3. Read data ---------------------------------------------------------------
df <- read_excel(data_path, sheet = 1)

# 4. Filter to ABP and HTGP -------------------------------------------------------
df_task <- df %>% filter(Group2 %in% c("ABP", "HTGP"))
df_task <- df_task %>%
  mutate(outcome = ifelse(Group2 == "ABP", 1, 0)) %>%
  filter(!is.na(outcome))

cat("ABP sample size:", sum(df_task$outcome), "\n")
cat("HTGP sample size:", nrow(df_task) - sum(df_task$outcome), "\n")
if(length(unique(df_task$outcome)) < 2) stop("Outcome variable has only one level")

# 5. Define clinical variables (based on actual column names; ensure exact match) ---------------------------
clinical_vars <- c("Age", "Onsettime", "BMI", "Temperature", "Triglyceride",
                   "Balthazar", "RansonCT", "CTSI", "MCTSI",
                   "FibrinogenDegradationProducts", "D_Dimer", "CRP", "PCT",
                   "ALT", "AST", "ALP", "GGT", "Totalbilirubin", "Directbilirubin",
                   "creatinin_e", "GCS")
# Check for missing columns
missing_clinical <- setdiff(clinical_vars, colnames(df_task))
if(length(missing_clinical) > 0) {
  stop("The following clinical variables do not exist: ", paste(missing_clinical, collapse = ", "))
}

# 6. Extract metabolite data (all numeric columns that are not clinical, grouping, or ID columns) -----------------------
exclude_cols <- c("ID", "SampleID", "Group", "Group2", "outcome", clinical_vars)
exclude_cols <- intersect(exclude_cols, colnames(df_task))
metabolite_cols <- setdiff(colnames(df_task), exclude_cols)
is_numeric <- sapply(df_task[, metabolite_cols], is.numeric)
metabolite_cols <- metabolite_cols[is_numeric]
if(length(metabolite_cols) == 0) stop("No numeric metabolite columns detected")
cat("Number of metabolites detected:", length(metabolite_cols), "\n")

metabolite_data <- df_task[, metabolite_cols] %>% mutate(across(everything(), as.numeric))

# 7. Clean clinical variables and drop samples with NA -------------------------------------------
df_task <- df_task %>%
  mutate(across(all_of(clinical_vars), as.numeric)) %>%
  drop_na(all_of(clinical_vars))
metabolite_data <- metabolite_data[rownames(df_task), , drop = FALSE]

# 8. Median-impute missing metabolite values ---------------------------------------------------
fill_na_median <- function(x) {
  if(is.numeric(x)) { med <- median(x, na.rm = TRUE); x[is.na(x)] <- med }
  return(x)
}
metabolite_data <- metabolite_data %>% mutate(across(everything(), fill_na_median))

# 9. Remove constant columns (zero variance) ---------------------------------------------------
zero_var <- sapply(metabolite_data, function(x) length(unique(x)) == 1)
if(sum(zero_var) > 0) {
  cat("Removing constant columns:", colnames(metabolite_data)[zero_var], "\n")
  metabolite_data <- metabolite_data[, !zero_var, drop = FALSE]
}
metabolite_cols <- colnames(metabolite_data)

# 10. Combine features and re-check for missing values ------------------------------------------------
X_raw <- cbind(df_task[, clinical_vars], metabolite_data)
y <- df_task$outcome
X_raw <- X_raw %>% mutate(across(everything(), fill_na_median))
df_model <- cbind(outcome = y, X_raw)
df_model <- df_model[complete.cases(df_model), ]
y <- df_model$outcome
X <- df_model[, -1, drop = FALSE]
cat("Final modeling sample size:", nrow(df_model), "\n")

# 11. LASSO variable selection --------------------------------------------------------
set.seed(123)
X_scaled <- scale(X)
cv_lasso <- cv.glmnet(x = X_scaled, y = y, family = "binomial",
                      alpha = 1, nfolds = min(10, nrow(df_model)))
pdf(file.path(output_dir, "LASSO_plot.pdf"), width = 8, height = 6)
plot(cv_lasso)
dev.off()

coef_lasso <- as.matrix(coef(cv_lasso, s = "lambda.min"))
selected_vars <- rownames(coef_lasso)[which(coef_lasso != 0)][-1]
if(length(selected_vars) == 0) stop("LASSO did not select any variables")
cat("Number of variables selected by LASSO:", length(selected_vars), "\n")
print(selected_vars)

# 12. Keep at most 3 variables (by absolute coefficient) ---------------------------------------
if(length(selected_vars) > 3) {
  coef_abs <- abs(coef_lasso[selected_vars, 1])
  best_vars <- names(sort(coef_abs, decreasing = TRUE))[1:3]
} else {
  best_vars <- selected_vars
}
cat("Final model variables:\n", paste(best_vars, collapse = "\n"), "\n")

# 13. Build the formula (add backticks around special variable names) --------------------------------------
safe_vars <- sapply(best_vars, function(v) {
  if(grepl("^[0-9]|[^a-zA-Z0-9_.]", v)) paste0("`", gsub("`", "", v), "`") else v
})
formula_str <- paste("outcome ~", paste(safe_vars, collapse = " + "))
formula_obj <- as.formula(formula_str)

# 14. Firth logistic regression (handles perfect separation) ------------------------------------------
final_model <- logistf(formula_obj, data = df_model, pl = TRUE)
pred <- 1 / (1 + exp(-final_model$linear.predictors))

# 15. Model evaluation: ROC and AUC --------------------------------------------------
roc_obj <- roc(y, pred)
auc_val <- auc(roc_obj)

# 16. Plot ROC curves (combined model + individual variables) -----------------------------------
clean_vars <- gsub("`", "", best_vars)
roc_curves <- list()
roc_curves[["Combined Model"]] <- roc(y, pred, quiet = TRUE)
for (v in clean_vars) {
  if(v %in% colnames(df_model)) {
    roc_curves[[v]] <- roc(y, df_model[[v]], quiet = TRUE)
  } else {
    warning("Variable does not exist: ", v)
  }
}

n_curves <- length(roc_curves)
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")[1:n_curves]
line_types <- c(1, 2, 3, 4, 5, 6)[1:n_curves]

pdf(file.path(output_dir, "ROC_curves_combined_and_individual.pdf"), width = 8, height = 8)
plot(roc_curves[[1]], col = colors[1], lwd = 2.5, lty = line_types[1],
     main = "ROC Curves: Combined Model vs Individual Variables",
     xlab = "1 - Specificity", ylab = "Sensitivity")
for (i in 2:n_curves) {
  lines(roc_curves[[i]], col = colors[i], lwd = 1.5, lty = line_types[i])
}
legend_auc <- sapply(roc_curves, auc)
legend_labels <- paste0(names(roc_curves), " (AUC = ", round(legend_auc, 3), ")")
legend("bottomright", legend = legend_labels, col = colors, lty = line_types, lwd = 2, cex = 0.8)
dev.off()
cat("ROC curve plot saved\n")

# 17. Boxplots (triglyceride + 4 metabolites) -----------------------------------------
plot_vars <- clean_vars
df_long <- df_model %>%
  select(outcome, all_of(plot_vars)) %>%
  pivot_longer(cols = -outcome, names_to = "Variable", values_to = "Value") %>%
  mutate(Group = ifelse(outcome == 1, "ABP", "HTGP")) %>%
  mutate(Variable = factor(Variable, levels = plot_vars))

# Calculate p-values
p_values <- df_long %>%
  group_by(Variable) %>%
  summarise(p = wilcox.test(Value ~ outcome)$p.value, .groups = "drop") %>%
  mutate(p_label = case_when(
    p < 0.001 ~ "p < 0.001",
    p < 0.01  ~ sprintf("p = %.003f", p),
    p < 0.05  ~ sprintf("p = %.03f", p),
    TRUE      ~ "ns"
  ))

# Sample-size annotation
n_abp <- sum(df_model$outcome == 1)
n_htgp <- sum(df_model$outcome == 0)
sample_label <- paste0("ABP (n=", n_abp, ")   HTGP (n=", n_htgp, ")")

p_box <- ggplot(df_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 1.2, alpha = 0.8, lwd = 0.4) +
  geom_text(data = p_values, aes(x = 1.5, y = Inf, label = p_label),
            inherit.aes = FALSE, vjust = 1.5, hjust = 0.5, size = 3.5) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Expression Level / Value", caption = sample_label) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#F0F0F0", color = "black", size = 0.5),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E0E0E0", size = 0.3),
    plot.caption = element_text(hjust = 0, face = "italic", size = 9)
  ) +
  scale_fill_manual(values = c("ABP" = "#F4A582", "HTGP" = "#92C5DE")) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.8)))

pdf_height <- max(5, length(plot_vars) * 1.3)
pdf(file.path(output_dir, "Boxplots_Triglyceride_and_Metabolites_light.pdf"),
    width = 8, height = pdf_height)
print(p_box)
dev.off()
cat("Boxplot saved\n")

# 18. Calibration curve (manual binning) ---------------------------------------------------
pdf(file.path(output_dir, "Calibration_curve.pdf"), width = 6, height = 6)
breaks <- unique(quantile(pred, probs = seq(0, 1, 0.1), na.rm = TRUE))
if(length(breaks) < 2) breaks <- seq(0, 1, 0.1)
pred_group <- cut(pred, breaks = breaks, include.lowest = TRUE)
obs_rate <- tapply(y, pred_group, mean)
pred_mean <- tapply(pred, pred_group, mean)
plot(pred_mean, obs_rate, xlim = c(0,1), ylim = c(0,1),
     xlab = "Predicted Probability", ylab = "Observed Probability",
     main = "Calibration Curve", pch = 16)
abline(0, 1, lty = 2, col = "gray")
lines(lowess(pred_mean, obs_rate), col = "red", lwd = 2)
dev.off()

# 19. Hosmer-Lemeshow test -------------------------------------------------
if(length(unique(pred)) > 1) {
  hl <- hoslem.test(y, pred, g = 10)
  hl_results <- data.frame(Chi_sq = hl$statistic, df = hl$parameter, P_value = hl$p.value)
} else {
  hl_results <- data.frame(Note = "Hosmer-Lemeshow test failed")
}

# 20. Bootstrap validation of AUC ----------------------------------------------------
set.seed(123)
boot_auc_firth <- function(data, indices, form) {
  d <- data[indices, ]
  model <- tryCatch(logistf(form, data = d, pl = FALSE), error = function(e) NULL)
  if(is.null(model)) return(NA)
  pred_boot <- 1 / (1 + exp(-model$linear.predictors))
  roc_obj <- tryCatch(roc(d$outcome, pred_boot, quiet = TRUE), error = function(e) NULL)
  if(is.null(roc_obj)) return(NA)
  return(auc(roc_obj))
}
boot_res <- boot(df_model, statistic = boot_auc_firth, R = 200, form = formula_obj)
boot_vals <- boot_res$t[!is.na(boot_res$t)]
if(length(boot_vals) > 0) {
  boot_summary <- data.frame(
    Original_AUC = auc_val,
    Bootstrap_Mean = mean(boot_vals),
    Optimism_Corrected = 2*auc_val - mean(boot_vals),
    CI_lower = quantile(boot_vals, 0.025),
    CI_upper = quantile(boot_vals, 0.975)
  )
} else {
  boot_summary <- data.frame(Note = "Bootstrap failed")
}

# 21. Extract model coefficients and save to Excel ---------------------------------------------
coef_vec <- coef(final_model)
se_vec <- sqrt(diag(vcov(final_model)))
z_vec <- coef_vec / se_vec
p_vec <- 2 * (1 - pnorm(abs(z_vec)))
or_vec <- exp(coef_vec)
or_lower <- exp(coef_vec - 1.96 * se_vec)
or_upper <- exp(coef_vec + 1.96 * se_vec)

results <- data.frame(
  Variable = names(coef_vec),
  Coefficient = coef_vec,
  SE = se_vec,
  Z = z_vec,
  P_value = p_vec,
  OR = or_vec,
  OR_CI_low = or_lower,
  OR_CI_high = or_upper
)

wb <- createWorkbook()
addWorksheet(wb, "Model_Results")
writeData(wb, "Model_Results", results)
addWorksheet(wb, "Bootstrap")
writeData(wb, "Bootstrap", boot_summary)
addWorksheet(wb, "Hosmer_Lemeshow")
writeData(wb, "Hosmer_Lemeshow", hl_results)
saveWorkbook(wb, file.path(output_dir, "Etiology_Results.xlsx"), overwrite = TRUE)

cat("Analysis complete! Results saved in:", output_dir, "\n")
