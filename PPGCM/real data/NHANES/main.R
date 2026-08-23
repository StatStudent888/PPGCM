# ============================================================
# NHANES real-data analysis for PPGCM
#
# Purpose:
#   This script implements the semi-supervised NHANES analysis.
#   It tests the conditional independence between waist circumference
#   X and two-hour OGTT plasma glucose Y given cardiometabolic
#   covariates Z, using HbA1c as the auxiliary variable S.
#
# Data structure:
#   Labeled observations contain (Y, X, S, Z), whereas unlabeled
#   observations contain only (X, S, Z). The labeled and unlabeled
#   samples arise from the actual NHANES data-collection design.
#
# Analysis steps:
#   1. load and standardize the labeled and unlabeled samples;
#   2. test X ⟂ S | Z using all observations for which X, S, and Z
#      are available;
#   3. assess the predictive contribution of S to Y;
#   4. apply GCM, PPGCM (w=1), and PPGCM (opt);
#   5. apply CRT and Maxway CRT for comparison.
#
# Variable definitions:
#   Y = two-hour OGTT plasma glucose
#   X = waist circumference
#   S = HbA1c
#   Z = BMI, age, systolic blood pressure, HDL cholesterol,
#       total cholesterol, and creatinine
#
# Main outputs:
#   - GCM result for testing X ⟂ S | Z;
#   - PPGCM and GCM p-values, test statistics, variance estimates,
#     and estimated weights;
#   - CRT and Maxway CRT p-values and test statistics.
#
# Required files:
#   - PPGCM.R
#   - maxway.R
#
# Required packages:
#   dplyr, isotree, randomForest
# ============================================================

# Clear the current R session, load the testing implementations,
# and initialize the required packages and random seed.
gc(); rm(list = ls()); cat("\014")
library(dplyr)
source('PPGCM/PPGCM.R')
library(isotree)
source("PPGCM/simulation/type-I_power/CRT_MaxwayCRT/maxway.R")
library(randomForest)
options(warn = -1)
set.seed(123)


# ============================================================
# Load the labeled and unlabeled NHANES samples
# ============================================================

# Define the outcome, exposure, auxiliary variable, and conditioning covariates.
y_col <- "Y_ogtt_2h_mgdl"
x_col <- "X_waist_cm"
s_col <- "S_hba1c_pct"
z_cols <- c("Z_bmi", "Z_age_years", "Z_sbp_mmhg", "Z_hdl_mgdl",
            "Z_total_chol_mgdl", "Z_creatinine_mgdl")

# Load the genuinely labeled and unlabeled samples.
label_raw <- read.csv("PPGCM/real data/NHANES/label.csv", check.names = FALSE)
unlabel_raw <- read.csv("PPGCM/real data/NHANES/unlabel.csv", check.names = FALSE)

# ------------------------------------------------------------
# Standardize the analysis variables.
#
# X, S, and Z are standardized using their pooled labeled and unlabeled
# distributions, since these variables are observed in both samples.
# Y is standardized using the labeled sample only because it is
# unavailable in the unlabeled data.
# ------------------------------------------------------------
standardize_inputs <- function(label, unlabel) {
  shared <- c(x_col, s_col, z_cols)
  # Pool labeled and unlabeled observations to estimate the means and
  # standard deviations of variables observed in both samples.
  combined <- rbind(label[, shared, drop = FALSE], unlabel[, shared, drop = FALSE])
  centers <- vapply(combined, mean, numeric(1))
  scales <- vapply(combined, sd, numeric(1))
  if (any(!is.finite(scales) | scales <= 0)) stop("A shared analysis variable has zero variance.")
  label_std <- label[, c(y_col, shared), drop = FALSE]
  unlabel_std <- unlabel[, shared, drop = FALSE]
  for (v in shared) {
    label_std[[v]] <- (label_std[[v]] - centers[[v]]) / scales[[v]]
    unlabel_std[[v]] <- (unlabel_std[[v]] - centers[[v]]) / scales[[v]]
  }
  # Standardize Y using the labeled observations only.
  y_center <- mean(label_std[[y_col]])
  y_scale <- sd(label_std[[y_col]])
  label_std[[y_col]] <- (label_std[[y_col]] - y_center) / y_scale
  list(label = label_std, unlabel = unlabel_std)
}

# Apply the same standardization to the labeled and unlabeled samples.
dat <- standardize_inputs(label_raw, unlabel_raw)
label_raw <- dat$label
unlabel_raw <- dat$unlabel

# ============================================================
# Test X ⟂ S | Z using all available observations
# ============================================================
# Pool the labeled and unlabeled samples because this diagnostic involves
# only X, S, and Z and therefore does not require Y.

# Use generalized random forests for nuisance-function and
# prediction-model estimation throughout the NHANES analysis.
method_fit = "grf"

# Use generalized random forests for nuisance-function and
# prediction-model estimation throughout the NHANES analysis.
dataXSZ <- rbind((label_raw %>% select(-y_col)),unlabel_raw)

# Treat S as the response and apply the cross-fitted GCM to test X ⟂ S | Z.
GCM <- ppicit_test(
  label_data = dataXSZ,
  unlabel_data = NULL,
  K = 20,
  fit_h = method_fit,
  fit_g = method_fit,
  fit_f = method_fit,
  x_col = x_col,
  y_col = s_col,
  s_col = NULL,
  z_cols = NULL,
  alpha = 0.05,
  seed = 123
)
cat("CIT of X and S given Z:","\n",
    "GCM p-value: ",GCM$gcm$p_value,"\n",
    "GCM reject: ",GCM$gcm$reject,"\n")

# ============================================================
# Assess the predictive contribution of the auxiliary variable S to Y
# ============================================================
# First report the marginal correlation between HbA1c and OGTT glucose.
# Then regress Y on X, Z, and S to assess whether S remains predictive
# of Y after adjusting for the exposure and conditioning covariates.
print(cor.test(label_raw[[s_col]], label_raw[[y_col]]))

# Fit the full linear regression Y ~ X + Z + S.
selected_vars <- c(x_col,z_cols,s_col)
df_lm <- data.frame(label_raw[, selected_vars])
Y <- as.matrix(label_raw[y_col])
lm_fit <- lm(Y ~ ., data = df_lm)
summary(lm_fit)  

# ============================================================
# Main semi-supervised analysis: GCM and PPGCM
# ============================================================
# Apply the cross-fitted GCM and the two PPGCM variants to the genuine
# labeled and unlabeled NHANES samples.

# PPGCM uses HbA1c through the prediction model f(X,Z,S), while the
# tested hypothesis remains X ⟂ Y | Z.
res_r <- ppicit_test(
  label_data = label_raw,
  unlabel_data = unlabel_raw,
  K = 20,
  fit_h = method_fit,
  fit_g = method_fit,
  fit_f = method_fit,
  x_col = x_col,
  y_col = y_col,
  s_col = s_col,
  z_cols = NULL,
  alpha = 0.05,
  seed = 123
)

# Extract the main inferential quantities returned by each procedure.
if (is.null(res_r)) {
  has_error  <- TRUE
} else {
  ppi_pvals  <- res_r$ppi$p_value
  ppi_reject  <- res_r$ppi$reject
  ppi_T  <- res_r$ppi$T_hat
  ppi_var  <- res_r$ppi$variance_hat
  ppi_weight  <- res_r$ppi$weight_hat_mean
  
  ppi_w1_pvals  <- res_r$ppi_w1$p_value
  ppi_w1_reject  <- res_r$ppi_w1$reject
  ppi_w1_T  <- res_r$ppi_w1$T_hat
  ppi_w1_var  <- res_r$ppi_w1$variance_hat
  ppi_w1_weight  <- res_r$ppi_w1$weight_hat_mean
  
  gcm_pvals  <- res_r$gcm$p_value
  gcm_reject  <- res_r$gcm$reject
  gcm_T  <- res_r$gcm$T_hat
  gcm_var  <- res_r$gcm$variance_hat
}

# Assemble a compact comparison of the two PPGCM variants and GCM.
summary <- data.frame(
  method = c("PPI-GCM", "PPI-GCM (w=1)", "GCM"),
  p_value = c(
    mean(ppi_pvals ),
    mean(ppi_w1_pvals ),
    mean(gcm_pvals )
  ),
  T_hat = c(
    mean(ppi_T ),
    mean(ppi_w1_T ),
    mean(gcm_T )
  ),
  variance_hat = c(
    mean(ppi_var ),
    mean(ppi_w1_var ),
    mean(gcm_var )
  ),
  weight_hat = c(
    mean(ppi_weight ),
    mean(ppi_w1_weight ),
    0
  ))
print(summary)

# Assemble a compact comparison of the two PPGCM variants and GCM.

# Remove S from the CRT analysis because the competing CRT procedures
# condition only on Z and do not use the auxiliary variable.
label_raw <- label_raw %>% select(-s_col)

# Use the nonlinear/random-forest implementation for the NHANES data.
res_r <-  maxwayCRT_insample_one(
  label_data = label_raw,
  unlabel_data = unlabel_raw,
  gbar_type = "nonlinear",
  x_col = x_col,
  y_col = y_col,
  z_cols = NULL,
  seed = 123, 
  M = 200
)
if (is.null(res_r)) {
  has_error  <- TRUE
} else {
  maxway_pvals  <- res_r$MaxwayCRT$p_value
  maxway_T  <- res_r$MaxwayCRT$T_signed
  
  crt_pvals  <- res_r$CRT$p_value
  crt_T  <- res_r$CRT$T_signed
}

# Assemble a compact comparison of CRT and Maxway CRT.
summaryCRT <- data.frame(
  method = c("Maxway CRT", "CRT"),
  mean_p_value = c(
    mean(maxway_pvals ),
    mean(crt_pvals )
  ),
  mean_T_signed = c(
    mean(maxway_T ),
    mean(crt_T )
  ))


print(summaryCRT)


