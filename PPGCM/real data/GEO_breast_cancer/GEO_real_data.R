# ============================================================
# GEO breast cancer real-data analysis for PPGCM
#
# Purpose:
#   This script implements the breast cancer analysis in the paper.
#   It preprocesses the GEO gene-expression data, removes anomalous
#   observations, and tests the conditional independence between
#   ESR1 (X) and CCND1 (Y) given the remaining PAM50 genes (Z),
#   using PPFIA1 as the auxiliary variable S.
#
# Analysis steps:
#   1. restrict the data to Luminal A/B breast cancer samples;
#   2. retain ESR1, CCND1, PPFIA1, and the PAM50 genes;
#   3. remove anomalous observations using isolation forest;
#   4. standardize all retained gene-expression variables;
#   5. assess whether PPFIA1 is predictive of CCND1;
#   6. apply the full-data GCM to X and Y given Z;
#   7. test X ⟂ S | Z using the full sample;
#   8. repeatedly mask Y to construct an artificial semi-supervised
#      setting and compare GCM, PPGCM (w=1), PPGCM (opt), CRT,
#      and Maxway CRT.
#
# Variable definitions:
#   Y = expression of CCND1
#   X = expression of ESR1
#   S = expression of PPFIA1
#   Z = the remaining PAM50 genes after excluding ESR1
#
# Main outputs:
#   - full-data GCM result;
#   - GCM result for testing X ⟂ S | Z;
#   - Monte Carlo summaries for GCM and the two PPGCM variants;
#   - Monte Carlo summaries for CRT and Maxway CRT.
#
# Required files:
#   - PPGCM.R
#   - maxway.R
#
# Required packages:
#   glmnet, limma, dplyr, isotree
# ============================================================

# Clear the current R session, load the testing implementations,
# and initialize the packages and random seed used in the analysis.
gc(); rm(list = ls()); cat("\014")
source('PPGCM/PPGCM.R')
source("PPGCM/simulation/type-I_power/CRT_MaxwayCRT/maxway.R")
library(glmnet)
library(limma)
library(dplyr)
library(isotree)
options(warn = -1)
seed = 123
set.seed(seed)

# Compute a standard deviation safely when at least two valid values are available.
sd_or_na <- function(x) {
  if (sum(!is.na(x)) <= 1) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

# ============================================================
# Load and preprocess the GEO breast cancer data
# ============================================================
# Load the combined gene-expression/clinical dataset and remove
# observations with missing values.
GEO_gene_clinical0 <- read.csv("PPGCM/real data/GEO_breast_cancer/GEO_data.csv", stringsAsFactors = FALSE)
GEO_gene_clinical0 <- na.omit(GEO_gene_clinical0)

GEO_gene_clinical <- GEO_gene_clinical0[,-c(1,3,4,5)]
# Restrict the analysis to the luminal breast cancer subtypes
# by excluding Basal and HER2-enriched samples.
GEO_gene_clinical <- subset(GEO_gene_clinical, subtype != "Basal")
GEO_gene_clinical <- subset(GEO_gene_clinical, subtype != "Her2")
GEO_gene_clinical <- GEO_gene_clinical[,-c(1)]

# List the PAM50 genes used to define the conditioning set Z.
genes <- c(
  "ACTR3B", "ANLN", "BAG1", "BCL2", "BIRC5", "BLVRA", "CCNB1", "CCNE1",
  "CDC20", "CDC6", "CDH3", "CENPF", "CEP55", "CXXC5", "EGFR", "ERBB2",
  "ESR1", "EXO1", "FGFR4", "FOXA1", "FOXC1", "GPR160", "GRB7", "KIF2C",
  "KRT14", "KRT17", "KRT5", "MAPT", "MDM2", "MELK", "MIA", "MKI67",
  "MLPH", "MMP11", "MYBL2", "MYC", "NAT1", "NDC80", "NUF2", "ORC6",
  "PGR", "PHGDH", "PTTG1", "RRM2", "SFRP1", "SLC39A6", "TMEM45B",
  "TYMS", "UBE2C", "UBE2T"
)

# Define the outcome, exposure, and auxiliary gene used in the analysis.
Y_var <- "CCND1"          
X_var <- "ESR1" 
S_var <- "PPFIA1"
# Retain the three focal genes together with the PAM50 adjustment genes.
genes <- c(Y_var,X_var,S_var,genes)

GEO_gene_clinical <- GEO_gene_clinical[, colnames(GEO_gene_clinical) %in% genes]

# ============================================================
# Remove anomalous observations using isolation forest
# ============================================================
# Fit an isolation forest to the retained gene-expression variables
# and remove observations whose anomaly scores exceed the 95th percentile.
model <- isolation.forest(
  GEO_gene_clinical,
  ntrees = 100,
  sample_size = 256,
  ndim = ncol(GEO_gene_clinical),
  nthreads = 1)
scores <- predict(model, GEO_gene_clinical, type = "score")
GEO_gene_clinical$anomaly_score <- scores
threshold <- quantile(scores, 0.95)
GEO_gene_clinical$is_anomaly <- GEO_gene_clinical$anomaly_score > threshold

# Identify and remove the observations flagged as anomalous.
anomalies <- GEO_gene_clinical %>% filter(is_anomaly) %>% arrange(desc(anomaly_score))
GEO_gene_clinical <- GEO_gene_clinical[!GEO_gene_clinical$is_anomaly, ]
GEO_gene_clinical <- GEO_gene_clinical %>% select(-"is_anomaly")
GEO_gene_clinical <- GEO_gene_clinical %>% select(-"anomaly_score")
GEO_gene_clinical <- as.matrix(GEO_gene_clinical)
GEO_gene_clinical[is.infinite(GEO_gene_clinical)] <- NA
GEO_gene_clinical <- na.omit(as.data.frame(GEO_gene_clinical))
GEO_gene_clinical <- na.omit(GEO_gene_clinical)

# Standardize all retained gene-expression variables before fitting the tests.
considered_data <- as.data.frame(scale(GEO_gene_clinical))

# ============================================================
# Assess the predictive contribution of the auxiliary variable S to Y
# ============================================================
# Use lasso to select predictors of Y from X, S, and Z, and then refit
# an ordinary linear regression on the selected variables to assess the
# contribution and statistical significance of the selected predictors.

dataSZ <- considered_data %>% select(-Y_var)
dataSZ <- as.matrix(dataSZ)
Y <- as.matrix(considered_data[Y_var])

# Select predictors using 10-fold cross-validated lasso.
cv_fit <- cv.glmnet(
  x = dataSZ,
  y = Y,
  alpha = 1,
  nfolds = 10
)
lambda_opt <- cv_fit$lambda.min

final_model <- glmnet(
  dataSZ,
  Y,
  alpha = 1,
  lambda = lambda_opt
)

# Identify variables with nonzero lasso coefficients.
coef_vec <- coef(final_model)
selected_vars <- which(coef_vec[-1] != 0)

# Refit an ordinary linear regression on the selected variables
# to obtain coefficient estimates and p-values.
if (length(selected_vars) == 0) {
  cat("No variables were selected by the lasso.\n")
} else {
  df_lm <- data.frame(dataSZ[, selected_vars, drop = FALSE])
  lm_fit <- lm(Y ~ ., data = df_lm)
  print(summary(lm_fit))
}

# ============================================================
# Regression method
# ============================================================
# Use lasso for nuisance-function and prediction-model estimation
# throughout the GEO analysis.
usedmethod <- "lasso"    

# ============================================================
# Full-data GCM test of X ⟂ Y | Z
# ============================================================
# Use all observations with Y available to assess the conditional
# association between ESR1 and CCND1 given the PAM50 covariates.

# Exclude S because the null hypothesis conditions only on Z.
dataXYZ <- considered_data %>% select(-S_var)
GCM_Oracle <- ppicit_test(
  label_data = dataXYZ,
  unlabel_data = NULL,
  K = 20,
  fit_h = usedmethod,
  fit_g = usedmethod,
  fit_f = usedmethod,
  x_col = X_var,
  y_col = Y_var,
  s_col = NULL,
  z_cols = NULL,
  alpha = 0.05,
  seed = seed
)
cat("Oracle GCM of X and Y given Z:","\n",
    "GCM p-value: ",GCM_Oracle$gcm$p_value,"\n",
    "GCM reject: ",GCM_Oracle$gcm$reject,"\n",
    "GCM variance: ",GCM_Oracle$gcm$variance_hat,"\n",
    "GCM statistic: ",GCM_Oracle$gcm$T_hat,"\n")

# ============================================================
# Full-data GCM test of X ⟂ S | Z
# ============================================================
# Check whether ESR1 and PPFIA1 remain conditionally associated
# after adjusting for Z. This diagnostic uses no Y information.

# Remove Y and treat S as the response in the conditional-independence test.
dataXSZ <- considered_data %>% select(-Y_var)
GCM <- ppicit_test(
  label_data = dataXSZ,
  unlabel_data = NULL,
  K = 20,
  fit_h = usedmethod,
  fit_g = usedmethod,
  fit_f = usedmethod,
  x_col = X_var,
  y_col = S_var,
  s_col = NULL,
  z_cols = NULL,
  alpha = 0.05,
  seed = seed
)
cat("CIT of X and S given Z:","\n",
    "GCM p-value: ",GCM$gcm$p_value,"\n",
    "GCM reject: ",GCM$gcm$reject,"\n")

# ============================================================
# Artificial semi-supervised experiment: GCM and PPGCM
# ============================================================
# Repeatedly draw 300 observations as the labeled sample and treat
# all remaining observations as unlabeled by masking Y. The same split
# is used to compare GCM, PPGCM (w=1), and PPGCM (opt).

# Monte Carlo configuration matching the real-data experiment in the paper.
n_rep <- 2000
sample_size <-  300

# Initialize storage for p-values, rejection indicators, test statistics,
# variance estimates, and PPGCM weights across random splits.
ppi_pvals <- rep(NA_real_, n_rep)
ppi_reject <- rep(NA, n_rep)
ppi_T <- rep(NA_real_, n_rep)
ppi_var <- rep(NA_real_, n_rep)
ppi_weight <- rep(NA_real_, n_rep)

gcm_pvals <- rep(NA_real_, n_rep)
gcm_reject <- rep(NA, n_rep)
gcm_T <- rep(NA_real_, n_rep)
gcm_var <- rep(NA_real_, n_rep)

ppi_w1_pvals <- rep(NA_real_, n_rep)
ppi_w1_reject <- rep(NA, n_rep)
ppi_w1_T <- rep(NA_real_, n_rep)
ppi_w1_var <- rep(NA_real_, n_rep)
ppi_w1_weight <- rep(NA_real_, n_rep)
has_error <- rep(FALSE, n_rep)

n_total <- nrow(considered_data)

# Repeat the random labeled/unlabeled split and testing procedure.
for (r in seq_len(n_rep)){
  set.seed(seed+r)
  # Randomly select the labeled sample and use the remaining observations
  # as the unlabeled sample.
  sampled_indices <- sample(1:n_total, size = sample_size, replace = FALSE)
  label_realdata <- considered_data[sampled_indices, ]
  unlabel_realdata_full <- considered_data[-sampled_indices, ]
  # Mask Y in the unlabeled sample while retaining X, S, and Z.
  unlabel_realdata <- unlabel_realdata_full%>% select(-Y_var)
  
  # Apply the cross-fitted GCM and the two PPGCM variants to this split.
  res_r <- ppicit_test(
    label_data = label_realdata,
    unlabel_data = unlabel_realdata,
    K = 20,
    fit_h = usedmethod,
    fit_g = usedmethod,
    fit_f = usedmethod,
    x_col = X_var,
    y_col = Y_var,
    s_col = S_var,
    z_cols = NULL,
    alpha = 0.05,
    seed = seed + r
  )
  # Store the inferential quantities for successful replications and
  # flag failed runs without terminating the Monte Carlo experiment.
  if (is.null(res_r)) {
    has_error[r] <- TRUE
  } else {
    ppi_pvals[r] <- res_r$ppi$p_value
    ppi_reject[r] <- res_r$ppi$reject
    ppi_T[r] <- res_r$ppi$T_hat
    ppi_var[r] <- res_r$ppi$variance_hat
    ppi_weight[r] <- res_r$ppi$weight_hat_mean
    
    ppi_w1_pvals[r] <- res_r$ppi_w1$p_value
    ppi_w1_reject[r] <- res_r$ppi_w1$reject
    ppi_w1_T[r] <- res_r$ppi_w1$T_hat
    ppi_w1_var[r] <- res_r$ppi_w1$variance_hat
    ppi_w1_weight[r] <- res_r$ppi_w1$weight_hat_mean
    
    gcm_pvals[r] <- res_r$gcm$p_value
    gcm_reject[r] <- res_r$gcm$reject
    gcm_T[r] <- res_r$gcm$T_hat
    gcm_var[r] <- res_r$gcm$variance_hat  
  }
}

# Restrict the Monte Carlo summary to splits completed successfully
# for all three procedures.
valid <- !has_error &
  !is.na(ppi_reject) &
  !is.na(ppi_w1_reject) &
  !is.na(gcm_reject)

# Summarize rejection rates, p-values, test statistics, estimated
# asymptotic variances, and PPGCM weights across random splits.
summary <- data.frame(
  method = c("PPI-GCM", "PPI-GCM (w=1)", "GCM"),
  rejection_rate = c(
    mean(ppi_reject[valid]),
    mean(ppi_w1_reject[valid]),
    mean(gcm_reject[valid])
  ),
  mean_p_value = c(
    mean(ppi_pvals[valid]),
    mean(ppi_w1_pvals[valid]),
    mean(gcm_pvals[valid])
  ),
  sd_p_value = c(
    sd_or_na(ppi_pvals[valid]),
    sd_or_na(ppi_w1_pvals[valid]),
    sd_or_na(gcm_pvals[valid])
  ),
  median_p_value = c(
    stats::median(ppi_pvals[valid]),
    stats::median(ppi_w1_pvals[valid]),
    stats::median(gcm_pvals[valid])
  ),
  mean_T_hat = c(
    mean(ppi_T[valid]),
    mean(ppi_w1_T[valid]),
    mean(gcm_T[valid])
  ),
  sd_T_hat = c(
    sd_or_na(ppi_T[valid]),
    sd_or_na(ppi_w1_T[valid]),
    sd_or_na(gcm_T[valid])
  ),
  mean_variance_hat = c(
    mean(ppi_var[valid]),
    mean(ppi_w1_var[valid]),
    mean(gcm_var[valid])
  ),
  sd_variance_hat = c(
    sd_or_na(ppi_var[valid]),
    sd_or_na(ppi_w1_var[valid]),
    sd_or_na(gcm_var[valid])
  ),
  mean_weight_hat = c(
    mean(ppi_weight[valid]),
    mean(ppi_w1_weight[valid]),
    0
  ),
  sd_weight_hat = c(
    sd_or_na(ppi_weight[valid]),
    sd_or_na(ppi_w1_weight[valid]),
    0
  ))

print(summary)


# ============================================================
# Artificial semi-supervised experiment: CRT and Maxway CRT
# ============================================================
# Use the same random labeled/unlabeled splits as above. Neither CRT
# nor Maxway CRT uses the auxiliary variable S in this analysis.
alpha = 0.05

# Initialize storage for rejection indicators, p-values, and signed
# residual-product statistics across random splits.
maxway_reject <- rep(NA_real_, n_rep)
crt_reject <- rep(NA_real_, n_rep)
maxway_pvals <- rep(NA_real_, n_rep)
crt_pvals <- rep(NA_real_, n_rep)
maxway_T <- rep(NA_real_, n_rep)
crt_T <- rep(NA_real_, n_rep)
has_error <- rep(FALSE, n_rep)

# Repeat the same labeled/unlabeled splitting scheme used for PPGCM.
for (r in seq_len(n_rep)){
  set.seed(seed+r)
  sampled_indices <- sample(1:n_total, size = sample_size, replace = FALSE)
  label_realdata <- considered_data[sampled_indices, ]
  
  # Remove S because the CRT-type competitors do not use the auxiliary
  # variable, and mask Y in the unlabeled sample.
  label_realdata <- label_realdata %>% select(-S_var)
  unlabel_realdata_full <- considered_data[-sampled_indices, ]
  unlabel_realdata <- unlabel_realdata_full %>% select(-Y_var)
  unlabel_realdata <- unlabel_realdata %>% select(-S_var)
  
  # Apply CRT and the in-sample Maxway CRT using the linear/lasso implementation.
  res_r <-  maxwayCRT_insample_one(
    label_data = label_realdata,
    unlabel_data = unlabel_realdata,
    gbar_type = "linear",
    x_col = X_var,
    y_col = Y_var,
    z_cols = NULL,
    seed = seed + r, 
    M = 200
  )
  if (is.null(res_r)) {
    has_error[r] <- TRUE
  } else {
    maxway_pvals[r] <- res_r$MaxwayCRT$p_value
    maxway_reject[r] <- res_r$MaxwayCRT$p_value <= alpha
    maxway_T[r] <- res_r$MaxwayCRT$T_signed
    
    crt_pvals[r] <- res_r$CRT$p_value
    crt_reject[r] <- res_r$CRT$p_value <= alpha
    crt_T[r] <- res_r$CRT$T_signed
  }
}

# Keep only splits completed successfully for both CRT procedures.
valid <- !has_error &
  !is.na(maxway_reject) &
  !is.na(crt_reject) 
summaryCRT <- data.frame(
  method = c("Maxway CRT", "CRT"),
  rejection_rate = c(
    mean(maxway_reject[valid]),
    mean(crt_reject[valid])
  ),
  mean_p_value = c(
    mean(maxway_pvals[valid]),
    mean(crt_pvals[valid])
  ),
  sd_p_value = c(
    sd_or_na(maxway_pvals[valid]),
    sd_or_na(crt_pvals[valid])
  ),
  median_p_value = c(
    stats::median(maxway_pvals[valid]),
    stats::median(crt_pvals[valid])
  ),
  mean_T_signed = c(
    mean(maxway_T[valid]),
    mean(crt_T[valid])
  ),
  sd_T_signed = c(
    sd_or_na(maxway_T[valid]),
    sd_or_na(crt_T[valid])
  ))

# Summarize rejection rates, p-values, and signed test statistics
# for CRT and Maxway CRT across the random splits.
print(summaryCRT)




