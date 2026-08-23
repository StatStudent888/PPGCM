# ============================================================
# Monte Carlo simulation driver for CRT and Maxway CRT
#
# Purpose:
#   This script runs the CRT and the in-sample/CV version of Maxway CRT
#   under the linear and nonlinear simulation designs used in the PPGCM
#   study. For each signal level c, it repeatedly generates labeled and
#   unlabeled samples, computes CRT and Maxway CRT p-values, and summarizes
#   their empirical rejection rates and randomization-statistic behavior.
#
# Data:
#   Labeled observations contain (Y, X, S, Z), whereas unlabeled
#   observations contain only (X, S, Z). CRT and Maxway CRT use X and Z
#   from the unlabeled sample to estimate X | Z; the auxiliary variable S
#   is not used by either method.
#
# Implementation:
#   - linear design: lasso is used for X | Z and Y | Z;
#   - nonlinear design: random forests are used for nuisance estimation
#     and Maxway calibration;
#   - in the nonlinear design, the original RF variable-importance
#     statistic is replaced by the residual-product statistic
#         |mean(r_X r_Y)|,
#     so that the test statistic is comparable with the GCM-based
#     statistic used in the PPGCM simulations;
#   - Monte Carlo replications can be run in parallel, while each
#     individual replication is forced to use a single computational thread.
#
# Main inputs:
#   n_rep      : number of Monte Carlo replications
#   c_coef     : signal coefficient c
#   n_label    : labeled sample size n
#   n_unlabel  : unlabeled sample size N
#   gbar_type  : "linear" or "nonlinear" data-generating mechanism
#   M          : number of CRT randomization draws
#
# Main outputs:
#   summary : method-level Monte Carlo summaries;
#   details : replication-level p-values and test statistics;
#   settings: complete simulation configuration.
#
# Required files:
#   - generate_data.R
#   - maxway.R
#
# Required packages:
#   glmnet, randomForest, MASS
#   CompQuadForm is additionally required when use_interaction = TRUE.
# ============================================================

# ------------------------------------------------------------
# Force all numerical libraries used within a replication to one thread.
# Parallelism is handled across Monte Carlo replications rather than
# within individual regression or matrix computations.
# ------------------------------------------------------------
.force_single_threading <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
  }
  invisible(NULL)
}

# Load a required package and stop with an informative message if unavailable.
.load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Please install it first.", pkg), call. = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Keep the master process single-threaded as well. The parallelism is across replications.
.force_single_threading()

.load_pkg("glmnet")
.load_pkg("randomForest")
.load_pkg("MASS")

# CompQuadForm is only needed if use_interaction = TRUE for the linear dICRT statistic.
if (requireNamespace("CompQuadForm", quietly = TRUE)) {
  suppressPackageStartupMessages(library(CompQuadForm))
}

# Set the project directory and load the data-generating mechanism
# together with the original Maxway helper functions.

setwd("PPGCM")
source("simulation/type-I_power/generate_data.R")
source("simulation/type-I_power/CRT_MaxwayCRT/maxway.R")

# ----------------------------
# 1. Small helper functions
# Utilities for generating fixed sparse coefficients, safely summarizing
# Monte Carlo output, and standardizing matrix inputs.
# ----------------------------

# Construct a sparse coefficient vector with Gaussian nonzero entries
# on a user-specified consecutive support.
make_sparse_beta <- function(p_z, beta_start, beta_nnz, z_names = paste0("Z", seq_len(p_z))) {
  if (beta_start < 1) stop("beta_start must be at least 1.")
  if (beta_start + beta_nnz - 1 > p_z) {
    stop("beta_start + beta_nnz - 1 cannot exceed p_z.")
  }
  beta <- rep(0, p_z)
  idx <- beta_start:(beta_start + beta_nnz - 1)
  beta[idx] <- stats::rnorm(beta_nnz)
  names(beta) <- z_names
  beta
}

# Summary functions that return NA rather than failing when no valid
# Monte Carlo observations are available.
safe_mean <- function(x) {
  if (sum(!is.na(x)) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (sum(!is.na(x)) == 0) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) <= 1) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

# Convert an input object to a numeric matrix for regression functions.
as_num_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

# ------------------------------------------------------------
# 2. Linear dCRT statistic with the observed test statistic returned.
#
# This is the same residual-based statistic used by dCRT_res() in
# maxway.R.
# ------------------------------------------------------------
dCRT_res_with_stat <- function(res_x, res_y, Z_sub = NULL, k = 0, d.interaction = FALSE) {
  res_x <- as.numeric(res_x)
  res_y <- as.numeric(res_y)
  
  # Standardize the fitted X residuals before forming the test statistic.
  sx <- stats::sd(res_x)
  if (!is.finite(sx) || sx <= 0) stop("sd(res_x) is zero or non-finite.")
  res_x <- res_x / sx

  n <- length(res_x)
  # Basic residual-product dCRT statistic.
  if (!d.interaction || k == 0) {
    T_signed <- mean(res_x * res_y)
    T_hat <- abs(T_signed)
    emp_var <- mean(res_y^2)
    p_value <- 2 * stats::pnorm(-sqrt(n) * T_hat / sqrt(emp_var))
    return(list(
      p_value = as.numeric(p_value),
      T_hat = as.numeric(T_hat),
      T_signed = as.numeric(T_signed),
      T_null_mean = NA_real_,
      T_null_sd = NA_real_,
      T_centered = NA_real_,
      T_z = NA_real_
    ))
  }
  
  # Interaction-augmented dCRT statistic using selected Z variables.
  if (!requireNamespace("CompQuadForm", quietly = TRUE)) {
    stop("Package 'CompQuadForm' is required when use_interaction = TRUE.")
  }

  Z_sub <- as_num_matrix(Z_sub)
  weight_inter <- 1 / sqrt(k)
  W_inter <- res_y
  for (l in seq_len(k)) {
    W_inter <- cbind(W_inter, Z_sub[, l] * res_y)
  }

  X_design <- cbind(1, Z_sub)
  XTX_inv <- solve(t(X_design) %*% X_design)

  Z_dI <- diag(c(1, rep(weight_inter, k))) %*%
    XTX_inv %*% t(W_inter) %*% res_x
  T_hat <- sum(Z_dI^2)

  WTW <- t(W_inter) %*% W_inter
  svd_WTW <- svd(WTW)
  root_WTW <- svd_WTW$u %*% diag(sqrt(svd_WTW$d)) %*% t(svd_WTW$v)
  lambda_W <- svd(root_WTW %*% XTX_inv %*%
                    diag(c(1, rep(weight_inter^2, k))) %*%
                    XTX_inv %*% root_WTW)$d

  p_value <- tryCatch({
    out <- CompQuadForm::imhof(T_hat, lambda_W, epsabs = 1e-8)
    as.vector(abs(out$Qq))
  }, error = function(e) 0)

  list(
    p_value = as.numeric(p_value),
    T_hat = as.numeric(T_hat),
    T_signed = NA_real_,
    T_null_mean = NA_real_,
    T_null_sd = NA_real_,
    T_centered = NA_real_,
    T_z = NA_real_
  )
}

# ------------------------------------------------------------
# 3. Residual-product CRT for the nonlinear design.
#
# Random forests are used only to estimate X | Z and Y | Z. The test
# statistic itself is the residual product
#
#   T = |mean(r_X r_Y)|,
#
# rather than the RF variable-importance statistic in the original
# Maxway implementation. The randomization distribution is generated by
# drawing standardized residuals r_X^* ~ N(0, I_n).
# ------------------------------------------------------------
dCRT_res_cov_randomization_with_stat <- function(res_x, res_y, M = 500) {
  res_x <- as.numeric(res_x)
  res_y <- as.numeric(res_y)

  sx <- stats::sd(res_x)
  if (!is.finite(sx) || sx <= 0) stop("sd(res_x) is zero or non-finite.")
  res_x <- res_x / sx

  n <- length(res_x)

  T_signed <- mean(res_x * res_y)
  T_hat <- abs(T_signed)

  res_x_sample <- MASS::mvrnorm(M, rep(0, n), diag(rep(1, n)))
  T_sample_signed <- as.vector((res_x_sample %*% res_y) / n)
  T_sample <- abs(T_sample_signed)

  p_value <- (1 + sum(T_sample >= T_hat)) / (1 + M)
  
  # Summarize the randomization-null distribution for diagnostic purposes.
  T_null_mean <- mean(T_sample, na.rm = TRUE)
  T_null_sd <- stats::sd(T_sample, na.rm = TRUE)
  T_centered <- T_hat - T_null_mean
  T_z <- if (is.finite(T_null_sd) && T_null_sd > 0) T_centered / T_null_sd else NA_real_

  list(
    p_value = as.numeric(p_value),
    T_hat = as.numeric(T_hat),
    T_signed = as.numeric(T_signed),
    T_null_mean = as.numeric(T_null_mean),
    T_null_sd = as.numeric(T_null_sd),
    T_centered = as.numeric(T_centered),
    T_z = as.numeric(T_z)
  )
}

# ------------------------------------------------------------
# 4. Apply CRT and in-sample/CV Maxway CRT to one simulated dataset.
#
# X | Z is estimated using the unlabeled observations, while Y | Z is
# estimated using the labeled sample. Maxway then recalibrates the fitted
# X residuals using a low-dimensional summary extracted from the Y model.
#
# The auxiliary variable S is deliberately excluded from both CRT methods.
# ------------------------------------------------------------
maxwayCRT_insample_one <- function(label_data,
                                   unlabel_data,
                                   gbar_type = c("linear", "nonlinear"),
                                   x_col = "X",
                                   y_col = "Y",
                                   z_cols = NULL,
                                   seed = 1,
                                   M = 200,
                                   RF.num.trees = c(100, 100, 50, 20),
                                   lambda.seq = NULL,
                                   k = NULL,
                                   use_interaction = FALSE) {
  gbar_type <- match.arg(gbar_type)

  .force_single_threading()
  set.seed(seed)

  if (is.null(unlabel_data)) {
    stop("unlabel_data is NULL. MaxwayCRT/CRT here is intended to learn X|Z from unlabeled data.")
  }
  if (is.null(z_cols)) {
    z_cols <- setdiff(colnames(label_data), c(y_col, x_col, "S"))
  }
  
  # Extract labeled and unlabeled X and Z variables in the notation
  # expected by the original Maxway implementation.
  A_label <- as.numeric(label_data[[x_col]])
  y_label <- as.numeric(label_data[[y_col]])
  Z_label <- as_num_matrix(label_data[, z_cols, drop = FALSE])

  A_add <- as.numeric(unlabel_data[[x_col]])
  Z_add <- as_num_matrix(unlabel_data[, z_cols, drop = FALSE])

  p <- ncol(Z_label)
  N <- length(A_add)
  
  # Linear design: estimate both X | Z and Y | Z using lasso.
  if (gbar_type == "linear") {
    # Same default k as Maxway_CRT() for Gaussian_lasso.
    if (is.null(k)) {
      k <- as.integer(min(2 * log(p), p - 1, N^(1 / 3)))
    }
    k <- max(1L, min(as.integer(k), p))

    # X | Z is learned from unlabeled data only.
    set.seed(seed)
    fit_x <- fit_gaussian(
      Z_train = Z_add,
      y_train = A_add,
      Z_test = Z_label,
      y_test = A_label,
      model = "linear",
      lambda.seq = lambda.seq
    )
    res_x_lab <- fit_x$res_test
    res_x_add <- fit_x$res_train

    # Y | Z is learned from labeled data.
    set.seed(seed)
    fit_y <- fit_gaussian(
      Z_train = Z_label,
      y_train = y_label,
      Z_test = Z_add,
      y_test = NULL,
      model = "linear",
      lambda.seq = lambda.seq
    )
    res_y_lab <- fit_y$res_train
    beta_y <- fit_y$coef
    gZ_label <- fit_y$pred_train
    gZ_add <- fit_y$pred_test
    
    # Select the k covariates with the largest absolute coefficients
    # in the fitted Y regression.
    beta_sort <- sort(abs(as.vector(beta_y)), decreasing = TRUE, index.return = TRUE)
    index_gZ <- beta_sort$ix[seq_len(k)]

    Z_label_sub <- Z_label[, index_gZ, drop = FALSE]
    Z_add_sub <- Z_add[, index_gZ, drop = FALSE]
    
    # Construct the low-dimensional outcome-informed summary used
    # for Maxway calibration.
    gZ_label <- orthogonalize(gZ_label, Z_label_sub)
    gZ_add <- orthogonalize(gZ_add, Z_add_sub)

    # Maxway calibration: adjust fitted X residuals against the low-dimensional g(Z).
    cal_fit <- cal_fit_gaussian(gZ_add, res_x_add, gZ_label, res_x_lab)
    res_x_lab_cal <- cal_fit$res_test
    
    # Compute the original CRT statistic and the calibrated Maxway statistic.
    crt <- dCRT_res_with_stat(
      res_x = res_x_lab,
      res_y = res_y_lab,
      Z_sub = Z_label_sub,
      k = k,
      d.interaction = use_interaction
    )
    maxway <- dCRT_res_with_stat(
      res_x = res_x_lab_cal,
      res_y = res_y_lab,
      Z_sub = Z_label_sub,
      k = k,
      d.interaction = use_interaction
    )

    return(list(
      CRT = crt,
      MaxwayCRT = maxway,
      model_x = "Gaussian_lasso",
      model_y = "Gaussian_lasso",
      k = k
    ))
  }

  if (gbar_type == "nonlinear") {
    # Nonlinear design: estimate X | Z and Y | Z using random forests.
    if (is.null(k)) {
      k <- as.integer(2 * log(p))
    }
    k <- max(1L, min(as.integer(k), p))

    # X | Z is learned from unlabeled data only.
    set.seed(seed)
    fit_x <- fit_gaussian(
      Z_train = Z_add,
      y_train = A_add,
      Z_test = Z_label,
      y_test = A_label,
      model = "RF",
      RF.num.trees = RF.num.trees[1],
      CV = TRUE
    )
    res_x_lab <- fit_x$res_test
    res_x_add <- fit_x$res_train

    # Y | Z is learned from labeled data.
    set.seed(seed)
    fit_y <- fit_gaussian(
      Z_train = Z_label,
      y_train = y_label,
      Z_test = Z_add,
      y_test = NULL,
      model = "RF",
      RF.num.trees = RF.num.trees[2],
      CV = TRUE
    )
    res_y_lab <- fit_y$res_train
    gZ_label <- fit_y$pred_train
    gZ_add <- fit_y$pred_test
    
    # Select the k most important Z variables according to the fitted
    # random-forest outcome model.
    RF_y_importance <- as.vector(fit_y$model$importance)
    imp_sort <- sort(RF_y_importance, decreasing = TRUE, index.return = TRUE)
    index_gZ <- imp_sort$ix[seq_len(k)]

    Z_label_sub <- Z_label[, index_gZ, drop = FALSE]
    Z_add_sub <- Z_add[, index_gZ, drop = FALSE]
    
    # Combine the fitted Y regression with the selected important covariates
    # to form the Maxway calibration features.
    gZ_label <- cbind(gZ_label, Z_label_sub)
    gZ_add <- cbind(gZ_add, Z_add_sub)

    # Maxway calibration: RF calibration by default, matching Maxway_CRT().
    cal_fit <- cal_fit_gaussian(
      gZ_train = gZ_add,
      y_train = res_x_add,
      gZ_test = gZ_label,
      y_test = res_x_lab,
      model = "RF",
      RF.num.trees = RF.num.trees[3]
    )
    res_x_lab_cal <- cal_fit$res_test

    # Modified nonlinear statistic: residual product, not RF variable importance.
    crt <- dCRT_res_cov_randomization_with_stat(
      res_x = res_x_lab,
      res_y = res_y_lab,
      M = M
    )
    maxway <- dCRT_res_cov_randomization_with_stat(
      res_x = res_x_lab_cal,
      res_y = res_y_lab,
      M = M
    )

    return(list(
      CRT = crt,
      MaxwayCRT = maxway,
      model_x = "Gaussian_RF",
      model_y = "Gaussian_RF",
      k = k
    ))
  }
}

# ------------------------------------------------------------
# 5. Run Monte Carlo replications for one signal level c.
#
# Sparse coefficient vectors are drawn once and shared across all
# replications. Each replication generates a new labeled/unlabeled
# dataset, applies CRT and Maxway CRT, and stores the resulting
# p-values and test-statistic diagnostics.
#
# Replications may be evaluated sequentially or in parallel.
# ------------------------------------------------------------
run_maxwayCRT_replications_one_c <- function(n_rep = 2000,
                                             c_coef = 0,
                                             n_label = 300,
                                             n_unlabel = 2000,
                                             p_z = 30,
                                             gbar_type = c("linear", "nonlinear"),
                                             rho = 0.5,
                                             sx_gamma = 0,
                                             sx_power = 1,
                                             alpha = 0.05,
                                             seed = 123,
                                             beta_start = 1,
                                             beta_nnz = 5,
                                             beta_seed = 999,
                                             M = 200,
                                             RF.num.trees = c(100, 100, 50, 20),
                                             lambda.seq = NULL,
                                             k = NULL,
                                             use_interaction = FALSE,
                                             show_progress = TRUE,
                                             use_progress_bar = TRUE,
                                             parallel = TRUE,
                                             n_cores = max(1, parallel::detectCores() - 1),
                                             save_prefix = NULL) {
  gbar_type <- match.arg(gbar_type)
  
  # Validate settings required by the nonlinear data-generating mechanism
  # and the random-forest implementation.
  if (gbar_type == "nonlinear" && p_z < 22) {
    stop("For the nonlinear generator in generate_data.R, please use p_z >= 22.")
  }
  if (length(RF.num.trees) != 4) {
    stop("RF.num.trees must be a vector of length 4: c(X|Z, Y|Z, Maxway calibration, unused_statistic_slot). In this modified script, RF.num.trees[4] is not used because the nonlinear statistic is residual product.")
  }

  z_names <- paste0("Z", seq_len(p_z))
  if (is.null(beta_seed)) beta_seed <- seed
  
  # Generate sparse coefficient vectors once and hold them fixed across
  # all Monte Carlo replications, matching the main simulation design.
  # beta_g corresponds to beta_y in the notation of the paper.
  set.seed(beta_seed)
  beta_x <- make_sparse_beta(p_z, beta_start, beta_nnz, z_names)
  beta_s <- make_sparse_beta(p_z, beta_start + 3, beta_nnz, z_names)
  beta_g <- make_sparse_beta(p_z, beta_start + 3, beta_nnz, z_names)
  
  # Run one complete data-generation and testing replication.
  one_replication <- function(r) {
    tryCatch({
      .force_single_threading()
      rep_seed <- seed + r
      
      # Generate a new labeled and unlabeled dataset from the selected DGP.
      dat <- simulate_ppicit_plm(
        n_label = n_label,
        n_unlabel = n_unlabel,
        p_z = p_z,
        gbar_type = gbar_type,
        c_coef = c_coef,
        rho = rho,
        sx_gamma = sx_gamma,
        sx_power = sx_power,
        beta_x = beta_x,
        beta_s = beta_s,
        beta_g = beta_g,
        seed = rep_seed
      )
      
      # Apply CRT and Maxway CRT to the generated dataset.
      ans <- maxwayCRT_insample_one(
        label_data = dat$label_data,
        unlabel_data = dat$unlabel_data,
        gbar_type = gbar_type,
        z_cols = dat$truth$z_cols,
        seed = rep_seed,
        M = M,
        RF.num.trees = RF.num.trees,
        lambda.seq = lambda.seq,
        k = k,
        use_interaction = use_interaction
      )
      
      # Retain p-values together with observed and randomization-null
      # statistic summaries for both procedures.
      list(
        error = FALSE,
        error_msg = NA_character_,
        k = ans$k,
        MaxwayCRT_p_value = ans$MaxwayCRT$p_value,
        MaxwayCRT_reject = ans$MaxwayCRT$p_value <= alpha,
        MaxwayCRT_T_hat = ans$MaxwayCRT$T_hat,
        MaxwayCRT_T_signed = ans$MaxwayCRT$T_signed,
        MaxwayCRT_T_null_mean = ans$MaxwayCRT$T_null_mean,
        MaxwayCRT_T_null_sd = ans$MaxwayCRT$T_null_sd,
        MaxwayCRT_T_centered = ans$MaxwayCRT$T_centered,
        MaxwayCRT_T_z = ans$MaxwayCRT$T_z,
        CRT_p_value = ans$CRT$p_value,
        CRT_reject = ans$CRT$p_value <= alpha,
        CRT_T_hat = ans$CRT$T_hat,
        CRT_T_signed = ans$CRT$T_signed,
        CRT_T_null_mean = ans$CRT$T_null_mean,
        CRT_T_null_sd = ans$CRT$T_null_sd,
        CRT_T_centered = ans$CRT$T_centered,
        CRT_T_z = ans$CRT$T_z
      )
    }, error = function(e) {
      list(
        error = TRUE,
        error_msg = paste0("Replication ", r, " failed: ", conditionMessage(e)),
        k = NA_integer_,
        MaxwayCRT_p_value = NA_real_,
        MaxwayCRT_reject = NA,
        MaxwayCRT_T_hat = NA_real_,
        MaxwayCRT_T_signed = NA_real_,
        MaxwayCRT_T_null_mean = NA_real_,
        MaxwayCRT_T_null_sd = NA_real_,
        MaxwayCRT_T_centered = NA_real_,
        MaxwayCRT_T_z = NA_real_,
        CRT_p_value = NA_real_,
        CRT_reject = NA,
        CRT_T_hat = NA_real_,
        CRT_T_signed = NA_real_,
        CRT_T_null_mean = NA_real_,
        CRT_T_null_sd = NA_real_,
        CRT_T_centered = NA_real_,
        CRT_T_z = NA_real_
      )
    })
  }

  r_grid <- seq_len(n_rep)
  
  # Run Monte Carlo replications in parallel when multiple cores are requested.
  if (parallel && n_cores > 1) {
    n_cores <- min(n_cores, n_rep)

    if (show_progress) {
      cat("Running", n_rep, "MaxwayCRT/CRT replications on", n_cores, "cores ...\n")
      cat("RF/BLAS/OpenMP threads are set to 1 inside each worker.\n")
      cat("Note: base R parallel does not show a real-time progress bar on Windows.\n")
    }

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # Initialize each worker with single-threaded numerical libraries
    # and the packages required by the simulation.
    parallel::clusterEvalQ(cl, {
      Sys.setenv(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS = "1",
        RCPP_PARALLEL_NUM_THREADS = "1"
      )
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
        try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
      }
      suppressPackageStartupMessages(library(glmnet))
      suppressPackageStartupMessages(library(randomForest))
      suppressPackageStartupMessages(library(MASS))
      if (requireNamespace("CompQuadForm", quietly = TRUE)) {
        suppressPackageStartupMessages(library(CompQuadForm))
      }
      NULL
    })
    
    # Export the data-generation and CRT/Maxway helper functions to the workers.
    export_env <- environment(run_maxwayCRT_replications_one_c)
    parallel::clusterExport(
      cl,
      varlist = c(
        ".force_single_threading",
        "simulate_ppicit_plm",
        "fit_gaussian",
        "fit_binary",
        "cal_fit_gaussian",
        "cal_fit_binary",
        "orthogonalize",
        "glogit",
        "logit",
        "as_num_matrix",
        "dCRT_res_with_stat",
        "dCRT_res_cov_randomization_with_stat",
        "maxwayCRT_insample_one"
      ),
      envir = export_env
    )

    # Load balancing is useful because nonlinear RF randomization can vary in runtime.
    # Reproducibility is preserved because each replication uses seed + r.
    rep_list <- parallel::parLapplyLB(cl, r_grid, one_replication)

  } else {
    if (show_progress && use_progress_bar) {
      pb <- utils::txtProgressBar(min = 0, max = n_rep, style = 3)
      on.exit(close(pb), add = TRUE)
    }

    rep_list <- vector("list", n_rep)
    for (r in r_grid) {
      if (show_progress && !use_progress_bar) {
        cat("Running replication", r, "of", n_rep, "\n")
      }
      rep_list[[r]] <- one_replication(r)
      if (show_progress && use_progress_bar) {
        utils::setTxtProgressBar(pb, r)
      }
    }
  }
  
  # Convert the replication-level result lists into typed vectors.
  get_col <- function(name, mode = c("numeric", "logical", "integer", "character")) {
    mode <- match.arg(mode)
    if (mode == "numeric") {
      as.numeric(vapply(rep_list, function(x) x[[name]], numeric(1)))
    } else if (mode == "logical") {
      as.logical(vapply(rep_list, function(x) x[[name]], logical(1)))
    } else if (mode == "integer") {
      as.integer(vapply(rep_list, function(x) x[[name]], integer(1)))
    } else {
      as.character(vapply(rep_list, function(x) x[[name]], character(1)))
    }
  }
  
  # Identify failed replications and report their error messages.
  has_error <- get_col("error", "logical")
  error_msg <- get_col("error_msg", "character")
  error_msg <- error_msg[!is.na(error_msg) & nzchar(error_msg)]
  if (length(error_msg) > 0) {
    message(paste(error_msg, collapse = "\n"))
  } 
  
  # Extract test results and randomization-statistic diagnostics
  # for Maxway CRT and the original CRT.
  k_used <- get_col("k", "integer")
  maxway_p <- get_col("MaxwayCRT_p_value", "numeric")
  maxway_reject <- get_col("MaxwayCRT_reject", "logical")
  maxway_T <- get_col("MaxwayCRT_T_hat", "numeric")
  maxway_T_signed <- get_col("MaxwayCRT_T_signed", "numeric")
  maxway_T_null_mean <- get_col("MaxwayCRT_T_null_mean", "numeric")
  maxway_T_null_sd <- get_col("MaxwayCRT_T_null_sd", "numeric")
  maxway_T_centered <- get_col("MaxwayCRT_T_centered", "numeric")
  maxway_T_z <- get_col("MaxwayCRT_T_z", "numeric")

  crt_p <- get_col("CRT_p_value", "numeric")
  crt_reject <- get_col("CRT_reject", "logical")
  crt_T <- get_col("CRT_T_hat", "numeric")
  crt_T_signed <- get_col("CRT_T_signed", "numeric")
  crt_T_null_mean <- get_col("CRT_T_null_mean", "numeric")
  crt_T_null_sd <- get_col("CRT_T_null_sd", "numeric")
  crt_T_centered <- get_col("CRT_T_centered", "numeric")
  crt_T_z <- get_col("CRT_T_z", "numeric")
  
  # Restrict Monte Carlo summaries to replications completed successfully
  # for both procedures.
  valid <- !has_error &
    !is.na(maxway_p) & !is.na(crt_p) &
    !is.na(maxway_T) & !is.na(crt_T) &
    !is.na(maxway_reject) & !is.na(crt_reject)
  
  # Summarize rejection rates, p-values, and test-statistic behavior
  # across successful Monte Carlo replications.
  summary <- data.frame(
    method = c("MaxwayCRT_insample", "CRT"),
    gbar_type = gbar_type,
    c_coef = c_coef,
    n_label = n_label,
    n_unlabel_used_for_X_given_Z = n_unlabel,
    p_z = p_z,
    sx_gamma = sx_gamma,
    sx_power = sx_power,
    sx_lambda = sx_gamma / (n_label^sx_power),
    alpha = alpha,
    n_rep = n_rep,
    n_success = sum(valid),
    n_error = sum(has_error),
    parallel = parallel,
    n_cores = ifelse(parallel && n_cores > 1, n_cores, 1),
    RF_threads_per_replication = 1,
    rejection_rate = c(
      safe_mean(maxway_reject[valid]),
      safe_mean(crt_reject[valid])
    ),
    mean_p_value = c(
      safe_mean(maxway_p[valid]),
      safe_mean(crt_p[valid])
    ),
    median_p_value = c(
      safe_median(maxway_p[valid]),
      safe_median(crt_p[valid])
    ),
    mean_T_hat = c(
      safe_mean(maxway_T[valid]),
      safe_mean(crt_T[valid])
    ),
    sd_T_hat = c(
      safe_sd(maxway_T[valid]),
      safe_sd(crt_T[valid])
    ),
    mean_T_signed = c(
      safe_mean(maxway_T_signed[valid]),
      safe_mean(crt_T_signed[valid])
    ),
    sd_T_signed = c(
      safe_sd(maxway_T_signed[valid]),
      safe_sd(crt_T_signed[valid])
    ),
    mean_T_null_mean = c(
      safe_mean(maxway_T_null_mean[valid]),
      safe_mean(crt_T_null_mean[valid])
    ),
    mean_T_null_sd = c(
      safe_mean(maxway_T_null_sd[valid]),
      safe_mean(crt_T_null_sd[valid])
    ),
    mean_T_centered = c(
      safe_mean(maxway_T_centered[valid]),
      safe_mean(crt_T_centered[valid])
    ),
    sd_T_centered = c(
      safe_sd(maxway_T_centered[valid]),
      safe_sd(crt_T_centered[valid])
    ),
    mean_T_z = c(
      safe_mean(maxway_T_z[valid]),
      safe_mean(crt_T_z[valid])
    ),
    sd_T_z = c(
      safe_sd(maxway_T_z[valid]),
      safe_sd(crt_T_z[valid])
    ),
    mean_k = c(safe_mean(k_used[valid]), safe_mean(k_used[valid])),
    stringsAsFactors = FALSE
  )
  
  # Retain the complete replication-level results for diagnostics
  # and optional external analysis.
  details <- data.frame(
    replication = seq_len(n_rep),
    valid = valid,
    error = has_error,
    error_message = ifelse(has_error, get_col("error_msg", "character"), NA_character_),
    k = k_used,
    MaxwayCRT_p_value = maxway_p,
    MaxwayCRT_reject = maxway_reject,
    MaxwayCRT_T_hat = maxway_T,
    MaxwayCRT_T_signed = maxway_T_signed,
    MaxwayCRT_T_null_mean = maxway_T_null_mean,
    MaxwayCRT_T_null_sd = maxway_T_null_sd,
    MaxwayCRT_T_centered = maxway_T_centered,
    MaxwayCRT_T_z = maxway_T_z,
    CRT_p_value = crt_p,
    CRT_reject = crt_reject,
    CRT_T_hat = crt_T,
    CRT_T_signed = crt_T_signed,
    CRT_T_null_mean = crt_T_null_mean,
    CRT_T_null_sd = crt_T_null_sd,
    CRT_T_centered = crt_T_centered,
    CRT_T_z = crt_T_z,
    stringsAsFactors = FALSE
  )
  
  # Optionally save the summary and replication-level results as CSV files.
  if (!is.null(save_prefix)) {
    utils::write.csv(summary, paste0(save_prefix, "_summary.csv"), row.names = FALSE)
    utils::write.csv(details, paste0(save_prefix, "_details.csv"), row.names = FALSE)
  }
  
  # Return the compact summary, replication-level details,
  # and complete simulation settings.
  out <- list(
    summary = summary,
    rejection_rate = summary[, c("method", "rejection_rate")],
    details = details,
    settings = list(
      n_rep = n_rep,
      c_coef = c_coef,
      n_label = n_label,
      n_unlabel = n_unlabel,
      p_z = p_z,
      gbar_type = gbar_type,
      rho = rho,
      sx_gamma = sx_gamma,
      sx_power = sx_power,
      sx_lambda = sx_gamma / (n_label^sx_power),
      alpha = alpha,
      seed = seed,
      beta_seed = beta_seed,
      beta_x = beta_x,
      beta_s = beta_s,
      beta_g = beta_g,
      M = M,
      RF.num.trees = RF.num.trees,
      lambda.seq = lambda.seq,
      k = k,
      use_interaction = use_interaction,
      parallel = parallel,
      n_cores = ifelse(parallel && n_cores > 1, n_cores, 1),
      RF_threads_per_replication = 1
    )
  )

  class(out) <- "maxwayCRT_sim_one_c"
  return(out)
}

# ------------------------------------------------------------
# 6. Run CRT and Maxway CRT over a grid of signal strengths.
#
# The same Monte Carlo configuration is used for each c value, and the
# method-level summaries are combined into one output table.
# ------------------------------------------------------------
run_maxwayCRT_grid <- function(c_grid,
                               gbar_type = c("linear", "nonlinear"),
                               save_prefix = NULL,
                               ...) {
  gbar_type <- match.arg(gbar_type)
  out <- vector("list", length(c_grid))
  names(out) <- paste0("c_", c_grid)
  
  # Run the one-c simulation separately for each signal strength.
  for (i in seq_along(c_grid)) {
    c0 <- c_grid[i]
    prefix_i <- NULL
    if (!is.null(save_prefix)) {
      prefix_i <- paste0(save_prefix, "_", gbar_type, "_c", gsub("\\.", "p", as.character(c0)))
    }

    out[[i]] <- run_maxwayCRT_replications_one_c(
      c_coef = c0,
      gbar_type = gbar_type,
      save_prefix = prefix_i,
      ...
    )
  }
  
  # Combine summaries across all signal levels.
  summary_all <- do.call(rbind, lapply(out, function(x) x$summary))
  if (!is.null(save_prefix)) {
    utils::write.csv(summary_all, paste0(save_prefix, "_", gbar_type, "_summary_all.csv"), row.names = FALSE)
  }

  list(summary_all = summary_all, results = out)
}

# ============================================================
# 7. Main simulation settings
# ============================================================

# Linear design: use lasso-based nuisance estimation for c in
# {0, 0.1, 0.15, 0.2}.
c_grid <- c(0, 0.1, 0.15, 0.2)
res_linear_grid <- run_maxwayCRT_grid(
  c_grid = c_grid,
  gbar_type = "linear",
  n_rep = 2000,
  n_label = 200, #200,300
  n_unlabel = 2000, #2000,10000
  p_z = 30,
  sx_gamma = 0,
  sx_power = 1,
  alpha = 0.05,
  seed = 123,
  beta_start = 1,
  beta_nnz = 5,
  beta_seed = 999,
  M = 200,
  parallel = TRUE,
  n_cores = 20
)
print(res_linear_grid$summary_all)

# ------------------------------------------------------------
# Alternative: nonlinear design
# ------------------------------------------------------------
# Uncomment this block to run the same signal grid using
# random-forest nuisance estimation and Maxway calibration.

# c_grid <- c(0, 0.1, 0.15, 0.2)
# res_linear_grid <- run_maxwayCRT_grid(
#   c_grid = c_grid,
#   gbar_type = "nonlinear",
#   n_rep = 2000,
#   n_label = 200, #200,300
#   n_unlabel = 2000, #2000,10000
#   p_z = 30,
#   sx_gamma = 0,
#   sx_power = 1,
#   alpha = 0.05,
#   seed = 123,
#   beta_start = 1,
#   beta_nnz = 5,
#   beta_seed = 999,
#   M = 200,
#   RF.num.trees = c(100, 100, 50, 20),
#   parallel = TRUE,
#   n_cores = 20
# )
# print(res_linear_grid$summary_all)

