# ============================================================
# Prediction-Powered Generalized Covariance Measure (PPGCM)
#
# Purpose:
#   This script implements the PPGCM test for testing the conditional
#   independence hypothesis X ⟂ Y | Z using labeled observations
#   (X, Y, Z, S) together with optional unlabeled observations (X, Z, S).
#   For comparison, it also computes the classical cross-fitted GCM
#   and the PPGCM with fixed weight w = 1.
#
# Main input:
#   label_data   : labeled data containing Y, X, S, and Z
#   unlabel_data : optional unlabeled data containing X, S, and Z
#   K            : number of folds used for cross-fitting
#   fit_h        : regression method for h(Z) = E[X | Z]
#   fit_g        : regression method for g(Z) = E[Y | Z]
#   fit_f        : regression method for f(X,S,Z) = E[Y | X,S,Z]
#
# Main output:
#   ppicit_test() returns a list containing results for:
#     1. PPGCM with the estimated optimal weight;
#     2. PPGCM with fixed weight w = 1;
#     3. the classical cross-fitted GCM.
#   For each method, the test statistic, estimated asymptotic variance,
#   standard error, z-score, p-value, and rejection decision are reported.
#
# Required packages:
#   fit_* = "lasso" : install.packages("glmnet")
#   fit_* = "grf"   : install.packages("grf")
# ============================================================

# ------------------------------------------------------------
# Fit a nuisance regression model.
# The function provides a common interface for fitting either
# Lasso or generalized random forests and returns a prediction function.
# A constant model is used when no predictors are available or the
# response is essentially constant.
# ------------------------------------------------------------
.fit_regression <- function(y, X,
                            method = c("grf", "lasso"),
                            standardize_lasso = TRUE,
                            grf_num_trees = 200,
                            eps = 1e-12) {
  method <- match.arg(method)
  
  y <- as.numeric(y)
  X <- as.data.frame(X)
  x_names <- make.names(names(X), unique = TRUE)
  names(X) <- x_names
  
  n_train <- length(y)
  
  .n_new <- function(newX) {
    newX <- as.data.frame(newX)
    nrow(newX)
  }
  
  # Constant model fallback
  if (ncol(X) == 0 || stats::sd(y) < eps) {
    mu <- mean(y)
    return(list(
      method = "constant",
      predict = function(newX) rep(mu, .n_new(newX))
    ))
  }
  
  .prep_new <- function(newX) {
    newX <- as.data.frame(newX)
    if (ncol(newX) != length(x_names)) {
      stop("newX has a different number of columns from training X.")
    }
    names(newX) <- x_names
    newX
  }
  
  if (method == "lasso") {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package 'glmnet' is required for method = 'lasso'.")
    }
    nfolds <- min(5, n_train)
    if (nfolds < 3) stop("Lasso needs at least 3 training observations in each training fold.")
    fit <- glmnet::cv.glmnet(
      x = as.matrix(X), y = y,
      family = "gaussian", alpha = 1, nfolds = nfolds,
      standardize = standardize_lasso
    )
    return(list(
      method = "lasso", fit = fit,
      predict = function(newX) {
        newX <- .prep_new(newX)
        as.numeric(stats::predict(fit, newx = as.matrix(newX), s = "lambda.min"))
      }
    ))
  }
  
  if (method == "grf") {
    if (!requireNamespace("grf", quietly = TRUE)) {
      stop("Package 'grf' is required for method = 'grf'.")
    }
    fit <- grf::regression_forest(
      X = as.matrix(X), Y = y, num.trees = grf_num_trees, num.threads = 1
    )
    return(list(
      method = "grf", fit = fit,
      predict = function(newX) {
        newX <- .prep_new(newX)
        as.numeric(stats::predict(fit, newdata = as.matrix(newX))$predictions)
      }
    ))
  }
}

# ------------------------------------------------------------
# Main PPGCM testing function.
#
# The labeled sample is split into K folds. For each fold, nuisance
# functions are trained on the remaining labeled observations and
# evaluated out of fold. The resulting residual products are then
# combined across folds to construct the GCM and PPGCM statistics.
# Unlabeled observations are used only by PPGCM.
# ------------------------------------------------------------
ppicit_test <- function(label_data,
                        unlabel_data = NULL,
                        K = 5,
                        fit_h = c("grf", "lasso"),
                        fit_g = c("grf", "lasso"),
                        fit_f = c("grf", "lasso"),
                        x_col = "X",
                        y_col = "Y",
                        s_col = "S",
                        z_cols = NULL,
                        alpha = 0.05,
                        seed = NULL,
                        drop_missing = TRUE,
                        standardize_lasso = TRUE,
                        grf_num_trees = 200,
                        eps = 1e-8) {
  
  fit_h <- match.arg(fit_h)
  fit_g <- match.arg(fit_g)
  fit_f <- match.arg(fit_f)
  
  # Prepare the labeled sample: identify the required variables,
  # remove or check missing values, and convert all variables to numeric form.
  label_data <- as.data.frame(label_data)
  
  if (is.null(z_cols)) {
    z_cols <- setdiff(names(label_data), c(y_col, x_col, s_col))
  }
  
  label_req <- c(y_col, x_col, s_col, z_cols)
  miss_label <- setdiff(label_req, names(label_data))
  if (length(miss_label) > 0) {
    stop("label_data is missing columns: ", paste(miss_label, collapse = ", "))
  }
  
  label_data <- label_data[, label_req, drop = FALSE]
  
  if (drop_missing) {
    label_data <- label_data[stats::complete.cases(label_data), , drop = FALSE]
  } else {
    if (anyNA(label_data)) stop("label_data contains missing values.")
  }
  
  label_data <- data.frame(
    lapply(label_data, function(v) as.numeric(v)),
    check.names = FALSE
  )
  
  if (!all(is.finite(as.matrix(label_data)))) {
    stop("label_data contains non-finite values after numeric conversion.")
  }
  
  n <- nrow(label_data)
  if (n < 2) stop("label_data must contain at least 2 complete observations.")
  if (K < 2) stop("K must be at least 2 for cross-fitting.")
  if (K > n) stop("K cannot be larger than the labeled sample size.")
  
  # Extract the outcome, exposure, conditioning variables, and predictors
  # used by the prediction model f(X, S, Z).
  Y <- label_data[[y_col]]
  X <- label_data[[x_col]]
  Z <- label_data[, z_cols, drop = FALSE]
  F_label <- label_data[, c(x_col, s_col, z_cols), drop = FALSE]
  
  # ============================================================
  # Handle unlabeled data
  # The unlabeled sample contains X, S, and Z but does not require Y.
  # It contributes to the prediction-powered correction term.
  # ============================================================
  
  use_unlabel <- !is.null(unlabel_data) && nrow(as.data.frame(unlabel_data)) > 0
  N <- 0
  
  if (use_unlabel) {
    unlabel_data <- as.data.frame(unlabel_data)
    unlabel_req <- c(x_col, s_col, z_cols)
    
    miss_unlabel <- setdiff(unlabel_req, names(unlabel_data))
    if (length(miss_unlabel) > 0) {
      stop("unlabel_data is missing columns: ", paste(miss_unlabel, collapse = ", "))
    }
    
    unlabel_data <- unlabel_data[, unlabel_req, drop = FALSE]
    
    if (drop_missing) {
      unlabel_data <- unlabel_data[stats::complete.cases(unlabel_data), , drop = FALSE]
    } else {
      if (anyNA(unlabel_data)) stop("unlabel_data contains missing values.")
    }
    
    unlabel_data <- data.frame(
      lapply(unlabel_data, function(v) as.numeric(v)),
      check.names = FALSE
    )
    
    if (!all(is.finite(as.matrix(unlabel_data)))) {
      stop("unlabel_data contains non-finite values after numeric conversion.")
    }
    
    N <- nrow(unlabel_data)
    
    if (N == 0) {
      use_unlabel <- FALSE
    } else {
      X_unlabel <- unlabel_data[[x_col]]
      Z_unlabel <- unlabel_data[, z_cols, drop = FALSE]
      F_unlabel <- unlabel_data[, c(x_col, s_col, z_cols), drop = FALSE]
    }
  }
  
  # ============================================================
  # Cross-fitting folds
  # Randomly partition the labeled observations into K folds and initialize
  # storage for fold-specific statistics and out-of-fold quantities.
  # ============================================================
  
  if (!is.null(seed)) set.seed(seed)
  fold_id <- sample(rep(seq_len(K), length.out = n))
  folds <- split(seq_len(n), fold_id)
  
  fold_details <- data.frame(
    fold = seq_len(K),
    n_fold = NA_integer_,
    
    Abar = NA_real_,
    var_A = NA_real_,
    
    Bbar_label = NA_real_,
    Bbar_unlabel = NA_real_,
    cov_AB = NA_real_,
    var_B = NA_real_,
    
    weight = 0,
    T_gcm_fold = NA_real_,
    T_ppi_fold = NA_real_,
    T_ppi_w1_fold = NA_real_
  )
  
  A_oof <- rep(NA_real_, n)
  B_oof <- rep(NA_real_, n)
  
  if (use_unlabel) {
    B_unlabel_mat <- matrix(NA_real_, nrow = N, ncol = K)
  } else {
    B_unlabel_mat <- NULL
  }
  
  # ============================================================
  # Main cross-fitting loop
  # For each fold, nuisance functions are fitted on the training folds and
  # evaluated on the held-out fold. This produces out-of-fold estimates of
  # the GCM component A and the prediction-powered component B.
  # ============================================================
  
  for (k in seq_len(K)) {
    idx_test <- folds[[k]]
    idx_train <- setdiff(seq_len(n), idx_test)
    
    if (length(idx_train) < 2) {
      stop("Each training fold must contain at least 2 observations.")
    }
    
    # h(z) = E[X | Z]
    h_fit <- .fit_regression(
      y = X[idx_train],
      X = Z[idx_train, , drop = FALSE],
      method = fit_h,
      standardize_lasso = standardize_lasso,
      grf_num_trees = grf_num_trees
    )
    
    # g(z) = E[Y | Z]
    g_fit <- .fit_regression(
      y = Y[idx_train],
      X = Z[idx_train, , drop = FALSE],
      method = fit_g,
      standardize_lasso = standardize_lasso,
      grf_num_trees = grf_num_trees
    )
    
    # Compute out-of-fold residuals for X and Y.
    hx_label <- h_fit$predict(Z[idx_test, , drop = FALSE])
    gy_label <- g_fit$predict(Z[idx_test, , drop = FALSE])
    
    eps_x_label <- X[idx_test] - hx_label
    eps_y_label <- Y[idx_test] - gy_label
    
    # A_i = eps_x * eps_y
    # This is the GCM component.
    A <- eps_x_label * eps_y_label
    Abar <- mean(A)
    var_A <- mean((A - Abar)^2)
    
    A_oof[idx_test] <- A
    
    # GCM fold statistic
    T_gcm_fold <- Abar
    
    if (use_unlabel) {
      # Construct the prediction-powered component B using both labeled
      # and unlabeled observations.
      # f(x, s, z) = E[Y | X, S, Z]
      f_fit <- .fit_regression(
        y = Y[idx_train],
        X = F_label[idx_train, , drop = FALSE],
        method = fit_f,
        standardize_lasso = standardize_lasso,
        grf_num_trees = grf_num_trees
      )
      
      # Compute B on the held-out labeled observations.
      f_label <- f_fit$predict(F_label[idx_test, , drop = FALSE])
      # misspecified prediction model
      # f_label <- 0.6*f_fit$predict(F_label[idx_test, , drop = FALSE])
      V_label <- f_label - gy_label
      B_label <- eps_x_label * V_label
      
      Bbar_label <- mean(B_label)
      B_oof[idx_test] <- B_label
      
      # Compute the corresponding B values on the unlabeled observations.
      hx_unlabel <- h_fit$predict(Z_unlabel)
      gy_unlabel <- g_fit$predict(Z_unlabel)
      f_unlabel <- f_fit$predict(F_unlabel)
      # misspecified prediction model
      # f_unlabel <- 0.6*f_fit$predict(F_unlabel)
      
      eps_x_unlabel <- X_unlabel - hx_unlabel
      V_unlabel <- f_unlabel - gy_unlabel
      B_unlabel <- eps_x_unlabel * V_unlabel
      
      Bbar_unlabel <- mean(B_unlabel)
      B_unlabel_mat[, k] <- B_unlabel
      
      # Combined mean for Var(B)
      Bbar_combined <- (sum(B_label) + sum(B_unlabel)) /
        (length(B_label) + length(B_unlabel))
      
      cov_AB <- mean((A - Abar) * (B_label - Bbar_label))
      
      var_B <- (
        sum((B_label - Bbar_combined)^2) +
          sum((B_unlabel - Bbar_combined)^2)
      ) / (length(B_label) + length(B_unlabel))
      
      # Estimate the fold-specific variance-minimizing PPGCM weight.
      # The weight is set to zero if Var(B) is numerically degenerate.
      if (!is.finite(var_B) || var_B <= eps) {
        w_hat <- 0
      } else {
        w_hat <- cov_AB / ((1 + n / N) * var_B)
      }
      
      # PPGCM fold statistic
      T_ppi_fold <- Abar + w_hat * (Bbar_unlabel - Bbar_label)
      # PPGCM fold statistic with fixed weight = 1
      T_ppi_w1_fold <- Abar + (Bbar_unlabel - Bbar_label)
      
    } else {
      Bbar_label <- NA_real_
      Bbar_unlabel <- NA_real_
      cov_AB <- NA_real_
      var_B <- NA_real_
      w_hat <- 0
      
      # No unlabeled data: PPGCM degenerates to GCM.
      T_ppi_fold <- T_gcm_fold
      T_ppi_w1_fold <- T_gcm_fold
    }
    
    fold_details[k, ] <- list(
      fold = k,
      n_fold = length(idx_test),
      
      Abar = Abar,
      var_A = var_A,
      
      Bbar_label = Bbar_label,
      Bbar_unlabel = Bbar_unlabel,
      cov_AB = cov_AB,
      var_B = var_B,
      
      weight = w_hat,
      T_gcm_fold = T_gcm_fold,
      T_ppi_fold = T_ppi_fold,
      T_ppi_w1_fold = T_ppi_w1_fold
    )
  }
  
  # ============================================================
  # Pool out-of-fold quantities across folds
  # ============================================================
  # The out-of-fold A and B values are pooled to estimate the asymptotic
  # variances and the common optimal PPGCM weight used in the final test.
  
  A_pool <- A_oof
  
  # GCM statistic remains the original cross-fitted statistic
  T_gcm <- mean(fold_details$T_gcm_fold)
  
  # GCM variance now uses pooled A
  V_gcm <- mean((A_pool - mean(A_pool))^2)
  
  if (!is.finite(V_gcm) || V_gcm <= eps) {
    V_gcm <- eps
  }
  
  # Form the two-sided GCM test using the asymptotic normal approximation.
  se_gcm <- sqrt(V_gcm / n)
  z_gcm <- sqrt(n) * T_gcm / sqrt(V_gcm)
  p_gcm <- 2 * (1 - stats::pnorm(abs(z_gcm)))
  reject_gcm <- as.logical(p_gcm < alpha)
  
  
  # ============================================================
  # PPGCM with estimated optimal weight
  # ============================================================
  # Estimate a common optimal weight from the pooled out-of-fold quantities
  # and use it to construct the final PPGCM statistic and variance estimator.
  
  if (use_unlabel) {
    B_label_pool <- B_oof
    
    # Pool all fold-specific unlabeled B values.
    # This matches the fold-based construction because each fold has its own nuisance fit.
    B_unlabel_pool <- as.vector(B_unlabel_mat)
    
    # Cov(A, B) can only be estimated from labeled out-of-fold pairs
    cov_AB_pool <- mean(
      (A_pool - mean(A_pool)) *
        (B_label_pool - mean(B_label_pool))
    )
    
    # Var(B) is estimated from pooled labeled and unlabeled B values
    B_pool <- c(B_label_pool, B_unlabel_pool)
    var_B_pool <- mean((B_pool - mean(B_pool))^2)
    
    # Pooled optimal weight
    # Estimate the fold-specific variance-minimizing PPGCM weight.
    # The weight is set to zero if Var(B) is numerically degenerate.
    if (!is.finite(var_B_pool) || var_B_pool <= eps) {
      w_hat <- 0
    } else {
      w_hat <- cov_AB_pool / ((1 + n / N) * var_B_pool)
    }
    
    # Recompute the fold statistics using the common pooled weight.
    fold_details$weight <- w_hat
    fold_details$T_ppi_fold <- fold_details$Abar +
      w_hat * (fold_details$Bbar_unlabel - fold_details$Bbar_label)
    
    T_ppi <- mean(fold_details$T_ppi_fold)
    
    # Pooled empirical version of:
    # V_opt = Var(eps_x eps_y)
    #         - Cov(eps_x eps_y, eps_x V)^2 /
    #           ((1 + n/N) Var(eps_x V))
    V_ppi <- V_gcm -
      cov_AB_pool^2 / ((1 + n / N) * var_B_pool)
    
    if (!is.finite(V_ppi) || V_ppi <= eps) {
      V_ppi <- V_gcm
    }
    
  } else {
    w_hat <- 0
    cov_AB_pool <- NA_real_
    var_B_pool <- NA_real_
    
    fold_details$weight <- 0
    fold_details$T_ppi_fold <- fold_details$T_gcm_fold
    
    T_ppi <- T_gcm
    V_ppi <- V_gcm
  }
  
  # Form the two-sided PPGCM test with the estimated optimal weight.
  se_ppi <- sqrt(V_ppi / n)
  z_ppi <- sqrt(n) * T_ppi / sqrt(V_ppi)
  p_ppi <- 2 * (1 - stats::pnorm(abs(z_ppi)))
  reject_ppi <- as.logical(p_ppi < alpha)
  
  # ============================================================
  # PPGCM result with fixed weight = 1
  # This version uses the same prediction-powered correction but fixes
  # the weight at one instead of estimating the variance-minimizing weight.
  # ============================================================
  
  T_ppi_w1 <- mean(fold_details$T_ppi_w1_fold)
  
  if (use_unlabel) {
    B_label_pool <- B_oof
    
    # Pool all fold-specific unlabeled B values.
    # This matches the fold-based construction because each fold has its own nuisance fit.
    B_unlabel_pool <- as.vector(B_unlabel_mat)
    
    # Cov(A, B) can only be estimated from labeled out-of-fold pairs
    cov_AB_pool <- mean(
      (A_pool - mean(A_pool)) *
        (B_label_pool - mean(B_label_pool))
    )
    
    # Var(B) is estimated from pooled labeled and unlabeled B values
    B_pool <- c(B_label_pool, B_unlabel_pool)
    var_B_pool <- mean((B_pool - mean(B_pool))^2)
    
    # For a general fixed weight w:
    # V(w) = Var(A - w B_label) + (n / N) w^2 Var(B_unlabel)
    #      = Var(A) - 2w Cov(A, B) + (1 + n / N) w^2 Var(B)
    # Here w = 1.
    V_ppi_w1 <- V_gcm - 2 * cov_AB_pool + (1 + n / N) * var_B_pool
    
    if (!is.finite(V_ppi_w1) || V_ppi_w1 <= eps) {
      V_ppi_w1 <- max(V_gcm, eps)
    }
    
  } else {
    V_ppi_w1 <- V_gcm
  }
  
  se_ppi_w1 <- sqrt(V_ppi_w1 / n)
  z_ppi_w1 <- sqrt(n) * T_ppi_w1 / sqrt(V_ppi_w1)
  p_ppi_w1 <- 2 * (1 - stats::pnorm(abs(z_ppi_w1)))
  reject_ppi_w1 <- as.logical(p_ppi_w1 < alpha)
  
  # ============================================================
  # Summary table
  # Collect the main inferential results for the three procedures.
  # ============================================================
  
  method_summary <- data.frame(
    method = c("PPI-GCM", "PPI-GCM (w=1)", "GCM"),
    n_label = c(n, n, n),
    n_unlabel_used = c(
      ifelse(use_unlabel, N, 0),
      ifelse(use_unlabel, N, 0),
      0
    ),
    T_hat = c(T_ppi, T_ppi_w1, T_gcm),
    variance_hat = c(V_ppi, V_ppi_w1, V_gcm),
    se_hat = c(se_ppi, se_ppi_w1, se_gcm),
    z_score = c(z_ppi, z_ppi_w1, z_gcm),
    p_value = c(p_ppi, p_ppi_w1, p_gcm),
    weight_hat_mean = c(mean(fold_details$weight), 1, 0),
    reject = c(reject_ppi, reject_ppi_w1, reject_gcm)
  )
  
  # Return both the compact summary table and method-specific results.
  out <- list(
    summary = method_summary,
    
    ppi = list(
      T_hat = T_ppi,
      variance_hat = V_ppi,
      se_hat = se_ppi,
      z_score = z_ppi,
      p_value = p_ppi,
      reject = reject_ppi,
      weights = rep(w_hat, K),
      weight_hat_mean = w_hat,
      cov_AB_pool = cov_AB_pool,
      var_B_pool = var_B_pool
    ),
    
    gcm = list(
      T_hat = T_gcm,
      variance_hat = V_gcm,
      se_hat = se_gcm,
      z_score = z_gcm,
      p_value = p_gcm,
      reject = reject_gcm,
      weights = rep(0, K),
      weight_hat_mean = 0
    ),
    
    ppi_w1 = list(
      T_hat = T_ppi_w1,
      variance_hat = V_ppi_w1,
      se_hat = se_ppi_w1,
      z_score = z_ppi_w1,
      p_value = p_ppi_w1,
      reject = reject_ppi_w1,
      weights = rep(1, K),
      weight_hat_mean = 1
    )
  )
  
  class(out) <- "ppicit"
  return(out)
}
