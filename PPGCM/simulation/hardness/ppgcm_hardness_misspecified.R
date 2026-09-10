# ============================================================
# Oracle hardness simulation for PPGCM under a misspecified
# prediction model
#
# Purpose:
#   This script implements the oracle simulation used to compare
#   PPGCM with the classical GCM when the prediction model is
#   deliberately misspecified. The true nuisance functions h(Z) and
#   g(Z) are known and evaluated without estimation error, while the
#   prediction function f(X,Z) used by PPGCM differs from the true
#   conditional mean E[Y | X,Z].
#
# Data-generating model:
#   X = Z^T beta_x + eps_x,
#   Y = c X + Z^T beta_g + eps_y,
#   where Z has an autoregressive covariance structure and eps_x,
#   eps_y are independent Gaussian errors.
#
# Main inputs:
#   n_label   : number of labeled observations (X, Y, Z)
#   n_unlabel : number of unlabeled observations (X, Z)
#   c_coef    : signal coefficient c
#   local_t   : optional local-alternative parameter, with
#               c = local_t / sqrt(n_label)
#   p_z       : dimension of Z
#   n_rep     : number of Monte Carlo replications
#
# Main outputs:
#   The simulation functions return empirical rejection rates, estimated
#   asymptotic variances and weights, their theoretical counterparts,
#   replication-level results, and the simulation configuration for:
#     1. PPGCM with the estimated optimal weight;
#     2. PPGCM with fixed weight w = 1;
#     3. the classical GCM.
#
# No external R packages are required.
# ============================================================

# ------------------------------------------------------------
# Generate a sparse coefficient vector for the linear model.
# Nonzero coefficients are placed in a consecutive block and drawn
# independently from a standard normal distribution.
# ------------------------------------------------------------
make_sparse_beta <- function(p_z,
                             beta_start = 1,
                             beta_nnz = 5,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (beta_start < 1) stop("beta_start must be at least 1.")
  if (beta_start + beta_nnz - 1 > p_z) {
    stop("beta_start + beta_nnz - 1 cannot exceed p_z.")
  }
  beta <- rep(0, p_z)
  idx <- beta_start:(beta_start + beta_nnz - 1)
  beta[idx] <- stats::rnorm(beta_nnz)
  names(beta) <- paste0("Z", seq_len(p_z))
  beta
}

# Construct the AR(1)-type covariance matrix for Z, with
# Cov(Z_j, Z_k) = rho^{|j-k|}.
make_ar_cov <- function(p_z, rho = 0.5) {
  outer(seq_len(p_z), seq_len(p_z), function(i, j) rho^abs(i - j))
}

# Compute the quadratic/bilinear form a^T Sigma b used in the
# analytical variance calculations.
quad_form <- function(a, Sigma, b = a) {
  as.numeric(crossprod(a, Sigma %*% b))
}

# Generate multivariate Gaussian covariates Z with the specified
# autoregressive covariance structure.
generate_Z <- function(n, p_z, rho = 0.5) {
  Sigma <- make_ar_cov(p_z, rho)
  Z <- matrix(stats::rnorm(n * p_z), nrow = n, ncol = p_z) %*% chol(Sigma)
  colnames(Z) <- paste0("Z", seq_len(p_z))
  Z
}

# ------------------------------------------------------------
# Generate one labeled and unlabeled dataset from the linear model.
#
# The labeled sample contains (Y, X, Z), whereas the unlabeled sample
# contains only (X, Z). The true data-generating parameters are returned
# together with the simulated observations for subsequent oracle
# calculations.
# ------------------------------------------------------------
simulate_linear_noS_oracle_fminus <- function(n_label = 300,
                                              n_unlabel = 2000,
                                              p_z = 30,
                                              c_coef = 0,
                                              rho = 0.5,
                                              beta_x,
                                              beta_g,
                                              sigma_x = 1,
                                              sigma_y = 1,
                                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (n_label < 2) stop("n_label must be at least 2.")
  if (n_unlabel < 1) stop("n_unlabel must be positive for PPGCM.")
  if (length(beta_x) != p_z) stop("length(beta_x) must equal p_z.")
  if (length(beta_g) != p_z) stop("length(beta_g) must equal p_z.")
  
  # Generate the labeled sample from the specified linear model.
  Z_l <- generate_Z(n_label, p_z, rho)
  eps_x_l <- stats::rnorm(n_label, sd = sigma_x)
  eps_y_l <- stats::rnorm(n_label, sd = sigma_y)

  X_l <- as.numeric(Z_l %*% beta_x + eps_x_l)
  Y_l <- as.numeric(c_coef * X_l + Z_l %*% beta_g + eps_y_l)
  
  # Generate an independent unlabeled sample containing X and Z only.
  Z_u <- generate_Z(n_unlabel, p_z, rho)
  eps_x_u <- stats::rnorm(n_unlabel, sd = sigma_x)
  X_u <- as.numeric(Z_u %*% beta_x + eps_x_u)

  list(
    label = list(Y = Y_l, X = X_l, Z = Z_l),
    unlabel = list(X = X_u, Z = Z_u),
    truth = list(
      c_coef = c_coef,
      beta_x = beta_x,
      beta_g = beta_g,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      rho = rho
    )
  )
}

# ------------------------------------------------------------
# Compute the oracle GCM and prediction-powered score components
# under the misspecified prediction model.
#
# The true functions h(Z) and g(Z) are used, while fminus(X,Z)
# deliberately differs from the true conditional mean E[Y | X,Z].
# The function returns
#   A   = eps_x * eps_y,
#   B_l = eps_x * {fminus(X,Z) - g(Z)} for labeled observations,
#   B_u = the corresponding quantity for unlabeled observations.
# ------------------------------------------------------------
oracle_scores_noS_fminus <- function(dat) {
  c_coef <- dat$truth$c_coef
  beta_x <- dat$truth$beta_x
  beta_g <- dat$truth$beta_g

  Y <- dat$label$Y
  X <- dat$label$X
  Z <- dat$label$Z

  Xu <- dat$unlabel$X
  Zu <- dat$unlabel$Z

  # Evaluate the true h(Z) and g(Z), together with the deliberately
  # misspecified prediction function, on the labeled observations.
  h_l <- as.numeric(Z %*% beta_x)
  gbar_l <- as.numeric(Z %*% beta_g)
  g_l <- as.numeric(c_coef * h_l + gbar_l)
  
  # Misspecified prediction function used by PPGCM.
  fminus_l <- as.numeric(-0.2 * (c_coef * X + gbar_l))
  
  # Construct the labeled residuals and prediction component.
  eps_x_l <- X - h_l
  eps_y_l <- Y - g_l
  V_l <- fminus_l - g_l
  
  # Form the GCM score A and the labeled prediction-powered score B.
  A <- eps_x_l * eps_y_l
  B_l <- eps_x_l * V_l

  # Evaluate the corresponding quantities on the unlabeled observations
  # to construct the unlabeled prediction-powered score.
  h_u <- as.numeric(Zu %*% beta_x)
  gbar_u <- as.numeric(Zu %*% beta_g)
  g_u <- as.numeric(c_coef * h_u + gbar_u)

  fminus_u <- as.numeric(-0.2*(c_coef * Xu + gbar_u))

  eps_x_u <- Xu - h_u
  V_u <- fminus_u - g_u
  B_u <- eps_x_u * V_u

  list(A = A, B_l = B_l, B_u = B_u)
}

# ------------------------------------------------------------
# Compute the theoretical asymptotic variances and optimal weight
# under the misspecified prediction model.
#
# The Gaussian linear model gives closed-form expressions for Var(A),
# Var(B), and Cov(A,B). These quantities are used to obtain the
# theoretical variance of GCM, PPGCM (w=1), and PPGCM (opt), together
# with their fixed-t local asymptotic limits.
# ------------------------------------------------------------
theory_variance_noS_fminus <- function(c_coef,
                                       n_label,
                                       n_unlabel,
                                       beta_x,
                                       beta_g,
                                       rho = 0.5,
                                       sigma_x = 1,
                                       sigma_y = 1,
                                       eps = 1e-12) {
  p_z <- length(beta_x)
  if (length(beta_g) != p_z) stop("beta_x and beta_g must have the same length.")
  
  # Compute the population variances and covariance of the linear
  # components h(Z) and Z^T beta_g.
  Sigma <- make_ar_cov(p_z, rho)

  var_h <- quad_form(beta_x, Sigma)
  var_G <- quad_form(beta_g, Sigma)
  cov_hG <- quad_form(beta_x, Sigma, beta_g)
  
  # Variance of the Z-dependent component entering the misspecified score.
  var_D <- c_coef^2 * var_h + 2 * c_coef * cov_hG + var_G
  
  # Ratio of labeled to unlabeled sample sizes.
  tau <- n_label / n_unlabel
  
  # Population variance and covariance terms entering the PPGCM variance.
  var_A <- sigma_x^2 * sigma_y^2 + 2 * c_coef^2 * sigma_x^4
  var_B <- 2 * 0.2^2 * c_coef^2 * sigma_x^4 +
    1.2^2 * sigma_x^2 * var_D
  cov_AB <- -2 * 0.2 * c_coef^2 * sigma_x^4
  
  # Theoretical asymptotic variance for the fixed-weight choice w = 1.
  var_w1 <- var_A - 2 * cov_AB + (1 + tau) * var_B
  
  # Compute the variance-minimizing weight and the corresponding
  # optimal asymptotic variance.
  if (!is.finite(var_B) || var_B <= eps) {
    w_star <- 0
    var_opt <- var_A
  } else {
    w_star <- cov_AB / ((1 + tau) * var_B)
    var_opt <- var_A - cov_AB^2 / ((1 + tau) * var_B)
  }
  
  # Limiting asymptotic variances under the fixed-t local alternative
  # c_n = t / sqrt(n).
  local_limit_gcm <- sigma_x^2 * sigma_y^2
  local_limit_opt <- sigma_x^2 * sigma_y^2
  local_limit_w1 <- sigma_x^2 * sigma_y^2 +
    (1 + tau) * 1.2^2 * sigma_x^2 * var_G

  data.frame(
    method = c("PPGCM (opt)", "PPGCM (w=1)", "GCM"),
    theory_weight = c(w_star, 1, 0),
    theory_variance = c(var_opt, var_w1, var_A),
    theory_variance_local_limit = c(local_limit_opt, local_limit_w1, local_limit_gcm),
    var_A = rep(var_A, 3),
    var_B = rep(var_B, 3),
    cov_AB = rep(cov_AB, 3),
    var_D = rep(var_D, 3),
    var_G = rep(var_G, 3),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Apply the three oracle tests to one simulated dataset.
#
# Empirical moments of A and B are used to construct the GCM,
# PPGCM (w=1), and PPGCM (opt) statistics and variance estimators
# under the misspecified prediction model. Each method is evaluated
# using a two-sided normal test.
# ------------------------------------------------------------
oracle_test_one_rep_fminus <- function(dat,
                                       alpha = 0.05,
                                       eps = 1e-10) {
  sc <- oracle_scores_noS_fminus(dat)

  A <- sc$A
  B_l <- sc$B_l
  B_u <- sc$B_u

  n <- length(A)
  N <- length(B_u)
  tau <- n / N
  
  # Compute the empirical means and second-moment quantities needed
  # for the test statistics and variance estimators.
  mean_A <- mean(A)
  mean_Bl <- mean(B_l)
  mean_Bu <- mean(B_u)

  var_A <- mean((A - mean_A)^2)
  cov_AB <- mean((A - mean_A) * (B_l - mean_Bl))

  B_pool <- c(B_l, B_u)
  var_B_pool <- mean((B_pool - mean(B_pool))^2)

  if (!is.finite(var_A) || var_A <= eps) var_A <- eps

  # Classical GCM statistic and variance estimator.
  T_gcm <- mean_A
  V_gcm <- var_A

  # PPGCM with fixed weight w = 1.
  T_w1 <- mean_A + (mean_Bu - mean_Bl)
  V_w1 <- var_A - 2 * cov_AB + (1 + tau) * var_B_pool
  if (!is.finite(V_w1) || V_w1 <= eps) V_w1 <- eps

  # PPGCM with the estimated variance-minimizing weight.
  if (!is.finite(var_B_pool) || var_B_pool <= eps) {
    w_hat <- 0
    V_opt <- var_A
  } else {
    w_hat <- cov_AB / ((1 + tau) * var_B_pool)
    V_opt <- var_A - cov_AB^2 / ((1 + tau) * var_B_pool)
    if (!is.finite(V_opt) || V_opt <= eps) V_opt <- eps
  }
  T_opt <- mean_A + w_hat * (mean_Bu - mean_Bl)
  
  # Convert a statistic and variance estimate into a two-sided z-test result.
  make_row <- function(method, T_hat, V_hat, weight_hat) {
    z <- sqrt(n) * T_hat / sqrt(V_hat)
    p_value <- 2 * (1 - stats::pnorm(abs(z)))
    data.frame(
      method = method,
      T_hat = T_hat,
      variance_hat = V_hat,
      weight_hat = weight_hat,
      z_score = z,
      p_value = p_value,
      reject = as.integer(p_value < alpha),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    make_row("PPGCM (opt)", T_opt, V_opt, w_hat),
    make_row("PPGCM (w=1)", T_w1, V_w1, 1),
    make_row("GCM", T_gcm, V_gcm, 0)
  )
}

# ------------------------------------------------------------
# Summarize Monte Carlo results across replications.
#
# For each method, report empirical means and standard deviations of
# the statistic, variance estimator, and estimated weight, together
# with the theoretical benchmarks, rejection rate, and mean p-value.
# Additional population variance components specific to the
# misspecified prediction model are also retained.
# ------------------------------------------------------------
summarize_oracle_results_fminus <- function(results_long,
                                            theory_df,
                                            c_coef,
                                            local_t,
                                            n_label,
                                            n_unlabel,
                                            n_rep) {
  methods <- c("PPGCM (opt)", "PPGCM (w=1)", "GCM")
  
  # Compute method-specific Monte Carlo summaries and attach the
  # corresponding theoretical benchmarks.
  out <- do.call(rbind, lapply(methods, function(m) {
    x <- results_long[results_long$method == m, , drop = FALSE]
    th <- theory_df[theory_df$method == m, , drop = FALSE]

    data.frame(
      method = m,
      c_coef = c_coef,
      local_t = ifelse(is.null(local_t), NA_real_, local_t),
      n_label = n_label,
      n_unlabel = n_unlabel,
      n_rep = n_rep,

      mean_T_hat = mean(x$T_hat, na.rm = TRUE),
      sd_T_hat = stats::sd(x$T_hat, na.rm = TRUE),

      mean_variance_hat = mean(x$variance_hat, na.rm = TRUE),
      sd_variance_hat = stats::sd(x$variance_hat, na.rm = TRUE),

      theory_variance = th$theory_variance,
      theory_variance_local_limit = th$theory_variance_local_limit,

      mean_weight_hat = mean(x$weight_hat, na.rm = TRUE),
      sd_weight_hat = stats::sd(x$weight_hat, na.rm = TRUE),
      theory_weight = th$theory_weight,

      rejection_rate = mean(x$reject, na.rm = TRUE),
      mean_p_value = mean(x$p_value, na.rm = TRUE),

      var_A = th$var_A,
      var_B = th$var_B,
      cov_AB = th$cov_AB,
      var_D = th$var_D,
      var_G = th$var_G,

      stringsAsFactors = FALSE
    )
  }))

  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------
# Run the misspecified-prediction simulation for one signal level.
#
# The signal can be specified directly through c_coef or through a
# local parameter local_t, in which case c = local_t / sqrt(n_label).
# The same sparse coefficient vectors are held fixed across Monte Carlo
# replications, while new labeled and unlabeled samples are generated.
# ------------------------------------------------------------
run_hardness_oracle_noS_fminus_one_c <- function(n_rep = 1000,
                                                 c_coef = 0,
                                                 local_t = NULL,
                                                 n_label = 300,
                                                 n_unlabel = 2000,
                                                 p_z = 30,
                                                 rho = 0.5,
                                                 sigma_x = 1,
                                                 sigma_y = 1,
                                                 alpha = 0.05,
                                                 seed = 123,
                                                 beta_start = 1,
                                                 beta_nnz = 5,
                                                 beta_seed = 999,
                                                 show_progress = TRUE) {
  # Convert the local-alternative parameter into the sample-size-dependent signal.
  if (!is.null(local_t)) {
    c_coef <- local_t / sqrt(n_label)
  }
  
  # Generate fixed sparse coefficient vectors for X and Y.
  beta_x <- make_sparse_beta(
    p_z = p_z,
    beta_start = beta_start,
    beta_nnz = beta_nnz,
    seed = beta_seed
  )

  beta_g <- make_sparse_beta(
    p_z = p_z,
    beta_start = beta_start + 5,
    beta_nnz = beta_nnz,
    seed = beta_seed + 1
  )

  all_results <- vector("list", n_rep)

  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = n_rep, style = 3)
  }
  
  # Repeat data generation and testing independently across Monte Carlo replications.
  for (r in seq_len(n_rep)) {
    dat <- simulate_linear_noS_oracle_fminus(
      n_label = n_label,
      n_unlabel = n_unlabel,
      p_z = p_z,
      c_coef = c_coef,
      rho = rho,
      beta_x = beta_x,
      beta_g = beta_g,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      seed = seed + r
    )

    res_r <- oracle_test_one_rep_fminus(dat, alpha = alpha)
    res_r$replication <- r
    all_results[[r]] <- res_r

    if (show_progress) utils::setTxtProgressBar(pb, r)
  }

  if (show_progress) close(pb)
  
  # Combine replication-level results and compute the theoretical benchmarks.
  results_long <- do.call(rbind, all_results)

  theory_df <- theory_variance_noS_fminus(
    c_coef = c_coef,
    n_label = n_label,
    n_unlabel = n_unlabel,
    beta_x = beta_x,
    beta_g = beta_g,
    rho = rho,
    sigma_x = sigma_x,
    sigma_y = sigma_y
  )
  
  # Aggregate the simulation results into a method-level summary.
  summary <- summarize_oracle_results_fminus(
    results_long = results_long,
    theory_df = theory_df,
    c_coef = c_coef,
    local_t = local_t,
    n_label = n_label,
    n_unlabel = n_unlabel,
    n_rep = n_rep
  )
  
  # Return summary results, replication-level results, theoretical values,
  # and the full simulation configuration.
  list(
    summary = summary,
    results = results_long,
    theory = theory_df,
    config = list(
      n_rep = n_rep,
      c_coef = c_coef,
      local_t = local_t,
      n_label = n_label,
      n_unlabel = n_unlabel,
      p_z = p_z,
      rho = rho,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      alpha = alpha,
      seed = seed,
      beta_seed = beta_seed,
      beta_x = beta_x,
      beta_g = beta_g
    )
  )
}

# ------------------------------------------------------------
# Run the misspecified-prediction simulation over a grid of signal strengths.
#
# The grid can be specified either by fixed c values or by local-alternative
# parameters t. Results from all settings are combined into one summary.
# ------------------------------------------------------------
run_hardness_oracle_noS_fminus_grid <- function(c_values = c(0, 0.05, 0.10, 0.15),
                                                local_t_values = NULL,
                                                ...) {
  # Use local alternatives when local_t_values are supplied; otherwise
  # run the simulation over the specified fixed c values.
  if (!is.null(local_t_values)) {
    res_list <- lapply(local_t_values, function(tt) {
      run_hardness_oracle_noS_fminus_one_c(local_t = tt, ...)
    })
  } else {
    res_list <- lapply(c_values, function(cc) {
      run_hardness_oracle_noS_fminus_one_c(c_coef = cc, ...)
    })
  }

  summary_all <- do.call(rbind, lapply(res_list, function(x) x$summary))
  rownames(summary_all) <- NULL

  list(summary = summary_all, by_setting = res_list)
}

# ------------------------------------------------------------
# Run the fixed-t simulation over labeled sample sizes under the
# misspecified prediction model.
#
# For each n, set
#   c_n = fixed_t / sqrt(n)
# and choose the unlabeled sample size so that n/N is approximately
# equal to n_over_N. This function studies how local power, asymptotic
# variance, and the estimated optimal weight behave as n increases.
# ------------------------------------------------------------
run_hardness_oracle_noS_fminus_by_n <- function(n_label_values = c(100, 200, 300, 500, 800, 1000),
                                                fixed_t = 1,
                                                n_over_N = 0.15,
                                                n_rep = 1000,
                                                p_z = 30,
                                                rho = 0.5,
                                                sigma_x = 1,
                                                sigma_y = 1,
                                                alpha = 0.05,
                                                seed = 123,
                                                beta_start = 1,
                                                beta_nnz = 5,
                                                beta_seed = 999,
                                                independent_by_n = TRUE,
                                                seed_gap = 100000,
                                                show_progress = TRUE) {
  if (length(n_label_values) < 1) stop("n_label_values must be non-empty.")
  if (any(n_label_values < 2)) stop("Every n_label must be at least 2.")
  if (!is.finite(fixed_t) || fixed_t < 0) stop("fixed_t must be a nonnegative finite number.")
  if (!is.finite(n_over_N) || n_over_N <= 0) stop("n_over_N must be positive.")
  
  # Choose N for each n to keep the labeled-to-unlabeled ratio approximately fixed.
  n_label_values <- as.integer(n_label_values)
  n_unlabel_values <- pmax(1L, as.integer(round(n_label_values / n_over_N)))

  res_list <- vector("list", length(n_label_values))
  
  # Run the fixed-t simulation separately for each labeled sample size.
  for (idx in seq_along(n_label_values)) {
    n_now <- n_label_values[idx]
    N_now <- n_unlabel_values[idx]
    # Optionally use separated random-number streams for different sample sizes.
    seed_now <- if (independent_by_n) seed + idx * seed_gap else seed

    if (show_progress) {
      cat("\n============================================================\n")
      cat("Running fixed-t oracle no-S fminus simulation:\n")
      cat("  f0(X,Z) = -0.2 E(Y|X,Z)\n")
      cat("  n_label =", n_now, "\n")
      cat("  n_unlabel =", N_now, "\n")
      cat("  actual n/N =", round(n_now / N_now, 6), "\n")
      cat("  fixed t =", fixed_t, "\n")
      cat("  c_n = t / sqrt(n) =", fixed_t / sqrt(n_now), "\n")
      cat("============================================================\n")
    }

    res_list[[idx]] <- run_hardness_oracle_noS_fminus_one_c(
      n_rep = n_rep,
      local_t = fixed_t,
      n_label = n_now,
      n_unlabel = N_now,
      p_z = p_z,
      rho = rho,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      alpha = alpha,
      seed = seed_now,
      beta_start = beta_start,
      beta_nnz = beta_nnz,
      beta_seed = beta_seed,
      show_progress = show_progress
    )
  }
  
  # Combine results across n and record the target and realized n/N ratios.
  summary_all <- do.call(rbind, lapply(res_list, function(x) x$summary))
  summary_all$n_over_N_target <- n_over_N
  summary_all$n_over_N_actual <- summary_all$n_label / summary_all$n_unlabel
  summary_all$fixed_t <- fixed_t
  summary_all$c_coef <- fixed_t / sqrt(summary_all$n_label)
  
  # Reorder the main simulation descriptors to the front of the summary table.
  front_cols <- c(
    "method", "fixed_t", "c_coef", "n_label", "n_unlabel",
    "n_over_N_target", "n_over_N_actual", "n_rep"
  )
  other_cols <- setdiff(names(summary_all), front_cols)
  summary_all <- summary_all[, c(front_cols, other_cols), drop = FALSE]
  rownames(summary_all) <- NULL
  
  # Return the combined summary, results for each sample size,
  # and the simulation configuration.
  list(
    summary = summary_all,
    by_n = res_list,
    config = list(
      n_label_values = n_label_values,
      n_unlabel_values = n_unlabel_values,
      fixed_t = fixed_t,
      n_over_N = n_over_N,
      n_rep = n_rep,
      p_z = p_z,
      rho = rho,
      sigma_x = sigma_x,
      sigma_y = sigma_y,
      alpha = alpha,
      seed = seed,
      beta_seed = beta_seed,
      independent_by_n = independent_by_n
    )
  )
}

# ------------------------------------------------------------
# Extract a compact table for the fixed-t misspecified-prediction simulation.
# The retained columns focus on rejection rates, estimated and theoretical
# variances, optimal weights, and the population variance components
# specific to the misspecified prediction model.
# ------------------------------------------------------------
make_hardness_by_n_table_fminus <- function(res_by_n) {
  s <- res_by_n$summary
  keep <- c(
    "method", "fixed_t", "n_label", "n_unlabel", "n_over_N_actual",
    "c_coef", "mean_T_hat", "sd_T_hat",
    "mean_variance_hat", "sd_variance_hat",
    "theory_variance", "theory_variance_local_limit",
    "mean_weight_hat", "theory_weight", "rejection_rate",
    "var_A", "var_B", "cov_AB", "var_D", "var_G"
  )
  s[, keep, drop = FALSE]
}

