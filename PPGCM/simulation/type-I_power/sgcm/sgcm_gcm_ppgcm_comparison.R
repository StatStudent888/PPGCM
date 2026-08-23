# ============================================================
# Simulation comparing GCM, SGCM, and PPGCM under H0
#
# Purpose:
#   This script compares the oracle GCM, cross-fitted GCM, oracle SGCM,
#   and optimal-weight PPGCM under a data-generating mechanism for which
#   X is conditionally independent of Y given Z, but the auxiliary variable
#   S is conditionally associated with both X and Y.
#
#   The simulation is designed to illustrate that directly adjusting Y
#   for S, as in SGCM, can destroy type-I error control even though H0 holds,
#   whereas GCM and PPGCM continue to target X ⟂ Y | Z.
#
# Data-generating mechanism:
#   Z ~ N_p(0, Sigma),  Sigma[j,k] = 0.5^|j-k|
#   X = Z^T beta_x + eps_x
#   Y = Z^T beta_y + eps_y
#   S = Z^T beta_s + rho_x eps_x + rho_y eps_y + eps_s
#
# where eps_x, eps_y, and eps_s are independent N(0,1).
#
# Under this model:
#   X ⟂ Y | Z,
# but S is conditionally associated with both X and Y when
# rho_x and rho_y are nonzero.
#
# Methods compared:
#   1. Oracle GCM;
#   2. cross-fitted GCM;
#   3. Oracle SGCM;
#   4. PPGCM with the estimated optimal weight.
#
# Main outputs:
#   raw_results     : replication-level test results;
#   summary_results : Monte Carlo rejection rates and summaries;
#   report_table    : formatted simulation table.
#
# Required:
#   - PPGCM.R
#   - glmnet
#
# Optional for parallel computation:
#   - future
#   - future.apply
# ============================================================

# ----------------------------
# 0. User settings
# Specify the sample sizes, Monte Carlo repetitions, cross-fitting settings,
# dependence parameters, and random seeds used throughout the simulation.
# ----------------------------

N_GRID <- c(300, 500, 1000) # Grid of labeled sample sizes n.
B <- 2000
K <- 20
ALPHA <- 0.05

# Strength of the residual dependence of S on X and Y, respectively.
RHO_X <- 0.4
RHO_Y <- 0.4

# Keep tau = n / N fixed, so N = round(n / tau).
TAU <- 0.15

P <- 30
COVARIANCE_BASE <- 0.5
COEFFICIENT_SEED <- 20260801
BASE_SEED <- 20260801

# Set to TRUE to use future.apply when available.
USE_PARALLEL <- FALSE
N_WORKERS <- max(1L, parallel::detectCores() - 1L)

# ----------------------------
# 1. Load PPGCM implementation
# Load the testing implementation and check the regression package
# required for the fitted GCM and PPGCM procedures.
# ----------------------------

source('PPGCM/PPGCM.R')

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Please install the R package 'glmnet'.")
}

# ----------------------------
# 2. Fixed simulation parameters
# ----------------------------

# ------------------------------------------------------------
# Generate the fixed covariance matrix and sparse coefficient vectors.
#
# The coefficient vectors are drawn once and then held fixed across all
# Monte Carlo replications and sample-size settings.
# ------------------------------------------------------------
make_simulation_parameters <- function(
    p = 30,
    covariance_base = 0.5,
    coefficient_seed = 20260801) {

  set.seed(coefficient_seed)
  
  # Construct the AR(1)-type covariance matrix for Z.
  Sigma <- outer(
    seq_len(p),
    seq_len(p),
    function(j, k) covariance_base^abs(j - k)
  )

  beta_x <- rep(0, p)
  beta_y <- rep(0, p)
  beta_s <- rep(0, p)
  
  # Generate sparse coefficient vectors with fixed supports.
  beta_x[1:5] <- stats::rnorm(5)
  beta_y[4:8] <- stats::rnorm(5)
  beta_s[4:8] <- stats::rnorm(5)

  list(
    p = p,
    Sigma = Sigma,
    chol_Sigma = chol(Sigma),
    beta_x = beta_x,
    beta_y = beta_y,
    beta_s = beta_s,
    z_cols = paste0("Z", seq_len(p))
  )
}

# ----------------------------
# 3. Data generation
# ----------------------------

# ------------------------------------------------------------
# Generate one labeled or unlabeled sample from the null model.
#
# X and Y have independent residuals given Z, so X ⟂ Y | Z.
# The auxiliary variable S shares residual components with both X and Y,
# which induces conditional dependence between S and each of them.
# ------------------------------------------------------------
generate_sample <- function(
    sample_size,
    parameters,
    rho_x = 0.4,
    rho_y = 0.4,
    labeled = TRUE) {

  p <- parameters$p
  
  # Generate correlated covariates and independent standard Gaussian errors.
  Z <- matrix(
    stats::rnorm(sample_size * p),
    nrow = sample_size,
    ncol = p
  ) %*% parameters$chol_Sigma
  colnames(Z) <- parameters$z_cols

  eps_x <- stats::rnorm(sample_size)
  eps_y <- stats::rnorm(sample_size)
  eps_s <- stats::rnorm(sample_size)
  
  # Generate X and Y under H0, then construct S so that it carries
  # residual information about both X and Y beyond Z.
  X <- drop(Z %*% parameters$beta_x) + eps_x
  Y <- drop(Z %*% parameters$beta_y) + eps_y

  S <- drop(Z %*% parameters$beta_s) +
    rho_x * eps_x +
    rho_y * eps_y +
    eps_s

  if (labeled) {
    data.frame(
      X = X,
      Y = Y,
      S = S,
      Z,
      check.names = FALSE
    )
  } else {
    # eps_y is generated because S depends on eps_y, but Y is masked.
    data.frame(
      X = X,
      S = S,
      Z,
      check.names = FALSE
    )
  }
}

# ----------------------------
# 4. Generic oracle score test
# ----------------------------

# ------------------------------------------------------------
# Generic two-sided oracle score test.
#
# Given a vector of oracle scores, compute its sample mean and variance
# and form the corresponding asymptotic normal test.
# ------------------------------------------------------------
oracle_score_test <- function(score, alpha = 0.05, eps = 1e-10) {
  n <- length(score)

  T_hat <- mean(score)
  variance_hat <- mean((score - T_hat)^2)

  if (!is.finite(variance_hat) || variance_hat <= eps) {
    variance_hat <- eps
  }

  z_score <- sqrt(n) * T_hat / sqrt(variance_hat)
  p_value <- 2 * (1 - stats::pnorm(abs(z_score)))
  reject <- as.numeric(p_value < alpha)

  list(
    T_hat = T_hat,
    variance_hat = variance_hat,
    z_score = z_score,
    p_value = p_value,
    reject = reject
  )
}

# ----------------------------
# 5. Oracle GCM
# ----------------------------

# ------------------------------------------------------------
# Oracle GCM.
#
# Use the true conditional means
#   h(Z) = E[X | Z] = Z^T beta_x
#   g(Z) = E[Y | Z] = Z^T beta_y
# to construct the oracle GCM score.
# ------------------------------------------------------------
oracle_gcm_test <- function(
    label_data,
    parameters,
    alpha = 0.05) {

  Z <- as.matrix(label_data[, parameters$z_cols, drop = FALSE])

  h_true <- drop(Z %*% parameters$beta_x)
  g_true <- drop(Z %*% parameters$beta_y)
  
  # Under H0 this score has population mean zero.
  score <- (label_data$X - h_true) *
    (label_data$Y - g_true)

  oracle_score_test(score = score, alpha = alpha)
}

# ------------------------------------------------------------
# 6. Oracle SGCM.
#
# SGCM keeps h(Z) = E[X | Z] but replaces g(Z) by
#   q(Z,S) = E[Y | Z,S].
# Under the present model, this produces a nonzero population score
# even though X ⟂ Y | Z, illustrating the type-I error failure of SGCM.
# ------------------------------------------------------------

oracle_sgcm_test <- function(
    label_data,
    parameters,
    rho_x = 0.4,
    rho_y = 0.4,
    alpha = 0.05) {

  Z <- as.matrix(label_data[, parameters$z_cols, drop = FALSE])

  h_true <- drop(Z %*% parameters$beta_x)
  
  # Since S - Z^T beta_s = rho_x eps_x + rho_y eps_y + eps_s,
  # Gaussian conditioning gives the exact conditional mean E[Y | Z,S].
  D <- 1 + rho_x^2 + rho_y^2
  q_true <- drop(Z %*% parameters$beta_y) +
    rho_y / D * (
      label_data$S - drop(Z %*% parameters$beta_s)
    )
  
  # The resulting SGCM score is generally not centered at zero under H0.
  score <- (label_data$X - h_true) *
    (label_data$Y - q_true)

  oracle_score_test(score = score, alpha = alpha)
}

# ------------------------------------------------------------
# 7. One Monte Carlo replication
#
# Generate independent labeled and unlabeled samples, compute the two
# oracle tests, and then apply the fitted cross-fitted GCM and PPGCM.
# ------------------------------------------------------------
run_one_replication <- function(
    n,
    N,
    parameters,
    rho_x = 0.4,
    rho_y = 0.4,
    K = 20,
    alpha = 0.05,
    replication_seed = 1) {

  set.seed(replication_seed)
  
  # Generate independent labeled and unlabeled samples from the same population.
  label_data <- generate_sample(
    sample_size = n,
    parameters = parameters,
    rho_x = rho_x,
    rho_y = rho_y,
    labeled = TRUE
  )

  unlabel_data <- generate_sample(
    sample_size = N,
    parameters = parameters,
    rho_x = rho_x,
    rho_y = rho_y,
    labeled = FALSE
  )

  # Oracle GCM
  oracle_gcm <- oracle_gcm_test(
    label_data = label_data,
    parameters = parameters,
    alpha = alpha
  )

  # Oracle SGCM
  oracle_sgcm <- oracle_sgcm_test(
    label_data = label_data,
    parameters = parameters,
    rho_x = rho_x,
    rho_y = rho_y,
    alpha = alpha
  )

  # Cross-fitted GCM and optimal-weight PPGCM.
  # h, g, and f are all estimated by lasso.
  # PPGCM uses S only through the prediction model f(X,Z,S), not in the null hypothesis.
  fitted <- ppicit_test(
    label_data = label_data,
    unlabel_data = unlabel_data,
    K = K,
    fit_h = "lasso",
    fit_g = "lasso",
    fit_f = "lasso",
    x_col = "X",
    y_col = "Y",
    s_col = "S",
    z_cols = parameters$z_cols,
    alpha = alpha,
    seed = replication_seed + 1000000L,
    standardize_lasso = TRUE
  )

  fitted_gcm <- fitted$gcm
  fitted_ppgcm <- fitted$ppi
  
  # Store the rejection indicator, test statistic, variance estimate,
  # and estimated PPGCM weight for this replication.
  data.frame(
    n = rep(n, 4),
    N = rep(N, 4),
    method = c(
      "Oracle GCM",
      "GCM",
      "Oracle SGCM",
      "PPGCM (opt)"
    ),
    reject = c(
      oracle_gcm$reject,
      as.numeric(fitted_gcm$reject),
      oracle_sgcm$reject,
      as.numeric(fitted_ppgcm$reject)
    ),
    TS = c(
      oracle_gcm$T_hat,
      fitted_gcm$T_hat,
      oracle_sgcm$T_hat,
      fitted_ppgcm$T_hat
    ),
    AV = c(
      oracle_gcm$variance_hat,
      fitted_gcm$variance_hat,
      oracle_sgcm$variance_hat,
      fitted_ppgcm$variance_hat
    ),
    weight = c(
      NA_real_,
      NA_real_,
      NA_real_,
      fitted_ppgcm$weight_hat_mean
    )
  )
}

# ----------------------------
# 8. Run one sample-size setting
# ----------------------------

# Return an NA-filled result when a replication fails, so that the
# remaining Monte Carlo replications can still be summarized.
failed_replication <- function(n, N) {
  data.frame(
    n = rep(n, 4),
    N = rep(N, 4),
    method = c(
      "Oracle GCM",
      "GCM",
      "Oracle SGCM",
      "PPGCM (opt)"
    ),
    reject = NA_real_,
    TS = NA_real_,
    AV = NA_real_,
    weight = NA_real_
  )
}

# ------------------------------------------------------------
# Run all Monte Carlo replications for one labeled sample size n.
#
# The unlabeled sample size is chosen as N = round(n / tau), so that
# the labeled-to-unlabeled ratio remains approximately fixed.
# ------------------------------------------------------------
run_one_sample_size <- function(
    n,
    B,
    parameters,
    tau,
    rho_x,
    rho_y,
    K,
    alpha,
    base_seed = 20260801,
    use_parallel = FALSE) {
  
  # Keep the ratio n/N approximately equal to the target tau.
  N <- round(n / tau)
  n_index <- match(n, N_GRID)
  
  # Wrap one replication with error handling and a reproducible seed.
  one_job <- function(b) {
    tryCatch(
      run_one_replication(
        n = n,
        N = N,
        parameters = parameters,
        rho_x = rho_x,
        rho_y = rho_y,
        K = K,
        alpha = alpha,
        replication_seed = base_seed + 100000L * n_index + b
      ),
      error = function(e) {
        message(
          "Replication ", b,
          " failed for n = ", n,
          ": ", conditionMessage(e)
        )
        failed_replication(n = n, N = N)
      }
    )
  }
  
  # Run replications in parallel when future.apply is available;
  # otherwise fall back to sequential evaluation.
  if (
    use_parallel &&
      requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)
  ) {
    out <- future.apply::future_lapply(
      seq_len(B),
      one_job,
      future.seed = TRUE
    )
  } else {
    out <- lapply(seq_len(B), one_job)
  }

  do.call(rbind, out)
}

# ------------------------------------------------------------
# 9. Monte Carlo summary.
#
# For each method, compute the rejection rate and the mean and standard
# deviation of the test statistic, variance estimate, and PPGCM weight.
# ------------------------------------------------------------
summarize_simulation <- function(raw_results) {

  method_order <- c(
    "Oracle GCM",
    "GCM",
    "Oracle SGCM",
    "PPGCM (opt)"
  )

  raw_results$method <- factor(
    raw_results$method,
    levels = method_order
  )

  split_results <- split(
    raw_results,
    list(raw_results$n, raw_results$method),
    drop = TRUE
  )

  summary_list <- lapply(split_results, function(dat) {
    valid <- stats::complete.cases(dat[, c("reject", "TS", "AV")])
    dat <- dat[valid, , drop = FALSE]

    finite_weights <- dat$weight[is.finite(dat$weight)]

    data.frame(
      n = dat$n[1],
      N = dat$N[1],
      method = as.character(dat$method[1]),
      replications = nrow(dat),
      RR = mean(dat$reject),
      TS_mean = mean(dat$TS),
      TS_sd = stats::sd(dat$TS),
      AV_mean = mean(dat$AV),
      AV_sd = stats::sd(dat$AV),
      weight_mean = if (length(finite_weights) > 0) {
        mean(finite_weights)
      } else {
        NA_real_
      },
      weight_sd = if (length(finite_weights) > 1) {
        stats::sd(finite_weights)
      } else {
        NA_real_
      }
    )
  })

  out <- do.call(rbind, summary_list)

  out$method <- factor(out$method, levels = method_order)
  out <- out[order(out$n, out$method), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------
# Format the Monte Carlo summary for reporting.
#
# Test statistics, asymptotic variances, and weights are shown as
# mean(standard deviation), while rejection rates are reported separately.
# ------------------------------------------------------------
format_simulation_table <- function(summary_results) {
  data.frame(
    n = summary_results$n,
    N = summary_results$N,
    method = as.character(summary_results$method),
    RR = sprintf("%.4f", summary_results$RR),
    TS = sprintf(
      "%.3f(%.3f)",
      summary_results$TS_mean,
      summary_results$TS_sd
    ),
    AV = sprintf(
      "%.3f(%.3f)",
      summary_results$AV_mean,
      summary_results$AV_sd
    ),
    weight = ifelse(
      is.finite(summary_results$weight_mean),
      sprintf(
        "%.3f(%.3f)",
        summary_results$weight_mean,
        summary_results$weight_sd
      ),
      "--"
    ),
    check.names = FALSE
  )
}

# ----------------------------
# 10. Run the simulation
# Generate the fixed model parameters and run the complete simulation
# over all labeled sample sizes in the grid.
# ----------------------------

parameters <- make_simulation_parameters(
  p = P,
  covariance_base = COVARIANCE_BASE,
  coefficient_seed = COEFFICIENT_SEED
)

if (
  USE_PARALLEL &&
    requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)
) {
  # Initialize the optional multisession parallel backend.
  future::plan(
    future::multisession,
    workers = N_WORKERS
  )
}

# Run all Monte Carlo replications for each sample-size setting.
raw_results <- do.call(
  rbind,
  lapply(
    N_GRID,
    run_one_sample_size,
    B = B,
    parameters = parameters,
    tau = TAU,
    rho_x = RHO_X,
    rho_y = RHO_Y,
    K = K,
    alpha = ALPHA,
    base_seed = BASE_SEED,
    use_parallel = USE_PARALLEL
  )
)

if (
  USE_PARALLEL &&
    requireNamespace("future", quietly = TRUE)
) {
  # Restore the default sequential execution plan after the simulation.
  future::plan(future::sequential)
}

# Aggregate the replication-level results and print the formatted table.
summary_results <- summarize_simulation(raw_results)
report_table <- format_simulation_table(summary_results)

print(report_table, row.names = FALSE)


