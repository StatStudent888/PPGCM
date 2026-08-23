# ============================================================
# Main Monte Carlo simulation for PPGCM type-I error and power
#
# Purpose:
#   This script runs the main simulation study in Section 4.2.1.
#   For a fixed signal level c, it repeatedly generates labeled and
#   unlabeled samples, applies GCM, PPGCM (w=1), and PPGCM (opt),
#   and summarizes their rejection rates, test statistics, estimated
#   asymptotic variances, p-values, and estimated optimal weights.
#
# Simulation designs:
#   Both the linear and nonlinear data-generating mechanisms from
#   Section 4.2.1 are supported. The nonlinear design uses generalized
#   random forests, whereas the linear design uses the lasso.
#
# Main inputs:
#   n_rep      : number of Monte Carlo replications
#   c_coef     : signal coefficient c; c = 0 evaluates type-I error,
#                while c > 0 evaluates power
#   n_label    : labeled sample size n
#   n_unlabel  : unlabeled sample size N
#   p_z        : dimension of Z
#   gbar_type  : "linear" or "nonlinear" simulation design
#   K          : number of cross-fitting folds
#   fit_h, fit_g, fit_f : regression methods for the nuisance and
#                         prediction models
#
# Main output:
#   run_ppicit_replications_one_c() returns Monte Carlo summaries,
#   replication-level results, and the complete simulation configuration.
#
# The coefficient vectors are generated once and held fixed across
# Monte Carlo replications, as specified in Section 4.2.1.
# ============================================================

# Set the project directory and load the PPGCM implementation
# together with the simulation data-generating functions.
setwd("PPGCM")

source("PPGCM.R")
source("simulation/type-I_power/generate_data.R")

# ------------------------------------------------------------
# Run Monte Carlo replications for one signal level c.
#
# The coefficient vectors are generated once before the replication
# loop and shared by all simulated datasets. Each replication then
# generates a new labeled/unlabeled sample and applies the three tests.
# Replications can be run sequentially or in parallel.
# ------------------------------------------------------------
run_ppicit_replications_one_c <- function(n_rep = 500,
                                          c_coef = 0,
                                          n_label = 300,
                                          n_unlabel = 2000,
                                          p_z = 30,
                                          gbar_type = c("linear", "nonlinear"),
                                          rho = 0.5,
                                          sx_gamma = 0,
                                          sx_power = 1,
                                          K = 5,
                                          fit_h = "lasso",
                                          fit_g = "lasso",
                                          fit_f = "lasso",
                                          alpha = 0.05,
                                          seed = 123,
                                          beta_start = 1,
                                          beta_nnz = 5,
                                          beta_seed = NULL,
                                          show_progress = TRUE,
                                          use_progress_bar = TRUE,
                                          parallel = TRUE,
                                          n_cores = max(1, parallel::detectCores() - 1)) {
  gbar_type <- match.arg(gbar_type)

  # ============================================================
  # Generate fixed coefficients shared by all replications
  # Following Section 4.2.1, beta_x is supported on {1,...,5},
  # while beta_s and beta_g (beta_y in the paper) are supported on
  # {4,...,8}. Their nonzero entries are drawn once from N(0,1).
  # ============================================================
  z_names <- paste0("Z", seq_len(p_z))
  
  # Construct a sparse coefficient vector with Gaussian nonzero entries.
  make_sparse_beta <- function(p_z, beta_start, beta_nnz) {
    if (beta_start < 1) {
      stop("beta_start must be at least 1.")
    }
    if (beta_start + beta_nnz - 1 > p_z) {
      stop("beta_start + beta_nnz - 1 cannot exceed p_z.")
    }

    beta <- rep(0, p_z)
    idx <- beta_start:(beta_start + beta_nnz - 1)
    beta[idx] <- stats::rnorm(beta_nnz)
    names(beta) <- z_names
    beta
  }

  if (is.null(beta_seed)) beta_seed <- seed

  set.seed(beta_seed)
  
  # Generate the fixed coefficient vectors used throughout the simulation.
  # beta_g corresponds to beta_y in the notation of the paper.
  beta_x <- make_sparse_beta(p_z, beta_start, beta_nnz)
  beta_s <- make_sparse_beta(p_z, beta_start + 3, beta_nnz)
  beta_g <- make_sparse_beta(p_z, beta_start + 3, beta_nnz)

  # ============================================================
  # Run one Monte Carlo replication
  # ============================================================
  # Generate a new dataset, apply GCM and both PPGCM variants,
  # and retain the quantities needed for the final Monte Carlo summary.
  # In parallel mode, this function is sent independently to each worker.
  one_replication <- function(r) {
    tryCatch({
      # Generate independent labeled and unlabeled samples from the selected DGP.
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
        seed = seed + r
      )
      
      # Apply the cross-fitted GCM, PPGCM (w=1), and PPGCM (opt).
      res <- ppicit_test(
        label_data = dat$label_data,
        unlabel_data = dat$unlabel_data,
        K = K,
        fit_h = fit_h,
        fit_g = fit_g,
        fit_f = fit_f,
        z_cols = dat$truth$z_cols,
        alpha = alpha,
        seed = seed + r
      )
      
      # Retain the inferential quantities reported in the simulation tables.
      list(
        error = FALSE,
        error_msg = NA_character_,

        ppi_p_value = res$ppi$p_value,
        ppi_reject = res$ppi$reject,
        ppi_T_hat = res$ppi$T_hat,
        ppi_variance_hat = res$ppi$variance_hat,
        ppi_weight_hat = res$ppi$weight_hat_mean,

        ppi_w1_p_value = res$ppi_w1$p_value,
        ppi_w1_reject = res$ppi_w1$reject,
        ppi_w1_T_hat = res$ppi_w1$T_hat,
        ppi_w1_variance_hat = res$ppi_w1$variance_hat,
        ppi_w1_weight_hat = res$ppi_w1$weight_hat_mean,

        gcm_p_value = res$gcm$p_value,
        gcm_reject = res$gcm$reject,
        gcm_T_hat = res$gcm$T_hat,
        gcm_variance_hat = res$gcm$variance_hat
      )
    }, error = function(e) {
      list(
        error = TRUE,
        error_msg = paste0("Replication ", r, " failed: ", e$message),

        ppi_p_value = NA_real_,
        ppi_reject = NA,
        ppi_T_hat = NA_real_,
        ppi_variance_hat = NA_real_,
        ppi_weight_hat = NA_real_,

        ppi_w1_p_value = NA_real_,
        ppi_w1_reject = NA,
        ppi_w1_T_hat = NA_real_,
        ppi_w1_variance_hat = NA_real_,
        ppi_w1_weight_hat = NA_real_,

        gcm_p_value = NA_real_,
        gcm_reject = NA,
        gcm_T_hat = NA_real_,
        gcm_variance_hat = NA_real_
      )
    })
  }

  # ============================================================
  # Run replications: sequential or parallel
  # Use PSOCK parallelization when requested; otherwise run the same
  # replication function sequentially.
  # ============================================================
  r_grid <- seq_len(n_rep)

  if (parallel && n_cores > 1) {
    n_cores <- min(n_cores, n_rep)

    if (show_progress) {
      cat("Running", n_rep, "replications on", n_cores, "cores ...\n")
      cat("Note: base R parallel does not show a real-time progress bar on Windows.\n")
    }
    
    # Create a worker cluster and export the functions required by each replication.
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export functions defined in sourced files to PSOCK workers.
    parallel::clusterExport(
      cl,
      varlist = c("simulate_ppicit_plm", "ppicit_test", ".fit_regression"),
      envir = .GlobalEnv
    )

    # Load package namespaces on workers only when needed.
    # requireNamespace() is enough because the functions use glmnet:: and grf::.
    if (fit_h == "lasso" || fit_g == "lasso" || fit_f == "lasso") {
      parallel::clusterEvalQ(cl, {
        if (!requireNamespace("glmnet", quietly = TRUE)) {
          stop("Package 'glmnet' is required for method = 'lasso'.")
        }
        NULL
      })
    }
    if (fit_h == "grf" || fit_g == "grf" || fit_f == "grf") {
      parallel::clusterEvalQ(cl, {
        if (!requireNamespace("grf", quietly = TRUE)) {
          stop("Package 'grf' is required for method = 'grf'.")
        }
        NULL
      })
    }

    # Load balancing is useful because some replications may take longer.
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

  # Print error messages, if any
  error_msg <- vapply(rep_list, function(x) x$error_msg, character(1))
  error_msg <- error_msg[!is.na(error_msg)]
  if (length(error_msg) > 0) {
    message(paste(error_msg, collapse = "\n"))
  }

  # ============================================================
  # Collect results
  # Convert the replication-level lists into vectors for Monte Carlo summaries.
  # ============================================================
  get_col <- function(name, mode = c("numeric", "logical")) {
    mode <- match.arg(mode)
    if (mode == "numeric") {
      as.numeric(vapply(rep_list, function(x) x[[name]], numeric(1)))
    } else {
      as.logical(vapply(rep_list, function(x) x[[name]], logical(1)))
    }
  }

  has_error <- get_col("error", "logical")

  ppi_pvals <- get_col("ppi_p_value", "numeric")
  ppi_reject <- get_col("ppi_reject", "logical")
  ppi_T <- get_col("ppi_T_hat", "numeric")
  ppi_var <- get_col("ppi_variance_hat", "numeric")
  ppi_weight <- get_col("ppi_weight_hat", "numeric")

  ppi_w1_pvals <- get_col("ppi_w1_p_value", "numeric")
  ppi_w1_reject <- get_col("ppi_w1_reject", "logical")
  ppi_w1_T <- get_col("ppi_w1_T_hat", "numeric")
  ppi_w1_var <- get_col("ppi_w1_variance_hat", "numeric")
  ppi_w1_weight <- get_col("ppi_w1_weight_hat", "numeric")

  gcm_pvals <- get_col("gcm_p_value", "numeric")
  gcm_reject <- get_col("gcm_reject", "logical")
  gcm_T <- get_col("gcm_T_hat", "numeric")
  gcm_var <- get_col("gcm_variance_hat", "numeric")
  
  # Restrict summary calculations to replications completed successfully
  # for all three testing procedures.
  valid <- !has_error &
    !is.na(ppi_reject) &
    !is.na(ppi_w1_reject) &
    !is.na(gcm_reject)

  sd_or_na <- function(x) {
    if (sum(!is.na(x)) <= 1) return(NA_real_)
    stats::sd(x, na.rm = TRUE)
  }
  
  # Summarize rejection rates and the distributions of the main estimated
  # quantities across successful Monte Carlo replications.
  summary <- data.frame(
    method = c("PPI-GCM", "PPI-GCM (w=1)", "GCM"),
    c_coef = c(c_coef, c_coef, c_coef),
    n_label = c(n_label, n_label, n_label),
    n_unlabel_used = c(n_unlabel, n_unlabel, 0),
    n_rep = c(n_rep, n_rep, n_rep),
    n_success = c(sum(valid), sum(valid), sum(valid)),
    n_error = c(sum(has_error), sum(has_error), sum(has_error)),

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
    )
  )
  
  # Return both the compact Monte Carlo summary and the full
  # replication-level results for subsequent tables or diagnostics.
  out <- list(
    summary = summary,

    rejection_rate = summary[, c("method", "rejection_rate")],

    results = data.frame(
      replication = seq_len(n_rep),

      ppi_p_value = ppi_pvals,
      ppi_reject = ppi_reject,
      ppi_T_hat = ppi_T,
      ppi_variance_hat = ppi_var,
      ppi_weight_hat = ppi_weight,

      ppi_w1_p_value = ppi_w1_pvals,
      ppi_w1_reject = ppi_w1_reject,
      ppi_w1_T_hat = ppi_w1_T,
      ppi_w1_variance_hat = ppi_w1_var,
      ppi_w1_weight_hat = ppi_w1_weight,

      gcm_p_value = gcm_pvals,
      gcm_reject = gcm_reject,
      gcm_T_hat = gcm_T,
      gcm_variance_hat = gcm_var,

      error = has_error
    ),

    config = list(
      n_rep = n_rep,
      c_coef = c_coef,
      n_label = n_label,
      n_unlabel = n_unlabel,
      p_z = p_z,
      gbar_type = gbar_type,
      rho = rho,
      sx_gamma = sx_gamma,
      sx_power = sx_power,
      K = K,
      fit_h = fit_h,
      fit_g = fit_g,
      fit_f = fit_f,
      alpha = alpha,
      seed = seed,
      beta_x = beta_x,
      beta_s = beta_s,
      beta_g = beta_g,
      parallel = parallel,
      n_cores = ifelse(parallel, n_cores, 1)
    )
  )

  class(out) <- "ppicit_sim_one_c"
  return(out)
}

# ============================================================
# Main simulation settings used in Section 4.2.1
# ============================================================

# Nonlinear design: use generalized random forests for h, g, and f.
# The setting below corresponds to one combination of
# (n, N, c) used in the main simulation study.
sim_res <- run_ppicit_replications_one_c(
  n_rep = 2000,
  c_coef = 0.2, #0,0.1,0.15,0.2
  n_label = 200, #200,300
  n_unlabel = 2000, #2000,10000
  p_z = 30,
  gbar_type = "nonlinear",
  sx_gamma = 0,
  sx_power = 1,
  K = 20, # Main setting; K = 10 and 30 are used in the sensitivity analysis.
  fit_h = "grf",
  fit_g = "grf",
  fit_f = "grf",
  alpha = 0.05,
  seed = 123,
  beta_start = 1,
  beta_nnz = 5,
  beta_seed = 999,
  parallel = TRUE,
  n_cores = 20
)

# ------------------------------------------------------------
# Alternative: linear design
# ------------------------------------------------------------
# Use the lasso for all three fitted regression functions, as in
# Section 4.2.1. Uncomment this block to run the linear setting.

# sim_res <- run_ppicit_replications_one_c(
#   n_rep = 2000,
#   c_coef = 0.2, #0,0.1,0.15,0.2
#   n_label = 300, #200,300
#   n_unlabel = 2000, #2000,10000
#   p_z = 30,
#   gbar_type = "linear",
#   sx_gamma = 0,
#   sx_power = 1,
#   K = 20, # Main setting; K = 10 and 30 are used in the sensitivity analysis.
#   fit_h = "lasso",
#   fit_g = "lasso",
#   fit_f = "lasso",
#   alpha = 0.05,
#   seed = 123,
#   beta_start = 1,
#   beta_nnz = 5,
#   beta_seed = 999,
#   parallel = TRUE,
#   n_cores = 20
# )

# Inspect the simulation configuration and the resulting Monte Carlo summary.
sim_res$config
sim_res$summary


