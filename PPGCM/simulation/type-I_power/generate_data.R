# ============================================================
# Data generation for the main PPGCM simulation
#
# Purpose:
#   This script generates labeled and unlabeled samples for the
#   linear and nonlinear simulation designs in Section 4.2.1.
#   The generated data are used to compare the cross-fitted GCM,
#   PPGCM, and competing conditional-independence tests.
#
# Data structure:
#   Labeled observations contain (Y, X, S, Z), whereas unlabeled
#   observations contain only (X, S, Z). Both samples are generated
#   independently from the same population.
#
# Common setup:
#   Z ~ N(0, Sigma), with Sigma_jk = rho^{|j-k|},
#   and eps_x, eps_s, eps_y are mutually independent N(0,1).
#
# Available designs:
#   gbar_type = "linear"    : linear model in Section 4.2.1
#   gbar_type = "nonlinear" : nonlinear model in Section 4.2.1
#
# Main inputs:
#   n_label   : labeled sample size
#   n_unlabel : unlabeled sample size
#   p_z       : dimension of Z
#   c_coef    : signal coefficient c
#   rho       : correlation parameter for Z
#   beta_x, beta_s, beta_g : coefficient vectors used in the linear model
#
# Main output:
#   A list containing the labeled data, unlabeled data, the full
#   simulated unlabeled sample, and the data-generating parameters.
# ============================================================

# ------------------------------------------------------------
# Generate one labeled/unlabeled dataset for the main simulation.
#
# The function supports both the linear and nonlinear designs.
# The same data-generating mechanism is used for labeled and unlabeled
# observations, with Y removed from the latter before returning.
# ------------------------------------------------------------
simulate_ppicit_plm <- function(n_label,
                                n_unlabel = 0,
                                p_z = 30,
                                gbar_type = c("linear", "nonlinear"),
                                c_coef = 0,
                                rho = 0.5,
                                sx_gamma = 0,
                                sx_power = 1,
                                beta_x = NULL,
                                beta_s = NULL,
                                beta_g = NULL,
                                seed = NULL) {
  gbar_type <- match.arg(gbar_type)
  
  if (!is.null(seed)) set.seed(seed)
  if (p_z < 1) stop("p_z must be at least 1.")
  if (n_label < 1) stop("n_label must be positive.")
  if (n_unlabel < 0) stop("n_unlabel cannot be negative.")
  
  z_names <- paste0("Z", seq_len(p_z))
  
  # Construct the AR(1)-type covariance matrix of Z.
  Sigma <- outer(seq_len(p_z), seq_len(p_z), function(i, j) rho^abs(i - j))
  chol_Sigma <- chol(Sigma)
  
  # Supply default coefficient vectors when none are provided by the caller.
  default_coef <- function(p, phase = 0) {
    v <- sin(seq_len(p) + phase) / (seq_len(p)^0.7)
    v / sqrt(sum(v^2))
  }
  
  if (is.null(beta_x)) beta_x <- default_coef(p_z, phase = 0.0)
  if (is.null(beta_s)) beta_s <- default_coef(p_z, phase = 1.3)
  if (is.null(beta_g)) beta_g <- default_coef(p_z, phase = 2.1)
  
  beta_x <- as.numeric(beta_x)
  beta_s <- as.numeric(beta_s)
  beta_g <- as.numeric(beta_g)
  
  if (length(beta_x) != p_z) stop("length(beta_x) must equal p_z.")
  if (length(beta_s) != p_z) stop("length(beta_s) must equal p_z.")
  if (length(beta_g) != p_z) stop("length(beta_g) must equal p_z.")
  
  sx_lambda <- sx_gamma / (n_label^sx_power)
  
  # Generate one sample of size n from the selected data-generating mechanism.
  .generate_one <- function(n) {
    # Generate the correlated covariates and independent Gaussian errors.
    Z <- matrix(stats::rnorm(n * p_z), nrow = n, ncol = p_z) %*% chol_Sigma
    colnames(Z) <- z_names
    
    eps_x <- stats::rnorm(n)
    eps_s <- stats::rnorm(n)
    eps_y <- stats::rnorm(n)
    
    
    # Linear design:
    #   X = Z^T beta_x + eps_x,
    #   S = Z^T beta_s + eps_s,
    #   Y = c X + 0.7 S + Z^T beta_g + eps_y.
    if (gbar_type == "linear") {
      X <- as.numeric(Z %*% beta_x + eps_x)
      S <- as.numeric(Z %*% beta_s + sx_lambda * eps_x + eps_s)
      # S <- as.numeric(Z %*% beta_s + eps_s)
      gbar <- 0.7 * S + as.numeric(Z %*% beta_g)
      Y <- c_coef * X + gbar + eps_y
    } else {
      # Nonlinear design based on thresholded functions and interactions of Z.
      # Indicator variables entering the nonlinear regression functions.
      I1 <- I(Z[,1] > 0)
      I2 <- I(Z[,2] > 0.5)
      I3 <- I(Z[,3] > -0.5)
      I4 <- I(abs(Z[,4]) > 1)

      I11 <- I(Z[,11] > 0)
      I12 <- I(Z[,12] > 0)
      I21 <- I(Z[,21] > 0)
      I22 <- I(Z[,22] > 0)
      # Conditional mean of X given Z.
      mean_x <- 0.5 * (-I1 + 0.8 * I2 + I3 - 0.8 * I4 - 0.8 * (I1 * I4 + I2 * I3)) +
        0.1 * (I11 + I12 - I21 * I22)
      X <- as.numeric(mean_x + eps_x)
      # Conditional mean of S given Z.
      mean_s <- 0.5 * (-I2 + 0.8 * I1 + I4 - 0.8 * I3 - 0.8 * (I1 * I2 + I3 * I4)) -
        0.1 * (I21 + I22 - I11 * I12 )
      S <- as.numeric(mean_s + sx_lambda * eps_x + eps_s)
      # Conditional mean of Y given X, S, and Z.
      mean_y <- c_coef*X + 0.6 * (S + 0.8*sin(S)) + 0.5 * (-0.8 * I1 + I2 + I3 - 0.8 * I4 + (I1 * I2 + I3 * I4)) +
        0.1 * (I11 + I12 + 2 * I21 - I11 * I12 + I21 * I22)
      Y <- as.numeric(mean_y + eps_y)
    }
    
    
    # Assemble the full simulated observation with Y included.
    data.frame(
      Y = Y,
      X = X,
      S = S,
      Z,
      check.names = FALSE
    )
  }
  
  # Generate the labeled sample.
  label_data <- .generate_one(n_label)
  
  # Generate an independent unlabeled sample and remove Y before returning it.
  if (n_unlabel > 0) {
    unlabel_full <- .generate_one(n_unlabel)
    
    # Direct input for ppicit_test(): no Y column
    unlabel_data <- unlabel_full[, c("X", "S", z_names), drop = FALSE]
  } else {
    unlabel_full <- NULL
    unlabel_data <- NULL
  }
  
  # Return the simulated samples together with the data-generating parameters.
  list(
    label_data = label_data,
    unlabel_data = unlabel_data,
    unlabel_full = unlabel_full,
    truth = list(
      c_coef = c_coef,
      gbar_type = gbar_type,
      rho = rho,
      sx_gamma = sx_gamma,
      sx_power = sx_power,
      sx_lambda = sx_lambda,
      beta_x = beta_x,
      beta_s = beta_s,
      beta_g = beta_g,
      z_cols = z_names
    )
  )
}