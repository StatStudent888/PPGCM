# ============================================================
# Main script for the fixed-t hardness simulation
#
# Purpose:
#   This script runs the oracle fixed-t simulation for PPGCM and GCM
#   across a range of labeled sample sizes, and then produces plots of
#   empirical power and estimated asymptotic variance.
#
# Simulation settings:
#   Two prediction-model settings are available:
#     1. correctly specified prediction model;
#     2. deliberately misspecified prediction model.
#   Activate the desired setting by sourcing the corresponding simulation
#   file and using the matching run_hardness_*_by_n() function below.
#
# Main inputs:
#   n_label_values : labeled sample sizes considered in the simulation
#   fixed_t        : local-alternative parameter, with c_n = t / sqrt(n)
#   n_over_N       : target labeled-to-unlabeled sample-size ratio n/N
#   n_rep          : number of Monte Carlo replications
#   p_z            : dimension of Z
#   rho            : correlation parameter in the covariance matrix of Z
#
# Main outputs:
#   res_by_n : simulation results across labeled sample sizes
#   p1       : empirical power versus labeled sample size
#   p2       : estimated asymptotic variance versus labeled sample size
#
# The simulation and plotting functions are defined in separate R scripts
# under simulation/hardness/.
# ============================================================

# Set the project working directory so that the relative source paths below work.
setwd('PPGCM')

# ============================================================
# Choose the prediction-model setting
# ============================================================
# By default, run the correctly specified prediction-model simulation.
# To run the misspecified setting instead, comment out the first line
# and uncomment the second line together with the corresponding function below.

source("simulation/hardness/ppgcm_hardness_correctly_specified.R")
# source("simulation/hardness/ppgcm_hardness_misspecified.R")

# ============================================================
# Run the fixed-t simulation: correctly specified prediction model
# ============================================================
# For each labeled sample size n, the signal is set to
#   c_n = fixed_t / sqrt(n),
# while N is chosen so that n/N remains approximately fixed.
res_by_n <- run_hardness_oracle_noS_by_n(
  n_label_values = c(100, 300, 500, 1000, 1500, 2000, 3000),
  fixed_t = 1.5,
  n_over_N = 0.15,   # Keep n/N fixed; for example, n = 300 gives N = 2000.
  n_rep = 2000,
  p_z = 30,
  rho = 0.5,
  alpha = 0.05,
  seed = 123,
  beta_start = 1,
  beta_nnz = 5,
  beta_seed = 999,
  show_progress = TRUE
)

# ============================================================
# Alternative: misspecified prediction-model simulation
# ============================================================
# Uncomment this block, together with the misspecified source file above,
# to run the same fixed-t experiment under the misspecified prediction model.

# res_by_n <- run_hardness_oracle_noS_fminus_by_n(
#   n_label_values = c(100, 300, 500, 1000, 1500, 2000, 3000),
#   fixed_t = 1.5,
#   n_over_N = 0.15,   # Keep n/N fixed; for example, n = 300 gives N = 2000.
#   n_rep = 2000,
#   p_z = 30,
#   rho = 0.5,
#   alpha = 0.05,
#   seed = 123,
#   beta_start = 1,
#   beta_nnz = 5,
#   beta_seed = 999,
#   show_progress = TRUE
# )

# ============================================================
# Load plotting functions
# ============================================================
source("simulation/hardness/ppigcm_hardness_plot_functions.R")

# Plot empirical power against labeled sample size.
# The y-axis range and legend position are chosen for this simulation setting.
p1 <- plot_rejection_by_n_pretty(res_by_n,y_limits=c(0.275,0.375),y_breaks=seq(0.275, 0.375, by = 0.025),legend_position = c(0.98, 0.98))
print(p1)

# Plot the mean estimated asymptotic variance against labeled sample size.
# Error bars represent one standard deviation across Monte Carlo replications.
p2 <- plot_variance_by_n_pretty(
  res_by_n,
  y_limits=c(0.6,1.4),y_breaks=seq(0.6, 1.4, by = 0.1),legend_position = c(0.98, 0.98),
  errorbar_type = "sd"
)
print(p2)

# ============================================================
# Alternative plotting ranges for the misspecified setting
# ============================================================

# p1 <- plot_rejection_by_n_pretty(res_by_n,y_limits=c(0.08,0.34),y_breaks=seq(0.08, 0.34, by = 0.04),legend_position = c(0.98, 0.68))
# print(p1)
# 
# p2 <- plot_variance_by_n_pretty(
#   res_by_n,
#   y_limits=c(0,8),y_breaks=seq(0, 8, by = 2),legend_position = c(0.98, 0.58),
#   errorbar_type = "sd"
# )
# 
# print(p2)

