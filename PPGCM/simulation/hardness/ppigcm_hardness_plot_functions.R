# ============================================================
# Plotting functions for the fixed-t hardness simulation
#
# Purpose:
#   This script provides plotting and saving functions for summarizing
#   simulation results across different labeled sample sizes n.
#   It produces:
#     1. empirical rejection rate / power versus n;
#     2. estimated asymptotic variance versus n.
#
# Main input:
#   res_by_n : simulation results containing a summary data frame with
#              one row for each combination of labeled sample size and
#              testing method.
#
# Expected variables in res_by_n$summary include:
#   n_label            : labeled sample size
#   method             : testing method
#   rejection_rate     : empirical rejection rate
#   mean_variance_hat  : mean estimated asymptotic variance
#   sd_variance_hat    : standard deviation of the variance estimates
#   n_rep              : number of simulation repetitions
#
# Main output:
#   plot_rejection_by_n_pretty() returns the power/rejection-rate plot;
#   plot_variance_by_n_pretty() returns the estimated-variance plot;
#   save_hardness_plots() saves both plots as PNG files and returns
#   the plot objects and output file paths invisibly.
#
# Required package:
#   ggplot2
# ============================================================

# ------------------------------------------------------------
# Plot empirical rejection rate / power against labeled sample size.
#
# The function compares PPGCM (opt), PPGCM (w=1), and GCM using
# method-specific colors, point shapes, and connecting lines.
# Axis ranges and legend placement can be customized by the user.
# ------------------------------------------------------------
plot_rejection_by_n_pretty <- function(res_by_n,
                                       method_colors = c(
                                         "PPGCM (opt)" = "#1f77b4",
                                         "PPGCM (w=1)" = "#d62728",
                                         "GCM" = "#2ca02c"
                                       ),
                                       method_shapes = c(
                                         "PPGCM (opt)" = 16,
                                         "PPGCM (w=1)" = 17,
                                         "GCM" = 15
                                       ),
                                       y_limits = NULL,
                                       y_breaks = NULL,
                                       point_size = 2.8,
                                       line_width = 1.0,
                                       base_size = 13,
                                       legend_position = c(0.98, 0.98)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please run install.packages('ggplot2').")
  }
  
  # Extract and order the simulation summary used for plotting.
  s <- res_by_n$summary
  s <- s[order(s$n_label, s$method), , drop = FALSE]

  # Only change the plotting labels; do not change the original result object.
  s$method_plot <- gsub("^PPI-GCM", "PPGCM", s$method)
  
  # Define the plotting order and display labels of the methods.
  method_levels <- names(method_colors)
  
  method_labels <- c(
    "PPGCM (opt)" = "Naive PP-GCM (opt)",
    "PPGCM (w=1)" = "Naive PP-GCM (w=1)",
    "GCM" = "GCM"
  )
  
  s$method_plot <- factor(s$method_plot, levels = method_levels)
  x_breaks <- sort(unique(s$n_label))
  
  # Construct the power curve with one line for each testing method.
  p <- ggplot2::ggplot(
    s,
    ggplot2::aes(
      x = n_label,
      y = rejection_rate,
      color = method_plot,
      shape = method_plot,
      group = method_plot
    )
  ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(
      values = method_colors,
      breaks = method_levels,
      labels = unname(method_labels[method_levels]),
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = method_shapes,
      breaks = method_levels,
      labels = unname(method_labels[method_levels]),
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::labs(
      x = "n (the sample size of labeled data)",
      y = "Power",
      color = "Method",
      shape = "Method"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#E5E5E5", linewidth = 0.35),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.7),
      legend.position = legend_position,
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(
        fill = ggplot2::alpha("white", 0.88),
        color = "grey70",
        linewidth = 0.3
      ),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.box.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          shape = unname(method_shapes[method_levels]),
          linetype = "solid",
          linewidth = line_width,
          size = point_size
        )
      ),
      shape = "none"
    )
  
  # Apply optional user-specified y-axis ticks and plotting range.
  if (!is.null(y_breaks)) {
    p <- p + ggplot2::scale_y_continuous(breaks = y_breaks)
  }

  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }

  return(p)
}

# ------------------------------------------------------------
# Plot estimated asymptotic variance against labeled sample size.
#
# The plotted point is the mean estimated asymptotic variance across
# simulation repetitions. Error bars can represent either one standard
# deviation ("sd") or one standard error ("se").
# ------------------------------------------------------------
plot_variance_by_n_pretty <- function(res_by_n,
                                      method_colors = c(
                                        "PPGCM (opt)" = "#1f77b4",
                                        "PPGCM (w=1)" = "#d62728",
                                        "GCM" = "#2ca02c"
                                      ),
                                      method_shapes = c(
                                        "PPGCM (opt)" = 16,
                                        "PPGCM (w=1)" = 17,
                                        "GCM" = 15
                                      ),
                                      y_limits = NULL,
                                      y_breaks = NULL,
                                      errorbar_type = c("sd", "se"),
                                      errorbar_width = 0.02,
                                      point_size = 2.8,
                                      line_width = 1.0,
                                      errorbar_linewidth = 0.65,
                                      base_size = 13,
                                      legend_position = c(0.98, 0.98)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please run install.packages('ggplot2').")
  }

  errorbar_type <- match.arg(errorbar_type)
  
  # Extract and order the simulation summary used for plotting.
  s <- res_by_n$summary
  s <- s[order(s$n_label, s$method), , drop = FALSE]

  # Only change the plotting labels; do not change the original result object.
  s$method_plot <- gsub("^PPI-GCM", "PPGCM", s$method)
  
  # Define the plotting order and display labels of the methods.
  method_levels <- names(method_colors)
  method_labels <- c(
    "PPGCM (opt)" = "Naive PP-GCM (opt)",
    "PPGCM (w=1)" = "Naive PP-GCM (w=1)",
    "GCM" = "GCM"
  )
  s$method_plot <- factor(s$method_plot, levels = method_levels)
  x_breaks <- sort(unique(s$n_label))
  
  # Construct uncertainty bars using either the across-repetition standard
  # deviation or the corresponding standard error of the mean.
  if (errorbar_type == "sd") {
    s$err <- s$sd_variance_hat
  } else {
    s$err <- s$sd_variance_hat / sqrt(s$n_rep)
  }

  s$ymin <- pmax(s$mean_variance_hat - s$err, 0)
  s$ymax <- s$mean_variance_hat + s$err
  
  # Scale the horizontal width of the error bars to the range of n values.
  if (length(x_breaks) >= 2) {
    eb_width <- errorbar_width * diff(range(x_breaks))
  } else {
    eb_width <- 0.1
  }
  
  # Construct the variance curve and add uncertainty bars for each method.
  p <- ggplot2::ggplot(
    s,
    ggplot2::aes(
      x = n_label,
      y = mean_variance_hat,
      color = method_plot,
      shape = method_plot,
      group = method_plot
    )
  ) +
    ggplot2::geom_line(linewidth = line_width, linetype = "solid") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ymin, ymax = ymax),
      width = eb_width,
      linewidth = errorbar_linewidth
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(
      values = method_colors,
      breaks = method_levels,
      labels = unname(method_labels[method_levels]),
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = method_shapes,
      breaks = method_levels,
      labels = unname(method_labels[method_levels]),
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::labs(
      x = "n (the sample size of labeled data)",
      y = "Estimated asymptotic variance",
      color = "Method",
      shape = "Method"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#E5E5E5", linewidth = 0.35),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.7),
      legend.position = legend_position,
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(
        fill = ggplot2::alpha("white", 0.88),
        color = "grey70",
        linewidth = 0.3
      ),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.box.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          shape = unname(method_shapes[method_levels]),
          linetype = "solid",
          linewidth = line_width,
          size = point_size
        )
      ),
      shape = "none"
    )
  
  # Apply optional user-specified y-axis ticks and plotting range.
  if (!is.null(y_breaks)) {
    p <- p + ggplot2::scale_y_continuous(breaks = y_breaks)
  }

  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }

  return(p)
}

# ------------------------------------------------------------
# Generate and save both hardness-simulation plots.
#
# The function calls the two plotting functions above, saves the power
# and variance plots as PNG files, and invisibly returns the plot objects
# together with the corresponding output file paths.
# ------------------------------------------------------------
save_hardness_plots <- function(res_by_n,
                                out_dir = getwd(),
                                prefix = "hardness_oracle_noS",
                                width = 8.5,
                                height = 5.5,
                                dpi = 320,
                                rejection_y_limits = NULL,
                                rejection_y_breaks = NULL,
                                variance_y_limits = NULL,
                                variance_y_breaks = NULL,
                                errorbar_type = "sd") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please run install.packages('ggplot2').")
  }
  
  # Generate the power/rejection-rate and estimated-variance plots.
  p1 <- plot_rejection_by_n_pretty(
    res_by_n,
    y_limits = rejection_y_limits,
    y_breaks = rejection_y_breaks
  )

  p2 <- plot_variance_by_n_pretty(
    res_by_n,
    y_limits = variance_y_limits,
    y_breaks = variance_y_breaks,
    errorbar_type = errorbar_type
  )
  
  # Construct the output file names using the user-specified prefix.
  f1 <- file.path(out_dir, paste0(prefix, "_rejection_vs_n.png"))
  f2 <- file.path(out_dir, paste0(prefix, "_variance_vs_n.png"))
  
  # Save both figures using the requested size and resolution.
  ggplot2::ggsave(filename = f1, plot = p1, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(filename = f2, plot = p2, width = width, height = height, dpi = dpi)
  
  # Return the plot objects and file paths without printing them automatically.
  invisible(list(
    rejection_plot = p1,
    variance_plot = p2,
    rejection_file = f1,
    variance_file = f2
  ))
}
