# ==============================================================================
# hestia: Vignette Figure Rapid Prototyping Harness
#
# This script enables rapid design, rendering, and visual inspection of
# proposed figures for vignettes without refitting Stan models.
# Pre-computed package data objects (sir_res, sir_cov_res, siir_res, sir)
# are used as inputs.
#
# Usage:
#   Rscript vignettes/prototype_figures.R [--out-dir <path>]
#   Or source() directly in an interactive R session.
# ==============================================================================

suppressPackageStartupMessages({
  library(hestia)
  library(posterior)
  library(ggplot2)
  library(dplyr)
  library(bayesplot)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# Default Plot Styling
# ------------------------------------------------------------------------------
theme_hestia <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(color = "grey40", size = rel(0.9)),
      plot.caption = element_text(color = "grey50", size = rel(0.8), hjust = 0),
      legend.position = "bottom",
      strip.text = element_text(face = "bold", hjust = 0)
    )
}

# ------------------------------------------------------------------------------
# Prototype 1: SIR Parameter Recovery (vignettes/SIR.Rmd)
# ------------------------------------------------------------------------------
proto_sir_params <- function() {
  data("sir_res", package = "hestia", envir = environment())

  # Summary of posterior draws
  draws_df <- as_draws_df(sir_res)

  # Parameter labels and true simulation values
  param_meta <- tibble::tribble(
    ~variable,  ~label,                                    ~true_val, ~category,
    "eh_prob",  "Extra-household infection prob (daily)",  0.01,      "Infection Probability",
    "ih_prob",  "Intra-household infection prob (daily)",  0.05,      "Infection Probability",
    "gamma",    "Recovery rate (per day)",                 0.20,      "Transition Rate"
  )

  draws_long <- draws_df |>
    tidyr::pivot_longer(cols = c("eh_prob", "ih_prob", "gamma"),
                        names_to = "variable", values_to = "value") |>
    left_join(param_meta, by = "variable")

  summary_stats <- draws_long |>
    group_by(variable, label, category, true_val) |>
    summarise(
      median = median(value),
      q025 = quantile(value, 0.025),
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75),
      q975 = quantile(value, 0.975),
      .groups = "drop"
    )

  p <- ggplot(summary_stats, aes(y = label)) +
    # 95% CrI
    geom_errorbar(aes(xmin = q025, xmax = q975), width = 0, color = "#2b5c8f", linewidth = 1) +
    # 50% CrI
    geom_errorbar(aes(xmin = q25, xmax = q75), width = 0, color = "#183e6b", linewidth = 2.5) +
    # Posterior median
    geom_point(aes(x = median), color = "#0f2540", size = 3.5, shape = 21, fill = "white", stroke = 1.8) +
    # True value reference line
    geom_point(aes(x = true_val, color = "True Simulation Value"), size = 3, shape = 4, stroke = 2) +
    facet_wrap(~category, scales = "free", ncol = 1) +
    scale_color_manual(name = "", values = c("True Simulation Value" = "#d95f02")) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
    labs(
      title = "Posterior Parameter Estimates (SIR Model)",
      subtitle = "Points show posterior median; thick/thin bars show 50% & 95% credible intervals",
      x = "Parameter Value (Natural Scale)",
      y = NULL
    ) +
    theme_hestia() +
    theme(
      legend.position = "top",
      legend.justification = "right",
      panel.grid.major.y = element_line(color = "grey90")
    )

  p
}

# ------------------------------------------------------------------------------
# Prototype 2: SIR Covariates Effects (vignettes/SIR.Rmd)
# ------------------------------------------------------------------------------
proto_sir_cov_params <- function() {
  data("sir_cov_res", package = "hestia", envir = environment())

  draws_df <- as_draws_df(sir_cov_res)

  cov_meta <- tibble::tribble(
    ~variable, ~label,                            ~true_val,  ~category,
    "eh_prob", "Extra-household baseline prob",   0.01,       "Baseline Probabilities",
    "ih_prob", "Intra-household baseline prob",   0.05,       "Baseline Probabilities",
    "gamma",   "Recovery rate (gamma)",           0.20,       "Recovery Rate",
    "x1_eh",   "x1 effect on extra-household",    exp(-0.4),  "Relative Risk / Odds Ratio",
    "x2_eh",   "x2 effect on extra-household",    exp(0.7),   "Relative Risk / Odds Ratio",
    "x1_ih",   "x1 effect on intra-household",    exp(0.8),   "Relative Risk / Odds Ratio",
    "x2_ih",   "x2 effect on intra-household",    exp(0.1),   "Relative Risk / Odds Ratio"
  )

  draws_long <- draws_df |>
    tidyr::pivot_longer(cols = all_of(cov_meta$variable),
                        names_to = "variable", values_to = "value") |>
    left_join(cov_meta, by = "variable")

  summary_stats <- draws_long |>
    group_by(variable, label, category, true_val) |>
    summarise(
      median = median(value),
      q025 = quantile(value, 0.025),
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75),
      q975 = quantile(value, 0.975),
      .groups = "drop"
    )

  p <- ggplot(summary_stats, aes(y = label)) +
    # Null effect reference line for relative risk
    geom_vline(data = filter(summary_stats, category == "Relative Risk / Odds Ratio"),
               aes(xintercept = 1), linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = q025, xmax = q975), width = 0, color = "#1b9e77", linewidth = 1) +
    geom_errorbar(aes(xmin = q25, xmax = q75), width = 0, color = "#0d664c", linewidth = 2.5) +
    geom_point(aes(x = median), color = "#0d664c", size = 3, shape = 21, fill = "white", stroke = 1.5) +
    geom_point(aes(x = true_val, color = "True Simulation Value"), size = 3, shape = 4, stroke = 2) +
    facet_wrap(~category, scales = "free", ncol = 1) +
    scale_color_manual(name = "", values = c("True Simulation Value" = "#d95f02")) +
    labs(
      title = "Posterior Estimates with Covariates (SIR Model)",
      subtitle = "Coefficients are on the exponentiated (natural) scale",
      x = "Estimate",
      y = NULL
    ) +
    theme_hestia() +
    theme(
      legend.position = "top",
      legend.justification = "right",
      panel.grid.major.y = element_line(color = "grey90")
    )

  p
}

# ------------------------------------------------------------------------------
# Prototype 3: Multiple Infection Compartments (vignettes/multiple_infection_compartments.Rmd)
# ------------------------------------------------------------------------------
proto_siir_params <- function() {
  data("siir_res", package = "hestia", envir = environment())

  draws_df <- as_draws_df(siir_res)

  siir_meta <- tibble::tribble(
    ~variable,     ~label,                                    ~true_val, ~category,
    "eh_prob",     "Extra-household infection prob",          0.01,      "Infection Probabilities",
    "ih_prob_Is",  "Intra-household (Symptomatic)",           0.05,      "Infection Probabilities",
    "ih_prob_Ia",  "Intra-household (Asymptomatic)",          0.025,     "Infection Probabilities",
    "gamma_s",     "Recovery rate (Symptomatic)",             0.20,      "Recovery Rates",
    "gamma_a",     "Recovery rate (Asymptomatic)",            1/3,       "Recovery Rates",
    "phi",         "Symptomatic proportion (split)",          0.70,      "Symptomatic Proportion"
  )

  draws_long <- draws_df |>
    tidyr::pivot_longer(cols = all_of(siir_meta$variable),
                        names_to = "variable", values_to = "value") |>
    left_join(siir_meta, by = "variable")

  summary_stats <- draws_long |>
    group_by(variable, label, category, true_val) |>
    summarise(
      median = median(value),
      q025 = quantile(value, 0.025),
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75),
      q975 = quantile(value, 0.975),
      .groups = "drop"
    )

  p <- ggplot(summary_stats, aes(y = label)) +
    geom_errorbar(aes(xmin = q025, xmax = q975), width = 0, color = "#7570b3", linewidth = 1) +
    geom_errorbar(aes(xmin = q25, xmax = q75), width = 0, color = "#4d4785", linewidth = 2.5) +
    geom_point(aes(x = median), color = "#4d4785", size = 3, shape = 21, fill = "white", stroke = 1.5) +
    geom_point(aes(x = true_val, color = "True Simulation Value"), size = 3, shape = 4, stroke = 2) +
    facet_wrap(~category, scales = "free", ncol = 1) +
    scale_color_manual(name = "", values = c("True Simulation Value" = "#d95f02")) +
    labs(
      title = "SIIR Model Estimates: Symptomatic vs. Asymptomatic",
      subtitle = "Intra-household transmission and recovery stratified by compartment",
      x = "Parameter Value",
      y = NULL
    ) +
    theme_hestia() +
    theme(
      legend.position = "top",
      legend.justification = "right",
      panel.grid.major.y = element_line(color = "grey90")
    )

  p
}

# ------------------------------------------------------------------------------
# Prototype 4: Simulated Data Trajectory (vignettes/SIR.Rmd)
# ------------------------------------------------------------------------------
proto_sir_data_trajectory <- function() {
  data("sir", package = "hestia", envir = environment())

  obs_trajectory <- sir |>
    group_by(t) |>
    summarise(
      pcr_pos = mean(pcr, na.rm = TRUE),
      igg_pos = mean(igg, na.rm = TRUE),
      n_obs = n(),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(cols = c(pcr_pos, igg_pos),
                        names_to = "test_type", values_to = "positivity") |>
    mutate(test_label = ifelse(test_type == "pcr_pos", "PCR Positive (Active)", "IgG Positive (Recovered)"))

  p <- ggplot(obs_trajectory, aes(x = t, y = positivity, color = test_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5, shape = 21, aes(fill = test_label), color = "white") +
    scale_color_manual(values = c("PCR Positive (Active)" = "#e41a1c", "IgG Positive (Recovered)" = "#377eb8")) +
    scale_fill_manual(values = c("PCR Positive (Active)" = "#e41a1c", "IgG Positive (Recovered)" = "#377eb8")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Simulated SIR Epidemic Dynamics Across Households",
      subtitle = "Proportion of enrolled individuals testing positive over observation days",
      x = "Day (t)",
      y = "Sample Positivity Rate",
      color = NULL,
      fill = NULL
    ) +
    theme_hestia() +
    theme(
      legend.position = "top",
      legend.justification = "left"
    )

  p
}

# ------------------------------------------------------------------------------
# Render and Save Helper
# ------------------------------------------------------------------------------
render_all_prototypes <- function(out_dirs = c("vignettes/figures")) {
  for (d in out_dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  plots <- list(
    "sir_param_recovery.png" = list(plot = proto_sir_params(), width = 7.5, height = 4.8),
    "sir_cov_effects.png"    = list(plot = proto_sir_cov_params(), width = 7.5, height = 5.8),
    "siir_stratified.png"    = list(plot = proto_siir_params(), width = 7.5, height = 5.5),
    "sir_trajectory.png"     = list(plot = proto_sir_data_trajectory(), width = 7.5, height = 4.2)
  )

  saved_paths <- c()
  for (fn in names(plots)) {
    item <- plots[[fn]]
    for (d in out_dirs) {
      target <- file.path(d, fn)
      ggsave(target, plot = item$plot, width = item$width, height = item$height, dpi = 300)
      saved_paths <- c(saved_paths, target)
      message("Rendered: ", target)
    }
  }

  invisible(saved_paths)
}

# ------------------------------------------------------------------------------
# Script Execution
# ------------------------------------------------------------------------------
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir_arg <- "vignettes/figures"
  if (length(args) >= 2 && args[1] == "--out-dir") {
    out_dir_arg <- args[2]
  }

  render_all_prototypes(out_dirs = c(out_dir_arg))
}
