#' Example: Replicating UK Manufacturing Analysis
#'
#' This script demonstrates how to use the gsf package to replicate
#' the analysis from Kumbhakar and Tsionas (2021).
#'
#' Note: Uses simulated data since the actual UK manufacturing data
#' is not publicly available.

library(gsf)

# =============================================================================
# Step 1: Generate Simulated UK Manufacturing Data
# =============================================================================

message("Generating simulated UK manufacturing data...")

set.seed(42)
uk_data <- simulate_uk_manufacturing(
  n_firms = 100,  # Paper uses 582 firms, but we use fewer for speed
  avg_periods = 9,
  seed = 42
)

message("Data generated: ", nrow(uk_data), " observations from ",
        length(unique(uk_data$firm_id)), " firms")

# Summary statistics
cat("\nSummary Statistics:\n")
print(summary(uk_data[, c("Y", "L", "K", "CR", "MKSH", "FP", "IMP")]))

# True efficiency (from simulation)
true_eff <- attr(uk_data, "true_efficiency")
cat("\nTrue Technical Inefficiency (%):\n")
cat("  Mean:", mean(true_eff$technical_inefficiency), "\n")
cat("  SD:", sd(true_eff$technical_inefficiency), "\n")

cat("\nTrue Labor Slack (%):\n")
cat("  Mean:", mean(true_eff$input_slacks[, 1]), "\n")
cat("  SD:", sd(true_eff$input_slacks[, 1]), "\n")

cat("\nTrue Capital Slack (%):\n")
cat("  Mean:", mean(true_eff$input_slacks[, 2]), "\n")
cat("  SD:", sd(true_eff$input_slacks[, 2]), "\n")

# =============================================================================
# Step 2: Fit the GSF Model
# =============================================================================

message("\nFitting GSF model (this may take 10-30 minutes)...")

fit <- gsf_fit(
  formula = log(Y) ~ log(L) + log(K),
  data = uk_data,
  z_vars = c("CR", "MKSH", "FP", "IMP", "TREND"),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = parallel::detectCores() - 1,
  seed = 123
)

# =============================================================================
# Step 3: Model Summary
# =============================================================================

cat("\n")
cat("=================================================================\n")
cat("                    MODEL SUMMARY                                \n")
cat("=================================================================\n")
summary(fit)

# =============================================================================
# Step 4: Efficiency Measures (Table 2 in paper)
# =============================================================================

cat("\n")
cat("=================================================================\n")
cat("            EFFICIENCY MEASURES (Table 2)                        \n")
cat("=================================================================\n")

eff_table <- efficiency_summary(fit)
print(eff_table)

# Compare with true values
cat("\nComparison with True Values:\n")
cat("----------------------------\n")

# Technical inefficiency
ti <- extract_technical_inefficiency(fit)
cat("Technical Inefficiency:\n")
cat("  Estimated Mean:", round(mean(ti$mean), 3), "%\n")
cat("  True Mean:", round(mean(true_eff$technical_inefficiency), 3), "%\n")

# Input slacks
slacks <- extract_input_slacks(fit)
labor_slack <- slacks[slacks$input == "log(L)", ]
capital_slack <- slacks[slacks$input == "log(K)", ]

cat("\nLabor Slack:\n")
cat("  Estimated Mean:", round(mean(labor_slack$mean), 3), "%\n")
cat("  True Mean:", round(mean(true_eff$input_slacks[, 1]), 3), "%\n")

cat("\nCapital Slack:\n")
cat("  Estimated Mean:", round(mean(capital_slack$mean), 3), "%\n")
cat("  True Mean:", round(mean(true_eff$input_slacks[, 2]), 3), "%\n")

# =============================================================================
# Step 5: Delta Coefficients (Table 3 in paper)
# =============================================================================

cat("\n")
cat("=================================================================\n")
cat("        SLACK DETERMINANTS (Delta Coefficients)                  \n")
cat("=================================================================\n")

delta_coef <- extract_delta_coefficients(fit)
print(delta_coef)

cat("\nInterpretation:\n")
cat("- Positive coefficients: Variable increases slack (more inefficient)\n")
cat("- Negative coefficients: Variable decreases slack (more efficient)\n")
cat("- Note: Slack values are negative, so positive Delta means more negative slack\n")

# =============================================================================
# Step 6: Returns to Scale
# =============================================================================

cat("\n")
cat("=================================================================\n")
cat("                   RETURNS TO SCALE                              \n")
cat("=================================================================\n")

rts <- extract_returns_to_scale(fit)
cat("Mean RTS:", round(mean(rts$mean), 3), "\n")
cat("SD RTS:", round(sd(rts$mean), 3), "\n")
cat("Percent with DRS:", round(100 * mean(rts$mean < 1), 1), "%\n")
cat("Percent with IRS:", round(100 * mean(rts$mean > 1), 1), "%\n")

# =============================================================================
# Step 7: Visualizations (Figures 2, 3, 4 in paper)
# =============================================================================

cat("\n")
cat("=================================================================\n")
cat("                   CREATING PLOTS                                \n")
cat("=================================================================\n")

# Figure 2a,c: Input slacks and technical inefficiency
p_eff <- plot(fit, type = "efficiency")
print(p_eff)
cat("Created: Efficiency distributions plot\n")

# Figure 2b,d: Output loss
p_loss <- plot(fit, type = "output_loss")
print(p_loss)
cat("Created: Output loss plot\n")

# Figure 4: Input elasticities
p_elas <- plot(fit, type = "elasticities")
print(p_elas)
cat("Created: Input elasticities plot\n")

# Returns to scale distribution
p_rts <- plot(fit, type = "rts")
print(p_rts)
cat("Created: Returns to scale plot\n")

# Marginal effects of z variables
p_marginal <- plot_marginal_effects(fit)
print(p_marginal)
cat("Created: Marginal effects plot\n")

# Figure 3: Best and worst firms
best_worst <- plot_best_worst_firms(fit, n_firms = 5)
print(best_worst$worst)
print(best_worst$best)
cat("Created: Best/worst firms comparison\n")

# =============================================================================
# Step 8: Model Comparison (Optional)
# =============================================================================

cat("\n")
cat("=================================================================\n")
cat("                MODEL COMPARISON METRICS                         \n")
cat("=================================================================\n")

# WAIC
waic_result <- waic(fit)
cat("WAIC:", round(waic_result$waic, 2), "\n")
cat("LPPD:", round(waic_result$lppd, 2), "\n")
cat("p_WAIC:", round(waic_result$p_waic, 2), "\n")

# LOO-CV (if loo package is available)
if (requireNamespace("loo", quietly = TRUE)) {
  loo_result <- loo_cv(fit)
  cat("\nLOO-CV:\n")
  print(loo_result)
}

# =============================================================================
# Step 9: Save Results
# =============================================================================

cat("\n")
cat("=================================================================\n")
cat("                   SAVING RESULTS                                \n")
cat("=================================================================\n")

# Save the fit object
saveRDS(fit, "gsf_fit_uk_manufacturing.rds")
cat("Saved model fit to: gsf_fit_uk_manufacturing.rds\n")

# Save efficiency estimates
results <- list(
  technical_inefficiency = ti,
  input_slacks = slacks,
  delta_coefficients = delta_coef,
  returns_to_scale = rts,
  efficiency_summary = eff_table
)
saveRDS(results, "gsf_results_uk_manufacturing.rds")
cat("Saved results to: gsf_results_uk_manufacturing.rds\n")

# Save plots
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggsave("efficiency_distributions.png", p_eff, width = 10, height = 6)
  ggplot2::ggsave("output_loss.png", p_loss, width = 10, height = 6)
  ggplot2::ggsave("elasticities.png", p_elas, width = 10, height = 6)
  ggplot2::ggsave("returns_to_scale.png", p_rts, width = 8, height = 6)
  ggplot2::ggsave("marginal_effects.png", p_marginal, width = 10, height = 6)
  cat("Saved plots as PNG files\n")
}

cat("\n")
cat("=================================================================\n")
cat("                     ANALYSIS COMPLETE                           \n")
cat("=================================================================\n")
