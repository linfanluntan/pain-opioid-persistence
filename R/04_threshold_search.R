# =============================================================================
# 04_threshold_search.R
# Pain AUC & Opioid Persistence — Grid Search for Clinical AUC Cutpoints
# ASTRO 2026 | Mahin, Nkuku, He, Fuller, Moreno, Javed
# =============================================================================
#
# Exploratory risk-stratification analysis:
#   For each candidate threshold, binarize AUC and fit a multivariable
#   logistic model to estimate the OR for above-threshold vs. below.
#
# NOTE: The primary inference uses continuous AUC (more powerful).
#       Thresholds are secondary — for clinical decision-making only.
# =============================================================================

library(logistf)
library(dplyr)

# ─── Grid search function ───────────────────────────────────────────────────
threshold_grid_search <- function(data, outcome_col, auc_col, adjust_cols,
                                   n_thresholds = 60,
                                   pctl_range = c(0.05, 0.95),
                                   use_firth = TRUE) {
  # data:         analytic data frame
  # outcome_col:  name of binary outcome (0/1)
  # auc_col:      name of continuous AUC variable (in %)
  # adjust_cols:  character vector of adjustment covariates
  # n_thresholds: number of evenly spaced thresholds to test
  # pctl_range:   percentile range for thresholds
  # use_firth:    if TRUE, use Firth; if FALSE, use standard GLM
  
  auc_values <- data[[auc_col]]
  outcome    <- data[[outcome_col]]
  
  # Define threshold grid
  lo <- quantile(auc_values, pctl_range[1], na.rm = TRUE)
  hi <- quantile(auc_values, pctl_range[2], na.rm = TRUE)
  thresholds <- seq(lo, hi, length.out = n_thresholds)
  
  cat("Testing", length(thresholds), "thresholds from",
      round(lo, 1), "to", round(hi, 1), "\n\n")
  
  results <- list()
  
  for (t in thresholds) {
    # Create binary above/below threshold
    data$AUC_above <- as.integer(auc_values >= t)
    
    n_low  <- sum(data$AUC_above == 0)
    n_high <- sum(data$AUC_above == 1)
    
    # Skip if either group too small
    if (n_low < 10 || n_high < 10) next
    
    # Build formula
    rhs <- paste(c("AUC_above", adjust_cols), collapse = " + ")
    formula <- as.formula(paste(outcome_col, "~", rhs))
    
    tryCatch({
      if (use_firth) {
        m <- logistf(formula, data = data)
        beta <- coef(m)["AUC_above"]
        p_val <- m$prob["AUC_above"]
      } else {
        m <- glm(formula, data = data, family = binomial())
        s <- summary(m)
        beta <- coef(m)["AUC_above"]
        p_val <- s$coefficients["AUC_above", "Pr(>|z|)"]
      }
      
      results <- c(results, list(tibble(
        threshold = round(t, 1),
        n_low     = n_low,
        n_high    = n_high,
        beta      = round(beta, 3),
        OR        = round(exp(beta), 2),
        p_value   = round(p_val, 4),
        significant = p_val < 0.05
      )))
    }, error = function(e) {
      # Skip thresholds that cause model failure
    })
  }
  
  result_df <- bind_rows(results)
  
  # Summary
  sig_results <- result_df %>% filter(significant)
  
  if (nrow(sig_results) > 0) {
    best <- sig_results %>% arrange(desc(abs(OR))) %>% slice(1)
    lowest_sig <- sig_results %>% arrange(threshold) %>% slice(1)
    
    cat("--- Threshold Search Results ---\n")
    cat("Total thresholds tested:", nrow(result_df), "\n")
    cat("Significant thresholds:", nrow(sig_results), "\n\n")
    
    cat("Lowest significant threshold:\n")
    cat("  AUC >=", lowest_sig$threshold, "% → OR =", lowest_sig$OR,
        "(p =", lowest_sig$p_value, ")\n")
    
    cat("\nStrongest significant association:\n")
    cat("  AUC >=", best$threshold, "% → OR =", best$OR,
        "(p =", best$p_value, ")\n")
    cat("  N below:", best$n_low, "| N above:", best$n_high, "\n")
  } else {
    cat("No significant thresholds found.\n")
  }
  
  result_df
}

# ─── Internal validation (bootstrap) ────────────────────────────────────────
bootstrap_threshold <- function(data, outcome_col, auc_col, adjust_cols,
                                 optimal_threshold, n_boot = 500,
                                 use_firth = TRUE) {
  # Bootstrap the optimal threshold's OR to assess stability
  
  cat("\n--- Bootstrap Validation (n =", n_boot, ") ---\n")
  cat("Testing stability of threshold:", optimal_threshold, "\n\n")
  
  boot_ors <- numeric(n_boot)
  
  for (b in seq_len(n_boot)) {
    idx <- sample(nrow(data), replace = TRUE)
    boot_data <- data[idx, ]
    
    boot_data$AUC_above <- as.integer(boot_data[[auc_col]] >= optimal_threshold)
    
    rhs <- paste(c("AUC_above", adjust_cols), collapse = " + ")
    formula <- as.formula(paste(outcome_col, "~", rhs))
    
    tryCatch({
      if (use_firth) {
        m <- logistf(formula, data = boot_data)
        boot_ors[b] <- exp(coef(m)["AUC_above"])
      } else {
        m <- glm(formula, data = boot_data, family = binomial())
        boot_ors[b] <- exp(coef(m)["AUC_above"])
      }
    }, error = function(e) {
      boot_ors[b] <- NA
    })
  }
  
  boot_ors <- boot_ors[!is.na(boot_ors)]
  
  cat("Bootstrap OR — Mean:", round(mean(boot_ors), 2), "\n")
  cat("Bootstrap OR — Median:", round(median(boot_ors), 2), "\n")
  cat("Bootstrap 95% CI:", round(quantile(boot_ors, 0.025), 2), "–",
      round(quantile(boot_ors, 0.975), 2), "\n")
  cat("Proportion OR > 1:", round(mean(boot_ors > 1), 3), "\n")
  
  invisible(boot_ors)
}

# =============================================================================
# Main execution
# =============================================================================
cat("=== AUC Threshold Search — Exploratory Risk Stratification ===\n\n")

# --- PLACEHOLDER: Load your data ---
# fu1 <- readRDS("data/processed/fu1_analytic.rds")
# fu2 <- readRDS("data/processed/fu2_analytic.rds")
#
# # FU1 threshold search
# cat(">>> FU1 Threshold Search <<<\n")
# fu1_thresholds <- threshold_grid_search(
#   data = fu1,
#   outcome_col = "Opioid_Active_At_FU1Date",
#   auc_col = "AUCpain_acute_pct",
#   adjust_cols = c("Smk_Ever", "Age60plus", "TxIntensity"),
#   use_firth = FALSE  # FU1 is large enough for standard GLM
# )
#
# # FU2 threshold search (use Firth due to small N)
# cat("\n>>> FU2 Threshold Search <<<\n")
# fu2_thresholds <- threshold_grid_search(
#   data = fu2,
#   outcome_col = "Opioid_Active_At_FU2Date",
#   auc_col = "AUCpain_acute_pct",
#   adjust_cols = c("Smk_Ever", "Age60plus", "TxIntensity"),
#   use_firth = TRUE
# )
#
# # Bootstrap validation of best FU1 threshold
# boot_fu1 <- bootstrap_threshold(
#   data = fu1,
#   outcome_col = "Opioid_Active_At_FU1Date",
#   auc_col = "AUCpain_acute_pct",
#   adjust_cols = c("Smk_Ever", "Age60plus", "TxIntensity"),
#   optimal_threshold = 53.2,  # from grid search
#   n_boot = 500
# )

# --- Demo ---
cat("--- Running demo with simulated data ---\n")
set.seed(42)
n <- 300
demo_df <- data.frame(
  Outcome     = rbinom(n, 1, 0.4),
  AUC_pct     = pmax(0, pmin(70, rnorm(n, 30, 15))),
  Smk_Ever    = factor(rbinom(n, 1, 0.4)),
  Age60plus   = factor(rbinom(n, 1, 0.5))
)
# Create slight association
demo_df$Outcome[demo_df$AUC_pct > 45] <-
  rbinom(sum(demo_df$AUC_pct > 45), 1, 0.65)

demo_thresholds <- threshold_grid_search(
  data = demo_df,
  outcome_col = "Outcome",
  auc_col = "AUC_pct",
  adjust_cols = c("Smk_Ever", "Age60plus"),
  n_thresholds = 30,
  use_firth = FALSE
)

cat("\n✓ Demo complete.\n")
