# =============================================================================
# 03_l2_sensitivity.R
# Pain AUC & Opioid Persistence — L2 Regularized Logistic Regression
# ASTRO 2026 | Mahin, Nkuku, He, Fuller, Moreno, Javed
# =============================================================================
#
# L2 (ridge) regularized logistic regression as SENSITIVITY analysis.
#
# Purpose: Confirm that the primary Pain AUC effect is directionally
# consistent under an alternative shrinkage approach.
#
# This is NOT the primary inference model. Firth is primary.
# L2 is used here for robustness checking only.
#
# REQUIRES: install.packages("glmnet")
# =============================================================================

library(glmnet)
library(dplyr)

# ─── Configuration ───────────────────────────────────────────────────────────
DATA_DIR   <- "data/processed/"
OUTPUT_DIR <- "results/"

# ─── L2 Sensitivity Model ───────────────────────────────────────────────────
run_l2_sensitivity <- function(data, outcome_col, predictor_cols,
                                model_name = "L2 Sensitivity") {
  cat("\n", strrep("=", 60), "\n")
  cat(" ", model_name, "\n")
  cat(strrep("=", 60), "\n\n")
  
  # Prepare design matrix (glmnet needs matrix + vector)
  formula <- as.formula(paste(outcome_col, "~",
                               paste(predictor_cols, collapse = " + ")))
  
  X <- model.matrix(formula, data = data)[, -1]  # drop intercept (glmnet adds it)
  y <- data[[outcome_col]]
  
  cat("Design matrix:", nrow(X), "x", ncol(X), "\n")
  cat("Events:", sum(y == 1), "| Non-events:", sum(y == 0), "\n\n")
  
  # Cross-validated L2 (alpha = 0 = ridge)
  set.seed(42)
  cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 0,
                       type.measure = "deviance", nfolds = 5)
  
  cat("Optimal lambda (min deviance):", round(cv_fit$lambda.min, 5), "\n")
  cat("Lambda 1SE:", round(cv_fit$lambda.1se, 5), "\n\n")
  
  # Coefficients at lambda.min
  coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
  
  results <- data.frame(
    Term = rownames(coefs),
    Beta_L2 = round(as.numeric(coefs), 4),
    OR_L2   = round(exp(as.numeric(coefs)), 3),
    row.names = NULL
  )
  
  cat("--- L2 Coefficients (lambda.min) ---\n")
  print(results, row.names = FALSE)
  
  # Predictions
  pred <- predict(cv_fit, X, s = "lambda.min", type = "response")
  pred_class <- ifelse(pred > 0.5, 1, 0)
  accuracy <- mean(pred_class == y) * 100
  cat("\nClassification accuracy:", round(accuracy, 1), "%\n")
  
  # AUC-ROC
  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(y, as.numeric(pred), quiet = TRUE)
    cat("AUC-ROC:", round(pROC::auc(roc_obj), 3), "\n")
  }
  
  list(cv_fit = cv_fit, results = results, accuracy = accuracy)
}

# ─── Compare Firth vs L2 ────────────────────────────────────────────────────
compare_firth_l2 <- function(firth_results, l2_results) {
  cat("\n", strrep("=", 60), "\n")
  cat("  Firth vs. L2 — Directional Consistency Check\n")
  cat(strrep("=", 60), "\n\n")
  
  # Merge by term name
  merged <- merge(
    firth_results %>% select(Term, OR_Firth = OR, p_Firth = p_value),
    l2_results   %>% select(Term, OR_L2),
    by = "Term",
    all = TRUE
  )
  
  merged$Direction_Firth <- ifelse(merged$OR_Firth > 1, "+", "-")
  merged$Direction_L2    <- ifelse(merged$OR_L2 > 1, "+", "-")
  merged$Consistent      <- merged$Direction_Firth == merged$Direction_L2
  
  print(merged, row.names = FALSE)
  
  n_consistent <- sum(merged$Consistent, na.rm = TRUE)
  n_total      <- sum(!is.na(merged$Consistent))
  cat("\nDirectional consistency:", n_consistent, "/", n_total, "terms\n")
  
  invisible(merged)
}

# =============================================================================
# Main execution
# =============================================================================
cat("=== L2 Regularized Logistic Regression — Sensitivity Analysis ===\n\n")

# --- PLACEHOLDER: Load your data ---
# fu2 <- readRDS(file.path(DATA_DIR, "fu2_analytic.rds"))
#
# l2_res <- run_l2_sensitivity(
#   data = fu2,
#   outcome_col = "Opioid_Active_At_FU2Date",
#   predictor_cols = c("AUC10", "Smk_Ever", "Age60plus", "TxIntensity"),
#   model_name = "FU2 L2 Sensitivity"
# )
#
# # Compare with Firth results (load from 02_firth_primary_models.R output)
# firth_res <- readRDS(file.path(OUTPUT_DIR, "fu2_firth_results.rds"))
# compare_firth_l2(firth_res$results, l2_res$results)

# --- Demo ---
cat("--- Running demo with simulated data ---\n")
set.seed(42)
n <- 200
demo_df <- data.frame(
  Outcome     = rbinom(n, 1, 0.4),
  AUC10       = rnorm(n, mean = 3, sd = 1.5),
  Smk_Ever    = factor(rbinom(n, 1, 0.4)),
  Age60plus   = factor(rbinom(n, 1, 0.5)),
  TxIntensity = factor(sample(1:3, n, replace = TRUE))
)

l2_demo <- run_l2_sensitivity(
  data = demo_df,
  outcome_col = "Outcome",
  predictor_cols = c("AUC10", "Smk_Ever", "Age60plus", "TxIntensity"),
  model_name = "Demo: Simulated L2 Model"
)

cat("\n✓ Demo complete.\n")
