# =============================================================================
# 02_firth_primary_models.R
# Pain AUC & Opioid Persistence — Firth Penalized Logistic Regression
# ASTRO 2026 | Mahin, Nkuku, He, Fuller, Moreno, Javed
# =============================================================================
#
# Primary inference models using Firth bias-reduced penalized LR.
#
# Model A (FU1): Full multivariable model — larger cohort supports more params
# Model B (FU2): Prespecified parsimonious model — ~4-6 effective parameters
# Model B2 (FU2 sensitivity): Model B + one behavioral covariate
#
# Firth regression adds a Jeffreys-prior penalty:
#   l*(β) = l(β) + ½ log|I(β)|
# Guarantees finite estimates under separation, corrects small-sample bias.
#
# REQUIRES: install.packages("logistf")
# =============================================================================

library(logistf)
library(dplyr)

# ─── Configuration ───────────────────────────────────────────────────────────
DATA_DIR   <- "data/processed/"
OUTPUT_DIR <- "results/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ─── Helper: Format Firth results ───────────────────────────────────────────
format_firth_results <- function(model, model_name = "Model") {
  beta <- coef(model)
  ci   <- confint(model)
  pval <- model$prob
  
  results <- data.frame(
    Term    = names(beta),
    Beta    = round(beta, 4),
    OR      = round(exp(beta), 3),
    CI_low  = round(exp(ci[, 1]), 3),
    CI_high = round(exp(ci[, 2]), 3),
    p_value = round(pval, 4),
    row.names = NULL
  )
  
  cat("\n===", model_name, "===\n")
  cat("N =", model$n, "\n")
  cat("Penalized log-likelihood:", round(model$loglik[2], 3), "\n\n")
  print(results, row.names = FALSE)
  
  invisible(results)
}

# ─── Helper: Compare MLE vs Firth ───────────────────────────────────────────
compare_mle_firth <- function(formula, data, model_name = "Model") {
  cat("\n", strrep("=", 60), "\n")
  cat(" ", model_name, "— MLE vs. Firth Comparison\n")
  cat(strrep("=", 60), "\n")
  
  # Standard MLE (may warn/fail)
  cat("\n--- Standard MLE ---\n")
  m_mle <- tryCatch(
    glm(formula, data = data, family = binomial()),
    warning = function(w) {
      cat("WARNING:", conditionMessage(w), "\n")
      suppressWarnings(glm(formula, data = data, family = binomial()))
    },
    error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(m_mle)) {
    mle_summary <- summary(m_mle)
    cat("Converged:", m_mle$converged, "\n")
    cat("Max |OR|:", round(max(abs(exp(coef(m_mle)))), 1), "\n")
    
    # Check for separation indicators
    max_se <- max(mle_summary$coefficients[, "Std. Error"])
    if (max_se > 10) {
      cat("⚠ SEPARATION DETECTED: Max SE =", round(max_se, 1), "\n")
      cat("  Standard MLE estimates are UNRELIABLE.\n")
    }
  }
  
  # Firth
  cat("\n--- Firth Penalized LR ---\n")
  m_firth <- logistf(formula, data = data)
  results <- format_firth_results(m_firth, model_name)
  
  # Classification accuracy
  pred_prob <- predict(m_firth, type = "response")
  pred_class <- ifelse(pred_prob > 0.5, 1, 0)
  actual <- data[[all.vars(formula)[1]]]
  accuracy <- mean(pred_class == actual, na.rm = TRUE) * 100
  cat("\nClassification accuracy:", round(accuracy, 1), "%\n")
  
  list(mle = m_mle, firth = m_firth, results = results)
}

# =============================================================================
# MODEL A: FU1 — Full Multivariable Model
# =============================================================================
run_fu1_model <- function(fu1) {
  cat("\n", strrep("#", 60), "\n")
  cat(" FU1 MODELS (Subacute Opioid Persistence)\n")
  cat(strrep("#", 60), "\n")
  
  cat("\nFU1 cohort: N =", nrow(fu1), "\n")
  cat("Events:", sum(fu1$Opioid_Active_At_FU1Date == 1), "\n")
  cat("Non-events:", sum(fu1$Opioid_Active_At_FU1Date == 0), "\n")
  
  # FU1 can support more predictors (larger N)
  formula_fu1 <- Opioid_Active_At_FU1Date ~
    AUC10 +
    Smk_Ever +
    Age60plus +
    TxIntensity +
    Gender +
    Race2 +
    Site2
  
  # Note: If Race has 4+ levels, use collapse_race(., n_levels=2) from prep script
  
  compare_mle_firth(formula_fu1, fu1, "Model A: FU1 Opioid Persistence")
}

# =============================================================================
# MODEL B: FU2 — Prespecified Parsimonious Model (PRIMARY)
# =============================================================================
run_fu2_model <- function(fu2) {
  cat("\n", strrep("#", 60), "\n")
  cat(" FU2 MODELS (Chronic Opioid Persistence)\n")
  cat(strrep("#", 60), "\n")
  
  cat("\nFU2 cohort: N =", nrow(fu2), "\n")
  cat("Events:", sum(fu2$Opioid_Active_At_FU2Date == 1), "\n")
  cat("Non-events:", sum(fu2$Opioid_Active_At_FU2Date == 0), "\n")
  
  n_events <- sum(fu2$Opioid_Active_At_FU2Date == 1)
  cat("\nEvents-per-variable budget:", floor(n_events / 10),
      "predictors (at 10 EPV)\n")
  
  # Primary parsimonious model: ~4-6 parameters
  formula_fu2 <- Opioid_Active_At_FU2Date ~
    AUC10 +
    Smk_Ever +
    Age60plus +
    TxIntensity
  
  result_b <- compare_mle_firth(formula_fu2, fu2, "Model B: FU2 Primary (Parsimonious)")
  
  # Sensitivity: add benzo use
  if ("Benzo_Use" %in% names(fu2)) {
    cat("\n--- Model B2: Sensitivity (+ Benzo Use) ---\n")
    formula_fu2_sens <- Opioid_Active_At_FU2Date ~
      AUC10 +
      Smk_Ever +
      Age60plus +
      TxIntensity +
      Benzo_Use
    
    result_b2 <- compare_mle_firth(formula_fu2_sens, fu2,
                                    "Model B2: FU2 + Benzo Sensitivity")
  }
  
  result_b
}

# =============================================================================
# OR per 10% AUC increase (clinical interpretation helper)
# =============================================================================
interpret_auc_effect <- function(firth_model) {
  beta_auc <- coef(firth_model)["AUC10"]
  or_10pct <- exp(beta_auc)
  
  cat("\n=== Clinical Interpretation ===\n")
  cat("AUC is scaled per 10% increase.\n")
  cat("OR per 10% increase in Pain AUC:", round(or_10pct, 3), "\n")
  cat("Meaning: For every 10% increase in acute pain burden,\n")
  cat("  odds of opioid persistence increase by",
      round((or_10pct - 1) * 100, 1), "%.\n")
  
  # Back-calculate per-1% for abstract
  or_1pct <- exp(beta_auc / 10)
  cat("\nOR per 1% increase:", round(or_1pct, 4), "\n")
  cat("(Same model — just different scaling for reporting.)\n")
}

# =============================================================================
# Main execution
# =============================================================================
cat("=== Firth Penalized Logistic Regression — Primary Inference ===\n")
cat("Using logistf package (Firth / Jeffreys-prior penalization)\n\n")

# --- PLACEHOLDER: Load your prepared data ---
# fu1 <- readRDS(file.path(DATA_DIR, "fu1_analytic.rds"))
# fu2 <- readRDS(file.path(DATA_DIR, "fu2_analytic.rds"))
#
# res_fu1 <- run_fu1_model(fu1)
# res_fu2 <- run_fu2_model(fu2)
#
# interpret_auc_effect(res_fu2$firth)
#
# # Save results
# saveRDS(res_fu1, file.path(OUTPUT_DIR, "fu1_firth_results.rds"))
# saveRDS(res_fu2, file.path(OUTPUT_DIR, "fu2_firth_results.rds"))

# --- Demo with simulated data ---
cat("\n--- Running demo with simulated data ---\n")
set.seed(42)
n <- 200
demo_df <- data.frame(
  Outcome    = rbinom(n, 1, 0.4),
  AUC10      = rnorm(n, mean = 3, sd = 1.5),
  Smk_Ever   = factor(rbinom(n, 1, 0.4), labels = c("Never", "Ever")),
  Age60plus  = factor(rbinom(n, 1, 0.5), labels = c("<60", ">=60")),
  TxIntensity = factor(sample(c("RT_only", "PORT", "Def_ChemoRT"), n,
                               replace = TRUE, prob = c(0.2, 0.3, 0.5)),
                        levels = c("RT_only", "PORT", "Def_ChemoRT"))
)

m_demo <- logistf(Outcome ~ AUC10 + Smk_Ever + Age60plus + TxIntensity,
                   data = demo_df)
format_firth_results(m_demo, "Demo: Simulated FU2-like Model")

cat("\n✓ Demo complete. Uncomment the pipeline above with your real data.\n")
