# =============================================================================
# 01_data_preparation.R
# Pain AUC & Opioid Persistence — Data Preparation
# ASTRO 2026 | Mahin, Nkuku, He, Fuller, Moreno, Javed
# =============================================================================
#
# This script:
#   1. Loads the analytic cohort (MOSAIC-eligible HNC patients)
#   2. Defines FU1 and FU2 cohorts based on time windows
#   3. Computes normalized Pain AUC (%) via trapezoidal rule
#   4. Collapses sparse categories for stable modeling
#   5. Defines the opioid persistence outcome
#
# INPUT:  Raw analytic dataset (one row per patient)
# OUTPUT: fu1_analytic.rds, fu2_analytic.rds — ready for modeling
# =============================================================================

library(dplyr)
library(tidyr)

# ─── Configuration ───────────────────────────────────────────────────────────
# Update these paths to match your environment
DATA_PATH   <- "data/raw_cohort.csv"   # Your raw data file
OUTPUT_DIR  <- "data/processed/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Time windows (days after RT end date)
FU1_WINDOW <- c(42, 105)   # 6–15 weeks
FU2_WINDOW <- c(150, 210)  # 5–7 months

# Minimum weekly visits required (first + last + at least 1-2 in between)
MIN_WEEKLY_VISITS <- 4  # WSV1 + (WSV6 or WSV7) + >=1 WSV2-5 + FUP

# ─── Helper: Trapezoidal AUC ────────────────────────────────────────────────
compute_pain_auc <- function(pain_scores, timepoints) {
  # pain_scores: numeric vector of pain values (0–10 scale)
  # timepoints:  numeric vector of corresponding time indices (e.g., week numbers)
  #
  # Returns: list(raw_auc, normalized_auc_pct)
  
  valid <- !is.na(pain_scores) & !is.na(timepoints)
  ps <- pain_scores[valid]
  tp <- timepoints[valid]
  
  if (length(ps) < 2) return(list(raw_auc = NA, auc_pct = NA))
  
  ord <- order(tp)
  ps <- ps[ord]
  tp <- tp[ord]
  
  # Trapezoidal rule: sum of (avg of consecutive pain) * (time difference)
  raw_auc <- sum(
    ((ps[-length(ps)] + ps[-1]) / 2) * diff(tp)
  )
  
  # Normalize: max possible AUC = 10 * (last_timepoint - first_timepoint)
  max_auc <- 10 * (tp[length(tp)] - tp[1])
  auc_pct <- ifelse(max_auc > 0, (raw_auc / max_auc) * 100, NA)
  
  list(raw_auc = raw_auc, auc_pct = auc_pct)
}

# ─── Helper: Impute missing weeks ───────────────────────────────────────────
impute_pain <- function(pain_vec) {
  # Linear interpolation between nearest non-NA neighbors
  # pain_vec: named numeric vector with names = week identifiers
  #
  # Logic from study protocol:
  #   If W7 missing: W7 = (W6 + FU1) / 2
  #   If W6 also missing: W7 = (W5 + FU1) / 2
  #   For missing W1: half of W2 (if W2 exists)
  #   For missing FU3: half of FU2 (if FU2 exists)
  
  n <- length(pain_vec)
  imputed <- pain_vec
  was_imputed <- rep(FALSE, n)
  
  for (i in seq_along(imputed)) {
    if (is.na(imputed[i])) {
      # Find nearest non-NA left and right
      left_val  <- NA
      right_val <- NA
      
      if (i > 1) {
        left_idx <- max(which(!is.na(imputed[1:(i-1)])), na.rm = FALSE)
        if (length(left_idx) > 0 && !is.na(left_idx)) left_val <- imputed[left_idx]
      }
      if (i < n) {
        right_candidates <- which(!is.na(imputed[(i+1):n]))
        if (length(right_candidates) > 0) {
          right_idx <- min(right_candidates) + i
          right_val <- imputed[right_idx]
        }
      }
      
      if (!is.na(left_val) && !is.na(right_val)) {
        imputed[i] <- (left_val + right_val) / 2
        was_imputed[i] <- TRUE
      } else if (!is.na(left_val)) {
        imputed[i] <- left_val / 2
        was_imputed[i] <- TRUE
      } else if (!is.na(right_val)) {
        imputed[i] <- right_val / 2
        was_imputed[i] <- TRUE
      }
    }
  }
  
  list(values = imputed, was_imputed = was_imputed)
}

# ─── Helper: Collapse categories ────────────────────────────────────────────
collapse_smoking <- function(smoking_raw) {
  # Never vs. Ever (current + former)
  ifelse(smoking_raw %in% c("Current", "Former", "Current/Former"),
         "Ever", "Never")
}

collapse_site <- function(site_raw, n_levels = 2) {
  if (n_levels == 2) {
    # Oropharynx vs. Other
    ifelse(site_raw == "Oropharynx", "Oropharynx", "Other")
  } else {
    # 3-level: Oropharynx / Oral Cavity / Everything else
    case_when(
      site_raw == "Oropharynx"        ~ "Oropharynx",
      site_raw == "Oral Cavity/Pharynx" ~ "OralCavity",
      TRUE                             ~ "Other"
    )
  }
}

collapse_treatment <- function(rt_only, surgery, concurrent, induction, n_levels = 3) {
  # 3-level treatment intensity:
  #   1. RT_only:      RT alone, no surgery/chemo/induction
  #   2. PORT:         Surgery + RT (+/- chemo)  [post-operative RT]
  #   3. Def_ChemoRT:  No surgery, any chemo (CRT/ICRT/IRT) [definitive chemoRT]
  
  tx <- case_when(
    rt_only == 1                                  ~ "RT_only",
    surgery == 1                                  ~ "PORT",
    (concurrent == 1 | induction == 1)            ~ "Def_ChemoRT",
    TRUE                                          ~ "RT_only"  # fallback
  )
  
  if (n_levels == 2) {
    # Collapse to: ChemoRT vs. RT_only_or_PORT
    tx <- ifelse(tx == "Def_ChemoRT", "ChemoRT", "RT_only_or_PORT")
  }
  
  factor(tx, levels = if (n_levels == 3) {
    c("RT_only", "PORT", "Def_ChemoRT")
  } else {
    c("RT_only_or_PORT", "ChemoRT")
  })
}

collapse_race <- function(race_raw, n_levels = 2) {
  if (n_levels == 2) {
    ifelse(race_raw == "White", "White", "Non_White")
  } else {
    # Keep White, collapse rest
    case_when(
      race_raw == "White" ~ "White",
      race_raw == "Black" ~ "Black",
      race_raw == "Asian" ~ "Asian",
      TRUE                ~ "Other"
    )
  }
}

# ─── Opioid persistence definition ──────────────────────────────────────────
define_opioid_outcome <- function(opioid_dates, fu_date, window_days = 7) {
  # Active opioid prescription within +/- window_days of the FU pain assessment
  # opioid_dates: vector of dates with active Rx
  # fu_date:      the patient's FU visit date
  # window_days:  tolerance window (default +/- 7 days)
  
  if (is.na(fu_date) || length(opioid_dates) == 0) return(0L)
  
  any_active <- any(
    abs(as.numeric(difftime(opioid_dates, fu_date, units = "days"))) <= window_days,
    na.rm = TRUE
  )
  
  as.integer(any_active)
}

# ─── Sparse factor safety check ─────────────────────────────────────────────
check_sparse_factors <- function(df, outcome_col, factor_cols, min_cell = 5) {
  # For each categorical variable, check if any level has <min_cell total

  # or 0 events / 0 non-events. Returns a tibble of warnings.
  
  warnings <- list()
  
  for (col in factor_cols) {
    tab <- table(df[[col]], df[[outcome_col]])
    
    for (level in rownames(tab)) {
      n_total  <- sum(tab[level, ])
      n_events <- tab[level, "1"]  # assuming 0/1 coding
      n_nonevt <- tab[level, "0"]
      
      if (n_total < min_cell || n_events == 0 || n_nonevt == 0) {
        warnings <- c(warnings, list(tibble(
          variable = col,
          level    = level,
          n_total  = n_total,
          n_events = n_events,
          n_nonevents = n_nonevt,
          issue    = case_when(
            n_events == 0   ~ "zero events (separation risk)",
            n_nonevt == 0   ~ "zero non-events (separation risk)",
            n_total < min_cell ~ paste0("sparse (<", min_cell, " total)")
          )
        )))
      }
    }
  }
  
  if (length(warnings) > 0) bind_rows(warnings) else tibble()
}

# ─── Main pipeline ──────────────────────────────────────────────────────────
cat("=== Pain AUC & Opioid Persistence — Data Preparation ===\n\n")

# --- PLACEHOLDER: Load your data here ---
# df <- read.csv(DATA_PATH, stringsAsFactors = FALSE)
#
# Expected columns (adapt to your naming):
#   MRN, Age, Gender, Race, Ethnicity, Smoking,
#   PrimarySite, Surgery, Concurrent, Induction, RT_only,
#   RT_EndDate, FU1_Date, FU2_Date,
#   Pain_W1..Pain_W7, Pain_FU1, Pain_FU2,
#   MOSAIC_flag, Opioid_Active_At_FU1Date, Opioid_Active_At_FU2Date,
#   Benzo_Use

cat(">> Load your data and uncomment the pipeline below.\n")
cat(">> The helper functions above are ready to use.\n\n")

# --- Example pipeline (uncomment and adapt) ---
#
# # 1. Filter MOSAIC-eligible
# df <- df %>% filter(MOSAIC_flag == 1)
# cat("MOSAIC-eligible:", nrow(df), "patients\n")
#
# # 2. Compute Pain AUC for each patient
# # ... (apply compute_pain_auc to each row's weekly pain + FU pain)
#
# # 3. Build FU1 cohort
# fu1 <- df %>%
#   filter(!is.na(Pain_FU1),
#          FU1_days_from_RT_end >= FU1_WINDOW[1],
#          FU1_days_from_RT_end <= FU1_WINDOW[2]) %>%
#   mutate(
#     AUC10        = AUCpain_acute_pct / 10,
#     Smk_Ever     = collapse_smoking(Smoking),
#     Age60plus    = as.integer(Age >= 60),
#     Site2        = collapse_site(PrimarySite, n_levels = 2),
#     TxIntensity  = collapse_treatment(RT_only, Surgery, Concurrent, Induction, n_levels = 3),
#     Race2        = collapse_race(Race, n_levels = 2)
#   )
#
# # 4. Build FU2 cohort
# fu2 <- df %>%
#   filter(!is.na(Pain_FU2),
#          FU2_days_from_RT_end >= FU2_WINDOW[1],
#          FU2_days_from_RT_end <= FU2_WINDOW[2]) %>%
#   mutate(
#     AUC10        = AUCpain_acute_pct / 10,
#     Smk_Ever     = collapse_smoking(Smoking),
#     Age60plus    = as.integer(Age >= 60),
#     TxIntensity  = collapse_treatment(RT_only, Surgery, Concurrent, Induction, n_levels = 3)
#   )
#
# # 5. Sparse factor check
# cat("\n--- FU2 Sparse Factor Check ---\n")
# sparse_warnings <- check_sparse_factors(
#   fu2, "Opioid_Active_At_FU2Date",
#   c("Smk_Ever", "TxIntensity")
# )
# if (nrow(sparse_warnings) > 0) print(sparse_warnings)
#
# # 6. Save
# saveRDS(fu1, file.path(OUTPUT_DIR, "fu1_analytic.rds"))
# saveRDS(fu2, file.path(OUTPUT_DIR, "fu2_analytic.rds"))
# cat("\nSaved analytic datasets to", OUTPUT_DIR, "\n")

cat("Done.\n")
