# Technical Report: Acute-on-Chronic Pain Trajectories and Opioid Persistence in Head and Neck Cancer

**Md Mahin, Oliver Nkuku, Renji He, Clifton Dave Fuller, Amy Moreno, Saba Javed**  
Department of Radiation Oncology, The University of Texas MD Anderson Cancer Center

*Internal Technical Document — April 2026*

---

## 1. Introduction and Scientific Context

Head and neck cancer (HNC) represents a heterogeneous group of malignancies arising from the mucosal surfaces of the oral cavity, oropharynx, hypopharynx, larynx, nasopharynx, and paranasal sinuses. Radiation therapy (RT), often combined with concurrent systemic therapy, remains a cornerstone of curative-intent treatment for locally advanced disease (Pfister et al., 2020). During a standard fractionation course of 6–7 weeks, patients frequently develop severe acute oral mucositis, dysphagia, and pain that significantly impair quality of life (Elting et al., 2007; Trotti et al., 2003).

Pain management during RT for HNC relies heavily on opioid analgesics, yet prescribing practices remain largely non-standardized (Moreno et al., 2023). A substantial proportion of patients initiated on opioids during the acute treatment phase continue to fill opioid prescriptions well beyond the expected resolution of treatment-related toxicity. This phenomenon — termed *opioid persistence* — is clinically distinct from opioid use disorder (OUD) as defined by the DSM-5 (American Psychiatric Association, 2013), but nonetheless carries significant risks including dose escalation, physiologic dependence, and adverse downstream health consequences (Dowell et al., 2022).

The oncology literature has increasingly recognized that the transition from acute to chronic opioid use may be predictable. Several studies have identified pre-treatment risk factors such as tobacco use, psychiatric comorbidity, and prior substance use history (Vitzthum et al., 2020; Patel et al., 2021). However, few studies have leveraged *dynamic, longitudinal symptom trajectories* collected during active treatment as predictors of post-treatment opioid exposure. Patient-reported outcome (PRO) instruments such as the MD Anderson Symptom Inventory for Head and Neck Cancer (MDASI-HN) provide granular, weekly pain assessments that can be summarized quantitatively (Rosenthal et al., 2007; Cleeland et al., 2000).

This report provides a comprehensive technical exposition of a study designed to evaluate whether the **total burden of acute pain during radiation therapy**, quantified as a normalized area under the curve (AUC), independently predicts opioid persistence at subacute (6–15 weeks post-RT) and chronic (21–30 weeks post-RT) follow-up.

---

## 2. Study Population and Data Architecture

Patients were drawn from a prospective clinical registry of HNC patients treated with curative-intent RT at MD Anderson Cancer Center between February 2021 and June 2024. Eligibility was determined by the MOSAIC quality flag (`flag_dose5000to7400_frac20to44_and_4weeklyvisits = 1`), which ensures that patients received a radiobiologically appropriate dose-fractionation scheme (50–74 Gy in 20–44 fractions) and completed at least four weekly symptom visits during treatment. From the initial MOSAIC-eligible pool of **1,118 patients**, analytic cohorts were derived based on the availability of follow-up pain assessments within defined temporal windows.

### 2.1 Demographic and Clinical Characteristics

The FU1 analytic cohort (N = 397) had a median age of 62 years (range 19–88), was 76.5% male, 89.7% White, and 44.1% current or former smokers. The predominant primary sites were oropharynx (40.9%) and oral cavity (17.3%), consistent with the epidemiologic profile of HPV-associated HNC at a tertiary referral center (Chaturvedi et al., 2011). Treatment distributions reflected contemporary patterns: 33.2% received concurrent chemoradiation without surgery (CRT), 24.5% received RT with prior surgery (PORT), 12.1% received RT alone, and the remainder received induction-based regimens (Mahal et al., 2019).

The FU2 cohort (N = 106) was a subset enriched for patients with sufficient longitudinal follow-up. Compared to FU1, this cohort had a higher proportion of male patients (85.0%) and oropharyngeal primaries (57.5%), reflecting the typical selection bias inherent in longitudinal PRO studies where patients with better functional outcomes are more likely to complete extended follow-up (Mercieca-Bebber et al., 2017).

---

## 3. Pain Area Under the Curve: Definition and Computation

Rather than relying on single time-point pain assessments, this study employs the **area under the curve (AUC)** of patient-reported pain scores to capture the cumulative pain burden over time. This approach has precedent in pharmacokinetic modeling (Rowland & Tozer, 2011) and has been adapted for symptom trajectory analysis in oncology (Kahn et al., 2012). The AUC integrates both the intensity and duration of pain, providing a more comprehensive measure than peak or average scores alone.

### 3.1 Trapezoidal Rule

Pain AUC was computed using the trapezoidal rule:

> AUC_raw = Σ [(p(t_i) + p(t_{i+1})) / 2] × (t_{i+1} − t_i)

The raw AUC was then normalized to a percentage:

> AUC_% = [AUC_raw / (10 × (t_last − t_first))] × 100

A score of 50% means the patient experienced, on average, half of the maximum possible pain burden. The MDASI pain item uses a 0–10 numeric rating scale (Cleeland et al., 2000).

### 3.2 Rationale for Continuous Modeling

Pain AUC was retained as a **continuous predictor** in all regression models. Dichotomizing continuous variables discards information and reduces statistical power (Royston et al., 2006; Altman & Royston, 2006). Results are reported per **10% increase** — a linear transformation that does not alter model fit or p-values (Vittinghoff et al., 2012). If the 1% coefficient is β₁, then the 10% coefficient is 10β₁ with SE = 10×SE₁ and identical Z-statistic.

---

## 4. Temporal Phases: Acute, Subacute, and Chronic

The study divides the post-RT trajectory into three biologically motivated phases, consistent with pain chronification models (Glare et al., 2014):

**Acute Phase (During RT):** The standard fractionation course (~6–7 weeks). Weekly MDASI-HN pain scores (WSV1–WSV7) capture the tissue injury response including mucositis, dermatitis, and pharyngeal inflammation (Trotti et al., 2003).

**Subacute Phase (FU1: 6–15 weeks post-RT):** MDASI-HN assessments sent between 42–105 days after RT end date. Captures early recovery. Among 1,118 eligible patients, 794 (71%) were sent an MDASI within this window; 397 (35.5%) met all analytic criteria (median 66 days post-RT).

**Chronic Phase (FU2: 21–30 weeks post-RT):** MDASI-HN assessments sent between 150–210 days after RT end date. Represents early survivorship. Only 106 patients (9.5%) met all FU2 criteria (median 180 days post-RT), reflecting substantial longitudinal PRO attrition (Bottomley et al., 2016).

These definitions were established *a priori* based on clinical consensus and are consistent with NCI-CTCAE late toxicity definitions (U.S. DHHS, 2017).

---

## 5. Outcome Definition: Opioid Persistence

The primary outcome was **active opioid prescription within ±7 days of the completed follow-up pain assessment**, anchoring the opioid exposure assessment to the same clinical encounter at which pain was evaluated.

### 5.1 Terminological Precision

The DSM-5 criteria for OUD require evidence of impaired control, social impairment, risky use, and pharmacological indicators assessed over 12 months (APA, 2013). Our endpoint measures **opioid persistence** — ongoing prescribing behavior — not addiction, misuse, or functional impairment. This distinction prevents overclaiming and shapes clinical implications: interventions to reduce persistence (e.g., earlier non-opioid transition) differ from OUD interventions (e.g., buprenorphine, behavioral therapy) (Dowell et al., 2022; Volkow & McLellan, 2016).

---

## 6. Multivariable Logistic Regression Framework

Logistic regression modeled the probability of opioid persistence:

> logit(P(Y = 1 | X)) = X^T β

where Y is the binary outcome and β the coefficient vector estimated by maximum likelihood (McCullagh & Nelder, 1989).

### 6.1 Predictor Selection

Covariates were selected on clinical relevance, not statistical screening: Pain AUC (continuous, per 10%), smoking (current/former vs. never), age, gender, race, ethnicity, primary cancer site, treatment regimen, and benzodiazepine use. Each was included based on prior literature (Vitzthum et al., 2020; Brummett et al., 2017).

### 6.2 Events-Per-Variable Considerations

The ~10 EPV rule (Peduzzi et al., 1996; Vittinghoff & McCulloch, 2007) constrains model complexity. FU1: 287 events / 26 parameters = EPV ≈ 11 (acceptable). FU2: 47 events / 23 parameters = EPV ≈ 2 (unacceptable → parsimonious model required).

### 6.3 Odds Ratio Interpretation

OR of 1.05 per 1% means OR₁₀% = (1.05)¹⁰ ≈ 1.63 — a 10% increase in acute pain AUC corresponds to 63% higher odds of opioid persistence. For context, established surgical risk factors typically show ORs of 1.5–3.0 (Brummett et al., 2017; Sun et al., 2016). The Pain AUC threshold-level effects (OR 7–8) are among the strongest reported predictors in oncology.

---

## 7. Threshold Search for Clinical Risk Stratification

An evenly spaced grid search between the 5th–95th percentiles identified clinically actionable cutpoints. At each threshold, the cohort was dichotomized and a multivariable model fitted.

**FU1:** Lowest significant threshold AUC ≥ 5% (OR 2.6, p = 0.02); strongest at 53.2% (OR 7.8). **FU2:** Lowest at 15.8% (OR 4.3, p = 0.01); strongest at 18.3% (OR 7.0).

The divergent thresholds are biologically coherent: very high acute burden predicts early persistence; even moderate burden predicts late persistence — consistent with postoperative opioid literature (Brummett et al., 2017; Sun et al., 2016).

Grid search carries data-dredging risk (Lausen & Schumacher, 1992). Mitigations include: framing as exploratory, bootstrap validation, optional Bayesian optimization with GP surrogate (Snoek et al., 2012), and retaining the continuous model as primary inference. Binary thresholding is inherently information-lossy (Royston et al., 2006).

---

## 8. Data Imputation Strategy

Missing weekly pain scores were imputed using linear interpolation between nearest available assessments. Among 628 patients with trajectory data, 508 (80.9%) had ≥1 imputed timepoint; Week 7 was imputed in 51.3% (the most frequent). FU1 was never imputed (required anchor).

The high rate means the trapezoidal AUC is partly synthetic. Sensitivity analyses comparing imputed vs. complete-case results are recommended. If effect direction persists, the imputation is unlikely to have introduced substantial bias (Sterne et al., 2009). More sophisticated approaches such as MICE could be considered (van Buuren & Groothuis-Oudshoorn, 2011).

---

## 9. Cohort Attrition and Statistical Power

| Stage | N | % of Eligible |
|-------|---|---------------|
| MOSAIC-eligible cohort | 1,118 | 100% |
| FU1 MDASI sent (42–105 days) | 794 | 71.0% |
| FU1 analytic (all criteria) | 397 | **35.5%** |
| FU2 MDASI sent (150–210 days) | 281 | 25.1% |
| FU2 analytic (all criteria) | 106 | **9.5%** |

This 90.5% attrition from eligible to FU2 is not unusual for longitudinal PRO studies (Mercieca-Bebber et al., 2017; Bottomley et al., 2016), but has three consequences: **(1) Reduced power** — 47 events support ~4–5 predictors at 10 EPV (Peduzzi et al., 1996); **(2) Potential selection bias** — completers may differ systematically; **(3) Subgroup instability** — real effects may appear non-significant, and spurious ones may appear significant (Button et al., 2013).

---

## 10. Key Findings and Clinical Interpretation

### 10.1 Pain AUC: The Dominant Predictor

AUC_pain-acute was the strongest predictor at both timepoints. FU1: OR 1.05 per 1% (95% CI 1.03–1.06, p < 0.001); FU2: OR 1.05 per 1% (95% CI 1.02–1.09, p = 0.006). A 10% increase ≈ 63% higher odds. The biological plausibility rests on central sensitization — prolonged nociceptive input induces maladaptive neuroplastic changes (Woolf, 2011; Latremoliere & Woolf, 2009).

### 10.2 Smoking: Temporal Amplification

FU1: OR 2.15 (95% CI 1.25–3.71, p = 0.006). FU2 (unstabilized): OR 5.27 (95% CI 1.54–18.03). Tobacco use is a known predictor of persistent pain and slower healing (Shi et al., 2010; Ditre et al., 2011). The temporal amplification suggests smokers may transition to chronic persistence rather than tapering.

### 10.3 Age: Late-Phase Protection

FU1: non-significant (OR 1.07). FU2: age ≥60 → OR 0.28 (95% CI 0.09–0.89, p = 0.031). Possible explanations include more cautious taper protocols for older patients or differential opioid metabolism (Chau et al., 2008).

### 10.4 Race and Tumor Site

FU1: Asian race OR 0.16 (p = 0.004); larynx/hypopharynx OR 0.37 (p = 0.03). Not reproduced in FU2 (reduced power, more homogeneous cohort). Racial disparities in opioid prescribing are well-documented (Meghani et al., 2012).

---

## 11. Model Instability, Separation, and Penalized Regression

### 11.1 What Is Separation?

Separation occurs when predictor values perfectly predict the outcome. The MLE attempts P(Y=1) = 1, requiring X^T β → +∞, so β → ±∞ (Albert & Anderson, 1984; Heinze & Schemper, 2002). The likelihood has a supremum but no finite maximum.

### 11.2 Sources in This Study

- Very small race groups (Black n = 1 in FU2)
- Very small site groups (Eye & Orbit n = 1)
- Small treatment subcategories (ICRT + Surgery Yes n = 1)
- 26 predictors with 110 non-events (FU1); 23 predictors with 59 non-events (FU2)

Pain AUC is **not** affected — it is continuous, broadly distributed (0–65%+), and shows outcome overlap.

### 11.3 Firth Penalized Logistic Regression

Firth (1993) adds a Jeffreys invariant prior:

> l*(β) = l(β) + ½ log|I(β)|

Benefits: finite estimates under separation, O(n⁻¹) bias reduction, valid likelihood-based inference. Implemented in R via `logistf` (Heinze et al., 2024), in Python via custom Newton-Raphson. Widely accepted for rare-event clinical research (King & Zeng, 2001). Unlike L2, requires no tuning hyperparameter and preserves profile likelihood ratio tests.

### 11.4 FU2 Parsimonious Model

| Variable | Coding | Parameters |
|----------|--------|------------|
| Pain AUC (acute) | Continuous, per 10% | 1 |
| Smoking | Ever vs. Never | 1 |
| Age | ≥60 vs. <60 | 1 |
| Treatment intensity | 2–3 levels | 1–2 |
| **Total** | | **4–5** |

Race, detailed site, and detailed regimen excluded from inference due to sparse strata. L2 sensitivity (Friedman et al., 2010) confirmed consistent directionality.

---

## 12. Limitations

1. **Attrition**: 90% drop from eligible to FU2 may introduce selection bias (Mercieca-Bebber et al., 2017)
2. **Imputation**: 80.9% with ≥1 imputed value; AUC partly synthetic (Little & Rubin, 2019)
3. **Missing confounders**: No prior opioid use, psychiatric history, or substance use data (Vitzthum et al., 2020)
4. **Outcome scope**: Measures persistence, not OUD or misuse
5. **FU2 power**: 47 events limits subgroup inference
6. **Single institution**: Generalizability unknown; external validation needed

---

## 13. References

See [`references.md`](references.md) for the complete bibliography (46 cited works).
