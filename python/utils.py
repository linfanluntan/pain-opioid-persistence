"""
utils.py
Utility functions for Pain AUC computation and variable collapsing.
ASTRO 2026 | Mahin, Nkuku, He, Fuller, Moreno, Javed
MD Anderson Cancer Center

KEY FUNCTIONS:
    compute_pain_auc()     - Trapezoidal rule AUC with normalization
                             [Rowland & Tozer, 2011; Kahn et al., 2012]
    impute_pain_scores()   - Linear interpolation for missing weekly PROs
                             [See docs/technical_report.md §8 for caveats]
    collapse_smoking()     - Ever vs. Never (prevents sparse cells in FU2)
    collapse_site()        - Oropharynx vs. Other (prevents n=1 separation)
    collapse_treatment()   - RT alone / PORT / Def ChemoRT intensity levels
    collapse_race()        - White vs. Non-White (for FU2 sparse strata)
    check_sparse_factors() - Identify cells with <5 patients or 0 events
                             that cause quasi-complete separation
                             [Albert & Anderson, 1984]
    build_design_matrix()  - Numeric matrix with dummy-coded categoricals

COLLAPSING RATIONALE (docs/statistical_notes.md §5):
    FU2 has 47 events → ~4-5 predictors at 10 EPV [Peduzzi et al., 1996].
    Detailed race (4 levels), site (8-10 levels), and regimen (8 levels)
    create sparse cells that cause separation. Collapsing eliminates
    infinite ORs while preserving the primary Pain AUC signal.

See docs/references.md for the complete bibliography.
"""

import numpy as np
import pandas as pd
from typing import List, Optional, Tuple


# ─── Pain AUC Computation ───────────────────────────────────────────────────

def compute_pain_auc(
    pain_scores: np.ndarray,
    timepoints: np.ndarray,
) -> Tuple[float, float]:
    """
    Compute normalized Pain AUC (%) via trapezoidal rule.

    Parameters
    ----------
    pain_scores : array-like
        Pain values on 0–10 scale at each timepoint.
    timepoints : array-like
        Corresponding time indices (e.g., week numbers: 1,2,...,7,8 for FU1).

    Returns
    -------
    raw_auc : float
        Raw trapezoidal AUC.
    auc_pct : float
        Normalized AUC as percentage (0–100).
        Formula: AUC_raw / (10 * (last_timepoint - first_timepoint)) * 100
    """
    pain_scores = np.asarray(pain_scores, dtype=float)
    timepoints = np.asarray(timepoints, dtype=float)

    # Remove NaN pairs
    valid = ~(np.isnan(pain_scores) | np.isnan(timepoints))
    ps = pain_scores[valid]
    tp = timepoints[valid]

    if len(ps) < 2:
        return np.nan, np.nan

    # Sort by timepoint
    order = np.argsort(tp)
    ps = ps[order]
    tp = tp[order]

    # Trapezoidal rule
    raw_auc = np.sum(((ps[:-1] + ps[1:]) / 2.0) * np.diff(tp))

    # Normalize
    max_auc = 10.0 * (tp[-1] - tp[0])
    auc_pct = (raw_auc / max_auc * 100.0) if max_auc > 0 else np.nan

    return raw_auc, auc_pct


def impute_pain_scores(pain_vec: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Impute missing pain scores using linear interpolation between neighbors.

    Protocol:
        If W7 missing → W7 = (W6 + FU1) / 2
        If W6 also missing → W7 = (W5 + FU1) / 2
        For endpoints: half of nearest neighbor

    Parameters
    ----------
    pain_vec : array-like
        Pain scores indexed by time (may contain NaN).

    Returns
    -------
    imputed : ndarray
        Pain scores with NaN values filled.
    was_imputed : ndarray of bool
        True for positions that were imputed.
    """
    pain_vec = np.asarray(pain_vec, dtype=float)
    imputed = pain_vec.copy()
    was_imputed = np.zeros(len(pain_vec), dtype=bool)

    for i in range(len(imputed)):
        if np.isnan(imputed[i]):
            left_val = np.nan
            right_val = np.nan

            # Find nearest non-NaN to the left
            for j in range(i - 1, -1, -1):
                if not np.isnan(imputed[j]):
                    left_val = imputed[j]
                    break

            # Find nearest non-NaN to the right
            for j in range(i + 1, len(imputed)):
                if not np.isnan(imputed[j]):
                    right_val = imputed[j]
                    break

            if not np.isnan(left_val) and not np.isnan(right_val):
                imputed[i] = (left_val + right_val) / 2.0
            elif not np.isnan(left_val):
                imputed[i] = left_val / 2.0
            elif not np.isnan(right_val):
                imputed[i] = right_val / 2.0

            if not np.isnan(imputed[i]):
                was_imputed[i] = True

    return imputed, was_imputed


# ─── Variable Collapsing ────────────────────────────────────────────────────

def collapse_smoking(smoking_series: pd.Series) -> pd.Series:
    """Collapse smoking to Never vs. Ever (current + former)."""
    return smoking_series.map(
        lambda x: "Ever" if x in ("Current", "Former", "Current/Former") else "Never"
    )


def collapse_site(site_series: pd.Series, n_levels: int = 2) -> pd.Series:
    """
    Collapse primary cancer site.
    n_levels=2: Oropharynx vs. Other
    n_levels=3: Oropharynx / Oral Cavity / Other
    """
    if n_levels == 2:
        return site_series.map(
            lambda x: "Oropharynx" if x == "Oropharynx" else "Other"
        )
    else:
        def _map(x):
            if x == "Oropharynx":
                return "Oropharynx"
            elif x in ("Oral Cavity/Pharynx", "Oral Cavity"):
                return "OralCavity"
            else:
                return "Other"
        return site_series.map(_map)


def collapse_treatment(
    df: pd.DataFrame,
    rt_only_col: str = "RT_only",
    surgery_col: str = "Surgery",
    concurrent_col: str = "Concurrent",
    induction_col: str = "Induction",
    n_levels: int = 3,
) -> pd.Series:
    """
    Collapse treatment regimen to intensity levels.

    3-level:
        RT_only:     RT alone
        PORT:        Surgery + RT (± chemo)
        Def_ChemoRT: No surgery, any chemo

    2-level:
        ChemoRT vs. RT_only_or_PORT
    """
    def _assign(row):
        if row[rt_only_col] == 1:
            return "RT_only"
        elif row[surgery_col] == 1:
            return "PORT"
        elif row[concurrent_col] == 1 or row[induction_col] == 1:
            return "Def_ChemoRT"
        else:
            return "RT_only"

    tx = df.apply(_assign, axis=1)

    if n_levels == 2:
        tx = tx.map(lambda x: "ChemoRT" if x == "Def_ChemoRT" else "RT_only_or_PORT")

    return pd.Categorical(tx)


def collapse_race(race_series: pd.Series, n_levels: int = 2) -> pd.Series:
    """Collapse race. n_levels=2: White vs. Non-White."""
    if n_levels == 2:
        return race_series.map(lambda x: "White" if x == "White" else "Non_White")
    else:
        def _map(x):
            if x in ("White", "Black", "Asian"):
                return x
            else:
                return "Other"
        return race_series.map(_map)


# ─── Sparse Factor Check ────────────────────────────────────────────────────

def check_sparse_factors(
    df: pd.DataFrame,
    outcome_col: str,
    factor_cols: List[str],
    min_cell: int = 5,
) -> pd.DataFrame:
    """
    Check for sparse categories that could cause separation.

    Returns a DataFrame of warnings for any factor level with:
        - Total n < min_cell
        - Zero events
        - Zero non-events
    """
    warnings = []

    for col in factor_cols:
        ct = pd.crosstab(df[col], df[outcome_col])

        for level in ct.index:
            n_total = ct.loc[level].sum()
            n_events = ct.loc[level].get(1, 0)
            n_nonevents = ct.loc[level].get(0, 0)

            issue = None
            if n_events == 0:
                issue = "zero events (separation risk)"
            elif n_nonevents == 0:
                issue = "zero non-events (separation risk)"
            elif n_total < min_cell:
                issue = f"sparse (<{min_cell} total)"

            if issue:
                warnings.append({
                    "variable": col,
                    "level": level,
                    "n_total": n_total,
                    "n_events": n_events,
                    "n_nonevents": n_nonevents,
                    "issue": issue,
                })

    return pd.DataFrame(warnings) if warnings else pd.DataFrame()


# ─── Design Matrix Builder ──────────────────────────────────────────────────

def build_design_matrix(
    df: pd.DataFrame,
    continuous_cols: List[str],
    categorical_cols: List[str],
    add_intercept: bool = True,
) -> Tuple[np.ndarray, List[str]]:
    """
    Build a numeric design matrix with dummy-coded categoricals.

    Returns
    -------
    X : ndarray of shape (n, p)
    term_names : list of str
    """
    parts = []
    names = []

    if add_intercept:
        parts.append(np.ones((len(df), 1)))
        names.append("Intercept")

    for col in continuous_cols:
        parts.append(df[col].values.reshape(-1, 1).astype(float))
        names.append(col)

    for col in categorical_cols:
        dummies = pd.get_dummies(df[col], prefix=col, drop_first=True, dtype=float)
        parts.append(dummies.values)
        names.extend(dummies.columns.tolist())

    X = np.hstack(parts)
    return X, names


# ─── Demo ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=== Utils Demo ===\n")

    # Pain AUC computation
    pain = np.array([2, 4, 5, 6, 7, 5, 3, 1])  # W1-W7 + FU1
    weeks = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    raw, pct = compute_pain_auc(pain, weeks)
    print(f"Pain AUC: raw={raw:.1f}, normalized={pct:.1f}%")

    # Imputation
    pain_missing = np.array([2, 4, np.nan, 6, np.nan, 5, np.nan, 1])
    imputed, flags = impute_pain_scores(pain_missing)
    print(f"\nOriginal:  {pain_missing}")
    print(f"Imputed:   {imputed}")
    print(f"Flags:     {flags}")

    # Design matrix
    demo = pd.DataFrame({
        "AUC10": [3.2, 1.5, 4.8, 2.1],
        "Smoking": ["Never", "Ever", "Ever", "Never"],
        "TxIntensity": ["RT_only", "PORT", "Def_ChemoRT", "PORT"],
    })
    X, names = build_design_matrix(demo, ["AUC10"], ["Smoking", "TxIntensity"])
    print(f"\nDesign matrix shape: {X.shape}")
    print(f"Terms: {names}")

    print("\n✓ Utils demo complete.")
