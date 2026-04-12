"""
firth_logistic.py
Pure NumPy implementation of Firth (Jeffreys-prior) penalized logistic regression.

Firth regression modifies the score function to prevent coefficient divergence
under separation and correct small-sample bias:

    l*(β) = l(β) + ½ log|I(β)|

where I(β) is the Fisher information matrix. This is equivalent to adding a
Jeffreys invariant prior to the likelihood, producing finite estimates even
under complete or quasi-complete separation.

THEORETICAL BACKGROUND (see docs/technical_report.md §11 and docs/statistical_notes.md §3-4):
    - Separation occurs when predictor values perfectly predict the outcome,
      causing the MLE to diverge to ±∞ (Albert & Anderson, Biometrika 1984)
    - Firth's penalty guarantees finite estimates and reduces O(n⁻¹) bias
      (Firth, Biometrika 1993; Heinze & Schemper, Stat Med 2002)
    - Unlike L2/ridge regularization, Firth requires no tuning hyperparameter
      and preserves likelihood-based inference (valid CIs, profile LR tests)
    - Standard for rare-event logistic regression in clinical research
      (King & Zeng, Polit Anal 2001)

WHY THIS IMPLEMENTATION EXISTS:
    In our FU2 cohort (N=106, 47 events, 23+ parameters), standard MLE
    produces ORs of 47,924 and 7,532,188 with SEs > 2,000. Firth regression
    resolves this by adding curvature to the likelihood surface via the
    Jeffreys prior. The primary Pain AUC predictor (continuous, broadly
    distributed) is unaffected by separation — it survives penalization
    because it has no sparse cells and overlaps between outcome groups.

R EQUIVALENT: logistf::logistf() (Heinze, Ploner & Jiricka, 2024)

References:
    Albert A, Anderson JA. Biometrika. 1984;71(1):1-10.
    Firth D. Biometrika. 1993;80(1):27-38.
    Heinze G, Schemper M. Stat Med. 2002;21(16):2409-2419.
    King G, Zeng L. Polit Anal. 2001;9(2):137-163.
    See docs/references.md for the complete bibliography.

Usage:
    from firth_logistic import FirthLogisticRegression

    model = FirthLogisticRegression()
    model.fit(X, y)
    print(model.summary())
"""

import numpy as np
from typing import Optional, Tuple, Dict


def _sigmoid(z: np.ndarray) -> np.ndarray:
    """Numerically stable sigmoid function."""
    z = np.clip(z, -500, 500)
    return 1.0 / (1.0 + np.exp(-z))


def firth_logistic(
    X: np.ndarray,
    y: np.ndarray,
    max_iter: int = 100,
    tol: float = 1e-8,
    ridge: float = 1e-10,
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Fit Firth penalized logistic regression via modified Newton-Raphson (IRLS).

    Parameters
    ----------
    X : ndarray of shape (n, p)
        Design matrix. Should include an intercept column if desired.
    y : ndarray of shape (n,)
        Binary outcome vector (0/1).
    max_iter : int
        Maximum Newton iterations.
    tol : float
        Convergence tolerance on max |Δβ|.
    ridge : float
        Tiny diagonal jitter for numerical stability of matrix inversion.

    Returns
    -------
    beta : ndarray of shape (p,)
        Firth-penalized coefficient estimates.
    cov : ndarray of shape (p, p)
        Approximate covariance matrix (inverse Fisher information at solution).
    info : dict
        Convergence diagnostics: iterations, converged, penalized_loglik.
    """
    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64).ravel()
    n, p = X.shape

    beta = np.zeros(p)
    converged = False

    for it in range(1, max_iter + 1):
        eta = X @ beta
        mu = _sigmoid(eta)

        # Weights: W = diag(mu * (1 - mu))
        w = mu * (1.0 - mu)
        w = np.maximum(w, 1e-12)  # prevent zero weights

        # Fisher information: I = X' W X
        I = X.T @ (X * w[:, None]) + ridge * np.eye(p)

        # Inverse information
        I_inv = np.linalg.inv(I)

        # Hat matrix diagonal: h_i = diag(W^{1/2} X I^{-1} X' W^{1/2})
        sqrt_w = np.sqrt(w)
        A = X * sqrt_w[:, None]        # (n, p)
        AIinv = A @ I_inv              # (n, p)
        h = np.sum(AIinv * A, axis=1)  # (n,)

        # Firth adjustment: a = ½ (1 - 2μ) h
        a = 0.5 * (1.0 - 2.0 * mu) * h

        # Modified score: U* = X' (y - μ + a)
        U_star = X.T @ (y - mu + a)

        # Newton step: β_new = β + I^{-1} U*
        step = I_inv @ U_star
        beta_new = beta + step

        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            converged = True
            break

        beta = beta_new

    # Final covariance at solution
    eta = X @ beta
    mu = _sigmoid(eta)
    w = mu * (1.0 - mu)
    w = np.maximum(w, 1e-12)
    I = X.T @ (X * w[:, None]) + ridge * np.eye(p)
    cov = np.linalg.inv(I)

    # Penalized log-likelihood
    loglik = np.sum(y * np.log(mu + 1e-15) + (1 - y) * np.log(1 - mu + 1e-15))
    pen_loglik = loglik + 0.5 * np.log(np.linalg.det(I))

    info = {
        "iterations": it,
        "converged": converged,
        "loglik": loglik,
        "penalized_loglik": pen_loglik,
    }

    return beta, cov, info


def odds_ratios_ci(
    beta: np.ndarray,
    cov: np.ndarray,
    alpha: float = 0.05,
    term_names: Optional[list] = None,
) -> Dict:
    """
    Compute odds ratios and Wald confidence intervals from Firth estimates.

    Parameters
    ----------
    beta : ndarray of shape (p,)
        Coefficient estimates.
    cov : ndarray of shape (p, p)
        Covariance matrix.
    alpha : float
        Significance level (default 0.05 for 95% CI).
    term_names : list of str, optional
        Names for each coefficient.

    Returns
    -------
    dict with keys: OR, CI_low, CI_high, SE, Z, p_value, terms
    """
    se = np.sqrt(np.diag(cov))
    z_crit = 1.959963984540054  # norm.ppf(1 - alpha/2)
    z_stat = beta / se
    p_values = 2.0 * (1.0 - _norm_cdf(np.abs(z_stat)))

    if term_names is None:
        term_names = [f"X{i}" for i in range(len(beta))]

    return {
        "terms": term_names,
        "beta": beta,
        "SE": se,
        "Z": z_stat,
        "p_value": p_values,
        "OR": np.exp(beta),
        "CI_low": np.exp(beta - z_crit * se),
        "CI_high": np.exp(beta + z_crit * se),
    }


def _norm_cdf(x):
    """Standard normal CDF using error function."""
    return 0.5 * (1.0 + _erf(x / np.sqrt(2.0)))


def _erf(x):
    """Approximation of error function (Abramowitz & Stegun)."""
    # Use numpy's built-in if available
    try:
        return np.vectorize(lambda z: np.math.erf(z))(x)
    except AttributeError:
        from scipy.special import erf
        return erf(x)


# ─── Convenience class ──────────────────────────────────────────────────────
class FirthLogisticRegression:
    """
    Firth penalized logistic regression with scikit-learn-like API.

    Example
    -------
    >>> model = FirthLogisticRegression()
    >>> model.fit(X, y, term_names=["intercept", "AUC10", "Smoking", "Age60"])
    >>> print(model.summary())
    >>> probs = model.predict_proba(X_new)
    """

    def __init__(self, max_iter=100, tol=1e-8, ridge=1e-10):
        self.max_iter = max_iter
        self.tol = tol
        self.ridge = ridge
        self.beta_ = None
        self.cov_ = None
        self.info_ = None
        self.term_names_ = None

    def fit(self, X, y, term_names=None):
        """Fit Firth logistic regression."""
        self.beta_, self.cov_, self.info_ = firth_logistic(
            X, y, self.max_iter, self.tol, self.ridge
        )
        self.term_names_ = term_names or [f"X{i}" for i in range(len(self.beta_))]
        return self

    def predict_proba(self, X):
        """Predicted probabilities."""
        X = np.asarray(X, dtype=np.float64)
        return _sigmoid(X @ self.beta_)

    def predict(self, X, threshold=0.5):
        """Binary predictions."""
        return (self.predict_proba(X) >= threshold).astype(int)

    def summary(self, alpha=0.05):
        """Print formatted summary table."""
        res = odds_ratios_ci(self.beta_, self.cov_, alpha, self.term_names_)

        lines = []
        lines.append(f"\nFirth Logistic Regression Results")
        lines.append(f"{'=' * 75}")
        lines.append(f"Converged: {self.info_['converged']} "
                      f"(iterations: {self.info_['iterations']})")
        lines.append(f"Penalized log-likelihood: {self.info_['penalized_loglik']:.3f}")
        lines.append(f"{'=' * 75}")
        lines.append(
            f"{'Term':<20} {'Beta':>8} {'SE':>8} {'OR':>8} "
            f"{'CI_low':>8} {'CI_high':>8} {'p':>8}"
        )
        lines.append(f"{'-' * 75}")

        for i in range(len(self.beta_)):
            sig = ""
            if res["p_value"][i] < 0.001:
                sig = "***"
            elif res["p_value"][i] < 0.01:
                sig = "**"
            elif res["p_value"][i] < 0.05:
                sig = "*"

            lines.append(
                f"{res['terms'][i]:<20} "
                f"{res['beta'][i]:>8.4f} "
                f"{res['SE'][i]:>8.4f} "
                f"{res['OR'][i]:>8.3f} "
                f"{res['CI_low'][i]:>8.3f} "
                f"{res['CI_high'][i]:>8.3f} "
                f"{res['p_value'][i]:>7.4f} {sig}"
            )

        lines.append(f"{'=' * 75}")
        lines.append("Signif: *** p<0.001, ** p<0.01, * p<0.05")

        output = "\n".join(lines)
        print(output)
        return res

    def get_results(self, alpha=0.05):
        """Return results as dictionary (for programmatic use)."""
        return odds_ratios_ci(self.beta_, self.cov_, alpha, self.term_names_)


# ─── Demo ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=== Firth Logistic Regression — Demo with Separation ===\n")

    np.random.seed(42)
    n = 200

    x1 = np.random.randn(n)
    x2 = (np.random.rand(n) < 0.1).astype(float)  # rare binary → separation risk

    linpred = -0.5 + 1.2 * x1 + 4.0 * x2
    p = _sigmoid(linpred)
    y = (np.random.rand(n) < p).astype(float)
    y[x2 == 1] = 1  # force quasi-complete separation

    # Design matrix with intercept
    X = np.column_stack([np.ones(n), x1, x2])
    term_names = ["Intercept", "x1 (continuous)", "x2 (rare binary)"]

    # Fit
    model = FirthLogisticRegression()
    model.fit(X, y, term_names=term_names)
    model.summary()

    # Predictions
    probs = model.predict_proba(X)
    print(f"\nMean predicted probability: {probs.mean():.3f}")
    print(f"Classification accuracy: {(model.predict(X) == y).mean() * 100:.1f}%")
