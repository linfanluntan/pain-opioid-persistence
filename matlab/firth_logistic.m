function [beta, cov_mat, info] = firth_logistic(X, y, maxIter, tol)
%FIRTH_LOGISTIC Firth (Jeffreys-prior penalized) logistic regression.
%
%   [beta, cov_mat, info] = firth_logistic(X, y)
%   [beta, cov_mat, info] = firth_logistic(X, y, maxIter, tol)
%
%   Fits logistic regression with Firth's bias-reducing penalty,
%   equivalent to adding a Jeffreys invariant prior:
%
%       l*(beta) = l(beta) + 0.5 * log|I(beta)|
%
%   This prevents coefficient divergence under separation and corrects
%   small-sample bias -- ideal for sparse clinical datasets.
%
%   WHEN TO USE:
%   In our FU2 cohort (N=106, 47 events, 23+ parameters), standard MLE
%   produces quasi-complete separation with ORs of 47,924 and 7,532,188.
%   Firth regression guarantees finite estimates without a tuning
%   hyperparameter. The primary Pain AUC predictor is unaffected by
%   separation (continuous, broadly distributed, outcome overlap).
%
%   HIERARCHY (docs/statistical_notes.md Section 4.5):
%     1. Firth penalized LR     -- primary inference (this function)
%     2. Parsimonious model     -- collapsed categories
%     3. L2 regularized LR      -- sensitivity only
%     4. Standard MLE           -- not reportable under separation
%
%   Inputs:
%       X       - (n x p) design matrix. Include intercept column if desired.
%       y       - (n x 1) binary outcomes (0 or 1).
%       maxIter - (optional) maximum Newton iterations, default 100.
%       tol     - (optional) convergence tolerance, default 1e-8.
%
%   Outputs:
%       beta    - (p x 1) Firth-penalized coefficient estimates.
%       cov_mat - (p x p) approximate covariance (inverse Fisher info).
%       info    - struct: .iterations, .converged, .loglik, .pen_loglik
%
%   References:
%       Firth D. Biometrika. 1993;80(1):27-38.
%       Albert A, Anderson JA. Biometrika. 1984;71(1):1-10.
%       Heinze G, Schemper M. Stat Med. 2002;21(16):2409-2419.
%       King G, Zeng L. Polit Anal. 2001;9(2):137-163.
%       See docs/references.md for complete bibliography.
%
%   R equivalent: logistf::logistf() (Heinze et al., 2024)
%   Python equivalent: python/firth_logistic.py
%
%   ASTRO 2026 | Mahin, Nkuku, He, Fuller, Moreno, Javed

    if nargin < 3 || isempty(maxIter), maxIter = 100; end
    if nargin < 4 || isempty(tol), tol = 1e-8; end

    X = double(X);
    y = double(y(:));
    [n, p] = size(X);

    beta = zeros(p, 1);
    ridge = 1e-10;
    converged = false;

    for it = 1:maxIter
        eta = X * beta;
        mu = 1 ./ (1 + exp(-eta));

        % Weights W = mu .* (1 - mu)
        w = mu .* (1 - mu);
        w = max(w, 1e-12);  % prevent zero weights

        % Fisher information: I = X' * diag(w) * X
        % Efficient: X' * (X .* w) avoids forming full diagonal matrix
        I_mat = X' * (X .* w) + ridge * eye(p);

        % Inverse information
        I_inv = inv(I_mat);

        % Hat matrix diagonal: h = diag(W^{1/2} X I^{-1} X' W^{1/2})
        sqrtw = sqrt(w);
        A = X .* sqrtw;             % (n x p)
        AIinv = A * I_inv;          % (n x p)
        h = sum(AIinv .* A, 2);     % (n x 1)

        % Firth adjustment: a = 0.5 * (1 - 2*mu) .* h
        a = 0.5 * (1 - 2*mu) .* h;

        % Modified score: U* = X' * (y - mu + a)
        Ustar = X' * (y - mu + a);

        % Newton step: beta_new = beta + I^{-1} * U*
        step = I_inv * Ustar;
        beta_new = beta + step;

        if max(abs(beta_new - beta)) < tol
            beta = beta_new;
            converged = true;
            break;
        end

        beta = beta_new;
    end

    % Final covariance at solution
    eta = X * beta;
    mu = 1 ./ (1 + exp(-eta));
    w = mu .* (1 - mu);
    w = max(w, 1e-12);
    I_mat = X' * (X .* w) + ridge * eye(p);
    cov_mat = inv(I_mat);

    % Log-likelihoods
    loglik = sum(y .* log(mu + 1e-15) + (1 - y) .* log(1 - mu + 1e-15));
    pen_loglik = loglik + 0.5 * log(det(I_mat));

    % Output info
    info.iterations = it;
    info.converged = converged;
    info.loglik = loglik;
    info.pen_loglik = pen_loglik;
end
