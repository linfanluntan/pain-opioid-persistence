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
%   small-sample bias — ideal for sparse clinical datasets.
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
%       info    - struct with fields:
%                   .iterations  - number of iterations used
%                   .converged   - logical, whether convergence was achieved
%                   .loglik      - standard log-likelihood at solution
%                   .pen_loglik  - penalized log-likelihood at solution
%
%   Reference:
%       Firth D. Bias reduction of maximum likelihood estimates.
%       Biometrika. 1993;80(1):27-38.
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
