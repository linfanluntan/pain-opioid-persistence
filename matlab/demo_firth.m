%% demo_firth.m
%  Demonstration of Firth penalized logistic regression with separation.
%  ASTRO 2026 | Mahin, Nkuku, He, Fuller, Moreno, Javed
%
%  This script:
%    1. Generates simulated data with quasi-complete separation
%    2. Shows that standard MLE fails (huge coefficients)
%    3. Shows that Firth produces finite, stable estimates
%    4. Demonstrates the FU2-like parsimonious model structure

clear; clc;
fprintf('=== Firth Logistic Regression — Demo ===\n\n');

%% 1) Simulate data with separation
rng(42);
n = 200;

x1 = randn(n, 1);                            % continuous predictor
x2 = double(rand(n, 1) < 0.1);              % rare binary (10% prevalence)

linpred = -0.5 + 1.2*x1 + 4.0*x2;
p = 1 ./ (1 + exp(-linpred));
y = double(rand(n, 1) < p);

% Force separation: all x2=1 cases are positive
y(x2 == 1) = 1;

fprintf('Simulated data: N = %d\n', n);
fprintf('Events = %d, Non-events = %d\n', sum(y), sum(1-y));
fprintf('x2 prevalence: %.1f%% (n=%d)\n', 100*mean(x2), sum(x2));
fprintf('All x2=1 patients have y=1 → quasi-complete separation\n\n');

%% 2) Standard MLE (for comparison — expect huge coefficients)
fprintf('--- Standard MLE (glmfit) ---\n');
X_pred = [x1, x2];

try
    [b_mle, ~, stats_mle] = glmfit(X_pred, y, 'binomial', 'logit');
    fprintf('Intercept: %.3f\n', b_mle(1));
    fprintf('x1:        %.3f (SE = %.3f)\n', b_mle(2), stats_mle.se(2));
    fprintf('x2:        %.3f (SE = %.3f)\n', b_mle(3), stats_mle.se(3));
    fprintf('x2 OR:     %.1f\n', exp(b_mle(3)));

    if stats_mle.se(3) > 10
        fprintf('\n⚠ SEPARATION DETECTED: x2 SE = %.1f\n', stats_mle.se(3));
        fprintf('  MLE estimates are UNRELIABLE.\n');
    end
catch ME
    fprintf('MLE failed: %s\n', ME.message);
end

%% 3) Firth penalized logistic regression
fprintf('\n--- Firth Penalized LR ---\n');
X = [ones(n, 1), x1, x2];   % include intercept
term_names = {'Intercept', 'x1 (continuous)', 'x2 (rare binary)'};

[beta, cov_mat, info] = firth_logistic(X, y);

fprintf('Converged: %d (iterations: %d)\n', info.converged, info.iterations);
fprintf('Penalized log-likelihood: %.3f\n\n', info.pen_loglik);

% OR and CI
se = sqrt(diag(cov_mat));
z = 1.96;
OR = exp(beta);
CI_lo = exp(beta - z*se);
CI_hi = exp(beta + z*se);
p_val = 2 * (1 - normcdf(abs(beta ./ se)));

fprintf('%-25s %8s %8s %10s %8s\n', 'Term', 'OR', 'CI_low', 'CI_high', 'p');
fprintf('%s\n', repmat('-', 1, 65));
for i = 1:length(beta)
    fprintf('%-25s %8.3f %8.3f %10.3f %8.4f\n', ...
        term_names{i}, OR(i), CI_lo(i), CI_hi(i), p_val(i));
end

%% 4) Predictions
pred_prob = 1 ./ (1 + exp(-X * beta));
pred_class = double(pred_prob >= 0.5);
accuracy = mean(pred_class == y) * 100;
fprintf('\nClassification accuracy: %.1f%%\n', accuracy);
fprintf('Mean predicted probability: %.3f\n', mean(pred_prob));

%% 5) FU2-like parsimonious model demo
fprintf('\n\n--- FU2-like Parsimonious Model Demo ---\n');

n2 = 106;
rng(7);

AUC10 = max(0, min(7, randn(n2, 1)*1.5 + 3));   % AUC scaled per 10%
Smk = double(rand(n2, 1) < 0.4);                  % smoking ever/never
Age60 = double(rand(n2, 1) < 0.5);                % age >=60
Tx_PORT = double(rand(n2, 1) < 0.3);              % PORT dummy
Tx_CRT = double(rand(n2, 1) < 0.5);               % Def_ChemoRT dummy

logit2 = -1.5 + 0.4*AUC10 + 0.6*Smk - 0.4*Age60 + 0.2*Tx_CRT;
p2 = 1 ./ (1 + exp(-logit2));
y2 = double(rand(n2, 1) < p2);

X2 = [ones(n2, 1), AUC10, Smk, Age60, Tx_PORT, Tx_CRT];
names2 = {'Intercept', 'AUC10 (per 10%)', 'Smk_Ever', 'Age60plus', ...
           'TxIntensity_PORT', 'TxIntensity_DefChemoRT'};

fprintf('N = %d, Events = %d, Non-events = %d\n', n2, sum(y2), sum(1-y2));
fprintf('Parameters: %d (EPV = %.1f)\n\n', size(X2,2)-1, sum(y2)/(size(X2,2)-1));

[beta2, cov2, info2] = firth_logistic(X2, y2);

se2 = sqrt(diag(cov2));
OR2 = exp(beta2);
CI2_lo = exp(beta2 - z*se2);
CI2_hi = exp(beta2 + z*se2);
p2_val = 2 * (1 - normcdf(abs(beta2 ./ se2)));

fprintf('%-30s %8s %8s %10s %8s\n', 'Term', 'OR', 'CI_low', 'CI_high', 'p');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:length(beta2)
    sig = '';
    if p2_val(i) < 0.001, sig = '***';
    elseif p2_val(i) < 0.01, sig = '**';
    elseif p2_val(i) < 0.05, sig = '*';
    end
    fprintf('%-30s %8.3f %8.3f %10.3f %7.4f %s\n', ...
        names2{i}, OR2(i), CI2_lo(i), CI2_hi(i), p2_val(i), sig);
end

fprintf('\nConverged: %d | Penalized LL: %.3f\n', info2.converged, info2.pen_loglik);
fprintf('\n✓ Demo complete.\n');
