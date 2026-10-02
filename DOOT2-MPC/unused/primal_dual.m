function [alpha, lambda] = primal_dual(kernels, integral, w_cont, c_sq, params, alpha0, lambda0)
%primal_dual  Projected primal-dual iterations for the discretized dual (12).
%
%   Solves
%       max_alpha  alpha'*w - (eps/2)*alpha'*K*alpha - (delta/2)*||alpha||^2
%       s.t.       |grad phi(x_i)|^2 <= 1   at every agent (and antenna)
%   with eps = params.epsilon_reg, K = integral.K_reg and the OPTIONAL ridge
%   delta = params.delta_ridge. If delta_ridge is absent or 0, the iteration
%   is exactly the original one.
%
%   Why the ridge. K is the Gram matrix of heavily overlapping kernels
%   (agent spacing << h), so its condition number is ~1e12. Fixed-step
%   gradient iterations converge along an eigenvector of eps*K at rate
%   ~etap*eps*lambda_k per iteration, i.e. never along the small-eigenvalue
%   (wiggly) directions, and the output then depends on pd_iters. With
%   delta > 0 every eigenvalue of eps*K + delta*I is >= delta, so the problem
%   is well conditioned and pd_iters iterations converge. Early stopping
%   after n iterations acts like a ridge with delta ~ 1/(etap*n); delta = 4e-3
%   (etap = 0.05, 3000 iterations) therefore reproduces the current fields,
%   but as the converged solution of a well-defined problem.
%
%   The name delta_ridge avoids a clash with params.delta (obstacle width).
%
%   [alpha, lambda] = primal_dual(kernels, integral, w_cont, c_sq, params)
%   [...] = primal_dual(..., alpha0)           warm-start alpha from alpha0
%   [...] = primal_dual(..., alpha0, lambda0)  ...and lambda from lambda0
%   c_sq is accepted for interface compatibility and is not used.

pd_iters    = params.pd_iters;
etap        = params.etap;
etad        = params.etad;
epsilon_reg = params.epsilon_reg;
if isfield(params, 'delta_ridge') && ~isempty(params.delta_ridge)
    delta_ridge = params.delta_ridge;
else
    delta_ridge = 0;
end

N_tot = size(kernels.D_x, 1) - params.M;
M = params.M;

if nargin < 6 || isempty(alpha0)
    alpha = zeros(N_tot + M, 1);
else
    alpha = alpha0(:);
end
if nargin < 7 || isempty(lambda0)
    lambda = zeros(N_tot + M, 1);
else
    lambda = lambda0(:);
end

for iter = 1:pd_iters
    v_x = kernels.D_x * alpha;
    v_y = kernels.D_y * alpha;
    violation = v_x.^2 + v_y.^2 - 1;
    lambda = max(0, lambda + etad * violation);
    grad_alpha = w_cont - kernels.D_x' * (lambda .* v_x) - kernels.D_y' * (lambda .* v_y) ...
        - epsilon_reg * (integral.K_reg * alpha) - delta_ridge * alpha;
    grad_alpha = max(min(grad_alpha, 1), -1);
    alpha = alpha + etap * grad_alpha;
end

end
