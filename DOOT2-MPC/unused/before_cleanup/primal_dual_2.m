function [alpha, lambda, info] = primal_dual_2(kernels, integral, w_cont, c_sq, params, alpha0, lambda0)
%primal_dual_2  Projected primal-dual iterations for the discretized dual.
%
%   Solves
%       max_alpha  alpha'*w - R(alpha) - (delta/2)*||alpha||^2
%       s.t.       |grad phi(x_i)|^2 <= 1   at every agent (and antenna)
%   with phi = sum_j alpha_j K_field(x - x_j) (+ antenna kernels),
%   w_j = int K_j (rho_N - rho*) dx and the OPTIONAL ridge delta = params.delta_ridge.
%
%   TERMINAL COST. F enters only through R(alpha) = int h(phi(x); rho*(x)) dx,
%   with h'(phi) = rho_T(phi) - rho*, where rho_T(phi) is the terminal state
%   the price phi asks for. params.terminal selects F:
%
%     'l2_relaxed'  F = ||rho - rho*||^2/(2 eps) over SIGNED rho (original)
%                   rho_T = rho* + eps*phi  (can be negative)
%                   R = (eps/2) alpha'*K_reg*alpha
%     'l2'          same F over rho >= 0 (clipped)
%                   rho_T = max(rho* + eps*phi, 0)
%     'kl'          F = kappa * KL(rho | rho_f)  (unnormalised KL)
%                   rho_f = max(rho*, kl_floor*max(rho*)) > 0
%                   rho_T = rho_f .* exp(phi/kappa)  (> 0 automatically)
%                   h(phi; r) = kappa*r_f*(exp(phi/kappa) - 1) - r*phi
%       In 'kl' the rho* terms of w and of h' cancel in the matching
%       agent-kernel projection, so the solver uses the corresponding
%       agent-KDE approximation to the dual of F = kappa*KL(.|rho_f).
%       Mass balance comes from the (approximate) constant mode of phi.
%       This cancellation needs params.w_gain = 1.
%   If params.terminal is absent, params.clip_positive picks 'l2' (true) or
%   'l2_relaxed' (false), exactly as in the previous version.
%
%   WELL-POSEDNESS ('kl'). h is convex and C^1 (exp is continued linearly
%   above phi/kappa = 30, so it cannot overflow), hence the objective is
%   concave, and strongly concave for delta_ridge > 0: unique maximiser.
%   The floor keeps rho_f > 0, so phi stays bounded even where the Lipschitz
%   constraint is inactive: phi/kappa <~ log(max rho_N / min rho_f) ~ 15.
%   The curvature of R is ~ (rho_T/kappa) * Gram, so etap = 0.05 is safe for
%   kappa >~ 0.1; reduce etap for smaller kappa. Nonlinear terminal modes
%   use the agent-kernel integrals from computeTargetint_w.m, rather than a
%   separate grid-by-basis matrix.
%
%   Why the ridge. K is the Gram matrix of heavily overlapping kernels
%   (agent spacing << h), so it is badly conditioned. With delta > 0 every
%   eigenvalue of the curvature is >= delta, the problem is well conditioned,
%   and the output does not depend on pd_iters once converged.
%   The name delta_ridge avoids a clash with params.delta (obstacle width).
%
%   [alpha, lambda] = primal_dual(kernels, integral, w_cont, c_sq, params)
%   [...] = primal_dual(..., alpha0)           warm-start alpha from alpha0
%   [...] = primal_dual(..., alpha0, lambda0)  ...and lambda from lambda0
%   [alpha, lambda, info] = primal_dual(...)   also returns, at the final alpha:
%       info.mass_T     agent-quadrature estimate of int rho_T dx
%       info.neg_mass   int max(-(rho* + eps*phi), 0)   (L2 modes; 0 for 'kl')
%       info.max_grad   max |grad phi| at the agents     (constraint: <= 1)
%       info.step_last  etap*||grad|| at the last iteration (small = converged)
%       info.terminal   terminal mode used
%   Needs computeTargetint_w.m providing agents_ag, target_ag and K_reg.
%   c_sq is accepted for interface compatibility and not used.

pd_iters    = params.pd_iters;
etap        = params.etap;
etad        = params.etad;
epsilon_reg = params.epsilon_reg;
if isfield(params, 'delta_ridge') && ~isempty(params.delta_ridge)
    delta_ridge = params.delta_ridge;
else
    delta_ridge = 0;
end

% ---- terminal cost ------------------------------------------------------
if isfield(params, 'terminal') && ~isempty(params.terminal)
    terminal = lower(char(params.terminal));
elseif isfield(params, 'clip_positive') && ~isempty(params.clip_positive) ...
        && params.clip_positive
    terminal = 'l2';
else
    terminal = 'l2_relaxed';
end
switch terminal
    case 'l2_relaxed', tmode = 0;
    case 'l2',         tmode = 1;
    case 'kl',         tmode = 2;
    otherwise
        error('primal_dual:terminal', ...
            'params.terminal must be ''l2_relaxed'', ''l2'' or ''kl''.');
end
if tmode == 2
    if ~isfield(params, 'kappa') || ~isscalar(params.kappa) || ~(params.kappa > 0)
        error('primal_dual:kappa', 'terminal = ''kl'' needs a scalar params.kappa > 0.');
    end
    kappa = params.kappa;
    if isfield(params, 'kl_floor') && ~isempty(params.kl_floor)
        kl_floor = params.kl_floor;
    else
        kl_floor = 1e-6;
    end
    exp_cap = 30;
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

% ---- agent-kernel integrals --------------------------------------------
% agents_ag(:,i) = int K_field_j(x) K_agent_i(x) dx.  Represent the
% nonlinear terminal density as the agent KDE with weights q_i/N, where q_i
% is its ratio to the current agent density at agent i.
if tmode > 0
    if ~isfield(integral, 'agents_ag') || ~isfield(integral, 'target_ag')
        error('primal_dual_2:missingAgentIntegrals', ...
            'The nonlinear terminal modes need agents_ag and target_ag from computeTargetint_w.');
    end
    if ~isfield(integral, 'rho_target_sm')
        error('primal_dual_2:missingTarget', ...
            'The nonlinear terminal modes need integral.rho_target_sm.');
    end
    field_at_agents = kernels.K_reg(1:N_tot, :);
    density_at_agents = sum(kernels.K_aa(1:N_tot, 1:params.N), 2) / params.N;
    density_floor = max(eps, 1e-12 * max(density_at_agents));
    if isfield(integral, 'rho_target_at_agents')
        rho_s_agents = integral.rho_target_at_agents(:);
    else
        % Backward-compatible fallback for integral structs created before
        % computeTargetint_w returned rho_target_at_agents.
        gram_scale = max(trace(integral.K_reg) / size(integral.K_reg, 1), eps);
        target_coeff = (integral.K_reg + 1e-8 * gram_scale * ...
            eye(size(integral.K_reg))) \ integral.target_ag;
        rho_s_agents = field_at_agents * target_coeff;
    end
    if tmode == 2
        rho_f = max(rho_s_agents, kl_floor * max(max(rho_s_agents), eps));
    end
end

grad_alpha = zeros(size(alpha));
for iter = 1:pd_iters
    v_x = kernels.D_x * alpha;
    v_y = kernels.D_y * alpha;
    violation = v_x.^2 + v_y.^2 - 1;
    lambda = max(0, lambda + etad * violation);

    switch tmode
        case 0   % relaxed L2: eps*K_reg*alpha
            reg_grad = epsilon_reg * (integral.K_reg * alpha);
        case 1   % clipped L2, projected through the agent KDE integrals
            phi_agents = field_at_agents * alpha;
            rho_T_agents = max(rho_s_agents + epsilon_reg * phi_agents, 0);
            q_agents = rho_T_agents ./ max(density_at_agents, density_floor);
            reg_grad = integral.agents_ag(:,1:params.N) * (q_agents(1:params.N) / params.N) ...
                - integral.target_ag;
        case 2   % KL, projected through the agent KDE integrals
            phi_agents = field_at_agents * alpha;
            rho_T_agents = rho_f .* safe_exp(phi_agents / kappa, exp_cap);
            q_agents = rho_T_agents ./ max(density_at_agents, density_floor);
            reg_grad = integral.agents_ag(:,1:params.N) * (q_agents(1:params.N) / params.N) ...
                - integral.target_ag;
    end

    grad_alpha = w_cont - kernels.D_x' * (lambda .* v_x) - kernels.D_y' * (lambda .* v_y) ...
        - reg_grad - delta_ridge * alpha;
    grad_alpha = max(min(grad_alpha, 1), -1);
    alpha = alpha + etap * grad_alpha;
end

if nargout > 2
    switch tmode
        case 2
            phi_agents = field_at_agents * alpha;
            rho_T_agents = rho_f .* safe_exp(phi_agents / kappa, exp_cap);
            info.neg_mass = 0;
        case 1
            phi_agents = field_at_agents * alpha;
            rho_relax_agents = rho_s_agents + epsilon_reg * phi_agents;
            rho_T_agents = max(rho_relax_agents, 0);
            info.neg_mass = mean(max(-rho_relax_agents, 0) ./ ...
                max(density_at_agents, density_floor));
        otherwise
            info.neg_mass = NaN;
    end
    if tmode > 0
        info.mass_T = mean(rho_T_agents(1:params.N) ./ ...
            max(density_at_agents(1:params.N), density_floor));
    else
        info.mass_T = NaN;
    end
    gx = kernels.D_x * alpha;
    gy = kernels.D_y * alpha;
    info.max_grad  = max(hypot(gx(1:N_tot), gy(1:N_tot)));
    info.step_last = etap * norm(grad_alpha);
    info.terminal  = terminal;
end

end


function y = safe_exp(u, cap)
%SAFE_EXP  exp(u) for u <= cap, continued linearly above cap.
%   C^1, convex and non-decreasing, so h stays convex; it never overflows.
y = exp(min(u, cap));
over = u > cap;
y(over) = y(over) .* (1 + u(over) - cap);
end
