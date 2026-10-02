function [alpha, lambda] = primal_dual(kernels, integral, w, params, alpha, lambda)
%primal_dual  Projected primal-dual iterations for the DOOT dual problem.
%
%   Solves
%       max_alpha  alpha'*w - R(alpha) - (delta/2)*||alpha||^2
%       s.t.       |grad phi(x_i)|^2 <= 1   at every agent
%   with phi = sum_j alpha_j K_field(x - x_j), w the density mismatch,
%   delta = params.delta_ridge and R the dual of the terminal cost F,
%   chosen by params.terminal:
%
%     'l2'  F = ||rho - rho*||^2 / (2 eps),  eps = params.epsilon_reg
%           R = (eps/2) alpha'*K_reg*alpha,  K_reg = Gram matrix of the field basis
%     'kl'  F = kappa * KL(rho | rho_f),  rho_f = max(rho*, kl_floor*max(rho*))
%           grad R = int K_field (rho_T - rho*) dx with the terminal density
%           rho_T = rho_f .* exp(phi/kappa), integrated with the agent KDE
%           (agents_ag). exp is continued linearly above phi/kappa = 30.
%
%   alpha, lambda: warm start in, solution out.

N = params.N;

if strcmp(params.terminal, 'kl')
    % KL terms that stay fixed during the iterations
    kappa = params.kappa;
    exp_cap = 30;
    field_at_agents = kernels.K_field;                  % phi(x_i) = field_at_agents(i,:) * alpha
    density_at_agents = sum(kernels.K_agent, 2) / N;    % rho_N(x_i)
    density_at_agents = max(density_at_agents, max(eps, 1e-12 * max(density_at_agents)));
    rho_s = integral.rho_target_at_agents;              % rho*(x_i)
    rho_f = max(rho_s, params.kl_floor * max(max(rho_s), eps));
end

for iter = 1:params.pd_iters
    % Dual step: lambda grows where the speed limit |grad phi| <= 1 is violated
    v_x = kernels.D_x * alpha;
    v_y = kernels.D_y * alpha;
    lambda = max(0, lambda + params.etad * (v_x.^2 + v_y.^2 - 1));

    % Gradient of the terminal cost
    switch params.terminal
        case 'l2'
            reg_grad = params.epsilon_reg * (integral.K_reg * alpha);
        case 'kl'
            rho_T = rho_f .* safe_exp(field_at_agents * alpha / kappa, exp_cap);
            reg_grad = integral.agents_ag * (rho_T ./ density_at_agents / N) - integral.target_ag;
    end

    % Primal step: clipped gradient ascent
    grad_alpha = w - kernels.D_x' * (lambda .* v_x) - kernels.D_y' * (lambda .* v_y) ...
        - reg_grad - params.delta_ridge * alpha;
    grad_alpha = max(min(grad_alpha, 1), -1);
    alpha = alpha + params.etap * grad_alpha;
end

end


function y = safe_exp(u, cap)
%SAFE_EXP  exp(u) for u <= cap, continued linearly above cap (never overflows).
y = exp(min(u, cap));
over = u > cap;
y(over) = y(over) .* (1 + u(over) - cap);
end
