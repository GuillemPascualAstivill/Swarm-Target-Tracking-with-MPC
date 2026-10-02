function [v, beta, lam, info] = plan_fixed(pos, beta, lam, nb, rho_s, params)
%PLAN_FIXED  One DOOT planning step with the potential phi on a FIXED node basis.
%
%   phi(x) = sum_a beta_a K_node(x - A_a), with the nodes A_a from build_nodes.m.
%   Projected primal-dual iterations solve the dual of
%       min_rho  W1(rho_N, rho) + F(rho),     F = terminal cost (params.terminal)
%   i.e.
%       max_beta  int phi rho_N dx - F*(phi) - (delta/2)*||beta||^2
%       s.t.      |grad phi|^2 <= 1   at every agent and at every node
%   with rho_N the Wendland KDE of the swarm (support h_agent), F* the convex
%   conjugate of F and delta = delta_ridge. F enters only through the terminal
%   density rho_T(phi) = grad F*(phi) that the price phi asks for:
%       'kl'          F = kappa * KL(rho | rho_f)     rho_T = rho_f .* exp(phi/kappa)  (> 0)
%                     rho_f = max(rho*, kl_floor*max(rho*)); exp is continued
%                     linearly above phi/kappa = 30, so it never overflows
%       'l2'          F = ||rho - rho*||^2 / (2 eps)  rho_T = max(rho* + eps*phi - mu, 0)  (>= 0)
%                     restricted to P = {rho >= 0, int rho dx = m_N}, m_N = int rho_N dx:
%                     rho_T is the exact L2 projection of rho* + eps*phi onto P, and the
%                     scalar mu = eps*lambda(phi) fixes the mass (found exactly by sorting)
%       'l2_relaxed'  F = ||rho - rho*||^2 / (2 eps)  rho_T = rho* + eps*phi  (can be < 0)
%   The beta-gradient is int K_node(x - A_a) (rho_N - rho_T(phi)) dx (+ constraint
%   and ridge terms), computed as w - int K_node (rho_T - rho*) dx with
%   w_a = int K_node(x - A_a) (rho_N(x) - rho*(x)) dx; the rho* terms cancel
%   exactly (same grid quadrature). The agents move with v_i = -grad phi(x_i).
%
%   Why a fixed basis: with kernels centred at the agents, phi can only be
%   built from bumps/dips centred at agents. Where the swarm is thinner than
%   the target (phi < 0) those dips make neighbouring agents ATTRACT each
%   other, and an agent with no neighbour within h gets exactly zero velocity.
%   A fixed basis removes both effects and lets phi reach empty regions.
%
%   Inputs : pos (N x 2) agent positions; beta (M x 1) and lam (N+M x 1) warm
%            start; nb from build_nodes; rho_s target density on the grid; params.
%   Outputs: v (N x 2) agent velocities; beta, lam the solution;
%            info.rhoN       rho_N on the grid
%            info.phi_grid   phi on the grid
%            info.mass_T     terminal mass int rho_T dx
%            info.mass_N     swarm mass int rho_N dx (the mass rho_T is projected to in 'l2')
%            info.step_last  etap*||grad|| at the last iteration (small = converged)

dA = params.dx * params.dy;         % grid cell area
etap  = params.etap;
etad  = params.etad;
delta = params.delta_ridge;
rho_s = rho_s(:);

% Swarm density on the grid, its mass and the linear term w
rhoN = kde_grid(params.grid_pts, pos, params.h_agent);
mN   = dA * sum(rhoN);              % = 1 up to KDE truncation at the boundary
w    = nb.KgnT * ((rhoN - rho_s) * dA);

% Terminal density rho_T(phi) for the chosen terminal cost F (see above)
switch params.terminal
    case 'kl'
        kappa = params.kappa;
        rho_f = max(rho_s, params.kl_floor * max(rho_s));    % target floored to stay > 0
        rho_T = @(phi) rho_f .* safe_exp(phi / kappa, 30);
    case 'l2'
        eps_reg = params.epsilon_reg;
        rho_T = @(phi) proj_mass(rho_s + eps_reg * phi, dA, mN);   % max(. - mu, 0), mass mN
    case 'l2_relaxed'
        eps_reg = params.epsilon_reg;
        rho_T = @(phi) rho_s + eps_reg * phi;
    otherwise
        error('params.terminal must be ''kl'', ''l2'' or ''l2_relaxed''.');
end

% Node-kernel gradients at the agents: grad phi(x_i) = [Dax(i,:); Day(i,:)] * beta.
% Stacked on the node rows, they give the speed limit at the agents and the nodes.
[~, Dfac] = wendland_kernel(pdist2(pos, nb.A).^2, nb.h);
Dax = sparse(Dfac .* (pos(:,1) - nb.A(:,1).'));
Day = sparse(Dfac .* (pos(:,2) - nb.A(:,2).'));
Dx  = [Dax; nb.Dnx];   DxT = Dx.';
Dy  = [Day; nb.Dny];   DyT = Dy.';

grad = zeros(size(beta));
for iter = 1:params.pd_iters
    % Dual step: lam grows where the speed limit |grad phi| <= 1 is violated
    vx  = Dx * beta;
    vy  = Dy * beta;
    lam = max(0, lam + etad * (vx.^2 + vy.^2 - 1));

    % Primal step: clipped gradient ascent on beta
    phi_g = nb.Kgn * beta;                                  % phi on the grid
    hp    = rho_T(phi_g) - rho_s;                           % rho_T - rho* on the grid
    grad  = w - DxT * (lam .* vx) - DyT * (lam .* vy) - dA * (nb.KgnT * hp) - delta * beta;
    grad  = max(min(grad, 1), -1);
    beta  = beta + etap * grad;
end

v = -[Dax * beta, Day * beta];      % v_i = -grad phi(x_i)

phi_g = nb.Kgn * beta;
info.rhoN      = rhoN;
info.phi_grid  = phi_g;
info.mass_T    = dA * sum(rho_T(phi_g));
info.mass_N    = mN;
info.step_last = etap * norm(grad);
end


function rho = kde_grid(G, P, h)
%KDE_GRID  Wendland KDE of the points P evaluated on the grid G (chunked to limit memory).
nG = size(G, 1);  rho = zeros(nG, 1);  chunk = 2000;
for s = 1:chunk:nG
    idx = s:min(s+chunk-1, nG);
    rho(idx) = sum(wendland_kernel(pdist2(G(idx,:), P).^2, h), 2);
end
rho = rho / size(P, 1);
end


function y = safe_exp(u, cap)
%SAFE_EXP  exp(u) for u <= cap, continued linearly above cap (C^1, convex).
y = exp(min(u, cap));
over = u > cap;
y(over) = y(over) .* (1 + u(over) - cap);
end


function s = proj_mass(u, dA, m)
%PROJ_MASS  L2 projection of u onto {s >= 0, dA*sum(s) = m}:  s = max(u - mu, 0),
%   with mu the unique scalar giving mass m. Exact, O(n log n): while the active
%   set is the k largest entries, the mass is linear in mu, giving mu_k below; the
%   true active set is the largest k with us(k) > mu_k.
us = sort(u, 'descend');
mu = (cumsum(us) - m / dA) ./ (1:numel(us)).';
k  = find(us > mu, 1, 'last');
s  = max(u - mu(k), 0);
end
