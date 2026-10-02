%% main_w.m  --  WENDLAND version of main2.m
% Run parameters_w.m first. Identical flow to main2.m but calls the Wendland
% kernel/target builders. Solver stays primal_dual (the fast projected one).

if ~exist('params', 'var')
    error('Run parameters_w.m before main_w.m.');
end

plot_mode = 'picture';

%% ------------------------ solver / diagnostics flags ---------------------
% warm_start: feed the previous time step's alpha back into primal_dual as its
%   starting point instead of restarting from 0 every step. With a cold start
%   the elementwise clip inside primal_dual caps |alpha_j| at etap*pd_iters, so
%   pd_iters doubles as a step-size gain; warm-starting removes that coupling.
% warm_start_lambda: also carry the dual multiplier over. OFF by default --
%   lambda only ever accumulates (max(0, lambda + etad*violation)), so carrying
%   it across steps can let it grow without bound. Turn on deliberately.
warm_start        = false;
warm_start_lambda = false;
show_diagnostics  = true;

N = params.N; S = params.S; T = params.T; dt_cont = params.dt_cont;
h_agent = params.h_agent; h_ant = params.h_ant;
R_agent = params.R_agent; R_ant = params.R_ant;
var_agent = params.var_agent; C_agent = params.C_agent;
var_ant = params.var_ant; C_ant = params.C_ant;
pd_iters = params.pd_iters; etap = params.etap; etad = params.etad;
epsilon_reg_initial = params.epsilon_reg;
k = 5 / T;
epsilon_reg_schedule = epsilon_reg_initial * exp(-k * (0:T-1));
antw = params.antw; ant_cords = params.ant_cords; M = params.M;
ghost_w = params.ghost_w; N_tot = params.N_tot;

grid_pts = params.grid_pts; dx = params.dx; dy = params.dy;
X = params.X; Y = params.Y; grid_res = params.grid_res;
rho_target = params.rho_target;

if ~exist('pos_cont', 'var')
    error('Run parameters_w.m first so pos_cont exists.');
end

targetParams = params;
targetParams.grid_pts = grid_pts;
targetParams.rho_target = rho_target;
targetParams.dx = dx; targetParams.dy = dy;

%% ------------------------------ moving target ---------------------------
% Rigid rotation of the star, set up in parameters_w.m. rho_target is rebuilt
% at the top of every step, so everything downstream that reads it -- the
% target integrals, the per-step error diagnostic, and panel 1 of the final
% plot -- automatically tracks the moving target.
if isfield(params, 'target_rot') && params.target_rot
    target_rot   = true;
    target_omega = params.target_omega;
else
    target_rot   = false;
    target_omega = 0;
end
if target_rot
    % the precomputed smoothed target is only valid at phase 0; blank it so
    % computeTargetint_w rebuilds it from the current rho_target (that branch
    % is only reached when no_smoothing is false)
    targetParams.rho_target_sm = [];
end

pos_hist = zeros(N_tot, 2, T);
alpha_hist = zeros(N_tot + M, T);

%% ------------------------- diagnostic time series ------------------------
% (named diagn, not diag -- diag is a MATLAB builtin)
diagn = struct();
diagn.d_alpha      = zeros(T, 1);   % ||alpha_t - alpha_{t-1}||  (outer convergence)
diagn.rel_d_alpha  = zeros(T, 1);   % same, relative to ||alpha_{t-1}||
diagn.alpha_inf    = zeros(T, 1);   % max|alpha_j|; compare to the cold-start cap
diagn.max_gradphi  = zeros(T, 1);   % max|grad phi| -- Kantorovich needs <= 1
diagn.frac_violate = zeros(T, 1);   % fraction of agents with |grad phi| > 1
diagn.mean_speed   = zeros(T, 1);
diagn.max_speed    = zeros(T, 1);
diagn.w_norm       = zeros(T, 1);   % ||w|| mismatch residual (drives the flow)
diagn.err_L1       = zeros(T, 1);   % ||rho_N - mu*||_1 at every step
diagn.err_L2       = zeros(T, 1);
diagn.warm_start        = warm_start;
diagn.warm_start_lambda = warm_start_lambda;
diagn.pd_iters          = pd_iters;
diagn.alpha_cap_cold    = etap * pd_iters;   % |alpha_j| ceiling if cold-started

fprintf('Starting Evolution (Wendland)...\n');
if warm_start && warm_start_lambda
    fprintf('  warm start: alpha AND lambda carried over between steps\n');
elseif warm_start
    fprintf('  warm start: alpha carried over, lambda reset each step\n');
else
    fprintf('  cold start: alpha reset to 0 each step (|alpha_j| <= %.3f)\n', ...
            etap * pd_iters);
end
if target_rot
    fprintf('  moving target: %.5f rad/step (%.2f revolutions over %d steps)\n', ...
            target_omega, target_omega*T/(2*pi), T);
else
    fprintf('  static target\n');
end
alpha_prev = zeros(N_tot + M, 1);
lambda_prev = zeros(N_tot + M, 1);
for t = 1:T
    params.epsilon_reg = epsilon_reg_schedule(1);
    targetParams.epsilon_reg = epsilon_reg_schedule(1);

    % ---- moving target: rebuild the star at this step's rotation phase ----
    % phase 0 at t = 1, so the run starts from exactly the static target.
    if target_rot
        rho_target = star_target_w(params, target_omega * (t - 1));
        targetParams.rho_target = rho_target;
    end

    kernels  = computeKernels_w(pos_cont, ant_cords, params);
    integral = computeTargetint_w(pos_cont, ant_cords, targetParams);

    int_agents_ag  = sum(integral.agents_ag(1:N, 1:N), 2);
    if ~isempty(integral.agents_ant)
        int_agents_ant = sum(integral.agents_ant, 1)'/N;
    else
        int_agents_ant = [];
    end
    int_target_ag  = integral.target_ag(1:N);
    int_target_ant = integral.target_ant;

    w_agent = zeros(N_tot, 1);
    w_agent(1:N) = (int_agents_ag / N) - int_target_ag;
    w_agent(N+1:N_tot) = sum(integral.agents_ag(N+1:N_tot, N+1:N_tot), 2) + ghost_w;
    if ~isempty(int_agents_ant)
        w_ant = antw * (int_agents_ant/N - int_target_ant);
    else
        w_ant = [];
    end
    w_cont = [w_agent; w_ant];
    w_cont = params.w_gain * w_cont;
    w_cont(N+1:N_tot) = w_cont(N+1:N_tot) + ghost_w;

    c_sq = conformal(pos_cont, ant_cords, params);
    if warm_start
        if warm_start_lambda
            [alpha, lambda] = primal_dual(kernels, integral, w_cont, c_sq, params, alpha_prev, lambda_prev);
        else
            [alpha, lambda] = primal_dual(kernels, integral, w_cont, c_sq, params, alpha_prev);
        end
    else
        [alpha, lambda] = primal_dual(kernels, integral, w_cont, c_sq, params);
    end

    % ---- outer convergence metrics (before alpha_prev is overwritten) ----
    d_alpha = alpha - alpha_prev;
    diagn.d_alpha(t)     = norm(d_alpha);
    diagn.rel_d_alpha(t) = norm(d_alpha) / max(norm(alpha_prev), eps);
    diagn.alpha_inf(t)   = max(abs(alpha));

    alpha_prev = alpha;
    lambda_prev = lambda;
    alpha_hist(:, t) = alpha;

    v_agents_x = -(kernels.D_x(1:N_tot, :) * alpha);
    v_agents_y = -(kernels.D_y(1:N_tot, :) * alpha);
    v_agents = [v_agents_x, v_agents_y];
    v_move = v_agents;

    % Congestion term: repulsion from density gradient (prevents collapse/clustering)
    grad_rhoN_x = sum(kernels.D_x_agent_aa(1:N, 1:N), 2) / N;
    grad_rhoN_y = sum(kernels.D_y_agent_aa(1:N, 1:N), 2) / N;
    grad_rhoN = [grad_rhoN_x, grad_rhoN_y];


    % Ghosts and antennas don't move
    v_move(N+1:N_tot, :) = 0;
    pos_cont = pos_cont + dt_cont * v_move;
    pos_hist(:, :, t) = pos_cont;

    % ---- constraint / residual / error metrics ---------------------------
    % |grad phi| for the real (moving) agents. The Kantorovich constraint the
    % dual step enforces is |grad phi| <= 1, so max_gradphi >> 1 means the
    % lambda updates are not keeping up.
    gmag = hypot(v_agents_x(1:N), v_agents_y(1:N));
    diagn.max_gradphi(t)  = max(gmag);
    diagn.frac_violate(t) = mean(gmag > 1);
    diagn.mean_speed(t)   = mean(gmag);
    diagn.max_speed(t)    = max(gmag);
    diagn.w_norm(t)       = norm(w_agent(1:N));

    % rho_N vs target, reusing the density kernel computeTargetint_w already
    % built this step (K_grid_agent is grid_pts x N_tot at h_agent) -- free.
    rho_N = sum(integral.K_grid_agent(:, 1:N), 2) / N;
    mass = sum(rho_N) * dx * dy;
    if mass > 0
        rho_N = rho_N / mass;
    end
    d_rho = rho_N - rho_target;
    diagn.err_L1(t) = sum(abs(d_rho)) * dx * dy;
    diagn.err_L2(t) = sqrt(sum(d_rho.^2) * dx * dy);
end
fprintf('Simulation complete (Wendland).\n');

% Continuous Wasserstein diagnostic for the final agent KDE and smoothed target.
target_final = rho_target;
target_final_sm = smooth_on_grid_wendland(target_final, params, h_agent);
diffs_final = permute(grid_pts, [1 3 2]) - permute(pos_cont(1:N, :), [3 1 2]);
sq_dists_final = sum(diffs_final.^2, 3);
rho_agents_final = sum(wendland_kernel(sq_dists_final, h_agent), 2) / N;
rho_agents_final = rho_agents_final / (sum(rho_agents_final) * dx * dy);
[W_final, W_info] = computeW2(rho_agents_final, target_final_sm, x, y, 'p', 1);
fprintf('Final continuous W1 = %.6f (source integral = %.6f, target integral = %.6f, ...\n', ...
    W_final, sum(rho_agents_final) * dx * dy, sum(target_final_sm) * dx * dy);
fprintf('  Sinkhorn marginal error = %.3e, converged = %d, epsilon = %.3e\n', ...
    W_info.marg_err, W_info.converged, W_info.epsilon);

%% ---------------------------- diagnostics summary ------------------------
tail_idx = max(1, round(0.9*T)):T;          % last 10% of the run
[best_L1, best_t] = min(diagn.err_L1);
fprintf('\n--- Diagnostics (pd_iters = %d, warm_start = %d) ---\n', pd_iters, warm_start);
fprintf('  L1 error   final = %.4f | best = %.4f @ t = %d | mean(last 10%%) = %.4f\n', ...
        diagn.err_L1(T), best_L1, best_t, mean(diagn.err_L1(tail_idx)));
fprintf('  L1 std over last 10%%      = %.4f   (large => orbiting, not settled)\n', ...
        std(diagn.err_L1(tail_idx)));
fprintf('  ||alpha_t - alpha_{t-1}||  final = %.3e | mean(last 10%%) = %.3e\n', ...
        diagn.d_alpha(T), mean(diagn.d_alpha(tail_idx)));
fprintf('  relative d_alpha           mean(last 10%%) = %.3e\n', ...
        mean(diagn.rel_d_alpha(tail_idx)));
fprintf('  max|alpha_j|               final = %.4f  (cold-start cap would be %.4f)\n', ...
        diagn.alpha_inf(T), etap * pd_iters);
fprintf('  max|grad phi|              final = %.4f | mean(last 10%%) = %.4f  (constraint: <= 1)\n', ...
        diagn.max_gradphi(T), mean(diagn.max_gradphi(tail_idx)));
fprintf('  agents violating |grad phi| <= 1, mean(last 10%%) = %.1f%%\n', ...
        100*mean(diagn.frac_violate(tail_idx)));
fprintf('  ||w|| residual             final = %.3e | mean(last 10%%) = %.3e\n', ...
        diagn.w_norm(T), mean(diagn.w_norm(tail_idx)));
fprintf('---------------------------------------------------\n\n');

plotData = struct( ...
    'X', X, 'Y', Y, ...
    'rho_target_mat', reshape(rho_target, grid_res, grid_res), ...
    'pos_hist', pos_hist, 'alpha_hist', alpha_hist, ...
    'grid_pts', grid_pts, 'grid_res', grid_res, 'N', N, 'dt_cont', dt_cont, ...
    'h_agent', h_agent, 'R_agent', R_agent, 'var_agent', var_agent, 'C_agent', C_agent, ...
    'ant_cords', ant_cords, 'R_ant', R_ant, 'var_ant', var_ant, 'C_ant', C_ant, ...
    'h_field', params.h_field, 'h_ant', params.h_ant);

plotresults_w(plot_mode, plotData);

if show_diagnostics
    plot_diagnostics_w(diagn, params, sprintf('main\\_w  ($h_{agent} = %.2f$, $h_{field} = %.2f$)', ...
                                              h_agent, params.h_field));
end
