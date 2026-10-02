clc; close all;
%% main_w_final.m  --  DOOT swarm control toward a rotating star target,
% optionally tracked by the low-level MPC. Run parameters_w.m first.

if ~exist('params', 'var') || ~isfield(params, 'pos0')
    error('Run parameters_w.m before main_w_final.m.');
end

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(fileparts(script_dir), 'WINDS'));

%% Settings
use_mpc = false;       % false: agents move with the DOOT velocity
                       % true:  agents track the DOOT plan with the low-level MPC
n_inner = 10;          % MPC steps per DOOT step
Q_x = 1;               % MPC position weight
Q_u = 3.5;             % MPC control weight
Nh = 10;               % MPC horizon
gamma = 1;             % MPC terminal-cost scale
u_max = [];            % MPC acceleration bound |u| <= u_max on each axis
                       % []:    unconstrained, constant gain from idare
                       % value: box-constrained QP over the horizon (quadprog)
wind_model = 'none';   % MPC only: 'none', 'uniform' or 'taylor_green'
wind_U = 0.01;         % wind speed
wind_theta = 0;        % 'uniform' wind direction [rad]
wind_L = 6;            % 'taylor_green' vortex period
drift_mode = 'measured';  % wind the MPC compensates for:
                       % 'measured':  the true local wind (wind sensor)
                       % 'estimated': estimated from the agent's own motion
                       % 'none':      no compensation
drift_alpha = 0.2;     % 'estimated': estimator gain in (0, 1], 1 = last step only

N = params.N;
T = params.T;
dt_outer = params.dt_cont;
dt_inner = dt_outer / n_inner;

if use_mpc
    mpc = low_level_mpc_setup(dt_inner, Q_x, Q_u, Nh, gamma, u_max);
    switch wind_model
        case 'none',         wind_fn = @(tt,p) zeros(size(p));
        case 'uniform',      wind_fn = @(tt,p) uniform_wind(tt, p, wind_U, wind_theta);
        case 'taylor_green', wind_fn = @(tt,p) taylor_green_wind(tt, p, wind_U, wind_L);
    end
end

%% Simulation
pos = params.pos0;          % agent positions
vel = zeros(N, 2);          % agent velocities (MPC only)
d_hat = zeros(N, 2);        % agents' wind estimates ('estimated' only)
n_active = 0;               % per-axis MPC solves with an active bound (u_max only)
alpha = zeros(N, 1);        % potential coefficients, warm start for the next step
lambda = zeros(N, 1);       % speed-limit multipliers, warm start for the next step
pos_hist = zeros(N, 2, T);
alpha_hist = zeros(N, T);
min_rho_terminal_hist = zeros(T, 1);

for t = 1:T
    % Target: the star rotated by target_omega per step
    rho_target = star_target_w(params, params.target_omega * (t - 1));

    % DOOT: velocity command for every agent
    [v_doot, alpha, lambda] = plan_doot(pos, rho_target, alpha, lambda, params);

    if use_mpc
        % Low-level MPC tracks the reference p_start + tau*v_doot
        p_start = pos;
        for s = 1:n_inner
            tau = (s - 1) * dt_inner;
            time_now = (t - 1) * dt_outer + tau;
            p_ref = p_start + tau * v_doot;
            d = wind_fn(time_now, pos);     % true wind at every agent

            % Wind the MPC compensates for; the agents always drift with the true d
            switch drift_mode
                case 'measured',  d_mpc = d;
                case 'estimated', d_mpc = d_hat;
                case 'none',      d_mpc = zeros(N, 2);
                otherwise, error('drift_mode must be ''measured'', ''estimated'' or ''none''.');
            end
            [u, n_act] = low_level_mpc(pos, vel, p_ref, v_doot, d_mpc, mpc);
            n_active = n_active + n_act;

            % Exact step of xdot = v + d, vdot = u (u and d held over the step)
            pos_pred = pos + dt_inner * vel + 0.5 * dt_inner^2 * u;   % without wind
            pos = pos + dt_inner * (vel + d) + 0.5 * dt_inner^2 * u;
            vel = vel + dt_inner * u;

            % Estimate the wind from the gap between the measured and the
            % predicted position: needs self-localization, not a wind sensor
            if strcmp(drift_mode, 'estimated')
                d_hat = (1 - drift_alpha) * d_hat + drift_alpha * (pos - pos_pred) / dt_inner;
            end
        end
    else
        pos = pos + dt_outer * v_doot;
    end

    pos_hist(:, :, t) = pos;
    alpha_hist(:, t) = alpha;

    % Check the signed terminal density on the complete spatial grid.
    sq_dists_grid = pdist2(params.grid_pts, pos).^2;
    phi_N = wendland_kernel(sq_dists_grid, params.h_field) * alpha;
    min_rho_terminal_hist(t) = min(rho_target + params.epsilon_reg * phi_N);
    fprintf('Step %d/%d: min(rho^* + epsilon phi_N) = %.6e%s\n', ...
        t, T, min_rho_terminal_hist(t), ...
        ternary(min_rho_terminal_hist(t) < 0, '  NEGATIVE', ''));

    if mod(t, 10) == 0 || t == T
        fprintf('Step %d/%d complete.\n', t, T);
    end
end

if use_mpc && ~isempty(u_max)
    % 0% means the bound never acted, so the run equals the unconstrained one
    fprintf('Acceleration bound active in %.2f%% of the per-axis MPC solves.\n', ...
        100 * n_active / (2 * N * n_inner * T));
end

%% Video of the run
plotresults_w(pos_hist, alpha_hist, rho_target, params, min_rho_terminal_hist);

%% DOOT planner: velocity v = -grad phi, phi = sum_j alpha_j K_field(x - x_j)
function [v, alpha, lambda] = plan_doot(pos, rho_target, alpha, lambda, params)
    kernels = computeKernels_w(pos, params);
    integral = computeTargetint_w(pos, rho_target, params);

    % Density mismatch w_j = int K_field(x - x_j) (rho_N(x) - rho_target(x)) dx
    w = sum(integral.agents_ag, 2) / params.N - integral.target_ag;

    [alpha, lambda] = primal_dual_2(kernels, integral, w, params, alpha, lambda);
    v = -[kernels.D_x * alpha, kernels.D_y * alpha];
end

function out = ternary(condition, true_value, false_value)
if condition
    out = true_value;
else
    out = false_value;
end
end
