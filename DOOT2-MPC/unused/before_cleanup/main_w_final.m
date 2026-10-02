clc; close all;
%% DOOT + low-level MPC only.
% Run parameters_w.m before this script.

if ~exist('params', 'var') || ~isfield(params, 'pos_cont0')
    error('Run parameters_w.m before main_w_final.m.');
end

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(fileparts(script_dir), 'MPC'));

%% Settings
wind_model = 'none';  % 'uniform', 'taylor_green', or 'none'
wind_U = 0.01;
wind_theta = 0;
wind_L = 6;

warm_start = true;
warm_start_lambda = true;
n_inner = 10;
Q_x = 1;
Q_u = 3.5;
Nh = 10;
gamma = 1;
u_max = [];

N = params.N;
T = params.T;
dt_outer = params.dt_cont;
dt_inner = dt_outer / n_inner;

assert(T >= 2 && dt_outer > 0 && n_inner >= 1 && mod(n_inner, 1) == 0, ...
    'Need T >= 2, positive dt_cont, and a positive integer n_inner.');

mpc = low_level_mpc_setup(dt_inner, Q_x, Q_u, Nh, gamma, u_max);
%mpc = low_level_mpc_setup_qp(dt_inner, Q_x, Q_u, Nh, gamma, u_max);

switch lower(wind_model)
    case 'none'
        wind_fn = @(tt,p) zeros(size(p));
    case 'uniform'
    addpath(fullfile(fileparts(script_dir), 'WINDS'));
        wind_fn = @(tt,p) uniform_wind(tt, p, wind_U, wind_theta);
    case 'taylor_green'
        addpath(fullfile(fileparts(script_dir), 'WINDS'));
        wind_fn = @(tt,p) taylor_green_wind(tt, p, wind_U, wind_L);
    otherwise
        error('dootMpc:windModel', ...
            'wind_model must be ''uniform'', ''taylor_green'', or ''none''.');
end

targetParams = params;
target_rot = isfield(params, 'target_rot') && params.target_rot;
if target_rot
    targetParams.rho_target_sm = [];
end

%% DOOT + MPC state
pos = params.pos_cont0;
vel = zeros(N, 2);
alpha = zeros(params.N_tot + params.M, 1);
lambda = zeros(params.N_tot + params.M, 1);

pos_hist = zeros(params.N_tot, 2, T + 1);
pos_hist(:, :, 1) = pos;
alpha_hist = zeros(params.N_tot + params.M, T);

fprintf('Running DOOT + low-level MPC only.\n');
fprintf('dt_outer=%g, dt_inner=%g, Nh=%d, Q_u=%g\n', ...
    dt_outer, dt_inner, Nh, Q_u);

for t = 1:T
    if target_rot
        targetParams.rho_target = ...
            star_target_w(params, params.target_omega * (t - 1));
    end

    %DOOT 
    [v_doot, alpha, lambda] = plan_doot( pos, alpha, lambda, targetParams,warm_start, warm_start_lambda);
    alpha_hist(:, t) = alpha;
    p_start = pos(1:N, :);

    %{
    for s = 1:n_inner
        tau = (s - 1) * dt_inner;
        time_now = (t - 1) * dt_outer + tau;
        p_ref = p_start + tau * v_doot;

        for m = 1:N
            p_old = pos(m, :);
            v_old = vel(m, :);
            [u, v_next, ~, info] = low_level_mpc(p_old, v_old, p_ref(m, :), v_doot(m, :), time_now, wind_fn, mpc);
            %[u, v_next, ~, info] = low_level_mpc_qp(p_old, v_old, p_ref(m, :), v_doot(m, :),time_now, wind_fn, mpc);

            vel(m, :) = v_next;
            pos(m, :) = p_old + dt_inner * (v_old + info.d) + ...
                0.5 * dt_inner^2 * u;
        end
    end
    %}

    pos(1:N, :) = pos(1:N, :) + dt_outer * v_doot;
    pos_hist(:, :, t + 1) = pos;

    if mod(t, 10) == 0 || t == T
        fprintf('Step %d/%d complete.\n', t, T);
    end
end

%% Plot final DOOT + MPC result
plotData_mpc = struct( ...
    'X', params.X, ...
    'Y', params.Y, ...
    'rho_target_mat', reshape(targetParams.rho_target, ...
                              params.grid_res, params.grid_res), ...
    'pos_hist', pos_hist(:, :, 2:end), ...
    'alpha_hist', alpha_hist, ...
    'grid_pts', params.grid_pts, ...
    'grid_res', params.grid_res, ...
    'N', params.N, ...
    'dt_cont', params.dt_cont, ...
    'h_agent', params.h_agent, ...
    'h_field', params.h_field, ...
    'h_ant', params.h_ant, ...
    'ant_cords', params.ant_cords);

plotresults_w('video', plotData_mpc);

%% Local helper
function [v, alpha, lambda] = plan_doot( ...
        pos, alpha_prev, lambda_prev, params, ws, wsl)
    N = params.N;
    N_tot = params.N_tot;

    kernels = computeKernels_w(pos, params.ant_cords, params);
    integral = computeTargetint_w(pos, params.ant_cords, params);

    w_agent = zeros(N_tot, 1);
    w_agent(1:N) = sum(integral.agents_ag(1:N, 1:N), 2) / N - ...
        integral.target_ag(1:N);
    w_agent(N+1:N_tot) = ...
        sum(integral.agents_ag(N+1:N_tot, N+1:N_tot), 2) + ...
        params.ghost_w;

    if ~isempty(integral.agents_ant)
        int_agents_ant = sum(integral.agents_ant, 1)' / N;
        w_ant = params.antw * ...
            (int_agents_ant / N - integral.target_ant);
    else
        w_ant = [];
    end

    w_cont = params.w_gain * [w_agent; w_ant];
    w_cont(N+1:N_tot) = w_cont(N+1:N_tot) + params.ghost_w;
    c_sq = conformal(pos, params.ant_cords, params);

    if ws && wsl
        [alpha, lambda] = primal_dual_2( ...
            kernels, integral, w_cont, c_sq, params, ...
            alpha_prev, lambda_prev);
    elseif ws
        [alpha, lambda] = primal_dual_2( ...
            kernels, integral, w_cont, c_sq, params, alpha_prev);
    else
        [alpha, lambda] = primal_dual_2( ...
            kernels, integral, w_cont, c_sq, params);
    end

    v = -[kernels.D_x(1:N, :) * alpha, ...
          kernels.D_y(1:N, :) * alpha];
end
