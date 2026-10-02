%% regularization_sweep_w.m  --  control effort over time for several values of
% the l2 regularization eps, each averaged over independent initial swarms
% (the regularization panel (a) of the paper figure). Runs parameters_w.m
% itself: it calls clear, so it must come first.

parameters_w;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(fileparts(script_dir), 'WINDS'), ...
    fullfile(fileparts(script_dir), 'MPC'));
results_dir = fullfile(script_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% Settings
eps_list = [0 0.05 0.25];  % l2 regularization eps (params.epsilon_reg)
eps_star = 0.25;           % eps drawn thicker (value used in the paper)
seeds = 1:6;               % initial swarms, N(0, 3 I) as in parameters_w.m
wind_U = 0.01;             % Taylor-Green wind, compensated by the MPC (0 = none)
wind_L = 6;
n_inner = 10;              % MPC steps per DOOT step
Q_x = 1;                   % MPC weights, as in main_w_final.m
Q_u = 3.5;
Nh = 10;
gamma = 1;
smooth_win = 10;           % moving average of the plotted curves [steps]
tail = 50;                 % last steps averaged for the steady-state values
quick_test = false;        % true: 5-step runs to check the script

params.terminal = 'l2';
if quick_test, params.T = 5; end
T = params.T;
mpc = low_level_mpc_setup(params.dt_cont / n_inner, Q_x, Q_u, Nh, gamma, []);

%% Runs: r = (k - 1) * ne + j is eps_list(j) with seeds(k)
ne = numel(eps_list);
nk = numel(seeds);
runs = cell(ne * nk, 1);
parfor r = 1:ne * nk
    [j, k] = ind2sub([ne, nk], r);
    p = params;
    p.epsilon_reg = eps_list(j);
    p.pos0 = sqrt(3) * randn(RandStream('mt19937ar', 'Seed', seeds(k)), p.N, 2);
    [effort, L1] = run_reg(p, wind_U, wind_L, n_inner, mpc);
    runs{r} = [effort, L1];
end
R = reshape(cell2mat(runs'), T, 2, ne, nk);
effort = squeeze(R(:, 1, :, :));           % T x ne x nk
L1 = squeeze(R(:, 2, :, :));
save(fullfile(results_dir, 'regularization_sweep_results.mat'), 'effort', 'L1', 'eps_list', 'seeds', ...
    'wind_U', 'Q_x', 'Q_u', 'params');

%% Steady state, mean over seeds
tail_idx = max(1, T - tail + 1):T;
fprintf('\nMean over %d seeds of the last %d steps:\n', nk, numel(tail_idx));
fprintf('   eps     ||u||^2      L1\n');
for j = 1:ne
    fprintf('%6.3f   %.3e   %.4f\n', eps_list(j), mean(effort(tail_idx, j, :), 'all'), ...
        mean(L1(tail_idx, j, :), 'all'));
end

%% Figure: seed mean, light band = seed min-max
cols = parula(ne + 1);
cols = cols(1:ne, :);
lw = 1.7;
t_axis = (1:T)' * params.dt_cont;

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 11 10]);
hold on; box on;
Y = movmean(effort, smooth_win, 1);
for j = 1:ne                                 % bands first, lines on top
    Yj = squeeze(Y(:, j, :));
    fill([t_axis; flipud(t_axis)], [min(Yj, [], 2); flipud(max(Yj, [], 2))], ...
        0.3 * cols(j, :) + 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
for j = 1:ne
    w = lw * (1 + (abs(eps_list(j) - eps_star) < 1e-12));
    plot(t_axis, mean(Y(:, j, :), 3), '-', 'Color', cols(j, :), 'LineWidth', w, ...
        'DisplayName', sprintf('$\\varepsilon = %g$', eps_list(j)));
end
set(gca, 'YScale', 'log', 'FontSize', 13);
xlim([0, t_axis(end)]);
xlabel('time');
ylabel('$\|u\|^2$');
title('(a) regularization');
legend('Location', 'east');
exportgraphics(fig, fullfile(results_dir, 'regularization_sweep_w.png'), 'Resolution', 200);

%% One run of main_w_final.m with the MPC in Taylor-Green wind
% effort(t): mean over agents of |u|^2, averaged over DOOT step t
% L1(t):     ||rho_N - rho*||_1 after step t (Wendland KDE at h_agent)
function [effort, L1] = run_reg(params, wind_U, wind_L, n_inner, mpc)
    t_run = tic;
    N = params.N;
    T = params.T;
    dt_inner = params.dt_cont / n_inner;
    dA = params.dx * params.dy;

    pos = params.pos0;
    vel = zeros(N, 2);
    alpha = zeros(N, 1);
    lambda = zeros(N, 1);
    effort = zeros(T, 1);
    L1 = zeros(T, 1);

    for t = 1:T
        rho_target = star_target_w(params, params.target_omega * (t - 1));
        [v_doot, alpha, lambda] = plan_doot(pos, rho_target, alpha, lambda, params);

        p_start = pos;
        for s = 1:n_inner
            d = taylor_green_wind(0, pos, wind_U, wind_L);   % true wind (steady field)
            p_ref = p_start + (s - 1) * dt_inner * v_doot;
            u = low_level_mpc(pos, vel, p_ref, v_doot, d, mpc);
            effort(t) = effort(t) + sum(u(:).^2) / (N * n_inner);

            % Exact step of xdot = v + d, vdot = u (u and d held over the step)
            pos = pos + dt_inner * (vel + d) + 0.5 * dt_inner^2 * u;
            vel = vel + dt_inner * u;
        end

        sq_dists = (params.grid_pts(:, 1) - pos(:, 1)').^2 + (params.grid_pts(:, 2) - pos(:, 2)').^2;
        rho_N = sum(wendland_kernel(sq_dists, params.h_agent), 2);
        rho_N = rho_N / (sum(rho_N) * dA);
        L1(t) = sum(abs(rho_N - rho_target)) * dA;
    end

    fprintf('eps = %.3f: mean ||u||^2 = %.3e, final L1 = %.4f (%.1f min)\n', ...
        params.epsilon_reg, mean(effort), L1(T), toc(t_run) / 60);
end

%% DOOT planner, as in main_w_final.m
function [v, alpha, lambda] = plan_doot(pos, rho_target, alpha, lambda, params)
    kernels = computeKernels_w(pos, params);
    integral = computeTargetint_w(pos, rho_target, params);
    w = sum(integral.agents_ag, 2) / params.N - integral.target_ag;
    [alpha, lambda] = primal_dual(kernels, integral, w, params, alpha, lambda);
    v = -[kernels.D_x * alpha, kernels.D_y * alpha];
end
