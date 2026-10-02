%% wind_sweep_w.m  --  L1 error under Taylor-Green wind, with and without the
% low-level MPC (the wind panel (b) of the paper figure). Runs parameters_w.m
% itself: it calls clear, so it must come first.

parameters_w;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(fileparts(script_dir), 'WINDS'), ...
    fullfile(fileparts(script_dir), 'MPC'));
results_dir = fullfile(script_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% Settings
winds = [0.03 0.05 0.1];   % Taylor-Green peak wind speeds U
wind_L = 6;                % vortex period
no_ll = 'vel';             % agents without the low level:
                           % 'vel':   same double integrator, velocity loop
                           %          u = -k_v (v - v_doot) with the MPC's k_v,
                           %          no position reference, blind to the wind
                           % 'first': first order, xdot = v_doot + d
n_inner = 10;              % low-level steps per DOOT step
Q_x = 1;                   % MPC weights, as in main_w_final.m; the MPC
Q_u = 3.5;                 % compensates the measured wind
Nh = 10;
gamma = 1;
smooth_win = 10;           % moving average of the plotted curves [steps]
tail = 50;                 % last steps averaged for the steady-state L1
quick_test = false;        % true: 5-step runs to check the script

if quick_test, params.T = 5; end
T = params.T;
mpc = low_level_mpc_setup(params.dt_cont / n_inner, Q_x, Q_u, Nh, gamma, []);

%% Runs: r = 1..nU with the MPC, nU+1..2nU without low level
nU = numel(winds);
runs = cell(2 * nU, 1);
parfor r = 1:2 * nU
    [k, c] = ind2sub([nU, 2], r);
    runs{r} = run_wind(params, c == 1, no_ll, winds(k), wind_L, n_inner, mpc);
end
L1 = reshape(cell2mat(runs'), T, nU, 2);     % L1(:, k, 1) MPC, L1(:, k, 2) no low level
save(fullfile(results_dir, 'wind_sweep_results.mat'), 'L1', 'winds', 'wind_L', 'no_ll', 'Q_x', 'Q_u', 'params');

%% Final distance
tail_idx = max(1, T - tail + 1):T;
fprintf('\n             final L1             mean L1, last %d steps\n', numel(tail_idx));
fprintf('   U      MPC   no low level      MPC   no low level\n');
for k = 1:nU
    fprintf('%6.2f   %.4f     %.4f       %.4f     %.4f\n', winds(k), L1(T, k, 1), L1(T, k, 2), ...
        mean(L1(tail_idx, k, 1)), mean(L1(tail_idx, k, 2)));
end

%% Figure: hue = controller, lightness = wind speed
c_mpc = [0 0.45 0.74];
c_nll = [0.85 0.33 0.1];
a = linspace(0.4, 1, nU);                    % colour strength per wind
tint = @(c, ak) ak * c + (1 - ak) * [1 1 1];
lw = 1.7;
t_axis = (1:T)' * params.dt_cont;

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 11 10]);
hold on; box on;
for k = 1:nU                                 % no low level below, MPC on top
    plot(t_axis, movmean(L1(:, k, 2), smooth_win), '--', 'Color', tint(c_nll, a(k)), ...
        'LineWidth', lw, 'HandleVisibility', 'off');
end
for k = 1:nU
    plot(t_axis, movmean(L1(:, k, 1), smooth_win), '-', 'Color', tint(c_mpc, a(k)), ...
        'LineWidth', lw, 'HandleVisibility', 'off');
end
% Legend: controllers by hue and line style, winds by lightness (grey)
plot(nan, nan, '-', 'Color', c_mpc, 'LineWidth', lw, 'DisplayName', 'MPC');
plot(nan, nan, '--', 'Color', c_nll, 'LineWidth', lw, 'DisplayName', 'no low level');
for k = 1:nU
    plot(nan, nan, '-', 'Color', tint([0 0 0], a(k)), 'LineWidth', lw, ...
        'DisplayName', sprintf('$U = %g$', winds(k)));
end
xlim([0, t_axis(end)]);
ylabel('$\|\rho_N-\rho^*\|_{L^1}$');
title('(b) wind disturbance');
legend('Location', 'northeast');
set(gca, 'FontSize', 13);
exportgraphics(fig, fullfile(results_dir, 'wind_sweep_w.png'), 'Resolution', 200);

%% One run of main_w_final.m in Taylor-Green wind of peak speed U
function L1 = run_wind(params, use_mpc, no_ll, U, wind_L, n_inner, mpc)
    t_run = tic;
    N = params.N;
    T = params.T;
    dt_inner = params.dt_cont / n_inner;
    dA = params.dx * params.dy;

    pos = params.pos0;
    vel = zeros(N, 2);
    alpha = zeros(N, 1);
    lambda = zeros(N, 1);
    L1 = zeros(T, 1);

    for t = 1:T
        rho_target = star_target_w(params, params.target_omega * (t - 1));
        [v_doot, alpha, lambda] = plan_doot(pos, rho_target, alpha, lambda, params);

        p_start = pos;
        for s = 1:n_inner
            d = taylor_green_wind(0, pos, U, wind_L);    % true wind (steady field)
            if use_mpc
                % Track p_start + tau*v_doot, compensating the measured wind
                p_ref = p_start + (s - 1) * dt_inner * v_doot;
                u = low_level_mpc(pos, vel, p_ref, v_doot, d, mpc);
            elseif strcmp(no_ll, 'vel')
                u = -mpc.K(2) * (vel - v_doot);
            else
                pos = pos + dt_inner * (v_doot + d);
                continue
            end
            % Exact step of xdot = v + d, vdot = u (u and d held over the step)
            pos = pos + dt_inner * (vel + d) + 0.5 * dt_inner^2 * u;
            vel = vel + dt_inner * u;
        end

        % L1 error after the move: Wendland KDE at h_agent, unit mass on the grid
        sq_dists = (params.grid_pts(:, 1) - pos(:, 1)').^2 + (params.grid_pts(:, 2) - pos(:, 2)').^2;
        rho_N = sum(wendland_kernel(sq_dists, params.h_agent), 2);
        rho_N = rho_N / (sum(rho_N) * dA);
        L1(t) = sum(abs(rho_N - rho_target)) * dA;
    end

    if use_mpc, name = 'MPC'; else, name = ['no low level (' no_ll ')']; end
    fprintf('U = %.2f, %s: final L1 = %.4f (%.1f min)\n', U, name, L1(T), toc(t_run) / 60);
end

%% DOOT planner, as in main_w_final.m
function [v, alpha, lambda] = plan_doot(pos, rho_target, alpha, lambda, params)
    kernels = computeKernels_w(pos, params);
    integral = computeTargetint_w(pos, rho_target, params);
    w = sum(integral.agents_ag, 2) / params.N - integral.target_ag;
    [alpha, lambda] = primal_dual(kernels, integral, w, params, alpha, lambda);
    v = -[kernels.D_x * alpha, kernels.D_y * alpha];
end
