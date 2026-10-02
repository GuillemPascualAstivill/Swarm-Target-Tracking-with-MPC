%% paper_figure_w.m  --  the 2x2 paper figure, exported as a square PDF:
%   (a) control effort over time for several l2 regularizations eps
%   (b) W1 distance to the target over time in Taylor-Green wind, with the
%       low-level MPC and without low level
%   (c), (d) final agent density with eps = eps_star and with eps = 0
% Curves are means over independent initial swarms, with a light min-max
% band. Runs parameters_w.m itself: it calls clear, so it must come first.
%
% W1 at every step of (b) dominates the cost (~6 s per W1 on one worker):
% ~2.5 h on 20 workers. The runs are saved to results_file before plotting;
% run_sims = false only redraws the figure from that file.

parameters_w;
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(fileparts(script_dir), 'WINDS'), ...
    fullfile(fileparts(script_dir), 'MPC'));
results_dir = fullfile(script_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% Settings
run_sims = false;          % false: redraw the figure from results_file
results_file = fullfile(results_dir, 'paper_figure_results.mat');
seeds = 1:6;               % initial swarms N(0, 3 I); seeds(1) is drawn in (c), (d)

% (a) regularization: MPC in weak Taylor-Green wind
eps_list = [0 0.05 0.25];  % l2 regularization eps
eps_star = 0.25;           % eps of (b) and (c), drawn thicker in (a)
U_reg = 0.01;              % wind speed of (a), (c), (d)

% (b) wind disturbance
winds = [0.03 0.05 0.1];   % Taylor-Green peak wind speeds U
W1_every = 1;              % W1 every k DOOT steps
wind_L = 6;                % vortex period

% Low level: MPC weights as in main_w_final.m. Without low level the agents
% are the same double integrators with the velocity loop u = -k_v (v - v_doot)
% (MPC's k_v, no position reference, blind to the wind)
n_inner = 10;
Q_x = 1;
Q_u = 3.5;
Nh = 10;
gamma = 1;

% Figure
smooth_win = 10;           % moving average of the curves [steps]
grid_cm = 13;              % width and height of the square PDF
font_pt = 10;
pdf_file = fullfile(results_dir, 'paper_figure_w');
quick_test = false;        % true: 5-step runs, 2 seeds, to check the script

%% Runs
if run_sims
    if quick_test, params.T = 5; seeds = 1:2; end
    mpc = low_level_mpc_setup(params.dt_cont / n_inner, Q_x, Q_u, Nh, gamma, []);

    % (b) first: the W1 runs are the long ones
    specs = struct('panel', {}, 'eps', {}, 'U', {}, 'use_mpc', {}, 'seed', {});
    for use_mpc = [true false]
        for U = winds
            for s = seeds
                specs(end + 1) = struct('panel', 'b', 'eps', eps_star, 'U', U, ...
                    'use_mpc', use_mpc, 'seed', s); %#ok<SAGROW>
            end
        end
    end
    for e = eps_list
        for s = seeds
            specs(end + 1) = struct('panel', 'a', 'eps', e, 'U', U_reg, ...
                'use_mpc', true, 'seed', s); %#ok<SAGROW>
        end
    end
    fprintf('%d runs.\n', numel(specs));

    t_all = tic;
    res = cell(numel(specs), 1);
    one_at_a_time = parforOptions(gcp, 'RangePartitionMethod', 'fixed', 'SubrangeSize', 1);
    parfor (r = 1:numel(specs), one_at_a_time)
        res{r} = run_fig(params, specs(r), wind_L, n_inner, mpc, W1_every);
    end
    fprintf('All runs done in %.1f h.\n', toc(t_all) / 3600);
    save(results_file, 'specs', 'res', 'seeds', 'eps_list', 'eps_star', 'U_reg', ...
        'winds', 'W1_every', 'wind_L', 'n_inner', 'Q_x', 'Q_u', 'Nh', 'gamma', 'params');
else
    load(results_file, 'specs', 'res', 'seeds', 'eps_list', 'eps_star', 'U_reg', ...
        'winds', 'W1_every', 'params');
end

%% Collect: effort(t, eps, seed), W1(t, U, seed, 1 = MPC / 2 = no low level)
T = params.T;
ne = numel(eps_list);
nU = numel(winds);
nk = numel(seeds);
effort = zeros(T, ne, nk);
W1 = nan(T, nU, nk, 2);
pos_snap = cell(ne, 1);
for r = 1:numel(specs)
    sp = specs(r);
    k = find(seeds == sp.seed);
    if sp.panel == 'a'
        j = find(eps_list == sp.eps);
        effort(:, j, k) = res{r}.effort;
        if k == 1, pos_snap{j} = res{r}.pos; end
    else
        W1(:, winds == sp.U, k, 2 - sp.use_mpc) = res{r}.W1;
    end
end
t_axis = (1:T)' * params.dt_cont;
t_W1 = (W1_every:W1_every:T)';
W1 = W1(t_W1, :, :, :);

tail_idx = max(1, T - 49):T;
fprintf('\n(a) mean over %d seeds of the last %d steps\n   eps     ||u||^2\n', nk, numel(tail_idx));
for j = 1:ne
    fprintf('%6.3f   %.3e\n', eps_list(j), mean(effort(tail_idx, j, :), 'all'));
end
tail_W1 = t_W1 >= tail_idx(1);
fprintf('\n(b) W1, mean over %d seeds\n', nk);
fprintf('            final            mean, last %d steps\n', numel(tail_idx));
fprintf('   U      MPC   no low level      MPC   no low level\n');
for k = 1:nU
    fprintf('%6.2f   %.4f     %.4f       %.4f     %.4f\n', winds(k), ...
        mean(W1(end, k, :, 1)), mean(W1(end, k, :, 2)), ...
        mean(W1(tail_W1, k, :, 1), 'all'), mean(W1(tail_W1, k, :, 2), 'all'));
end

%% Figure: four equal square panels
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 grid_cm grid_cm]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
lw = 1.4;
axs = gobjects(4, 1);
titles = {'(a) regularization', '(b) wind disturbance', ...
    sprintf('(c) $\\varepsilon = %g$', eps_star), '(d) $\varepsilon = 0$'};

% (a) effort per eps: hue = eps, eps_star thicker
ax = nexttile(tl); axs(1) = ax;
hold(ax, 'on'); box(ax, 'on');
cols = parula(ne + 1);
cols = cols(1:ne, :);
Y = movmean(effort, smooth_win, 1);
for j = 1:ne                                 % bands first, lines on top
    band(ax, t_axis, squeeze(Y(:, j, :)), 0.3 * cols(j, :) + 0.7);
end
for j = 1:ne
    w = lw * (1 + (eps_list(j) == eps_star));
    plot(ax, t_axis, mean(Y(:, j, :), 3), '-', 'Color', cols(j, :), 'LineWidth', w, ...
        'DisplayName', sprintf('$\\varepsilon = %g$', eps_list(j)));
end
set(ax, 'YScale', 'log');
xlim(ax, [0, t_axis(end)]);
xlabel(ax, 'time');
ylabel(ax, '$\|u\|^2$');
lg = legend(ax, 'Location', 'east');
lg.FontSize = font_pt - 2;
lg.ItemTokenSize = [8 6];

% (b) W1 per wind: hue = controller, lightness = wind speed
ax = nexttile(tl); axs(2) = ax;
hold(ax, 'on'); box(ax, 'on');
c_mpc = [0 0.45 0.74];
c_nll = [0.85 0.33 0.1];
a = linspace(0.4, 1, nU);                    % colour strength per wind
tint = @(c, ak) ak * c + (1 - ak) * [1 1 1];
Y = movmean(W1, max(1, round(smooth_win / W1_every)), 1);
for k = 1:nU
    band(ax, t_W1 * params.dt_cont, squeeze(Y(:, k, :, 2)), tint(c_nll, 0.25));
    band(ax, t_W1 * params.dt_cont, squeeze(Y(:, k, :, 1)), tint(c_mpc, 0.25));
end
for k = 1:nU                                 % no low level below, MPC on top
    plot(ax, t_W1 * params.dt_cont, mean(Y(:, k, :, 2), 3), '--', 'Color', tint(c_nll, a(k)), ...
        'LineWidth', 1.2 * lw, 'HandleVisibility', 'off');
end
for k = 1:nU
    plot(ax, t_W1 * params.dt_cont, mean(Y(:, k, :, 1), 3), '-', 'Color', tint(c_mpc, a(k)), ...
        'LineWidth', 1.2 * lw, 'HandleVisibility', 'off');
end
% Legend: controllers by hue and line style, winds by lightness (grey)
plot(ax, nan, nan, '-', 'Color', c_mpc, 'LineWidth', 1.2 * lw, 'DisplayName', 'MPC');
plot(ax, nan, nan, '--', 'Color', c_nll, 'LineWidth', 1.2 * lw, 'DisplayName', 'no low level');
for k = 1:nU
    plot(ax, nan, nan, '-', 'Color', tint([0 0 0], a(k)), 'LineWidth', 1.2 * lw, ...
        'DisplayName', sprintf('$U = %g$', winds(k)));
end
xlim(ax, [0, t_axis(end)]);
xlabel(ax, 'time');
ylabel(ax, '$W_1(\rho_N, \rho^*)$');
lg = legend(ax, 'Location', 'northeast');
lg.FontSize = font_pt - 2;
lg.ItemTokenSize = [10 6];

% (c), (d) final agent density of seeds(1), one colour scale capped at the
% peak of the final target smoothed with the same kernel (clumps saturate)
tgt = star_target_w(params, params.target_omega * (T - 1));
cmax = max(smooth_on_grid_wendland(tgt, params, params.h_agent));
j_snap = [find(eps_list == eps_star), find(eps_list == 0)];
for q = 1:2
    axs(2 + q) = nexttile(tl);
    draw_snap(axs(2 + q), pos_snap{j_snap(q)}, params, cmax);
end

for q = 1:4
    if q <= 2                    % (c), (d) are square already from axis equal;
        axis(axs(q), 'square');  % square on top breaks their images in the PDF
    end
    set(axs(q), 'FontSize', font_pt, 'TickLabelInterpreter', 'latex');
    title(axs(q), titles{q}, 'FontSize', font_pt + (q >= 3));
end

% Vector PDF of the whole figure, grid_cm x grid_cm: fonts keep their point
% size. Not print -vector: it renders the density images of (c), (d) flat
exportgraphics(fig, [pdf_file '.pdf'], 'ContentType', 'vector', 'Padding', 'figure');
exportgraphics(fig, [pdf_file '.png'], 'Resolution', 200);
fprintf('Saved %s.pdf (%g x %g cm).\n', pdf_file, grid_cm, grid_cm);

%% One run of main_w_final.m in Taylor-Green wind of peak speed sp.U
%   out.effort(t)  mean over agents of |u|^2 over DOOT step t
%   out.W1(t)      entropic W1 between the agent KDE and the target, both
%                  smoothed at h_agent, after step t (panel 'b' only)
%   out.pos        final agent positions
function out = run_fig(params, sp, wind_L, n_inner, mpc, W1_every)
    t_run = tic;
    N = params.N;
    T = params.T;
    dt_inner = params.dt_cont / n_inner;
    params.epsilon_reg = sp.eps;
    x = params.X(1, :);
    y = params.Y(:, 1);

    pos = sqrt(3) * randn(RandStream('mt19937ar', 'Seed', sp.seed), N, 2);
    vel = zeros(N, 2);
    d_hat = zeros(N, 2);
    alpha = zeros(N, 1);
    lambda = zeros(N, 1);
    out.effort = zeros(T, 1);
    out.W1 = nan(T, 1);

    for t = 1:T
        rho_target = star_target_w(params, params.target_omega * (t - 1));
        [v_doot, alpha, lambda] = plan_doot(pos, rho_target, alpha, lambda, params);

        p_start = pos;
        for s = 1:n_inner
            d = taylor_green_wind(0, pos, sp.U, wind_L);   % true wind (steady field)
            if sp.use_mpc
                % Track p_start + tau*v_doot, compensating the estimated wind
                p_ref = p_start + (s - 1) * dt_inner * v_doot;
                u = low_level_mpc(pos, vel, p_ref, v_doot, d_hat, mpc);
            else
                u = -mpc.K(2) * (vel - v_doot);
            end
            out.effort(t) = out.effort(t) + sum(u(:).^2) / (N * n_inner);

            % Exact step of xdot = v + d, vdot = u (u and d held over the step)
            pos_pred = pos + dt_inner * vel + 0.5 * dt_inner^2 * u;
            pos = pos + dt_inner * (vel + d) + 0.5 * dt_inner^2 * u;
            vel = vel + dt_inner * u;
            if sp.use_mpc
                d_hat = 0.8 * d_hat + 0.2 * (pos - pos_pred) / dt_inner;
            end
        end

        if sp.panel == 'b' && mod(t, W1_every) == 0
            sq_dists = (params.grid_pts(:, 1) - pos(:, 1)').^2 + (params.grid_pts(:, 2) - pos(:, 2)').^2;
            rho_N = sum(wendland_kernel(sq_dists, params.h_agent), 2);
            tgt_sm = smooth_on_grid_wendland(rho_target, params, params.h_agent);
            out.W1(t) = computeW2(rho_N, tgt_sm, x, y, 'p', 1);
        end
    end
    out.pos = pos;

    if sp.use_mpc, name = 'MPC'; else, name = 'no low level'; end
    fprintf('(%s) eps = %g, U = %g, %s, seed %d: mean ||u||^2 = %.3e, final W1 = %.4f (%.1f min)\n', ...
        sp.panel, sp.eps, sp.U, name, sp.seed, mean(out.effort), out.W1(T), toc(t_run) / 60);
end

%% DOOT planner, as in main_w_final.m
function [v, alpha, lambda] = plan_doot(pos, rho_target, alpha, lambda, params)
    kernels = computeKernels_w(pos, params);
    integral = computeTargetint_w(pos, rho_target, params);
    w = sum(integral.agents_ag, 2) / params.N - integral.target_ag;
    [alpha, lambda] = primal_dual(kernels, integral, w, params, alpha, lambda);
    v = -[kernels.D_x * alpha, kernels.D_y * alpha];
end

%% Min-max band of the columns of Y over t
function band(ax, t, Y, col)
    fill(ax, [t; flipud(t)], [min(Y, [], 2); flipud(max(Y, [], 2))], col, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end

%% Agent density (Wendland KDE at h_agent) as a smooth image, agents on top
function draw_snap(ax, P, params, cmax)
    hold(ax, 'on');
    x = params.X(1, :);
    y = params.Y(:, 1)';
    sq_dists = (params.grid_pts(:, 1) - P(:, 1)').^2 + (params.grid_pts(:, 2) - P(:, 2)').^2;
    kde = sum(wendland_kernel(sq_dists, params.h_agent), 2) / size(P, 1);
    % An image, not contourf: vector PDFs show seams between filled contour
    % bands, an image is embedded as one bitmap
    Z = reshape(min(kde, cmax), params.grid_res, params.grid_res);
    xf = linspace(x(1), x(end), 300);
    yf = linspace(y(1), y(end), 300);
    [XF, YF] = meshgrid(xf, yf);
    Zf = min(max(interp2(params.X, params.Y, Z, XF, YF, 'cubic'), 0), cmax);
    imagesc(ax, xf, yf, Zf);
    set(ax, 'YDir', 'normal');
    colormap(ax, parula);
    clim(ax, [0 cmax]);
    scatter(ax, P(:, 1), P(:, 2), 1, [0.9 0.9 0.9], 'filled');
    axis(ax, 'equal');
    set(ax, 'XLim', [x(1) x(end)], 'YLim', [y(1) y(end)], 'XTick', [-5 0 5], 'YTick', [-5 0 5]);
    box(ax, 'on');
end
