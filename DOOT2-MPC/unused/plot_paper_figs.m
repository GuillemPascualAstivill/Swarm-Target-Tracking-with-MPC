function figs = plot_paper_figs(R, params, export_pdf)
%PLOT_PAPER_FIGS  Figures A-C from sweep_paper_figs.m, plus a summary table.
%
%   A  final distance vs wind: (a) W1 at the final time, (b) steady L1
%      (mean +- std over the last 40% of the run), with vs without low level
%   B  Lyapunov functional V(t) for every wind, with vs without low level
%   C  final agent density rho_N with and without the regularization
%      (same wind, same color scale in both panels)

if nargin < 3, export_pdf = false; end
U     = R.winds(:);
nU    = numel(U);
c_mpc = [0 0.45 0.74];
c_vel = [0.85 0.33 0.1];
pick  = @(C, f) cellfun(@(r) r.(f), C(:));
tailsd = @(C) cellfun(@(r) std(r.L1(floor(0.6 * numel(r.L1)) + 1:end)), C(:));
names = {'with low-level MPC', 'without low level'};
sets  = {R.mpc, R.vel};
figs  = gobjects(0);

% ---------------- A: final distances vs wind -------------------------------
figs(end + 1) = figure('Name', 'A: final distance vs wind', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 17 6.5]);
subplot(1, 2, 1); hold on; box on;
plot(U, pick(R.mpc, 'W1_final'), '-o', 'Color', c_mpc, 'LineWidth', 1.4, ...
    'MarkerFaceColor', c_mpc, 'DisplayName', names{1});
plot(U, pick(R.vel, 'W1_final'), '-s', 'Color', c_vel, 'LineWidth', 1.4, ...
    'MarkerFaceColor', c_vel, 'DisplayName', names{2});
xlabel('wind amplitude $U$');
ylabel('$W_1(\rho_T,\rho^*_T)$');
legend('Location', 'northwest');
title('(a) final Wasserstein distance');
subplot(1, 2, 2); hold on; box on;
errorbar(U, pick(R.mpc, 'L1_ss'), tailsd(R.mpc), '-o', 'Color', c_mpc, ...
    'LineWidth', 1.4, 'MarkerFaceColor', c_mpc, 'DisplayName', names{1});
errorbar(U, pick(R.vel, 'L1_ss'), tailsd(R.vel), '-s', 'Color', c_vel, ...
    'LineWidth', 1.4, 'MarkerFaceColor', c_vel, 'DisplayName', names{2});
xlabel('wind amplitude $U$');
ylabel('steady $\|\rho_N-\rho^*\|_{L^1}$');
legend('Location', 'northwest');
title('(b) steady density error (mean $\pm$ std)');
save_pdf(figs(end), 'fig_final_distance', export_pdf);

% ---------------- B: Lyapunov functional over time -------------------------
figs(end + 1) = figure('Name', 'B: Lyapunov functional', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 10 17 6.5]);
cols = parula(nU + 1);
ax = gobjects(2, 1);
for j = 1:2
    ax(j) = subplot(1, 2, j); hold on; box on;
    for k = 1:nU
        res = sets{j}{k};
        plot(res.t, res.V, 'Color', cols(k, :), 'LineWidth', 1.2, ...
            'DisplayName', sprintf('$U = %g$', U(k)));
    end
    set(gca, 'YScale', 'log');
    xlabel('time');
    ylabel('$V(\rho_t)$');
    title(sprintf('(%s) %s', char('a' + j - 1), names{j}));
end
legend(ax(1), 'Location', 'northeast');
linkaxes(ax, 'y');
save_pdf(figs(end), 'fig_lyapunov', export_pdf);

% ---------------- C: with / without regularization -------------------------
figs(end + 1) = figure('Name', 'C: with / without regularization', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [20 2 17 8]);
cases  = {R.reg, R.noreg};
labels = {sprintf('$\\varepsilon = %g$', R.reg.settings.epsilon_reg), ...
          sprintf('$\\varepsilon = %g$ (no regularization)', R.noreg.settings.epsilon_reg)};
% agent density rho_N (Wendland KDE at h_agent), same color scale in both panels
kdes = cell(1, 2);
for j = 1:2
    P = cases{j}.snap.pos(:, :, end);
    kdes{j} = sum(wendland_kernel(pdist2(params.grid_pts, P).^2, params.h_agent), 2) / size(P, 1);
end
cmax = max(cellfun(@max, kdes));
for j = 1:2
    res = cases{j};
    P = res.snap.pos(:, :, end);
    subplot(1, 2, j); cla; hold on;
    contourf(params.X, params.Y, reshape(kdes{j}, params.grid_res, params.grid_res), ...
        20, 'LineColor', 'none');
    colormap(gca, parula);
    caxis([0 cmax]);
    scatter(P(:, 1), P(:, 2), 2, 'w', 'MarkerEdgeAlpha', 0.2);
    axis equal;
    set(gca, 'XLim', [params.x(1) params.x(end)], 'YLim', [params.y(1) params.y(end)]);
    title(sprintf('%s, $W_1 = %.3f$', labels{j}, res.W1_final), 'Interpreter', 'latex');
    box on;
end
save_pdf(figs(end), 'fig_regularization', export_pdf);

% ---------------- summary table ---------------------------------------------
fprintf('\n%-6s %-20s %8s %8s %8s %9s %8s %8s %9s\n', 'U', 'controller', 'L1_ss', ...
    'L1_end', 'W1_end', 'W1_marg', 'offgrid', 'nn_drop', 'track');
for k = 1:nU
    for j = 1:2
        r = sets{j}{k};
        fprintf('%-6g %-20s %8.3f %8.3f %8.3f %9.1e %8.3f %8.3f %9.2e\n', U(k), names{j}, ...
            r.L1_ss, r.L1(end), r.W1_final, r.W1_marg_err, r.off_grid, r.nn_drop, r.track_ss);
    end
end
for j = 1:2
    r = cases{j};
    fprintf('regularization eps = %-5g: L1_ss = %.3f, W1_end = %.3f, nn_drop = %.3f, off-grid = %.3f\n', ...
        r.settings.epsilon_reg, r.L1_ss, r.W1_final, r.nn_drop, r.off_grid);
end
end

function save_pdf(fig, name, flag)
if flag
    exportgraphics(fig, [name '.pdf'], 'ContentType', 'vector');
end
end
