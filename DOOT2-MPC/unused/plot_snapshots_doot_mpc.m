function fig = plot_snapshots_doot_mpc(res, params, export_pdf)
%PLOT_SNAPSHOTS_DOOT_MPC  Agent density rho_N (Wendland KDE at h_agent) at
%   the stored steps, same style as panel 2 of plotresults_w and one common
%   color scale for all panels. res comes from run_doot_mpc.m.

if nargin < 3, export_pdf = false; end
K = numel(res.snap.step);
kdes = cell(1, K);
for k = 1:K
    P = res.snap.pos(:, :, k);
    kdes{k} = sum(wendland_kernel(pdist2(params.grid_pts, P).^2, params.h_agent), 2) / size(P, 1);
end
cmax = max(cellfun(@max, kdes));

fig = figure('Name', 'Snapshots', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.05 0.35 0.9 0.28]);
for k = 1:K
    P = res.snap.pos(:, :, k);
    subplot(1, K, k); cla; hold on;
    contourf(params.X, params.Y, reshape(kdes{k}, params.grid_res, params.grid_res), ...
        20, 'LineColor', 'none');
    colormap(gca, parula);
    caxis([0 cmax]);
    scatter(P(:, 1), P(:, 2), 2, 'w', 'MarkerEdgeAlpha', 0.2);
    axis equal;
    set(gca, 'XLim', [params.x(1) params.x(end)], 'YLim', [params.y(1) params.y(end)]);
    title(sprintf('$t = %g$', res.snap.step(k) * params.dt_cont), 'Interpreter', 'latex');
    box on;
end
if export_pdf
    exportgraphics(fig, 'snapshots.pdf', 'ContentType', 'vector');
end
end
