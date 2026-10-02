function plotresults_w(pos_hist, alpha_hist, rho_target, params, min_rho_terminal_hist)
%PLOTRESULTS_W  Save the run as the video OT_Wendland.avi.
%   Panels: target, agent density rho_N, potential phi, and either agent
%   speed or the active wind vector field.
%   pos_hist(:,:,t) and alpha_hist(:,t) are the positions and potential
%   coefficients after DOOT step t. rho_target is shown in panel 1.

if nargin < 5
    min_rho_terminal_hist = [];
end

results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
fn = fullfile(results_dir, 'OT_Wendland.avi');
frame_step = 5;             % draw every 5th step
w = VideoWriter(fn); w.FrameRate = 10; open(w);
fig = figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
T = size(pos_hist, 3);
for t = unique([1:frame_step:T, T])
    draw_frame(fig, pos_hist, alpha_hist, rho_target, params, t); drawnow limitrate;
    writeVideo(w, getframe(fig));
end
close(w); fprintf('Saved %s\n', fn);

if ~isempty(min_rho_terminal_hist)
    t_axis = (1:numel(min_rho_terminal_hist)) * params.dt_cont;
    negative = min_rho_terminal_hist < 0;
    figure('Color', 'w'); hold on; box on;
    plot(t_axis, min_rho_terminal_hist, 'k-', 'LineWidth', 1.5, ...
        'DisplayName', '$\min(\rho^* + \varepsilon\phi_N)$');
    yline(0, 'r--', 'LineWidth', 1.2, 'DisplayName', 'zero threshold');
    if any(negative)
        plot(t_axis(negative), min_rho_terminal_hist(negative), 'ro', ...
            'MarkerFaceColor', 'r', 'DisplayName', 'negative');
    end
    xlabel('time');
    ylabel('$\min(\rho^* + \varepsilon\phi_N)$');
    title('Terminal density positivity check');
    legend('Location', 'best');
    xlim([t_axis(1), t_axis(end)]);
    if any(negative)
        warning('main_w_final:negativeTerminalDensity', ...
            'rho^* + epsilon phi_N became negative at %d step(s).', nnz(negative));
    end
end
end

function draw_frame(fig, pos_hist, alpha_hist, rho_target, params, t)
figure(fig);
X = params.X; Y = params.Y; gr = params.grid_res;
pos = pos_hist(:, :, t);
sq_dists = pdist2(params.grid_pts, pos).^2;     % grid points x agents

% Panel 1: target
subplot(2,2,1); cla; hold on;
contourf(X, Y, reshape(rho_target, gr, gr), 20, 'LineColor','none');
colormap(gca, parula); axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title('Target Distribution','FontSize',14,'Interpreter','latex'); box on;

% Panel 2: agent density rho_N (Wendland KDE at h_agent)
kde = sum(wendland_kernel(sq_dists, params.h_agent), 2);
subplot(2,2,2); cla; hold on;
contourf(X, Y, reshape(kde, gr, gr), 20, 'LineColor','none');
colormap(gca, parula);
scatter(pos(:,1), pos(:,2), 2, 'w', 'MarkerEdgeAlpha',0.2);
axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title('Agent density $\rho_N$','FontSize',14,'Interpreter','latex'); box on;

% Panel 3: potential phi = sum_j alpha_j K_field(x - x_j)
phi = wendland_kernel(sq_dists, params.h_field) * alpha_hist(:, t);
subplot(2,2,3); cla; hold on;
contourf(X, Y, reshape(phi, gr, gr), 20, 'LineColor','none');
colorbar;
scatter(pos(:,1), pos(:,2), 5, 'w', 'filled', 'MarkerFaceAlpha',0.5);
axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title('Potential field','FontSize',14,'Interpreter','latex'); box on;

subplot(2,2,4); cla; hold on;
if isfield(params, 'wind_model') && ~strcmp(params.wind_model, 'none')
    wind_pts = [X(:), Y(:)];
    wind = wind_at_points(wind_pts, params);
    skip = max(1, ceil(gr / 20));
    idx = 1:skip:size(wind_pts, 1);
    quiver(wind_pts(idx,1), wind_pts(idx,2), wind(idx,1), wind(idx,2), ...
        'Color', [0.1 0.2 0.7], 'LineWidth', 1);
    scatter(pos(:,1), pos(:,2), 12, 'k', 'filled');
    title('Wind vector field','FontSize',14,'Interpreter','latex');
else
    prev = pos_hist(:, :, max(t - 1, 1));
    vel = (pos - prev) / params.dt_cont;
    speed = sqrt(sum(vel.^2, 2));
    scatter(pos(:,1), pos(:,2), 25, speed, 'filled', 'MarkerEdgeColor','k');
    quiver(pos(:,1), pos(:,2), -vel(:,1), -vel(:,2), 1.5, 'w', 'LineWidth',1);
    colormap(gca, turbo); colorbar;
    caxis([0, max(0.1, max(speed))]);
    title('Motion Magnitude ','FontSize',14,'Interpreter','latex');
end
axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
box on;
end

function d = wind_at_points(points, params)
switch params.wind_model
    case 'uniform'
        d = uniform_wind(0, points, params.wind_U, params.wind_theta);
    case 'taylor_green'
        d = taylor_green_wind(0, points, params.wind_U, params.wind_L);
    otherwise
        error('plotresults_w:windModel', ...
            'Unsupported active wind model ''%s''.', params.wind_model);
end
end
