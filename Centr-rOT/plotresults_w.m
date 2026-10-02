function plotresults_w(pos_hist, phi_hist, target_hist, params)
%PLOTRESULTS_W  Save the run as the video OT_Wendland.avi.
%   Panels: target, agent density rho_N, potential phi, agent speed.
%   For DOOT step t: pos_hist(:,:,t) agent positions after the step,
%   phi_hist(:,t) phi on the grid, target_hist(:,:,t) target on the grid.

results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
fn = fullfile(results_dir, 'OT_Wendland.avi');
frame_step = 5;             % draw every 5th step
w = VideoWriter(fn); w.FrameRate = 10; open(w);
fig = figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
T = size(pos_hist, 3);
for t = unique([1:frame_step:T, T])
    draw_frame(fig, pos_hist, phi_hist, target_hist, params, t); drawnow limitrate;
    writeVideo(w, getframe(fig));
end
close(w); fprintf('Saved %s\n', fn);
end

function draw_frame(fig, pos_hist, phi_hist, target_hist, params, t)
figure(fig);
X = params.X; Y = params.Y; gr = params.grid_res;
pos = pos_hist(:, :, t);

% Panel 1: target
subplot(2,2,1); cla; hold on;
contourf(X, Y, double(target_hist(:, :, t)), 20, 'LineColor','none');
colormap(gca, parula); axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title('Target Distribution','FontSize',14,'Interpreter','latex'); box on;

% Panel 2: agent density rho_N (Wendland KDE at h_agent)
kde = sum(wendland_kernel(pdist2(params.grid_pts, pos).^2, params.h_agent), 2);
subplot(2,2,2); cla; hold on;
contourf(X, Y, reshape(kde, gr, gr), 20, 'LineColor','none');
colormap(gca, parula);
scatter(pos(:,1), pos(:,2), 2, 'w', 'MarkerEdgeAlpha',0.2);
axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title('Agent density $\rho_N$','FontSize',14,'Interpreter','latex'); box on;

% Panel 3: potential phi
subplot(2,2,3); cla; hold on;
contourf(X, Y, reshape(double(phi_hist(:, t)), gr, gr), 20, 'LineColor','none');
colorbar;
scatter(pos(:,1), pos(:,2), 5, 'w', 'filled', 'MarkerFaceAlpha',0.5);
axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title('Potential field','FontSize',14,'Interpreter','latex'); box on;

% Panel 4: agent speed over the last step, arrows along the motion
prev = pos_hist(:, :, max(t - 1, 1));
vel = (pos - prev) / params.dt_cont;
speed = sqrt(sum(vel.^2, 2));
subplot(2,2,4); cla; hold on;
scatter(pos(:,1), pos(:,2), 25, speed, 'filled', 'MarkerEdgeColor','k');
quiver(pos(:,1), pos(:,2), vel(:,1), vel(:,2), 1.5, 'w', 'LineWidth',1);
colormap(gca, turbo); colorbar;
caxis([0, max(0.1, max(speed))]);
axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title('Motion Magnitude ','FontSize',14,'Interpreter','latex'); box on;
end
