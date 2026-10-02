function output_file = plotresults_w(mode, plotData)
%PLOTRESULTS_W  WENDLAND version of plotresults2.
%   Panel 2 (density) reconstructs rho_N with a Wendland kernel at h_agent.
%   Panel 3 (potential) reconstructs phi with Wendland at h_field (+ h_ant).

if nargin < 1 || isempty(mode), mode = 'picture'; end
mode = lower(string(mode));
output_file = '';
switch mode
    case "video",   output_file = run_video(plotData);
    case "picture", output_file = run_picture(plotData);
    otherwise, error('Unknown mode %s', mode);
end
end

function fn = run_video(plotData)
fn = 'OT_Wendland.avi';
if isfield(plotData, 'video_speed') && isscalar(plotData.video_speed) && ...
        isfinite(plotData.video_speed) && plotData.video_speed > 0
    video_speed = plotData.video_speed;
else
    video_speed = 5;
end
frame_step = max(1, round(video_speed));
w = VideoWriter(fn); w.FrameRate = 10; open(w);
fig = figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
frame_idx = unique([1:frame_step:size(plotData.pos_hist,3), size(plotData.pos_hist,3)]);
for t = frame_idx
    draw_frame(fig, plotData, t); drawnow limitrate;
    writeVideo(w, getframe(fig));
end
close(w); fprintf('Saved %s\n', fn);
end

function fn = run_picture(plotData)
fn = 'OT_Wendland_final.png';
fig = figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
draw_frame(fig, plotData, size(plotData.pos_hist,3));
drawnow; print(fig, fn, '-dpng', '-r300'); fprintf('Saved %s\n', fn);
end

function draw_frame(fig, plotData, t)
figure(fig);
h_agent = plotData.h_agent;
if isfield(plotData,'h_field'), h_field = plotData.h_field; else, h_field = h_agent; end
if isfield(plotData,'h_ant'),   h_ant   = plotData.h_ant;   else, h_ant   = h_agent; end

gp = plotData.grid_pts; gr = plotData.grid_res;
current_pos = plotData.pos_hist(:, :, t);

% Panel 1: target
subplot(2,2,1); cla; hold on;
if isfield(plotData, 'target_hist')
    target_mat = plotData.target_hist(:, :, t);
else
    target_mat = plotData.rho_target_mat;
end
contourf(plotData.X, plotData.Y, target_mat, 20, 'LineColor','none');
colormap(gca, parula); axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title('Target Distribution','FontSize',14,'Interpreter','latex'); box on;

% Panel 2: DENSITY rho_N via Wendland at h_agent
sqd = pdist2(gp, current_pos).^2;
kde = sum(wendland_kernel(sqd, h_agent), 2);
subplot(2,2,2); cla; hold on;
contourf(plotData.X, plotData.Y, reshape(kde, gr, gr), 20, 'LineColor','none');
colormap(gca, parula);
scatter(current_pos(1:plotData.N,1), current_pos(1:plotData.N,2), 2, 'w', 'MarkerEdgeAlpha',0.2);
axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title(sprintf('Agent density $\\rho_N$'), ...
    'FontSize',14,'Interpreter','latex'); box on;

% Panel 3: POTENTIAL phi via Wendland field (+ antenna)
alpha = plotData.alpha_hist(:, t);
Kf = wendland_kernel(pdist2(gp, current_pos).^2, h_field);
if ~isempty(plotData.ant_cords)
    Ka = wendland_kernel(pdist2(gp, plotData.ant_cords).^2, h_ant);
    Kphi = [Kf, Ka];
else
    Kphi = Kf;
end
phi = Kphi * alpha;
subplot(2,2,3); cla; hold on;
contourf(plotData.X, plotData.Y, reshape(phi, gr, gr), 20, 'LineColor','none');
colorbar;
scatter(current_pos(1:plotData.N,1), current_pos(1:plotData.N,2), 5, 'w', 'filled', 'MarkerFaceAlpha',0.5);
if ~isempty(plotData.ant_cords)
    scatter(plotData.ant_cords(:,1), plotData.ant_cords(:,2), 20, 'y', 'filled', 'MarkerEdgeColor','k');
end
axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title('Potential field', ...
    'FontSize',14,'Interpreter','latex'); box on;

% Panel 4: motion magnitude
subplot(2,2,4); cla; hold on;
if t == 1, prev = plotData.pos_hist(:,:,1); else, prev = plotData.pos_hist(:,:,t-1); end
vel = (current_pos - prev) / plotData.dt_cont;
gmag = sqrt(sum(vel.^2, 2));
scatter(current_pos(1:plotData.N,1), current_pos(1:plotData.N,2), 25, gmag(1:plotData.N), 'filled', 'MarkerEdgeColor','k');
quiver(current_pos(:,1), current_pos(:,2), -vel(:,1), -vel(:,2), 1.5, 'w', 'LineWidth',1);
colormap(gca, turbo); cbar = colorbar;
caxis([0, max(0.1, max(gmag))]);
axis equal; set(gca,'XLim',[-8 8],'YLim',[-8 8]);
title('Motion Magnitude ', ...
    'FontSize',14,'Interpreter','latex'); box on;
end
