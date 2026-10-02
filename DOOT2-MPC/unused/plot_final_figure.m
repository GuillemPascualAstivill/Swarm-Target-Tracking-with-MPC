function figs = plot_final_figure(R, Rrot, params, varargin)
%PLOT_FINAL_FIGURE  Paper figure with four EQUAL square panels.
%   (a) 'map'  : regime map over (eps, target speed) from Rmap: colour =
%                inter-agent spacing (clustering), black lines = control
%                effort in decades (needs 'Rmap', sweep_regime_map.m)
%       'lines': metric vs target speed, with / without regularization (Rrot)
%   (b) 'time' : L1 error over time for every wind amplitude, with (solid)
%                and without (dashed) low level (R)
%       'final': metric vs wind amplitude, with / without low level (R)
%   (c) final agent density with regularization     (R.reg)
%   (d) final agent density without regularization  (R.noreg)
%
%   plot_final_figure(R, Rrot, params, Name, Value, ...)
%     'panel_a'   'map' (default if Rmap is given) | 'time' | 'lines'
%                 'time': like panel (b), one curve per eps over time at one
%                 target speed, from the Rmap runs (no new simulations)
%     'a_var'     'effort' (default) | 'spacing'   variable of panel (a) 'time'
%     'a_revs'    1       target speed of panel (a) 'time' (rev per 500 steps)
%     'a_eps'     []      eps values shown (default: 0 0.05 0.25)
%     'b_winds'   [0 0.03 0.1]  wind amplitudes shown in panel (b) 'time'
%     'a_ylabel_x' -0.10  horizontal position of the y label of panel (a)
%                 'time' (fraction of the axes width; closer to 0 = closer
%                 to the axis; [] = MATLAB's automatic position)
%     'a_legend'  'east'  legend location of panel (a) 'time' (any legend
%                 Location, e.g. 'northeast', 'southwest')
%     'panel_b'   'time' (default) | 'final'
%     'Rmap'      struct from sweep_regime_map.m (regime_map_results.mat);
%                 seeds (3rd dimension of Rmap.res) are averaged
%     'Rwind'     struct from sweep_wind_seeds.m (wind_seeds_results.mat):
%                 panel (b) then shows the seed mean with a min-max band;
%                 without it, panel (b) uses the single runs in R
%     'eps_star'  0.25   eps marked in the map (value used in the paper)
%     'metric_a'  'effort' (default) | 'close' | 'nn' | 'W1' | 'L1'
%     'metric_b'  'W1' (default) | 'L1'
%     'layout'    'grid'     : 2x2 figure, all four axes square and equal,
%                              exported as fig_final.pdf (grid_cm wide)
%                 'separate' : one figure per panel, identical page size and
%                              axes box, exported as fig_a.pdf ... fig_d.pdf
%                              (no titles: use \subfloat captions in LaTeX)
%     'export'    false (default) | true
%     'file'      'paper-fig'  name of the 2x2 PDF ('grid' layout); it is
%                 printed at exactly grid_cm x grid_cm with font_pt fonts, so
%                 \includegraphics[width=\columnwidth]{paper-fig} is scale 1
%     'panel_cm'  [4.3 4.3]  page size of each separate panel: fits
%                 0.24\textwidth in a figure* row, or 0.49\columnwidth in a
%                 2x2 single-column figure, at scale ~1 (fonts stay font_pt)
%     'font_pt'   8
%     'grid_cm'   13         width and height of the square 2x2 figure
%                 \includegraphics[width=\columnwidth] keeps fonts at font_pt
%   Density panels share one color scale capped at the peak of the target
%   smoothed with the same kernel; denser clumps saturate.

ip = inputParser;
ip.addParameter('metric_a', 'effort');
ip.addParameter('metric_b', 'W1');
ip.addParameter('layout', 'grid');
ip.addParameter('export', false);
ip.addParameter('panel_cm', [4.3 4.3]);
ip.addParameter('font_pt', 10);
ip.addParameter('grid_cm', 13);
ip.addParameter('file', 'paper-fig');
ip.addParameter('panel_a', '');
ip.addParameter('panel_b', 'time');
ip.addParameter('Rmap', []);
ip.addParameter('Rwind', []);
ip.addParameter('eps_star', 0.25);
ip.addParameter('a_var', 'effort');
ip.addParameter('a_revs', 1);
ip.addParameter('a_eps', []);
ip.addParameter('a_legend', 'east');
ip.addParameter('b_winds', [0.03 0.05 0.1]);
ip.addParameter('a_ylabel_x', -0.10);
ip.parse(varargin{:});
o = ip.Results;
o.h_eval = 1;
o.compact = strcmp(o.layout, 'separate') || o.grid_cm < 12;
if isempty(o.panel_a)
    if isempty(o.Rmap), o.panel_a = 'lines'; else, o.panel_a = 'map'; end
end

% common color scale for the density panels
tgt  = star_target_w(params, R.reg.snap.phase(end));
cmax = max(smooth_on_grid_wendland(tgt, params, params.h_agent));

switch o.panel_a
    case 'map',   fa = @(ax) draw_map(ax, o.Rmap, params, o);  ta = '(a)';
    case 'time',  fa = @(ax) draw_a_time(ax, o.Rmap, params, o); ta = '(a)';
    case 'lines', fa = @(ax) draw_a(ax, Rrot, params, o);      ta = '(a)';
    otherwise,    error('panel_a must be ''map'', ''time'' or ''lines''.');
end
switch o.panel_b
    case 'time',  fb = @(ax) draw_b_time(ax, R, o);            tb = '(b)';
        if ~isempty(o.Rwind), fb = @(ax) draw_b_time(ax, o.Rwind, o); end
    case 'final', fb = @(ax) draw_b(ax, R, params, o);         tb = '(b)';
    otherwise,    error('panel_b must be ''time'' or ''final''.');
end
draw = {fa, fb, ...
        @(ax) draw_snap(ax, R.reg, params, cmax, o), ...
        @(ax) draw_snap(ax, R.noreg, params, cmax, o)};
titles = {ta, tb, ...
          sprintf('(c) $\\varepsilon = %g$', R.reg.settings.epsilon_reg), ...
          sprintf('(d) $\\varepsilon = %g$', R.noreg.settings.epsilon_reg)};
figs = gobjects(0);

switch o.layout
    case 'grid'
        fig = figure('Name', 'final figure', 'Color', 'w', ...
            'Units', 'centimeters', 'Position', [2 2 o.grid_cm o.grid_cm]);
        tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        for k = 1:4
            ax = nexttile(tl);
            draw{k}(ax);
            axis(ax, 'square');                 % identical square boxes
            set(ax, 'FontSize', o.font_pt, 'TickLabelInterpreter', 'latex');
            if k == 1 && strcmp(o.metric_a, 'effort')
                yl = ylim(ax);
                set(ax, 'YTick', yl);
            end
            title_size = o.font_pt + double(k >= 3);
            title(ax, titles{k}, 'Interpreter', 'latex', 'FontSize', title_size);
        end
        figs = fig;
        if o.export
            save_exact(fig, o.file, o.grid_cm, o.grid_cm);
            fprintf('Saved %s.pdf (%.1f x %.1f cm, %g pt fonts).\n', o.file, ...
                o.grid_cm, o.grid_cm, o.font_pt);
        end

    case 'separate'
        W = o.panel_cm(1);  H = o.panel_cm(2);
        box_pos = [1.15, 0.95, W - 1.35, H - 1.25];   % same axes box in every panel
        names = {'fig_a', 'fig_b', 'fig_c', 'fig_d'};
        for k = 1:4
            fig = figure('Name', names{k}, 'Color', 'w', 'Units', 'centimeters', ...
                'Position', [2 + 5 * (k - 1), 4, W, H]);
            ax = axes(fig, 'Units', 'centimeters', 'Position', box_pos);
            draw{k}(ax);
            set(ax, 'FontSize', o.font_pt, 'TickLabelInterpreter', 'latex');
            figs(end + 1) = fig; %#ok<AGROW>
            if o.export
                save_exact(fig, names{k}, W, H);
            end
        end
        fprintf('Subcaption values: (c) eps = %g, W1 = %.3f | (d) eps = %g, W1 = %.3f\n', ...
            R.reg.settings.epsilon_reg, R.reg.W1_final, ...
            R.noreg.settings.epsilon_reg, R.noreg.W1_final);
    otherwise
        error('layout must be ''grid'' or ''separate''.');
end
end

% =============================================================================
function draw_a(ax, Rrot, params, o)
hold(ax, 'on'); box(ax, 'on');
[c1, c2, lw, ms] = style(o);
spd = cellfun(@(r) r.target_speed, Rrot.reg(:));
f   = @(res) panel_value(res, params, o.metric_a, o.h_eval);
ya  = cellfun(f, Rrot.reg(:));
yn  = cellfun(f, Rrot.noreg(:));
plot(ax, spd, ya, '-o', 'Color', c1, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', c1, 'DisplayName', sprintf('$\\varepsilon = %g$', Rrot.eps(1)));
plot(ax, spd, yn, '-s', 'Color', c2, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', c2, 'DisplayName', sprintf('$\\varepsilon = %g$', Rrot.eps(2)));
if strcmp(o.metric_a, 'nn')
    hy = yline(ax, reference_nn(params, Rrot.reg{1}, 20), 'k--');
    hy.HandleVisibility = 'off';
end
if any(strcmp(o.metric_a, {'nn', 'effort'}))
    set(ax, 'YScale', 'log');
end
xlim(ax, [0, max(spd) * 1.05]);
xlabel(ax, pick_label(o, 'target speed $\omega r_0$', '$\omega r_0$'), 'Interpreter', 'latex');
hy = ylabel(ax, axis_label(o.metric_a, o), 'Interpreter', 'latex');
if strcmp(o.metric_a, 'effort') && ~isempty(o.a_ylabel_x)
    hy.Units = 'normalized';
    hy.Position(1) = o.a_ylabel_x;
    hy.Clipping = 'off';
end
lg = legend(ax, 'Location', 'best', 'Interpreter', 'latex');
if o.compact, lg.FontSize = o.font_pt - 1; lg.ItemTokenSize = [12 8]; end
fprintf('(a) %s vs target speed:\n', o.metric_a);
fprintf('   %8.3f  %10.4g  %10.4g\n', [spd, ya, yn]');
end

function draw_b(ax, R, params, o)
hold(ax, 'on'); box(ax, 'on');
[c1, c2, lw, ms] = style(o);
U  = R.winds(:);
f  = @(res) panel_value(res, params, o.metric_b, o.h_eval);
ym = cellfun(f, R.mpc(:));
yv = cellfun(f, R.vel(:));
plot(ax, U, ym, '-o', 'Color', c1, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', c1, 'DisplayName', pick_label(o, 'with low-level MPC', 'MPC'));
plot(ax, U, yv, '-s', 'Color', c2, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', c2, 'DisplayName', pick_label(o, 'without low level', 'no low level'));
xlim(ax, [0, max(U) * 1.05]);
xlabel(ax, pick_label(o, 'wind amplitude $U$', '$U$'), 'Interpreter', 'latex');
hy = ylabel(ax, axis_label(o.metric_b, o), 'Interpreter', 'latex');
if strcmp(o.metric_b, 'effort') && ~isempty(o.a_ylabel_x)
    hy.Units = 'normalized';
    hy.Position(1) = o.a_ylabel_x;
end
lg = legend(ax, 'Location', 'northwest', 'Interpreter', 'latex');
if o.compact, lg.FontSize = o.font_pt - 1; lg.ItemTokenSize = [12 8]; end
fprintf('(b) %s vs wind:\n', o.metric_b);
fprintf('   %8.3f  %10.4g  %10.4g\n', [U, ym, yv]');
end

function draw_map(ax, Rmap, params, o)
% Regime map over (eps, target speed), seed-averaged:
%   colour : inter-agent spacing = median nearest-neighbour distance divided
%            by that of independent samples of rho* (1 = as regular as ideal
%            samples, 0 = collapsed into clusters); labelled contours
%   black  : control effort per agent, one line per decade (10^k)
%   dots   : simulated settings; dotted line: eps used in the paper
hold(ax, 'on'); box(ax, 'on');
E = Rmap.eps(:)';
S = cellfun(@(r) r.target_speed, Rmap.res(:, 1, 1))';
ref = reference_nn(params, Rmap.res{1, 1, 1}, 20);
Sp  = mean(cellfun(@(r) panel_value(r, params, 'nn', o.h_eval), Rmap.res), 3) / ref;
Sp  = min(Sp, 1);
Ef  = mean(cellfun(@(r) panel_value(r, params, 'effort', o.h_eval), Rmap.res), 3);
[EE, SS] = meshgrid(E, S);
fs = max(o.font_pt - 1, 6);

% colour: spacing
lev = linspace(0, 1, 11);
[cm, hc] = contourf(ax, EE, SS, Sp, lev, 'LineColor', [1 1 1] * 0.55, 'LineWidth', 0.3);
colormap(ax, parula);
caxis(ax, [0 1]);
clabel(cm, hc, lev([3 6 9]), 'FontSize', fs, 'Color', 'w', 'LabelSpacing', 500);

% black lines: control effort in decades, labelled 10^k at the longest piece
Lg = log10(Ef);
dec = ceil(min(Lg(:))):floor(max(Lg(:)));
if ~isempty(dec)
    if isscalar(dec), dec = [dec dec]; end
    Cc = contour(ax, EE, SS, Lg, dec, 'k-', 'LineWidth', 0.9);
    best = containers.Map('KeyType', 'double', 'ValueType', 'any');
    k = 1;
    while k < size(Cc, 2)
        lv = Cc(1, k);  np = Cc(2, k);
        seg = Cc(:, k + 1:k + np);
        if ~isKey(best, lv) || size(best(lv), 2) < np
            best(lv) = seg;
        end
        k = k + np + 1;
    end
    for lv = cell2mat(keys(best))
        seg = best(lv);
        m = seg(:, max(1, round(size(seg, 2) / 2)));
        text(ax, m(1), m(2), sprintf('$10^{%d}$', round(lv)), 'Interpreter', 'latex', ...
            'FontSize', fs, 'BackgroundColor', 'w', 'Margin', 0.5, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

plot(ax, EE(:), SS(:), 'k.', 'MarkerSize', 3);
hx = xline(ax, o.eps_star, 'k:', 'LineWidth', 1.2);
hx.HandleVisibility = 'off';
xlim(ax, [E(1) E(end)]);  ylim(ax, [S(1) S(end)]);
xlabel(ax, '$\varepsilon$', 'Interpreter', 'latex');
ylabel(ax, pick_label(o, 'target speed $\omega r_0$', '$\omega r_0$'), 'Interpreter', 'latex');
fprintf('(a) rows: speed %s; cols: eps %s\n', mat2str(S, 3), mat2str(E, 3));
fprintf('    spacing (median NN / ideal, 1 = regular):\n');  disp(round(Sp, 2));
fprintf('    control effort per agent:\n');                  disp(Ef);
end

function draw_a_time(ax, Rmap, params, o)
% Same style as panel (b): one curve per regularization eps (colour), at one
% target speed, seed mean with a light min-max band, 10-step moving average.
%   o.a_var 'effort'  : control effort per agent and unit time (log scale)
%           'spacing' : median NN distance / that of samples of rho* (1 = regular)
hold(ax, 'on'); box(ax, 'on');
[~, ~, lw] = style(o);
i = find(abs(Rmap.revs - o.a_revs) < 1e-9, 1);
if isempty(i)
    error('a_revs = %g is not in Rmap.revs = %s', o.a_revs, mat2str(Rmap.revs));
end
if isempty(o.a_eps)
    want = [0 0.05 0.25];
    J = find(any(abs(Rmap.eps(:) - want) < 1e-12, 2))';
    if isempty(J), J = 1:numel(Rmap.eps); end
else
    J = arrayfun(@(e) find(abs(Rmap.eps - e) < 1e-12, 1), o.a_eps);
end
n = numel(J);
cols = parula(n + 1);
cols = cols(1:n, :);
r0 = Rmap.res{i, J(1), 1};
t  = r0.t;
N  = size(r0.pos_final, 1);
dt = t(1);
switch o.a_var
    case 'effort'
        f = @(r) movmean(r.effort(:) / (N * dt), 10);
        ylab = '$\|u\|^2$';
        logy = true;
    case 'spacing'
        ref = reference_nn(params, r0, 20);
        f = @(r) movmean(r.nn_med(:) / ref, 10);
        ylab = pick_label(o, 'inter-agent spacing (1 = regular)', 'spacing');
        logy = false;
    otherwise
        error('a_var must be ''effort'' or ''spacing''.');
end
Y = cell(n, 1);
for q = 1:n
    Y{q} = cell2mat(cellfun(f, reshape(Rmap.res(i, J(q), :), 1, []), 'UniformOutput', false));
end
if size(Y{1}, 2) > 1                          % bands first, lines on top
    for q = 1:n
        fill(ax, [t; flipud(t)], [min(Y{q}, [], 2); flipud(max(Y{q}, [], 2))], ...
            0.3 * cols(q, :) + 0.7, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
h = gobjects(n, 1);
for q = 1:n
    e  = Rmap.eps(J(q));
    w  = lw * (1 + (abs(e - o.eps_star) < 1e-12));      % eps used in the paper: thicker
    h(q) = plot(ax, t, mean(Y{q}, 2), '-', 'Color', cols(q, :), 'LineWidth', w, ...
        'DisplayName', sprintf('$\\varepsilon=%g$', e));
end
lg = legend(ax, h, 'Location', o.a_legend, 'Interpreter', 'latex');
lg.FontSize = max(o.font_pt - 2, 6);
lg.ItemTokenSize = [8 6];
if logy, set(ax, 'YScale', 'log'); end
xlim(ax, [0, t(end)]);
xlabel(ax, '', 'Interpreter', 'latex');
hy = ylabel(ax, ylab, 'Interpreter', 'latex');
if logy && ~isempty(o.a_ylabel_x)       % MATLAB leaves a wide gap on log axes
    hy.Units = 'normalized';
    hy.Position(1) = o.a_ylabel_x;
    hy.Clipping = 'off';
end
fprintf('(a) %s over time at %.2f rev/500 steps, eps = %s\n', o.a_var, Rmap.revs(i), ...
    mat2str(Rmap.eps(J)));
end

function draw_b_time(ax, R, o)
% L1 error over time for a few wind amplitudes, with the low-level MPC and
% without low level; 10-step moving average.
%   hue       = controller: blue solid = MPC, orange dashed = no low level
%   lightness = wind amplitude: light = weak, dark = strong
% R.mpc / R.vel are nU x nSeeds cells (Rwind) or nU x 1 (single runs):
% lines are the seed mean, light bands the seed min-max.
hold(ax, 'on'); box(ax, 'on');
[c_mpc, c_vel, lw] = style(o);
lw = 1.2 * lw;
K = find(any(abs(R.winds(:) - o.b_winds(:)') < 1e-12, 2))';
if isempty(K), K = 1:numel(R.winds); end
U  = R.winds(K);
nU = numel(U);
a  = linspace(0.4, 1, nU);                          % strength of the colour per wind
tint = @(c, ak) ak * c + (1 - ak) * [1 1 1];
t = R.mpc{1}.t;
curves = @(C, k) cell2mat(cellfun(@(r) movmean(r.L1(:), 10), C(k, :), 'UniformOutput', false));
Ym = cell(nU, 1);  Yv = cell(nU, 1);
for k = 1:nU
    Ym{k} = curves(R.mpc, K(k));
    Yv{k} = curves(R.vel, K(k));
end
if size(Ym{1}, 2) > 1                               % bands first, lines on top
    for k = 1:nU
        band(ax, t, Yv{k}, tint(c_vel, 0.25));
        band(ax, t, Ym{k}, tint(c_mpc, 0.25));
    end
end
for k = 1:nU                                        % no low level below, MPC on top
    plot(ax, t, mean(Yv{k}, 2), '--', 'Color', tint(c_vel, a(k)), 'LineWidth', lw, ...
        'HandleVisibility', 'off');
end
for k = 1:nU
    plot(ax, t, mean(Ym{k}, 2), '-', 'Color', tint(c_mpc, a(k)), 'LineWidth', lw, ...
        'HandleVisibility', 'off');
end
% legend: controllers by hue and line style, winds by lightness (grey)
hs = [plot(ax, nan, nan, '-',  'Color', c_mpc, 'LineWidth', lw, ...
          'DisplayName', pick_label(o, 'with low-level MPC', 'MPC')); ...
      plot(ax, nan, nan, '--', 'Color', c_vel, 'LineWidth', lw, ...
          'DisplayName', pick_label(o, 'without low level', 'no low level'))];
hU = gobjects(nU, 1);
for k = 1:nU
    hU(k) = plot(ax, nan, nan, '-', 'Color', tint([0 0 0], a(k)), 'LineWidth', lw, ...
        'DisplayName', sprintf('$U=%g$', U(k)));
end
xlim(ax, [0, t(end)]);
xlabel(ax, '', 'Interpreter', 'latex');
ylabel(ax, '$\|\rho_N-\rho^*\|_{L^1}$', 'Interpreter', 'latex');
lg = legend(ax, [hs; hU], 'Location', 'northeast', 'Interpreter', 'latex');
lg.FontSize = max(o.font_pt - 2, 6);
lg.ItemTokenSize = [10 6];
end

function band(ax, t, Y, col)
fill(ax, [t; flipud(t)], [min(Y, [], 2); flipud(max(Y, [], 2))], col, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function draw_snap(ax, res, params, cmax, o)
hold(ax, 'on');
P = res.snap.pos(:, :, end);
kde = sum(wendland_kernel(pdist2(params.grid_pts, P).^2, params.h_agent), 2) / size(P, 1);
% smooth image (not contourf): vector PDFs show seams between filled contour
% bands; an image is embedded as one bitmap, so there are no lines
Z  = reshape(min(kde, cmax), params.grid_res, params.grid_res);
xf = linspace(params.x(1), params.x(end), 300);
yf = linspace(params.y(1), params.y(end), 300);
[XF, YF] = meshgrid(xf, yf);
Zf = min(max(interp2(params.X, params.Y, Z, XF, YF, 'cubic'), 0), cmax);
imagesc(ax, xf, yf, Zf);
set(ax, 'YDir', 'normal');
colormap(ax, parula);
caxis(ax, [0 cmax]);
scatter(ax, P(:, 1), P(:, 2), 1, [0.9 0.9 0.9], 'filled');
axis(ax, 'equal');
set(ax, 'XLim', [params.x(1) params.x(end)], 'YLim', [params.y(1) params.y(end)]);
if o.compact
    set(ax, 'XTick', [-5 0 5], 'YTick', [-5 0 5]);
end
box(ax, 'on');
end

% =============================================================================
function [c1, c2, lw, ms] = style(o)
c1 = [0 0.45 0.74];
c2 = [0.85 0.33 0.1];
if o.compact, lw = 1.0; ms = 3.5; else, lw = 1.4; ms = 6; end
end

function s = pick_label(o, long, short)
if o.compact, s = short; else, s = long; end
end

function s = axis_label(metric, o)
switch metric
    case 'effort', s = '$\|u\|^2$';
    case 'close',  s = pick_label(o, 'agents with a neighbour closer than 0.05 (\%)', 'near-collisions (\%)');
    case 'nn',     s = pick_label(o, 'median nearest-neighbour distance', 'median NN distance');
    case 'W1',     s = pick_label(o, 'final $W_1(\rho_T,\rho^*_T)$', '$W_1(\rho_T,\rho^*_T)$');
    case 'L1',     s = sprintf('$\\|\\rho_N-\\rho^*\\|_{L^1}$ ($h = %g$)', o.h_eval);
end
end

function y = panel_value(res, params, metric, h_eval)
switch metric
    case 'effort'
        y = res.effort_total / (size(res.pos_final, 1) * res.t(end));
    case 'close'
        if isfield(res, 'close_frac_ss')
            y = 100 * res.close_frac_ss;
        else
            y = 100 * mean(nn_distances(res.pos_final) < 0.05);
        end
    case 'nn'
        if isfield(res, 'nn_med_ss')
            y = res.nn_med_ss;
        else
            y = median(nn_distances(res.pos_final));
        end
    case 'W1'
        y = res.W1_final;
    case 'L1'
        P   = res.pos_final;
        rho = sum(wendland_kernel(pdist2(params.grid_pts, P).^2, h_eval), 2);
        rho = rho / (sum(rho) * params.dx * params.dy);
        tgt = star_target_w(params, res.snap.phase(end));
        y   = sum(abs(rho - tgt)) * params.dx * params.dy;
    otherwise
        error('metric must be ''effort'', ''close'', ''nn'', ''W1'' or ''L1''.');
end
end

function d = nn_distances(P)
N = size(P, 1);
d2 = (P(:, 1) - P(:, 1)').^2 + (P(:, 2) - P(:, 2)').^2;
d2(1:N+1:end) = Inf;
d = sqrt(min(d2, [], 2));
end

function m = reference_nn(params, res, ndraw)
% Median NN distance of N independent samples of the target.
N = size(res.pos_final, 1);
w = star_target_w(params, res.snap.phase(end));
c = cumsum(w) / sum(w);
s = rng;  rng(1);
m = 0;
for k = 1:ndraw
    idx = discretize(rand(N, 1), [0; c]);
    P = params.grid_pts(idx, :) + (rand(N, 2) - 0.5) * params.dx;
    m = m + median(nn_distances(P)) / ndraw;
end
rng(s);
end

function save_exact(fig, name, W, H)
% Vector PDF of exactly W x H cm (fonts keep their point size), independent
% of the screen resolution and of Windows display scaling.
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [W H], ...
    'PaperPositionMode', 'manual', 'PaperPosition', [0 0 W H]);
try
    print(fig, name, '-dpdf', '-vector');       % R2022a and later
catch
    print(fig, name, '-dpdf', '-painters');     % older releases
end
end
