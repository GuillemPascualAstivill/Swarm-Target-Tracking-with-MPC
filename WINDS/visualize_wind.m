function visualize_wind(wind_fn, name, dom, t_eval, n_grid, v_ref)
% VISUALIZE_WIND  Plot a drift field: direction (unit quiver) over magnitude.
%
%   VISUALIZE_WIND(wind_fn) plots the field over the working domain of
%   doot_mpc.m. wind_fn is the same handle the simulation uses,
%   d = wind_fn(t, x), so what is drawn is exactly what the agents feel --
%   no second definition to drift out of sync.
%
%   Arrows are NORMALIZED to unit length, so direction and magnitude are read
%   from separate channels: every arrow is the same size and shows only where
%   the flow points, while the background colour carries |d|. A raw quiver
%   conflates the two and goes unreadable near stagnation points, where the
%   direction is still well defined but the arrow vanishes.
%
%   Inputs
%     wind_fn : handle d = wind_fn(t, x), x an Mx2 stack (both wind files
%               accept this)
%     name    : title string (optional; default 'wind field')
%     dom     : half-width of the square domain (optional; default 9, the
%               plotting range doot_mpc.m uses)
%     t_eval  : time at which to evaluate (optional; default 0). Only matters
%               for a time-modulated field.
%     n_grid  : arrows per axis (optional; default 25)
%     v_ref   : reference speed for the strength annotation (optional; default
%               0.3, the typical |v_nom| in doot_mpc.m). The colourbar is
%               labelled with |d|/v_ref so the field can be judged against
%               what the swarm is actually doing -- |d|/v_ref near 1 means the
%               wind is as strong as the transport command.
%
%   Examples
%     visualize_wind(@(t,x) uniform_wind(t,x,0.10,0), 'uniform')
%     visualize_wind(@(t,x) taylor_green_wind(t,x,0.10,6), 'Taylor-Green')
%
%   Reads the max spatial gradient off the sampled grid and prints it, since
%   ||grad d|| is what survives the -d feedforward in low_level_mpc.m: the
%   standing tracking error goes like ||grad d||*|v_nom| / sqrt(Q_x/Q_u).
%   A uniform field reports ~0 there, which is the point of running it first.

    % No arguments (e.g. the editor's Run button): show both fields the
    % simulation ships with, so the file is directly runnable.
    if nargin < 1 || isempty(wind_fn)
        visualize_wind(@(t,x) uniform_wind(t,x,0.10,0), 'uniform');
        visualize_wind(@(t,x) taylor_green_wind(t,x,0.10,6), 'Taylor-Green');
        return;
    end

    if nargin < 2 || isempty(name),   name   = 'wind field'; end
    if nargin < 3 || isempty(dom),    dom    = 9;    end
    if nargin < 4 || isempty(t_eval), t_eval = 0;    end
    if nargin < 5 || isempty(n_grid), n_grid = 25;   end
    if nargin < 6 || isempty(v_ref),  v_ref  = 0.3;  end

    % ---- Coarse grid for the arrows -------------------------------------
    q = linspace(-dom, dom, n_grid);
    [Xq, Yq] = meshgrid(q, q);
    Dq = wind_fn(t_eval, [Xq(:), Yq(:)]);
    Dqx = reshape(Dq(:,1), size(Xq));
    Dqy = reshape(Dq(:,2), size(Xq));
    Sq  = hypot(Dqx, Dqy);

    % Unit arrows. Stagnation points (|d| -> 0) have no defined direction, so
    % they are dropped rather than drawn as noise.
    tol = 1e-9 * max(Sq(:), [], 'omitnan');
    Uq = Dqx ./ max(Sq, eps);
    Vq = Dqy ./ max(Sq, eps);
    dead = Sq <= max(tol, eps);
    Uq(dead) = NaN;  Vq(dead) = NaN;

    % ---- Fine grid for the magnitude background -------------------------
    p = linspace(-dom, dom, 8*n_grid);
    [Xp, Yp] = meshgrid(p, p);
    Dp = wind_fn(t_eval, [Xp(:), Yp(:)]);
    Sp = reshape(hypot(Dp(:,1), Dp(:,2)), size(Xp));

    Smax = max(Sp(:));

    figure('Units','pixels','Position',[100 100 660 580]);

    imagesc(p, p, Sp / max(Smax, eps));
    set(gca, 'YDir', 'normal');
    colormap(parula);
    caxis([0 1]);
    cb = colorbar;
    cb.Label.String = sprintf('$\\|d\\| / \\|d\\|_{\\max}$ (max $= %.3g$)', Smax);
    cb.Label.Interpreter = 'latex';
    hold on;

    % Scale so arrows span most of a grid cell without overlapping.
    cell_w = 2*dom / (n_grid - 1);
    quiver(Xq, Yq, Uq*cell_w*0.45, Vq*cell_w*0.45, 0, 'k', ...
           'LineWidth', 0.8, 'MaxHeadSize', 0.6);

    axis equal tight; grid on;
    xlabel('$x$'); ylabel('$y$');
    title(sprintf('%s at $t = %.3g$ (arrows unit-length)', ...
                  strrep(name, '_', '\_'), t_eval), 'Interpreter', 'latex');
    hold off;

    % ---- Numbers worth having next to the picture -----------------------
    % ||grad d|| from central differences on the fine grid. This is the
    % quantity the feedforward cannot cancel.
    h = p(2) - p(1);
    [dxdx, dxdy] = gradient(reshape(Dp(:,1), size(Xp)), h, h);
    [dydx, dydy] = gradient(reshape(Dp(:,2), size(Xp)), h, h);
    % Frobenius norm bounds the spectral norm and needs no per-point eig.
    gnorm = sqrt(dxdx.^2 + dxdy.^2 + dydx.^2 + dydy.^2);
    % Trim the border, where central differences are one-sided.
    gnorm = gnorm(2:end-1, 2:end-1);
    div   = dxdx + dydy;
    div   = div(2:end-1, 2:end-1);

    fprintf('\n---- %s (t = %.3g) ----\n', name, t_eval);
    fprintf('  max |d|                 %.4f\n', Smax);
    fprintf('  max |d| / v_ref         %.3f   (v_ref = %.3g)\n', ...
            Smax / max(v_ref, eps), v_ref);
    fprintf('  max |grad d| (Frob)     %.4f\n', max(gnorm(:)));
    fprintf('  max |div d|             %.3e   (0 => incompressible)\n', ...
            max(abs(div(:))));
    fprintf(['  Standing tracking error under this field is about\n' ...
             '  |grad d|*|v_nom| / sqrt(Q_x/Q_u); at Q_u = 35, Q_x = 1,\n' ...
             '  |v_nom| = %.3g that is %.3f. Compare against bw = 0.3.\n\n'], ...
            v_ref, max(gnorm(:)) * v_ref / sqrt(1/35));
end
