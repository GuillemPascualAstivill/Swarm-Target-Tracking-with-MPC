function figs = plot_sweep_doot_mpc(results, meta, sweeps, params, export_pdf)
%PLOT_SWEEP_DOOT_MPC  One figure per sweep run by sweep_doot_mpc.m.
%
%   (a) L1 density error vs time, color = swept value, line style = variant
%   (b) Lyapunov functional V(rho_t) vs time
%   (c) steady-state L1 error vs the swept value (one curve per variant);
%       for the target-speed sweep also the steady max commanded speed
%       against the target speed (below the identity line = cannot keep up)
%   (d) sweep-specific summary
%       target_omega : steady-state V vs Delta with the tightest affine
%                      bound V_0 + b*Delta (form of Theorem 6.2, b = 1/(2c))
%       wind_U       : low-level tracking RMS vs U
%       Q_u          : tracking RMS and control effort vs 1/lambda
%       epsilon_reg  : convergence time and degradation after the best fit
%   For the epsilon sweep, panel (b) shows the nearest-neighbour distance
%   instead of V (clumping detector; V changes definition with eps).
%   A summary table is printed to the command window.

if nargin < 5, export_pdf = false; end
styles = {'-', '--', ':', '-.'};
marks  = {'o', 's', '^', 'd'};
figs   = gobjects(0);

for si = unique(meta(:, 1))'
    sw   = sweeps{si};
    rows = find(meta(:, 1) == si);
    nval = numel(sw.values);
    nvar = numel(sw.variants);
    cols = parula(nval + 1);
    cols = cols(1:nval, :);                  % drop the pale yellow end

    fig = figure('Name', ['Sweep: ' sw.name], 'Color', 'w', ...
        'Units', 'normalized', 'Position', [0.05 0.08 0.8 0.75]);
    figs(end + 1) = fig; %#ok<AGROW>

    % ---------------- (a) L1 vs time, (b) V or saturation vs time -------
    h_val = gobjects(nval, 1);
    for panel = 1:2
        subplot(2, 2, panel); hold on; box on;
        for r = rows'
            res = results{r};
            k = meta(r, 2);
            v = meta(r, 3);
            if panel == 1
                y = res.L1;
            elseif strcmp(sw.field, 'epsilon_reg')
                y = res.nn_dist;
            else
                y = res.V;
            end
            h = plot(res.t, y, 'Color', cols(k, :), ...
                'LineStyle', styles{1 + mod(v - 1, 4)}, 'LineWidth', 1.2);
            if panel == 1 && v == 1
                h_val(k) = h;
            end
        end
        xlabel('time');
        if panel == 1
            ylabel('$\|\rho_N-\rho^*\|_{L^1}$');
            set(gca, 'YScale', 'log');
            title('(a) density error');
        elseif strcmp(sw.field, 'epsilon_reg')
            ylabel('mean nearest-neighbour distance');
            title('(b) clumping (falling = crumbling)');
        else
            ylabel('$V(\rho_t)$');
            set(gca, 'YScale', 'log');
            title('(b) Lyapunov functional');
        end
    end
    subplot(2, 2, 1);
    labels = arrayfun(@(x) sprintf('$%s = %.3g$', sw.symbol, x), sw.values, ...
        'UniformOutput', false);
    ok_h = isgraphics(h_val);
    hv = h_val(ok_h);
    lv = labels(ok_h);
    if nvar > 1
        hs = gobjects(nvar, 1);
        for v = 1:nvar
            hs(v) = plot(nan, nan, 'k', 'LineStyle', styles{1 + mod(v - 1, 4)});
        end
        var_labels = cellfun(@(s) s.label, sw.variants, 'UniformOutput', false);
        legend([hv(:); hs(:)], [lv(:); var_labels(:)], 'Location', 'best');
    else
        legend(hv(:), lv(:), 'Location', 'best');
    end

    % ---------------- x values for the summaries ------------------------
    x = nan(numel(rows), 1);
    for j = 1:numel(rows)
        r = rows(j);
        switch sw.field
            case 'Q_u'
                x(j) = results{r}.inv_lambda;
            case 'target_omega'
                x(j) = results{r}.target_speed;
            otherwise
                x(j) = sw.values(meta(r, 2));
        end
    end
    switch sw.field
        case 'Q_u',          xlab = 'low-level time constant $1/\lambda$';
        case 'target_omega', xlab = 'target speed $\omega r_0$';
        otherwise,           xlab = sw.xlabel;
    end

    % ---------------- (c) steady-state L1 vs swept value ----------------
    subplot(2, 2, 3); hold on; box on;
    for v = 1:nvar
        sel = meta(rows, 3) == v;
        xs  = x(sel);
        ys  = cellfun(@(res) res.L1_ss, results(rows(sel)));
        [xs, o] = sort(xs);
        plot(xs, ys(o), ['-' marks{1 + mod(v - 1, 4)}], 'LineWidth', 1.3, ...
            'DisplayName', sw.variants{v}.label);
    end
    ylabel('steady-state $L^1$ error');
    if strcmp(sw.field, 'Q_u')
        set(gca, 'XScale', 'log');
        hx = xline(params.dt_cont, 'k--', '$\Delta t$');
        hx.Interpreter = 'latex';
        hx.HandleVisibility = 'off';
    elseif strcmp(sw.field, 'epsilon_reg')
        set(gca, 'XScale', 'log');
    elseif strcmp(sw.field, 'target_omega')
        % Can the swarm keep up? Steady max commanded speed vs target speed:
        % below the dotted identity line the target outruns the swarm.
        yyaxis right;
        ms = cellfun(@(res) res.max_gradphi_ss, results(rows));
        [xs, o] = sort(x);
        plot(xs, ms(o), '-s', 'LineWidth', 1.1, 'DisplayName', 'max commanded speed');
        plot(xs, xs, ':', 'Color', [0.3 0.3 0.3], 'DisplayName', 'target speed');
        ylabel('speed');
        yyaxis left;
    end
    xlabel(xlab);
    legend('Location', 'best');
    title('(c) steady state');

    % ---------------- (d) sweep-specific summary ------------------------
    subplot(2, 2, 4); hold on; box on;
    switch sw.field
        case 'target_omega'
            % Theorem 6.2 at steady state: V_ss <= (Delta_ss + Delta_0)/(2c),
            % an AFFINE bound in the measured input. Delta_0 is the part of the
            % input the measurement cannot see (discretization, finite N); the
            % static-target run exposes it as the floor V_0. The dashed line
            % is the tightest such bound through that floor; its slope is an
            % empirical 1/(2c). Expect the data to bend below it (the bound
            % is linear, the lag saturates once the target outruns the swarm).
            Vss = cellfun(@(res) res.V_ss, results(rows));
            Dss = cellfun(@(res) res.Delta_ss, results(rows));
            plot(Dss, Vss, 'o', 'MarkerFaceColor', [0 .45 .74], 'DisplayName', 'runs');
            [D0, i0] = min(Dss);
            V0 = Vss(i0);
            above = Dss > D0;
            if any(above)
                b  = max((Vss(above) - V0) ./ (Dss(above) - D0));
                dd = linspace(D0, max(Dss), 50);
                plot(dd, V0 + b * (dd - D0), 'k--', ...
                    'DisplayName', sprintf('affine bound, slope %.3g', b));
                fprintf('[%s] tightest affine bound V_ss <= V_0 + b*Delta_ss: b = %.3g, i.e. c = %.3g\n', ...
                    sw.name, b, 1 / (2 * b));
            end
            xlabel('steady-state input $\bar\Delta$');
            ylabel('steady-state $\bar V$');
            legend('Location', 'best');
            title('(d) ISS estimate: $\bar V \le (\bar\Delta+\Delta_0)/(2c)$');
        case 'wind_U'
            for v = 1:nvar
                sel = meta(rows, 3) == v;
                xs  = x(sel);
                ys  = cellfun(@(res) res.track_ss, results(rows(sel)));
                [xs, o] = sort(xs);
                plot(xs, ys(o), ['-' marks{1 + mod(v - 1, 4)}], 'LineWidth', 1.3, ...
                    'DisplayName', sw.variants{v}.label);
            end
            xlabel(xlab);
            ylabel('ground-velocity tracking RMS');
            legend('Location', 'best');
            title('(d) tracking error of the commanded velocity');
        case 'Q_u'
            trk = cellfun(@(res) res.track_ss, results(rows));
            eff = cellfun(@(res) res.effort_total, results(rows));
            [xs, o] = sort(x);
            yyaxis left;
            plot(xs, trk(o), '-o', 'LineWidth', 1.3);
            ylabel('tracking RMS $e^{\mathrm{tr}}$');
            yyaxis right;
            plot(xs, eff(o), '-s', 'LineWidth', 1.3);
            ylabel('control effort $\int\|u\|^2 dt$');
            set(gca, 'XScale', 'log');
            hx = xline(params.dt_cont, 'k--');
            hx.HandleVisibility = 'off';
            xlabel(xlab);
            title('(d) tracking vs effort');
        case 'epsilon_reg'
            tc = cellfun(@(res) res.t_conv, results(rows));
            lr = cellfun(@(res) res.L1_rise, results(rows));
            [xs, o] = sort(x);
            yyaxis left;
            plot(xs, tc(o), '-o', 'LineWidth', 1.3);
            ylabel('time to 10\% of initial error');
            yyaxis right;
            plot(xs, 100 * lr(o), '-s', 'LineWidth', 1.3);
            ylabel('steady $L^1$ above its best (\%)');
            set(gca, 'XScale', 'log');
            xlabel(xlab);
            title('(d) speed vs degradation');
        otherwise
            ys = cellfun(@(res) res.track_ss, results(rows));
            [xs, o] = sort(x);
            plot(xs, ys(o), '-o', 'LineWidth', 1.3);
            xlabel(xlab);
            ylabel('tracking RMS');
            title('(d) tracking error');
    end

    sgtitle(fig, ['sweep: ' sw.name], 'Interpreter', 'none');
    if export_pdf
        exportgraphics(fig, sprintf('sweep_%s.pdf', sw.name), 'ContentType', 'vector');
    end
end

% ---------------- summary table ------------------------------------------
fprintf('\n%-13s %-9s %9s %9s %8s %8s %10s %10s %9s %9s %8s %9s %8s\n', 'sweep', 'variant', ...
    'value', 'L1_ss', 'L1_rise', 'nn_drop', 'V_ss', 'Delta_ss', 'track', 'maxgphi', 'sat', ...
    't_conv', 'W1_end');
for r = 1:size(meta, 1)
    sw  = sweeps{meta(r, 1)};
    res = results{r};
    fprintf('%-13s %-9s %9.3g %9.4f %8.3f %8.3f %10.3e %10.3e %9.4f %9.3f %8.3f %9.1f %8.4f\n', ...
        sw.name, sw.variants{meta(r, 3)}.label, sw.values(meta(r, 2)), res.L1_ss, ...
        res.L1_rise, res.nn_drop, res.V_ss, res.Delta_ss, res.track_ss, ...
        res.max_gradphi_max, res.frac_sat_mean, res.t_conv, res.W1_final);
end
end
