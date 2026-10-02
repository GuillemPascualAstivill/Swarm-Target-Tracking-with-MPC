function fh = plot_diagnostics_w(diagn, params, tag)
%PLOT_DIAGNOSTICS_W  Convergence / constraint diagnostics for one Wendland run.
%
%   fh = plot_diagnostics_w(diagn, params, tag)
%
%   diagn is the struct main_w.m (or sweep_w.m's run_sim) fills in per time
%   step. The panels are meant to answer three questions:
%
%     1. Is the swarm actually settling, or orbiting the target?
%        -> panel 1 (L1/L2 vs t). A flat tail = settled; a ripple = orbiting,
%           and then a single end-of-run snapshot is a luck-of-the-draw metric.
%
%     2. Is alpha converging across time steps, or still marching?
%        -> panels 2 and 3. With a COLD start the elementwise clip inside
%           primal_dual caps max|alpha_j| at etap*pd_iters (dashed line in
%           panel 3): if alpha_inf sits on that line, pd_iters is acting as a
%           step-size gain, not a convergence knob, and comparing pd_iters = 30
%           against 50 compares two different step sizes.
%
%     3. Is the Kantorovich constraint |grad phi| <= 1 being enforced?
%        -> panel 4. Sustained values above the dashed line at 1 mean the
%           lambda updates never catch up within pd_iters.

if nargin < 3 || isempty(tag), tag = ''; end

T = numel(diagn.err_L1);
tt = 1:T;
tail_idx = max(1, round(0.9*T)):T;

fh = figure('Units','normalized','Position',[0.06 0.08 0.80 0.78], ...
            'Name', sprintf('Diagnostics  (pd\\_iters = %d, warm\\_start = %d)', ...
                            diagn.pd_iters, diagn.warm_start), ...
            'NumberTitle','off');

% ---- 1. distribution error ---------------------------------------------
subplot(2,3,1); hold on; box on;
plot(tt, diagn.err_L1, 'LineWidth', 1.4);
plot(tt, diagn.err_L2, 'LineWidth', 1.4);
xlabel('step $t$'); ylabel('error');
legend({'$L^1$','$L^2$'}, 'Location','northeast');
title({'Error vs target', ...
       sprintf('final $L^1$ = %.4f, tail mean = %.4f, tail std = %.4f', ...
               diagn.err_L1(T), mean(diagn.err_L1(tail_idx)), std(diagn.err_L1(tail_idx)))});

% ---- 2. outer convergence of alpha --------------------------------------
subplot(2,3,2); hold on; box on;
semilogy(tt, max(diagn.d_alpha, eps), 'LineWidth', 1.2);
semilogy(tt, max(diagn.rel_d_alpha, eps), 'LineWidth', 1.2);
set(gca, 'YScale', 'log');
xlabel('step $t$'); ylabel('change per step');
legend({'$\|\alpha_t-\alpha_{t-1}\|$','relative'}, 'Location','northeast');
title('Outer convergence of $\alpha$');

% ---- 3. alpha magnitude vs the cold-start cap ---------------------------
subplot(2,3,3); hold on; box on;
plot(tt, diagn.alpha_inf, 'LineWidth', 1.4);
cap = diagn.alpha_cap_cold;
plot([1 T], [cap cap], 'r--', 'LineWidth', 1.2);
xlabel('step $t$'); ylabel('$\max_j |\alpha_j|$');
legend({'$\|\alpha\|_\infty$', ...
        sprintf('cold-start cap $\\eta_p \\cdot$ pd\\_iters = %.3f', cap)}, ...
       'Location','southeast');
title({'$\alpha$ magnitude', ...
       'sitting on the cap $\Rightarrow$ pd\_iters is a gain'});

% ---- 4. constraint satisfaction -----------------------------------------
subplot(2,3,4); hold on; box on;
yyaxis left;
plot(tt, diagn.max_gradphi, 'LineWidth', 1.4);
plot([1 T], [1 1], 'k--', 'LineWidth', 1.2);
ylabel('$\max |\nabla \phi|$');
yyaxis right;
plot(tt, 100*diagn.frac_violate, 'LineWidth', 1.2);
ylabel('\% agents with $|\nabla \phi| > 1$');
xlabel('step $t$');
title({'Kantorovich constraint', ...
       sprintf('tail mean $\\max|\\nabla\\phi|$ = %.3f', ...
               mean(diagn.max_gradphi(tail_idx)))});

% ---- 5. speeds -----------------------------------------------------------
subplot(2,3,5); hold on; box on;
plot(tt, diagn.mean_speed, 'LineWidth', 1.4);
plot(tt, diagn.max_speed,  'LineWidth', 1.0);
xlabel('step $t$'); ylabel('speed');
legend({'mean','max'}, 'Location','northeast');
title('Agent speed (per step)');

% ---- 6. mismatch residual ------------------------------------------------
subplot(2,3,6); hold on; box on;
plot(tt, diagn.w_norm, 'LineWidth', 1.4);
xlabel('step $t$'); ylabel('$\|w\|$');
title({'Mismatch residual', sprintf('final = %.3e', diagn.w_norm(T))});

% ---- super-title ---------------------------------------------------------
if diagn.warm_start && diagn.warm_start_lambda
    ws = 'warm start ($\alpha$ + $\lambda$)';
elseif diagn.warm_start
    ws = 'warm start ($\alpha$ only)';
else
    ws = 'cold start';
end
sup = sprintf('Diagnostics --- %s, pd\\_iters $= %d$, $\\eta_p = %g$, $\\eta_d = %g$, $T = %d$', ...
              ws, diagn.pd_iters, params.etap, params.etad, T);
if ~isempty(tag)
    sup = sprintf('%s \\quad [%s]', sup, tag);
end
if exist('sgtitle', 'file') == 2 || exist('sgtitle', 'builtin') == 5
    sgtitle(fh, sup, 'Interpreter', 'latex', 'FontSize', 14);
else
    annotation(fh, 'textbox', [0 0.955 1 0.045], 'String', sup, ...
        'Interpreter', 'latex', 'FontSize', 14, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none');
end
drawnow;
end
