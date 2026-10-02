%% sweep_doot_mpc.m  --  parameter sweeps for the paper figures
%
% Runs run_doot_mpc.m (function version of main_w_final.m) over one-at-a-time
% sweeps, saves everything to a .mat file and plots one figure per sweep plus
% a snapshot figure.
%
% Files needed next to main_w_final.m: this script, run_doot_mpc.m,
% plot_sweep_doot_mpc.m, plot_snapshots_doot_mpc.m and the updated
% primal_dual.m (identical to the old one when delta_ridge = 0).
%
% 1) Run with quick_test = true (short runs) to check paths and plots.
% 2) Set quick_test = false for the real sweeps. Each run costs T planner
%    solves; use n_workers > 0 (Parallel Computing Toolbox) or run the
%    sweeps separately with run_sweeps.
% 3) To re-plot without re-running:  load(out_file) and call
%    plot_sweep_doot_mpc(results, meta, sweeps, params, export_pdf)

parameters_w;          % NOTE: parameters_w.m calls clear, so it must come first
close all;

%% ------------------------------ settings --------------------------------
quick_test = false;             % short runs to check that everything works
n_workers  = 0;                % 0 = serial; >0 = parfor with that many workers
run_sweeps = 2:4;                % stage 1: eps sweep; stage 2: 2:4 (5 = ridge check)
out_file   = sprintf('sweep_results_%s.mat', sprintf('%d', run_sweeps));
export_pdf = false;            % one vector PDF per figure

base = struct();               % shared defaults (= main_w_final.m + ridge)
base.T            = params.T;
base.epsilon_reg  = 0.25;       % replica: tracks better than 0.25, no crumbling (0.05 crumbles)
base.delta_ridge  = 4e-3;      % ridge on kernel weights (0 = original solver)
base.pd_iters     = params.pd_iters;
base.target_omega = 2*pi/500;  % one revolution per 500 steps
base.llc          = 'lqr';
base.n_inner      = 10;
base.Q_x = 1;  base.Q_u = 35;  base.Nh = 10;  base.gamma = 1;  base.u_max = [];
base.wind_model   = 'taylor_green';
base.wind_U       = 0.01;
base.wind_L       = 6;
base.W1_every     = 0;         % Sinkhorn W1 only at the last step (slow)

if quick_test                   % only checks that every code path runs
    base.T = 20;  base.pd_iters = 500;  base.W1_every = -1;
    fprintf('*** QUICK TEST MODE: T = %d, pd_iters = %d, no W1 ***\n', base.T, base.pd_iters);
end

%% ------------------------------ sweeps ----------------------------------
% Each sweep varies ONE field over 'values'. 'fixed' overrides base for the
% whole sweep; each entry of 'variants' repeats the sweep (one curve each).
mpc_variant  = struct('label', 'DOOT+MPC', 'llc', 'lqr');
doot_variant = struct('label', 'DOOT',     'llc', 'none');

sweeps = {};
% 1) Regularization -> tracking vs crumbling. RUN THIS FIRST, pick eps, set
%    base.epsilon_reg, then run sweeps 2-4. Smaller eps = higher gain 1/eps:
%    faster and closer tracking, but with the slow low level (Q_u = 35) the
%    swarm fits best early and then clumps. Pick the smallest eps whose
%    L1_rise and nn_drop (table, panels b and d) stay near zero. The
%    saturated transport regime (eps <= 0.02) also needs a smaller dt_cont.
sweeps{end+1} = make_sweep('epsilon', 'epsilon_reg', [0.3 0.25 0.2 0.15 0.1 0.075 0.05], ...
    '\varepsilon', 'regularization $\varepsilon$', struct('wind_U', 0), {mpc_variant});
% 2) Target speed -> lag of Remark 6.1 and the ISS estimate of Theorem 6.2
%    (wind off). Rigid rotation is divergence-free: the nu = 0 case.
%    Values = revolutions per 500 steps.
sweeps{end+1} = make_sweep('target_speed', 'target_omega', 2*pi*[0 0.25 0.5 1 2 4]/500, ...
    '\omega', 'target rotation $\omega$ (rad/step)', struct('wind_U', 0), {mpc_variant});
% 3) Wind amplitude -> disturbance rejection of the low level.
sweeps{end+1} = make_sweep('wind', 'wind_U', [0 0.01 0.03 0.1], ...
    'U', 'wind amplitude $U$', struct(), {doot_variant, mpc_variant});
% 4) Low-level bandwidth -> time-scale separation (plotted against 1/lambda).
sweeps{end+1} = make_sweep('lowlevel', 'Q_u', [0.0035 0.035 0.35 3.5 35], ...
    'Q_u', 'low-level weight $Q_u$', struct('wind_U', 0.03), {mpc_variant});
% 5) (optional) ridge sanity check: results should barely change in [1e-3, 1e-2].
sweeps{end+1} = make_sweep('ridge', 'delta_ridge', [0 1e-3 4e-3 1e-2], ...
    '\delta', 'ridge $\delta$', struct(), {mpc_variant});

%% ------------------------------ run list --------------------------------
runs = {};
meta = zeros(0, 3);            % [sweep index, value index, variant index]
for si = run_sweeps
    sw = sweeps{si};
    for vi = 1:numel(sw.variants)
        for k = 1:numel(sw.values)
            o = merge_struct(base, sw.fixed);
            o = merge_struct(o, rmfield(sw.variants{vi}, 'label'));
            o.(sw.field) = sw.values(k);
            runs{end+1} = o;                 %#ok<SAGROW>
            meta(end+1, :) = [si, k, vi];    %#ok<SAGROW>
        end
    end
end
n_runs = numel(runs);
fprintf('%d runs in %d sweeps.\n', n_runs, numel(run_sweeps));

%% ------------------------------ execute ---------------------------------
results = cell(n_runs, 1);
t_all = tic;
parfor (r = 1:n_runs, n_workers)
    res = run_doot_mpc(params, runs{r});
    results{r} = res;
    fprintf('run %d/%d done in %.0f s  (L1_ss = %.4f, max|grad phi| = %.3f)\n', ...
        r, n_runs, res.runtime, res.L1_ss, res.max_gradphi_max);
end
fprintf('All runs done in %.1f min.\n', toc(t_all) / 60);
save(out_file, 'results', 'runs', 'meta', 'sweeps', 'base', 'params', '-v7.3');

%% ------------------------------ plots -----------------------------------
plot_sweep_doot_mpc(results, meta, sweeps, params, export_pdf);

% Snapshot figure from the DOOT+MPC run with the base configuration
r_snap = find_run(runs, struct('llc', 'lqr', 'wind_U', base.wind_U, ...
    'target_omega', base.target_omega, 'epsilon_reg', base.epsilon_reg, ...
    'Q_u', base.Q_u, 'delta_ridge', base.delta_ridge));
if isempty(r_snap)
    fprintf(['No run matches the base configuration (it is part of sweep 3); ' ...
        'snapshots skipped. For them alone: res = run_doot_mpc(params, base);\n' ...
        'plot_snapshots_doot_mpc(res, params, true)\n']);
else
    plot_snapshots_doot_mpc(results{r_snap}, params, export_pdf);
end

%% ------------------------------ helpers ---------------------------------
function sw = make_sweep(name, field, values, symbol, xlab, fixed, variants)
sw = struct('name', name, 'field', field, 'values', values, ...
    'symbol', symbol, 'xlabel', xlab, 'fixed', fixed);
sw.variants = variants;      % assigned separately: a cell in struct() makes an array
end

function a = merge_struct(a, b)
f = fieldnames(b);
for k = 1:numel(f)
    a.(f{k}) = b.(f{k});
end
end

function r = find_run(runs, crit)
r = [];
f = fieldnames(crit);
for k = 1:numel(runs)
    ok = true;
    for j = 1:numel(f)
        if ~isfield(runs{k}, f{j}) || ~isequal(runs{k}.(f{j}), crit.(f{j}))
            ok = false;
            break;
        end
    end
    if ok
        r = k;
        return;
    end
end
end
