%% sweep_paper_figs.m  --  runs for the paper figures (eps = 0.25)
%
%  Fig. A  final distance to the moving target vs wind intensity, with the
%          low-level MPC vs without low level (same double-integrator agents;
%          without low level each agent only tracks the commanded velocity
%          with its own velocity loop and does not know the wind)
%  Fig. B  Lyapunov functional V(t) of the same runs: exponential transient
%          plus a plateau set by the inputs -- the structure of Theorem 6.2
%  Fig. C  final snapshot with (eps = 0.25) and without (eps = 0) the
%          terminal-cost regularization
%
%  Needs run_doot_mpc.m (with the 'vel' option) and the updated primal_dual.m.
%  11 runs. Re-plot later without re-running:
%      parameters_w; load(out_file); plot_paper_figs(R, params, true)

parameters_w;          % NOTE: parameters_w.m calls clear, so it must come first
close all;

%% ------------------------------ settings --------------------------------
quick_test = false;              % 20-step runs to check paths; then set false
n_workers  = 0;                 % 0 = serial; >0 = parfor with that many workers
out_file   = 'paper_figs_results.mat';
export_pdf = true;

winds     = [0 0.01 0.03 0.05 0.1];   % wind amplitudes for Fig. A/B
U_snap    = 0.01;                      % wind of the Fig. C pair (must be in winds)
eps_noreg = 0;                         % "without regularization" (ridge kept)

base = struct();
base.T            = params.T;
base.epsilon_reg  = 0.25;
base.delta_ridge  = 4e-3;       % numerical ridge only; keeps eps = 0 well posed
base.pd_iters     = params.pd_iters;
base.target_omega = 2*pi/500;   % moving target: one revolution per 500 steps
base.n_inner      = 10;
base.Q_x = 1;  base.Q_u = 0.035;  base.Nh = 10;  base.gamma = 1;  base.u_max = [];
base.wind_model   = 'taylor_green';
base.wind_L       = 6;
base.W1_every     = 0;          % W1 at the last step only
base.W1_eps       = 0.05;       % smaller Sinkhorn blur than the default
base.W1_iters     = 2000;       % check W1_marg_err in the table (< 1e-3)

if quick_test
    base.T = 20;  base.pd_iters = 500;  base.W1_every = -1;
    fprintf('*** QUICK TEST MODE: T = %d, no W1 ***\n', base.T);
end

%% ------------------------------ run list --------------------------------
runs = {};  tags = {};
for U = winds
    o = base;  o.llc = 'lqr';  o.wind_U = U;
    runs{end+1} = o;  tags{end+1} = sprintf('MPC  U=%g', U);  %#ok<SAGROW>
    o = base;  o.llc = 'vel';  o.wind_U = U;
    runs{end+1} = o;  tags{end+1} = sprintf('noLL U=%g', U);  %#ok<SAGROW>
end
o = base;  o.llc = 'lqr';  o.wind_U = U_snap;  o.epsilon_reg = eps_noreg;
runs{end+1} = o;  tags{end+1} = 'MPC  eps=0';
fprintf('%d runs.\n', numel(runs));

%% ------------------------------ execute ---------------------------------
res = cell(numel(runs), 1);
t_all = tic;
parfor (r = 1:numel(runs), n_workers)
    out = run_doot_mpc(params, runs{r});
    res{r} = out;
    fprintf('%-12s done in %.1f min  (L1_ss = %.3f, W1 = %.3f, off-grid = %.3f)\n', ...
        tags{r}, out.runtime / 60, out.L1_ss, out.W1_final, out.off_grid);
end
fprintf('All runs done in %.1f min.\n', toc(t_all) / 60);

nU = numel(winds);
R = struct();
R.winds  = winds;
R.mpc    = res(1:2:2*nU);             % with low-level MPC, one per wind
R.vel    = res(2:2:2*nU);             % without low level, one per wind
R.reg    = R.mpc{find(winds == U_snap, 1)};
R.noreg  = res{end};
R.base   = base;
save(out_file, 'R', 'params', '-v7.3');

%% ------------------------------ plots -----------------------------------
plot_paper_figs(R, params, export_pdf);
