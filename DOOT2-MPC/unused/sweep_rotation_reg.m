%% sweep_rotation_reg.m  --  final distance vs target rotation speed, with
%  (eps = 0.25) and without (eps = 0) the terminal-cost regularization.
%
%  Same settings as the eps pair of sweep_paper_figs.m (low-level MPC,
%  Taylor-Green wind U = 0.01), so the 1 rev/500 steps points reproduce the
%  snapshot runs of that sweep. 12 runs.
%
%  Figure (with paper_figs_results.mat from sweep_paper_figs.m):
%      parameters_w; load('paper_figs_results.mat'); load('rotation_reg_results.mat');
%      plot_final_figure(R, Rrot, params, 'L1', true)

parameters_w;          % NOTE: parameters_w.m calls clear, so it must come first
close all;

%% ------------------------------ settings --------------------------------
quick_test = false;              % 20-step runs to check paths; then set false
n_workers  = 3;                 % 0 = serial; >0 = parfor with that many workers
out_file   = 'rotation_reg_results.mat';
revs       = [0 0.25 0.5 1 2 4];     % target revolutions per 500 steps
eps_pair   = [0.225 0];               % with / without regularization

base = struct();
base.T            = params.T;
base.delta_ridge  = 4e-3;       % numerical ridge only; keeps eps = 0 well posed
base.pd_iters     = params.pd_iters;
base.llc          = 'lqr';
base.n_inner      = 10;
base.Q_x = 1;  base.Q_u = 0.035;  base.Nh = 10;  base.gamma = 1;  base.u_max = [];
base.wind_model   = 'taylor_green';
base.wind_U       = 0.01;
base.wind_L       = 6;
base.W1_every     = 0;
base.W1_eps       = 0.05;
base.W1_iters     = 2000;

if quick_test
    base.T = 20;  base.pd_iters = 500;  base.W1_every = -1;
    fprintf('*** QUICK TEST MODE: T = %d, no W1 ***\n', base.T);
end

%% ------------------------------ run list --------------------------------
runs = {};  tags = {};
for e = eps_pair
    for rv = revs
        o = base;  o.epsilon_reg = e;  o.target_omega = 2*pi*rv/500;
        runs{end+1} = o;                                   %#ok<SAGROW>
        tags{end+1} = sprintf('eps=%g rev=%g', e, rv);     %#ok<SAGROW>
    end
end
fprintf('%d runs.\n', numel(runs));

%% ------------------------------ execute ---------------------------------
res = cell(numel(runs), 1);
t_all = tic;
parfor (r = 1:numel(runs), n_workers)
    out = run_doot_mpc(params, runs{r});
    res{r} = out;
    fprintf('%-16s done in %.1f min  (L1_ss = %.3f, W1 = %.3f, off-grid = %.3f)\n', ...
        tags{r}, out.runtime / 60, out.L1_ss, out.W1_final, out.off_grid);
end
fprintf('All runs done in %.1f min.\n', toc(t_all) / 60);

nr = numel(revs);
Rrot = struct();
Rrot.revs  = revs;
Rrot.eps   = eps_pair;
Rrot.reg   = res(1:nr);            % eps = eps_pair(1), one per rotation speed
Rrot.noreg = res(nr+1:2*nr);       % eps = eps_pair(2)
Rrot.base  = base;
save(out_file, 'Rrot', 'params', '-v7.3');
