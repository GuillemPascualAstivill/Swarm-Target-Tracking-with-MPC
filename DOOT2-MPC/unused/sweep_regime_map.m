%% sweep_regime_map.m  --  data for the regime map of panel (a)
%
%  Steady tracking error and clustering over a grid of regularization eps
%  and target rotation speed, averaged over independent initial swarms
%  (low-level MPC, Taylor-Green wind U = 0.01, same settings as the other
%  paper runs). W1 is skipped: the map uses the L1 error and nearest-
%  neighbour distances, so each run takes ~2 min.
%
%  Full grid: 10 eps x 8 speeds x 3 seeds = 240 runs, ~9 h serial (divide
%  by n_workers with the Parallel Computing Toolbox). Results are saved
%  after every row of speeds, so an interrupted run keeps what finished.
%
%  Figure:
%      parameters_w; load('paper_figs_results.mat'); load('regime_map_results.mat');
%      load('wind_seeds_results.mat');
%      plot_final_figure(R, [], params, 'Rmap', Rmap, 'Rwind', Rwind, 'export', true)

parameters_w;          % NOTE: parameters_w.m calls clear, so it must come first
close all;

%% ------------------------------ settings --------------------------------
quick_test = false;              % tiny grid, 20-step runs; then set false
n_workers  = 5;                 % 0 = serial; >0 = parfor with that many workers
out_file   = 'regime_map_results.mat';
eps_list   = [0 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4];   % regularization
revs       = [0 0.25 0.5 1 1.5 2 3 4];                        % rev per 500 steps
seeds      = 1:3;                                             % initial swarms

base = struct();
base.T            = params.T;
base.delta_ridge  = 4e-3;
base.pd_iters     = params.pd_iters;
base.llc          = 'lqr';
base.n_inner      = 10;
base.Q_x = 1;  base.Q_u = 0.035;  base.Nh = 10;  base.gamma = 1;  base.u_max = [];
base.wind_model   = 'taylor_green';
base.wind_U       = 0.01;
base.wind_L       = 6;
base.W1_every     = -1;         % no Sinkhorn W1 (not used by the map)

if quick_test
    base.T = 20;  base.pd_iters = 500;
    eps_list = [0 0.25];  revs = [0 1];  seeds = 1:2;
    fprintf('*** QUICK TEST MODE: T = %d, 2x2x2 grid ***\n', base.T);
end

ne = numel(eps_list);  ns = numel(revs);  nk = numel(seeds);
Rmap = struct('eps', eps_list, 'revs', revs, 'seeds', seeds, 'base', base);
Rmap.res = cell(ns, ne, nk);    % rows: speed, cols: eps, pages: seed
fprintf('%d runs.\n', ns * ne * nk);
t_all = tic;

%% ------------------------------ execute (row by row) ---------------------
for i = 1:ns
    runs = cell(ne * nk, 1);
    for j = 1:ne
        for k = 1:nk
            o = base;
            o.epsilon_reg  = eps_list(j);
            o.target_omega = 2*pi*revs(i)/500;
            o.seed         = seeds(k);
            runs{(k - 1) * ne + j} = o;
        end
    end
    res = cell(ne * nk, 1);
    parfor (r = 1:ne * nk, n_workers)
        res{r} = run_doot_mpc(params, runs{r});
    end
    Rmap.res(i, :, :) = reshape(res, 1, ne, nk);
    save(out_file, 'Rmap', 'params', '-v7.3');          % checkpoint
    fprintf('speed %d/%d (%.2f rev/500 steps) done, %.1f h elapsed\n', ...
        i, ns, revs(i), toc(t_all) / 3600);
end
fprintf('All runs done in %.1f h.\n', toc(t_all) / 3600);
