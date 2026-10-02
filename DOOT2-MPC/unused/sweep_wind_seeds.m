%% sweep_wind_seeds.m  --  data for the wind panel (b), over several seeds
%
%  L1 error over time for several wind amplitudes, with the low-level MPC
%  and without low level (same double-integrator agents), each repeated for
%  independent initial swarms so the figure can show the mean and spread.
%  eps = 0.25, rotating target (1 rev per 500 steps). W1 skipped.
%  7 winds x 2 controllers x 3 seeds = 42 runs, ~1.5 h serial.

parameters_w;          % NOTE: parameters_w.m calls clear, so it must come first
close all;

%% ------------------------------ settings --------------------------------
quick_test = false;              % 20-step runs to check paths; then set false
n_workers  = 5;
out_file   = 'wind_seeds_results.mat';
winds      = [0 0.01 0.02 0.03 0.05 0.07 0.1];
seeds      = 1:3;

base = struct();
base.T            = params.T;
base.epsilon_reg  = 0.25;
base.delta_ridge  = 4e-3;
base.pd_iters     = params.pd_iters;
base.target_omega = 2*pi/500;
base.n_inner      = 10;
base.Q_x = 1;  base.Q_u = 0.035;  base.Nh = 10;  base.gamma = 1;  base.u_max = [];
base.wind_model   = 'taylor_green';
base.wind_L       = 6;
base.W1_every     = -1;

if quick_test
    base.T = 20;  base.pd_iters = 500;  winds = [0 0.1];  seeds = 1:2;
    fprintf('*** QUICK TEST MODE: T = %d ***\n', base.T);
end

nU = numel(winds);  nk = numel(seeds);
runs = cell(2 * nU * nk, 1);
llcs = {'lqr', 'vel'};
for c = 1:2
    for k = 1:nU
        for s = 1:nk
            o = base;  o.llc = llcs{c};  o.wind_U = winds(k);  o.seed = seeds(s);
            runs{sub2ind([nk, nU, 2], s, k, c)} = o;
        end
    end
end
fprintf('%d runs.\n', numel(runs));

res = cell(numel(runs), 1);
t_all = tic;
parfor (r = 1:numel(runs), n_workers)
    out = run_doot_mpc(params, runs{r});
    res{r} = out;
    fprintf('run %2d done in %.1f min  (%s, U = %g, seed = %d, L1_ss = %.3f)\n', r, ...
        out.runtime / 60, out.settings.llc, out.settings.wind_U, out.settings.seed, out.L1_ss);
end
fprintf('All runs done in %.1f h.\n', toc(t_all) / 3600);

res   = reshape(res, nk, nU, 2);
Rwind = struct('winds', winds, 'seeds', seeds, 'base', base);
Rwind.mpc = permute(res(:, :, 1), [2 1]);     % nU x nk
Rwind.vel = permute(res(:, :, 2), [2 1]);
save(out_file, 'Rwind', 'params', '-v7.3');
