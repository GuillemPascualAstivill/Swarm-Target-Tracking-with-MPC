function out = run_doot_mpc(params, opts)
%RUN_DOOT_MPC  One closed-loop run of the two-level controller (function
%   version of main_w_final.m) returning per-step metrics for sweeps.
%
%   out = run_doot_mpc(params, opts)
%
%   params : struct built by parameters_w.m (initial swarm params.pos_cont0,
%            grid, star geometry, bandwidths, primal-dual step sizes).
%   opts   : optional overrides. Omitted fields use the defaults below,
%            which reproduce main_w_final.m (plus the optional ridge).
%
%   Planner
%     .T             outer steps                                params.T
%     .epsilon_reg   regularization eps in (12)                 params.epsilon_reg
%     .delta_ridge   ridge (delta/2)*||alpha||^2 in (12)        params.delta_ridge or 0
%     .pd_iters      primal-dual iterations per step            params.pd_iters
%   Target
%     .target_omega  star rotation per outer step [rad], 0 = static
%                                                               params.target_omega
%   Low level
%     .llc  'lqr'  : the low_level_mpc law (constant gain from the Riccati
%                    recursion with terminal cost gamma*P*), vectorized over
%                    agents -- same control law and same integrator as
%                    main_w_final.m, much faster
%           'qp'   : low_level_mpc_qp per agent (input limits via u_max)
%           'vel'  : NO low level, same double-integrator agents: each agent
%                    tracks the commanded velocity with its own velocity loop
%                    u = -k_vel*(v - v_doot), with no position reference and no
%                    knowledge of the wind (ground velocity = v + d)
%           'none' : first-order DOOT, agents move with v_doot (+ wind)
%     .n_inner, .Q_x, .Q_u, .Nh, .gamma, .u_max           10, 1, 35, 10, 1, []
%     .k_vel         velocity gain for 'vel' (default: the velocity gain of
%                    the low-level MPC with the same Q_x, Q_u, so both
%                    controllers have the same velocity bandwidth)
%           (u_max: scalar or 1x2 row; [] = unconstrained)
%   Disturbance
%     .wind_model    'none' | 'uniform' | 'taylor_green'        'taylor_green'
%     .wind_U, .wind_theta, .wind_L                             0.01, 0, 6
%   Output
%     .W1_every      Sinkhorn W1 every k steps (0 = last step only,
%                    -1 = never; each evaluation takes ~1-2 min)  0
%     .W1_eps        Sinkhorn regularization (computeW2 default ~0.11 adds a
%                    blur of ~0.2 to every W1; 0.05 halves it)  computeW2 default
%     .W1_iters      Sinkhorn iterations (raise when W1_eps is small)  500
%     .snap_steps    steps whose agent positions are stored     5 steps
%     .r_close       near-collision distance for close_frac      0.05
%   Initial swarm
%     .seed          [] = params.pos_cont0; an integer draws an independent
%                    initial swarm from N(init_mean, init_std^2 I), the same
%                    distribution as parameters_w.m (local RandStream, the
%                    global random state is untouched)
%     .init_mean, .init_std                                     [0 0], sqrt(3)
%     .sat_tol       |grad phi| >= sat_tol counts as saturated  0.98
%     .tail_frac     final fraction of the run = steady state   0.4
%     .verbose       progress every 50 steps                    false
%
%   Per-step series (T x 1) in out:
%     V            planner optimum = Lyapunov functional V(rho_t) of Thm 6.2
%     Delta        mean_i |v_i - v*(x_i)|^2, v_i = ground velocity - v_doot.
%                  This is Delta_t of Thm 6.2 with the identity coupling,
%                  exact while no agent saturates (sigma = rho); check frac_sat.
%     L1           ||rho_N - rho*||_1 after the move (grid-normalized KDE)
%     W1           Sinkhorn W1 against the smoothed target (NaN if not computed)
%     max_gradphi, frac_sat, mean_cmd_speed, track_rms (= e^tr), effort
%     nn_dist      mean nearest-neighbour distance (clumping detector)
%     nn_med       median nearest-neighbour distance (robust to stray agents)
%     close_frac   fraction of agents with a neighbour closer than r_close
%   Summaries: V_ss, Delta_ss, L1_ss, track_ss, max_gradphi_ss (tail
%     means; max_gradphi_ss = steady max commanded speed), effort_total,
%     L1_rise = (L1_ss - min L1)/min L1 and nn_drop = relative fall of the
%     nearest-neighbour distance after the best fit (both > 0 = crumbling),
%     t_conv (time to 10% of the initial excess L1), inv_lambda (low-level
%     time constant), target_speed (omega*r0 per unit time), off_grid, snap.

if nargin < 2 || isempty(opts), opts = struct(); end
t_start = tic;
add_local_paths();

% ------------------------------ options ------------------------------------
T        = getopt(opts, 'T', params.T);
dt_outer = params.dt_cont;
N        = params.N;
if isfield(params, 'target_rot') && params.target_rot
    omega0 = params.target_omega;
else
    omega0 = 0;
end
target_omega = getopt(opts, 'target_omega', omega0);

params.epsilon_reg = getopt(opts, 'epsilon_reg', params.epsilon_reg);
params.delta_ridge = getopt(opts, 'delta_ridge', getopt(params, 'delta_ridge', 0));
params.pd_iters    = getopt(opts, 'pd_iters', params.pd_iters);

llc      = lower(getopt(opts, 'llc', 'lqr'));
n_inner  = getopt(opts, 'n_inner', 10);
dt_inner = dt_outer / n_inner;
Q_x      = getopt(opts, 'Q_x', 1);
Q_u      = getopt(opts, 'Q_u', 35);
Nh       = getopt(opts, 'Nh', 10);
gam      = getopt(opts, 'gamma', 1);
u_max    = getopt(opts, 'u_max', []);

wind_model = lower(getopt(opts, 'wind_model', 'taylor_green'));
wind_U     = getopt(opts, 'wind_U', 0.01);
wind_theta = getopt(opts, 'wind_theta', 0);
wind_L     = getopt(opts, 'wind_L', 6);

W1_every   = getopt(opts, 'W1_every', 0);
W1_eps     = getopt(opts, 'W1_eps', []);
W1_iters   = getopt(opts, 'W1_iters', 500);
W1_marg    = NaN;
r_close    = getopt(opts, 'r_close', 0.05);
snap_steps = getopt(opts, 'snap_steps', unique(max(1, round([1, T/20, T/5, T/2, T]))));
sat_tol    = getopt(opts, 'sat_tol', 0.98);
tail_frac  = getopt(opts, 'tail_frac', 0.4);
verbose    = getopt(opts, 'verbose', false);

% ------------------------------ low level ----------------------------------
switch llc
    case 'lqr'
        mpc = low_level_mpc_setup(dt_inner, Q_x, Q_u, Nh, gam, u_max);
        lam_llc = -log(mpc.rho_cl) / dt_inner;
    case 'qp'
        mpc = low_level_mpc_setup_qp(dt_inner, Q_x, Q_u, Nh, gam, u_max);
        lam_llc = -log(mpc.rho_cl) / dt_inner;
    case 'vel'
        mpc = low_level_mpc_setup(dt_inner, Q_x, Q_u, Nh, gam, u_max);   % for the gain only
        k_vel = getopt(opts, 'k_vel', mpc.K(2));
        lam_llc = -log(abs(1 - dt_inner * k_vel)) / dt_inner;
    case 'none'
        mpc = [];
        lam_llc = Inf;
    otherwise
        error('run_doot_mpc:llc', 'llc must be ''lqr'', ''qp'', ''vel'' or ''none''.');
end

% ------------------------------ wind ---------------------------------------
if wind_U == 0, wind_model = 'none'; end
switch wind_model
    case 'none'
        wind_fn = @(tt, p) zeros(size(p, 1), 2);
    case 'uniform'
        need_function('uniform_wind');
        wind_fn = @(tt, p) uniform_wind(tt, p, wind_U, wind_theta);
    case 'taylor_green'
        need_function('taylor_green_wind');
        wind_fn = @(tt, p) taylor_green_wind(tt, p, wind_U, wind_L);
    otherwise
        error('run_doot_mpc:wind', 'wind_model must be ''none'', ''uniform'' or ''taylor_green''.');
end

% ------------------------------ storage ------------------------------------
targetParams = params;              % carries epsilon, ridge and pd_iters
targetParams.rho_target_sm = [];    % force a rebuild if smoothing is ever on

V = nan(T, 1);  Delta = nan(T, 1);  L1 = nan(T, 1);  W1 = nan(T, 1);
max_gradphi = nan(T, 1);  frac_sat = nan(T, 1);  mean_cmd_speed = nan(T, 1);
track_rms = nan(T, 1);  effort = zeros(T, 1);  nn_dist = nan(T, 1);
nn_med = nan(T, 1);  close_frac = nan(T, 1);
n_snap = numel(snap_steps);
snap_pos = zeros(N, 2, n_snap);  snap_phase = zeros(1, n_snap);

pos   = params.pos_cont0;
seed  = getopt(opts, 'seed', []);
if ~isempty(seed)
    rs = RandStream('mt19937ar', 'Seed', seed);
    pos(1:N, :) = getopt(opts, 'init_mean', [0 0]) + ...
        getopt(opts, 'init_std', sqrt(3)) * randn(rs, N, 2);
end
vel   = zeros(N, 2);
dxdy  = params.dx * params.dy;
Omega = target_omega / dt_outer;    % angular velocity per unit time

% ------------------------------ closed loop --------------------------------
for t = 1:T
    phase = target_omega * (t - 1);
    targetParams.rho_target = star_target_w(params, phase);

    % Top level: regularized OT planner, target frozen over the interval
    [v_doot, V(t), gmag] = plan_step(pos, targetParams);
    max_gradphi(t)    = max(gmag);
    frac_sat(t)       = mean(gmag >= sat_tol);
    mean_cmd_speed(t) = mean(gmag);

    % Low level: track the ramp p_start + tau*v_doot over n_inner substeps
    p_start  = pos(1:N, :);
    track_sq = 0;
    for s = 1:n_inner
        tau      = (s - 1) * dt_inner;
        time_now = (t - 1) * dt_outer + tau;
        P = pos(1:N, :);
        d = wind_fn(time_now, P);
        d = d(:, 1:2);
        switch llc
            case 'none'
                ground = v_doot + d;
                pos(1:N, :) = P + dt_inner * ground;
                vel = v_doot;
            case 'lqr'
                ground = vel + d;
                p_ref  = p_start + tau * v_doot;
                u = -mpc.K(1) * (P - p_ref) - mpc.K(2) * (vel - (v_doot - d));
                if ~isempty(mpc.u_max)
                    u = max(min(u, mpc.u_max), -mpc.u_max);
                end
                pos(1:N, :) = P + dt_inner * (vel + d) + 0.5 * dt_inner^2 * u;
                vel = vel + dt_inner * u;
                effort(t) = effort(t) + sum(u(:).^2) * dt_inner;
            case 'vel'
                ground = vel + d;
                u = -k_vel * (vel - v_doot);           % no reference, no wind info
                if ~isempty(u_max)
                    u = max(min(u, u_max), -u_max);
                end
                pos(1:N, :) = P + dt_inner * (vel + d) + 0.5 * dt_inner^2 * u;
                vel = vel + dt_inner * u;
                effort(t) = effort(t) + sum(u(:).^2) * dt_inner;
            case 'qp'
                ground = vel + d;
                p_ref  = p_start + tau * v_doot;
                for m = 1:N
                    [u_m, v_next, ~, info] = low_level_mpc_qp(P(m, :), vel(m, :), ...
                        p_ref(m, :), v_doot(m, :), time_now, wind_fn, mpc);
                    pos(m, :) = P(m, :) + dt_inner * (vel(m, :) + info.d) + 0.5 * dt_inner^2 * u_m;
                    vel(m, :) = v_next;
                    effort(t) = effort(t) + sum(u_m.^2) * dt_inner;
                end
        end
        err = ground - v_doot;                       % ground velocity - command
        track_sq = track_sq + sum(err(:).^2) * dt_inner;
    end
    track_rms(t) = sqrt(track_sq / (N * dt_outer));

    % Theorem 6.2 input: perturbation v of the nominal velocity vs target v*
    v_pert = (pos(1:N, :) - p_start) / dt_outer - v_doot;
    v_star = Omega * [-p_start(:, 2), p_start(:, 1)];      % rigid rotation
    Delta(t) = mean(sum((v_pert - v_star).^2, 2));

    % Density error after the move, against this interval's target
    rho = agent_density(pos(1:N, :), params);
    L1(t) = sum(abs(rho - targetParams.rho_target)) * dxdy;

    % Clumping detector: mean nearest-neighbour distance (the L1 above is
    % computed after smoothing at h_agent and cannot see clumps below that)
    Pn = pos(1:N, :);
    d2 = (Pn(:, 1) - Pn(:, 1)').^2 + (Pn(:, 2) - Pn(:, 2)').^2;
    d2(1:N+1:end) = Inf;
    nn_i = sqrt(min(d2, [], 2));
    nn_dist(t)    = mean(nn_i);
    nn_med(t)     = median(nn_i);          % robust to a few stray agents
    close_frac(t) = mean(nn_i < r_close);  % agents with a near-collision
    if (W1_every > 0 && mod(t, W1_every) == 0) || (t == T && W1_every >= 0)
        tgt_sm = smooth_on_grid_wendland(targetParams.rho_target, params, params.h_agent);
        if isempty(W1_eps)
            [W1(t), w1info] = computeW2(rho, tgt_sm, params.x, params.y, 'p', 1, ...
                'iters', W1_iters);
        else
            [W1(t), w1info] = computeW2(rho, tgt_sm, params.x, params.y, 'p', 1, ...
                'epsilon', W1_eps, 'iters', W1_iters);
        end
        W1_marg = w1info.marg_err;
    end

    k_snap = find(snap_steps == t, 1);
    if ~isempty(k_snap)
        snap_pos(:, :, k_snap) = pos(1:N, :);
        snap_phase(k_snap) = phase;
    end

    if verbose && (mod(t, 50) == 0 || t == T)
        fprintf('  step %d/%d  L1 = %.4f  V = %.3e  max|grad phi| = %.3f\n', ...
            t, T, L1(t), V(t), max_gradphi(t));
    end
end

% ------------------------------ summaries ----------------------------------
tail  = max(1, floor((1 - tail_frac) * T) + 1):T;
L1_ss = mean(L1(tail));
k_conv = find(L1 <= L1_ss + 0.1 * (L1(1) - L1_ss), 1);
if isempty(k_conv) || L1(1) <= L1_ss
    t_conv = NaN;
else
    t_conv = k_conv * dt_outer;
end

out = struct();
out.opts  = opts;
out.t     = (1:T)' * dt_outer;
out.V     = V;
out.Delta = Delta;
out.L1    = L1;
out.W1    = W1;
out.max_gradphi    = max_gradphi;
out.frac_sat       = frac_sat;
out.mean_cmd_speed = mean_cmd_speed;
out.track_rms      = track_rms;
out.effort         = effort;
out.lambda_llc     = lam_llc;
out.inv_lambda     = 1 / lam_llc;
out.target_speed   = abs(Omega) * params.star_r0;
out.V_ss      = mean(V(tail));
out.Delta_ss  = mean(Delta(tail));
out.L1_ss     = L1_ss;
out.track_ss  = sqrt(mean(track_rms(tail).^2));
out.effort_total    = sum(effort);
[L1_min, k_min] = min(L1);
out.nn_dist  = nn_dist;
out.nn_ss    = mean(nn_dist(tail));
out.nn_med   = nn_med;
out.close_frac    = close_frac;
out.nn_med_ss     = mean(nn_med(tail));
out.close_frac_ss = mean(close_frac(tail));
out.L1_min   = L1_min;
out.L1_rise  = (L1_ss - L1_min) / L1_min;                   % > 0: worse than its best fit
out.nn_drop  = (nn_dist(k_min) - out.nn_ss) / nn_dist(k_min); % > 0: clumping after best fit
out.max_gradphi_max = max(max_gradphi);
out.max_gradphi_ss  = mean(max_gradphi(tail));   % steady max commanded speed
out.frac_sat_mean   = mean(frac_sat);
out.t_conv    = t_conv;
out.W1_final  = W1(T);
out.W1_marg_err = W1_marg;   % Sinkhorn marginal error of the last W1 (want < 1e-3)
out.off_grid  = mean(pos(1:N, 1) < params.x(1) | pos(1:N, 1) > params.x(end) | ...
                     pos(1:N, 2) < params.y(1) | pos(1:N, 2) > params.y(end));
out.snap      = struct('step', snap_steps, 'pos', snap_pos, 'phase', snap_phase);
out.pos_final = pos(1:N, :);
out.settings  = struct('T', T, 'epsilon_reg', params.epsilon_reg, ...
    'delta_ridge', params.delta_ridge, 'pd_iters', params.pd_iters, ...
    'target_omega', target_omega, 'llc', llc, 'n_inner', n_inner, ...
    'Q_x', Q_x, 'Q_u', Q_u, 'Nh', Nh, 'gamma', gam, 'wind_model', wind_model, ...
    'wind_U', wind_U, 'wind_L', wind_L);
out.settings.u_max = u_max;
out.settings.seed  = seed;
out.runtime   = toc(t_start);
end

% =========================================================================
function [v, Vt, gmag] = plan_step(pos, P)
% Same planner algebra as plan_doot in main_w_final.m, plus the objective.
N = P.N;  N_tot = P.N_tot;
kernels  = computeKernels_w(pos, P.ant_cords, P);
integral = computeTargetint_w(pos, P.ant_cords, P);

w_agent = zeros(N_tot, 1);
w_agent(1:N) = sum(integral.agents_ag(1:N, 1:N), 2) / N - integral.target_ag(1:N);
w_agent(N+1:N_tot) = sum(integral.agents_ag(N+1:N_tot, N+1:N_tot), 2) + P.ghost_w;
if ~isempty(integral.agents_ant)
    % single 1/N (main_w_final.m divides twice; only matters with antennas)
    w_ant = P.antw * (sum(integral.agents_ant, 1)' / N - integral.target_ant);
else
    w_ant = [];
end
w_cont = P.w_gain * [w_agent; w_ant];
w_cont(N+1:N_tot) = w_cont(N+1:N_tot) + P.ghost_w;

c_sq  = conformal(pos, P.ant_cords, P);
alpha = primal_dual(kernels, integral, w_cont, c_sq, P);

% Optimal value of the (ridge-)regularized dual = Lyapunov functional V
Vt = alpha' * w_cont - 0.5 * P.epsilon_reg * (alpha' * (integral.K_reg * alpha)) ...
     - 0.5 * P.delta_ridge * (alpha' * alpha);

G = [kernels.D_x(1:N, :) * alpha, kernels.D_y(1:N, :) * alpha];   % grad phi
v = -G;
gmag = sqrt(sum(G.^2, 2));
end

function rho = agent_density(pos, params)
% Wendland KDE on the grid at h_agent, normalized to unit mass on the grid.
h   = params.h_agent;
gp  = params.grid_pts;
rho = zeros(size(gp, 1), 1);
for first = 1:2000:numel(rho)
    idx = first:min(first + 1999, numel(rho));
    d2 = (gp(idx, 1) - pos(:, 1)').^2 + (gp(idx, 2) - pos(:, 2)').^2;
    rho(idx) = mean(wendland_kernel(d2, h), 2);
end
mass = sum(rho) * params.dx * params.dy;
if mass > 0
    rho = rho / mass;
end
end

function v = getopt(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function add_local_paths()
% main_w_final.m keeps the MPC and WINDS folders next to the script folder.
here   = fileparts(mfilename('fullpath'));
parent = fileparts(here);
subs   = {'MPC', 'WINDS'};
for k = 1:numel(subs)
    p = fullfile(parent, subs{k});
    if exist(p, 'dir') == 7
        addpath(p);
    end
end
end

function need_function(name)
if exist(name, 'file') ~= 2
    error('run_doot_mpc:path', ...
        '%s.m not found: add the WINDS folder to the MATLAB path.', name);
end
end
