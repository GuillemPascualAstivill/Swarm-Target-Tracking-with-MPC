%% main_nodes.m  --  DOOT swarm control toward a rotating star target.
% The potential phi lives on a FIXED grid of nodes (build_nodes.m). Every DOOT
% step solves for its coefficients (plan_fixed.m) and gives each agent the
% velocity v = -grad phi, which the agent follows through the per-agent
% low-level MPC under wind (or directly). Run parameters_w.m first.

clc; close all;

if ~exist('params', 'var') || ~isfield(params, 'pos0')
    error('Run parameters_w.m before main_nodes.m.');
end

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, fullfile(fileparts(script_dir), 'MPC'), ...
    fullfile(fileparts(script_dir), 'WINDS'));

N = params.N;
T = params.T;
dt_outer = params.dt_cont;
dt_inner = dt_outer / params.n_inner;
dA = params.dx * params.dy;  % grid cell area
gr = params.grid_res;

switch params.wind_model
    case 'none',         wind_fn = @(tt,p) zeros(size(p));
    case 'uniform',      wind_fn = @(tt,p) uniform_wind(tt, p, params.wind_U, params.wind_theta);
    case 'taylor_green', wind_fn = @(tt,p) taylor_green_wind(tt, p, params.wind_U, params.wind_L);
    otherwise, error('main_nodes:windModel', ...
        'wind_model must be ''none'', ''uniform'' or ''taylor_green''.');
end

if params.use_mpc
    mpc = low_level_mpc_setup(dt_inner, params.Q_x, params.Q_u, params.Nh, ...
        params.gamma, params.u_max);
end

%% Fixed node basis for phi (built once)
nb = build_nodes(params);
M = size(nb.A, 1);
fprintf('Node basis: %d nodes, spacing %.2f, h_node = %.2f\n', M, ...
    (params.x(end) - params.x(1)) / (params.num_nodes_axis - 1), params.h_node);

%% Simulation
pos = params.pos0;           % agent positions
vel = zeros(N, 2);           % agent velocities (MPC only)
d_hat = zeros(N, 2);         % agents' wind estimates ('estimated' only)
% beta and lam are carried over between steps as the warm start (exact: the basis is fixed)
beta = zeros(M, 1);          % node coefficients, phi(x) = sum_a beta_a K_node(x - A_a)
lam  = zeros(N + M, 1);      % multipliers of |grad phi| <= 1, at the agents then the nodes

pos_hist    = zeros(N, 2, T);               % agent positions after each step
phi_hist    = zeros(gr*gr, T, 'single');    % phi on the grid
target_hist = zeros(gr, gr, T, 'single');   % target density on the grid
diagn = struct( ...
    'err_L1',    zeros(T,1), ...   % L1 distance between the swarm density and the target
    'frac_sat',  zeros(T,1), ...   % fraction of agents at the speed limit (on W1 rays)
    'max_speed', zeros(T,1), ...
    'lam_max',   zeros(T,1), ...   % largest agent multiplier
    'mass_T',    zeros(T,1), ...   % terminal mass
    'step_last', zeros(T,1), ...   % last primal step of the solver (small = converged)
    'collapsed', zeros(T,1));      % fraction of agents closer than 0.05 to another agent

for t = 1:T
    % Target: the star rotated by target_omega per step
    rho_s = star_target_w(params, params.target_omega * (t - 1));

    % DOOT: velocity v = -grad phi for every agent
    [v, beta, lam, info] = plan_fixed(pos, beta, lam, nb, rho_s, params);

    % Diagnostics, before the agents move
    speed = hypot(v(:,1), v(:,2));
    rho_N = info.rhoN / (sum(info.rhoN) * dA);    % swarm density with unit mass
    D = pdist2(pos, pos);  D(1:N+1:end) = inf;    % agent-agent distances
    diagn.err_L1(t)    = sum(abs(rho_N - rho_s)) * dA;
    diagn.frac_sat(t)  = mean(speed > 0.95);
    diagn.max_speed(t) = max(speed);
    diagn.lam_max(t)   = max(lam(1:N));
    diagn.mass_T(t)    = info.mass_T;
    diagn.step_last(t) = info.step_last;
    diagn.collapsed(t) = mean(min(D, [], 2) < 0.05);
    phi_hist(:, t)       = info.phi_grid;
    target_hist(:, :, t) = reshape(rho_s, gr, gr);

    % Motion
    if params.use_mpc
        % Each agent tracks the reference p_start + tau*v over n_inner MPC steps
        p_start = pos;
        for s = 1:params.n_inner
            tau = (s - 1) * dt_inner;
            time_now = (t - 1) * dt_outer + tau;
            p_ref = p_start + tau * v;
            d = wind_fn(time_now, pos);
            switch params.drift_mode
                case 'measured',  d_mpc = d;
                case 'estimated', d_mpc = d_hat;
                case 'none',      d_mpc = zeros(N, 2);
                otherwise, error('main_nodes:driftMode', ...
                    'drift_mode must be ''measured'', ''estimated'' or ''none''.');
            end
            [u, ~] = low_level_mpc(pos, vel, p_ref, v, d_mpc, mpc);
            pos_pred = pos + dt_inner * vel + 0.5 * dt_inner^2 * u;
            pos = pos + dt_inner * (vel + d) + 0.5 * dt_inner^2 * u;
            vel = vel + dt_inner * u;
            if strcmp(params.drift_mode, 'estimated')
                d_hat = (1 - params.drift_alpha) * d_hat + ...
                    params.drift_alpha * (pos - pos_pred) / dt_inner;
            end
        end
    else
        d = wind_fn((t - 1) * dt_outer, pos);
        pos = pos + dt_outer * (v + d);
    end
    pos_hist(:, :, t) = pos;

    if mod(t, 10) == 0 || t == T
        fprintf('Step %d/%d  L1 = %.3f  max|v| = %.2f  sat = %.2f  collapsed = %.2f\n', ...
            t, T, diagn.err_L1(t), diagn.max_speed(t), diagn.frac_sat(t), diagn.collapsed(t));
    end
end

%% Diagnostics
figure; tiledlayout(3, 2);
nexttile; plot(diagn.err_L1);    ylabel('$L^1$ error');
nexttile; plot(diagn.frac_sat);  ylabel('saturated fraction');
nexttile; plot(diagn.max_speed); yline(1, '--'); ylabel('max speed');
nexttile; semilogy(max(diagn.lam_max, eps)); ylabel('max $\lambda$');
nexttile; plot(diagn.mass_T);    ylabel('terminal mass');
nexttile; plot(diagn.collapsed); ylabel('frac. with nn $<0.05$'); xlabel('step');

%% Video of the run
plotresults_w(pos_hist, phi_hist, target_hist, params);
