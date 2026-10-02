%% parameters_w.m  --  parameters for main_nodes.m (run this first)
% Wendland kernels have compact support = bandwidth h. If agents get stranded
% (isolated, zero velocity) or the density looks jagged, increase h.

clc; close all; clear;
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelinterpreter','latex')
set(groot,'defaultLegendinterpreter','latex')

rng(10,'twister');

%% Simulation
params.N = 500;             % number of agents
params.T = 1000;            % number of DOOT steps
params.dt_cont = 0.5;       % DOOT time step ('kl': keep dt_cont <= 0.5/kappa, else it may be unstable)

%% Low-level MPC and wind
params.use_mpc = false;     % false: DOOT velocity directly; true: use low-level MPC
params.n_inner = 10;        % MPC steps per DOOT step
params.Q_x = 1;             % MPC position weight
params.Q_u = 3.5e-3;        % MPC control weight
params.Nh = 10;             % MPC horizon
params.gamma = 1;            % MPC terminal-cost scale
params.u_max = [];           % acceleration bound; [] means unconstrained
params.wind_model = 'none';  % 'none', 'uniform' or 'taylor_green'
params.wind_U = 0.01;        % wind speed
params.wind_theta = 0;       % uniform wind direction [rad]
params.wind_L = 6;            % Taylor-Green vortex period
params.drift_mode = 'estimated'; % 'measured', 'estimated' or 'none'
params.drift_alpha = 0.2;    % estimated-drift filter gain in (0, 1]

%% Kernels (Wendland support radii)
params.h_agent = 2;       % agent density rho_N
params.num_nodes_axis = 33; % potential phi: 33 x 33 fixed nodes, spacing 16/32 = 0.5
params.h_node = 1.5;        % node kernel support, about 3x the node spacing

%% Primal-dual solver (plan_fixed.m)
params.pd_iters = 500;      % iterations per DOOT step
params.etap = 0.05;         % primal step size
params.etad = 0.25;         % dual step size
params.delta_ridge = 1e-3;  % ridge (delta/2)*||beta||^2

%% Terminal cost F (regularization toward the target)
params.terminal = 'kl';     % 'kl':         F = kappa * KL(rho | rho*)
                            % 'l2':         F = ||rho - rho*||^2 / (2 eps), rho kept >= 0 (clipped)
                            % 'l2_relaxed': F = ||rho - rho*||^2 / (2 eps), rho may go < 0
params.kappa = 1;           % 'kl': weight (larger = closer to the target)
params.kl_floor = 1e-6;     % 'kl': floor on rho*, relative to max(rho*)
params.epsilon_reg = 0.0025; % 'l2', 'l2_relaxed': eps (smaller = closer to the target)

%% Grid on [-8, 8]^2
grid_res = 100;
x = linspace(-8, 8, grid_res);
y = linspace(-8, 8, grid_res);
[X, Y] = meshgrid(x, y);
params.grid_res = grid_res;
params.x = x;
params.dx = x(2) - x(1);
params.dy = y(2) - y(1);
params.X = X;
params.Y = Y;
params.grid_pts = [X(:), Y(:)];

%% Target: rotating star, radius R(theta) = r0 + A cos(m theta)
params.star_r0 = 3.5;
params.star_A = 1.25;
params.star_m = 5;
params.star_sigma = 1.5;    % radial width
target_revs = 1;            % revolutions over the whole run (0 = static target)
params.target_omega = 2*pi*target_revs / params.T;   % rotation per step [rad]

%% Initial swarm: Gaussian, mean [0 0], variance 3 per axis
gm_init = gmdistribution([0, 0], [3, 3]);
params.pos0 = random(gm_init, params.N);
