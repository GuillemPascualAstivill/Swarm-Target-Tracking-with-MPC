%% parameters_w.m  --  parameters for main_w_final.m (run this first)
% Wendland kernels have compact support = bandwidth h. If agents get stranded
% (isolated, zero velocity) or the potential field looks jagged, increase h.

clc; close all; clear;
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelinterpreter','latex')
set(groot,'defaultLegendinterpreter','latex')

rng(10,'twister');

%% Simulation
params.N = 500;             % number of agents
params.T = 500;             % number of DOOT steps
params.dt_cont = 0.85;      % DOOT time step

%% Low-level MPC and wind
params.use_mpc = true;     % false: DOOT velocity directly; true: use low-level MPC
params.n_inner = 10;        % MPC steps per DOOT step
params.Q_x = 1;             % MPC position weight
params.Q_u = 3.5;           % MPC control weight
params.Nh = 10;             % MPC horizon
params.gamma = 1;            % MPC terminal-cost scale
params.u_max = [];           % acceleration bound; [] means unconstrained
params.wind_model = 'taylor_green';  % 'none', 'uniform' or 'taylor_green'
params.wind_U = 0.1;        % wind speed
params.wind_theta = 0;       % uniform wind direction [rad]
params.wind_L = 6;           % Taylor-Green vortex period
params.drift_mode = 'estimated'; % 'measured', 'estimated' or 'none'
params.drift_alpha = 0.2;    % estimated-drift filter gain in (0, 1]

%% Kernel bandwidths (Wendland support radii)
params.h_agent = 2;         % agent density rho_N
params.h_field = 2;         % potential field phi

%% Primal-dual solver (primal_dual.m)
params.pd_iters = 500;      % iterations per DOOT step
params.etap = 0.05;         % primal step size
params.etad = 0.25;         % dual step size
params.delta_ridge = 4e-3;  % ridge (delta/2)*||alpha||^2

%% Terminal cost
params.terminal = 'l2';     % 'kl' or 'l2'
params.kappa = 0.2;         % 'kl': weight of the KL terminal cost
params.kl_floor = 1e-6;     % 'kl': floor on the target, relative to its max
params.epsilon_reg = 0.15;  % 'l2': eps in the terminal cost ||rho - rho*||^2 / (2 eps)

%% Grid on [-8, 8]^2
grid_res = 100;
x = linspace(-8, 8, grid_res);
y = linspace(-8, 8, grid_res);
[X, Y] = meshgrid(x, y);
params.grid_res = grid_res;
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
params.target_floor = 0; % additive uniform background before normalization
target_revs = 1;            % revolutions over the whole run (0 = static target)
params.target_omega = 2*pi*target_revs / params.T;   % rotation per step [rad]

%% Initial swarm: Gaussian, mean [0 0], variance 3 per axis
gm_init = gmdistribution([0, 0], [6, 6]);
params.pos0 = random(gm_init, params.N);
