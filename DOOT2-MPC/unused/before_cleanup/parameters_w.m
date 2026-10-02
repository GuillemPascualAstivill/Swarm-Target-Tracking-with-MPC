%% parameters_w.m  --  WENDLAND version of parameters2.m
% Run this first, then main_w.m.
%
% KEY DIFFERENCE FROM THE GAUSSIAN VERSION:
%   Wendland support = bandwidth (HARD cutoff at r=h). The Gaussian version
%   had sigma=h but reached ~20*sigma before truncation, so its EFFECTIVE
%   reach was far larger than h. To keep the same neighbour connectivity you
%   must set the Wendland h to roughly the Gaussian's EFFECTIVE reach, i.e.
%   a few times the old sigma -- NOT the old sigma. Rule of thumb: Wendland h
%   ~ 2.5-3x the Gaussian sigma gives comparable smoothing/neighbours.
%   Starting values below are scaled up accordingly. If you see stranded
%   agents (isolated, zero velocity) or a jagged field, INCREASE h.

clc; close all; clear;
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelinterpreter','latex')
set(groot,'defaultLegendinterpreter','latex')

rng(10,'twister');

%% Main Parameters
N = 500;
S = 1000;
T = 500;
dt_cont = 0.85;
neigh_floor = 100;

% ------------------------------------------------------------------------
% WENDLAND bandwidths (= support radii, HARD). Scaled up vs the Gaussian
% sigmas (Gaussian used h_agent=1.0, h_field=1.2; here ~2.5x).
% ------------------------------------------------------------------------
h_agent = 2;          % DENSITY support  (was Gaussian sigma 1.0)
h_field = 2;          % FIELD support    (was Gaussian sigma 1.2)
h_ant   = 3.8;          % ANTENNA support  (was Gaussian sigma 1.5)

% Wendland needs NO truncation constants; keep var_* only for the plotter
% fallback / any legacy reference. They are NOT used by the Wendland kernels.
var_agent = h_agent^2;  C_agent = 2*pi*var_agent;  C_agent_r = C_agent;
var_field = h_field^2;  C_field = 2*pi*var_field;  C_field_r = C_field;
var_ant   = h_ant^2;    C_ant   = 2*pi*var_ant;    C_ant_r   = C_ant;

% Truncation radii are irrelevant for Wendland (support is h). Keep fields
% present for compatibility; set them to h so any leftover mask uses h.
R_agent = h_agent;
R_field = h_field;
R_ant   = h_ant;

% Primal-dual
pd_iters =500;
etap = 0.05;
etad = 0.25;
epsilon_reg = 0.15;
w_gain = 1;
delta_ridge = 4e-3; 
terminal = 'kl';       % primal_dual_2: 'l2_relaxed', 'l2', or 'kl'
kappa = 0.2;            % used by the KL terminal mode
kl_floor = 1e-6;        % used by the KL terminal mode

% Antennas (0 -> off)
num_antennas_axis = 0;
antw = 3;
if num_antennas_axis > 0
    [X_ant, Y_ant] = meshgrid(linspace(-8, 8, num_antennas_axis), ...
                              linspace(-8, 8, num_antennas_axis));
    ant_cords = [X_ant(:), Y_ant(:)];
else
    X_ant = []; Y_ant = []; ant_cords = zeros(0,2);
end
M = size(ant_cords, 1);

% Obstacle
obs_center = [1, 1]; obs_radius = 1.5; c_min = 1; c_max = 100; delta = 0.1;

% Ghosts (disabled)
ghost_pos = []; N_ghosts = 0; N_tot = N + N_ghosts; ghost_w = 1;
num_ghosts_ring = 5; theta_ghosts = []; ring_ghost_pos = [];

it_cons = 10; gamma_max_multiplier = 1.0;

params = struct();
params.N=N; params.S=S; params.T=T; params.dt_cont=dt_cont;
params.neigh_floor=neigh_floor;
params.h_agent=h_agent; params.h_field=h_field; params.h_ant=h_ant;
params.R_agent=R_agent; params.R_field=R_field; params.R_ant=R_ant;
params.var_agent=var_agent; params.C_agent=C_agent; params.C_agent_r=C_agent_r;
params.var_field=var_field; params.C_field=C_field; params.C_field_r=C_field_r;
params.var_ant=var_ant; params.C_ant=C_ant; params.C_ant_r=C_ant_r;
params.pd_iters=pd_iters; params.etap=etap; params.etad=etad;
params.epsilon_reg=epsilon_reg;
params.delta_ridge=delta_ridge;
params.terminal=terminal; params.kappa=kappa; params.kl_floor=kl_floor;
params.num_antennas_axis=num_antennas_axis; params.antw=antw;
params.X_ant=X_ant; params.Y_ant=Y_ant; params.ant_cords=ant_cords; params.M=M;
params.obs_center=obs_center; params.obs_radius=obs_radius;
params.c_min=c_min; params.c_max=c_max; params.delta=delta;
params.num_ghosts_ring=num_ghosts_ring; params.theta_ghosts=theta_ghosts;
params.ring_ghost_pos=ring_ghost_pos; params.ghost_pos=ghost_pos;
params.N_ghosts=N_ghosts; params.N_tot=N_tot; params.ghost_w=ghost_w;
params.it_cons=it_cons; params.gamma_max_multiplier=gamma_max_multiplier;
params.w_gain=w_gain;
params.no_smoothing = true;

%grid
grid_res = 100;
x = linspace(-8, 8, grid_res); y = linspace(-8, 8, grid_res);
dx = x(2)-x(1); dy = y(2)-y(1);
[X, Y] = meshgrid(x, y);
grid_pts = [X(:), Y(:)];
params.grid_res=grid_res; params.x=x; params.y=y; params.dx=dx; params.dy=dy;
params.X=X; params.Y=Y; params.grid_pts=grid_pts;

%% Target (star)
% Geometry lives in params so star_target_w.m can rebuild the star at any
% rotation phase; main_w.m does that every step when target_rot is on.
r0 = 3.5; A = 1.25; m = 5; sigma = 1.5;
params.star_r0=r0; params.star_A=A; params.star_m=m; params.star_sigma=sigma;

% polar grids are phase-independent -- build once, reuse every step
params.grid_theta = atan2(Y, X);
params.grid_R     = hypot(X, Y);

% ------------------------------------------------------------------------
% MOVING TARGET: uniform rigid rotation of the star.
%   target_revs = full revolutions over the whole run (T steps).
%   The star has m-fold symmetry, so one revolution passes through the same
%   shape m times; target_revs = 1/m rotates by exactly one petal.
%   Set target_rot = false for the original static target.
% ------------------------------------------------------------------------
target_rot   = true;
target_revs  = 1;
target_omega = 2*pi*target_revs / T;      % rad per step (set directly if you prefer)
params.target_rot   = target_rot;
params.target_omega = target_omega;

[rho_target, rho_target_mat] = star_target_w(params, 0);   % phase 0 at t = 1
params.rho_target = rho_target;
params.rho_target_mat = rho_target_mat;
params.target_sh = rho_target / sum(rho_target);

%% Precompute matched-smoothed target (Wendland), used only if no_smoothing=false
params.rho_target_sm = smooth_on_grid_wendland(rho_target, params, h_agent);

%% Initial swarm
mean_init = [0, 0]; var_init = [3, 3];
gm_init = gmdistribution(mean_init, var_init);
pos_init = random(gm_init, N);
pos_cont = [pos_init; ghost_pos];
params.pos_init = pos_init;
params.pos_cont0 = pos_cont;

fprintf('parameters_w.m (WENDLAND) loaded. h_agent=%.2f h_field=%.2f h_ant=%.2f\n', ...
        h_agent, h_field, h_ant);
fprintf('Wendland support = bandwidth (hard). If agents strand or field is jagged, raise h.\n');
fprintf('Now run main_w.m\n');
