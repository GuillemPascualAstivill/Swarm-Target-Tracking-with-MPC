function integral = computeTargetint_w(pos_cont, pos_ant, params)
%computeTargetint_w  WENDLAND version of computeTargetint2.
%   Same decoupled-bandwidth structure and matched target smoothing, but all
%   grid kernels are Wendland C2 (compact support = bandwidth), no C_*_r.
%
%   Mismatch:  w_j = int K_field(x-x_j) ( rho_N(x) - rho*_sm(x) ) dx
%     rho_N built with h_agent (density), projected with h_field (field basis).
%     target pre-smoothed by h_agent (matched filter) unless no_smoothing=true.

N = params.N;
N_tot = params.N_tot;
grid_pts = params.grid_pts;
rho_target = params.rho_target;
dx = params.dx;
dy = params.dy;

h_agent = params.h_agent;
if isfield(params,'h_field') && ~isempty(params.h_field), h_field = params.h_field; else, h_field = h_agent; end
if isfield(params,'h_ant')   && ~isempty(params.h_ant),   h_ant   = params.h_ant;   else, h_ant   = h_agent; end

% -------------------------------------------------------------------------
% DENSITY kernel on grid (h_agent) -> rho_N
% -------------------------------------------------------------------------
diffs_grid_agent = permute(grid_pts, [1 3 2]) - permute(pos_cont, [3 1 2]);
sq_dists_grid_agent = sum(diffs_grid_agent.^2, 3);
K_grid_agent = wendland_kernel(sq_dists_grid_agent, h_agent);
integral.K_grid_agent = K_grid_agent;

% -------------------------------------------------------------------------
% FIELD kernel on grid (h_field) -> basis projection
% -------------------------------------------------------------------------
K_grid_field = wendland_kernel(sq_dists_grid_agent, h_field);
integral.K_grid_field = K_grid_field;

% -------------------------------------------------------------------------
% MATCHED-SMOOTHED target (raw if no_smoothing)
% -------------------------------------------------------------------------
if isfield(params,'no_smoothing') && params.no_smoothing
    rho_target_sm = rho_target;
elseif isfield(params, 'rho_target_sm') && ~isempty(params.rho_target_sm)
    rho_target_sm = params.rho_target_sm;
else
    rho_target_sm = smooth_on_grid_wendland(rho_target, params, h_agent);
end
integral.rho_target_sm = rho_target_sm;

% Target density at the current agent locations.  The KL/L2 terminal modes
% use these values directly, avoiding a repeated Gram-system reconstruction.
rho_target_sm_mat = reshape(rho_target_sm, size(params.X));
integral.rho_target_at_agents = interp2(params.X, params.Y, ...
    rho_target_sm_mat, pos_cont(:,1), pos_cont(:,2), 'linear', 0);

% Target projected onto FIELD basis
integral.target_ag = K_grid_field' * (rho_target_sm * dx * dy);
integral.target_ag(N+1:N_tot) = 0;

% Density projected onto FIELD basis (cross term):
%   sum_i agents_ag(j,i)/N = int K_field(x-x_j) rho_N dx
integral.agents_ag = (K_grid_field' * K_grid_agent) * dx * dy;

% ghost agents (infinite tail). Wendland has no analytic infinite tail; use a
% Wendland self-convolution proxy at sqrt(2)*h_field support. Ghosts disabled
% in current runs, so this only matters if you enable them.
if N_tot > N
    hh = sqrt(2)*h_field;
    sq_dists_aa = sum((permute(pos_cont, [1 3 2]) - permute(pos_cont, [3 1 2])).^2, 3);
    K_tail = wendland_kernel(sq_dists_aa, hh);
    integral.agents_ag(N+1:end, 1:N)     = K_tail(N+1:end, 1:N);
    integral.agents_ag(1:N, N+1:end)     = K_tail(1:N, N+1:end);
    integral.agents_ag(N+1:end, N+1:end) = K_tail(N+1:end, N+1:end);
end

% FIELD Gram for the regulariser
integral.field_gram = (K_grid_field' * K_grid_field) * dx * dy;

% -------------------------------------------------------------------------
% Regularisation matrix (K_reg)
% -------------------------------------------------------------------------
if ~isempty(pos_ant)
    diffs_grid_ant = permute(grid_pts, [1 3 2]) - permute(pos_ant, [3 1 2]);
    sq_dists_grid_ant = sum(diffs_grid_ant.^2, 3);
    K_grid_ant = wendland_kernel(sq_dists_grid_ant, h_ant);

    integral.target_ant = K_grid_ant' * (rho_target_sm * dx * dy);
    integral.agents_ant = (K_grid_field' * K_grid_ant) * dx * dy;
    integral.ant_ant    = (K_grid_ant'  * K_grid_ant) * dx * dy;

    base_K_reg = [
        integral.field_gram,   integral.agents_ant;
        integral.agents_ant',  integral.ant_ant
    ];
    delta = 1;
    integral.K_reg = base_K_reg + delta * eye(size(base_K_reg,1));
else
    integral.target_ant = [];
    integral.agents_ant = [];
    integral.ant_ant    = [];
    delta = 0;
    integral.K_reg = integral.field_gram + delta * eye(N_tot);
end

% -------------------------------------------------------------------------
% Target integral (Voronoi) -- density assignment, raw target
% -------------------------------------------------------------------------
[~, closest_agent_idx] = min(sq_dists_grid_agent, [], 2);
N_grid = size(sq_dists_grid_agent, 1);
N_agents_total = size(sq_dists_grid_agent, 2);
is_closest_mask = false(N_grid, N_agents_total);
linear_indices = sub2ind([N_grid, N_agents_total], (1:N_grid)', closest_agent_idx);
is_closest_mask(linear_indices) = true;
mask_grid_agent = sq_dists_grid_agent <= h_agent^2;
mask_grid_agent_exclusive = mask_grid_agent & is_closest_mask;
integral.target_vor = double(mask_grid_agent_exclusive)' * (rho_target * dx * dy);
integral.target_vor(N+1:N_tot) = 0;

end


% =====================================================================
%  Wendland grid smoothing (matched filter), normalised to integrate to 1
% =====================================================================
function f_sm = smooth_on_grid_wendland(f, params, h)
    n  = params.grid_res;
    dg = params.dx;                       % assumes dx == dy
    F  = reshape(f, n, n);
    rad = ceil(h / dg);                   % Wendland support is exactly h
    [ix, iy] = meshgrid(-rad:rad, -rad:rad);
    d2 = (ix*dg).^2 + (iy*dg).^2;
    Kst = wendland_kernel(d2, h);
    Kst = Kst / (sum(Kst(:)) * dg^2);     % normalise
    F_sm = conv2(F, Kst, 'same') * dg^2;
    f_sm = F_sm(:);
end
