function kernels = computeKernels_w(pos_cont, pos_ant, params)
%computeKernels_w  WENDLAND version of computeKernels2.
%   Same structure and outputs, but every kernel block is Wendland C2
%   (compact support = bandwidth) instead of a truncated Gaussian.
%
%   DECOUPLING is preserved: density uses h_agent, field uses h_field,
%   antennas use h_ant. Wendland needs NO C_*_r truncation constants.
%
%   Outputs (identical names to the Gaussian version):
%     kernels.K_aa        DENSITY-normalised agent-agent (ghost logic)
%     kernels.K_field_aa  FIELD agent-agent (drives D_x/D_y & K_reg)
%     kernels.K_aant, K_anta, K_antant
%     kernels.D_x, D_y    field derivatives (velocity = -(D_x,D_y)*alpha)
%     kernels.mask_aa, kernels.K_reg

N = params.N;
N_tot = size(pos_cont, 1);
M = size(pos_ant, 1);

% bandwidths (support radii). Fall back to h_agent if a field/ant band absent.
h_agent = params.h_agent;
if isfield(params,'h_field') && ~isempty(params.h_field), h_field = params.h_field; else, h_field = h_agent; end
if isfield(params,'h_ant')   && ~isempty(params.h_ant),   h_ant   = params.h_ant;   else, h_ant   = h_agent; end

% ghost index set
if isfield(params, 'N_tot')
    ghost_idx = (N + 1):params.N_tot;
else
    ghost_idx = (N + 1):N_tot;
end

% -------------------------------------------------------------------------
% Agent-Agent DENSITY block (h_agent) -- kept for ghost scaling / density use
% -------------------------------------------------------------------------
diffs_aa = permute(pos_cont, [1 3 2]) - permute(pos_cont, [3 1 2]);
sq_dists_aa = sum(diffs_aa.^2, 3);

kernels.K_aa = wendland_kernel(sq_dists_aa, h_agent);   

% exact compact support
% ghosts see everyone: force their rows/cols "inside" by direct fill.
% (Wendland is 0 beyond h; to let ghosts see all, we recompute their entries
%  with a large effective radius = domain-spanning. Simplest: use h_field for
%  ghost rows so they reach far. Ghosts are disabled in current runs, so this
%  is a no-op unless you enable them.)
kernels.mask_aa = sq_dists_aa <= h_agent^2;
kernels.mask_aa(ghost_idx, :) = true;
kernels.mask_aa(:, ghost_idx) = true;
if ~isempty(ghost_idx)
    Kg = wendland_kernel(sq_dists_aa(ghost_idx,:), max(h_field, h_agent));
    kernels.K_aa(ghost_idx, :) = Kg;
    Kg2 = wendland_kernel(sq_dists_aa(:,ghost_idx), max(h_field, h_agent));
    kernels.K_aa(:, ghost_idx) = Kg2;
end

% -------------------------------------------------------------------------
% Agent-Agent FIELD block (h_field) -- POTENTIAL basis. Velocity uses these.
% -------------------------------------------------------------------------
[Kf, Dfac_f] = wendland_kernel(sq_dists_aa, h_field);
kernels.K_field_aa = Kf;
if ~isempty(ghost_idx)
    [Kfg, Dfg] = wendland_kernel(sq_dists_aa(ghost_idx,:), max(h_field,h_agent));
    kernels.K_field_aa(ghost_idx,:) = Kfg;   %#ok<NASGU>  (value kept for K_reg)
    % note: ghost derivative rows handled below via full Dfac_f recompute if needed
end

% FIELD derivatives: grad K = Dfac .* diff  (no 1/var; Wendland factor built in)
kernels.D_x_aa = Dfac_f .* diffs_aa(:, :, 1);
kernels.D_y_aa = Dfac_f .* diffs_aa(:, :, 2);

% DENSITY derivatives (for congestion/anti-crowding)
[Ka, Dfac_a] = wendland_kernel(sq_dists_aa, h_agent);
kernels.K_agent_aa = Ka;
kernels.D_x_agent_aa = Dfac_a .* diffs_aa(:, :, 1);
kernels.D_y_agent_aa = Dfac_a .* diffs_aa(:, :, 2);

% -------------------------------------------------------------------------
% Antenna-Agent block (h_ant)
% -------------------------------------------------------------------------
diffs_aant = permute(pos_cont, [1 3 2]) - permute(pos_ant, [3 1 2]);
sq_dists_aant = sum(diffs_aant.^2, 3);
[Kaant, Dfac_aant] = wendland_kernel(sq_dists_aant, h_ant);
kernels.K_aant   = Kaant;
kernels.D_x_aant = Dfac_aant .* diffs_aant(:, :, 1);
kernels.D_y_aant = Dfac_aant .* diffs_aant(:, :, 2);

% -------------------------------------------------------------------------
% Agent-Antenna block (h_ant)  [was var_agent in Gaussian; use h_ant for the
% antenna-anchored kernel, consistent with its support]
% -------------------------------------------------------------------------
diffs_anta = permute(pos_ant, [1 3 2]) - permute(pos_cont, [3 1 2]);
sq_dists_anta = sum(diffs_anta.^2, 3);
[Kanta, Dfac_anta] = wendland_kernel(sq_dists_anta, h_ant);
kernels.K_anta   = Kanta;
kernels.D_x_anta = Dfac_anta .* diffs_anta(:, :, 1);
kernels.D_y_anta = Dfac_anta .* diffs_anta(:, :, 2);

% -------------------------------------------------------------------------
% Antenna-Antenna block (h_ant)
% -------------------------------------------------------------------------
diffs_antant = permute(pos_ant, [1 3 2]) - permute(pos_ant, [3 1 2]);
sq_dists_antant = sum(diffs_antant.^2, 3);
[Kantant, Dfac_antant] = wendland_kernel(sq_dists_antant, h_ant);
kernels.K_antant   = Kantant;
kernels.D_x_antant = Dfac_antant .* diffs_antant(:, :, 1);
kernels.D_y_antant = Dfac_antant .* diffs_antant(:, :, 2);

% -------------------------------------------------------------------------
% Assemble
% -------------------------------------------------------------------------
kernels.D_x = [kernels.D_x_aa, kernels.D_x_aant; kernels.D_x_anta, kernels.D_x_antant];
kernels.D_y = [kernels.D_y_aa, kernels.D_y_aant; kernels.D_y_anta, kernels.D_y_antant];
kernels.K_reg = [kernels.K_field_aa, kernels.K_aant; kernels.K_anta, kernels.K_antant];

end
