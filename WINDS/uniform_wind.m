function d = uniform_wind(t, x, U, theta)
% UNIFORM_WIND  Spatially constant drift field (the wind_fn contract).
%
%   d = UNIFORM_WIND(t, x, U, theta) returns the world-frame drift
%
%       d = U * [cos(theta), sin(theta)]
%
%   at every point. Bind it for LOW_LEVEL_MPC with
%
%       wind_fn = @(t, x) uniform_wind(t, x, 0.10, 0);
%
%   d enters the plant as xdot = v + d, so it has the units of a VELOCITY --
%   size it against |v_nom| (~0.3 in doot_cbf_mpc.m), not against u.
%
%   PURPOSE: this is the NULL TEST for the -d feedforward in low_level_mpc.m.
%   grad d = 0, so v_ref = v_tgt - d makes zeta = 0 an exact equilibrium of
%   the plant and the standing position offset must come out numerically zero
%   FOR ANY Q_u. If it does not, the bug is in the controller, not the field,
%   and there is no point running taylor_green_wind yet.
%
%   For the same reason it cannot stress the MPC: no Q_u dependence exists to
%   measure. It is still a fair cross-branch experiment once the DOOT and
%   B-CBF branches are advected too -- they carry no feedforward, so they eat
%   the drift as a persistent bias while the MPC branch cancels it.
%
%   SIZING: keep U <~ 0.15. DOOT is closed-loop in position (v_nom is replanned
%   from the actual pos_nom every outer step), so it tolerates a bias while
%   |d| < |v_nom| ~ 0.3. Above that the swarm is blown off the target support
%   and the metrics measure a runaway, not a bias.
%
%   Inputs
%     t     : time (unused; present for the wind_fn signature)
%     x     : position, 1x2 row / 2x1 column / Mx2 stack
%     U     : drift speed                              (default 0.10)
%     theta : drift direction, radians                 (default 0)
%
%   Output
%     d     : Mx2 drift, one row per input position

    if nargin < 3 || isempty(U),     U     = 0.10; end
    if nargin < 4 || isempty(theta), theta = 0;    end

    if isvector(x), x = x(:).'; end
    d = U * repmat([cos(theta), sin(theta)], size(x,1), 1);
end
