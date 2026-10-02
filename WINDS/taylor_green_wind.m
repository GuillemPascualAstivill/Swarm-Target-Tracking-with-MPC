function d = taylor_green_wind(t, x, U, L, x0)
% TAYLOR_GREEN_WIND  Taylor-Green vortex drift field (the wind_fn contract).
%
%   d = TAYLOR_GREEN_WIND(t, x, U, L, x0) returns the world-frame drift
%
%       d = U * [  sin(k*X) .* cos(k*Y) ,
%                 -cos(k*X) .* sin(k*Y) ],     k = 2*pi/L,  [X Y] = x - x0
%
%   Bind it for LOW_LEVEL_MPC with
%
%       wind_fn = @(t, x) taylor_green_wind(t, x, 0.10, 6);
%
%   d enters the plant as xdot = v + d, so it has the units of a VELOCITY --
%   size it against |v_nom| (~0.3 in doot_cbf_mpc.m), not against u.
%
%   WHY THIS FIELD
%     - Bounded: |d| <= U everywhere, so no runaway advection.
%     - DIVERGENCE FREE: div d = U*k*cos(kX)*cos(kY) - U*k*cos(kX)*cos(kY) = 0.
%       This is the decisive property here. The B-CBF branch and every
%       epsilon / epsilon_min metric read a KDE density; a compressible field
%       would move that density on its own, and "wind compressed the swarm"
%       could not be separated from "controller failed".
%     - Zero mean over a cell, so it perturbs transport without biasing the
%       swarm off the target support.
%     - Analytic gradient, hence the closed-form error prediction below.
%
%   Cell corners (multiples of L/2 from x0) are stagnation points, cell
%   centers are vortex cores. At L = 6 the agents cross ~3 cells per axis over
%   the [-9,9] working domain.
%
%   WHAT IT STRESSES (this is TODO [5] in doot_cbf_mpc.m)
%   The -d feedforward cancels a CONSTANT drift exactly, so only advection of
%   the spatial variation survives. Per axis, with e = p - x_nom under
%   u = -k_p*e - k_v*(v - v_ref):
%
%       edotdot + k_v*edot + k_p*e = (grad d) * xdot
%
%   With |grad d| <= U*k (tight for Taylor-Green) and |xdot| ~ |v_nom|:
%
%       |e_inf| ~ U*k*|v_nom| / k_p ,        k_p = sqrt(Q_x/Q_u)
%
%   which GROWS with the smoothness knob Q_u -- the opposite of the constant-
%   drift case, where the offset is zero for any tuning. At Q_u = 35
%   (k_p = 0.169) and |v_nom| ~ 0.3:
%
%       U = 0.10, L = 6  ->  |e_inf| ~ 0.19   (< bw = 0.3, safe default)
%       U = 0.15, L = 4  ->  |e_inf| ~ 0.42   (> bw, visible density damage)
%
%   Sweep U up from the default until the W2 metric degrades; that crossing is
%   the reportable number, and the estimate says where to expect it.
%
%   Inputs
%     t  : time (unused -- steady field; wrap the call to modulate, e.g.
%          @(t,x) min(t/20,1) * taylor_green_wind(t,x,U,L) to ramp the wind in
%          rather than hitting agents that start at rest, or
%          @(t,x) exp(-2*nu*(2*pi/L)^2*t) * taylor_green_wind(...) for the
%          exact decaying Navier-Stokes solution)
%     x  : position, 1x2 row / 2x1 column / Mx2 stack
%     U  : peak drift speed                            (default 0.10)
%     L  : spatial period of the vortex lattice        (default 6.0)
%     x0 : lattice origin offset, 1x2                  (default [0 0])
%
%   Output
%     d  : Mx2 drift, one row per input position

    if nargin < 3 || isempty(U),  U  = 0.10; end
    if nargin < 4 || isempty(L),  L  = 6.0;  end
    if nargin < 5 || isempty(x0), x0 = [0 0]; end

    if isvector(x), x = x(:).'; end
    k  = 2*pi / L;
    X  = k * (x(:,1) - x0(1));
    Y  = k * (x(:,2) - x0(2));
    d  = U * [ sin(X).*cos(Y), -cos(X).*sin(Y) ];
end
