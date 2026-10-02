function [rho, rho_mat] = star_target_w(params, phase)
%STAR_TARGET_W  Star-shaped target density on the grid, rotated by PHASE.
%
%   rho = star_target_w(params, phase) returns the normalised star density as
%   a column vector, with the star rotated rigidly by PHASE radians:
%
%       R_star(theta) = r0 + A cos( m (theta - phase) )
%       rho(x)        = exp( -(|x| - R_star(angle(x)))^2 / (2 sigma^2) )
%
%   phase = 0 reproduces the static target parameters_w.m builds at t = 1.
%   Because the star has m-fold symmetry, a rotation of 2*pi/m maps it onto
%   itself -- one full revolution passes through the same shape m times.
%
%   Needs params.star_r0/_A/_m/_sigma, params.dx/dy/grid_res, and the grid.
%   Uses the cached polar grids params.grid_theta / params.grid_R when present
%   (they do not depend on phase, so parameters_w.m builds them once).

if isfield(params, 'grid_theta') && ~isempty(params.grid_theta)
    theta = params.grid_theta;
    R     = params.grid_R;
else
    theta = atan2(params.Y, params.X);
    R     = hypot(params.X, params.Y);
end

R_star = params.star_r0 + params.star_A * cos(params.star_m * (theta - phase));
rho = exp(-(R - R_star).^2 / (2 * params.star_sigma^2));
rho = rho(:);
rho = rho / (sum(rho) * params.dx * params.dy);

if nargout > 1
    rho_mat = reshape(rho, params.grid_res, params.grid_res);
end
end
