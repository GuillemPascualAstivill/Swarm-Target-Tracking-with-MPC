function rho = star_target_w(params, phase)
%STAR_TARGET_W  Star-shaped target density on the grid, rotated by PHASE [rad].
%       R_star(theta) = r0 + A cos( m (theta - phase) )
%       rho(x)        = exp( -(|x| - R_star(angle(x)))^2 / (2 sigma^2) )
%   normalised to unit mass on the grid, returned as a column vector.

theta = atan2(params.Y, params.X);
R     = hypot(params.X, params.Y);

R_star = params.star_r0 + params.star_A * cos(params.star_m * (theta - phase));
rho = exp(-(R - R_star).^2 / (2 * params.star_sigma^2));
if isfield(params, 'target_floor') && ~isempty(params.target_floor)
	rho = rho + params.target_floor;
end
rho = rho(:);
rho = rho / (sum(rho) * params.dx * params.dy);
end
