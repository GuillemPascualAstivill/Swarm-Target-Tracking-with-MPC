function c_sq = conformal(pos_cont, pos_ant, params)
%conformal Compute the obstacle eikonal factor squared for all coordinates.

obs_center = params.obs_center;
obs_radius = params.obs_radius;
c_min = params.c_min;
c_max = params.c_max;
delta = params.delta;

all_coords = [pos_cont; pos_ant];
r = sqrt(sum((all_coords - obs_center).^2, 2));
d = r - obs_radius;
c_factor = c_min + (c_max - c_min) ./ (1 + exp(d / delta));
c_sq = c_factor.^2;

end