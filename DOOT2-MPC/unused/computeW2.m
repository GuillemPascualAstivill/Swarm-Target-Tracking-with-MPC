function [W, info] = computeW2(source_density, target_density, x, y, varargin)
%COMPUTEW2 Entropic continuous Wasserstein distance on a Cartesian grid.
%   The densities are interpreted as continuous densities and normalised with
%   dx*dy. The Sinkhorn kernel is applied by separable convolutions, so the
%   full grid-by-grid transport matrix is never formed.
%
%   source_density and target_density may be vectors or matrices matching
%   meshgrid(x,y). Optional name/value pairs:
%     'epsilon'  Entropic regularisation for the selected ground cost.
%                Default: 1 percent of an estimated grid cost.
%     'p'        Ground-cost exponent: 1 for Euclidean W1, 2 for W2.
%     'iters'    Sinkhorn iterations (default 500).
%     'tol'      Marginal tolerance (default 1e-4).

ip = inputParser;
ip.addParameter('epsilon', [], @(v) isempty(v) || (isscalar(v) && v > 0));
ip.addParameter('p', 2, @(v) isscalar(v) && (v == 1 || v == 2));
ip.addParameter('iters', 500, @(v) isscalar(v) && v >= 1);
ip.addParameter('tol', 1e-4, @(v) isscalar(v) && v > 0);
ip.parse(varargin{:});
opt = ip.Results;

x = x(:)';
y = y(:);
if numel(x) < 2 || numel(y) < 2
    error('computeW2:badGrid', 'x and y must each contain at least two points.');
end
dx = x(2) - x(1);
dy = y(2) - y(1);
if any(abs(diff(x) - dx) > 1e-12) || any(abs(diff(y) - dy) > 1e-12)
    error('computeW2:nonuniformGrid', 'x and y must be uniformly spaced.');
end

nx = numel(x);
ny = numel(y);
source = reshape(double(source_density), ny, nx);
target = reshape(double(target_density), ny, nx);
if ~all(isfinite(source(:))) || ~all(isfinite(target(:))) || ...
        any(source(:) < 0) || any(target(:) < 0)
    error('computeW2:badDensity', 'Densities must be finite and nonnegative.');
end

source = source * (dx * dy);
target = target * (dx * dy);
source = source / sum(source(:));
target = target / sum(target(:));

if isempty(opt.epsilon)
    if opt.p == 1
        median_cost = hypot(x(end) - x(1), y(end) - y(1)) / 2;
    else
        median_cost = (x(end) - x(1))^2 / 6 + (y(end) - y(1))^2 / 6;
    end
    epsilon = 0.01 * median_cost;
else
    epsilon = opt.epsilon;
end

% Build the radial ground-cost kernel. W2 can use separable convolutions,
% but W1's Euclidean cost requires the full radial 2-D kernel.
ix = (-(nx - 1):(nx - 1)) * dx;
iy = (-(ny - 1):(ny - 1)) * dy;
[DX, DY] = meshgrid(ix, iy);
if opt.p == 1
    ground_cost = hypot(DX, DY);
else
    ground_cost = DX.^2 + DY.^2;
end
K = exp(-ground_cost / epsilon);
applyK = @(z) conv2(z, K, 'same');

u = ones(ny, nx);
v = ones(ny, nx);
for it = 1:opt.iters
    Kv = applyK(v);
    u = source ./ max(Kv, realmin);
    Ku = applyK(u);
    v = target ./ max(Ku, realmin);
end

row_mass = u .* applyK(v);
col_mass = v .* applyK(u);
marg_err = max(sum(abs(row_mass(:) - source(:))), ...
               sum(abs(col_mass(:) - target(:))));

% Apply the ground-cost kernel to v, then evaluate <C, Pi>.
cost_kernel = ground_cost .* K;
cost_field = conv2(v, cost_kernel, 'same');
transport_cost = sum(u(:) .* cost_field(:));
if opt.p == 1
    W = max(transport_cost, 0);
else
    W = sqrt(max(transport_cost, 0));
end

info = struct('converged', marg_err < opt.tol, 'marg_err', marg_err, ...
              'epsilon', epsilon, 'iters', opt.iters, 'p', opt.p, ...
              'n_grid', nx * ny);
end
