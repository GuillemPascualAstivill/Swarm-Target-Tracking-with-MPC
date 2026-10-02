function [K, Dfac] = wendland_kernel(sq_dists, h)
%WENDLAND_KERNEL  Wendland C2 kernel in 2D with compact support radius h.
%       K_h(x)      = (7/(pi h^2)) (1-r)_+^4 (4r+1),   r = |x|/h   (integrates to 1)
%       grad K_h(x) = Dfac .* x,   Dfac = -(140/(pi h^4)) (1-r)_+^3
%   sq_dists holds squared distances |x|^2 (any shape); K and Dfac match it.

    r = sqrt(sq_dists) / h;               % normalised radius
    inside = r < 1;
    q = zeros(size(r));
    q(inside) = 1 - r(inside);            % (1-r)_+

    K = zeros(size(r));
    K(inside) = (7 / (pi * h^2)) * (q(inside).^4) .* (4*r(inside) + 1);

    if nargout > 1
        Dfac = zeros(size(r));
        Dfac(inside) = -(140 / (pi * h^4)) * (q(inside).^3);
    end
end
