function [K, Dfac] = wendland_kernel(sq_dists, h)
%WENDLAND_KERNEL  Wendland C2 RBF value and gradient factor from squared distances.
%
%   [K, Dfac] = wendland_kernel(sq_dists, h)
%
%   Wendland C2 in 2D, compact support radius h, normalised so int_{R^2} K = 1:
%       K_h(x)      = (7/(pi h^2)) (1-r)_+^4 (4r+1),   r = |x|/h
%       grad K_h(x) = Dfac(r) * x,   Dfac(r) = -(140/(pi h^4)) (1-r)_+^3
%
%   So for a difference matrix diffs (.,.,1)=dx, (.,.,2)=dy:
%       D_x = Dfac .* dx,   D_y = Dfac .* dy
%   which mirrors the Gaussian code's  D_x = -(dx/var).*K  but with the
%   Wendland factor instead of -K/var.
%
%   INPUTS
%     sq_dists : matrix of SQUARED distances |x_i - x_j|^2 (any shape)
%     h        : bandwidth = COMPACT SUPPORT RADIUS (hard cutoff at r=h)
%   OUTPUTS
%     K    : kernel values (same shape); exactly 0 where |x| >= h
%     Dfac : scalar gradient factor (same shape); 0 where |x| >= h.
%            grad K(x) = Dfac .* x componentwise.
%
%   NOTE. Unlike the truncated Gaussian, Wendland needs NO renormalisation
%   constant (C_agent_r): the (7/pi h^2) prefactor already gives unit mass.
%   Support is EXACTLY h, so pick h ~ the Gaussian's effective reach (a few
%   sigma), NOT the Gaussian's sigma, or agents lose neighbours.

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
