function f_sm = smooth_on_grid_wendland(f, params, h)
%SMOOTH_ON_GRID_WENDLAND  Wendland C2 grid smoothing (matched filter).
%   Standalone copy so parameters_w.m can precompute the smoothed target.
%   Convolves grid field f with a Wendland C2 kernel of support radius h,
%   normalised to integrate to 1.

    n  = params.grid_res;
    dg = params.dx;                       % assumes dx == dy
    F  = reshape(f, n, n);
    rad = ceil(h / dg);
    [ix, iy] = meshgrid(-rad:rad, -rad:rad);
    d2 = (ix*dg).^2 + (iy*dg).^2;
    Kst = wendland_kernel(d2, h);
    Kst = Kst / (sum(Kst(:)) * dg^2);
    F_sm = conv2(F, Kst, 'same') * dg^2;
    f_sm = F_sm(:);
end
