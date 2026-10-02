function nb = build_nodes(params)
%BUILD_NODES  Fixed node basis for the potential phi (build ONCE, before the loop).
%
%   nb = build_nodes(params) places a regular grid of num_nodes_axis^2 nodes
%   over the domain and precomputes everything that does not depend on the
%   agents:
%     nb.A         M x 2 node positions
%     nb.h         node kernel support (params.h_node)
%     nb.Kgn       G x M sparse, Kgn(g,a) = K_node(grid_g - A_a)
%     nb.KgnT      its transpose
%     nb.Dnx/Dny   M x M sparse, gradient of every node kernel at every node
%                  (rows of the Lipschitz constraint at the nodes)
%
%   phi(x) = sum_a beta_a K_node(x - A_a) then lives on the whole domain,
%   independently of where the agents are.

n = params.num_nodes_axis;
h = params.h_node;
a = linspace(params.x(1), params.x(end), n);
[AX, AY] = meshgrid(a, a);
A = [AX(:), AY(:)];
M = size(A, 1);

G  = params.grid_pts;
nG = size(G, 1);
I = []; J = []; V = [];
chunk = 2000;                                  % limits memory to chunk x M
for s = 1:chunk:nG
    idx = (s:min(s+chunk-1, nG)).';
    K = wendland_kernel(pdist2(G(idx,:), A).^2, h);
    [i, j, v] = find(K);
    I = [I; idx(i)]; J = [J; j]; V = [V; v];  %#ok<AGROW>
end
nb.Kgn  = sparse(I, J, V, nG, M);
nb.KgnT = nb.Kgn.';

[~, Dfac] = wendland_kernel(pdist2(A, A).^2, h);
nb.Dnx = sparse(Dfac .* (A(:,1) - A(:,1).'));
nb.Dny = sparse(Dfac .* (A(:,2) - A(:,2).'));

nb.A = A;
nb.h = h;
end
