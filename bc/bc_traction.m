function f = bc_traction(f, mesh, edge_nodes, tx, ty)
% Add edge traction loads to force vector f (continuum elements only).
%   f          - [nDof x 1] existing force vector (modified in place)
%   mesh       - mesh struct
%   edge_nodes - [nEdges x 2] node index pairs defining loaded edges
%   tx, ty     - traction components (scalar for uniform)
%                or [nEdges x 1] per-edge values
%
% Uses consistent (linear) load integration along each edge.
% Equivalent nodal forces: f_i = tx * L/2,  f_j = tx * L/2  (uniform case)
if ~ismember(mesh.type, {'T3', 'Q4'})
    error('bc_traction: only for continuum elements (T3, Q4).');
end
t  = mesh.t;
nEdges = size(edge_nodes, 1);

if isscalar(tx), tx = repmat(tx, nEdges, 1); end
if isscalar(ty), ty = repmat(ty, nEdges, 1); end

for k = 1:nEdges
    ni = edge_nodes(k, 1);
    nj = edge_nodes(k, 2);
    xi = mesh.nodes(ni, :);
    xj = mesh.nodes(nj, :);
    L  = norm(xj - xi);
    fk = (L * t / 2) * [tx(k); ty(k)];  % force at each endpoint

    % DOF indices
    dof_i = [2*ni-1, 2*ni];
    dof_j = [2*nj-1, 2*nj];
    f(dof_i) = f(dof_i) + fk';
    f(dof_j) = f(dof_j) + fk';
end
end
