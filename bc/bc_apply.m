function [K_free, f_free, free_dofs, u] = bc_apply(K, f, mesh, bc)
% Apply Dirichlet boundary conditions by DOF elimination (partition method).
%   K    - [nDof x nDof] global stiffness (sparse)
%   f    - [nDof x 1]    global force vector
%   mesh - mesh struct (used for DOF numbering)
%   bc   - boundary condition struct:
%          bc.fixed_nodes  - [nFixed x 1] node indices
%          bc.fixed_dofs   - [nFixed x 1] local DOF per node (1,2 or 1,2,3)
%          bc.fixed_vals   - [nFixed x 1] prescribed displacement values
%
%   K_free   - reduced stiffness matrix (free DOFs only)
%   f_free   - reduced force vector (free DOFs only)
%   free_dofs - indices of free DOFs in global system
%   u         - [nDof x 1] displacement vector with prescribed values filled in

nDof = size(K, 1);

% Convert (node, local_dof) pairs to global DOF indices
constrained = zeros(length(bc.fixed_nodes), 1);
for i = 1:length(bc.fixed_nodes)
    node = bc.fixed_nodes(i);
    ldof = bc.fixed_dofs(i);
    switch mesh.type
        case {'T3', 'Q4', 'truss2d'}
            constrained(i) = 2*(node-1) + ldof;
        case {'beam2d', 'mixed'}
            % mixed uses 3 DOFs/node universally: 1=u, 2=v, 3=theta
            constrained(i) = 3*(node-1) + ldof;
    end
end

vals = bc.fixed_vals(:);
free_dofs = setdiff((1:nDof)', constrained);

% Adjust RHS for prescribed values: f_free -= K_fc * u_prescribed
f_mod = f - K(:, constrained) * vals;

K_free = K(free_dofs, free_dofs);
f_free = f_mod(free_dofs);

% Pre-fill prescribed values in solution vector
u = zeros(nDof, 1);
u(constrained) = vals;
end
