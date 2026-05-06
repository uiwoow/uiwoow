function u = solve_linear(K, f, mesh, bc)
% Solve the linear static FEM system K*u = f with boundary conditions.
%   K    - [nDof x nDof] global stiffness (sparse)
%   f    - [nDof x 1]    global force vector
%   mesh - mesh struct
%   bc   - BC struct (fixed_nodes, fixed_dofs, fixed_vals)
%
%   u    - [nDof x 1] full displacement vector
%
% Applies Dirichlet BCs by elimination, then solves the reduced system
% using MATLAB's sparse direct solver (backslash / UMFPACK).

[K_free, f_free, free_dofs, u] = bc_apply(K, f, mesh, bc);

% Check for singularity (unfixed rigid body modes)
if condest(K_free) > 1e14
    warning('solve_linear: K_free appears near-singular (condest > 1e14). Check BCs.');
end

u(free_dofs) = K_free \ f_free;
end
