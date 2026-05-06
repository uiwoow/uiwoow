function pass = test_beam_cantilever()
% Cantilever beam: compare FEM tip deflection to delta = P*L^3 / (3*E*I).
startup_fem2d;

E = 1; Iz = 1; A = 1; L = 1; P = 1;
delta_theory = P * L^3 / (3 * E * Iz);

nodes = [0 0; L 0];
conn  = [1 2];
mesh  = mesh_frame(nodes, conn, 'beam2d');
mat.E = E; mat.A = A; mat.Iz = Iz;

% Fully fixed at node 1 (u, v, theta)
bc.fixed_nodes = [1; 1; 1];
bc.fixed_dofs  = [1; 2; 3];
bc.fixed_vals  = [0; 0; 0];

% Transverse load P at node 2, dof 2 (v)
loads.point_loads = [2, 2, P];

[K, f] = assemble(mesh, mat, loads);
u = solve_linear(K, f, mesh, bc);

delta_fem = u(3*2 - 1);  % uy at node 2
err = abs(delta_fem - delta_theory) / abs(delta_theory);

fprintf('test_beam_cantilever: delta_fem=%.6e  delta_theory=%.6e  rel_err=%.2e\n', ...
        delta_fem, delta_theory, err);
pass = err < 1e-10;
if pass, disp('  PASS'); else, disp('  FAIL'); end
end
