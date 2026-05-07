function pass = test_truss_bar()
% Single-bar extension test: F=EA*delta/L.
% EA=1, L=1 -> delta = F*L/(EA) = F for unit inputs.
startup_fem2d;

E = 1; A = 1; L = 1; F = 1;

nodes = [0 0; L 0];
conn  = [1 2];
mesh  = mesh_frame(nodes, conn, 'truss2d');
mat.E = E; mat.A = A;

bc.fixed_nodes = [1; 1; 2];
bc.fixed_dofs  = [1; 2; 2];
bc.fixed_vals  = [0; 0; 0];

loads.point_loads = [2, 1, F];

[K, f] = assemble(mesh, mat, loads);
u = solve_linear(K, f, mesh, bc);

delta_fem    = u(3);  % ux at node 2
delta_theory = F * L / (E * A);
err = abs(delta_fem - delta_theory);

fprintf('test_truss_bar: delta_fem=%.6e  delta_theory=%.6e  err=%.2e\n', ...
        delta_fem, delta_theory, err);
pass = err < 1e-12;
if pass, disp('  PASS'); else, disp('  FAIL'); end
end
