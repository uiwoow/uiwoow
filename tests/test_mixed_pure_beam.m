function pass = test_mixed_pure_beam()
% Regression: a single beam element expressed as a mixed mesh must give
% the same tip deflection as the same model using the native beam2d type.
%
% Cantilever: length L=1, E=1, A=1, Iz=1/12, tip load P=1.
% Analytical (Euler-Bernoulli): delta = P*L^3 / (3*E*I) = 4.
startup_fem2d;

L = 1;  E = 1;  A = 1;  Iz = 1/12;  P = 1;

% --- reference: native beam2d mesh ---
nodes = [0 0; L 0];
conn  = [1 2];
ref_mesh = mesh_frame(nodes, conn, 'beam2d');
mat.E = E; mat.A = A; mat.Iz = Iz;

bc_ref.fixed_nodes = [1; 1; 1];
bc_ref.fixed_dofs  = [1; 2; 3];
bc_ref.fixed_vals  = zeros(3, 1);

loads_ref.point_loads = [2, 2, -P];  % downward at tip
[K_ref, f_ref] = assemble(ref_mesh, mat, loads_ref);
u_ref = solve_linear(K_ref, f_ref, ref_mesh, bc_ref);
tip_ref = u_ref(end-1);  % v at node 2

% --- mixed mesh with identical beam element ---
m = mesh_mixed(nodes);
mat_b.E = E; mat_b.A = A; mat_b.Iz = Iz;
m = mesh_add_elements(m, 'beam2d', conn, mat_b);

bc_mix.fixed_nodes = [1; 1; 1];
bc_mix.fixed_dofs  = [1; 2; 3];
bc_mix.fixed_vals  = zeros(3, 1);

loads_mix.point_loads = [2, 2, -P];  % local_dof 2 = v in 3-DOF system
[K_mix, f_mix] = assemble(m, [], loads_mix);
u_mix = solve_linear(K_mix, f_mix, m, bc_mix);
% In the 3-DOF system, node 2 DOFs are [3*2-2, 3*2-1, 3*2] = [4, 5, 6]
% DOF 5 = v at node 2
tip_mix = u_mix(5);

err = abs(tip_mix - tip_ref);
analytical = -P*L^3 / (3*E*Iz);
err_vs_analytical = abs(tip_mix - analytical);

fprintf('test_mixed_pure_beam:\n');
fprintf('  tip_ref = %.6f   tip_mix = %.6f   analytical = %.6f\n', ...
        tip_ref, tip_mix, analytical);
fprintf('  diff vs reference = %.2e\n', err);
fprintf('  diff vs analytical= %.2e\n', err_vs_analytical);

pass = err < 1e-12 && err_vs_analytical < 1e-10;
if pass, disp('  PASS'); else, disp('  FAIL'); end
end
