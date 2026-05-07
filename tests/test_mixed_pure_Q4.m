function pass = test_mixed_pure_Q4()
% Regression: a Q4 patch expressed as a mixed mesh must recover the same
% uniform stress field as the native Q4 patch test.
%
% Applies a uniaxial stress sigma_xx=1 via prescribed displacements and
% checks that all element stresses match the imposed field.
startup_fem2d;

E = 1e6;  nu = 0.3;
sigma = 1;
eps_xx =  sigma / E;
eps_yy = -nu * sigma / E;

% Generate Q4 mesh and pull its connectivity into a mixed mesh.
plate = mesh_rect_quad(1.0, 1.0, 3, 3);

m = mesh_mixed(plate.nodes);
m.t = 1;
mat_p.E = E; mat_p.nu = nu; mat_p.formulation = 'plane_stress';
m = mesh_add_elements(m, 'Q4', plate.conn, mat_p);

% Prescribe linear displacement field on all nodes (patch test setup).
% In the 3-DOF system, u-DOF = 3*(n-1)+1, v-DOF = 3*(n-1)+2.
nN = m.nNodes;
nodes_all = (1:nN)';
bc.fixed_nodes = [nodes_all; nodes_all];
bc.fixed_dofs  = [ones(nN,1); 2*ones(nN,1)];
bc.fixed_vals  = [eps_xx * m.nodes(:,1); eps_yy * m.nodes(:,2)];

[K, f] = assemble(m, [], []);
u      = solve_linear(K, f, m, bc);
stress_cells = postproc_stress(m, [], u);

% Collect sxx, syy, txy across all elements
sxx = cellfun(@(s) s(1), stress_cells);
syy = cellfun(@(s) s(2), stress_cells);
txy = cellfun(@(s) s(3), stress_cells);

err = max([max(abs(sxx - sigma)), max(abs(syy)), max(abs(txy))]);
fprintf('test_mixed_pure_Q4: max stress error = %.3e\n', err);

pass = err < 1e-6;
if pass, disp('  PASS'); else, disp('  FAIL'); end
end
