function pass = test_Q4_patch()
% Q4 patch test: uniform uniaxial stress field must be recovered exactly.
startup_fem2d;

E = 1e6; nu = 0.3;
mat.E = E; mat.nu = nu; mat.formulation = 'plane_stress';
sigma = 1;
eps_xx = sigma / E; eps_yy = -nu * sigma / E;

mesh = mesh_rect_quad(1.0, 1.0, 4, 2);
mesh.t = 1;

nNodes = mesh.nNodes;
bc.fixed_nodes = repmat((1:nNodes)', 2, 1);
bc.fixed_dofs  = [ones(nNodes,1); 2*ones(nNodes,1)];
bc.fixed_vals  = [eps_xx * mesh.nodes(:,1); eps_yy * mesh.nodes(:,2)];

[K, f] = assemble(mesh, mat);
u = solve_linear(K, f, mesh, bc);
[stress, ~] = postproc_stress(mesh, mat, u);

err = max([max(abs(stress(:,1) - sigma)), max(abs(stress(:,2))), max(abs(stress(:,3)))]);
fprintf('test_Q4_patch: max stress error = %.3e\n', err);
pass = err < 1e-8;
if pass, disp('  PASS'); else, disp('  FAIL'); end
end
