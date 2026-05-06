% Patch test for Q4 element.
startup_fem2d;

E = 1e6; nu = 0.3;
mat.E = E; mat.nu = nu; mat.formulation = 'plane_stress';

sigma = 1;
eps_xx = sigma / E;
eps_yy = -nu * sigma / E;

mesh = mesh_rect_quad(1.0, 1.0, 4, 2);
mesh.t = 1;

nNodes = mesh.nNodes;
bc.fixed_nodes = [repmat((1:nNodes)', 2, 1)];
bc.fixed_dofs  = [ones(nNodes,1); 2*ones(nNodes,1)];
bc.fixed_vals  = [eps_xx * mesh.nodes(:,1); eps_yy * mesh.nodes(:,2)];

[K, f] = assemble(mesh, mat);
u = solve_linear(K, f, mesh, bc);
[stress, ~] = postproc_stress(mesh, mat, u);

err_sxx = max(abs(stress(:,1) - sigma));
err_syy = max(abs(stress(:,2)));
err_txy = max(abs(stress(:,3)));

fprintf('Q4 Patch Test Results:\n');
fprintf('  max |sigma_xx - 1| = %.3e  (should be < 1e-10)\n', err_sxx);
fprintf('  max |sigma_yy|     = %.3e  (should be < 1e-10)\n', err_syy);
fprintf('  max |tau_xy|       = %.3e  (should be < 1e-10)\n', err_txy);

if max([err_sxx, err_syy, err_txy]) < 1e-8
    disp('PASS');
else
    disp('FAIL');
end
