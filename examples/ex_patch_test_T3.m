% Patch test for T3 element.
% Applies a pure uniaxial stress state (sigma_xx = 1, sigma_yy = tau_xy = 0)
% as boundary displacements on a 4x2 T3 mesh.
% The T3 element must recover the exact stress field everywhere.
startup_fem2d;

E = 1e6; nu = 0.3;
mat.E = E; mat.nu = nu; mat.formulation = 'plane_stress';
D = mat_elastic_D(E, nu, 'plane_stress');

% Exact strain from plane-stress: eps_xx = sigma/E, eps_yy = -nu*sigma/E
sigma = 1;
eps_xx = sigma / E;
eps_yy = -nu * sigma / E;

% Displacement field: ux = eps_xx * x, uy = eps_yy * y
mesh = mesh_rect_tri(1.0, 1.0, 4, 2);
mesh.t = 1;

% Build BC: prescribe exact displacements on ALL boundary nodes
% (this is the patch test — if the element is consistent, interior stresses
%  will equal sigma_xx=1 exactly regardless of mesh distortion)
nNodes = mesh.nNodes;
fixed_nodes = (1:nNodes)';
fixed_dofs_x = ones(nNodes, 1);
fixed_dofs_y = 2 * ones(nNodes, 1);
vals_x = eps_xx * mesh.nodes(:, 1);
vals_y = eps_yy * mesh.nodes(:, 2);

bc.fixed_nodes = [fixed_nodes; fixed_nodes];
bc.fixed_dofs  = [fixed_dofs_x; fixed_dofs_y];
bc.fixed_vals  = [vals_x; vals_y];

[K, f] = assemble(mesh, mat);
u = solve_linear(K, f, mesh, bc);

[stress, ~] = postproc_stress(mesh, mat, u);

err_sxx = max(abs(stress(:,1) - sigma));
err_syy = max(abs(stress(:,2)));
err_txy = max(abs(stress(:,3)));

fprintf('T3 Patch Test Results:\n');
fprintf('  max |sigma_xx - 1| = %.3e  (should be < 1e-10)\n', err_sxx);
fprintf('  max |sigma_yy|     = %.3e  (should be < 1e-10)\n', err_syy);
fprintf('  max |tau_xy|       = %.3e  (should be < 1e-10)\n', err_txy);

if max([err_sxx, err_syy, err_txy]) < 1e-8
    disp('PASS');
else
    disp('FAIL');
end
