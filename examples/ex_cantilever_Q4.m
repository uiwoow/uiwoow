% Cantilever beam under tip point load — Q4 mesh.
% Compares FEM tip deflection to Euler-Bernoulli beam theory.
%
%   Fixed left edge (x=0), tip load P downward at (L, h/2).
%   Expected tip deflection: delta = P*L^3 / (3*E*I)  (slender beam)
startup_fem2d;

L  = 10;    % beam length
h  = 1;     % beam height
t  = 1;     % thickness
E  = 1e7;
nu = 0.3;
P  = -100;  % tip load (negative = downward in y-direction)

I_beam = t * h^3 / 12;
delta_theory = P * L^3 / (3 * E * I_beam);

nx = 20; ny = 4;
mesh = mesh_rect_quad(L, h, nx, ny);
mesh.t = t;
mat.E = E; mat.nu = nu; mat.formulation = 'plane_stress';

% Fix all nodes on left edge
left_nodes = mesh_find_nodes(mesh, 1e-9, 'x', 0);
bc.fixed_nodes = [left_nodes; left_nodes];
bc.fixed_dofs  = [ones(length(left_nodes),1); 2*ones(length(left_nodes),1)];
bc.fixed_vals  = zeros(2*length(left_nodes), 1);

% Distribute tip load across right-edge nodes
right_nodes = mesh_find_nodes(mesh, 1e-9, 'x', L);
n_right = length(right_nodes);
p_per_node = P / n_right;
loads.point_loads = [right_nodes, 2*ones(n_right,1), p_per_node*ones(n_right,1)];

[K, f] = assemble(mesh, mat, loads);
u = solve_linear(K, f, mesh, bc);

% Tip deflection: average uy at right edge mid-height
tip_node = mesh_find_nodes(mesh, 1e-9, 'xy', L, h/2);
dof_uy   = 2 * tip_node(1);
delta_fem = u(dof_uy);

err = abs((delta_fem - delta_theory) / delta_theory) * 100;
fprintf('Cantilever Q4:\n');
fprintf('  Theory tip deflection   = %.6f\n', delta_theory);
fprintf('  FEM    tip deflection   = %.6f\n', delta_fem);
fprintf('  Error                   = %.2f%%\n', err);

% Plots
[stress, ~] = postproc_stress(mesh, mat, u);
postproc_plot_disp(mesh, u);
vm = postproc_von_mises(stress);
postproc_plot_stress(mesh, [vm, zeros(mesh.nElems,2)], 1, 'Von Mises stress');
