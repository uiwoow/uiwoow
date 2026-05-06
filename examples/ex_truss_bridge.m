% Simple 2D Pratt truss bridge under mid-span point load.
%
%  Node layout (6 top + 6 bottom = 12 nodes):
%  Top chord: y=H, Bottom chord: y=0, Span = 5 panels of width W each
%
%       3---4---5---6---7    (top chord, y=H)
%      /|\ /|\ /|\ /|\ /|\
%     1--2--8--9-10-11-12   (bottom chord, y=0)  <- fixed at 1 and 12
%
%  Loads: P downward at nodes 3 and 7 (third points)

startup_fem2d;

W = 1;     % panel width
H = 1;     % truss height
E = 2e11;  % steel Young's modulus [Pa]
A = 1e-3;  % cross-section area [m^2]
P = -1e4;  % load per loaded node [N]

% Node coordinates
%  Bottom chord: 0..5W at y=0  (nodes 1..6)
%  Top chord:    0..5W at y=H  (nodes 7..12)
x_bot = (0:5) * W;  y_bot = zeros(1,6);
x_top = (0:5) * W;  y_top = H * ones(1,6);
nodes = [x_bot', y_bot'; x_top', y_top'];
% Nodes: 1-6 bottom, 7-12 top

% Connectivity: bottom chord + top chord + verticals + diagonals
bottom = [(1:5)', (2:6)'];
top    = [(7:11)', (8:12)'];
vert   = [(1:6)', (7:12)'];
% Diagonals (Pratt: tension diagonals slope toward mid-span)
diag_left  = [(2:4)', (8:10)'];   % left half: lean right
diag_right = [(4:5)', (10:11)']; % right half: lean left
diag_left2 = [(8:10)', (3:5)'];  % crossing diagonals

conn = [bottom; top; vert; diag_left; diag_right];

mesh = mesh_frame(nodes, conn, 'truss2d');
mat.E = E; mat.A = A;

% Fix supports: node 1 (pin: ux=uy=0), node 6 (roller: uy=0)
bc.fixed_nodes = [1; 1; 6];
bc.fixed_dofs  = [1; 2; 2];
bc.fixed_vals  = [0; 0; 0];

% Load: downward force at mid-bottom nodes (nodes 3 and 4)
loads.point_loads = [3, 2, P; 4, 2, P];

[K, f] = assemble(mesh, mat, loads);
u = solve_linear(K, f, mesh, bc);

[forces, ~] = postproc_stress(mesh, mat, u);
axial = forces(:, 1);

fprintf('Truss Bridge Results:\n');
fprintf('  Mid-bottom node 3 deflection: uy = %.4e m\n', u(2*3));
fprintf('  Mid-bottom node 4 deflection: uy = %.4e m\n', u(2*4));
fprintf('  Max tensile  member force:    N  = %.4e N\n', max(axial));
fprintf('  Max compressive member force: N  = %.4e N\n', min(axial));

figure;
mesh_plot(mesh, 'deformed', u, 'scale', 50);
title('Truss bridge deformed shape (50x scale)');

postproc_plot_stress(mesh, [axial, zeros(mesh.nElems,1)], 1, 'Axial force N [N]');
