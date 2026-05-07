% Jib crane — upgraded with rigid link for proper moment transfer.
%
% Extends ex_jib_crane_mixed.m by tying all right-edge plate nodes to the
% beam-end (master) node via a rigid-link constraint.  This transmits the
% beam's end moment into the plate, capturing the stress concentration that
% would be missed by a simple pin connection.
%
% Geometry (metres):
%
%   Wall (fixed)
%   |
%   | Q4 plate 0.3 x 0.3 m  (3 x 2 elements)
%   | left edge clamped
%   |
%   +----[right-edge centre = master node (0.3, 0.15)]
%        |  <- rigid link to right-edge nodes above/below master
%        |
%     beam elem 1   (0.6 m)
%        |
%     beam elem 2   (0.6 m)
%        |
%      tip (1.5, 0.15)   <- P downward
%
% Run alongside ex_jib_crane_mixed.m to compare pin vs rigid-link results.
startup_fem2d;

%% Parameters
E     = 200e9;
nu    = 0.3;
t_p   = 0.02;
A_b   = 4e-4;
Iz_b  = 1.333e-8;
P     = 20e3;

%% Q4 plate (identical to ex_jib_crane_mixed)
Lp_x = 0.3;  nx_p = 3;
Lp_y = 0.3;  ny_p = 2;
plate = mesh_rect_quad(Lp_x, Lp_y, nx_p, ny_p);

n_conn = mesh_find_nodes(plate, 1e-9, 'xy', Lp_x, Lp_y/2);
n_conn = n_conn(1);

%% Append beam nodes
beam_nodes = [0.9, Lp_y/2; 1.5, Lp_y/2];
all_nodes  = [plate.nodes; beam_nodes];
n_mid = plate.nNodes + 1;
n_tip = plate.nNodes + 2;

%% Mixed mesh
m = mesh_mixed(all_nodes);
m.t = t_p;

mat_plate.E = E; mat_plate.nu = nu; mat_plate.formulation = 'plane_stress';
m = mesh_add_elements(m, 'Q4', plate.conn, mat_plate);

beam_conn = [n_conn, n_mid; n_mid, n_tip];
mat_beam.E = E; mat_beam.A = A_b; mat_beam.Iz = Iz_b;
m = mesh_add_elements(m, 'beam2d', beam_conn, mat_beam);

%% Boundary conditions
left_nodes = mesh_find_nodes(plate, 1e-9, 'x', 0);
nL = numel(left_nodes);
bc.fixed_nodes = [left_nodes; left_nodes];
bc.fixed_dofs  = [ones(nL,1); 2*ones(nL,1)];
bc.fixed_vals  = zeros(2*nL, 1);

loads.point_loads = [n_tip, 2, -P];

%% Assemble
[K, f] = assemble(m, [], loads);

%% Apply rigid link: tie all right-edge plate nodes to the beam-end master
right_nodes  = mesh_find_nodes(plate, 1e-9, 'x', Lp_x);
slave_nodes  = setdiff(right_nodes, n_conn);   % exclude master itself

[K, f] = constraint_rigid_link_apply(K, f, m, n_conn, slave_nodes);

%% Solve
u = solve_linear(K, f, m, bc);

%% Results
u_conn_v  = u(3*(n_conn-1) + 2);
u_conn_th = u(3*(n_conn-1) + 3);
u_tip_v   = u(3*(n_tip-1)  + 2);

fprintf('\n=== Jib Crane — Rigid Link Results ===\n');
fprintf('  Connection vertical disp    : %+.4f mm\n', u_conn_v*1e3);
fprintf('  Connection rotation (theta) : %+.4f mrad\n', u_conn_th*1e3);
fprintf('  Beam tip vertical disp      : %+.4f mm\n', u_tip_v*1e3);

stress_cells = postproc_stress(m, [], u);
n_q4 = plate.nElems;
fprintf('\n  Beam internal forces:\n');
for e = 1:size(beam_conn,1)
    s = stress_cells{n_q4 + e};
    fprintf('    Beam elem %d:  N=%+.2f kN  V=%+.2f kN  M1=%+.2f kN.m  M2=%+.2f kN.m\n', ...
            e, s(1)/1e3, s(2)/1e3, s(3)/1e3, s(4)/1e3);
end

sxx = cellfun(@(s) s(1), stress_cells(1:n_q4));
syy = cellfun(@(s) s(2), stress_cells(1:n_q4));
txy = cellfun(@(s) s(3), stress_cells(1:n_q4));
vm  = sqrt(sxx.^2 - sxx.*syy + syy.^2 + 3*txy.^2);
fprintf('\n  Plate max von-Mises : %.2f MPa\n', max(vm)/1e6);

%% Analytical pure-beam tip deflection
L_beam  = 1.2;
delta_EB = P * L_beam^3 / (3*E*Iz_b);
fprintf('\n  Pure EB tip (no plate) : %.4f mm\n', delta_EB*1e3);
fprintf('  Rigid-link model tip   : %.4f mm\n',   abs(u_tip_v)*1e3);
fprintf('  Stiffening from plate  : %.4f mm less deflection than pure EB\n', ...
        delta_EB*1e3 - abs(u_tip_v)*1e3);

%% Visualisation
figure('Name', 'Jib crane — rigid link');
postproc_plot_mixed(m, u, stress_cells);
sgtitle('Jib crane  (rigid-link moment transfer)');
