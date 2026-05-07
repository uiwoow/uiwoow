% Jib crane with mixed-element model: Q4 wall-mount plate + beam arm.
%
% Geometry (all dimensions in metres):
%
%   Wall (fixed)
%   |
%   |  Q4 plate (0.3 x 0.3 m, 3x2 elements)
%   |  left edge clamped
%   |
%   +----[plate right-edge centre (0.3, 0.15)]
%                  |
%               beam element 1 (L=0.6 m)
%                  |
%               (0.9, 0.15)
%                  |
%               beam element 2 (L=0.6 m)
%                  |
%               tip (1.5, 0.15)  <-- P downward
%
% The plate models a short, stocky gusset/bracket where beam theory
% is not applicable (aspect ratio < 5:1).  The beam models the slender
% jib arm where Euler-Bernoulli assumptions hold.
%
% Connection physics: because the beam end meets a single Q4 node,
% forces (axial and transverse) transfer but rotational stiffness does
% not.  The connection is effectively pinned.  For full moment transfer
% use a multi-node rigid-link — see the plan for future work.
startup_fem2d;

%% Parameters
E     = 200e9;    % steel [Pa]
nu    = 0.3;
t_p   = 0.02;     % plate thickness  [m]
A_b   = 4e-4;     % beam cross-section area [m^2]  (20 x 20 mm)
Iz_b  = 1.333e-8; % beam second moment [m^4]
P     = 20e3;     % tip load [N]

%% Build Q4 wall-mount plate (nx=3, ny=2 -> 4x3 = 12 nodes)
Lp_x = 0.3;  nx_p = 3;
Lp_y = 0.3;  ny_p = 2;
plate = mesh_rect_quad(Lp_x, Lp_y, nx_p, ny_p);

% Locate the right-edge centre node (x=Lp_x, y=Lp_y/2)
n_conn = mesh_find_nodes(plate, 1e-9, 'xy', Lp_x, Lp_y/2);
if isempty(n_conn)
    error('Connection node not found. Check that ny_p is even so y=Lp_y/2 is a node.');
end
n_conn = n_conn(1);
fprintf('Connection node: %d  at (%.3f, %.3f)\n', ...
        n_conn, plate.nodes(n_conn,1), plate.nodes(n_conn,2));

%% Append beam nodes (two extra nodes, two elements)
beam_nodes = [0.9, Lp_y/2;
              1.5, Lp_y/2];
all_nodes = [plate.nodes; beam_nodes];
n_mid = plate.nNodes + 1;
n_tip = plate.nNodes + 2;

%% Assemble mixed mesh
m = mesh_mixed(all_nodes);
m.t = t_p;

mat_plate.E = E;  mat_plate.nu = nu;  mat_plate.formulation = 'plane_stress';
m = mesh_add_elements(m, 'Q4', plate.conn, mat_plate);

beam_conn = [n_conn, n_mid;
             n_mid,  n_tip];
mat_beam.E = E;  mat_beam.A = A_b;  mat_beam.Iz = Iz_b;
m = mesh_add_elements(m, 'beam2d', beam_conn, mat_beam);

%% Boundary conditions: clamp left edge of plate (x=0)
left_nodes = mesh_find_nodes(plate, 1e-9, 'x', 0);
nL = numel(left_nodes);
bc.fixed_nodes = [left_nodes; left_nodes];
bc.fixed_dofs  = [ones(nL,1); 2*ones(nL,1)];
bc.fixed_vals  = zeros(2*nL, 1);

%% Load: downward force at beam tip (local DOF 2 = v in 3-DOF system)
loads.point_loads = [n_tip, 2, -P];

%% Solve
[K, f] = assemble(m, [], loads);
u = solve_linear(K, f, m, bc);

%% Extract key results
u_conn_v = u(3*(n_conn-1) + 2);
u_tip_v  = u(3*(n_tip-1)  + 2);

fprintf('\n=== Jib Crane Results ===\n');
fprintf('  Connection node vertical displacement : %+.4f mm\n', u_conn_v*1e3);
fprintf('  Beam tip vertical displacement        : %+.4f mm\n', u_tip_v*1e3);

%% Beam internal forces
stress_cells = postproc_stress(m, [], u);
n_q4 = plate.nElems;
fprintf('\n  Beam internal forces:\n');
for e = 1:numel(beam_conn)
    s = stress_cells{n_q4 + e};
    fprintf('    Beam elem %d:  N=%+.2f kN   V=%+.2f kN   M1=%+.2f kN.m   M2=%+.2f kN.m\n', ...
            e, s(1)/1e3, s(2)/1e3, s(3)/1e3, s(4)/1e3);
end

%% Plate max von-Mises stress
sxx = cellfun(@(s) s(1), stress_cells(1:n_q4));
syy = cellfun(@(s) s(2), stress_cells(1:n_q4));
txy = cellfun(@(s) s(3), stress_cells(1:n_q4));
vm  = sqrt(sxx.^2 - sxx.*syy + syy.^2 + 3*txy.^2);
fprintf('\n  Plate max von-Mises stress : %.2f MPa\n', max(vm)/1e6);

%% Beam tip deflection vs. pure Euler-Bernoulli (no plate compliance)
L_beam = 1.2;   % total beam length
delta_EB = P * L_beam^3 / (3 * E * Iz_b);
fprintf('\n  Pure EB tip deflection (no plate): %.4f mm\n', delta_EB*1e3);
fprintf('  Mixed model tip deflection       : %.4f mm\n', abs(u_tip_v)*1e3);
fprintf('  Plate compliance added           : %.4f mm\n', ...
        abs(u_tip_v)*1e3 - delta_EB*1e3);

%% Visualisation
figure('Name','Jib crane deformed shape');
mesh_plot(m, 'deformed', u, 'scale', 200);
title(sprintf('Jib crane — deformed x200  (tip disp = %.2f mm)', abs(u_tip_v)*1e3));
