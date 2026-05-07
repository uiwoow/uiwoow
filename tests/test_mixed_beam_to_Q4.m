function pass = test_mixed_beam_to_Q4()
% Verify force transfer between a beam element and a Q4 plate.
%
% Geometry:
%   Q4 plate  : 0.2 x 0.2 m, 2x2 elements, left edge clamped.
%   Beam arm  : 2 elements, 0.6 m long, attached to the center of the
%               plate's right edge at (0.2, 0.1).
%   Load      : P = 1000 N downward at beam tip.
%
% Checks:
%   1. Tip and connection nodes both displace downward.
%   2. |u_conn_v| < |u_tip_v|  (connection stiffer due to plate support).
%   3. Global vertical equilibrium: sum of wall reactions == P.
startup_fem2d;

E  = 200e9;  nu = 0.3;
Lp = 0.2;    np = 2;      % plate side length and element count per side
Lb = 0.6;    nb = 2;      % beam length and number of elements
t  = 0.01;                 % plate thickness [m]
A  = 0.001;  Iz = 8.33e-9; % beam section (10x10 mm square)
P  = 1000;                  % tip load [N]

% --- Build Q4 plate mesh ---
plate = mesh_rect_quad(Lp, Lp, np, np);

% --- Find the right-edge center node (x=Lp, y=Lp/2) ---
n_conn = mesh_find_nodes(plate, 1e-9, 'xy', Lp, Lp/2);
if isempty(n_conn)
    error('test_mixed_beam_to_Q4: connection node not found in plate mesh.');
end
n_conn = n_conn(1);

% --- Append beam nodes to the shared node list ---
dx = Lb / nb;
beam_extra = zeros(nb, 2);
for k = 1:nb
    beam_extra(k,:) = [Lp + k*dx, Lp/2];
end
all_nodes = [plate.nodes; beam_extra];
n_tip = plate.nNodes + nb;  % last beam node

% --- Build mixed mesh ---
m = mesh_mixed(all_nodes);
m.t = t;

mat_plate.E = E; mat_plate.nu = nu; mat_plate.formulation = 'plane_stress';
m = mesh_add_elements(m, 'Q4', plate.conn, mat_plate);

beam_node_seq = [n_conn, plate.nNodes + (1:nb)];
beam_conn = [beam_node_seq(1:end-1)', beam_node_seq(2:end)'];
mat_beam.E = E; mat_beam.A = A; mat_beam.Iz = Iz;
m = mesh_add_elements(m, 'beam2d', beam_conn, mat_beam);

% --- Boundary conditions: clamp left edge of plate (x=0) ---
left_nodes = mesh_find_nodes(plate, 1e-9, 'x', 0);
nL = numel(left_nodes);
% Also fix theta at n_conn: the Q4 plate has no rotational DOFs so it
% provides no moment resistance at the connection node.  Without this
% constraint the beam can rotate as a rigid body about n_conn (zero-energy
% mechanism), making K_free singular.  Fixing theta=0 models a clamped
% (welded) beam-to-plate joint.
bc.fixed_nodes = [left_nodes; left_nodes; n_conn];
bc.fixed_dofs  = [ones(nL,1); 2*ones(nL,1); 3];
bc.fixed_vals  = zeros(2*nL+1, 1);

% --- Tip load: downward (v-DOF = local dof 2 in 3-DOF system) ---
loads.point_loads = [n_tip, 2, -P];

% --- Solve ---
[K, f] = assemble(m, [], loads);
u      = solve_linear(K, f, m, bc);

% --- Extract displacements ---
% DOF layout: node n -> [3n-2, 3n-1, 3n] = [u, v, theta]
u_conn_v = u(3*(n_conn-1) + 2);
u_tip_v  = u(3*(n_tip-1)  + 2);

% --- Global vertical equilibrium via reaction forces ---
% Reactions at clamped DOFs: R = K*u - f_ext  (at constrained DOFs)
R_all = K * u - f;
v_dofs_left = 3*(left_nodes-1) + 2;
R_vert = sum(R_all(v_dofs_left));

fprintf('test_mixed_beam_to_Q4:\n');
fprintf('  u_conn_v = %.4e m  (connection node vertical disp)\n', u_conn_v);
fprintf('  u_tip_v  = %.4e m  (beam tip vertical disp)\n',       u_tip_v);
fprintf('  wall vertical reaction = %.4f N  (should be ~%.1f)\n', R_vert, P);

check1 = u_tip_v  < 0;
check2 = u_conn_v < 0;
check3 = abs(u_conn_v) < abs(u_tip_v);
check4 = abs(R_vert - P) / P < 0.01;  % reaction within 1% of applied load

pass = check1 && check2 && check3 && check4;
if pass, disp('  PASS'); else
    fprintf('  FAIL (checks: tip<0=%d  conn<0=%d  |conn|<|tip|=%d  equil=%d)\n', ...
            check1, check2, check3, check4);
end
end
