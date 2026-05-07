function pass = test_rigid_link()
% Verify that a rigid link transfers moment from a beam to a Q4 plate.
%
% Setup:
%   Very stiff Q4 plate (E=1e15), 0.2x0.2 m, 4x4 elements, left edge clamped.
%   Single beam element (E=200e9) from plate right-edge centre to tip, L=1 m.
%   Tip load P=1000 N downward.
%
% Two models:
%   (a) Pin connection  - single node, no rigid link
%   (b) Rigid link      - all right-edge nodes tied to beam-end node
%
% Expected behaviour:
%   The rigid link constrains beam-end rotation (plate resists it via
%   bending).  This makes the effective boundary condition closer to a
%   clamp, so:
%     |theta_conn| with rigid link  <  |theta_conn| without
%     |u_tip|      with rigid link  <  |u_tip|      without
%
%   Analytical reference (fully clamped beam):
%     delta_clamp = P*L^3 / (3*E_beam*Iz)
startup_fem2d;

% --- parameters ---
E_plate = 1e15;  nu = 0.3;  t = 0.01;
Lp = 0.2;  np = 4;
E_beam = 200e9;  A_beam = 1e-3;  Iz_beam = 8.33e-9;
L_beam = 1.0;
P = 1000;

% --- shared plate mesh ---
plate = mesh_rect_quad(Lp, Lp, np, np);
n_conn = mesh_find_nodes(plate, 1e-9, 'xy', Lp, Lp/2);
n_conn = n_conn(1);

% --- shared extra beam tip node ---
tip_xy = [Lp + L_beam, Lp/2];
all_nodes = [plate.nodes; tip_xy];
n_tip = plate.nNodes + 1;

% === Model A: pin connection (no rigid link) ===
mA = mesh_mixed(all_nodes);
mA.t = t;
mat_plate.E = E_plate; mat_plate.nu = nu; mat_plate.formulation = 'plane_stress';
mA = mesh_add_elements(mA, 'Q4',     plate.conn,         mat_plate);
mat_beam.E = E_beam; mat_beam.A = A_beam; mat_beam.Iz = Iz_beam;
mA = mesh_add_elements(mA, 'beam2d', [n_conn, n_tip],    mat_beam);

left_nodes = mesh_find_nodes(plate, 1e-9, 'x', 0);
nL = numel(left_nodes);
bc.fixed_nodes = [left_nodes; left_nodes];
bc.fixed_dofs  = [ones(nL,1); 2*ones(nL,1)];
bc.fixed_vals  = zeros(2*nL, 1);
loads.point_loads = [n_tip, 2, -P];

[KA, fA] = assemble(mA, [], loads);
uA = solve_linear(KA, fA, mA, bc);

theta_conn_A = uA(3*n_conn);          % rotation at connection (pin)
u_tip_A      = uA(3*(n_tip-1) + 2);  % vertical tip displacement

% === Model B: rigid link (right-edge plate nodes tied to beam end) ===
mB = mA;  % identical mesh

right_nodes = mesh_find_nodes(plate, 1e-9, 'x', Lp);
slave_nodes = setdiff(right_nodes, n_conn);  % exclude master itself

[KB, fB] = assemble(mB, [], loads);
[KB, fB] = constraint_rigid_link_apply(KB, fB, mB, n_conn, slave_nodes);
uB = solve_linear(KB, fB, mB, bc);

theta_conn_B = uB(3*n_conn);
u_tip_B      = uB(3*(n_tip-1) + 2);

% === Analytical reference: perfectly clamped beam ===
delta_clamp = -P * L_beam^3 / (3 * E_beam * Iz_beam);

fprintf('test_rigid_link:\n');
fprintf('  Connection rotation (pin)        : theta = %.4e rad\n', theta_conn_A);
fprintf('  Connection rotation (rigid link) : theta = %.4e rad\n', theta_conn_B);
fprintf('  Tip disp (pin)        : %.4e m\n', u_tip_A);
fprintf('  Tip disp (rigid link) : %.4e m\n', u_tip_B);
fprintf('  Tip disp (clamp ref)  : %.4e m\n', delta_clamp);

% Checks
chk1 = abs(theta_conn_B) < abs(theta_conn_A);   % rigid link reduces rotation
chk2 = abs(u_tip_B)      < abs(u_tip_A);        % stiffer -> smaller deflection
chk3 = abs(u_tip_B)      > 0;                   % model solved (non-trivial)
% Rigid-link model should approach the clamped solution as plate stiffens:
chk4 = abs(u_tip_B - delta_clamp) / abs(delta_clamp) < 0.05;  % within 5%

pass = chk1 && chk2 && chk3 && chk4;
if pass
    disp('  PASS');
else
    fprintf('  FAIL (chk1=%d chk2=%d chk3=%d chk4=%d)\n', chk1, chk2, chk3, chk4);
end
end
