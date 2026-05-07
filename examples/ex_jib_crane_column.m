% ex_jib_crane_column.m  —  Free-standing column jib crane.
%
%   Structure:
%     Baseplate  →  column tube  →  top plate / jib root
%                                   support bracket (diagonal brace) attached
%                                   to column FACE  →  jib arm  →  tip load
%
%   Key modelling feature:
%     The diagonal brace welds to the column FACE, not the centreline.
%     A rigid-link MPC ties the face attachment node (slave, node 6) to the
%     column centreline node at the same height (master, node 2).  This means
%     that the column's bending rotation (theta_2) drives a vertical offset
%     at node 6, transferring a moment into the brace end.  Without the
%     rigid link the connection acts as a pin and under-estimates the moment
%     carried by the column.
%
%   Geometry sketch (x to the right, y upward, not to scale):
%
%     y
%     |
%   H |  3────────────4──────────5   ← jib arm
%     |  │           /
%  Hb |  2──────────6             ← rigid link 2↔6  (col. face node)
%     |  │
%   0 |  1                        ← column base, clamped (u=v=θ=0)
%     └────────────────────────── x
%        0         xb         Lj
%
%   Node list:
%     1  (0,   0  )  column base  [clamped]
%     2  (0,   Hb )  column CL at brace height  [rigid-link master]
%     3  (0,   H  )  column top / jib root
%     4  (xb,  H  )  jib at brace-to-jib junction
%     5  (Lj,  H  )  jib tip  [load applied here]
%     6  (cw,  Hb )  column face at brace height  [rigid-link slave]
%                    cw = column half-width
%
%   Element list:
%     1  beam2d  1→2   column lower
%     2  beam2d  2→3   column upper
%     3  beam2d  3→4   jib inner section
%     4  beam2d  4→5   jib outer section (to tip)
%     5  beam2d  6→4   diagonal brace
%
%   Adjust the geometry and section parameters below, then run:
%       startup_fem2d;  ex_jib_crane_column

startup_fem2d;

%% ── Geometry (edit these) ────────────────────────────────────────────────
H   = 4.0;   % total column height [m]
Hb  = 2.5;   % brace attachment height on column [m]  (Hb < H)
Lj  = 3.0;   % jib total length from column CL [m]
xb  = 2.0;   % jib x-coord of brace-to-jib junction [m]  (xb < Lj)
cw  = 0.10;  % column half-width [m]  (200 mm RHS → 0.10 m)

P   = 5000;  % tip load [N], applied downward at jib tip

%% ── Section properties ───────────────────────────────────────────────────
% All dimensions in mm, converted to SI at end of each block.
E = 200e9;   % Young's modulus, structural steel [Pa]

% Column:  200 × 200 × 8 mm RHS
h_c = 200; t_c = 8;
A_col  = (h_c^2 - (h_c - 2*t_c)^2) * 1e-6;              % [m²]
Iz_col = (h_c^4 - (h_c - 2*t_c)^4) / 12 * 1e-12;        % [m⁴]

% Jib arm:  150 × 100 × 6 mm RHS, strong axis vertical (h=150 mm)
h_j = 150; b_j = 100; t_j = 6;
A_jib  = (h_j*b_j - (h_j-2*t_j)*(b_j-2*t_j)) * 1e-6;
Iz_jib = (b_j*h_j^3 - (b_j-2*t_j)*(h_j-2*t_j)^3) / 12 * 1e-12;

% Diagonal brace:  80 × 80 × 5 mm RHS
h_b = 80; t_b = 5;
A_br   = (h_b^2 - (h_b - 2*t_b)^2) * 1e-6;
Iz_br  = (h_b^4 - (h_b - 2*t_b)^4) / 12 * 1e-12;

mat_col = struct('E', E, 'A', A_col,  'Iz', Iz_col);
mat_jib = struct('E', E, 'A', A_jib,  'Iz', Iz_jib);
mat_br  = struct('E', E, 'A', A_br,   'Iz', Iz_br);

fprintf('Section properties:\n');
fprintf('  Column  A=%.2e m²  Iz=%.2e m⁴\n', A_col,  Iz_col);
fprintf('  Jib     A=%.2e m²  Iz=%.2e m⁴\n', A_jib,  Iz_jib);
fprintf('  Brace   A=%.2e m²  Iz=%.2e m⁴\n', A_br,   Iz_br);

%% ── Nodes ────────────────────────────────────────────────────────────────
nodes = [
    0,   0;    % 1  column base (clamped)
    0,   Hb;   % 2  column CL at brace height  [rigid-link master]
    0,   H;    % 3  column top / jib root
    xb,  H;    % 4  jib at brace junction
    Lj,  H;    % 5  jib tip
    cw,  Hb;   % 6  column face at brace height  [rigid-link slave]
];

%% ── Mixed mesh ───────────────────────────────────────────────────────────
m = mesh_mixed(nodes);

m = mesh_add_elements(m, 'beam2d', [1 2; 2 3], mat_col);  % column (2 elem)
m = mesh_add_elements(m, 'beam2d', [3 4; 4 5], mat_jib);  % jib    (2 elem)
m = mesh_add_elements(m, 'beam2d', [6 4],      mat_br);   % brace  (1 elem)

%% ── Boundary conditions ──────────────────────────────────────────────────
% Fully clamp the column base (node 1): u=0, v=0, θ=0
bc.fixed_nodes = [1; 1; 1];
bc.fixed_dofs  = [1; 2; 3];
bc.fixed_vals  = zeros(3, 1);

%% ── Applied load ─────────────────────────────────────────────────────────
% Downward point load at jib tip (node 5, v-DOF = DOF index 2)
loads.point_loads = [5, 2, -P];

%% ── Assemble global K and f ──────────────────────────────────────────────
[K, f] = assemble(m, [], loads);

%% ── Rigid link: column face (slave=6) is rigidly attached to CL (master=2)
% Kinematics enforced by penalty:
%   u_6  =  u_2                       (ry = y6 - y2 = 0, no vertical offset)
%   v_6  =  v_2 + theta_2 * cw        (column rotation × face offset)
% When theta_2 ≠ 0 the face node gets extra vertical displacement, which
% applies a real moment couple into the brace — this is the effect that
% disappears if you skip the rigid link and merge nodes 2 and 6.
[K, f] = constraint_rigid_link_apply(K, f, m, 2, 6);

%% ── Solve ────────────────────────────────────────────────────────────────
u = solve_linear(K, f, m, bc);

%% ── Print results ────────────────────────────────────────────────────────
% DOF layout in mixed mesh: node n → [3n-2, 3n-1, 3n] = [u, v, θ]
dof = @(n) 3*(n-1) + (1:3);

u1 = u(dof(1));   u2 = u(dof(2));   u3 = u(dof(3));
u4 = u(dof(4));   u5 = u(dof(5));   u6 = u(dof(6));

fprintf('\n=== Jib crane — column with rigid-link brace attachment ===\n');
fprintf('Geometry:  H=%.1f m  Hb=%.1f m  Lj=%.1f m  xb=%.1f m  cw=%.3f m\n', ...
        H, Hb, Lj, xb, cw);

brace_ang = atan2d(nodes(4,2)-nodes(6,2), nodes(4,1)-nodes(6,1));
fprintf('Brace angle from horizontal: %.1f°\n', brace_ang);
fprintf('Tip load: P = %.0f N\n', P);

fprintf('\nDisplacements:\n');
fprintf('  Node 1  col. base       u=%+.3f mm  v=%+.3f mm  θ=%+.4f mrad\n', u1(1)*1e3, u1(2)*1e3, u1(3)*1e3);
fprintf('  Node 2  col. CL @ Hb   u=%+.3f mm  v=%+.3f mm  θ=%+.4f mrad\n', u2(1)*1e3, u2(2)*1e3, u2(3)*1e3);
fprintf('  Node 3  col. top        u=%+.3f mm  v=%+.3f mm  θ=%+.4f mrad\n', u3(1)*1e3, u3(2)*1e3, u3(3)*1e3);
fprintf('  Node 4  brace junction  u=%+.3f mm  v=%+.3f mm  θ=%+.4f mrad\n', u4(1)*1e3, u4(2)*1e3, u4(3)*1e3);
fprintf('  Node 5  jib tip         u=%+.3f mm  v=%+.3f mm  θ=%+.4f mrad\n', u5(1)*1e3, u5(2)*1e3, u5(3)*1e3);
fprintf('  Node 6  col. face @ Hb  u=%+.3f mm  v=%+.3f mm  θ=%+.4f mrad\n', u6(1)*1e3, u6(2)*1e3, u6(3)*1e3);

% Base reactions: R = K*u - f  (non-zero only at constrained DOFs)
R_all  = K * u - f;
R_base = R_all(dof(1));
fprintf('\nBase reactions (node 1):\n');
fprintf('  Horizontal force:  Rx = %+.1f N\n',   R_base(1));
fprintf('  Vertical force:    Ry = %+.1f N   (equilibrium check: should be +%.1f)\n', R_base(2), P);
fprintf('  Bending moment:    Mz = %+.1f N·m\n', R_base(3));

% Internal forces at each beam element
stress_cells = postproc_stress(m, [], u);
elem_names = {'Col lower (1-2)', 'Col upper (2-3)', ...
              'Jib inner (3-4)', 'Jib outer (4-5)', 'Brace   (6-4) '};
fprintf('\nInternal forces  [N = axial (+tension), V = shear, M = moment at each end]:\n');
for e = 1:m.nElems
    s = stress_cells{e};
    fprintf('  %-18s  N=%+7.2f kN  V=%+6.2f kN  M1=%+7.2f kN·m  M2=%+7.2f kN·m\n', ...
            elem_names{e}, s(1)/1e3, s(2)/1e3, s(3)/1e3, s(4)/1e3);
end

% Quick utilisation check (bending stress at extreme fibre, column base)
s_base = stress_cells{1};       % element 1: col lower, forces at node-1 end
M_col  = abs(s_base(3));        % M at node 1 end (base)
N_col  = s_base(1);
c_col  = h_c/2 * 1e-3;         % extreme fibre distance [m]
sig_bend = M_col * c_col / Iz_col;
sig_ax   = abs(N_col) / A_col;
sig_tot  = sig_bend + sig_ax;
fprintf('\nColumn base (node 1) — combined stress:\n');
fprintf('  Bending:  %.1f MPa\n',  sig_bend * 1e-6);
fprintf('  Axial:    %.1f MPa\n',  sig_ax   * 1e-6);
fprintf('  Total:    %.1f MPa  (f_y = 355 MPa for S355)\n', sig_tot * 1e-6);

%% ── Visualisation ────────────────────────────────────────────────────────
figure('Name', 'Jib crane — column mount');
scale = 100;   % deformation magnification for visibility
mesh_plot(m, 'color', [0.6 0.6 0.6], 'node_numbers', true, 'elem_numbers', true);
mesh_plot(m, 'deformed', u, 'scale', scale, 'color', 'b');
title(sprintf('Jib crane column  |  tip disp = %.2f mm  (deformation ×%d)', ...
              abs(u5(2))*1e3, scale));
xlabel('x [m]');  ylabel('y [m]');  axis equal;  grid on;
