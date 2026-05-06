% 2D Jib crane analysis — beam2d elements.
%
% Geometry (wall-mounted jib crane):
%
%     2---------3   <- jib arm at height H
%    /|
%   / |            Node 1: wall mount (bottom, pinned)
%  /  |            Node 2: wall mount (top, pinned)
% 1   |            Node 3: jib tip (tip load P downward)
%   (wall)
%
%  Member 1-3: jib arm (horizontal boom, from wall to tip)
%  Member 2-3: tie rod (diagonal brace from upper mount to tip)
%
%  All members are beam elements (carry axial + bending).
%  For a pure truss analysis, replace 'beam2d' with 'truss2d'.
%
% Parameters
% ----------
%   Wall mount height:  H = 2 m (node 2 above node 1)
%   Jib arm length:     Lj = 4 m (horizontal, from wall to tip)
%   Tip load:           P  = 5000 N (downward)

startup_fem2d;

% ---- Geometry ----
H  = 2;     % height between wall mounts [m]
Lj = 4;     % jib arm length [m]

% ---- Material (steel I-beam) ----
E  = 2.1e11;    % Young's modulus [Pa]
A  = 6.0e-3;    % cross-section area [m^2]
Iz = 2.0e-5;    % second moment of area [m^4]

% ---- Load ----
P  = -5000;     % tip load [N], negative = downward

% ---- Nodes ----
%  1: lower wall mount (0, 0)
%  2: upper wall mount (0, H)
%  3: jib tip          (Lj, H)  <- arm and tie rod meet here
nodes = [0   0;     % node 1
         0   H;     % node 2
         Lj  H];    % node 3

% ---- Elements ----
%  1: jib arm  (1 -> 3, horizontal boom at height H... but node 1 is at y=0)
%     Actually model as: lower mount(1) -> tip(3), this is the inclined arm
%  2: tie rod  (2 -> 3, upper mount to tip)
%
% A more realistic layout: the arm goes from node 1 at y=0 to node 3 at y=H
% (inclined arm), with a horizontal tie at the top connecting node 2 to 3.
conn = [1  3;   % jib arm (inclined)
        2  3];  % tie rod (horizontal)

mesh = mesh_frame(nodes, conn, 'beam2d');
mat.E = E; mat.A = A; mat.Iz = Iz;

% ---- Boundary conditions ----
% Both wall mounts fully fixed (welded to wall): u, v, theta = 0
fix_nodes = [1; 1; 1; 2; 2; 2];
fix_dofs  = [1; 2; 3; 1; 2; 3];
bc.fixed_nodes = fix_nodes;
bc.fixed_dofs  = fix_dofs;
bc.fixed_vals  = zeros(6, 1);

% ---- Loads ----
% Tip load P at node 3, dof 2 (vertical)
loads.point_loads = [3, 2, P];

% ---- Solve ----
[K, f] = assemble(mesh, mat, loads);
u      = solve_linear(K, f, mesh, bc);

% ---- Results ----
[forces, ~] = postproc_stress(mesh, mat, u);

% Node 3 displacements (dofs 7,8,9 for beam2d with 3 dofs/node)
ux3 = u(3*3 - 2);
uy3 = u(3*3 - 1);

fprintf('\n=== Jib Crane Results ===\n');
fprintf('Tip deflection:\n');
fprintf('  Horizontal (ux): %.4f mm\n', ux3 * 1e3);
fprintf('  Vertical   (uy): %.4f mm\n', uy3 * 1e3);
fprintf('\nMember internal forces:\n');
member_names = {'Jib arm (1-3)', 'Tie rod (2-3)'};
for e = 1:mesh.nElems
    N  = forces(e,1);
    V  = forces(e,2);
    M1 = forces(e,3);
    M2 = forces(e,4);
    fprintf('  %s:\n', member_names{e});
    fprintf('    Axial N = %8.1f N  (+ tension, - compression)\n', N);
    fprintf('    Shear V = %8.1f N\n', V);
    fprintf('    Moment at node 1 of element = %8.1f N.m\n', M1);
    fprintf('    Moment at node 2 of element = %8.1f N.m\n', M2);
end

% ---- Stress check at critical section ----
% Maximum bending stress at fixed end: sigma = M * c / I  where c = section depth/2
c_section = 0.1;  % half-depth of section [m] — adjust to your section
for e = 1:mesh.nElems
    M_max = max(abs([forces(e,3), forces(e,4)]));
    sig_bending = M_max * c_section / Iz;
    sig_axial   = forces(e,1) / A;
    sig_total   = abs(sig_axial) + sig_bending;  % conservative combination
    fprintf('\n  %s: max combined stress = %.2f MPa\n', ...
            member_names{e}, sig_total * 1e-6);
end

% ---- Plots ----
figure;
subplot(1,2,1);
mesh_plot(mesh, 'node_numbers', true);
title('Jib crane — original geometry');

subplot(1,2,2);
scale = 200;
mesh_plot(mesh, 'deformed', u, 'scale', scale, 'color', 'b');
mesh_plot(mesh, 'color', [0.7 0.7 0.7]);
title(sprintf('Deformed shape (x%d)', scale));
