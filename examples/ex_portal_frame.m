% Portal frame under horizontal (wind) load — beam2d elements.
%
%  Geometry:
%      2-------3
%      |       |
%      |       |
%      1       4   <- pinned supports
%
%  Column height H, beam span L, horizontal load P at node 2.
startup_fem2d;

H = 3;      % column height [m]
L = 4;      % beam span [m]
E = 2e11;   % steel [Pa]
A = 1e-2;   % cross-section area [m^2]
Iz = 1e-4;  % second moment of area [m^4]
P = 1e4;    % horizontal load [N]

nodes = [0  0;    % 1 - bottom left
         0  H;    % 2 - top left
         L  H;    % 3 - top right
         L  0];   % 4 - bottom right

conn = [1 2;   % left column
        2 3;   % beam
        3 4];  % right column

mesh = mesh_frame(nodes, conn, 'beam2d');
mat.E = E; mat.A = A; mat.Iz = Iz;

% Pinned at nodes 1 and 4 (u and v fixed, rotation free)
bc.fixed_nodes = [1; 1; 4; 4];
bc.fixed_dofs  = [1; 2; 1; 2];
bc.fixed_vals  = zeros(4, 1);

% Horizontal load at node 2, DOF 1 (u direction)
loads.point_loads = [2, 1, P];

[K, f] = assemble(mesh, mat, loads);
u = solve_linear(K, f, mesh, bc);

[forces, ~] = postproc_stress(mesh, mat, u);
N  = forces(:,1);  V  = forces(:,2);
M1 = forces(:,3);  M2 = forces(:,4);

fprintf('Portal Frame Results:\n');
fprintf('  Node 2 horizontal disp: ux = %.4e m\n', u(3*1+1));
fprintf('  Node 2 vertical disp:   uy = %.4e m\n', u(3*1+2));
fprintf('  Element moments (M1, M2) [N.m]:\n');
for e = 1:mesh.nElems
    fprintf('    Elem %d: M1=%.3e  M2=%.3e  N=%.3e  V=%.3e\n', ...
            e, M1(e), M2(e), N(e), V(e));
end

figure;
mesh_plot(mesh, 'deformed', u, 'scale', 50);
title('Portal frame deformed (50x)');

postproc_plot_stress(mesh, forces, 3, 'Bending moment M1 [N.m]');
