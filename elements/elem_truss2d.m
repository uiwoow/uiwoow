function [Ke, dofs_local] = elem_truss2d(xy, E, A)
% 2D truss (bar) element stiffness matrix in global coordinates.
%   xy  - [2x2] node coordinates [x1 y1; x2 y2]
%   E   - Young's modulus
%   A   - cross-sectional area
%
%   Ke          - [4x4] global stiffness matrix
%   dofs_local  - dof layout: [u1 v1 u2 v2] (2 dofs per node)
dx = xy(2,1) - xy(1,1);
dy = xy(2,2) - xy(1,2);
L  = sqrt(dx^2 + dy^2);
if L < eps
    error('elem_truss2d: zero-length element.');
end

c = dx / L;  s = dy / L;

% Local stiffness (axial only): EA/L * [1 -1; -1 1]
% Rotated to global 4x4
T = [c  s  0  0;
    -s  c  0  0;
     0  0  c  s;
     0  0 -s  c];

Ke_local = (E*A/L) * [1  0 -1  0;
                       0  0  0  0;
                      -1  0  1  0;
                       0  0  0  0];

Ke = T' * Ke_local * T;
dofs_local = [];  % layout info: [u1 v1 u2 v2]
end
