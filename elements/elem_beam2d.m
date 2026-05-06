function Ke = elem_beam2d(xy, E, A, Iz)
% 2D Euler-Bernoulli beam element stiffness in global coordinates.
%   xy  - [2x2] node coordinates [x1 y1; x2 y2]
%   E   - Young's modulus
%   A   - cross-sectional area
%   Iz  - second moment of area about z-axis (bending)
%
%   Ke  - [6x6] global stiffness matrix
%   DOF order per node: [u, v, theta] (axial, transverse, rotation)
%   Full element DOF vector: [u1 v1 th1 u2 v2 th2]
dx = xy(2,1) - xy(1,1);
dy = xy(2,2) - xy(1,2);
L  = sqrt(dx^2 + dy^2);
if L < eps
    error('elem_beam2d: zero-length element.');
end

c = dx / L;  s = dy / L;

% Axial stiffness terms
EA_L = E * A / L;

% Bending stiffness terms (Euler-Bernoulli)
EI  = E * Iz;
k1  = 12*EI / L^3;
k2  =  6*EI / L^2;
k3  =  4*EI / L;
k4  =  2*EI / L;

% Local stiffness [6x6]: dof = [u1 v1 th1 u2 v2 th2]
Ke_local = [ EA_L   0     0    -EA_L   0     0   ;
              0      k1    k2    0     -k1    k2  ;
              0      k2    k3    0     -k2    k4  ;
             -EA_L   0     0     EA_L   0     0   ;
              0     -k1   -k2    0      k1   -k2  ;
              0      k2    k4    0     -k2    k3  ];

% Transformation matrix from local to global (rotation about z)
T = [ c   s  0   0   0  0;
     -s   c  0   0   0  0;
      0   0  1   0   0  0;
      0   0  0   c   s  0;
      0   0  0  -s   c  0;
      0   0  0   0   0  1];

Ke = T' * Ke_local * T;
end
