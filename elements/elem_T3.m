function [Ke, fe, B, A] = elem_T3(xy, D, t, body_force)
% T3 Constant Strain Triangle element stiffness matrix.
%   xy          - [3x2] node coordinates [x1 y1; x2 y2; x3 y3]
%   D           - [3x3] constitutive matrix from mat_elastic_D
%   t           - thickness (scalar)
%   body_force  - [2x1] [bx; by] body force per unit volume (optional)
%
%   Ke - [6x6] element stiffness matrix
%   fe - [6x1] consistent body-force vector
%   B  - [3x6] strain-displacement matrix (constant for T3)
%   A  - element area (scalar)
if nargin < 4, body_force = [0; 0]; end

x = xy(:,1); y = xy(:,2);

% Area (positive for CCW node ordering)
A = 0.5 * det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
if A <= 0
    error('elem_T3: element has non-positive area (check node ordering, must be CCW).');
end

b1 = y(2) - y(3);  b2 = y(3) - y(1);  b3 = y(1) - y(2);
c1 = x(3) - x(2);  c2 = x(1) - x(3);  c3 = x(2) - x(1);

B = (1/(2*A)) * [b1  0  b2  0  b3  0 ;
                  0  c1   0  c2   0  c3;
                 c1  b1  c2  b2  c3  b3];

Ke = t * A * (B' * D * B);

% Consistent body-force vector: fe_i = integral(N_i * b) dV = (A*t/3)*b per node
fe = (A * t / 3) * [body_force; body_force; body_force];
end
