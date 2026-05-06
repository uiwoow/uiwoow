function [Ke, fe] = elem_Q4(xy, D, t, body_force)
% Q4 bilinear quadrilateral element stiffness matrix (2x2 Gauss).
%   xy          - [4x2] node coordinates (CCW: SW, SE, NE, NW)
%   D           - [3x3] constitutive matrix
%   t           - thickness
%   body_force  - [2x1] [bx; by] body force per unit volume (optional)
%
%   Ke - [8x8] element stiffness matrix
%   fe - [8x1] consistent body-force vector
if nargin < 4, body_force = [0; 0]; end

[xi_g, eta_g, w_g] = gauss_quad(2);
nGP = length(w_g);

Ke = zeros(8, 8);
fe = zeros(8, 1);

for g = 1:nGP
    [B, detJ, N] = elem_Q4_B(xy, xi_g(g), eta_g(g));
    fac = detJ * t * w_g(g);
    Ke  = Ke + fac * (B' * D * B);
    % Consistent body force: N expanded to [2x8] for 2 dofs/node
    Nmat = [N(1) 0 N(2) 0 N(3) 0 N(4) 0;
             0 N(1)  0 N(2)  0 N(3)  0 N(4)];
    fe = fe + fac * (Nmat' * body_force);
end
end

function [B, detJ, N] = elem_Q4_B(xy, xi, eta)
% B matrix [3x8] and detJ at isoparametric coordinates (xi, eta).
dN_dxi  = [-(1-eta)  (1-eta)  (1+eta) -(1+eta)] / 4;
dN_deta = [-(1-xi)  -(1+xi)   (1+xi)   (1-xi) ] / 4;
N       = [(1-xi)*(1-eta)  (1+xi)*(1-eta)  (1+xi)*(1+eta)  (1-xi)*(1+eta)] / 4;

J    = [dN_dxi; dN_deta] * xy;  % 2x2 Jacobian
detJ = det(J);
if detJ <= 0
    error('elem_Q4: non-positive Jacobian (distorted element or wrong node ordering).');
end
invJ = inv(J);

dN_dx = invJ(1,1)*dN_dxi + invJ(1,2)*dN_deta;  % 1x4
dN_dy = invJ(2,1)*dN_dxi + invJ(2,2)*dN_deta;  % 1x4

B = zeros(3, 8);
B(1, 1:2:end) = dN_dx;
B(2, 2:2:end) = dN_dy;
B(3, 1:2:end) = dN_dy;
B(3, 2:2:end) = dN_dx;
end
