function [L1, L2, L3, w] = gauss_tri(order)
% Gauss quadrature points on reference triangle (L1+L2+L3=1).
%   order 1 - 1-point  (exact for linears)
%   order 2 - 3-point  (exact for quadratics)
%   order 3 - 4-point  (exact for cubics)
%   Weights sum to 0.5 (area of reference triangle).
switch order
    case 1
        L1 = 1/3; L2 = 1/3; L3 = 1/3;
        w  = 0.5;
    case 2
        L1 = [1/6; 2/3; 1/6];
        L2 = [1/6; 1/6; 2/3];
        L3 = [2/3; 1/6; 1/6];
        w  = [1/6; 1/6; 1/6];
    case 3
        L1 = [1/3;  1/5;  3/5;  1/5];
        L2 = [1/3;  1/5;  1/5;  3/5];
        L3 = [1/3;  3/5;  1/5;  1/5];
        w  = [-27/96; 25/96; 25/96; 25/96];
    otherwise
        error('gauss_tri: order must be 1, 2, or 3.');
end
end
