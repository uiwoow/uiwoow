function [xi, eta, w] = gauss_quad(n)
% Gauss-Legendre quadrature points and weights on [-1,1]^2.
%   n   - integration order (1, 2, or 3)
%   xi  - [n^2 x 1] xi coordinates
%   eta - [n^2 x 1] eta coordinates
%   w   - [n^2 x 1] weights
switch n
    case 1
        pts = 0; wts = 2;
    case 2
        pts = [-1 1] / sqrt(3); wts = [1 1];
    case 3
        pts = [-sqrt(3/5) 0 sqrt(3/5)]; wts = [5/9 8/9 5/9];
    otherwise
        error('gauss_quad: order must be 1, 2, or 3.');
end
[XI, ETA] = meshgrid(pts, pts);
[WX, WY]  = meshgrid(wts, wts);
xi  = XI(:);
eta = ETA(:);
w   = WX(:) .* WY(:);
end
