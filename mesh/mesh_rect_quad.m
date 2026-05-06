function mesh = mesh_rect_quad(Lx, Ly, nx, ny)
% Structured rectangular mesh with Q4 (bilinear quad) elements.
%   Domain: [0, Lx] x [0, Ly]
%   nx, ny  - number of elements in x and y directions
%
%   mesh.nodes  - [(nx+1)*(ny+1) x 2]
%   mesh.conn   - [nx*ny x 4]  (CCW node ordering: SW, SE, NE, NW)
%   mesh.type   - 'Q4'

x = linspace(0, Lx, nx+1);
y = linspace(0, Ly, ny+1);
[X, Y] = meshgrid(x, y);
nodes  = [X(:), Y(:)];

nNodes = (nx+1) * (ny+1);
nElems = nx * ny;
conn   = zeros(nElems, 4);

nodeIdx = @(i,j) i*(nx+1) + j + 1;  % i=row(0-based), j=col(0-based)

e = 1;
for i = 0:ny-1
    for j = 0:nx-1
        n_sw = nodeIdx(i,   j);
        n_se = nodeIdx(i,   j+1);
        n_ne = nodeIdx(i+1, j+1);
        n_nw = nodeIdx(i+1, j);
        conn(e, :) = [n_sw, n_se, n_ne, n_nw];
        e = e + 1;
    end
end

mesh.nodes  = nodes;
mesh.conn   = conn;
mesh.type   = 'Q4';
mesh.nNodes = nNodes;
mesh.nElems = nElems;
mesh.t      = 1;
end
