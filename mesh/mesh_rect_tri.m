function mesh = mesh_rect_tri(Lx, Ly, nx, ny)
% Structured rectangular mesh with T3 (CST) triangles.
%   Domain: [0, Lx] x [0, Ly]
%   nx, ny  - number of elements in x and y directions
%             (each quad cell is split into 2 triangles)
%
%   mesh.nodes  - [(nx+1)*(ny+1) x 2]
%   mesh.conn   - [2*nx*ny x 3]
%   mesh.type   - 'T3'

x = linspace(0, Lx, nx+1);
y = linspace(0, Ly, ny+1);
[X, Y] = meshgrid(x, y);
nodes  = [X(:), Y(:)];

nNodes = (nx+1) * (ny+1);
nElems = 2 * nx * ny;
conn   = zeros(nElems, 3);

% meshgrid outputs column-major (x-outer): node at (ix=j, iy=i) -> j*(ny+1)+i+1
nodeIdx = @(i,j) j*(ny+1) + i + 1;

e = 1;
for i = 0:ny-1
    for j = 0:nx-1
        n1 = nodeIdx(i,   j);
        n2 = nodeIdx(i,   j+1);
        n3 = nodeIdx(i+1, j+1);
        n4 = nodeIdx(i+1, j);
        % Split each quad into 2 CCW triangles
        conn(e,   :) = [n1, n2, n3];
        conn(e+1, :) = [n1, n3, n4];
        e = e + 2;
    end
end

mesh.nodes  = nodes;
mesh.conn   = conn;
mesh.type   = 'T3';
mesh.nNodes = nNodes;
mesh.nElems = nElems;
mesh.t      = 1;
end
