function mesh = mesh_frame(nodes, conn, elem_type)
% Create a frame or truss mesh from explicit node/element data.
%   nodes     - [nNodes x 2] node coordinates (x, y)
%   conn      - [nElems x 2] element connectivity (1-based node indices)
%   elem_type - 'beam2d' or 'truss2d'
%
% Each element connects two nodes. Nodes are shared between elements.
% Example (simple L-frame):
%   nodes = [0 0; 0 3; 2 3];
%   conn  = [1 2; 2 3];
%   mesh  = mesh_frame(nodes, conn, 'beam2d');

if ~ismember(elem_type, {'beam2d', 'truss2d'})
    error('mesh_frame: elem_type must be ''beam2d'' or ''truss2d''.');
end
if size(conn, 2) ~= 2
    error('mesh_frame: conn must be [nElems x 2].');
end

mesh.nodes  = nodes;
mesh.conn   = conn;
mesh.type   = elem_type;
mesh.nNodes = size(nodes, 1);
mesh.nElems = size(conn,  1);
mesh.t      = 1;  % unused for frame elements, kept for struct consistency
end
