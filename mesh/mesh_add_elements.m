function mesh = mesh_add_elements(mesh, elem_type, conn, mat)
% Append a group of same-type elements to a mixed mesh.
%   mesh      - mixed mesh struct (from mesh_mixed)
%   elem_type - string: 'T3' | 'Q4' | 'truss2d' | 'beam2d'
%   conn      - [nNew x nNodesPerElem] connectivity (1-based node indices)
%   mat       - material struct for these elements
%
% Each call may use a different material; the material is stored in
% mesh.materials and referenced by index.  Call mesh_add_elements
% multiple times to mix element types.

if ~strcmp(mesh.type, 'mixed')
    error('mesh_add_elements: mesh.type must be ''mixed''.');
end

valid = {'T3','Q4','truss2d','beam2d'};
if ~any(strcmp(elem_type, valid))
    error('mesh_add_elements: unknown element type ''%s''.', elem_type);
end

% Register material
mat_idx = numel(mesh.materials) + 1;
mesh.materials{mat_idx} = mat;

nNew = size(conn, 1);
for e = 1:nNew
    mesh.nElems = mesh.nElems + 1;
    mesh.elem_types{mesh.nElems} = elem_type;
    mesh.elem_conn{mesh.nElems}  = conn(e, :);
    mesh.elem_mat(mesh.nElems)   = mat_idx;
end
end
