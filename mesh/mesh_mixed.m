function mesh = mesh_mixed(nodes)
% Create an empty mixed-element mesh from a shared node list.
%   nodes  - [nNodes x 2] nodal (x,y) coordinates
%
% Use mesh_add_elements to populate elements.
% In a mixed mesh, all nodes carry 3 DOFs [u, v, theta].
% Continuum/truss elements only use [u,v]; beam elements use all three.
% Unstiffened theta DOFs at continuum-only nodes are auto-constrained
% by solve_linear before factorisation.
%
% Example:
%   mesh = mesh_mixed(nodes);
%   mesh = mesh_add_elements(mesh, 'Q4',     q4_conn,   mat_plate);
%   mesh = mesh_add_elements(mesh, 'beam2d', beam_conn, mat_beam);

mesh.type      = 'mixed';
mesh.nodes     = nodes;
mesh.nNodes    = size(nodes, 1);
mesh.elem_types = {};   % {nElems x 1} cell of type strings
mesh.elem_conn  = {};   % {nElems x 1} cell of [1 x nNodesPerElem]
mesh.elem_mat   = [];   % [nElems x 1] index into mesh.materials
mesh.materials  = {};   % {nMat x 1} cell of material structs
mesh.nElems     = 0;
mesh.t          = 1;    % thickness for continuum elements
end
