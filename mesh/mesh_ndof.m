function n = mesh_ndof(mesh)
% Return total number of DOFs for a mesh.
%   'T3', 'Q4'        -> 2 dofs/node
%   'truss2d'         -> 2 dofs/node
%   'beam2d'          -> 3 dofs/node
switch mesh.type
    case {'T3', 'Q4'}
        n = 2 * mesh.nNodes;
    case 'truss2d'
        n = 2 * mesh.nNodes;
    case 'beam2d'
        n = 3 * mesh.nNodes;
    otherwise
        error('mesh_ndof: unknown element type ''%s''.', mesh.type);
end
end
