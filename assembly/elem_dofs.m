function dofs = elem_dofs(mesh, e)
% Return global DOF indices for element e.
%   Continuum (T3, Q4):  2 dofs/node -> [2n-1, 2n]
%   Truss (truss2d):     2 dofs/node -> [2n-1, 2n]
%   Beam  (beam2d):      3 dofs/node -> [3n-2, 3n-1, 3n]
if strcmp(mesh.type, 'mixed')
    % Per-element type and connectivity stored in cell arrays.
    % Universal 3-DOF/node system; continuum/truss skip the theta slot.
    etype   = mesh.elem_types{e};
    nodes_e = mesh.elem_conn{e};
    switch etype
        case {'T3', 'Q4', 'truss2d'}
            dofs = reshape([3*nodes_e - 2; 3*nodes_e - 1], 1, []);
        case 'beam2d'
            dofs = reshape([3*nodes_e - 2; 3*nodes_e - 1; 3*nodes_e], 1, []);
        otherwise
            error('elem_dofs: unknown mixed element type ''%s''.', etype);
    end
    return;
end

nodes_e = mesh.conn(e, :);

switch mesh.type
    case {'T3', 'Q4', 'truss2d'}
        dofs = reshape([2*nodes_e - 1; 2*nodes_e], 1, []);
    case 'beam2d'
        dofs = reshape([3*nodes_e - 2; 3*nodes_e - 1; 3*nodes_e], 1, []);
    otherwise
        error('elem_dofs: unknown element type ''%s''.', mesh.type);
end
end
