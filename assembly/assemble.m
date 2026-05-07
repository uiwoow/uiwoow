function [K, f] = assemble(mesh, mat, loads)
% Assemble global sparse stiffness matrix K and force vector f.
%   mesh  - mesh struct (type, nodes, conn, t)
%   mat   - material struct
%         For continuum: mat.E, mat.nu, mat.formulation (or mat.D directly)
%         For beam2d:    mat.E, mat.A, mat.Iz
%         For truss2d:   mat.E, mat.A
%   loads - (optional) load struct from bc_loads.m
%         loads.point_loads  - [nLoads x 3]: [node, local_dof, value]
%         loads.body_force   - [2x1] uniform body force [bx; by]
%
%   K  - [nDof x nDof] sparse global stiffness matrix
%   f  - [nDof x 1]    global force vector
if nargin < 3, loads = []; end

% Delegate mixed-element meshes to a dedicated path and return immediately.
if strcmp(mesh.type, 'mixed')
    [K, f] = assemble_mixed(mesh, loads);
    return;
end

nDof = mesh_ndof(mesh);

% Pre-allocate triplets for sparse assembly
switch mesh.type
    case 'T3',      nEDof = 6;
    case 'Q4',      nEDof = 8;
    case 'truss2d', nEDof = 4;
    case 'beam2d',  nEDof = 6;
    otherwise, error('assemble: unknown element type ''%s''.', mesh.type);
end

ntriplets = mesh.nElems * nEDof^2;
Iidx  = zeros(ntriplets, 1);
Jidx  = zeros(ntriplets, 1);
Kvals = zeros(ntriplets, 1);
f     = zeros(nDof, 1);

% Body force
if ~isempty(loads) && isfield(loads, 'body_force') && ~isempty(loads.body_force)
    bf = loads.body_force;
else
    bf = [0; 0];
end

% Material setup
if ismember(mesh.type, {'T3', 'Q4'})
    if isfield(mat, 'D')
        D = mat.D;
    else
        D = mat_elastic_D(mat.E, mat.nu, mat.formulation);
    end
    t = mesh.t;
end

ptr = 0;
for e = 1:mesh.nElems
    nodes_e = mesh.conn(e, :);
    xy_e    = mesh.nodes(nodes_e, :);
    dofs_e  = elem_dofs(mesh, e);

    switch mesh.type
        case 'T3'
            [Ke, fe] = elem_T3(xy_e, D, t, bf);
        case 'Q4'
            [Ke, fe] = elem_Q4(xy_e, D, t, bf);
        case 'truss2d'
            Ke = elem_truss2d(xy_e, mat.E, mat.A);
            fe = zeros(4, 1);
        case 'beam2d'
            Ke = elem_beam2d(xy_e, mat.E, mat.A, mat.Iz);
            fe = zeros(6, 1);
    end

    % Scatter Ke into triplets
    [ii, jj] = ndgrid(dofs_e, dofs_e);
    idx = ptr + (1:nEDof^2);
    Iidx(idx)  = ii(:);
    Jidx(idx)  = jj(:);
    Kvals(idx) = Ke(:);
    ptr = ptr + nEDof^2;

    % Scatter fe into f
    f(dofs_e) = f(dofs_e) + fe;
end

K = sparse(Iidx, Jidx, Kvals, nDof, nDof);

% Apply point loads
if ~isempty(loads) && isfield(loads, 'point_loads') && ~isempty(loads.point_loads)
    pl = loads.point_loads;  % [nLoads x 3]: [node, local_dof, value]
    for i = 1:size(pl, 1)
        node      = pl(i, 1);
        local_dof = pl(i, 2);
        val       = pl(i, 3);
        switch mesh.type
            case {'T3', 'Q4', 'truss2d'}
                gdof = 2*(node-1) + local_dof;
            case 'beam2d'
                gdof = 3*(node-1) + local_dof;
        end
        f(gdof) = f(gdof) + val;
    end
end
end

% -----------------------------------------------------------------------
function [K, f] = assemble_mixed(mesh, loads)
% Sparse assembly for a mixed-element mesh (universal 3-DOF/node system).

nDof = 3 * mesh.nNodes;

% Pre-count triplets so we can pre-allocate without growing arrays.
nedof_of = @(t) deal_nedof(t);
ntriplets = 0;
for e = 1:mesh.nElems
    ne = nedof_of(mesh.elem_types{e});
    ntriplets = ntriplets + ne^2;
end

Iidx  = zeros(ntriplets, 1);
Jidx  = zeros(ntriplets, 1);
Kvals = zeros(ntriplets, 1);
f     = zeros(nDof, 1);

% Body force (applied only to continuum elements)
if ~isempty(loads) && isfield(loads, 'body_force') && ~isempty(loads.body_force)
    bf = loads.body_force;
else
    bf = [0; 0];
end

ptr = 0;
for e = 1:mesh.nElems
    etype   = mesh.elem_types{e};
    nodes_e = mesh.elem_conn{e};
    xy_e    = mesh.nodes(nodes_e, :);
    mat_e   = mesh.materials{mesh.elem_mat(e)};
    dofs_e  = elem_dofs(mesh, e);
    nEDof   = numel(dofs_e);

    switch etype
        case 'T3'
            D = get_D(mat_e);
            [Ke, fe] = elem_T3(xy_e, D, mesh.t, bf);
        case 'Q4'
            D = get_D(mat_e);
            [Ke, fe] = elem_Q4(xy_e, D, mesh.t, bf);
        case 'truss2d'
            Ke = elem_truss2d(xy_e, mat_e.E, mat_e.A);
            fe = zeros(4, 1);
        case 'beam2d'
            Ke = elem_beam2d(xy_e, mat_e.E, mat_e.A, mat_e.Iz);
            fe = zeros(6, 1);
        otherwise
            error('assemble_mixed: unknown element type ''%s''.', etype);
    end

    [ii, jj] = ndgrid(dofs_e, dofs_e);
    idx = ptr + (1:nEDof^2);
    Iidx(idx)  = ii(:);
    Jidx(idx)  = jj(:);
    Kvals(idx) = Ke(:);
    ptr = ptr + nEDof^2;

    f(dofs_e) = f(dofs_e) + fe;
end

K = sparse(Iidx, Jidx, Kvals, nDof, nDof);

% Point loads: local_dof in 3-DOF convention (1=u, 2=v, 3=theta)
if ~isempty(loads) && isfield(loads, 'point_loads') && ~isempty(loads.point_loads)
    pl = loads.point_loads;
    for i = 1:size(pl, 1)
        gdof = 3*(pl(i,1)-1) + pl(i,2);
        f(gdof) = f(gdof) + pl(i,3);
    end
end
end

function ne = deal_nedof(etype)
switch etype
    case 'T3',      ne = 6;
    case 'Q4',      ne = 8;
    case 'truss2d', ne = 4;
    case 'beam2d',  ne = 6;
    otherwise, error('deal_nedof: unknown type ''%s''.', etype);
end
end

function D = get_D(mat_e)
if isfield(mat_e, 'D')
    D = mat_e.D;
else
    D = mat_elastic_D(mat_e.E, mat_e.nu, mat_e.formulation);
end
end
