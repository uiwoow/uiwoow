function mesh = mesh_import_gmsh(filename, varargin)
% Import a 2D mesh from a Gmsh v2 ASCII (.msh) file.
%
%   mesh = mesh_import_gmsh(filename)
%   mesh = mesh_import_gmsh(filename, 'elem_type', 'Q4')
%
% Supported element types imported as 2D elements:
%   Gmsh type 2  -> 'T3' (3-node triangle)
%   Gmsh type 3  -> 'Q4' (4-node quadrilateral)
% Line and point elements are parsed for physical-group bookkeeping but
% are not added to mesh.conn.
%
% Options (name-value pairs):
%   'elem_type'  - force import of only 'T3' or 'Q4'; skip the other.
%                  Default: auto-detect (error if both found and not specified).
%   'dim'        - spatial dimensions to keep (default 2; z is discarded).
%
% Output mesh struct (same layout as mesh_rect_quad / mesh_rect_tri):
%   mesh.nodes          - [nNodes x 2]  (x,y) coordinates
%   mesh.conn           - [nElems x nNodesPerElem]
%   mesh.type           - 'T3' or 'Q4'
%   mesh.nNodes, nElems
%   mesh.t              - 1  (default thickness; set by caller)
%   mesh.gmsh_groups    - struct array: .name, .dim, .tag, .node_ids
%                         (physical groups, useful for applying BCs)
%
% The Gmsh v2 ASCII format spec:
%   https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2

% --- parse options ---
p = inputParser;
addParameter(p, 'elem_type', '', @ischar);
parse(p, varargin{:});
force_type = p.Results.elem_type;

% Gmsh element-type -> (number of nodes, our type string or 'skip')
GMSH_ETYPE = containers.Map( ...
    {1,  2,   3,   15}, ...
    {'line', 'T3', 'Q4', 'skip'});

fid = fopen(filename, 'r');
if fid < 0
    error('mesh_import_gmsh: cannot open file ''%s''.', filename);
end
lines = {};
while ~feof(fid)
    ln = strtrim(fgetl(fid));
    if ischar(ln)
        lines{end+1} = ln; %#ok<AGROW>
    end
end
fclose(fid);

% --- locate section boundaries ---
sec = struct();
for i = 1:numel(lines)
    ln = lines{i};
    switch ln
        case '$MeshFormat',      sec.MeshFormat_start   = i;
        case '$EndMeshFormat',   sec.MeshFormat_end     = i;
        case '$PhysicalNames',   sec.PhysNames_start    = i;
        case '$EndPhysicalNames',sec.PhysNames_end      = i;
        case '$Nodes',           sec.Nodes_start        = i;
        case '$EndNodes',        sec.Nodes_end          = i;
        case '$Elements',        sec.Elements_start     = i;
        case '$EndElements',     sec.Elements_end       = i;
    end
end

% --- validate format ---
if isfield(sec, 'MeshFormat_start')
    fmt_line = lines{sec.MeshFormat_start + 1};
    ver = sscanf(fmt_line, '%f', 1);
    if ver < 2 || ver >= 3
        error('mesh_import_gmsh: only Gmsh v2 ASCII format is supported (found v%.1f).', ver);
    end
end

% --- parse physical names (optional) ---
gmsh_groups = struct('name', {}, 'dim', {}, 'tag', {}, 'node_ids', {});
if isfield(sec, 'PhysNames_start')
    n_phys = str2double(lines{sec.PhysNames_start + 1});
    for i = 1:n_phys
        parts = strsplit(lines{sec.PhysNames_start + 1 + i});
        dim   = str2double(parts{1});
        tag   = str2double(parts{2});
        name  = strrep(strjoin(parts(3:end)), '"', '');
        gmsh_groups(end+1).name    = name; %#ok<AGROW>
        gmsh_groups(end).dim      = dim;
        gmsh_groups(end).tag      = tag;
        gmsh_groups(end).node_ids = [];
    end
end

% --- parse nodes ---
n_nodes = str2double(lines{sec.Nodes_start + 1});
nodes_raw = zeros(n_nodes, 4);  % [id x y z]
for i = 1:n_nodes
    nums = sscanf(lines{sec.Nodes_start + 1 + i}, '%f')';
    nodes_raw(i, :) = nums(1:4);
end
% Build index remapping: Gmsh node IDs may not be 1-based contiguous
max_id = max(nodes_raw(:,1));
id2idx = zeros(max_id, 1);
for i = 1:n_nodes
    id2idx(nodes_raw(i,1)) = i;
end
nodes = nodes_raw(:, 2:3);  % keep x,y; drop z

% --- parse elements ---
n_elem_total = str2double(lines{sec.Elements_start + 1});
conn_T3  = [];
conn_Q4  = [];
line_elems = struct('phys_tag', {}, 'conn', {});  % for physical group assignment

for i = 1:n_elem_total
    nums  = sscanf(lines{sec.Elements_start + 1 + i}, '%d')';
    % nums = [elemID, elemType, nTags, tag1, ..., tagN, nodeID1, ...]
    if numel(nums) < 3, continue; end
    etype  = nums(2);
    n_tags = nums(3);
    phys_tag = 0;
    if n_tags >= 1
        phys_tag = nums(4);
    end
    node_ids = nums(3 + n_tags + 1 : end);
    local_ids = id2idx(node_ids);

    if ~isKey(GMSH_ETYPE, etype), continue; end

    switch GMSH_ETYPE(etype)
        case 'T3'
            conn_T3(end+1, :) = local_ids(1:3); %#ok<AGROW>
        case 'Q4'
            conn_Q4(end+1, :) = local_ids(1:4); %#ok<AGROW>
        case 'line'
            line_elems(end+1).phys_tag = phys_tag; %#ok<AGROW>
            line_elems(end).conn = local_ids(1:2);
    end
end

% --- determine element type to use ---
has_T3 = ~isempty(conn_T3);
has_Q4 = ~isempty(conn_Q4);

if ~has_T3 && ~has_Q4
    error('mesh_import_gmsh: no T3 or Q4 elements found in ''%s''.', filename);
end

if ~isempty(force_type)
    switch upper(force_type)
        case 'T3'
            if ~has_T3
                error('mesh_import_gmsh: no T3 elements found (force_type=''T3'').');
            end
            conn = conn_T3;  elem_type = 'T3';
            if has_Q4
                warning('mesh_import_gmsh: %d Q4 elements skipped (force_type=T3).', size(conn_Q4,1));
            end
        case 'Q4'
            if ~has_Q4
                error('mesh_import_gmsh: no Q4 elements found (force_type=''Q4'').');
            end
            conn = conn_Q4;  elem_type = 'Q4';
            if has_T3
                warning('mesh_import_gmsh: %d T3 elements skipped (force_type=Q4).', size(conn_T3,1));
            end
        otherwise
            error('mesh_import_gmsh: force_type must be ''T3'' or ''Q4''.');
    end
elseif has_T3 && has_Q4
    % Both found: prefer the more abundant type; warn about the other.
    if size(conn_Q4,1) >= size(conn_T3,1)
        conn = conn_Q4;  elem_type = 'Q4';
        warning('mesh_import_gmsh: found both T3 (%d) and Q4 (%d) elements. Importing Q4; use ''elem_type'' option to override.', ...
                size(conn_T3,1), size(conn_Q4,1));
    else
        conn = conn_T3;  elem_type = 'T3';
        warning('mesh_import_gmsh: found both T3 (%d) and Q4 (%d) elements. Importing T3; use ''elem_type'' option to override.', ...
                size(conn_T3,1), size(conn_Q4,1));
    end
elseif has_T3
    conn = conn_T3;  elem_type = 'T3';
else
    conn = conn_Q4;  elem_type = 'Q4';
end

% --- remove unused nodes (keep only those referenced by selected elements) ---
used_nodes = unique(conn(:));
% Remap to contiguous indices
new_id = zeros(n_nodes, 1);
new_id(used_nodes) = 1:numel(used_nodes);
nodes  = nodes(used_nodes, :);
conn   = new_id(conn);

% --- populate physical group node lists ---
% Associate each physical group with the nodes from its line elements
for k = 1:numel(line_elems)
    pt = line_elems(k).phys_tag;
    for g = 1:numel(gmsh_groups)
        if gmsh_groups(g).tag == pt
            % remap line element nodes to new contiguous IDs
            ln = id2idx(line_elems(k).conn);
            mapped = new_id(ln);
            mapped = mapped(mapped > 0);  % keep only nodes that survived pruning
            gmsh_groups(g).node_ids = unique([gmsh_groups(g).node_ids(:); mapped(:)]);
            break;
        end
    end
end

% --- build output mesh struct ---
mesh.type        = elem_type;
mesh.nodes       = nodes;
mesh.conn        = conn;
mesh.nNodes      = size(nodes,  1);
mesh.nElems      = size(conn,   1);
mesh.t           = 1;
mesh.gmsh_groups = gmsh_groups;

fprintf('mesh_import_gmsh: imported %d %s elements, %d nodes from ''%s''.\n', ...
        mesh.nElems, mesh.type, mesh.nNodes, filename);
if ~isempty(gmsh_groups)
    fprintf('  Physical groups: ');
    for g = 1:numel(gmsh_groups)
        fprintf('''%s'' (%d nodes)  ', gmsh_groups(g).name, numel(gmsh_groups(g).node_ids));
    end
    fprintf('\n');
end
end
