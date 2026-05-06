function node_ids = mesh_find_nodes(mesh, tol, varargin)
% Find nodes matching a geometric criterion.
%   mesh     - mesh struct
%   tol      - coordinate tolerance
%   Criteria (one of):
%     'x', x0           - nodes where x == x0
%     'y', y0           - nodes where y == y0
%     'xy', x0, y0      - single node closest to (x0, y0)
%     'box', [x1 x2 y1 y2] - nodes inside bounding box
%
% Returns node_ids as column vector of 1-based node indices.
if isempty(varargin)
    error('mesh_find_nodes: specify a criterion (''x'', ''y'', ''xy'', or ''box'').');
end

x = mesh.nodes(:,1);
y = mesh.nodes(:,2);
criterion = varargin{1};

switch lower(criterion)
    case 'x'
        x0 = varargin{2};
        node_ids = find(abs(x - x0) <= tol);
    case 'y'
        y0 = varargin{2};
        node_ids = find(abs(y - y0) <= tol);
    case 'xy'
        x0 = varargin{2};  y0 = varargin{3};
        d = sqrt((x - x0).^2 + (y - y0).^2);
        node_ids = find(d <= tol);
        if isempty(node_ids)
            [~, node_ids] = min(d);  % fallback: nearest node
        end
    case 'box'
        bb = varargin{2};  % [x1 x2 y1 y2]
        node_ids = find(x >= bb(1)-tol & x <= bb(2)+tol & ...
                        y >= bb(3)-tol & y <= bb(4)+tol);
    otherwise
        error('mesh_find_nodes: unknown criterion ''%s''.', criterion);
end
node_ids = node_ids(:);
end
