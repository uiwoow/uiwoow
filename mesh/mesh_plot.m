function mesh_plot(mesh, varargin)
% Plot a 2D mesh.
%   mesh - mesh struct
%   Optional name-value pairs:
%     'node_numbers', true/false  (default false)
%     'elem_numbers', true/false  (default false)
%     'deformed',     u_vec       (displacement vector to overlay)
%     'scale',        factor      (deformation scale, default 1)
%     'color',        'b'         (element edge color, default 'k')

p = inputParser;
addOptional(p, 'node_numbers', false);
addOptional(p, 'elem_numbers', false);
addOptional(p, 'deformed',     []);
addOptional(p, 'scale',        1.0);
addOptional(p, 'color',        'k');
parse(p, varargin{:});
opts = p.Results;

nodes = mesh.nodes;
if ~isempty(opts.deformed)
    u = opts.deformed;
    dx = u(1:2:end) * opts.scale;
    dy = u(2:2:end) * opts.scale;
    nodes_def = nodes + [dx(:), dy(:)];
else
    nodes_def = nodes;
end

hold on;
switch mesh.type
    case {'T3', 'Q4'}
        for e = 1:mesh.nElems
            nd = mesh.conn(e, :);
            xp = [nodes_def(nd, 1); nodes_def(nd(1), 1)];
            yp = [nodes_def(nd, 2); nodes_def(nd(1), 2)];
            plot(xp, yp, opts.color, 'LineWidth', 0.5);
        end
    case {'beam2d', 'truss2d'}
        for e = 1:mesh.nElems
            nd = mesh.conn(e, :);
            xp = nodes_def(nd, 1);
            yp = nodes_def(nd, 2);
            plot(xp, yp, opts.color, 'LineWidth', 2);
        end
end

if opts.node_numbers
    for n = 1:mesh.nNodes
        text(nodes(n,1), nodes(n,2), sprintf(' %d', n), ...
             'FontSize', 8, 'Color', 'b');
    end
end
if opts.elem_numbers
    for e = 1:mesh.nElems
        nd  = mesh.conn(e, :);
        cx  = mean(nodes(nd, 1));
        cy  = mean(nodes(nd, 2));
        text(cx, cy, sprintf('%d', e), 'FontSize', 7, ...
             'HorizontalAlignment', 'center', 'Color', 'r');
    end
end

axis equal; box on; grid on;
xlabel('x'); ylabel('y');
title(sprintf('%s mesh (%d nodes, %d elements)', mesh.type, mesh.nNodes, mesh.nElems));
hold off;
end
