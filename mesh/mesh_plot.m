function mesh_plot(mesh, varargin)
% Plot a 2D mesh (original or deformed).
%   mesh - mesh struct (any type, including 'mixed')
%   Optional name-value pairs:
%     'node_numbers', true/false  (default false)
%     'elem_numbers', true/false  (default false)
%     'deformed',     u_vec       (displacement vector to overlay)
%     'scale',        factor      (deformation scale, default 1)
%     'color',        spec        (element edge color, default 'k')

p = inputParser;
addParameter(p, 'node_numbers', false);
addParameter(p, 'elem_numbers', false);
addParameter(p, 'deformed',     []);
addParameter(p, 'scale',        1.0);
addParameter(p, 'color',        'k');
parse(p, varargin{:});
opts = p.Results;

nodes = mesh.nodes;

% --- deformed node positions ---
% DOF layout depends on mesh type:
%   2-DOF/node (T3, Q4, truss2d): u = [u1,v1, u2,v2, ...]
%   3-DOF/node (beam2d, mixed):    u = [u1,v1,th1, u2,v2,th2, ...]
if ~isempty(opts.deformed)
    u = opts.deformed;
    if strcmp(mesh.type, 'beam2d') || strcmp(mesh.type, 'mixed')
        dx = u(1:3:end) * opts.scale;
        dy = u(2:3:end) * opts.scale;
    else
        dx = u(1:2:end) * opts.scale;
        dy = u(2:2:end) * opts.scale;
    end
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
        lw = 2.0 + strcmp(mesh.type, 'beam2d') * 0.5;
        for e = 1:mesh.nElems
            nd = mesh.conn(e, :);
            plot(nodes_def(nd, 1), nodes_def(nd, 2), opts.color, 'LineWidth', lw);
        end

    case 'mixed'
        for e = 1:mesh.nElems
            etype = mesh.elem_types{e};
            nd    = mesh.elem_conn{e};
            xp    = nodes_def(nd, 1);
            yp    = nodes_def(nd, 2);
            switch etype
                case {'T3', 'Q4'}
                    % Close the polygon
                    patch([xp; xp(1)], [yp; yp(1)], 'w', ...
                          'EdgeColor', opts.color, 'LineWidth', 0.5, ...
                          'FaceAlpha', 0);
                case 'beam2d'
                    plot(xp, yp, opts.color, 'LineWidth', 2.5);
                case 'truss2d'
                    plot(xp, yp, opts.color, 'LineWidth', 1.5, 'LineStyle', '--');
            end
        end
end

% --- node labels ---
if opts.node_numbers
    for n = 1:mesh.nNodes
        text(nodes(n,1), nodes(n,2), sprintf(' %d', n), ...
             'FontSize', 7, 'Color', 'b');
    end
end

% --- element labels ---
if opts.elem_numbers
    if strcmp(mesh.type, 'mixed')
        for e = 1:mesh.nElems
            nd = mesh.elem_conn{e};
            cx = mean(mesh.nodes(nd,1));
            cy = mean(mesh.nodes(nd,2));
            text(cx, cy, sprintf('%d', e), 'FontSize', 6, ...
                 'HorizontalAlignment', 'center', 'Color', 'r');
        end
    else
        for e = 1:mesh.nElems
            nd = mesh.conn(e, :);
            cx = mean(nodes(nd,1));
            cy = mean(nodes(nd,2));
            text(cx, cy, sprintf('%d', e), 'FontSize', 7, ...
                 'HorizontalAlignment', 'center', 'Color', 'r');
        end
    end
end

axis equal; box on; grid on;
xlabel('x'); ylabel('y');
title(sprintf('%s mesh  (%d nodes, %d elements)', ...
      mesh.type, mesh.nNodes, mesh.nElems));
hold off;
end
