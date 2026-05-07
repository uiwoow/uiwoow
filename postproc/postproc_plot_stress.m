function postproc_plot_stress(mesh, stress, component, label)
% Filled contour/patch plot of a stress component over the mesh.
%   mesh      - mesh struct
%   stress    - [nElems x nComp] stress array from postproc_stress
%   component - column index to plot (1=sxx, 2=syy, 3=txy for continuum)
%               or 1=N, 2=V, 3=M1, 4=M2 for beam
%   label     - string label for colorbar (optional)
if nargin < 4, label = sprintf('Component %d', component); end

if strcmp(mesh.type, 'mixed')
    error(['postproc_plot_stress: use postproc_plot_mixed() for mixed-element meshes.\n' ...
           'It renders continuum von Mises and beam forces in a single combined figure.']);
end

vals = stress(:, component);

figure;
hold on;

switch mesh.type
    case {'T3', 'Q4'}
        % Use patch with flat shading (one colour per element)
        for e = 1:mesh.nElems
            nd = mesh.conn(e, :);
            xp = mesh.nodes(nd, 1);
            yp = mesh.nodes(nd, 2);
            patch(xp, yp, vals(e), 'EdgeColor', 'none');
        end
        colormap(jet(256));
        colorbar;
        clim([min(vals) max(vals)]);

    case {'beam2d', 'truss2d'}
        % Draw coloured lines with thickness proportional to value
        cmap = jet(256);
        vmin = min(vals); vmax = max(vals);
        rng  = vmax - vmin;
        if rng < eps, rng = 1; end
        for e = 1:mesh.nElems
            nd  = mesh.conn(e, :);
            xp  = mesh.nodes(nd, 1);
            yp  = mesh.nodes(nd, 2);
            idx = round(1 + 255 * (vals(e) - vmin) / rng);
            idx = max(1, min(256, idx));
            plot(xp, yp, 'Color', cmap(idx,:), 'LineWidth', 4);
        end
        colormap(jet(256));
        cb = colorbar;
        clim([vmin vmax]);
        ylabel(cb, label);
end

axis equal; box on; grid on;
xlabel('x'); ylabel('y');
title(label);
hold off;
end
