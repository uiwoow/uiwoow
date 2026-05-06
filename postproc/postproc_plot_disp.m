function postproc_plot_disp(mesh, u, scale)
% Plot deformed shape overlaid on original mesh.
%   mesh  - mesh struct
%   u     - [nDof x 1] displacement vector
%   scale - deformation amplification factor (default: auto)
if nargin < 3 || isempty(scale)
    max_disp = max(abs(u));
    if max_disp > 0
        domain_size = max(max(mesh.nodes) - min(mesh.nodes));
        scale = 0.1 * domain_size / max_disp;
    else
        scale = 1;
    end
    fprintf('postproc_plot_disp: auto scale = %.4g\n', scale);
end

figure;
hold on;
% Original mesh (grey dashed)
mesh_plot(mesh, 'color', [0.7 0.7 0.7]);
% Deformed mesh (black)
mesh_plot(mesh, 'deformed', u, 'scale', scale, 'color', 'k');

title(sprintf('Deformed shape (scale x%.4g)', scale));
end
