function postproc_plot_mixed(mesh, u, stress_cells, scale)
% Visualise a mixed-element model in a single figure.
%
%   Top panel  : deformed shape (continuum elements grey, beams/trusses blue)
%   Bottom panel: von Mises stress on continuum; axial force on beams/trusses
%
%   mesh         - mixed mesh struct
%   u            - [3*nNodes x 1] displacement vector
%   stress_cells - cell array from postproc_stress (mixed output)
%   scale        - (optional) deformation amplification; auto if omitted
%
% For continuum elements, von Mises is computed from [sxx,syy,txy].
% For beam2d elements, the plotted scalar is shear force V (column 2).
% For truss2d elements, the plotted scalar is axial force N (column 1).

if ~strcmp(mesh.type, 'mixed')
    error('postproc_plot_mixed: mesh must be of type ''mixed''.');
end

% --- auto scale ---
if nargin < 4 || isempty(scale)
    u_trans = [u(1:3:end); u(2:3:end)];
    max_d = max(abs(u_trans));
    if max_d > 0
        domain = max(max(mesh.nodes) - min(mesh.nodes));
        scale = 0.1 * domain / max_d;
    else
        scale = 1;
    end
    fprintf('postproc_plot_mixed: auto scale = %.4g\n', scale);
end

% --- deformed node positions (3-DOF layout: u at 3n-2, v at 3n-1) ---
dx = u(1:3:end) * scale;
dy = u(2:3:end) * scale;
nodes_def = mesh.nodes + [dx(:), dy(:)];

% --- separate elements by type ---
idx_cont  = [];   % T3 / Q4
idx_beam  = [];   % beam2d
idx_truss = [];   % truss2d
for e = 1:mesh.nElems
    switch mesh.elem_types{e}
        case {'T3','Q4'},    idx_cont(end+1)  = e; %#ok<AGROW>
        case 'beam2d',       idx_beam(end+1)  = e; %#ok<AGROW>
        case 'truss2d',      idx_truss(end+1) = e; %#ok<AGROW>
    end
end

% --- compute scalar fields ---
% Continuum: von Mises
vm_cont = zeros(numel(idx_cont), 1);
for k = 1:numel(idx_cont)
    s = stress_cells{idx_cont(k)};
    sxx = s(1); syy = s(2); txy = s(3);
    vm_cont(k) = sqrt(sxx^2 - sxx*syy + syy^2 + 3*txy^2);
end

% Beam: shear V (col 2)
v_beam = zeros(numel(idx_beam), 1);
for k = 1:numel(idx_beam)
    v_beam(k) = stress_cells{idx_beam(k)}(2);
end

% Truss: axial N (col 1)
n_truss = zeros(numel(idx_truss), 1);
for k = 1:numel(idx_truss)
    n_truss(k) = stress_cells{idx_truss(k)}(1);
end

% --- figure layout ---
figure('Name', 'Mixed model visualisation', 'Position', [100 100 1000 700]);

%% -- Panel 1: Deformed shape --
subplot(1,2,1);
hold on;
% Continuum: undeformed (ghost) + deformed
for k = 1:numel(idx_cont)
    nd = mesh.elem_conn{idx_cont(k)};
    patch(mesh.nodes(nd,1), mesh.nodes(nd,2), 'w', ...
          'EdgeColor', [0.75 0.75 0.75], 'LineWidth', 0.5, 'FaceAlpha', 0);
    patch(nodes_def(nd,1), nodes_def(nd,2), 'w', ...
          'EdgeColor', 'k', 'LineWidth', 0.8, 'FaceAlpha', 0);
end
% Beams: deformed
for k = 1:numel(idx_beam)
    nd = mesh.elem_conn{idx_beam(k)};
    plot(nodes_def(nd,1), nodes_def(nd,2), 'b-', 'LineWidth', 2.5);
end
% Trusses: deformed
for k = 1:numel(idx_truss)
    nd = mesh.elem_conn{idx_truss(k)};
    plot(nodes_def(nd,1), nodes_def(nd,2), 'r-', 'LineWidth', 2);
end
axis equal; box on; grid on;
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Deformed shape  (\\times%.4g)', scale));
hold off;

%% -- Panel 2: Stress / force field --
subplot(1,2,2);
hold on;

% Continuum: von Mises filled patches
if ~isempty(vm_cont)
    vm_min = min(vm_cont);  vm_max = max(vm_cont);
    vm_range = vm_max - vm_min;
    if vm_range < eps, vm_range = 1; end
    cmap = jet(256);
    for k = 1:numel(idx_cont)
        nd  = mesh.elem_conn{idx_cont(k)};
        xp  = mesh.nodes(nd,1);
        yp  = mesh.nodes(nd,2);
        ci  = round(1 + 255*(vm_cont(k)-vm_min)/vm_range);
        ci  = max(1, min(256, ci));
        patch(xp, yp, vm_cont(k), ...
              'FaceColor', cmap(ci,:), 'EdgeColor', 'k', 'LineWidth', 0.3);
    end
    colormap(jet(256));
    cb = colorbar;
    clim([vm_min, vm_max]);
    ylabel(cb, 'von Mises stress [Pa]');
end

% Beams: coloured lines by shear V
if ~isempty(v_beam)
    vmin = min(v_beam); vmax = max(v_beam);
    vrng = vmax - vmin;  if vrng < eps, vrng = 1; end
    cmap = cool(256);
    for k = 1:numel(idx_beam)
        nd  = mesh.elem_conn{idx_beam(k)};
        xp  = mesh.nodes(nd,1);
        yp  = mesh.nodes(nd,2);
        ci  = round(1 + 255*(v_beam(k)-vmin)/vrng);
        ci  = max(1, min(256, ci));
        plot(xp, yp, 'Color', cmap(ci,:), 'LineWidth', 4);
    end
end

% Trusses: coloured lines by axial N
if ~isempty(n_truss)
    nmin = min(n_truss); nmax = max(n_truss);
    nrng = nmax - nmin;  if nrng < eps, nrng = 1; end
    for k = 1:numel(idx_truss)
        nd  = mesh.elem_conn{idx_truss(k)};
        xp  = mesh.nodes(nd,1);
        yp  = mesh.nodes(nd,2);
        ci  = round(1 + 255*(n_truss(k)-nmin)/nrng);
        ci  = max(1, min(256, ci));
        plot(xp, yp, 'Color', cool(256), 'LineWidth', 3.5);
    end
end

axis equal; box on; grid on;
xlabel('x [m]'); ylabel('y [m]');
title('von Mises [Pa]  |  Beam shear V [N]');
hold off;
end
