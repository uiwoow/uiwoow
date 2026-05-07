function bc = bc_auto_orphan_rot(mesh, bc)
% For a mixed mesh, auto-constrain theta=0 at nodes that have no beam
% element attached (continuum and truss elements carry no rotational
% stiffness, leaving those DOFs unstiffened and the system singular).
%
% Only DOFs not already present in bc are added, so calling this
% function on a fully-clamped wall (which already fixes theta) is safe.
%
% Called automatically by solve_linear for mixed meshes.

if ~strcmp(mesh.type, 'mixed')
    return;
end

% Collect nodes touched by at least one beam2d element
beam_nodes = [];
for e = 1:mesh.nElems
    if strcmp(mesh.elem_types{e}, 'beam2d')
        beam_nodes = [beam_nodes, mesh.elem_conn{e}]; %#ok<AGROW>
    end
end
beam_nodes = unique(beam_nodes);

% Nodes that have no beam element attached
non_beam = setdiff(1:mesh.nNodes, beam_nodes);
if isempty(non_beam)
    return;
end

% Build set of global DOFs already constrained for theta (local dof = 3)
existing_theta_nodes = bc.fixed_nodes(bc.fixed_dofs == 3);

% Restrict to non-beam nodes not yet constrained
to_constrain = setdiff(non_beam(:), existing_theta_nodes(:));
if isempty(to_constrain)
    return;
end

n_new = numel(to_constrain);
bc.fixed_nodes = [bc.fixed_nodes(:); to_constrain(:)];
bc.fixed_dofs  = [bc.fixed_dofs(:);  3*ones(n_new, 1)];
bc.fixed_vals  = [bc.fixed_vals(:);  zeros(n_new, 1)];
end
