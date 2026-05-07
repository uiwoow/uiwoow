function [K, f] = constraint_rigid_link_apply(K, f, mesh, master_node, slave_nodes, alpha)
% Enforce rigid-link kinematics between a beam master node and a set of
% Q4/T3 slave nodes using the penalty method.
%
% A rigid link transmits BOTH force AND moment from the beam end to the
% plate boundary.  Without it, a single-node beam-to-continuum connection
% is effectively pinned (force only, no moment).
%
% Kinematic constraint per slave node s, master at node m:
%   u_s = u_m - theta_m * (y_s - y_m)        [u constraint]
%   v_s = v_m + theta_m * (x_s - x_m)        [v constraint]
%
% Penalty: alpha * residual^2 added to potential energy.
% Rule of thumb for alpha: 1e6 to 1e8 times the largest diagonal entry
% of K.  Larger alpha = tighter constraint but worse conditioning.
%
% Usage:
%   slave_nodes = mesh_find_nodes(mesh, tol, 'x', x_plate_right);
%   [K, f] = constraint_rigid_link_apply(K, f, mesh, master_node, slave_nodes);
%   u = solve_linear(K, f, mesh, bc);
%
% Note: only valid for mixed meshes (3 DOFs/node: u, v, theta).
%
%   K, f       - in/out: global stiffness and force (sparse K modified)
%   mesh       - mixed mesh struct
%   master_node- scalar node index; must be a beam2d node (has theta DOF)
%   slave_nodes- [n x 1] node indices on the continuum boundary
%   alpha      - (optional) penalty factor; auto-selected if omitted

if ~strcmp(mesh.type, 'mixed')
    error('constraint_rigid_link_apply: only valid for mixed-element meshes.');
end
if isscalar(slave_nodes) && slave_nodes == master_node
    error('constraint_rigid_link_apply: master and slave cannot be the same node.');
end

if nargin < 6 || isempty(alpha)
    d = abs(diag(K));
    d_max = max(d(d > 0));
    if isempty(d_max) || d_max == 0
        d_max = 1;
    end
    alpha = 1e7 * d_max;
end

x_m = mesh.nodes(master_node, 1);
y_m = mesh.nodes(master_node, 2);

dof_um  = 3*(master_node-1) + 1;
dof_vm  = 3*(master_node-1) + 2;
dof_thm = 3*(master_node-1) + 3;

for k = 1:numel(slave_nodes)
    ns = slave_nodes(k);
    if ns == master_node, continue; end

    rx = mesh.nodes(ns,1) - x_m;   % x_s - x_m
    ry = mesh.nodes(ns,2) - y_m;   % y_s - y_m

    dof_us = 3*(ns-1) + 1;
    dof_vs = 3*(ns-1) + 2;

    % u constraint: u_s - u_m + theta_m*ry = 0
    %   c^T = [1, -1, ry]  mapped to dofs [dof_us, dof_um, dof_thm]
    du = [dof_us; dof_um; dof_thm];
    cu = [1; -1; ry];
    K(du, du) = K(du, du) + alpha * (cu * cu');

    % v constraint: v_s - v_m - theta_m*rx = 0
    %   c^T = [1, -1, -rx]  mapped to dofs [dof_vs, dof_vm, dof_thm]
    dv = [dof_vs; dof_vm; dof_thm];
    cv = [1; -1; -rx];
    K(dv, dv) = K(dv, dv) + alpha * (cv * cv');
end
end
