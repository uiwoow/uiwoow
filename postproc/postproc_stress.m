function [stress, strain] = postproc_stress(mesh, mat, u)
% Recover element stresses and strains from displacement vector.
%   mesh   - mesh struct
%   mat    - material struct (with mat.D for continuum, or mat.E/A/Iz for frame)
%   u      - [nDof x 1] global displacement vector
%
%   stress - [nElems x nComp]  stress components per element
%   strain - [nElems x nComp]  strain components per element
%
%   For T3/Q4:    columns = [sxx, syy, txy]
%   For truss2d:  columns = [axial_force]
%   For beam2d:   columns = [axial_force, shear_force, moment_node1, moment_node2]

switch mesh.type
    case {'T3', 'Q4'}
        if isfield(mat, 'D')
            D = mat.D;
        else
            D = mat_elastic_D(mat.E, mat.nu, mat.formulation);
        end
        stress = zeros(mesh.nElems, 3);
        strain = zeros(mesh.nElems, 3);
        for e = 1:mesh.nElems
            nodes_e = mesh.conn(e, :);
            xy_e    = mesh.nodes(nodes_e, :);
            dofs_e  = elem_dofs(mesh, e);
            u_e     = u(dofs_e);
            switch mesh.type
                case 'T3'
                    [~, ~, B] = elem_T3(xy_e, D, mesh.t);
                case 'Q4'
                    % Evaluate B at centroid (xi=0, eta=0)
                    B = q4_B_centroid(xy_e);
            end
            eps_e       = B * u_e;
            sig_e       = D * eps_e;
            strain(e,:) = eps_e';
            stress(e,:) = sig_e';
        end

    case 'truss2d'
        stress = zeros(mesh.nElems, 1);
        strain = zeros(mesh.nElems, 1);
        for e = 1:mesh.nElems
            nodes_e = mesh.conn(e, :);
            xy_e    = mesh.nodes(nodes_e, :);
            dofs_e  = elem_dofs(mesh, e);
            u_e     = u(dofs_e);
            dx = xy_e(2,1) - xy_e(1,1);
            dy = xy_e(2,2) - xy_e(1,2);
            L  = sqrt(dx^2 + dy^2);
            c  = dx/L; s = dy/L;
            % Axial elongation: delta = (u2-u1)*c + (v2-v1)*s
            delta      = (u_e(3)-u_e(1))*c + (u_e(4)-u_e(2))*s;
            strain(e)  = delta / L;
            stress(e)  = mat.E * strain(e) * mat.A;  % axial force N
        end

    case 'beam2d'
        % Returns [N, V, M1, M2] per element (internal forces at nodes)
        stress = zeros(mesh.nElems, 4);
        strain = zeros(mesh.nElems, 4);
        for e = 1:mesh.nElems
            nodes_e = mesh.conn(e, :);
            xy_e    = mesh.nodes(nodes_e, :);
            dofs_e  = elem_dofs(mesh, e);
            u_e     = u(dofs_e);
            dx = xy_e(2,1) - xy_e(1,1);
            dy = xy_e(2,2) - xy_e(1,2);
            L  = sqrt(dx^2 + dy^2);
            c  = dx/L; s = dy/L;
            T6 = [c  s  0  0  0  0;
                 -s  c  0  0  0  0;
                  0  0  1  0  0  0;
                  0  0  0  c  s  0;
                  0  0  0 -s  c  0;
                  0  0  0  0  0  1];
            u_local = T6 * u_e;
            EA = mat.E * mat.A;
            EI = mat.E * mat.Iz;
            N  =  EA/L  * (u_local(4) - u_local(1));
            V  =  12*EI/L^3 * (u_local(2) - u_local(5)) + ...
                   6*EI/L^2 * (u_local(3) + u_local(6));
            M1 =  6*EI/L^2 * (u_local(2) - u_local(5)) + ...
                   EI*(4/L*u_local(3) + 2/L*u_local(6));
            M2 = -6*EI/L^2 * (u_local(2) - u_local(5)) + ...
                   EI*(2/L*u_local(3) + 4/L*u_local(6));
            stress(e,:) = [N, V, M1, M2];
        end
        strain = stress;  % not meaningful separately; use stress
end
end

function B = q4_B_centroid(xy)
% B matrix for Q4 at centroid (xi=0, eta=0).
xi = 0; eta = 0;
dN_dxi  = [-(1-eta)  (1-eta)  (1+eta) -(1+eta)] / 4;
dN_deta = [-(1-xi)  -(1+xi)   (1+xi)   (1-xi) ] / 4;
J       = [dN_dxi; dN_deta] * xy;
invJ    = inv(J);
dN_dx   = invJ(1,1)*dN_dxi + invJ(1,2)*dN_deta;
dN_dy   = invJ(2,1)*dN_dxi + invJ(2,2)*dN_deta;
B = zeros(3, 8);
B(1, 1:2:end) = dN_dx;
B(2, 2:2:end) = dN_dy;
B(3, 1:2:end) = dN_dy;
B(3, 2:2:end) = dN_dx;
end
