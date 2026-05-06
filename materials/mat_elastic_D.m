function D = mat_elastic_D(E, nu, formulation)
% Plane-stress or plane-strain constitutive matrix D [3x3].
%   E           - Young's modulus
%   nu          - Poisson's ratio
%   formulation - 'plane_stress' or 'plane_strain'
switch lower(formulation)
    case 'plane_stress'
        c = E / (1 - nu^2);
        D = c * [1   nu  0;
                 nu  1   0;
                 0   0   (1-nu)/2];
    case 'plane_strain'
        c = E / ((1 + nu) * (1 - 2*nu));
        D = c * [1-nu  nu    0;
                 nu    1-nu  0;
                 0     0     (1-2*nu)/2];
    otherwise
        error('mat_elastic_D: formulation must be ''plane_stress'' or ''plane_strain''.');
end
end
