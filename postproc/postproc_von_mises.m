function vm = postproc_von_mises(stress)
% Compute von Mises stress from continuum stress components.
%   stress - [nElems x 3]: [sxx, syy, txy]
%   vm     - [nElems x 1] von Mises stress
sxx = stress(:,1);
syy = stress(:,2);
txy = stress(:,3);
vm  = sqrt(sxx.^2 - sxx.*syy + syy.^2 + 3*txy.^2);
end
