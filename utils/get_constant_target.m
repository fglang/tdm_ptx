function target = get_constant_target(targFA, mask, refphase)
% generate target vector for pulse design with constant target flip angle
% INPUTS:
%   targFA:     target flip angle [deg]
%   mask:       spatial mask, Nx x Ny x Nz
%   refphase:   target phase, optional
% OUTPUT:
%   target:     vector of targets N_r x 1, complex
if nargin < 3; refphase = []; end
target = ones(sum(mask(:)>0),1) .* deg2rad(targFA);
target = toSpatial(target,mask);
target = target .*mask;
if ~isempty(refphase)
    target = target .* refphase;
end
target = target(~isnan(target));