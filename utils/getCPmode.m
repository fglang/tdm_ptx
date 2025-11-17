function p_CP = getCPmode(nc, F, sgn)
% get pulse vector of the circularly-polarized (CP) mode
% e.g. for an circular 8ch coil array, this means phases of
% 0°, -45°, -90°,... and equal magnitude across channels.
% INPUT:
%   nc:     number of channels
%   F:      order of the CP mode (1 for standard CP, 2 for CP2+ etc.
%   sgn:    phase sign (+-45° for 2nd element?)
% OUTPUT:
%   p_CP:   complex pulse vector
if nargin < 3; sgn=-1; end
if nargin < 2; F=1; end
p_CP = exp(sgn.*1j.*(0:F.*2.*pi/nc:F*(2*pi-2.*pi/nc))).';