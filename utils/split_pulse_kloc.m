function varargout = split_pulse_kloc(x, nc, nt, outflag, ndk)
% convention: x is 1d flattened list of pulses (nt x nc x 2 (Re/Im)) and k-space
% locations (nt x 3)
% for p (and thus 2nd axis of Amatrix), time is fast and coil is slow index

% ndk: 2 or 3 for 2D or 3D k-locations

if nargin < 5; ndk=3; end
if nargin < 4; outflag='both'; end
    np = 2*nc*nt;

    p = x(1:np);
    k = x(np+1:end);

%     p = reshape(p, nt, nc);
    k = reshape(k, nt, ndk);

    if strcmp(outflag, 'both')
        varargout{1} = p;
        varargout{2} = k;
    elseif strcmp(outflag, 'k')
        varargout{1} = k;
    else
        varargout{1} = p;
    end
end