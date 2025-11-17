function [f, g] = loss_simult_MLS(x, b, s, dist, stepdur, nc, nt, lbd_power, dB0, swmask)
%%%
% magnitude least squares loss function
% f = || |Ax| - |b| ||_2^2 
% optional: lbd_power regularizes L2 norm of x
%
%   xri:    pulse vector in [Re, Im] format
%   A:      system matrix
%   b:      RHS (target flip angle)
%%%

arguments
    x % real flattened vector of optimizable variables: pulse (Re/Im) and kloc
    b % [rad], target flip angle vector
    s % [T/V] voxel- and channel-wise B1+ maps (complex), N_r x N_c
    dist % [m], spatial coordinates of voxels withing mask, N_r x 3
    stepdur % [s], sub-pulse duration
    nc % number of channels
    nt % number of time steps
    lbd_power = 0 % optional power regularization (L2)
    dB0 = [] % [Hz], off-resonance map, N_r x 1
    swmask = [] % switch state vector of length N_t, 1 for top row and 2 for bottom row
end

[pri,k] = split_pulse_kloc(x, nc, nt, 'both');

% create system matrix
A = createAmat(s,dist,k,stepdur,swmask,dB0);

% compute MLS loss
p = reIm2Complex(pri);
Ap = A*p;
y = abs(Ap) - abs(b);

% optional L2 regularization
if lbd_power > 0
    pulse_power = sum(abs(p).^2./50);
    f = sum(y.^2, 'all') + lbd_power .* pulse_power;
else
    f = sum(y.^2, 'all');
end

% analytical Jacobian
if nargout > 1 
    %%% gradient wrt pulse weights, see doi:10.1002/mrm.25902
    temp1 = y.*exp(1j.*angle(Ap)); % residual * phase term (from magnitude operation)
    temp = A'*temp1;
    g_p_Re = 2*real(temp);
    g_p_Im = 2*imag(temp);
    if lbd_power > 0
        g_p = [g_p_Re; g_p_Im] + 2.*lbd_power.*pri ./ 50;
    else
        g_p = [g_p_Re; g_p_Im];
    end

    %%% gradient wrt k-space locations, 
    % chain rule, combined known analytical derivatives for MLS objective (doi:10.1002/mrm.25902) and Fourier matrix (doi:10.1002/mrm.21262)
    g_k = 4.*pi.*conj(p).*A' * (dist .* temp1); 
    g_k = squeeze(imag(sum(reshape(g_k, nt, nc, []), 2))); % sum over coils, nt x 3 remains, [kx,ky,kz] on last axis
    
    %%% total Jacobian for pulse samples and k-space locations
    g = [g_p; g_k(:)];
    
end
end