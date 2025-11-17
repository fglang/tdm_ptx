function Afull = createAmat(s,dist,k,stepdur,swmask,b0,blipdur)
% system matrix for the spatial domain method
arguments
    s % [T/V] voxel- and channel-wise B1+ maps (complex), N_r x N_c
    dist % [m], spatial coordinates of voxels withing mask, N_r x 3
    k % [1/m], k-space locations, nkt x 3 (kx, ky, kz)
    stepdur % [s], sub-pulse duration
    swmask = [] % switch state vector of length N_t, 1 for top row and 2 for bottom row
    b0 = [] % [Hz], off-resonance map, N_r x 1
    blipdur = 0 % [s], gradient blip duration
end

nTx = size(s,2);
ntimes = size(k,1);
gamma = 2 * pi * 42576384.74; % Hz/T

% off-resonance effect
if isempty(b0)
    b0eff = 0;
else
    T = ntimes * (stepdur + blipdur);
    remainingTime = [stepdur/2];
    for kk=1:ntimes-1
        remainingTime(end+1) = remainingTime(end) + blipdur + stepdur;
    end
    remainingTime = remainingTime - T;
    b0eff = 2 * pi * b0 * remainingTime; % phase evolution due to dB0
end

% flatten dist if necessary
if ndims(dist) > 2
    dist = reshape(dist,[],3); % nvox x 3
end

% phase evolution due to k-space trajectory
temp = 2.*pi.* dist * k.';

% Fourier part of system matrix
A = 1j * gamma * stepdur * exp(1j * b0eff + 1j * temp);

% Tx sensitivity part of system matrix
if isempty(swmask) % switched case
    A1 = repmat(A,1,nTx);
    A2 = repelem(s,1,size(A,2));
    Afull = A1 .* A2;
else % static case
    s_sw = s(:,:,swmask); % vox, cha, time
    s_sw = permute(s_sw,[1,3,2]); % vox, time, cha

    Afull = A .* s_sw;
    Afull = reshape(Afull,[],nTx*ntimes);
end

end