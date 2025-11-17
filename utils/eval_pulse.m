function [FAmap, pnrmse, ppower, SAR_VOP] = eval_pulse(p, A, mask, targFA, VOPs, swmask, stepdur, TR)
% evaluate pulse
arguments
    p % complex pulse vector, N_t x N_c
    A % spatial domain system matrix
    mask % spatial mask, Nx x Ny x Nz
    targFA % [deg], target flip angle
    VOPs % VOP SAR matrices, Nc x Nc x NVOP
    swmask = [1] % switch state vector of length N_t, 1 for top row and 2 for bottom row
    stepdur = 0.2e-3 % [s], time step duration
    TR = 10e-3; % [s], sequence repetition time
end

nc = size(VOPs,1);

FAmap = toSpatial(calcFA(p,A), mask);
pnrmse = nrmse2(targFA, abs(FAmap));
ppower = calcPowerForPulse(p, nc, stepdur, TR);
[~, SAR_VOP] = calcSarForPulse(p, VOPs, stepdur, TR);

end