function [c, ceq, gc, gceq] = nonlinfcn_jointopti(x, nc, nt, VOPs, SARLIM, stepdur, TR)
% non-linear constraint function for SAR constraints
arguments
    x % real flattened vector of optimizable variables: pulse (Re/Im) and kloc
    nc % number of channels
    nt % number of time steps
    VOPs % VOP SAR matrices, Nc x Nc x NVOP
    SARLIM % [W/kg], local SAR limit
    stepdur = 0.2e-3; % [s], sub-pulse duration
    TR = 0.2e-3; % [s], sequence repetition time
end

% inequality constraint function
c = constraint_SAR(split_pulse_kloc(x,nc,nt,'p'), VOPs, SARLIM, stepdur, TR, 0);
if nargout > 1
    % equality constraint function
    ceq = [];
end
if nargout > 2
    % gradient of inequality constraint function
    gc = cat(1, constraint_SAR(split_pulse_kloc(x,nc,nt,'p'),VOPs, SARLIM, stepdur, TR, 1), zeros(nt*3, size(VOPs,3)));
end
if nargout > 3
    % gradient of equality constraint function
    gceq = [];
end