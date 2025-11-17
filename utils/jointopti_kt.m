function [p, k_opt, nrmse_simult, FAmap_simult, SAR_simult_VOP] = ...
    jointopti_kt(nkt, s, dist, VOPs, dB0, mask, target, k0, p0, SARLIM, switched, swmask, optonlyp, algo, stepdur, TR)

% joint optimization of pulse and kT points under VOP SAR constraints
% OUTPUTS:
%   p:              optimized pulse weights, [V], N_t x N_c
%   k_opt:          optimized kT locations [1/m], N_t x 3
%   nrmse_simult:   NRMSE of obtained FA w.r.t. target
%   FAmap_simult:   flip angle map for obtained pulse [deg], Nx x Ny x Nz

arguments
    nkt % number of kT points
    s % [T/V] voxel- and channel-wise B1+ maps (complex), N_r x N_c
    dist % [m], spatial coordinates of voxels withing mask, N_r x 3
    VOPs % VOP SAR matrices, Nc x Nc x NVOP
    dB0 % [Hz], off-resonance map, N_r x 1
    mask % spatial mask, Nx x Ny x Nz
    target % [rad], target FA vector, Nr x 1
    k0 % initial kT locations, nkt x 3
    p0 = [] % initial pulse weights, nkt x N_c
    SARLIM = 1e6 % [W/kg] local SAR constraint
    switched = 0 % flag, static or multiplexed?
    swmask = [] % switch state vector of length N_t, 1 for top row and 2 for bottom row
    optonlyp = 0 % flag: option to freeze k-space locations during optimization
    algo = 'interior-point' % algorithm for fmincon
    stepdur = 0.2e-3 % [s], sub-pulse duration
    TR = 10e-3 % [s], sequence repetition time, needed for duty cycle scaling of SAR
end


nc = size(VOPs,1);

if switched==1
    nt = 2*nkt; % switching -> two time steps per k-space location
else
    nt = nkt;
end


% initialize pulse
if isempty(p0)
    p0 = zeros(nt, nc);
end

% concatenated optimization vector: Re/Im of pulse + kloc
x0 = [complex2ReIm(vec(p0)); vec(k0)]; % 2*nc*nt + 3*nt


% boundaries for p and k
lb_p = -500.*ones(2*nc*nt,1);
ub_p = 500.*ones(2*nc*nt,1);

lb_k = -15.*ones(nt,3);
ub_k = +15.*ones(nt,3);

% fix last kT point to center
lb_k(end,:) = 0; 
ub_k(end,:) = 0;

if optonlyp % option to freeze k-space locations
    lb_k = k0;
    ub_k = k0;
end

lb = [lb_p; lb_k(:)];
ub = [ub_p; ub_k(:)];

% equality constraints for switching (need to tie two neighboring klocs together for multiplexing
if switched==1
    Aeq_k1 = kron(speye(nt/2), [1,-1]); % tie together single component of e.g. kx: neighboring have to be identical (k1-k2=0, k3-k4=0, ...)
    Aeq_k3 = blkdiag(Aeq_k1,Aeq_k1,Aeq_k1); % all 3 comps (x,y,z)
    Aeq_p = sparse(2*nc*nt,2*nc*nt); % pulse part: no constraint
    Aeq = blkdiag(Aeq_p, Aeq_k3); % should be 3*nt/2+2*length(vec(p0)) times length(x0) matrix (first is number of equality constraints: all 0's for p-part, and 3*nt/2 for k-loc component coupling part)
    beq = sparse(size(Aeq,1),1);
else
    Aeq = [];
    beq = [];
end

% PERFORM OPTIMIZATION
objfun = @(x) loss_simult_MLS(x, target, s, dist, stepdur, nc, nt, 0, dB0, swmask);
nonlinfcn = @(x) nonlinfcn_jointopti(x, nc, nt, VOPs, SARLIM, stepdur, TR);

options = optimoptions("fmincon",...
    'Algorithm',algo,...
    "Display","iter-detailed", ...
    "MaxIterations",200, ...
    "MaxFunctionEvaluations", 5000,...
    "SpecifyObjectiveGradient",true,...
    "SpecifyConstraintGradient",true,...
    "PlotFcn", {'optimplotfval'}...
    );

[xopt,fval,exitflag,output,lambda,grad,hessian] = fmincon(objfun,x0,[],[],Aeq,beq,lb,ub,nonlinfcn,options);

% evaluate obtained pulse
[p_opt, k_opt] = split_pulse_kloc(xopt,nc,nt,'both');
A = createAmat(s,dist,k_opt,stepdur,swmask,dB0);
p_simult = reIm2Complex(p_opt);
FAmap_simult = toSpatial(calcFA(p_simult,A), mask);
nrmse_simult = nrmse2(abs(rad2deg(target)),abs(FAmap_simult(~isnan(FAmap_simult))));
[~,SAR_simult_VOP] = calcSarForPulse(p_simult, VOPs, stepdur, TR, [], 1);
p = reshape(p_simult,nt,[]);