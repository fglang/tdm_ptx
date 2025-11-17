function out = constraint_SAR(xri, Qavg, SARLIM, stepdur, TR, return_grad)
% constraint function for local SAR

arguments
    xri % real flattened pulse vector (Re/Im)
    Qavg % SAR matrices, Nc x Nc x NVOP
    SARLIM % [W/kg], local SAR limit
    stepdur = 0.2e-3; % [s], sub-pulse duration
    TR = 10e-3; % [s], sequence repetition time
    return_grad = false; % flag: false: return function value, true: return Jacobian
end

x = reIm2Complex(xri);
if ~return_grad % function
    out = calcSarForPulse(x, Qavg, stepdur, TR, [], true) - SARLIM; % negative means not violated
else % Jacobian, see doi:10.1109/TMI.2013.2295465
    % basically this is the analytical derivative of a quadratic form
    [nc,~,nQ] = size(Qavg);
    x = reshape(x, [], nc);
    Qt = permute(Qavg,[3,1,2]); % vox,c,c
    Qt = reshape(Qt,[],nc); % vox*c, c
    T = 2 .* Qt * x.';
    T = permute(reshape(T, nQ,nc,[]),[3,2,1]);
    T = reshape(T,[],nQ);
    out = cat(1, real(T), imag(T)) .* stepdur ./ TR;
end
end