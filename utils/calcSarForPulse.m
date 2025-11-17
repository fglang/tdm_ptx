function [SAR,maxLocalSAR] = calcSarForPulse(pulse,Qmat,stepdur,TR,inputPower,quiet, doAtt)
% calculate SAR for given pulse and Q matrices
% OUTPUTS:
%   - SAR: [W/kg], SAR values from all provided matrices (3rd axis of Qmat)
%   - maxLocalSAR: [W/kg], maximum local SAR across provided matrices
arguments
    pulse % single time step: [Nc x 1] vector of pulse weights, multiple time steps: [Nt x Nc] pulse matrix
    Qmat % SAR matrices, Nc x Nc x NQ
    stepdur = 0.2e-3; % [s] duration of pulse time step
    TR = 10e-3; % [s] sequence repetition time
    inputPower = [] % optional: scale for given input power
    quiet = 1 % verbosity flag
    doAtt = 0 % doAtt: if 1 or true, consider pulse attenuation of -2 dB
end


R = 50; % [Ohm]
[nc,~,nm] = size(Qmat);

if doAtt % consider attenuation
    Att = -2.0; % dB
    pulse = pulse .* 10.^(Att/20);
end

nt = length(pulse(:)) ./ nc;

if nt == 1
    pulse = vec(pulse);
end

power_applied = sum(abs(pulse).^2/R) ./ nt;

if nt == 1 % single time point
    SAR = pagemtimes(pagemtimes(pulse',Qmat),pulse);
else % multiple time points
    pulse_r = reshape(pulse,[],nc);
    SAR = sum(Qmat .* repmat(pulse_r'*pulse_r, 1,1,nm), [1,2]); % see doi:10.1002/mrm.24138, Eq.5
end

SAR = abs(squeeze(SAR)); % small imag will be there due to numerical error

SAR = SAR .* stepdur./TR; % duty cylce rescaling

% optional: scale SAR to desired input power
if exist('inputPower','var') && ~isempty(inputPower)
    SAR = SAR./power_applied .* inputPower;
    power_applied = inputPower;
end

maxLocalSAR = max(abs(SAR));

if ~exist('quiet','var') || quiet ~= true
    fprintf('----------------------------------\n')
    fprintf('Total applied power: % 2.3f W\n',power_applied)
    fprintf('Max local SAR (10g): % 2.3f W/kg\n',maxLocalSAR);
end