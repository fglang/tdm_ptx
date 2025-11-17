function power_applied = calcPowerForPulse(pulse, nc, stepdur, TR)
arguments
    pulse
    nc = 16;
    stepdur = 0.2e-3;
    TR = 10e-3;
end

R = 50; % [Ohm]
power_applied = sum(abs(pulse(:)).^2/R) .* stepdur ./ TR;