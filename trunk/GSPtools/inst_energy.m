function [dvec] = inst_energy(mag, phase)
    % Takes cumm_energy.m phasor and returns it in addition to a modified
    % phase.

    % mag(1) = []; phase(1) = [];
    L = 3; upper = length(mag);
    x = min(max(1, (1:upper) - 1), upper - 1);

    a = fakeslope([mag(x); mag(x+1)]);
    b = fakeslope([phase(x); phase(x+1)]);
    D = exp(j*phase) .* (a + j*b.*mag);
    dvec = [abs(D); angle(D)];
end

function [m] = fakeslope(y)
    m = (y(2, :) - y(1, :))/2; % 3 - 1 = 2
end