% Takes cumm_energy.m phasor and returns it in addition to a modified
% phase.
function [dvec] = inst_energy(mag, phase)
    L = 1; upper = length(mag);
    x = min(max(1, (1:upper) - L), upper - L);
    y = [x(2:end) x(end)];

    dvec = diffvec(mag(x), mag(y), phase(x), phase(y));
end

function [dvec] = diffvec(m1,m2,p1,p2)
    % m1, p1 is the first vector in polar; m2, p2 the second
    dx = m2 .* cos(p2) - m1 .* cos(p1);
    dy = m2 .* sin(p2) - m1 .* sin(p1);
    dvec = [sqrt(dx .^ 2 + dy .^ 2); atan2(dy, dx)];
end