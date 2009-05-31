% Takes cumm_energy.m's three phasors and returns the instantaneous
% energy at each codon, suitable for calculating displacement.
function [dvec] = inst_energy(mag, phase)
    upper = length(mag(1, :));
    x = max(1, (1:upper) - 1);
    y = x([2:end end]);

    [m, p] = diffvec(mag(1:3, x), mag(1:3, y), ...
                     phase(1:3, x), phase(1:3, y));
    dvec(:, :, 1) = m;
    dvec(:, :, 2) = p;
end

function [mag, phase] = diffvec(m1,m2,p1,p2)
    % m1, p1 is the first vector in polar; m2, p2 the second
    dx    = m2 .* cos(p2) - m1 .* cos(p1);
    dy    = m2 .* sin(p2) - m1 .* sin(p1);
    mag   = sqrt(dx .^ 2 + dy .^ 2);
    phase = atan2(dy, dx);
end
