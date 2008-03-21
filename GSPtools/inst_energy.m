% Takes cumm_energy.m phasor and returns it in addition to a modified
% phase.
function [dvec] = inst_energy(mag, phase)
    L = 1; upper = length(mag);
    x = min(max(1, (1:upper) - L), upper - L);

    a = (mag(x + L) - mag(x))/L;
    b = (phase(x + L) - phase(x))/L;
    d = exp(j*phase) .* (a + j*b.*mag);
    dvec = [abs(d); angle(d)];
end