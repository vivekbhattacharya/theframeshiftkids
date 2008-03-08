function [Dvec, Phase] = diff_vectors(Mag, Phase)

L = 3; upper = length(Mag);
x = min(max(1, (1:upper) - 1), upper - 1);

a = fakeslope([Mag(x); Mag(x+1)]);
b = fakeslope([Phase(x); Phase(x+1)]);
D = exp(j*Phase) .* (a + j*b.*Mag);
Dvec = [abs(D); angle(D)];

cond = find(Phase < 0);
Phase(cond) = Phase(cond) + 2*pi;
