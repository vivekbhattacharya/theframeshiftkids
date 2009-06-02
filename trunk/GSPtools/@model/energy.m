% An input into the loop() protocol: returns a vector of energies at
% this given codon. Will later be passed back into force(). Energies
% are computed according to a fourth-degree polynomial fit model of a
% moving window on the free-energy signal.
function [energy] = energy(model, codon)
    power  = model.power;
    lower  = floor(-power/2);
    upper  = floor(power/2);
    center = 3*codon + 5 + model.shift - 2;

    % Our window.
    x      = 2*(lower:upper);
    y      = model.psignal(center+lower:center+upper);
    energy = polyfit(x, y, power);
end
