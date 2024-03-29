% A function that takes a codon number and returns the energy polynomial
% and function to calculate the force from the codon number and the
% position
function [energy, polyforce] = polyenergy(codon)
    global store;

    % Pad the signal with zeroes in case we use larger powers.
    new_signal = [zeros(1, 6) store.signal zeros(1, 6)];
    center = 3*codon + 5 + store.shift - 2;

    % Defines power, bounds for the codons, and the limits of
    % regression.  Takes into account that a distance of one
    % nucleotide is actually 2 in terms of displacement.
    power = 4;
    lower = floor(-power/2);
    upper = floor(power/2);
    x = 2*(lower:upper);
    y = new_signal(center+lower:center+upper);
    energy = polyfit(x, y, power);

    function [dx, force] = helper(codon, x)
        % F = -dU/dx; need to drop the last element of the force
        % vector after differentiation because Matlab starts the
        % vector with the coefficient of the highest power of
        % the polynomial.
        force = -energy .* (power:-1:0);
        force = force(1:power);
        dx    = polyval(force, x - 2*store.shift);
    end
    polyforce = @helper;
end

