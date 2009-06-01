% A method for the class model, part of the loop() protocol. Returns
% the force at the given position, given the codon number and energy
% vector.
function dx = force(model, energy, codon, x)
    % F = -dU/dx. We need to drop the last element of the force vector
    % after differentiation because Matlab starts the vector with the
    % coefficient of the highest power of the polynomial.
    power = model.power;
    force = -energy .* (power:-1:0);
    force = force(1:power);
    dx    = polyval(force, x - 2*model.shift);
end
