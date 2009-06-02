% A function that takes a codon number and returns the energy vector
% with magnitude and phase, along with a function to calculate the
% force from the codon number and position
function [energy, waveforce] = waveenergy(codon)
    global store Config;

    % Pad the signal with zeroes in accordance with polyenergy.
    new_signal = [zeros(1, 6) store.signal zeros(1, 6)];
    center = 3*codon + 5 + store.shift - 2;

    % Fits the energy through the function y = M cos(pi x/3 + theta).
    % Note that the positions for x are -2, 0, and 2.  Also, the
    % DC offset of the free energy is irrelevant, and so the values
    % for y are demeaned.
    x = -2:2:2;
    y = new_signal(center-1:center+1);
    y = y - mean(y);
    
    M = sqrt(4*(y(1)^2 + y(2)^2 + y(1)*y(2))/3);
    theta = acos(y(1)/M);
    energy = [M, theta];

    function [dx, force] = helper(codon, x)
        % F = -dU/dx; however, M and theta do not change, although
        % force is in the form y \propto M sin(pi x/3 + theta).
        % Also, force must have phi_sp subtracted off the argument of
        % the sin.
        force = energy;
        dx    = force(1) * sin((pi/3) * (x - 2 * store.shift) + force(2) - Config.phi_sp);
    end
    waveforce = @helper;
end

