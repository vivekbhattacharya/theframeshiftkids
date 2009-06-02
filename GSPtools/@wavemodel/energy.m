% An input into the loop() protocol: returns a vector of energies at
% this given codon. Will later be passed back into force(). Energies
% are computed according to a non-register sine fit model of a moving
% window on the free-energy signal.
function [energy] = energy(self, codon)
    shift   = get(self.model, 'shift');
    psignal = get(self.model, 'psignal');

    center = 3*codon + 5 + shift - 2;
    lim    = 1;

    % Our window.
    x = -2:2:2;
    y = psignal(center-lim:center+lim);
    y = y - mean(y);

    M = sqrt(4*(y(2)^2 + y(3)^2 + y(2)*y(3))/3);
    if M == 0
        % We never use theta in any significant way if M is zero.
        theta = 1281;
    else
        % Thanks to floating-point numbers, y(2)/M can be *slightly*
        % less than -1 or *slightly* greater than +1, leading Matlab
        % to return pi +/- 0.000[snip]01i for acos. Kludge:
        theta = real(acos(y(2)/M));
        if abs(y(1) - M * cos(-2*pi/3 + theta)) > 0.01
            theta = 2*pi - theta;
        end
    end
    energy = [M theta];
end
