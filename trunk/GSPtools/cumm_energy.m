% This function calculates the cummulative magnitude
% and phase of the input signal.
%
% USAGE:
%   [mag, phase] = cumm_energy(signal)
%   signal is an array of free energy values from Kidnap
%   phase is in radians, adjusted to be greater than zero.

function [mag, phase] = cumm_energy(signal)
    % Round signal off to a codon multiple
    L = length(signal);

    limit = floor(L/3);
    mag = zeros(1, limit); phase = zeros(1, limit);
    registers = zeros(1,3);

    % Calculate the magnitude and phase  of the signal
    % using sine-wave interpolation.
    for i=1:limit
        index = (i-1)*3 + 1;
        registers = registers + signal(index:index+2);

        % The averaging calculation works only for the case when
        % length is a multiple of 3. Let's assume that is true.
        % Then calculate magnitude and phase.
        %
        % atan2(...) translates into atan(tan(theta)). Trust us. We
        % worked it out. (Given r.1 = M * sin(theta), r.2 = M *
        % sin(theta + 2*pi/3), r.3 = M * sin(theta + 4*pi/3)).
        %
        % We center the free energy wave at zero first.
        M = registers - mean(registers);

        if M(1) ~= 0
            phase(i) = atan2(M(1)*sqrt(3), M(1) + 2*M(2));
            mag(i) = M(1)/sin(phase(i));
        else mag(i) = 0; phase(i) = 0;
        end
    end
end
