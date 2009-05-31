% This function calculates the magnitude and phase of the cummulative
% energy that arises from the input signal (cf. Kidnap).
%
% USAGE:
%   [mag, phase] = cumm_energy(signal)
%   signal is an array of free energy values from Kidnap
%   phase is in radians, adjusted to be greater than zero.

function [mag, phase] = cumm_energy(signal)
    % Round signal off to a codon multiple
    limit = floor(length(signal)/3);
    mag   = zeros(3, limit);
    phase = zeros(3, limit);
    for i = 1:3
        % Pad the first element of the signal with zero if cumm_energy_frame
        % rounds off a codon.
        s = [signal(i:end) zeros(1, i-1)];
        [mag(i, :) phase(i, :)] = cumm_energy_frame(s);
    end
end

function [mag, phase] = cumm_energy_frame(signal)
    % Round signal off to a codon multiple
    limit = floor(length(signal)/3);
    % Get a matrix of [1 2 3; 4 5 6; 7 8 9; et cetera].
    indices = reshape(1:limit*3, 3, limit);
    % mem stands for "memory" as represented by three registers.
    mem = signal(indices);

    for i = 2:limit
        % Accumulate and center the free energy wave at zero.
        x = mem(:, i) + mem(:, i - 1);
        mem(:, i) = x - mean(x);
    end

    % The atan2 expression simplifies into atan(tan(theta)) = theta, given
    % r.1 = M*sin(theta), r.2 = M*sin(theta + 2*pi/3), r.3 =
    % M*sin(theta + 4*pi/3). We worked it out.
    phase = atan2(mem(1, :)*sqrt(3), mem(1, :) + 2*mem(2, :));

    warning off MATLAB:divideByZero;
    mag = mem(1, :) ./ sin(phase);
    mag(isnan(mag)) = 0;
    warning on MATLAB:divideByZero;
end
