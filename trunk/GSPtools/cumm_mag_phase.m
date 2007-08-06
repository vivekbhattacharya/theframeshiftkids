% This function calculates the cummulative magnitude
% and phase of the input signal.
%
% USAGE:
%   [Mag, Phase, numcodons] = cumm_mag_phase(x)
%   x is a row vector that contains the free energy values
%   Phase is in radians

function [mag, phase] = cumm_mag_phase(signal)
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
        M = registers - mean(registers);
        grass = M(3) + M(2);
        if M(1) ~= 0
            phase(i) = atan2(M(1)*sqrt(3), M(1) + 2*M(2)); 
            mag(i) = M(1)/sin(phase(i));
        elseif grass ~= 0
            phase(i) = atan2(grass*sqrt(3), M(3) - M(2));
            mag(i) = -grass/sin(phase(i));
        else mag(i) = 0; phase(i) = 0;
        end
    end
end