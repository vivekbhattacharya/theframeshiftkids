% This function calculates the cummulative magnitude
% and phase of the input signal.
%
% USAGE:
%   [Mag, Phase, numcodons] = cumm_mag_phase(x)
%   x is a row vector that contains the free energy values
%   Phase is in radians

function [Mag, Phase] = cumm_mag_phase(x)
    % Signal rounded off to a codon multiple
    L = length(x);
    if rem(L,3) ~= 0, x = x(1:L-rem(L,3)); end;
    
    Mag = zeros(1,round(L/3)); Phase = zeros(1,round(L/3));
    for j=1:L/3
        [Mag(j), Phase(j)] = helper(x(1:j*3));
    end
end

function [A, theta] = helper(x)
    % This function calculates the magnitude and phase
    % of a given length of signal using sine-wave
    % interpolation.
    % 
    % INPUTS: 
    %   x = signal upto a certain number of codons
    %   The length of the signal must be a codon multiple
    % 
    % USAGE:
    %   [A,theta,Err] = helper(x, avg_choice);
    %   theta = phase angle in radians

    M = zeros(1,3); % Fill memory registers
    for j = 1:3:length(x)
        M(1:3) = M(1:3) + x(j:j+2);
    end    
    
    % The averaging calculation works only for the case when
    % length(x) is a multiple of 3. Let's assume that is true.
    % Then calculate magnitude and phase.
    M = M - mean(M); grass = M(3) + M(2);
    if M(1) ~= 0
        theta = atan2(M(1)*sqrt(3), M(1) + 2*M(2)); 
        A = M(1)/sin(theta);
    elseif grass ~= 0
        theta = atan2(grass*sqrt(3), M(3) - M(2));
        A = -grass/sin(theta);
    else A = 0; theta = 0;
    end
end