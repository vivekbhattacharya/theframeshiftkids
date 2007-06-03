function [Mag, Phase, numcodons] = calc_cumm_mag_phase(signal)

% ------------------------------------------------------------
% CALCULATE CUMMULATIVE MAGNITUDE AND PHASE
% ------------------------------------------------------------    
[Mag, Phase, Err] = cumm_mag_phase(signal);
I1 = find(Err(1,:));
if ~isempty(I1)
    fprintf(1,'Magnitude negative at %d indices:\n',length(I1));
end

I2 = find(Err(2,:)); 
if ~isempty(I2)
    fprintf(1,'Equations not satisfied at %d indices:\n',length(I2));
end
%fprintf(1,'\nNumber of codons: %d',length(Mag));
numcodons = length(Mag);