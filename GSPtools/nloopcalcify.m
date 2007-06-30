function N = nloopcalcify(x)
% NLOOPCALC : "Calculate number of loops"
% This function calculates the looping number for a specified input codon. 
% Contains an additional input parameter (Nstop) compared to nloopcalc
% 
% USAGE: 
% y = nloopcalc(x,Case,Type,Names,TAV,Nstop)
% TAV = codon abundance ratio table
% Names = cell array or character array of the codons in lower-case RNA format
% Nstop = number of loops for the stop codon (set using prfB)
% Note: Elements of TAV should have a one-to-one match with those of Names
global TAV Codon2Index; Nstop = 1000;
if length(x) ~= 3
    error('Codon is not of length 3.');
end

global maxTAV minTAV;
if isempty(maxTAV), maxTAV = max(TAV); end;
if isempty(minTAV), minTAV = min(TAV(find(TAV))); end;

% Stop codons should have high wait times.
% I manually set values for E.coli since the TAV will
% be zero for these codons.
I = getfield(Codon2Index, x);
if I >= 37 && I <= 39
    N = Nstop;
else
    % Extract abundance ratio, but first find out
    % where codon is in TAV. We don't have a hashmap
    % in Matlab, so here's how we fudge it.
    if length(I) ~= 1
        error('No unique match for the codon');
    end
    N = (maxTAV/minTAV) - floor(TAV(I)/minTAV);
end