function N = nloopcalc(x)
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

TAV2 = TAV(find(TAV));
% Stop codons should have high wait times, Manually set values for E.coli,
% since the TAV will be zero for these codons
if strmatch(x, ['uaa'; 'uga'; 'uga']) >= 1
    N = Nstop;        
else
    % Compare x with Names and extract abundance ratio
    I = getfield(Codon2Index, x);
    if length(I) ~= 1
        error('No unique match for the codon');
    end
    N = (max(TAV2)/min(TAV2))-floor(TAV(I)/min(TAV2));
end