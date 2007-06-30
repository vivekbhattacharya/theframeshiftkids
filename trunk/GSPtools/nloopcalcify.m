function N = nloopcalcify(x, TAV)
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
Nstop = 1000;
if length(x) ~= 3
    error('Codon is not of length 3.');
end

% Stop codons should have high wait times.
% I manually set values for E.coli since the TAV will
% be zero for these codons.
I = getfield(TAV, x);
if I == 0
    N = Nstop;
else
    % Extract abundance ratio, but first find out
    % where codon is in TAV. We don't have a hashmap
    % in Matlab, so here's how we fudge it.
    N = (TAV.max/TAV.min) - floor(I/TAV.min);
end