
function myTravel = get_travel(names, tav)

limit = length(names);
str = 'myTravel = struct(';

for i = 1:limit
    k = nloopcalc(names{i},0,1, names, tav, 1000);
    str = [str '''' names{i} ''',' num2str(k) ','];
end

str = str(1:end-1);
str = [str ');'];
eval(str);

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

function N = nloopcalc(x,Case,Type,Names,TAV,Nstop)

if length(x)~=3
    error('Codon is not of length 3 !!');
end

TAV2 = TAV(find(TAV));

% Stop codons should have high wait times, Manually set values for E.coli,
% since the TAV will be zero for these codons
if strcmp(x,'uaa')==1
    N = Nstop;        
elseif strcmp(x,'uag')==1
    N = Nstop;
elseif strcmp(x,'uga')==1
    N = Nstop;
else
    % Compare x with Names and extract abundance ratio
    I=strmatch(x,Names,'exact');
    if length(I)~=1
        error('No unique match for the codon');
    else
        N = (max(TAV2)/min(TAV2))-floor(TAV(I)/min(TAV2));
    end
end

