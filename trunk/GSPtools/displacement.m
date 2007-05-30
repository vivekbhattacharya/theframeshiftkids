function [Phase,x,diffx] = displacement(seq,Phase,numcodons,Dvec,love)
% ------------------------------------------------------------
% CALCULATE DISPLACEMENT 
% ------------------------------------------------------------    

% Magic numbers
phi_sp=-30*(pi/180); initialx = 0.1;
C1 = 0.005; C2 = initialx; spc = 1;

Nloop = []; x = [0 C2]; codon = []; InstPhase = [];
ants = [];

shift = 0;
for k=2:numcodons-1
    % Choose appropriate codon, depending on the specified spacing, and
    % calculate nloop accordingly
    initial = 3*(k-1) + 3*spc;
    try
        codon=seq(initial+1+shift:initial+3+shift);
        other_codon=seq(initial+2+shift:initial+4+shift);
    catch; break; end;

    Nloop(k) = nloopcalcify(codon);
    here_loops = real_loops(codon, k);
    there_loops = real_loops(other_codon, k);
    
    phi_signal(1,k) = Dvec(k,2); x0 = x(1,k);    
    phi_dx = ((pi/3)*x0)-phi_sp;

    here_fail = 1; there_fail = 1;
    for wt=1:1000
        a = (x0 - 2*shift)*pi/4; % Window function follows     
        [here_fail, here] = probabilities(here_loops, cos(a)^10, here_fail);
        [there_fail, there] = probabilities(there_loops, sin(a)^10, there_fail);        
        reloop = 1 - here - there;

        r = rand; % Mersenne Twister
        if (reloop < here) || (reloop < there)
            if(r < here), break;
            elseif (r < here + there)
                shift = shift + 1;
                ants =  strvcat(ants, [codon ',' num2str(k)]);
                break;
            end;
        end;

        phi_dx = ((pi/3)*x0) - phi_sp;
        dx = -C1*Dvec(k,1)*sin(phi_signal(1,k) + phi_dx);
        x0 = x0 + dx;
    end
    
    InstPhase(k) = (180/pi)*(phi_signal(1,k) + phi_dx);
    if InstPhase(k)<0, InstPhase(k) = InstPhase(k) + 360; end
    x(1,k+1) = x0;
end

for k=1:length(Phase)
    if Phase(k)<0; Phase(k)=Phase(k)+(2*pi); end
end
for k=1:length(x)-1; diffx(k)=x(k+1)-x(k); end;

%%% Counters for megaunity, a treatise into the
%%% failings of Matlab's strings
global shoals sands;
sands = sands + 1;

% Handles edge case (which is quite often) where
% both `ants` and `love` is [].
if strcmp(strvcat(ants, 'o'), strvcat(love, 'o'))
   shoals = shoals + 1;
end; 
   
if size(ants) ~= size([])
   pigs = cellstr(ants);
   for i=1:length(pigs), fprintf([pigs{i} ';']); end;
end

function [n] = real_loops(codon, index)
n = ceil(nloopcalcify(codon));
n = 2^(1/n);
n = n / (n - 1);

%% parameters bling
% loops: from real_loops
% weight: cos/sin factor
% sofar: acc fed back to us
function [acc, p] = probabilities(loops, weight, sofar)
acc = sofar * (1 - (weight/loops));
p = 1 - acc;