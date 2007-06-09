function [Phase,x,diffx] = displacement(seq,Phase,numcodons,Dvec,frontshifts,backshifts)
% ------------------------------------------------------------
% CALCULATE DISPLACEMENT 
% ------------------------------------------------------------    

% Magic numbers
phi_sp=-30*(pi/180); initialx = 0.1;
C1 = 0.005; C2 = initialx; spc = 1;

% Initiate InstPhase array if so felt
x = [0 C2]; ants = {}; termites = {};

shift = 0;
for k=2:numcodons-1
    % Choose appropriate codon, depending on the specified spacing, and
    % calculate nloop accordingly
    initial = 3*(k-1) + 3*spc;
    
    if(initial+4+shift > size(seq)), break; end;
    back_codon=seq(initial+shift:initial+2+shift);
    codon = seq(initial+1+shift:initial+3+shift);
    other_codon = seq(initial+2+shift:initial+4+shift);

    back_loops = real_loops(back_codon);
    here_loops = real_loops(codon);
    there_loops = real_loops(other_codon);
    
    phi_signal(1,k) = Dvec(k,2); x0 = x(1,k);    
    phi_dx = (pi/3)*x0 - phi_sp;

    here_fail = 1; there_fail = 1; back_fail=1;
    for wt=1:1000
        a = (x0 - 2*shift)*pi/4; % Window function follows
        [back_fail, back] = probabilities(back_loops, exsin(a, '<')^10, back_fail);
        [here_fail, here] = probabilities(here_loops, excos(a)^10, here_fail);
        [there_fail, there] = probabilities(there_loops, exsin(a, '>')^10, there_fail);        
        reloop = 1 - here - there - back;

        r = rand; % Mersenne Twister
        if (reloop < here) || (reloop < there) || (reloop < back)
            if(r < here), break;
            elseif (r < here + there)
                shift = shift + 1;
                ants(length(ants)+1) = {[codon ',' num2str(k)]};
                break;
            elseif (r < here + there + back)
                shift = shift - 1;
                termites(length(termites) + 1) = {[codon ',' num2str(k)]};
                break;
            end;
        end;

        phi_dx = ((pi/3)*x0) - phi_sp;
        dx = -C1*Dvec(k,1)*sin(phi_signal(1,k) + phi_dx);
        x0 = x0 + dx;
    end
    
    %InstPhase(k) = (180/pi)*(phi_signal(1,k) + phi_dx);
    %if InstPhase(k)<0, InstPhase(k) = InstPhase(k) + 360; end
    x(1,k+1) = x0;
end

for k=1:length(Phase)
    if Phase(k)<0; Phase(k)=Phase(k)+(2*pi); end
end

diffx = zeros(length(x));
for k=1:length(x)-1; diffx(k) = x(k+1) - x(k); end;

%%% Counters for megaunity, a treatise into the
%%% failings of Matlab's strings
global shoals sands beached_whale;
sands = sands + 1;

% Handles edge case (which is quite often) where
% both `ants` and `frontshifts` and `backshifts` is [].
if strcmp(char(ants), char(frontshifts))
    if strcmp(char(termites), char(backshifts)), shoals = shoals + 1; end;
end

% Is verbosity disabled?
if (beached_whale), return; end;
fprintf('> '); disp_shifts(ants, '+1 frameshifts');
fprintf('< '); disp_shifts(termites, '-1 frameshifts');

% Iterate over the cell array and print the results
% sanely. If only Matlab had a join function like
% every other modern language in the world.
function disp_shifts(insects, id)
if size(insects) ~= size({})
    for i=1:length(insects), fprintf([insects{i} ';']); end;
end
fprintf('\n');

function [n] = real_loops(codon)
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