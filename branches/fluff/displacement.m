function [Phase,x,diffx] = displacement(seq,Phase,numcodons,Dvec,love)
% ------------------------------------------------------------
% CALCULATE DISPLACEMENT 
% ------------------------------------------------------------    

% Magic numbers
phi_sp=-30*(pi/180); initialx = 0.1;
C1 = 0.005; C2 = initialx; spc = 1;
ftws = []; lols = []; wtfs = [];

% Initiate InstPhase array if so felt
x = [0 C2]; ants = {};

shift = 0;
for k=2:numcodons-1
    % Choose appropriate codon, depending on the specified spacing, and
    % calculate nloop accordingly
    initial = 3*(k-1) + 3*spc;
    
    if(initial+4+shift > size(seq)), break; end;
    codon = seq(initial+1+shift:initial+3+shift);
    other_codon = seq(initial+2+shift:initial+4+shift);

    here_loops = real_loops(codon);
    there_loops = real_loops(other_codon);
    
    phi_signal(1,k) = Dvec(k,2); x0 = x(1,k);    
    phi_dx = (pi/3)*x0 - phi_sp;

    here_fail = 1; there_fail = 1;
    for wt=1:1000
        a = (x0 - 2*shift)*pi/4; % Window function follows     
        [here_fail, here] = probabilities(here_loops, cos(a)^10, here_fail);
        [there_fail, there] = probabilities(there_loops, sin(a)^10, there_fail);        
        reloop = 1 - here - there;
		
%          ftws(wt) = here;
%          lols(wt) = there;
%          wtfs(wt) = reloop;

        r = rand; % Mersenne Twister
        if (reloop < here) || (reloop < there)
            if(r < here), break;
            elseif (r < here + there)
                shift = shift + 1;
                ants(length(ants)+1) = {[codon ',' num2str(k)]};
                break;
            end;
        end;

        phi_dx = ((pi/3)*x0) - phi_sp;
        dx = -C1*Dvec(k,1)*sin(phi_signal(1,k) + phi_dx);
        x0 = x0 + dx;
    end
    
%      figure(2);
%          subplot(211); plot(0,0);plot(1:length(ftws), ftws, 1:length(lols), lols, 1:length(wtfs), wtfs);
%      pause(0.1);
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
% both `ants` and `love` is [].
if strcmp(char(ants), char(love)), shoals = shoals + 1; end;

% Is verbosity disabled?
if (beached_whale), return; end;
if size(ants) ~= size([])
   for i=1:length(ants), fprintf([ants{i} ';']); end;
end

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