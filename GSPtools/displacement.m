% The function takes a sequence (without the 12-leader
% sequence), the number of codons, the differential vector
% from `diff_vector`, two lists of +1/-1 frameshifts against
% whcih to match.
% 
% This is the equivalent of an OOP private method. Its primary
% purpose is to aid in the development of unities.
function [x] = displacement(seq,Dvec,frontshifts,backshifts)
    global Travel;
    if isempty(Travel), load Travel.mat; end
    
    % ants: List of +1 frameshifts encountered.
    % termites: List of -1 frameshifts encountered.
    global ants termites;
    x = [0 0.1]; ants = {}; termites = {};

    % Length of codon sequence
    upper = length(Dvec) - 1;
    shift = 0; wts = [];
    for k=2:upper
        % Choose appropriate codon, depending on the specified spacing, and
        % calculate nloop accordingly
        index = 3*k + shift;
        
        if(index + 4 > size(seq)), break; end;
        [wt, x0, shift, overaged] = loop(shift, x(k), seq(index:index+4), k, Dvec(k, :));
        wts = [wts wt];
        if (overaged == 1)
            fprintf('   %s at %g found Wichita\n', seq(index+1:index+3), k);
            ants = {'witch'}; termites = {'coven'};
            break;
        elseif (overaged == -1)
            fprintf('   %s at %g found a stop codon\n', seq(index+1:index+3), k);
            break;
        end
        
        x(k+1) = x0;
    end
    % wts
     
    % Tally the number of times the gene sequence
    % correctly frameshifted (shoals) in addition to the
    % total number of displacement.m calls (sands).
    global shoals sands;
    if isempty(sands), shoals = 0; sands = 0; end
    sands = sands + 1;
    if strcmp(char(ants), char(frontshifts))
        if strcmp(char(termites), char(backshifts)), shoals = shoals + 1; end;
    end
end

function [dead] = is_stopper(codon)
    global Travel ants termites;
    if (Travel.(codon) == 1000)
        msg = ['died at ' codon];
        ants(end+1) = {msg};
        termites(end+1) = {msg};
        dead = true;
    else dead = false; end;
end

function [wt, x0, shift, overaged] = loop(shift, x0, fragment, k, diff)
    % Refer to papers published by Dr. Bitzer, Dr. Ponalla, et al.
    % for meanings and derivations. C1 chosen specifically to
    % make prfB work, cf. Lalit et al.
    C1 = 0.005; overaged = 0;
    age_limit = 200; power = 10; phi_sp = -30*(pi/180);
    
    wt = 0; here_fail = 1; back_fail = 1; there_fail = 1;
    
    back_codon = fragment(1:3);
    codon = fragment(2:4);
    there_codon = fragment(3:5);
    
    back_loops = real_loops(back_codon);
    here_loops = real_loops(codon);
    there_loops = real_loops(there_codon);
    
    global ants termites;
    for wt=1:age_limit + 1
        a = (x0 - 2*shift); % Window function follows
        [back_fail, back] = probabilities(back_loops, exsin(a, 771)^power, back_fail);
        [here_fail, here] = probabilities(here_loops, excos(a)^10, here_fail);
        [there_fail, there] = probabilities(there_loops, exsin(a, 117)^power, there_fail);
        reloop = 1 * (1 - (here + there + back));

        r = rand; % Mersenne Twister
        if (reloop < here) || (reloop < there) || (reloop < back)
            if(r < here)
                if (is_stopper(codon)), overaged = -1; return; end;
                break;
            elseif (r < here + there)
                shift = shift + 1;
                ants(end+1) = {[codon ',' num2str(k)]};
                if (is_stopper(there_codon)), overaged = -1; return; end;
                break;
            elseif (r < here + there + back)
                shift = shift - 1;
                termites(end+1) = {[codon ',' num2str(k)]};
                if (is_stopper(back_codon)), overaged = -1; return; end;
                break;
            end
        end
        phi_dx = ((pi/3)*x0) - phi_sp;
        
        % This follows from phi_signal(1,k) = Dvec(k,2)
        dx = -C1*diff(1)*sin(diff(2) + phi_dx);
        x0 = x0 + dx;
    end
    if (wt > age_limit), overaged = true; end;
end

% Calculates Nloops per <http://code.google.com/p/
% theframeshiftkids/wiki/MathBehindTheModel> for
% the given codon.
function [n] = real_loops(codon)
    global Travel;
    n = ceil(Travel.(codon));
    n = 2^(1/n);
    n = n / (n - 1);
end

% Calculates probability per <http://code.google.com/p/
% theframeshiftkids/wiki/MathBehindTheModel> for
% any given "thinking" time slice.
% Parameters:
%   loops: from real_loops
%   weight: cos/sin factor
function [acc, p] = probabilities(loops, weight, acc)
    acc = acc*(1 - weight/loops);
    p = 1 - acc;
end