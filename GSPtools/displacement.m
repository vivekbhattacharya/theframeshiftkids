% The function takes a sequence (without the 12-leader
% sequence), the number of codons, the differential vector
% from `diff_vector`, two lists of +1/-1 frameshifts against
% whcih to match.
% 
% This is the equivalent of an OOP private method. Its primary
% purpose is to aid in the development of unities.
function [x] = displacement(seq,Dvec,frontshifts,backshifts)
    % Refer to papers published by Dr. Bitzer, Dr. Ponalla, et al.
    % for meanings and derivations. C1 chosen specifically to
    % make prfB work, cf. Lalit et al.
    phi_sp = -30*(pi/180); initialx = 0.1;
    C1 = 0.005; C2 = initialx; spc = 1;
    
    global Travel;
    if isempty(Travel), load Travel.mat; end
    
    % ants: List of +1 frameshifts encountered.
    % termites: List of -1 frameshifts encountered.
    global ants termites;
    x = [0 C2]; ants = {}; termites = {};
    
    shift = 0; power = 10; wts = [];
    upper = length(Dvec) - 1; age_limit = 130;
    for k=2:upper
        % Choose appropriate codon, depending on the specified spacing, and
        % calculate nloop accordingly
        index = 3*(k-1) + 3*spc + shift;
        
        if(index + 4 > size(seq)), break; end;
        back_codon=seq(index:index+2);
        codon = seq(index+1:index+3);
        other_codon = seq(index+2:index+4);
    
        back_loops = real_loops(back_codon);
        here_loops = real_loops(codon);
        there_loops = real_loops(other_codon);
        
        % Again, refer to papers.
        phi_signal(1,k) = Dvec(k,2); x0 = x(1,k);    
        phi_dx = (pi/3)*x0 - phi_sp;
        
        % Refer to <http://code.google.com/p/theframeshiftkids/wiki/
        % MathBehindTheModel> for explanation.
        here_fail = 1; there_fail = 1; back_fail = 1; wt = 0;
        for wt=1:age_limit + 1
            a = (x0 - 2*shift); % Window function follows
            [back_fail, back] = probabilities(back_loops, exsin(a, 771)^power, back_fail);
            [here_fail, here] = probabilities(here_loops, excos(a)^10, here_fail);
            [there_fail, there] = probabilities(there_loops, exsin(a, 117)^power, there_fail);
            reloop = 1 * (1 - (here + there + back));
    
            r = rand; % Mersenne Twister
            if (reloop < here) || (reloop < there) || (reloop < back)
                if(r < here), break;
                elseif (r < here + there)
                    shift = shift + 1;
                    ants(end+1) = {[codon ',' num2str(k)]};
                    break;
                elseif (r < here + there + back)
                     shift = shift - 1;
                     termites(end+1) = {[codon ',' num2str(k)]};
                     break;
                end
            end
    
            phi_dx = ((pi/3)*x0) - phi_sp;
            dx = -C1*Dvec(k,1)*sin(phi_signal(1,k) + phi_dx);
            x0 = x0 + dx;
        end
        wts = [wts wt];
        if (wt > age_limit)
            fprintf('   %s at %g found Wichita\n', codon, k);
            ants = {'balls'}; termites = {'balls'};
            break;
        end;
        
        x(1,k+1) = x0;
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