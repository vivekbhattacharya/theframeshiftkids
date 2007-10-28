% The function takes a sequence (without the 12-leader
% sequence), the number of codons, the differential vector
% from `diff_vector`, two lists of +1/-1 frameshifts against
% which to match.
function [x, waits] = displacement(seq, Dvec, fs)
    global Travel Config;
    
    % ants: List of +1 frameshifts encountered.
    % termites: List of -1 frameshifts encountered.
    global ants anthill termites store;
    ants = {}; termites = {}; anthill = [];

    upper = length(Dvec) - 1;
    % Fudging factor: initial displacement
    store = struct('x', [0 Config.init_disp], 'shift', 0, 'wts', zeros(1, upper-2));
    
    % For the common case of no actual frameshifts,
    % avoid computing the actual shift ever.
    away = 0; no_frameshifting = isempty(fs);
    if ~no_frameshifting, criticals = find_criticals(fs); end
    for k = 2:upper
        % Take the codon and calculate nloop
        index = 3*k + store.shift;
        
        if(index + 4 > length(seq)), break; end;
        overaged = loop(seq(index:index+4), k, Dvec(k, :));
        switch overaged
          case 1
            fprintf(': %s at %g found ribosomal hyperpause\n', seq(index+1:index+3), k);
            break;
          case -1
            fprintf('%s at %g found a stop codon\n', seq(index+1:index+3), k);
            break;
          case 3
            if Config.dire, break; end;
          case 4
            if Config.dire
                if anthill(end) ~= fs(1), break; end;
            end
        end
        if Config.dire
            if k > fs(1)
                if length(anthill) == 0, break; end;
            end
        end
        
        if ~no_frameshifting
            away = away + (store.x(k+1) - criticals(k))^2;
        else, away = away + store.x(k+1)^2; end;
    end
    x = store.x;
    waits = store.wts;
     
    % Tally the number of times the gene sequence
    % correctly frameshifted (foals) in addition to the
    % total number of displacement.m calls (sands) in
    % addition to the total deviation (coals).
    global sands shoals;
    if isempty(sands), sands = 0; end;
    if isempty(shoals), shoals = 0; end;
    
    sands = sands + 1;
    if Config.yield == 0
        shoals = shoals + sqrt(away/upper);
    elseif Config.yield == 1
        if size(fs) == size(anthill)
            if isempty(fs), shoals = shoals + 1; end;
            if fs == anthill, shoals = shoals + 1; end;
        end
    end
end

% Returns a function that retrieves the actual
% shift given a codon number.
function [get] = find_criticals(f)
    % Use vectorized assignment. Cascade the shifts.
    c(f) = 2;
    c = cumsum(c);

    % Embodies the simple logic needed to read a value
    % from the array `find_criticals` returns. Take the
    % end because shifts cascade.
    function [delta] = helper(k)
        if k > length(c), delta = c(end);
        else delta = c(k);
        end
    end
    
    get = @helper;
end

function [dead] = is_stopper(codon)
    global Travel ants termites;
    dead = false;
    
    % termites doubles as a disp_shifts flag for now.
    if (Travel.(codon) == 1000)
        termites{end+1} = 0;
        dead = true;
    end
end

% Refer to papers published by Dr. Bitzer, Dr. Ponalla, et al.
% for meanings and derivations. C1 chosen specifically to
% make prfB work, cf. Lalit et al. This is a heavily optimized
% function. Do not refer to it for the math.
function [overaged] = loop(piece, k, diff)
    overaged = 0; age_limit = 1000; power = 10;
    
    % [back_fail, here_fail, there_fail]
    fails = [1 1 1]; wt = 0;
    back_codon = piece(1:3); codon = piece(2:4); there_codon = piece(3:5);
    loops = real_loops(back_codon, codon, there_codon);
    
    global store ants anthill termites Config;
    x0 = store.x(k);
    for wt=1:age_limit + 1
        % Window function
        w = realpow(weights(x0 - 2*store.shift), power);
        fails = fails .* (1 - w ./ loops);
        probs = 1 - fails;
        
        reloop = 1 - sum(probs);
        back = probs(1); here = probs(2); there = probs(3);
        r = rand;
        % Config checks are here to minimize function call overhead.
        if (reloop < here) || (reloop < there) || (reloop < back)
            if r < here
                if Config.detect_stops
                    if is_stopper(codon), overaged = -1; end;
                end
                break;
            elseif r < here + there
                store.shift = store.shift + 1;
                anthill(end+1) = k;
                ants{end+1} = codon;
                
                if Config.detect_stops
                    if is_stopper(there_codon), overaged = -1; end;
                end
                if Config.dire, overaged = 4; end;
                break;
            elseif r < here + there + back
                store.shift = store.shift - 1;
                termites{end+1} = [codon ' ' num2str(k)];

                if Config.detect_stops
                    if is_stopper(back_codon), overaged = -1; end;
                end
                if Config.dire, overaged = 3; end;
                break;
            end
        end
        
        % This follows from phi_signal(1,k) = Dvec(k,2)
        % "A model for +1 frameshifts in eubacteria" by Ponnala, et al.
        phi_dx = ((pi/3)*x0) - Config.phi_sp;
        dx = -0.005 * diff(1) * sin(diff(2) + phi_dx);
        x0 = x0 + dx;
    end
    if Config.detect_pauses
        if (wt > age_limit), overaged = 1; end;
    end
    store.x(k+1) = x0;
    store.wts(k-1) = wt;
end

% Calculates Nloops per <http://code.google.com/p/
% theframeshiftkids/wiki/MathBehindTheModel> for
% the given codon.
function [loops] = real_loops(a, b, c)
    global Travel;
    loops = [Travel.(a) Travel.(b) Travel.(c)];
    loops = 2 .^ (1 ./ ceil(loops));
    loops = loops ./ (loops - 1);
end