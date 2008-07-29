% The function takes a sequence (without the 12-leader
% sequence), the number of codons, the differential vector
% from `diff_vector`, two lists of +1/-1 frameshifts against
% which to match.
function [disp, waits] = displacement(seq, dvec, fs)
    global Travel Config;

    % ants: List of +1 frameshifts encountered.
    % termites: List of -1 frameshifts encountered.
    global ants anthill termites store;
    ants = {}; termites = {}; anthill = [];

    upper = length(dvec) - 1;
    % Fudging factor: initial displacement
    store = struct('x', [Config.init_disp], 'shift', 0, 'wts', zeros(1, upper-2));

    % For the common case of no actual frameshifts,
    % avoid computing displacement deviation.
    for k = 1:upper
        index = 3*k + store.shift;
        if(index + 4 > length(seq)), break; end;
        overaged = loop(seq(index:index+4), k, dvec(:, k));

        if overaged ~= 0,
            if handle_aging(overaged, anthill, fs) == 1, break; end;
        end;

        % Check if frameshift occurred too late.
        if Config.dire
            if (length(fs) > 0) && (k > fs(1)) && (length(anthill) == 0), break; end;
        end
    end
    disp = store.x;
    waits = store.wts;

    % Tally times displacement called. shoals can be the total deviation
    % or times antill == fs.
    global sands shoals;
    if isempty(sands), sands = 0; end;
    if isempty(shoals), shoals = 0; end;

    sands = sands + 1;
    if Config.yield == 0
        criticals = zeros(1, length(disp));
        criticals(fs+1) = 2;
        criticals = cumsum(criticals);
        deviation = mean((disp - criticals) .^ 2);
        shoals = shoals + sqrt(deviation);
    % Barring no backframeshifts, check for antill-fs equality.
    elseif Config.yield == 1
        if length(termites) > 0, return;
        elseif size(fs) == size(anthill)
            if isempty(fs) || all(fs == anthill), shoals = shoals + 1; end;
        end
    end
end

function [break_p] = handle_aging(code, anthill, fs)
    break_p = 0;
    switch code
      case 1
        fprintf(': %s at %g found ribosomal hyperpause\n', seq(index+1:index+3), k);
        break_p = 1;
      case -1
        fprintf(': %s at %g found a stop codon\n', seq(index+1:index+3), k);
        break_p = 1;

      % These imply Config.dire.
      case 3, break_p = 1; % Backframeshift => automatic error
      case 4 % Frameshift => check if valid
        if anthill(end) ~= fs(1), break_p = 1; end;
    end
end

% Refer to papers published by Dr. Bitzer, Dr. Ponalla, et al.
% Config.phi_sp chosen specifically to make prfB work, cf. Lalit et
% al. This is heavily optimized. Do not refer to it for the math.
function [overaged] = loop(piece, k, diff)
    global Config; config;
    overaged = 0; age_limit = 1000; power = Config.power; c1 = Config.c1;

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
                break;
            elseif r < here + there
                store.shift = store.shift + 1;
                anthill(end+1) = k;
                ants{end+1} = codon;

                % If there's a frameshift, check if it's the wrong one under dire.
                if Config.dire, overaged = 4; end;
                break;
            elseif r < here + there + back
                store.shift = store.shift - 1;
                termites{end+1} = [codon ' ' num2str(k)];

                % There's a backframeshift, so stop automatically.
                if Config.dire, overaged = 3; end;
                break;
            end
        end

        % This follows from phi_signal(1,k) = dvec(2, k); see
        % "A model for +1 frameshifts in eubacteria" by Ponnala, et al.
        phi_dx = ((pi/3)*x0) - Config.phi_sp;
        dx = -c1 * diff(1) * sin(diff(2) + phi_dx);
        x0 = x0 + dx;
    end
    if Config.detect_pauses
        if (wt > age_limit), overaged = 1; end;
    end
    store.x(k+1) = x0;
    store.wts(k) = wt;
end

% Calculate normalized wait cycles.
function [loops] = real_loops(a, b, c)
    global Travel;
    loops = [Travel.(a) Travel.(b) Travel.(c)];
    loops = 2 .^ (1 ./ ceil(loops));
    loops = loops ./ (loops - 1);
end
