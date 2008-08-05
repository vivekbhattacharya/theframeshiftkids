% The function takes a sequence (without the 12-leader sequence), the
% number of codons, the force vector from `force`, and a vector of
% codons indices that should frameshift. (For prfB, this is [25].)
function [disp, waits] = displacement(seq, force, fs)
    global Travel Config store;

    upper = length(force) - 1;

    % ants: List of +1 frameshifts encountered.
    % termites: List of -1 frameshifts encountered.
    % Initial displacement is a fudging factor.
    store = struct('x', [Config.init_disp], 'shift', 0, 'wts', ...
                   zeros(1, upper-2), 'ants', [], 'termites', [], ...
                   'anthill', [], 'force', force);

    % For the common case of no actual frameshifts,
    % avoid computing displacement deviation.
    for k = 1:upper
        index = 3*k + store.shift;
        if index + 4 > length(seq), break; end;
        should_break = loop(seq(index:index+4), k);

        % Check if frameshift occurred too late.
        if Config.dire
            switch should_break
              case 1, break; % Backframeshift => automatic error
              case 2 % Frameshift => check if valid
                if anthill(end) ~= fs(1), break; end;
            end
            % I expected a frameshift. None occurred. Abort.
            if (length(fs) > 0) && (k > fs(1)) && (length(store.anthill) == 0), break; end;
        end
    end
    disp = store.x;
    waits = store.wts;
    update_globals(fs, disp);
end

function update_globals(fs, disp)
    % Tally times displacement called. shoals can be the total deviation
    % or times antill == fs.
    global sands shoals store Config;
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
        if length(store.termites) > 0, return;
        elseif size(fs) == size(store.anthill)
            if isempty(fs) || all(fs == store.anthill), shoals = shoals + 1; end;
        end
    end
end

% Refer to papers published by Dr. Bitzer, Dr. Ponalla, et al. This is
% heavily optimized. Do not refer to it for the math.
function [should_break] = loop(piece, k)
    config; global Config store;
    should_break = 0;

    % [back_fail, here_fail, there_fail]
    fails = [1 1 1];
    back_codon = piece(1:3); codon = piece(2:4); there_codon = piece(3:5);
    loops = real_loops(back_codon, codon, there_codon);

    displace = store.x(k);
    wt = 0;
    % This is where Nloop used to be.
    for wt = 1:1000 + 1
        % Window function
        w = realpow(weights(displace - 2*store.shift), Config.power);
        fails = fails .* (1 - w ./ loops);
        probs = 1 - fails;

        r = rand;
            if r < probs(2), break;
            elseif r < probs(2) + probs(3)
                store.shift = store.shift + 1;
                store.anthill(end+1) = k;
                store.ants{end+1} = codon;

                % Dire Mode: If there's a frameshift, check if it's the wrong one.
                should_break = 2;
                break;
            elseif r < sum(probs)
                store.shift = store.shift - 1;
                store.termites{end+1} = [codon ' ' num2str(k)];

                % Dire Mode: There's a backframeshift, so stop automatically.
                should_break = 1;
                break;
            end

        n_base = 3*k - 0.5 + displace/2;
        n_base_index = round(n_base);

        previous = store.force(n_base_index - 1);
        next = store.force(n_base_index);

        % Slope is rise over run, and run is one here. The code condenses the
        % following three lines into one line.

        % slope = next - previous;
        % intercept = next - slope * (n_base_index - 0.5);
        % dx = slope * n_base + intercept;

        dx = (next - previous) * (n_base - n_base_index + 0.5) + next;
        displace = displace + Config.c1 * dx;
    end
    store.x(k+1) = displace;
    store.wts(k) = wt;
end

% Calculate normalized wait cycles.
function [loops] = real_loops(a, b, c)
    global Travel;
    loops = [Travel.(a) Travel.(b) Travel.(c)];
    loops = 2 .^ (1 ./ ceil(loops));
    loops = loops ./ (loops - 1);
end
