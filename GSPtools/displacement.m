% The function takes a sequence (without the 12-leader
% sequence), the number of codons, the differential vector
% from `diff_vector`, two lists of +1/-1 frameshifts against
% which to match.
function [d, waits] = displacement(seq, dvec_, fs)
    global Config dvec store;
    dvec = dvec_;

    % ants: List of +1 frameshifts encountered.
    % termites: List of -1 frameshifts encountered.
    upper = length(dvec) - 1;

    % Fudging factor: initial displacement.
    store = struct('x', [Config.init_disp], ...
                   'shift', 0, 'wts', zeros(1, upper-2), ...
                   'ants', [], 'termites', [], 'anthill', []);

    % For the common case of no actual frameshifts,
    % avoid computing displacement deviation.
    for k = 1:upper
        index = 3*k + store.shift;
        if index + 4 > length(seq), break; end;
        overaged = loop(seq(index:index+4), k);

        % Check if frameshift occurred too late.
        if Config.dire
            if handle_aging(overaged, store.anthill, fs), break; end;
            % I expected a frameshift. None occurred. Abort.
            if (length(fs) > 0) && (k > fs(1)) && (length(store.anthill) == 0), break; end;
        end
    end
    d = store.x;
    waits = store.wts;

    % Tally times displacement called. shoals can be the total deviation
    % or times antill == fs.
    global sands shoals;
    if isempty(sands), sands = 0; end;
    if isempty(shoals), shoals = 0; end;

    sands = sands + 1;
    if Config.yield == 0
        criticals = zeros(1, length(d));
        criticals(fs+1) = 2;
        criticals = cumsum(criticals);
        deviation = mean((d - criticals) .^ 2);
        shoals = shoals + sqrt(deviation);
    % Barring no backframeshifts, check for antill-fs equality.
    elseif Config.yield == 1
        if length(store.termites) > 0, return;
        elseif size(fs) == size(store.anthill)
            if isempty(fs) || all(fs == store.anthill), shoals = shoals + 1; end;
        end
    end
end

function [break_p] = handle_aging(code, anthill, fs)
    break_p = 0;
    % These imply Config.dire.
    switch code
      case 3, break_p = 1; % Backframeshift => automatic error
      case 4 % Frameshift => check if valid
        if length(fs) == 0, break_p = 1;
        elseif anthill(end) ~= fs(1), break_p = 1;
        end
    end
end

% Refer to papers published by Dr. Bitzer, Dr. Ponalla, et al.
% Config.phi_sp chosen specifically to make prfB work, cf. Lalit et
% al. This is heavily optimized. Do not refer to it for the math.
function [overaged] = loop(piece, k)
    global store Config dvec;
    overaged = 0;

    % [back_fail, here_fail, there_fail]
    fails = [1 1 1];
    back_codon = piece(1:3); codon = piece(2:4); there_codon = piece(3:5);
    loops = real_loops(back_codon, codon, there_codon);

    x0 = store.x(k);
    wt = 0;
    % This is where Nloop used to be.
    for wt = 1:1000 + 1
        % Window function
        w = realpow(weights(x0 - 2*store.shift), Config.power);
        fails = fails .* (1 - w ./ loops);
        probs = 1 - fails;

        % [here, here + there, here + there + back]
        % 1 - cumprobs(3) = reloop.
        cumprobs = cumsum(probs); r = rand;
        if any(cumprobs(3) > 1 - probs)
            if r < cumprobs(1), break;
            elseif r < cumprobs(2)
                store.shift = store.shift + 1;
                store.anthill(end+1) = k;
                store.ants{end+1} = codon;

                % If there's a frameshift, tell handle_aging to check if it's the
                % wrong one under dire.
                overaged = 4; break;
            elseif r < cumprobs(3)
                store.shift = store.shift - 1;
                store.termites{end+1} = [codon ' ' num2str(k)];

                % There's a backframeshift, tell handle_aging to stop automatically.
                overaged = 3; break;
            end
        end

        % This follows from phi_signal(1,k) = dvec(2, k); see
        % "A model for +1 frameshifts in eubacteria" by Ponnala, et al.
        phi_dx = pi*x0/3 - Config.phi_sp;
        dx = dvec(1, k) * sin(dvec(2, k) + phi_dx);
        x0 = x0 + -Config.c1 * dx;
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
