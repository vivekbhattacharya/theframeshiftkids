% The function takes a sequence (without the 12-leader
% sequence), the number of codons, the differential vector
% from `diff_vector`, two lists of +1/-1 frameshifts against
% which to match.
function [disp, waits] = displacement(seq, dvec, fs, varargin)
    global Travel Config store;

    upper = length(dvec) - 1;

    function do_nothing(x0, probs, k, diff), end
    if nargin > 3
        chunky_closure = varargin{1};
    else
        chunky_closure = @do_nothing;
    end

    % ants: List of +1 frameshifts encountered.
    % termites: List of -1 frameshifts encountered.
    % Initial displacement is a fudging factor.
    store = struct('x', [Config.init_disp], ...
                   'shift', 0, ...
                   'wts', zeros(1, upper-2), ...
                   'ants', [], ...
                   'termites', [], ...
                   'anthill', [], ...
                   'chunky_closure', chunky_closure);

    % For the common case of no actual frameshifts,
    % avoid computing displacement deviation.
    for k = 1:upper
        frame      = mod(store.shift, 3) + 1;
        %frame      = 1;
        x(1, :)    = dvec(frame, :, 1);
        x(2, :)    = dvec(frame, :, 2);
        dvec_frame = x;

        index = 3*k + store.shift;
        codon = k + floor(store.shift/3);
        %codon = k;

        if index + 4 > length(seq)
            break;
        end

        overaged = loop(seq(index:index+4), k, dvec_frame(:, codon));

        % Check if frameshift occurred too late.
        if Config.dire
            if handle_aging(overaged, store.anthill, fs), break; end;
            % I expected a frameshift. None occurred. Abort.
            if (length(fs) > 0) && (k > fs(1)) && (length(store.anthill) == 0)
                break;
            end
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

% Refer to papers published by Dr. Bitzer, Dr. Ponalla, et al.
% Config.phi_sp chosen specifically to make prfB work, cf. Lalit et
% al. This is heavily optimized. Do not refer to it for the math.
function [overaged] = loop(piece, k, diff)
    config; global Config;
    overaged = 0;

    % [back_fail, here_fail, there_fail]
    fails = [1 1 1];
    back_codon = piece(1:3); codon = piece(2:4); there_codon = piece(3:5);
    loops = real_loops(back_codon, codon, there_codon);

    global store Config;
    x0 = store.x(k);

    % This is where Nloop used to be.
    wt = 0;
    for wt = 1:1000 + 1
        % Window function
        w = realpow(weights(x0 - 2*store.shift), Config.power);
        fails = fails .* (1 - w ./ loops);
        probs = 1 - fails;

        r = rand;
        % Config checks are here to minimize function call overhead.
        if any(sum(probs) > 1 - probs)
            if r < probs(2)
                % Stay where we are.
                break;
            elseif r < probs(3)
                % Forward frameshift.
                store.shift = store.shift + 1;
                store.anthill(end+1) = k;
                store.ants{end+1} = codon;

                % If there's a frameshift, check if it's the wrong one under dire.
                overaged = 4;
                break;
            elseif r < probs(1)
                % Backward frameshift.
                store.shift = store.shift - 1;
                store.termites{end+1} = [codon ' ' num2str(k)];

                % There's a backframeshift, so stop automatically.
                overaged = 3;
                break;
            end
        end

        % This follows from phi_signal(1,k) = dvec(2, k); see
        % "A model for +1 frameshifts in eubacteria" by Ponnala, et al.

        % We have to subtract 2*store.shift from the position because
        % inst_energy now accounts for frameshifts, not us. (It's a
        % shift in responsibility, you see.)
        phi_dx = (pi/3)*(x0 - 2*store.shift) - Config.phi_sp;
        dx = diff(1) * sin(diff(2) + phi_dx);
        x0 = x0 + -Config.c1 * dx;

        store.chunky_closure(x0 - 2*store.shift, probs, k, diff);
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
