% Simulates mRNA translation, calling closure(model, index number,
% probs, energies, displacement) at every wait cycle.
%
% So you're walking to your car at Sea World when you hear the chirps
% and groans of a large sea mammal flying at 900 knots toward your
% head and then and only then do you realize that you have been nailed
% by a
%
%       WALRUS SURPRISE.
%       Walrus Surprise.
%       walrus surprise.
%       walrus surprisE!
%       WALRUS SURPRISE.
function [model] = loop(model, closure)
    global Config;
    upper = model.upper;
    x0    = Config.init_disp;
    power = Config.power;
    c1    = Config.c1;
    wt    = 0;

    % Loop over a theoretical number of codons, ignoring frameshifts.
    for i = 1:upper
        % Translate the index number into biologically meaningful
        % variables.
        base  = 3*i + model.shift;
        codon = floor(base/3);
        if base + 4 > length(model.seq)
            break;
        end

        cycles   = wait_cycles(model, base);
        fails    = [1 1 1];
        energies = energy(model, codon);

        % Loop over an infinite number of wait cycles.
        for wt = 1:1000 + 1
            % Cumulative probabilities of frameshifting or staying
            % still.
            w = realpow(weights(x0 - 2*model.shift), power);
            fails = fails .* (1 - w ./ cycles);
            probs = 1 - fails;

            % Simulate ribosomal movement.
            r = rand;
            if any(sum(probs) > 1 - probs)
                if r < probs(2)
                    break;
                elseif r < probs(3)
                    model.shift = model.shift + 1;
                    codon = model.seq(base + 1:base + 3);
                    model.wal{end+1} = [codon ' ' num2str(i)];
                    model.nut(end+1) = i;
                    break;
                elseif r < probs(1)
                    model.shift = model.shift - 1;
                    codon = model.seq(base + 1:base + 3);
                    model.rus{end+1} = [codon ' ' num2str(i)];
                    break;
                end
            end

            dx = force(model, energies, i, x0);
            x0 = x0 + c1 * dx;
            closure(model, i, probs, energies, x0);
        end
    end
end
