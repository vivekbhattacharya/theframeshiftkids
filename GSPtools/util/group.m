% Groups vectors together.
function [grouped] = group(vec)
    if length(vec) == 0
        grouped = []; return;
    end

    if length(vec) == 1
        grouped = [vec(1); 1]; return;
    end

    % The starts of new groups will have non-zero values.
    raw_starts = [vec(1) (vec(2:end) - vec(1:end - 1))];
    starts = vec(find(raw_starts));
    % Where are these starts?
    start_indices = [find(raw_starts) length(vec) + 1];
    % How long is each group?
    lengths = start_indices(2:end) - start_indices(1:end - 1);
    grouped = [starts; lengths];
end
