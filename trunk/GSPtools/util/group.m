% Groups vectors together. Completely vectorized.
%
% USAGE:
%  >>> group([1 2 2 2 3 4 4 5 6 6 1])
%  [1 2 3 4 5 6 1; 1 3 1 2 1 2 1]
function [grouped] = group(vec)
    if length(vec) == 0
        grouped = []; return;
    end

    if length(vec) == 1
        grouped = [vec(1); 1]; return;
    end

    % Obfuscate zeroes.
    vec(find(vec == 0)) = -42;
    % The starts of new groups will have non-zero values.
    raw_starts = [vec(1) (vec(2:end) - vec(1:end - 1))];
    starts = vec(find(raw_starts));
    % Where are these starts?
    start_indices = [find(raw_starts) length(vec) + 1];
    % How long is each group?
    lengths = start_indices(2:end) - start_indices(1:end - 1);
    % Unobfuscate zeroes.
    starts(find(starts == -42)) = 0;
    grouped = [starts; lengths];
end
