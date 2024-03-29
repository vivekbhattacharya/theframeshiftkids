% Takes a TAV vector and returns a matrix has slight variations of the
% initial TAV.
function tavs = generation_zero(init, times)
    L = length(init);
    tavs = init(ones(1, times), :);

    % How many changes per mutation?
    N = 6;

    for row = 1:times
        changes = rand_int(L, N);
        % Mutate within interval [0.5, 2].
        r = -1 + 2 * rand(1, N);
        r = 2.^r;
        tavs(row, changes) = tavs(row, changes) .* r;
    end
end