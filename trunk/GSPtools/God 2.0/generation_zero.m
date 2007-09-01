function tavs = generation_zero(init, times)

% Takes a TAV vector and returns a matrix has slight variations of the initial.

L = length(init);
tavs = zeros(times, L);

% How many changes per mutation?
N = 6;

for row = 1:times
    changes = rand_int(1, L, N);
    tavs(row, :) = init;
    for j = 1:N
        r = 0.5 + 1.5 * rand;
        tavs(row, changes(j)) = tavs(row, changes(j)) * r;
    end
end

function [numbers] = rand_int(min, max, times)
delta = max - min
numbers = (round(delta .* rand(times, 1)) + min)';