function tavs = generation_zero(init, times)

% Takes a TAV vector and returns a matrix has slight variations of the initial.

L = length(init);
tavs = zeros(times, L);

% How many changes per mutation?
N = 6;

for row = 1:times
    changes = unique(rand_int(1, L, N));
    tavs(row, :) = init;
    for j = 1:length(changes)
        r = 0.5 + 1.5 * rand;
        tavs(row, changes(j)) = tavs(row, changes(j)) * r;
    end
end

function [numbers] = rand_int(min, max, times)
delta = max - min;

% The augend creates numbers in the set [0, delta].
% The addend morphs it into [min, delta + min] = [min, max].
% Do not use ceiling!
numbers = (round(delta .* rand(times, 1)) + min)';