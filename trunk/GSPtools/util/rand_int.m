% `count` number of unique random integers in the interval 1:max.
function [numbers] = rand_int(max, count)
    [ignore, a] = sort(rand(1, max)); % randperm
    numbers = a(1:count);
end