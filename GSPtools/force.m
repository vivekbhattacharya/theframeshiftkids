% Calculates the inst. energy vector.
function [force] = force(signal)
    force = -dumbdiff(signal);
end
