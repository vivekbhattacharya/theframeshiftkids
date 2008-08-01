% Calculates and smooths the inst. energy vector and produces a matrix
% of slopes and intercepts so you can determine force just by flooring
% and looking it up.
function [force] = force(signal)
    dsignal = dumbdiff(signal);

    % Test length checks
    % dsignal = [1 1 1];
    % Test correctness
    % dsignal = [1 1 1 5 5 5 5 5 5 5 9 9 9 11];
    % Test correctness
    % dsignal = [1 2 3 4];
    % Test evenness
    % dsignal = [1 1 5 5 3 3 6 6 6 6 6 6 6];

    force = smooth(dsignal);
end

function [force] = smooth(raw_force)
    if length(raw_force) <= 2
        error('Unable to smooth inst. energy vector: too short.');
    end
    g = group(raw_force);
    midpoints = cumsum(g(2, :)) - (g(2, :)-1)/2;
    slopes = dumbdiff(g(1, :)) ./ dumbdiff(midpoints);

    % Have a zero slope for both the initial midpoint, which I added
    % here, and the last midpoint.
    midpoints = [0 midpoints];
    slopes = [0 slopes 0];
    % intercept = y - x * slope
    intercepts = [g(1, 1) g(1, :)] - slopes .* midpoints;

    % force = [expanded slopes; expanded intercepts].
    force = zeros(2, length(raw_force));
    % Append length(raw_force) to take care of base pairs in the interval
    % [midpoints(end), index of last base].
    midpoints = [midpoints length(raw_force)];
    for i = 1:length(midpoints)-1
        lower = ceil(midpoints(i)) + 1;
        upper = floor(midpoints(i+1)) + 1;
        force(1, lower:upper) = slopes(i);
        force(2, lower:upper) = intercepts(i);
    end
end
