% Calculates weights (we also call it exposures) of a given "window".
% x = back. y = current. z = next. The period is four, but we've
% optimized that away.
function [w] = weights(x)
    w = [0 0 0];
    if (x > -4 && x < 0), w(1) = x * pi/4 + pi; end;
    % This is cosine.
    if (x > -2 && x < 2), w(2) = x * pi/4 + pi/2; end;
    if (x > 0 && x < 4), w(3) = x * pi/4; end;
    % Conveniently, sin(0) is 0, which avoids a find() for us.
    w = sin(w);
end
