% Calculates weights (we also call it exposures) of a given "window".
% x = back. y = current. z = next. The period is four, but we've
% optimized that away.
function [w] = weights(x)
    a = 0; b = 0; c = 0;
    if (x > -4 && x < 0), a = sin(x * pi/4 + pi); end;
    if (x > -2) && (x < 2), b = cos(x * pi/4); end;
    if (x > 0) && (x < 4), c = sin(x * pi/4); end;

    w = [a b c];
end
