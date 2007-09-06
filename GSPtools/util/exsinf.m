function [n] = exsinf(x)

% period = 4;
% phase = pi/2 - 2*pi/period; => 0

% Draw a graph of sine at x = 0 to see why we're doing this.
% if (x > (2-period/2)) && (x < (2+period/2)), n = sin(x * pi/period + phase);

if (x > 0) && (x < 4), n = sin(x * pi/4);
else n = 0;
end