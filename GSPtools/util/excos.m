function [n] = excos(x)
period = 4;

% Draw a graph of cosine at x = 0 to see why we're doing this.
if (x > -period/2) && (x < period/2), n = cos(x * pi/period);
else n = 0;
end