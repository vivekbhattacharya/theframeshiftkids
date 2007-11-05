function [n] = exsinb(x)

% period = 4;
% phase = pi/2 + 2*pi/period; => pi

% Draw a graph of sine at x = 0 to see why we're doing this.
% if (x > (-2-period/2)) && (x < (-2+period/2)), n = sin(x * pi/period + phase);

if (x > -4 && x < 0), n = sin(x * pi/4 + pi);
else n = 0;
end