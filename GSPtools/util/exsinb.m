function [n] = exsinb(x)

period = 4;
phase = pi/2 + 2*pi/period;

% Draw a graph of sine at x = 0 to see why we're doing this.
if (x > (-2-period/2)) && (x < (-2+period/2)), n = sin(x * pi/period + phase);
else n = 0;
end