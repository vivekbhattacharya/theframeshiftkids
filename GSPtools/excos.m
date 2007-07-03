function [n] = excos(x)

period = 4.5;
x = x * pi/period;

if (x > -period/2) && (x < period/2), n = cos(x);
else n = 0;
end