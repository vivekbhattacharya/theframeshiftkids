function [n] = excos(x)

if (x > -2) && (x < 2), n = cos(x);
else n = 0;
end