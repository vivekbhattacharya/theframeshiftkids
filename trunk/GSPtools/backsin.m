function [n] = backsin(x)

if (x > -4) && (x < 0)
    n = sin(x);
else n = 0;
end