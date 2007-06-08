function [n] = exsin(x)

if (x > 0) && (x < 4)
    n = sin(x);
else n = 0;
end