function [n] = mysin(x)

if x < 0
    n = 0;
else n = sin(x);
end