function [n] = exsin(x, dir)

if dir == 117
    if (x > 0) && (x < 4), n = sin(x);
    else n = 0;
    end
elseif dir == 771
    if (x > -4) && (x < 0), n = sin(x);
    else n = 0;
    end
end