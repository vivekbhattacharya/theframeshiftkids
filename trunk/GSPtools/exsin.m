function [n] = exsin(x, dir)

if strcmp(dir, '>')
    if (x > 0) && (x < 4), n = sin(x);
    else n = 0;
    end
elseif strcmp(dir, '<')
    if (x > -4) && (x < 0), n = sin(x);
    else n = 0;
    end
end