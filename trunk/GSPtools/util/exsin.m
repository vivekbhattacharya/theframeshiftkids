function [n] = exsin(x, dir)

period = 4;
x = x * pi/period;
phase1 = pi/2 - 2*pi/period;
phase2 = pi/2 + 2*pi/period;

if dir == 117
    if (x > (2-period/2)) && (x < (2+period/2)), n = sin(x + phase1);
    else n = 0;
    end
elseif dir == 771
    if (x > (-2-period/2)) && (x < (-2+period/2)), n = sin(x + phase2);
    else n = 0;
    end
end