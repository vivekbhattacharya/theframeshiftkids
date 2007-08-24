function [n] = exsin(x, dir)

period = 4;
phase1 = pi/2 - 2*pi/period;
phase2 = pi/2 + 2*pi/period;

% Draw a graph of sine at x = 0 to see why we're doing this.
if dir == 117
    if (x > (2-period/2)) && (x < (2+period/2)), n = sin(x * pi/period + phase1);
    else n = 0;
    end
elseif dir == 771
    if (x > (-2-period/2)) && (x < (-2+period/2)), n = sin(x * pi/period + phase2);
    else n = 0;
    end
end