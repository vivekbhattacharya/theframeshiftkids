function [dy, dx] = diffvec(m1,m2,p1,p2)
    % m1, p1 is the first vector in polar; m2, p2 the second
    dx = m2 .* cos(p2) - m1 .* cos(p1);
    dy = m2 .* sin(p2) - m1 .* sin(p1);
end
