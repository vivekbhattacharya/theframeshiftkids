% Assumptions: x and y are monotone increasing.
% Linear regression the 5th grade way.
function [m] = fakeslope(x, y, p)
m = (y(end) - y(1))/(x(end) - x(1));