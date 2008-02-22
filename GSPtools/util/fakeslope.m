% Assumptions: x and y are monotone increasing.
% Linear regression the 5th grade way.
% This assumes L = 3
function [m] = fakeslope(y)
m = (y(end) - y(1))/1; % 2 - 1 = 1