% Linear regression the 5th grade way, assuming the matrix is a filled
% with two points and I have to calculate the slope of the line
% connecting them.
function [m] = fakeslope(ys)
m = ys(2, :) - ys(1, :);