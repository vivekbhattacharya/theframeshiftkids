% Like `mean` and `std`, `exmean` computes mean
% for each column of a matrix but ignores the value
% of -3.14 as if it never existed, properly modifying
% the final mean and standard deviation to account
% for a smaller set of values for that particular row.
% This is useful when a program like jejunity needs
% the average of a jagged matrix but can also stuff
% that matrix into a normal matrix, filling in -3.14
% for non-defined values.
function [avg, stddev] = exmean(matrix)

% Matrix of zeros and ones, zeros corresponding
% to where -3.14s are.
neg_board = ~(matrix + 3.14);
totals = sum(~neg_board);

% Convert -3.14 to zeros.
actual = matrix + neg_board*3.14;
avg = sum(actual) ./ totals;

% E(X^2) - E(X)^2
stddev = (sum(actual .^ 2) ./ totals) - avg .^ 2;