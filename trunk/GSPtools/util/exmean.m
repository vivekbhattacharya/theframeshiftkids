% Compute the mean of a matrix along its
% columns, ignoring any value equal to
% -3.14, which represents that there was
% no x value at that row's (iteration's)
% column. Thus, exmean helps compute the
% mean of several arrays that are not of the
% same length.
%
% Also, it spits out the standard deviation.
function [avg, stddev] = exmean(matrix)
avg = zeros(1, size(matrix, 2));
stddev = avg;
for i=1:size(matrix,2)
	total = matrix(:,i); sum = [0 0];
	for j=1:size(total,1)
		if total(j) == -3.14, continue; end;
        sum(1) = sum(1) + total(j);
        sum(2) = sum(2) + 1;
    end
    
    % All the values were -3.14.
    if sum(2) == 0
        avg(i) = 0; stddev(i) = 0;
        continue;
    end
	avg(i) = sum(1)/sum(2);
	
    % Now calculate the standard deviation
    % according to "s.x = E(x^2) - E(x)^2."
    % `sum` is the sum of squares now.
    sum = [0 0];
    for j=1:size(total,1)
        if total(j) == -3.14, continue; end;
        sum(1) = sum(1) + total(j)^2;
        sum(2) = sum(2) + 1;
    end
    stddev(i) = sum(1)/sum(2) - avg(i)^2;
end