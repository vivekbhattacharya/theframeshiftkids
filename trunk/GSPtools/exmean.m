function [avg,stddev] = exmean(matrix)
avg = zeros(1, size(matrix, 2));
stddev = avg;
for i=1:size(matrix,2)
	total = matrix(:,i);
	sum = [0 0];
	for j=1:size(total,1)
		if total(j) ~= -3.14
			sum(1) = sum(1) + total(j);
			sum(2) = sum(2) + 1;
		end
    end
    
    if sum(2) == 0
        avg(i) = 0;
        stddev(i) = 0;
        continue;
    end
	avg(i) = sum(1)/sum(2);
	
	% s.x = E(x^2) - E(x)^2
	% `sum` is the sum of squares.
	sum = [0 0];
	for j=1:size(total,1)
		if total(j) ~= -3.14
			sum(1) = sum(1) + total(j)^2;
			sum(2) = sum(2) + 1;
		end
	end
	stddev(i) = sum(1)/sum(2) - avg(i)^2;
end