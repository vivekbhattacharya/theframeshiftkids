function [choice, prob] = probe(x, n1, n2)

k = 0.4;
n=2;

%x = x - 2*floor(abs(x/2));
disp(x);

% p1 = 1 - (sin(x*pi/4-k*n1)^(2*n));%*(k*n2);
% p2 = 1 - (cos(x*pi/4-k*n2)^(2*n));%*(k*n1);
s1 = (1/n1)^k;
s2 = (1/n2)^k;

p1 = exp((-(x-0)^2)/(2*s1^2));
p2 = exp((-(x-2)^2)/(2*s2^2));

p3 = 1 - (p1+p2);
p = [p1 p2 p3];
disp(p);

if (max(p)==p(1,3))
	choice = 2;
    prob = p3;
elseif (max(p)==p(1,1)); choice = 0; prob = p1;
elseif (max(p)==p(1,2)); choice = 1; prob = p2;

end