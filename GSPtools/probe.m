function [choice, prob] = probe(x, n1, n2)

k = 1;
pi = 3.14159;

p1 = 1 - (cos(x*pi/4)^2)/(k*n2);
p2 = 1 - (sin(x*pi/4)^2)/(k*n1);

if ((p1 + p2) < 0.666); choice = 0; prob = p1;
elseif (p1 > p2); choice = 0; prob = p1;
elseif (p2 > p1); choice = 1; prob = p2;
elseif (p1 == p2); choice = 0; prob = p1; end