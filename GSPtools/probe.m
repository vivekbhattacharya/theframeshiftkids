function [choice, prob] = probe(x, this, next)
% Normal distribution to calculate probabilities.
% By The Frameshift Kids.
% Parameters: `this` is the current codon.
%   `next` is the next codon.
%
% Return values: `choice` is 0, 1, or 2, representing whether
%   the probability of p1, p2, or failure is highest.
%   `prob` is the corresponding probability.
k = exp(-1);
n1 = nloopcalcify(this); n2 = nloopcalcify(next);

% p1 = 1 - (sin(x*pi/4-k*n1)^(2*n));%*(k*n2);
% p2 = 1 - (cos(x*pi/4-k*n2)^(2*n));%*(k*n1);
s1 = (1/n1)^k;
s2 = (1/n2)^k;

p1 = exp((-(x-0)^2)/(2*s1^2));
p2 = exp((-(x-2)^2)/(2*s2^2));
p3 = 1 - (p1+p2);

if (p3 > p2 && p3 > p1); choice = 2; prob = p3;
elseif (p1 > p2); choice = 0; prob = p1;
elseif (p2 > p1); choice = 1; prob = p2;
end