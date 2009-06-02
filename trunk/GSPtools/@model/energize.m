% Computes the energy value at displacement x.
function y = energize(self, energy, codon, x)
    shift = get(self, 'shift');
    y = polyval(energy, x - 2*shift);
end
