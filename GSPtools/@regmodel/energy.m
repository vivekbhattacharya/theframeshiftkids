% Does nothing except assign the codon number to the energy because
% cumulative energy and instanteous energy are both calculated once
% during construction.
function [energy] = energy(self, codon)
    energy = codon;
end
