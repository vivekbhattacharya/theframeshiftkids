import numpy, ghost, os, sequence
from nloops import nloops
from Kidnap.scan import scan

def signal(file):
    file = os.path.abspath(file)
    seq = sequence.get(file)
    signal = scan('auuccuccacuag', file)
    
    # float() can't handle empty lines. numpy.array() introduces
    # twofold speedup.
    signal = numpy.array([line for line in signal])
    return (signal, seq)

#################################################################

# I do not compute error (a Boolean value!), because nobody ever uses it.
# See mechanics paper for details of and ideas behind magnitude and
# phase calculations.
from numpy import zeros, mean
from math import atan2, sin, sqrt
def magphase(signal):
    # Round signal to a multiple of 3.
    L = len(signal)
    
    limit = L/3
    mag, phase = zeros(limit), zeros(limit)
    registers = zeros(3)
    
    for i in xrange(0, limit):
        registers += signal[3*i : 3*(i+1)]
        
        # Assume avg_choice in the old code is 0.
        M = registers - mean(registers)
        grass = M[1] + M[2]
        
        if not M[0] == 0:
            phase[i] = atan2(M[0]*sqrt(3), M[0] + 2*M[1])
            mag[i] = M[0]/sin(phase[i])
        elif not grass == 0:
            phase[i] = atan2(grass*sqrt(3), M[2] - M[1])
            mag[i] = -grass/sin(phase[i])
        else:
            mag[i], phase[i] = 0, 0
    return mag, phase, limit

#################################################################

from ghost import fakeslope
from cmath import exp, log
def diff_vectors(mag, phase, count):
    vec = []
    for i in xrange(0, count-1):
        # God knows what this does.
        x = min(max(1, i), count-2)
        
        # Find the codon, remembering that
        # the array index starts at zero.
        magic, phaser = [fakeslope(data[x-1], data[x+1]) for data in (mag, phase)]
        D = exp(1j*phase[i]) * (magic + 1j*mag[i]*phaser)
        
        # angle(z) is equivalent to log(z).imag
        vec += [[abs(D), log(D).imag]]
    return vec
        
#################################################################

# http://code.google.com/p/theframeshiftkids/wiki/MathBehindTheModel
# Imagine `birth`, `calcloops`, and `nudge` in a class except, for
# performance reasons, the `self` object is an array and there is no class.
def b(codon, trig):
    """ Returns [codon, loops, fail probability, trig function] """
    return [codon, calcloops(codon), 1.0, trig]

from math import ceil
def calcloops(codon):
    n = ceil(nloops[codon])
    n = 2. ** (1./n)
    return n/(n-1)

def nudge(self, weight):
    weight = self[3](weight) ** 10
    self[2] = self[2] * (1 - weight/self[1])
    return 1 - self[2]

from ghost import fxsin, xcos, bxsin
from random import random
from math import sin, pi
def displacement(seq, diffs, fshifts=(), bshifts=()):
    def cheese(self, *args): self.append('%s,%s' % args)
    
    ants, termites = [], []
    species = -30.*(pi/180.)
    x = [0.0, 0.1]
    
    shift = 0; maximus = len(seq)
    for k in xrange(1, len(diffs)):
        # 3(k+1) + shift - 1 expanded because i+1:i+4
        # must be the +0 codon, starting with the THIRD.
        i = 3*k + shift + 2
        if i + 5 > maximus: break
        
        codons = [b(seq[i:i+3], bxsin), b(seq[i+1:i+4], xcos), b(seq[i+2:i+5], fxsin)]
        x.append(x[k])
        for persian in xrange(1, 1000):
            # Window function: two units of shift equals on frameshift
            back, here, there = [nudge(c, x[k+1] - 2*shift) for c in codons]
            reloop = 1 - (back + here + there)
            
            # Throw a stone in a rectangle divided into thirds
            # and find which third it landed on via elifs and math.
            r = random()
            if (reloop < here) or (reloop < there) or (reloop < back):
                if (r < here): break
                elif (r < here + there):
                    shift += 1; cheese(ants, codons[1][0], k+1)
                    break
                elif (r < here + there + back):
                    shift -= 1; cheese(termites, codons[1][0], k+1)
                    break
            
            # This follows from phi_signal(1,k) = Dvec(k,2)
            # "A model for +1 frameshifts in eubacteria" by Ponnala, et al.
            phi_dx = (pi/3)*x[k+1] - species + diffs[k][1]
            x[k+1] += -0.005 * diffs[k][0] * sin(phi_dx)
    print ants
    print termites
    return x

if __name__ == '__main__':
    import sys
    help(sys.modules[__name__])