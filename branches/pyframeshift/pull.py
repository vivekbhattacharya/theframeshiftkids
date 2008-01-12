""" Produces the data needed to generate a polar and displacement plot
    representing ribosomal translation and especially frameshifts """

from __future__ import division
import numpy, ghost, os.path, sequence
from nloops import nloops
from Kidnap.scan import scan

def signal(file):
    """ Returns free energy signal and the mRNA sequence """
    file = os.path.abspath(file)
    seq = sequence.get(file)
    signal = scan('auuccuccacuag', file)
    
    # float() can't handle empty lines. numpy.array() introduces
    # twofold speedup.
    signal = numpy.array([line for line in signal])
    return (signal, seq)

from numpy import zeros, mean
from math import atan2, sin, sqrt
def magphase(signal):
    """ Calculates magnitude and phase using the register-memory model
        outlined in the mechanics paper and ignores error """
    # Round signal to a multiple of 3.
    L = len(signal)
    
    limit = L//3
    mag, phase = zeros(limit), zeros(limit)
    registers = zeros(3)
    
    for i in xrange(0, limit):
        # Assume avg_choice in the old code is 0.
        registers += signal[3*i : 3*(i+1)]
        M = registers - mean(registers)
        
        if not M[0] == 0:
            phase[i] = atan2(M[0]*sqrt(3), M[0] + 2*M[1])
            mag[i] = M[0]/sin(phase[i])
        else:
            mag[i], phase[i] = 0, 0
    return mag, phase, limit

from ghost import fakeslope
from cmath import exp, log
def diff_vectors(mag, phase, count):
    vec = []
    for i in xrange(0, count):
        # We want (x-1, x+1) to be (0, 2), (0, 2), (1, 3),
        # (2, 4), ..., (n-4, n-2), (n-3, n-1). Therefore, we
        # need x to be 1, 1, 2, 3, ..., n-3, n-2.
        x = min(max(1, i), count-2)
        
        # Find the codon, remembering that
        # the array index starts at zero.
        magic, phaser = [fakeslope(data[x-1], data[x+1]) for data in (mag, phase)]
        D = exp(1j*phase[i]) * (magic + 1j*mag[i]*phaser)
        
        # angle(z) is equivalent to log(z).imag
        yield [abs(D), log(D).imag]
        
#################################################################

# http://code.google.com/p/theframeshiftkids/wiki/MathBehindTheModel
# Imagine `b`, `calcloops`, and `nudge` in a class except, for
# performance reasons, the `self` object is an array and there is no class.
def b(codon, trig):
    """ Returns (codon, loops, fail probability, trig function) tuple """
    return [codon, calcloops(codon), 1.0, trig]

from math import ceil
def calcloops(codon):
    """ Calculates the window function using tRNA availabilities and
        computed waiting times """
    n = ceil(nloops[codon])
    n = 2. ** (1./n)
    return n/(n-1)

def nudge(self, weight):
    """ Adjusts the (cumulative) probability of failure to choose
        using `calcloops` as the window function of the trig model """
    weight = self[3](weight) ** 10
    self[2] = self[2] * (1 - weight/self[1])
    return 1 - self[2]

from ghost import fxsin, xcos, bxsin
from random import random
from math import sin, pi
def displacement(seq, diffs, fs = []):
    """ Uses differential vectors to stochastically generate a displacement
        plot, detecting frameshifts, reloops, stop codons, and pausing """
    def cheese(self, *args): self.append('%s,%s' % args)
    
    ants, termites = [], []
    species = -30.*(pi/180.)
    x = [0.0, 0.1]
    away = 0
    
    shift = 0; maximus = len(seq); diffs.next()
    for k in xrange(1, maximus//3):
        # 3(k+1) + shift - 1 expanded because i+1:i+4
        # must be the +0 codon, starting with the THIRD.
        i = 3*k + shift + 2
        if i + 5 > maximus: break
        
        codons = [b(seq[i:i+3], bxsin), b(seq[i+1:i+4], xcos), b(seq[i+2:i+5], fxsin)]
        x.append(x[k]); diff = diffs.next()
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
            phi_dx = (pi/3)*x[k+1] - species + diff[1]
            x[k+1] += -0.005 * diff[0] * sin(phi_dx)
    
    # Calculate deviation yield (by parts if necessary)
    x = numpy.array(x)
    away = 0
    if not fs == []:
        # fs, the cumulative actuary, stops early in its tabulation of
        # the predicted shifts e.g. [0 0 0 ... 2] vs. [0 0 0 ... 2 2 2 ...]
        away = sum((x[0:len(fs)] - fs)**2) + sum((x[len(fs):] - fs[-1])**2)
    else: away  = sum((x - 0)**2)
    
    print ants
    print termites
    print sqrt(away / (maximus//3))
    return x

if __name__ == '__main__':
    import sys
    help(sys.modules[__name__])