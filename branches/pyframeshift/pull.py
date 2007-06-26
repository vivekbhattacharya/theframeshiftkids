""" Data retrieval module that mainly calls
    Perl scripts for heavy lifting """

import numpy

seq_cmd = r'perl.exe getseq.pl "%s"'
signal_cmd = r'perl.exe free_scan.pl -e -q -p FREIER auuccuccacuag "%s.fasta"'
def signal(file):
    import os, ghost    
    file = os.path.abspath(file)
    os.chdir('pearls')
    
    # Todo: write2fasta
    seq = ghost.steal_out(seq_cmd % file).read()
    signal = ghost.steal_out(signal_cmd % file).readlines()
    
    # float can't handle empty lines.
    signal = [float(line.strip()) for line in signal if line.strip()]
    signal = numpy.matrix(signal)
    return (signal, seq)

from numpy import polyfit, polyder, polyval
class PolarPoint(object):
    def __init__(num, index):
        L, P = 3, 1
        self.num = num
        # P = 1 implies a linear regression.
        self.poly = polyfit(range(1,L), num, P)
        self.der = polyval(polyder(self.poly), index)

from numpy import exp, complex, angle, abs
def diff_vectors(mag, phase, count):
    for i in range(0, count-1):
        # God knows what this does.
        x = min(max(1, i), count-2)
        # Find the codon, remembering that
        # the array index starts at zero.
        index = slice(x - 1, x + 2)
        magic, phaser = [PolarPoint(data[index], i+1) for data in (mag, phase)]
        D = exp(1j*phase[i]) * (magic.der + 1j*mag[i]*phaser.der)
        

# http://code.google.com/p/theframeshiftkids/wiki/MathBehindTheModel
class Codon(object):
    def __init__(self, codon):
        self.codon = codon
        self.loops = Codon.calcloops(codon)
        self.fail = 1.
    
    def calcloops(codon):
        from nloops import *
        n = ceil(nloops(codon))
        n = pow(2., 1./n)
        return n/(n-1)
    
    def nudge(weight):
        self.fail = self.fail * (1 - weight/self.loops)
        return 1 - self.fail

from ghost import fxsin, xcos, bxsin
from random import random
from math import sin, pi
def displacement(seq, phase, count, diffs, fshifts=(), bshifts=()):
    ants, termites = [], []
    species = -30*(pi/180)
    x, C1 = [0 0.1], 0.005
    
    def cheese(self, *args): self.append('%s,%s' % args)
    def delphi(x0): (pi/3)*x0 - species
    for k in xrange(2, count - 1):
        i = 3*k
        if i + 4 > len(seq): break
        codons = [Codon(j) for j in (seq[i:i+3], seq[i+1:i+4], seq[i+2:i+5])]
        
        phi_signal[k-1] = diffs[k,1]; x0 = x[k-1]
        phi_dx = delphi(x0)
        for persian in xrange(1, 1000):
            a = x0*pi/4;
            # Window function
            weights = [func(a) for func in (fxsin, xcos, bxsin)]
            back, here, there = [c.nudge(w) for (c, w) in zip(codons, weights)]
            reloop = 1 - (back + here + there)
            
            r = random() # Mersenne Twister
            if all([reloop < p for p in (back, here, there)]):
                if r < here: break
                else if r < (here + there):
                    shift += 1; cheese(ants, codons[1].codon, k)
                    break
                else if r < (here + there + back):
                    shift -= 1; cheese(termites, codons[1].codon, k)
                    break
            
            phi_dx = delphi(x0)
            dx = -C1 * diffs[k-1,0] * sin(phi_signal[k-1] + phi_dx)
            x0 += dx
    x[k] = x0

if __name__ == '__main__':
    import sys
    (signal,seq) = signal(sys.argv[1])
    print signal.T
    print seq