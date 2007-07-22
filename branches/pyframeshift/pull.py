""" Data retrieval module that mainly calls
    Perl scripts for heavy lifting """

import numpy, ghost, os
from nloops import nloops
from tempfile import mkstemp as tempfile

_, fasta = tempfile()
seq_cmd = r'perl.exe getseq.pl "%s" "%s"'
signal_cmd = r'perl.exe free_scan.pl -e -q -p FREIER auuccuccacuag "%s"'
def signal(file):
    file = os.path.abspath(file)
    os.chdir('pearls')
    
    # Todo: write2fasta
    seq = ghost.steal_out(seq_cmd % (file, fasta)).read()
    signal = ghost.steal_out(signal_cmd % fasta).readlines()
    
    # float can't handle empty lines.
    signal = [float(line.strip()) for line in signal if line.strip()]
    return (signal, seq)

#################################################################

# I do not compute error, because nobody ever uses it.
from numpy import zeros, arctan2, sin, sqrt, mean, round
def magphase(signal):
    L = len(signal)
    # Round signal to a multiple of 3.
    if L % 3: signal = signal[0:L - (L % 3)]
    
    limit = L/3
    Mag, Phase = zeros(limit), zeros(limit)
    for i in xrange(0, limit):
        # +1 comes from Python's indexing, which starts at zero,
        # that conflates with `i` in that it's also used for
        # creating a sequence of endpoints here.
        #
        # We want 0:3, 0:6, 0:9, 0:12, ...
        y = signal[0:((i+1)*3)]
        M = zeros(3)
        for j in xrange(0, len(y), 3): M[0:3] += y[j:j+3]
        # Assume avg_choice in the old code is 0.
        M -= mean(M)
        
        grass = M[1] + M[2]
        if not M[0] == 0:
            Phase[i] = arctan2(M[0]*sqrt(3), M[0] + 2*M[1])
            Mag[i] = M[0]/sin(Phase[i])
        elif not grass == 0:
            Phase[i] = arctan2(grass*sqrt(3), M[2] - M[1])
            Mag[i] = -grass/sin(Phase[i])
        else: Mag[i], Phase[i] = 0, 0
    return Mag, Phase, len(Mag)

#################################################################

from numpy import polyfit, polyder, polyval
class PolarPoint(object):
    def __init__(self, num, index):
        L, P = 3, 1
        self.num = num
        # P = 1 implies a linear regression.
        self.poly = polyfit(range(1,L+1), num, P)
        self.der = polyval(polyder(self.poly), index)

from numpy import exp, complex, angle, abs
def diff_vectors(mag, phase, count):
    vec = []
    for i in range(0, count-1):
        # God knows what this does.
        x = min(max(1, i), count-2)
        # Find the codon, remembering that
        # the array index starts at zero.
        index = slice(x - 1, x + 2)
        magic, phaser = [PolarPoint(data[index], i+1) for data in (mag, phase)]
        D = exp(1j*phase[i]) * (magic.der + 1j*mag[i]*phaser.der)
        
        vec += [[abs(D), angle(D)]]
    return vec
        
#################################################################

# http://code.google.com/p/theframeshiftkids/wiki/MathBehindTheModel
from numpy import ceil, round
class Codon(object):
    def __init__(self, codon):
        self.codon = codon
        self.loops = Codon.calcloops(codon)
        self.fail = 1.0
    
    @staticmethod
    def calcloops(codon):
        # Implement nloops
        n = ceil(nloops[codon])
        n = 2. ** (1./n)
        return n/(n-1)
    
    def nudge(self, weight):
        print round([self.fail, weight, self.loops], 2)
        self.fail = self.fail * (1 - weight/self.loops)
        return 1 - self.fail

from ghost import fxsin, xcos, bxsin
from random import random
from numpy import sin, pi
def displacement(seq, phase, count, diffs, fshifts=(), bshifts=()):
    ants, termites = [], []
    species = -30.*(pi/180.)
    x, C1 = [0.0, 0.1], 0.005
    
    def cheese(self, *args): self.append('%s,%s' % args)
    def delphi(x0): return (pi/3)*x0 - species
    phi_signal = [0]
    shift = 0
    for k in xrange(1, 3):
        i = 3*k + 2 + shift
        if i + 4 > len(seq): break
        codons = [Codon(j) for j in (seq[i:i+3], seq[i+1:i+4], seq[i+2:i+5])]
        
        phi_signal += [diffs[k][1]]
        x0 = x[k]
        phi_dx = delphi(x0)
        for persian in xrange(1, 1000):
            a = x0*pi/(4.0)
            # Window function
            weights = [func(a)**10 for func in (fxsin, xcos, bxsin)]
            back, here, there = [c.nudge(w) for (c, w) in zip(codons, weights)]
            reloop = 1 - (back + here + there)
            
            # r = random() # Mersenne Twister
            if (reloop < here) or (reloop < there) or (reloop < back):
                if (here > there) and (here > back): break
                elif (there > here) and (there > back):
                    shift += 1; cheese(ants, codons[1].codon, k)
                    break
                elif (back > here) and (back > there):
                    shift -= 1; cheese(termites, codons[1].codon, k)
                    break
            
            
            phi_dx = delphi(x0)
            dx = -C1 * diffs[k][0] * sin(phi_signal[k] + phi_dx)
            x0 += dx
        x += [x0]
    return x

if __name__ == '__main__':
    import sys
    (signal,seq) = signal(sys.argv[1])
    print signal.T
    print seq