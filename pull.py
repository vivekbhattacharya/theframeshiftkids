""" Data retrieval module that mainly calls
    Perl scripts for heavy lifting """

import numpy
def signal(file):
    import os
    seq = get_seq(file)
    signal = os.system('perl get_signal.pl')
    return (signal, seq)

def seq(file):
    return os.system('perl get_seq.pl')

def diff_vectors(mag, phase, count):
    (L, P) = (3, 1)
    for i in range(1, count):
        # God knows what this does.
        x = min(max(1, i - 1), count - 2)
        # Find the codon.
        index = slice(x, x + 2)
        
        # P = 1 implies a linear regression.
        poly_mag = numpy.polyfit(arange(1, L), Mag[index], P)
        poly_phase = numpy.polyfit(arange(1, L), Phase[index], P)

def displacement(seq, phase, count, diff_vectors, fshifts, bshifts):
    pass