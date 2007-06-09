""" Data retrieval module that mainly calls
    Perl scripts for heavy lifting """

import numpy
def signal(file):
    import os, popen2
    
    file = os.path.abspath(file);
    os.chdir('pearls')
    
    seq_cmd = 'perl.exe getseq.pl "%s"'
    signal_cmd = 'perl.exe free_scan.pl -e -q -p FREIER auuccuccacuag "%s"'
    seq = popen2.popen2(seq_cmd % file)[0].read()
    # Todo: write2fasta
    signal = popen2.popen2(signal_cmd % (file + '.fasta'))[0].readlines()
    # float can't handle empty lines.
    signal = [float(line.strip()) for line in signal if line.strip()]
    signal = numpy.matrix(signal)
    return (signal, seq)

def diff_vectors(mag, phase, count):
    (L, P) = (3, 1)
    for i in range(1, count):
        # God knows what this does.
        x = min(max(1, i - 1), count - 2)
        # Find the codon.
        index = slice(x, x + 2)
        
        # P = 1 implies a linear regression.
        poly_mag = numpy.polyfit(arange(1, L), mag[index], P)
        poly_phase = numpy.polyfit(arange(1, L), phase[index], P)

def displacement(seq, phase, count, diff_vectors, fshifts, bshifts):
    pass

if __name__ == '__main__':
    import sys
    (signal,seq) = signal(sys.argv[1])
    print signal.T