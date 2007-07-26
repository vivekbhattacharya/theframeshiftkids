from pylab import *
from numpy import min
import pull, ghost

def unity(file):
    (signal, seq) = pull.signal(file)
    mag, phase, codon_count = pull.magphase(signal)
    diff_vectors = pull.diff_vectors(mag, phase, codon_count)
    
    places = pull.displacement(seq[12:], diff_vectors)
    
    return
    figure(1)
    plot(arange(0, len(places)), places, linewidth=1.0)
    grid(True)
    axis([0, len(places), min(0, min(places)), max(3, max(places))])
    show()

def megaunity(file):
    (signal, seq) = pull.signal(file)
    mag, phase, codon_count = pull.magphase(signal)
    diff_vectors = pull.diff_vectors(mag, phase, codon_count)
    seq = seq[12:]
    
    while True:
        places = pull.displacement(seq, diff_vectors)

if __name__ == '__main__':
    import sys
    import cProfile
    cProfile.run('unity(sys.argv[1])', 'foo')
    # unity(sys.argv[1])