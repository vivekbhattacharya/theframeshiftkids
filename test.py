from pylab import *
from numpy import min
import pull, ghost

tav = ghost.read('tav.txt')
codons = ghost.read('codons.txt')

def unity(file):
    (signal, seq) = pull.signal(file)
    
    # Synthesize all the data.
    mag, phase, codon_count = pull.magphase(signal)
    diff_vectors = pull.diff_vectors(mag, phase, codon_count)
    places = pull.displacement(seq[12:], codon_count, diff_vectors)
    
    figure(1)
    plot(arange(0, len(places)), places, linewidth=1.0)
    grid(True)
    axis([0, len(places), min(0, min(places)), max(3, max(places))]);
    show()

if __name__ == '__main__':
    import sys
    unity(sys.argv[1])