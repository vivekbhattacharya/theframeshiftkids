from pylab import *
import pull, ghost

tav = ghost.read('tav.txt')
codons = ghost.read('codons.txt')

def unity(file):
    (signal, seq) = pull.signal(file)
    
    # Synthesize all the data.
    (mag, phase, codon_count) = pull.polarity(signal)
    diff_vectors = pull.diff_vectors(mag, phase, codon_count)
    (thetas, placements, forces) = pull.displacement(seq[13:], \
        phase, codon_count, diff_vectors, [], [])
    
    figure(1);
    plot(arange(1, len(placements)), placements)

if __name__ == '__main__':
    import sys
    unity(sys.argv[1])