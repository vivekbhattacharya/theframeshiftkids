from matplotlib.pylab import *
import pull, ghost

def unity(file):
    signal, seq = pull.signal(file)
    mag, phase, n_codons = pull.magphase(signal)
    diff_vectors = pull.diff_vectors(mag, phase, n_codons)
    
    places = pull.displacement(seq[12:], diff_vectors)
    
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
    
    figure(1); hold(False)
    while True:
        places = pull.displacement(seq, diff_vectors)
        
        plot(arange(0, len(places)), places, linewidth=1.0)
        grid(True)
        axis([0, len(places), min(0, min(places)), max(3, max(places))])
        show()

if __name__ == '__main__':
    import sys
    import cProfile
    cProfile.run('unity(sys.argv[1])', 'foo')
    #megaunity(sys.argv[1])