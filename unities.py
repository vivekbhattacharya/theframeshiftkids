from matplotlib.pylab import *
import pull, ghost

def unity(file, fs = None):
    displacement = walrus_surprise(file, fs)
    places = displacement()

    figure(1)
    plot(range(0, len(places)), places, linewidth=1.0)
    grid(True)
    axis([0, len(places), min(0, min(places)), max(3, max(places))])
    show()

def megaunity(file, fs = None):
    d = walrus_surprise(file, fs)

    figure(1); hold(False)
    while True:
        places = d()

        plot(arange(0, len(places)), places, linewidth=1.0)
        grid(True)
        axis([0, len(places), min(0, min(places)), max(3, max(places))])
        show()

def walrus_surprise(file, fs):
    from pull import signal, magphase, diff_vectors, displacement
    from numpy import array

    signal, seq = signal(file)
    mag, phase, n_codons = magphase(signal)
    diff_vectors = diff_vectors(mag, phase, n_codons)

    def walrus():
        return displacement(seq[12:], diff_vectors, fs)
    return walrus

if __name__ == '__main__':
    unity('c:/up/prfB.txt', 25)
