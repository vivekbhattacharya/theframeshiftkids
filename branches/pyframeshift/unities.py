from matplotlib.pylab import *
import pull, ghost

def unity(file, fs = []):
    displacement = walrus_surprise(file, fs)
    places = displacement()

    return
    figure(1)
    plot(range(0, len(places)), places, linewidth=1.0)
    grid(True)
    axis([0, len(places), min(0, min(places)), max(3, max(places))])
    show()

def megaunity(file, fs = []):
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

    actuary = []
    if not fs == []:
        fs = eval(fs)
        actuary = [0] * (fs[-1] + 1)
        for i in fs: actuary[i] = 2
        actuary = array(actuary).cumsum()

    signal, seq = signal(file)
    mag, phase, n_codons = magphase(signal)
    diff_vectors = diff_vectors(mag, phase, n_codons)

    def walrus():
        return displacement(seq[12:], diff_vectors, actuary)
    return walrus

if __name__ == '__main__':
    import sys
    import cProfile

    cProfile.run('unity(*sys.argv[1:])', 'foo')
    #megaunity(sys.argv[1])
