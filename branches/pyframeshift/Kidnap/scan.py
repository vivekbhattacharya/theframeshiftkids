import sys
sys.path = ['..'] + sys.path

from bind import Bind
from sequence import dna2rna, get

def scan(rna, seq_file):
    """ Yields a generator of free energy. How exciting. """
    
    bind = Bind()
    rna = dna2rna(rna).lower()
    seq = dna2rna(get(seq_file)).lower()
    
    n_seq = len(seq); n_rna = len(rna)
    for i in xrange(0, n_seq - n_rna + 1):
        segment = seq[i:i + n_rna]
        
        # Free energy values greater than zero represent binding
        # that would cannot take place without added energy,
        # equivalent to as if no binding had taken place
        e = bind.energy(segment, rna)
        if e > 0: e = 0
        yield e