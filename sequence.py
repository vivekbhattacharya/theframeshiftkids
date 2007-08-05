def read(file):
    f = open(file, 'r')
    a = f.readlines()
    f.close()
    return a

def get(file):
    """ Removes FASTA comments if present """
    lines = read(file)
    if lines[0].find('>') > -1: del lines[0]
    
    lines = [l.strip() for l in lines]
    return sanitize(''.join(lines))

from string import maketrans
def sanitize(seq):
    import re
    seq = re.sub('[\s0-9]', '', seq).lower().strip()
    return seq.translate(maketrans('t', 'u'))
    
def dna2rna(dna):
    return dna.translate(maketrans('tT', 'uU'))

def valid_pair(a,b):
    return (a + b).upper() in ('GC', 'CG', 'AU', 'UA', 'GU', 'UG')

def wc_pair(a,b):
    return (a + b).upper() in ('GC', 'CG', 'AU', 'UA')