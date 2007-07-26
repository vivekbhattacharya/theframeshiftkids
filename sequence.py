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
    return ''.join(lines)