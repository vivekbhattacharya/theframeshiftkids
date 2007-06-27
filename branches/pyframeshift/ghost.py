def read(file):
    return open(file, 'r').readlines()

def steal_out(cmd):
    from subprocess import Popen, PIPE
    return Popen(cmd, stdout=PIPE).stdout

def fxsin(): pass
def xcos(): pass
def bxsin(): pass