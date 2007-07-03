def read(file):
    return open(file, 'r').readlines()

def steal_out(cmd):
    from subprocess import Popen, PIPE
    return Popen(cmd, stdout=PIPE).stdout

from numpy import pi, sin, cos
period = 4.0
def fxsin(x):
    x *= pi/period
    phase = pi/2 - 2*pi/period
    
    if x > (2-period/2) and x < (2+period/2):
        return sin(x + phase)
    return 0.0

def xcos(x):
    x *= pi/period
    
    if x > (-period/2.0) and x < (period/2.0):
        return cos(x)
    return 0

def bxsin(x):
    x *= pi/period
    phase = pi/2.0 + 2.0*pi/period
    
    if x > (-2.0-period/2.0) and x < (-2.0+period/2.0):
        return sin(x + phase)
    return 0.0
