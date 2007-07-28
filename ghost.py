def steal_out(cmd):
    from subprocess import Popen, PIPE
    return Popen(cmd, stdout=PIPE).stdout

def fakeslope(y):
    return (y[-1] - y[0])/2;

from numpy import pi, sin, cos
def fxsin(x):
    """ A piecewise function that computes the window
        function on +1 frameshifts """
    period = 4
    magic = 2*pi/period
    phase = pi/2 - magic
    
    if x > 4 - period and x < 4 + period:
        return sin(x*magic/2 + phase)
    return 0.0

def xcos(x):
    """ A piecewise function that computes the window
        function on the current frame """
    period = 4
    magic = 2*pi/period

    if x > -period and x < period:
        return cos(x*magic/2)
    return 0

def bxsin(x):
    """ A piecewise function that computes the window
        function on -1 frameshifts """
    period = 4
    magic = 2*pi/period
    phase = pi/2 + magic
    
    if x > -4 - period and x < -4 + period:
        return sin(x*magic/2 + phase)
    return 0.0
