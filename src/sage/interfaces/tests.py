from all import *
from sage.misc.misc import cputime, walltime
import sys

def manyvars(s, num=70000, inlen=1, step=2000):
    """
    Test that > 65,000 variable names works in each system.
    """
    print "Testing -- %s"%s
    t = '"%s"'%('9'*int(inlen))
    try:
        t = cputime()
        w = walltime()
        v = []
        for i in range(num):
            if i%step==0:
                sys.stdout.write('%s '%i)
                sys.stdout.flush()
            v.append(s(t))
        print '\nsuccess -- time = cpu: %s, wall: %s'%(cputime(t), walltime(w))
    except:
        print "%s -- failed!"%s

def manyvars_all(num=70000):
    #for s in [gap, gp, singular, kash, magma, octave, maxima, mathematica]:
    for s in [kash, magma, octave, maxima, mathematica]:
        manyvars(s, num)

# bad: maple -- infinite loop -- exception pexpect.EOF: <pexpect.EOF instance at 0xb091250c> in  ignored

def manyvars_all2(num=70000):
    for s in [singular, maxima, mathematica, octave]:
        manyvars(s, num)
