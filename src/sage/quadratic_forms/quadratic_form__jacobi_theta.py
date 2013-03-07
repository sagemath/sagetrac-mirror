"""
Jacobi Theta Series of Quadratic Forms

Andrew Fiori

"""

from sage.modular.jacobi.jacobi_theta import __jacobi_theta
from sage.misc.misc import cputime, verbose

def bilinear_form(self,x,y):
    return self(x+y)-self(x)-self(y)

def jacobi_theta(self, Vec, Max=10):
    """
    Compute the theta series as a power series in the variable given
    in var_str (which defaults to '`q`'), up to the specified precision
    `O(q^max)`.

    """
    return __jacobi_theta(self,Vec,Max)



