from sage.all import *

class pAdicLseries(SageObject):
    r"""
    The `p`-adic `L`-series associated to an overconvergent eigenform.
    """

    def __init__(self, M, p):
        r"""
        """
        self._M = M
        self._p = p

    def prime(self):
        """
        """
        return self._p

    def _repr_(self):
        """
        Return print representation.
        """
        s = "%s-adic L-series of $s"%(self._p, self._M)
        return s

    def series(self, n, quadratic_twist, prec):
        """
        """
        pass
