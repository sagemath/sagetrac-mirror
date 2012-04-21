from sage.all import *

class pAdicLseries(SageObject):
    r"""
    The `p`-adic `L`-series associated to an overconvergent eigensymbol.
    """

    def __init__(self, symb):
        r"""

        INPUT:
            - ``symb`` -- overconvergent eigensymbol
        """
        self._symb = symb

    def prime(self):
        """
        """
        return self._symb.parent().prime()

    def _repr_(self):
        """
        Return print representation.
        """
        s = "%s-adic L-series of $s"%(self._p, self._symb)
        return s

    def series(self, n, quadratic_twist, prec):
        """
        """
        pass

