"""
Multiple zeta values
"""
from __future__ import print_function

from sage.symbolic.function import BuiltinFunction
from sage.rings.all import Integer
from sage.misc.latex import latex
from sage.symbolic.ring import SR
from sage.libs.all import pari


class MultipleZetaValue(BuiltinFunction):
    def __init__(self):
        """
        Initialize class.

        EXAMPLES::

            sage: maxima(hypergeometric)
            hypergeometric
        """
        BuiltinFunction.__init__(self, 'multizeta', nargs=1,
                                 conversions={'pari': 'zetamult'})

    def __call__(self, a):
        r"""
        Return the multiple zeta value of `n_1,n_2,...,n_k`.

        INPUT:

        - ``a`` -- a list or tuple of parameters

        EXAMPLES::

            sage: from sage.functions.multizeta import multizeta
            sage: multizeta((4,2))
            multizeta((4, 2))
        """
        return BuiltinFunction.__call__(self, SR._force_pyobject(a))

    def _print_latex_(self, a):
        r"""
        TESTS::

            sage: from sage.functions.multizeta import multizeta
            sage: latex(multizeta([2, 2]))
            \zeta(2,2)
        """
        aa = ",".join(latex(c) for c in a)
        return r"\zeta(" + aa + ")"

    def _eval_(self, a):
        """
        Evaluate the multiple zeta value.

        EXAMPLES::

            sage: from sage.functions.multizeta import multizeta
            sage: multizeta((3,3,2))
            multizeta((3, 3, 2))
        """
        if not isinstance(a, tuple):
            raise TypeError("The parameters must be of type tuple")
        if not all(isinstance(y, Integer) for y in a):
            raise TypeError('arguments must be integers')
        if not all(y >= 1 for y in a):
            raise TypeError('arguments must be positive integers')
        if a[0] == 1:
            raise ValueError('divergence')

    def _evalf_(self, a, parent, algorithm=None):
        """
        TESTS::

            sage: from sage.functions.multizeta import multizeta
            sage: M = multizeta((3,))
            sage: pari(M)
            1.20205690315959
            sage: M = multizeta((2,1))
            sage: pari(M)
            1.20205690315959
        """
        if not isinstance(a, tuple):
            raise TypeError("The first parameters must be of type tuple")

        if algorithm is None:
            algorithm = "pari"
        if algorithm == "pari":
            return pari(a).zetamult()

        raise ValueError("unknown algorithm")


multizeta = MultipleZetaValue()
