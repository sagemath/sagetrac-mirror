from sage.structure.unique_representation import UniqueRepresentation

from .ring import CommutativeRing
from .lazy_laurent_series import LazyLaurentSeries

class LazyLaurentSeriesRing(UniqueRepresentation, CommutativeRing):
    """
    A power series ring.
    """

    def __init__(self, base_ring, names=None, category=None):

        CommutativeRing.__init__(self, base_ring, names=names)

    def _repr_(self):
        """
        Print out a power series ring.

        EXAMPLES::

            sage: R = GF(17)[['y']]
            sage: R
            Laurent Series Ring in y over Finite Field of size 17
            sage: R.__repr__()
            'Laurent Series Ring in y over Finite Field of size 17'
            sage: R.rename('my power series ring')
            sage: R
            my power series ring
        """
        s = "Lazy Laurent Series Ring in %s over %s"%(self.variable_name(), self.base_ring())
        return s

    def gen(self, n=0):
        """
        Return the generator of this power series ring.

        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: R.gen()
            t
            sage: R.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: generator n>0 not defined
        """
        if n != 0:
            raise IndexError("generator n > 0 not defined")

        one = self.base_ring().one()
        zero = self.base_ring().zero()

        def f(s, m):
            return one if m == 1 else zero

        return LazyLaurentSeries(self, coefficient=f, approximate_valuation=1, constant=(zero,2))

    def ngens(self):
        """
        Return the number of generators of this power series ring.

        This is always 1.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: R.ngens()
            1
        """
        return 1

    def _coerce_map_from_(self, S):
        """
        A coercion from `S` exists, if `S` coerces into ``self``'s base ring,
        or if `S` is a univariate polynomial or power series ring with the
        same variable name as self, defined over a base ring that coerces into
        ``self``'s base ring.

        EXAMPLES::

            sage: A = GF(17)[['x']]
            sage: A.has_coerce_map_from(ZZ)  # indirect doctest
            True
            sage: A.has_coerce_map_from(ZZ['x'])
            True
            sage: A.has_coerce_map_from(ZZ['y'])
            False
            sage: A.has_coerce_map_from(ZZ[['x']])
            True

        """
        if self.base_ring().has_coerce_map_from(S):
            return True

        return False

    def _element_constructor_(self, x):
        """
        Coerce object to this power series ring.

        Returns a new instance unless the parent of f is self, in which
        case f is returned (since f is immutable).

        INPUT:


        -  ``f`` - object, e.g., a power series ring element

        -  ``prec`` - (default: infinity); truncation precision
           for coercion

        -  ``check`` - bool (default: True), whether to verify
           that the coefficients, etc., coerce in correctly.


        EXAMPLES::

            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: R(t+O(t^5))    # indirect doctest

        """
        zero = self.base_ring().zero()
        def f(s, m):
            return self.base_ring()(x) if m == 0 else zero

        return LazyLaurentSeries(self, coefficient=f, constant=(zero,1))



