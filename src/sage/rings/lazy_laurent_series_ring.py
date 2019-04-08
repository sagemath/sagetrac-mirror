r"""
Lazy Laurent Series Rings

The ring of lazy Laurent series over a ring has usual arithmetic operations,
but it is actually not a ring in the usual sense since every
arithmetic operation gives a new series.

EXAMPLES:

The definition of Laurent series rings is not initially imported into the
global namespace. You need to import it explicitly to use it::

    sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
    sage: L = LazyLaurentSeriesRing(QQ, 'z')
    sage: L.category()
    Category of magmas and additive magmas
    sage: z = L.gen()
    sage: 1/(1 - z)
    1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
    sage: 1/(1 - z) == 1/(1 - z)
    True

Lazy Laurent series ring over a finite field::

    sage: L = LazyLaurentSeriesRing(GF(3), 'z'); L
    Lazy Laurent Series Ring in z over Finite Field of size 3
    sage: z = L.gen()
    sage: e = 1/(1 + z)
    sage: e.coefficient(100)
    1
    sage: e.coefficient(100).parent()
    Finite Field of size 3

Generating functions of integer sequences are Laurent series over the integer
ring::

    sage: L = LazyLaurentSeriesRing(ZZ, 'z'); L
    Lazy Laurent Series Ring in z over Integer Ring
    sage: z = L.gen()
    sage: 1/(1 - 2*z)^3
    1 + 6*z + 24*z^2 + 80*z^3 + 240*z^4 + 672*z^5 + 1792*z^6 + ...

AUTHORS:

- Kwankyu Lee (2019-02-24): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Kwankyu Lee <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import random

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.categories.magmas_and_additive_magmas import MagmasAndAdditiveMagmas

from sage.misc.cachefunc import cached_method

from .lazy_laurent_series import (LazyLaurentSeries,
                                  LazyLaurentSeriesOperator_gen,
                                  LazyLaurentSeriesOperator_constant)


class LazyLaurentSeriesRing(UniqueRepresentation, Parent):
    """
    Lazy Laurent series ring.

    INPUT:

    - ``base_ring`` -- base ring of this Laurent series ring

    - ``names`` -- name of the generator of this Laurent series ring

    EXAMPLES::

        sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
        sage: LazyLaurentSeriesRing(ZZ, 't')
        Lazy Laurent Series Ring in t over Integer Ring
    """
    Element = LazyLaurentSeries

    def __init__(self, base_ring, names, category=None):
        """
        Initialize.

        TESTS::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 't')
            sage: TestSuite(L).run(skip = ['_test_pickling', '_test_elements'])
        """
        Parent.__init__(self, base=base_ring, names=names,
                        category=MagmasAndAdditiveMagmas().or_subcategory(category))

    def _repr_(self):
        """
        String representation of this Laurent series ring.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: LazyLaurentSeriesRing(GF(2), 'z')
            Lazy Laurent Series Ring in z over Finite Field of size 2
        """
        return "Lazy Laurent Series Ring in %s over %s"%(self.variable_name(), self.base_ring())

    @cached_method
    def gen(self, n=0):
        """
        Return the generator of this Laurent series ring.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.gen()
            z
            sage: L.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: there is only one generator
        """
        if n != 0:
            raise IndexError("there is only one generator")

        op = LazyLaurentSeriesOperator_gen(self)
        c = (self.base_ring().zero(), 2)
        return self.element_class(self, coefficient=op, valuation=1, constant=c)

    def ngens(self):
        """
        Return the number of generators of this Laurent series ring.

        This is always 1.

        EXAMPLES::

            sage: R.<t> = ZZ[[]]
            sage: R.ngens()
            1
        """
        return 1

    def _coerce_map_from_(self, S):
        """
        Return ``True`` if a coercion from ``S`` exists.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L.has_coerce_map_from(ZZ)
            True
            sage: L.has_coerce_map_from(GF(2))
            True
        """
        if self.base_ring().has_coerce_map_from(S):
            return True

        return False

    def _element_constructor_(self, x):
        """
        Construct a Laurent series from ``x``.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L(2)
            0
            sage: L(3)
            1
        """
        R = self.base_ring()

        op = LazyLaurentSeriesOperator_constant(self, R(x))

        return self.element_class(self, coefficient=op, constant=(R.zero(), 1))

    def _an_element_(self):
        """
        Return a Laurent series in this ring.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.an_element()  # random
            z^-10 + z^-9 + z^-8 + z^-7 + z^-6 + z^-5 + z^-4 + z^-3 + z^-2 + z^-1 + 1 + ...
        """
        N = 10

        e = self.base_ring().an_element()

        def r(s, i):
            return self.base_ring().an_element()

        n = random.randint(-N,N)
        m = random.randint(0,N)

        return self.element_class(self, coefficient=r, valuation=n, constant=(e,n+m))

    def series(self, coefficient, valuation, constant=None):
        """
        Return a lazy Laurent series.

        INPUT:

        - ``coefficient`` -- Python function that computes coefficients

        - ``valuation`` -- integer; approximate valuation of the series

        - ``constant`` -- either ``None`` or pair of an element of the base ring and an integer

        Let the coefficient of index `i` mean the coefficient of the term of the
        series with exponent `i`.

        Python function ``coefficient`` returns the value of the coefficient of
        index `i` from input `s` and `i` where `s` is the series itself.

        Let ``valuation`` be `n`. All coefficients of index below `n` are zero.  If
        ``constant`` is ``None``, then the ``coefficient`` function is responsible
        to compute the values of all coefficients of index `\ge n`. If ``constant``
        is a pair `(c,m)`, then the ``coefficient`` function is responsible to
        compute the values of all coefficients of index `\ge n` and `< m` and all
        the coefficients of index `\ge m` is the constant `c`.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.series(lambda s, i: i, 5, (1,10))
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + ...

            sage: def g(s, i):
            ....:     if i < 0:
            ....:         return 1
            ....:     else:
            ....:         return s.coefficient(i - 1) + i
            ....:
            sage: e = L.series(g, -5); e
            z^-5 + z^-4 + z^-3 + z^-2 + z^-1 + 1 + 2*z + ...
            sage: f = e^-1; f
            z^5 - z^6 - z^11 + ...
            sage: f.coefficient(10)
            0
            sage: f.coefficient(20)
            9
            sage: f.coefficient(30)
            -219
        """
        if constant is not None:
            try:
                c,m = constant
            except:
                raise TypeError('not a tuple')

            if valuation > m and c: # weird case
                raise ValueError('inappropriate valuation')

            constant = (self.base_ring()(c), m)

        return self.element_class(self, coefficient=coefficient, valuation=valuation, constant=constant)
