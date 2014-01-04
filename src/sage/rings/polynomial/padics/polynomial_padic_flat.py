"""
`p`-adic polynomials with a flat precision

This module provides classes for `p`-adic polynomials with a flat precision,
i.e., the absolute precision is the same for all coefficients of the
polynomial.

AUTHORS:

- Julian Rueth (2013- 09-12): added some docstrings and cleanup of the inheritance
"""

#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                     2008 John Cremona <john.cremona@gmail.com>
#                     2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import sage.rings.padics.misc
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.libs.all import pari_gen
from sage.structure.factorization import Factorization
from polynomial_padic_generic import Polynomial_padic_generic_ring

class Polynomial_padic_flat(Polynomial_padic_generic_ring):
    r"""
    A polynomial with a flat precision over a `p`-adic base ring.

    INPUT:

    - ``parent`` -- a polynomial ring over a `p`-adic ring

    - ``x`` -- data to construct a polynomial (default: ``None``)

    - ``check`` -- a boolean (default: ``True``), whether to make sure that all
      coefficients are in the `p`-adic base ring

    - ``is_gen`` -- a boolean (default: ``False``), whether this is the
      generator of the polynomial ring

    - ``construct`` -- a boolean (default: ``False``), unused

    - ``absprec`` -- an integer or ``None`` (default: ``None``), if given,
      every coefficient is capped to this precision

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.polynomial_padic_flat import Polynomial_padic_flat
        sage: R.<x> = ZpFM(3)[] # indirect doctest
        sage: isinstance(x, Polynomial_padic_flat)
        True

    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False, absprec=None):
        """
        Initialization.

        TESTS::

            sage: R.<x> = ZpFM(3)[]
            sage: type(x)
            <class 'sage.rings.polynomial.padics.polynomial_padic_flat.Polynomial_padic_flat'>

        """
        if x is None:
            Polynomial_padic_generic_ring.__init__(self, parent, x = None, is_gen = is_gen)
            return
        R = parent.base_ring()
        if sage.rings.fraction_field_element.is_FractionFieldElement(x):
            if x.denominator() != 1:
                raise TypeError, "denominator must be 1"
            else:
                x = x.numerator()
        from sage.rings.polynomial.polynomial_element import Polynomial
        if isinstance(x, Polynomial):
            if x.base_ring() is R:
                x = list(x.list())
            else:
                x = [R(a) for a in x.list()]
        elif isinstance(x, list):
            if check:
                x = [R(a) for a in x]
        elif isinstance(x, dict):
            if check:
                m = infinity
                zero = R(0)
                n = max(x.keys())
                v = [zero for _ in xrange(n+1)]
                for i, z in x.iteritems():
                    v[i] = R(z)
                    m = min(m, v[i].precision_absolute())
                x = v
            else:
                m = sage.rings.padics.misc.min(a.precision_absolute() for a in x.values())
            if not absprec is None:
                m = min(m, absprec)
            Polynomial_padic_generic_ring.__init__(self, parent, x, absprec = m)
            return
        elif isinstance(x, pari_gen):
            x = [R(w) for w in x.list()]
        else:
            x = [R(x)]
        if absprec is None:
            absprec = infinity
        m = min([a.precision_absolute() for a in x] + [absprec])
        Polynomial_padic_generic_ring.__init__(self, parent, x)
