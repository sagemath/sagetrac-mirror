r"""
Common base class for all polynomials over `p`-adic rings

AUTHORS:

- Julian Rueth (2013-09-12): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_domain

class Polynomial_padic(Polynomial_generic_domain):
    r"""
    A polynomial over a `p`-adic ring.

    INPUT:

    - ``parent`` -- a polynomial ring over a `p`-adic ring

    - ``is_gen`` -- whether this is the generator of the polynomial ring
      (default: ``False``)

    .. NOTE::

        In contrast to :class:`polynomial_padic.Polynomial_padic_generic`
        (which inherits from this class), this class is meant as a base class
        for implementations which provide their own handling of the polynomial
        data.

    EXAMPLES::

        sage: R.<x> = Zp(3)[] # indirect doctest
        sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
        sage: isinstance(x, Polynomial_padic)
        True

    """
    def __init__(self, parent, is_gen=False):
        r"""
        Initialization.

        TESTS::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R.<x> = Zp(3)[]
            sage: type(Polynomial_padic(R))
            <class 'sage.rings.polynomial.padics.polynomial_padic.Polynomial_padic'>

        """
        Polynomial_generic_domain.__init__(self, parent, is_gen=is_gen)
