r"""
Common base classes for all polynomial rings over `p`-adic rings

AUTHORS:

- Julian Rueth (2013-09-26): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import sage.rings.polynomial.polynomial_ring

class PolynomialRing_padic(sage.rings.polynomial.polynomial_ring.PolynomialRing_integral_domain):
    r"""
    A polynomial ring over a `p`-adic ring.

    TESTS::

        sage: R.<x> = Zp(3)[]
        sage: from sage.rings.polynomial.padics.polynomial_ring_padic import PolynomialRing_padic
        sage: isinstance(R, PolynomialRing_padic)
        True

    """
    pass

class PolynomialRing_dense_padic_ring_generic(PolynomialRing_padic):
    r"""
    A polynomial ring over a `p`-adic ring which is not a field.

    TESTS::

        sage: R.<x> = Zp(3)[]
        sage: from sage.rings.polynomial.padics.polynomial_ring_padic import PolynomialRing_dense_padic_ring_generic
        sage: isinstance(R, PolynomialRing_dense_padic_ring_generic)
        True

        sage: TestSuite(Zp(3)['x']).run()
        sage: TestSuite(ZpFM(3)['x']).run()
        sage: TestSuite(ZpCA(3)['x']).run()

    """
    pass

class PolynomialRing_dense_padic_field_generic(sage.rings.polynomial.polynomial_ring.PolynomialRing_field, PolynomialRing_padic):
    r"""
    A polynomial ring over a `p`-adic field.

    TESTS::

        sage: R.<x> = Qp(3)[]
        sage: from sage.rings.polynomial.padics.polynomial_ring_padic import PolynomialRing_dense_padic_field_generic
        sage: isinstance(R, PolynomialRing_dense_padic_field_generic)
        True

        sage: TestSuite(R).run()

    """
    pass
