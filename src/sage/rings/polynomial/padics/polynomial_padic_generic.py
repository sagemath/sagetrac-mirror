r"""
Common implementation for polynomials over `p`-adic rings

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
from sage.rings.polynomial.polynomial_element import Polynomial_generic_dense
from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_dense_field, Polynomial_generic_domain
from polynomial_padic import Polynomial_padic

class Polynomial_padic_generic_ring(Polynomial_padic, Polynomial_generic_domain, Polynomial_generic_dense):
    r"""
    A polynomial over a `p`-adic ring which is not a field.

    INPUT:

    - ``parent`` -- a polynomial ring over a `p`-adic ring

    - ``x`` -- data which can be used to construct a polynomial

    - ``check`` -- a boolean (default: ``True``), unused

    - ``is_gen`` -- whether this is the generator of the polynomial ring
      (default: ``False``)

    - ``construct`` -- a boolean (default: ``False``), unused

    - ``absprec`` -- an integer or ``None`` (default: ``None``), unused

    .. NOTE::

        In contrast to :class:`polynomial_padic.Polynomial_padic`, this class
        is meant as a base class for implementations which want to use the
        generic handling of the polynomial data from
        :class:`sage.rings.polynomial.polynomial_element.Polynomial_generic_dense`.

    EXAMPLES::

        sage: R.<x> = ZpFM(3)[] # indirect doctest
        sage: from sage.rings.polynomial.padics.polynomial_padic_generic import Polynomial_padic_generic_ring
        sage: isinstance(x, Polynomial_padic_generic_ring)
        True

    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False, absprec=None):
        r"""
        Initialization.

        TESTS::

            sage: from sage.rings.polynomial.padics.polynomial_padic_generic import Polynomial_padic_generic_ring
            sage: R.<x> = ZpFM(3)[]
            sage: type(Polynomial_padic_generic_ring(R, None))
            <class 'sage.rings.polynomial.padics.polynomial_padic_generic.Polynomial_padic_generic_ring'>

        """
        Polynomial_padic.__init__(self, parent, is_gen=is_gen)
        Polynomial_generic_domain.__init__(self, parent, is_gen=is_gen)
        Polynomial_generic_dense.__init__(self, parent, x)

class Polynomial_padic_generic_field(Polynomial_padic, Polynomial_generic_dense_field):
    r"""
    A polynomial over a `p`-adic field.

    INPUT:

    - ``parent`` -- a polynomial ring over a `p`-adic ring

    - ``x`` -- data which can be used to construct a polynomial

    - ``check`` -- a boolean (default: ``True``), unused

    - ``is_gen`` -- whether this is the generator of the polynomial ring
      (default: ``False``)

    - ``construct`` -- a boolean (default: ``False``), unused

    - ``absprec`` -- an integer or ``None`` (default: ``None``), unused

    .. NOTE::

        In contrast to :class:`polynomial_padic.Polynomial_padic`, this class
        is meant as a base class for implementations which want to use the
        generic handling of the polynomial data from
        :class:`sage.rings.polynomial.polynomial_element.Polynomial_generic_dense`.

    EXAMPLES::

        sage: R.<a> = Qp(3)[]
        sage: L.<a> = Qp(3).extension(a^2 - 3)
        sage: R.<x> = L[] # indirect doctest
        sage: from sage.rings.polynomial.padics.polynomial_padic_generic import Polynomial_padic_generic_field
        sage: isinstance(x, Polynomial_padic_generic_field)
        True

    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False, absprec=None):
        r"""
        Initialization.

        TESTS::

            sage: from sage.rings.polynomial.padics.polynomial_padic_generic import Polynomial_padic_generic_field
            sage: R.<a> = Qp(3)[]
            sage: L.<a> = Qp(3).extension(a^2 - 3)
            sage: R.<x> = L[]
            sage: type(Polynomial_padic_generic_field(R, None))
            <class 'sage.rings.polynomial.padics.polynomial_padic_generic.Polynomial_padic_generic_field'>

        """
        Polynomial_padic.__init__(self, parent, is_gen=is_gen)
        Polynomial_generic_dense_field.__init__(self, parent, x)
