"""
Finite fields implemented via PARI's FFELT type

AUTHORS:

- Peter Bruin (June 2013): initial version, based on
  finite_field_ext_pari.py by William Stein et al.
"""

#*****************************************************************************
#       Copyright (C) 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from .element_pari_ffelt import FiniteFieldElement_pari_ffelt
from .finite_field_base import FiniteField
from .finite_field_constructor import GF

class FiniteField_pari_ffelt(FiniteField):
    """
    Finite fields whose cardinality is a prime power (not a prime),
    implemented using PARI's ``FFELT`` type.

    INPUT:

    - ``p`` -- prime number

    - ``modulus`` -- an irreducible polynomial of degree at least 2
      over the field of `p` elements

    - ``name`` -- string: name of the distinguished generator
      (default: variable name of ``modulus``)

    OUTPUT:

    A finite field of order `q = p^n`, generated by a distinguished
    element with minimal polynomial ``modulus``.  Elements are
    represented as polynomials in ``name`` of degree less than `n`.

    .. NOTE::

        Direct construction of :class:`FiniteField_pari_ffelt` objects
        requires specifying a characteristic and a modulus.  To
        construct a finite field by specifying a cardinality and an
        algorithm for finding an irreducible polynomial, use the
        ``FiniteField`` constructor with ``impl='pari_ffelt'``.

    EXAMPLES:

    Some computations with a finite field of order 9::

        sage: k = FiniteField(9, 'a', impl='pari_ffelt')
        sage: k
        Finite Field in a of size 3^2
        sage: k.is_field()
        True
        sage: k.characteristic()
        3
        sage: a = k.gen()
        sage: a
        a
        sage: a.parent()
        Finite Field in a of size 3^2
        sage: a.charpoly('x')
        x^2 + 2*x + 2
        sage: [a^i for i in range(8)]
        [1, a, a + 1, 2*a + 1, 2, 2*a, 2*a + 2, a + 2]
        sage: TestSuite(k).run()

    Next we compute with a finite field of order 16::

        sage: k16 = FiniteField(16, 'b', impl='pari_ffelt')
        sage: z = k16.gen()
        sage: z
        b
        sage: z.charpoly('x')
        x^4 + x + 1
        sage: k16.is_field()
        True
        sage: k16.characteristic()
        2
        sage: z.multiplicative_order()
        15

    Illustration of dumping and loading::

        sage: K = FiniteField(7^10, 'b', impl='pari_ffelt')
        sage: loads(K.dumps()) == K
        True

        sage: K = FiniteField(10007^10, 'a', impl='pari_ffelt')
        sage: loads(K.dumps()) == K
        True
    """
    def __init__(self, q, modulus, name=None):
        """
        Create a finite field of ``q`` elements defined by the polynomial
        ``modulus``, with distinguished generator called ``name``.

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_pari_ffelt import FiniteField_pari_ffelt
            sage: R.<x> = PolynomialRing(GF(3))
            sage: k = FiniteField_pari_ffelt(3, x^2 + 2*x + 2, 'a'); k
            Finite Field in a of size 3^2
        """
        n = modulus.degree()
        if n < 2:
            raise ValueError("the degree must be at least 2")

        p, n = q.is_prime_power(get_data=True)
        FiniteField.__init__(self, base=GF(p), names=name, normalize=True)

        self._modulus = modulus
        self._degree = n

        self._gen_pari = modulus._pari_with_name(self._names[0]).ffgen()
        self._zero_element = self.element_class(self, 0)
        self._one_element = self.element_class(self, 1)
        self._gen = self.element_class(self, self._gen_pari)

    Element = FiniteFieldElement_pari_ffelt

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: k.<b> = FiniteField(5^20, impl='pari_ffelt')
            sage: type(k)
            <class 'sage.rings.finite_rings.finite_field_pari_ffelt.FiniteField_pari_ffelt_with_category'>
            sage: k is loads(dumps(k))
            True
        """
        return self._factory_data[0].reduce_data(self)

    def gen(self, n=0):
        """
        Return a generator of ``self`` over its prime field, which is a
        root of ``self.modulus()``.

        INPUT:

        - ``n`` -- must be 0

        OUTPUT:

        An element `a` of ``self`` such that ``self.modulus()(a) == 0``.

        .. WARNING::

            This generator is not guaranteed to be a generator for the
            multiplicative group.  To obtain the latter, use
            :meth:`~sage.rings.finite_rings.finite_field_base.FiniteFields.multiplicative_generator()`
            or use the ``modulus="primitive"`` option when constructing
            the field.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(GF(2))
            sage: FiniteField(2^4, 'b', impl='pari_ffelt').gen()
            b
            sage: k = FiniteField(3^4, 'alpha', impl='pari_ffelt')
            sage: a = k.gen()
            sage: a
            alpha
            sage: a^4
            alpha^3 + 1
        """
        if n:
            raise IndexError("only one generator")
        return self._gen

    def characteristic(self):
        """
        Return the characteristic of ``self``.

        EXAMPLES::

            sage: F = FiniteField(3^4, 'a', impl='pari_ffelt')
            sage: F.characteristic()
            3
        """
        # This works since self is not its own prime field.
        return self.base_ring().characteristic()

    def degree(self):
        """
        Returns the degree of ``self`` over its prime field.

        EXAMPLES::

            sage: F = FiniteField(3^20, 'a', impl='pari_ffelt')
            sage: F.degree()
            20
        """
        return self._degree
