"""
Finite fields implemented via PARI's FFELT type

AUTHORS:

- Peter Bruin (June 2013): initial version, based on
  finite_field_ext_pari.py by William Stein et al.
"""
from __future__ import absolute_import

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
    def __init__(self, p, modulus, name=None):
        """
        Create a finite field of characteristic `p` defined by the
        polynomial ``modulus``, with distinguished generator called
        ``name``.

        EXAMPLE::

            sage: from sage.rings.finite_rings.finite_field_pari_ffelt import FiniteField_pari_ffelt
            sage: R.<x> = PolynomialRing(GF(3))
            sage: k = FiniteField_pari_ffelt(3, x^2 + 2*x + 2, 'a'); k
            Finite Field in a of size 3^2
        """
        n = modulus.degree()
        if n < 2:
            raise ValueError("the degree must be at least 2")

        FiniteField.__init__(self, base=GF(p), names=name, normalize=True)

        self._modulus = modulus
        self._degree = n
        self._kwargs = {}

        self._gen_pari = modulus._pari_with_name(self._names[0]).ffgen()
        self._zero_element = self.element_class(self, 0)
        self._one_element = self.element_class(self, 1)
        self._gen = self.element_class(self, self._gen_pari)

    Element = FiniteFieldElement_pari_ffelt

    def __reduce__(self):
        """
        For pickling.

        EXAMPLE::

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

        EXAMPLE::

            sage: F = FiniteField(3^4, 'a', impl='pari_ffelt')
            sage: F.characteristic()
            3
        """
        # This works since self is not its own prime field.
        return self.base_ring().characteristic()

    def degree(self):
        """
        Returns the degree of ``self`` over its prime field.

        EXAMPLE::

            sage: F = FiniteField(3^20, 'a', impl='pari_ffelt')
            sage: F.degree()
            20
        """
        return self._degree

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        INPUT:

        - ``x`` -- object

        OUTPUT:

        A finite field element generated from `x`, if possible.

        .. NOTE::

            If `x` is a list or an element of the underlying vector
            space of the finite field, then it is interpreted as the
            list of coefficients of a polynomial over the prime field,
            and that polynomial is interpreted as an element of the
            finite field.

        EXAMPLES::

            sage: k = FiniteField(3^4, 'a', impl='pari_ffelt')
            sage: b = k(5) # indirect doctest
            sage: b.parent()
            Finite Field in a of size 3^4
            sage: a = k.gen()
            sage: k(a + 2)
            a + 2

        Univariate polynomials coerce into finite fields by evaluating
        the polynomial at the field's generator::

            sage: R.<x> = QQ[]
            sage: k.<a> = FiniteField(5^2, 'a', impl='pari_ffelt')
            sage: k(R(2/3))
            4
            sage: k(x^2)
            a + 3

            sage: R.<x> = GF(5)[]
            sage: k(x^3-2*x+1)
            2*a + 4

            sage: x = polygen(QQ)
            sage: k(x^25)
            a

            sage: Q.<q> = FiniteField(5^7, 'q', impl='pari_ffelt')
            sage: L = GF(5)
            sage: LL.<xx> = L[]
            sage: Q(xx^2 + 2*xx + 4)
            q^2 + 2*q + 4

            sage: k = FiniteField(3^11, 't', impl='pari_ffelt')
            sage: k.polynomial()
            t^11 + 2*t^2 + 1
            sage: P = k.polynomial_ring()
            sage: k(P.0^11)
            t^2 + 2

        An element can be specified by its vector of coordinates with
        respect to the basis consisting of powers of the generator:

            sage: k = FiniteField(3^11, 't', impl='pari_ffelt')
            sage: V = k.vector_space()
            sage: V
            Vector space of dimension 11 over Finite Field of size 3
            sage: v = V([0,1,2,0,1,2,0,1,2,0,1])
            sage: k(v)
            t^10 + 2*t^8 + t^7 + 2*t^5 + t^4 + 2*t^2 + t

        Multivariate polynomials only coerce if constant::

            sage: k = FiniteField(5^2, 'a', impl='pari_ffelt')
            sage: R = k['x,y,z']; R
            Multivariate Polynomial Ring in x, y, z over Finite Field in a of size 5^2
            sage: k(R(2))
            2
            sage: R = QQ['x,y,z']
            sage: k(R(1/5))
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

        Gap elements can also be coerced into finite fields::

            sage: F = FiniteField(2^3, 'a', impl='pari_ffelt')
            sage: a = F.multiplicative_generator(); a
            a
            sage: b = gap(a^3); b
            Z(2^3)^3
            sage: F(b)
            a + 1
            sage: a^3
            a + 1

            sage: a = GF(13)(gap('0*Z(13)')); a
            0
            sage: a.parent()
            Finite Field of size 13

            sage: F = FiniteField(2^4, 'a', impl='pari_ffelt')
            sage: F(gap('Z(16)^3'))
            a^3
            sage: F(gap('Z(16)^2'))
            a^2

        You can also call a finite extension field with a string
        to produce an element of that field, like this::

            sage: k = GF(2^8, 'a')
            sage: k('a^200')
            a^4 + a^3 + a^2

        This is especially useful for conversion from Singular etc.

        TESTS::

            sage: k = FiniteField(3^2, 'a', impl='pari_ffelt')
            sage: a = k(11); a
            2
            sage: a.parent()
            Finite Field in a of size 3^2
            sage: V = k.vector_space(); v = V((1,2))
            sage: k(v)
            2*a + 1

        We create elements using a list and verify that :trac:`10486` has
        been fixed::

            sage: k = FiniteField(3^11, 't', impl='pari_ffelt')
            sage: x = k([1,0,2,1]); x
            t^3 + 2*t^2 + 1
            sage: x + x + x
            0
            sage: pari(x)
            t^3 + 2*t^2 + 1

        If the list is longer than the degree, we just get the result
        modulo the modulus::

            sage: from sage.rings.finite_rings.finite_field_pari_ffelt import FiniteField_pari_ffelt
            sage: R.<a> = PolynomialRing(GF(5))
            sage: k = FiniteField_pari_ffelt(5, a^2 - 2, 't')
            sage: x = k([0,0,0,1]); x
            2*t
            sage: pari(x)
            2*t

        When initializing from a list, the elements are first coerced
        to the prime field (:trac:`11685`)::

            sage: k = FiniteField(3^11, 't', impl='pari_ffelt')
            sage: k([ 0, 1/2 ])
            2*t
            sage: k([ k(0), k(1) ])
            t
            sage: k([ GF(3)(2), GF(3^5,'u')(1) ])
            t + 2
            sage: R.<x> = PolynomialRing(k)
            sage: k([ R(-1), x/x ])
            t + 2

        Check that zeros are created correctly (:trac:`11685`)::

            sage: K = FiniteField(3^11, 't', impl='pari_ffelt'); a = K.0
            sage: v = 0; pari(K(v))
            0
            sage: v = Mod(0,3); pari(K(v))
            0
            sage: v = pari(0); pari(K(v))
            0
            sage: v = pari("Mod(0,3)"); pari(K(v))
            0
            sage: v = []; pari(K(v))
            0
            sage: v = [0]; pari(K(v))
            0
            sage: v = [0,0]; pari(K(v))
            0
            sage: v = pari("Pol(0)"); pari(K(v))
            0
            sage: v = pari("Mod(0, %s)"%K.modulus()); pari(K(v))
            0
            sage: v = pari("Mod(Pol(0), %s)"%K.modulus()); pari(K(v))
            0
            sage: v = K(1) - K(1); pari(K(v))
            0
            sage: v = K([1]) - K([1]); pari(K(v))
            0
            sage: v = a - a; pari(K(v))
            0
            sage: v = K(1)*0; pari(K(v))
            0
            sage: v = K([1])*K([0]); pari(K(v))
            0
            sage: v = a*0; pari(K(v))
            0
        """
        if isinstance(x, self.element_class) and x.parent() is self:
            return x
        else:
            return self.element_class(self, x)
