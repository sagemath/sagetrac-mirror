"""
Finite fields implemented via FLINT's fq_t type.

AUTHORS:

- Jean-Pierre Flori (July 2014): initial version, based on
  finite_field_pari_ffelt.py by Peter Bruin.
"""

#*****************************************************************************
#       Copyright (C) 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/stdsage.pxi"
from sage.libs.flint.types cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.fmpz_mod_poly cimport *
from sage.libs.flint.fq cimport *

import sage.rings.integer
from sage.rings.integer cimport Integer
from sage.rings.finite_rings.element_flint_fq cimport FiniteFieldElement_flint_fq

cdef class FiniteField_flint_fq(FiniteField):
    """
    Finite fields whose cardinality is a prime power (not a prime),
    implemented using FLINT's ``fq_t`` type.

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

        - Direct construction of FiniteField_flint_fq objects
          requires specifying a characteristic and a modulus.
          To construct a finite field by specifying a cardinality and
          an algorithm for finding an irreducible polynomial, use the
          FiniteField constructor with ``impl='flint_fq'``.

        - Two finite fields are considered equal if and only if they
          have the same cardinality, variable name, and modulus.

    EXAMPLES:

    Some computations with a finite field of order 9::

        sage: k = FiniteField(9, 'a', impl='flint_fq')
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

        sage: k16 = FiniteField(16, 'b', impl='flint_fq')
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

        sage: K = FiniteField(7^10, 'b', impl='flint_fq')
        sage: loads(K.dumps()) == K
        True

        sage: K = FiniteField(10007^10, 'a', impl='flint_fq')
        sage: loads(K.dumps()) == K
        True
    """
    def __init__(self, p, modulus, name=None):
        """
        Create a finite field of characteristic `p` defined by the
        polynomial ``modulus``, with distinguished generator called
        ``name``.

        EXAMPLE::

            sage: from sage.rings.finite_rings.finite_field_flint_fq import FiniteField_flint_fq
            sage: R.<x> = PolynomialRing(GF(3))
            sage: k = FiniteField_flint_fq(3, x^2 + 2*x + 2, 'a'); k
            Finite Field in a of size 3^2
        """
        from sage.rings.finite_rings.constructor import GF
        #from sage.rings.integer import Integer
        from sage.structure.proof.all import arithmetic
        proof = arithmetic()

        p = Integer(p)
        if ((p < 2)
            or (proof and not p.is_prime())
            or (not proof and not p.is_pseudoprime())):
            raise ArithmeticError("p must be a prime number")
        Fp = GF(p)

        if name is None:
            name = modulus.variable_name()

        FiniteField.__init__(self, base=Fp, names=name, normalize=True)

        modulus = self.polynomial_ring()(modulus)
        n = modulus.degree()
        if n < 2:
            raise ValueError("the degree must be at least 2")

        self._modulus = modulus
        self._degree = n
        self._kwargs = {}

        self.__hash = hash((p**n, name, modulus))

        # Cannot be done in cinit as we need modulus
        cdef Integer ci
        cdef fmpz_t cflint, pflint
        cdef fmpz_mod_poly_t modulus_flint
        fmpz_init_set_readonly(pflint, (<Integer>p).value)
        fmpz_mod_poly_init(modulus_flint, pflint)
        fmpz_clear_readonly(pflint)
        for i in xrange(n+1):
            ci = Integer(modulus[i])
            fmpz_init_set_readonly(cflint, ci.value)
            fmpz_mod_poly_set_coeff_fmpz(modulus_flint, i, cflint)
            fmpz_clear_readonly(cflint)
        fq_ctx_init_modulus(self._ctx, modulus_flint, <char *> (name[0]))
        self._ctx_initialized = 1
        fmpz_mod_poly_clear(modulus_flint)

        self._zero_element = self.element_class(self, 0)
        self._one_element = self.element_class(self, 1)
        cdef fq_t gen_flint
        cdef FiniteFieldElement_flint_fq gen
        fq_init(gen_flint, self._ctx)
        fq_gen(gen_flint, self._ctx)
        gen = (<FiniteFieldElement_flint_fq>self._zero_element)._new()
        gen.set_from_fq(gen_flint)
        fq_clear(gen_flint, self._ctx)
        self._gen = gen

    def __cinit__(FiniteField_flint_fq self):
        self._ctx = <fq_ctx_struct *>sage_malloc(sizeof(fq_ctx_t))

    def __dealloc__(FiniteField_flint_fq self):
        if self._ctx_initialized:
            fq_ctx_clear(self._ctx)
        if self._ctx:
            sage_free(self._ctx)

    Element = FiniteFieldElement_flint_fq

    def __hash__(self):
        """
        Return the hash of this field.

        EXAMPLE::

            sage: {GF(9, 'a', impl='flint_fq'): 1} # indirect doctest
            {Finite Field in a of size 3^2: 1}
            sage: {GF(9, 'b', impl='flint_fq'): 1} # indirect doctest
            {Finite Field in b of size 3^2: 1}
        """
        return self.__hash

    def __reduce__(self):
        """
        For pickling.

        EXAMPLE::

            sage: k.<b> = FiniteField(5^20, impl='flint_fq')
            sage: type(k)
            <type 'sage.rings.finite_rings.finite_field_flint_fq.FiniteField_flint_fq'>
            sage: k is loads(dumps(k))
            True
        """
        return self._factory_data[0].reduce_data(self)

    def _cmp_(self, other):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: k = FiniteField(7^20, 'a', impl='flint_fq')
            sage: k == k
            True
            sage: k2 = FiniteField(7^20, 'a', impl='flint_fq')
            sage: k2 == k
            True
            sage: kb = FiniteField(7^20, 'b', impl='flint_fq')
            sage: kb == k
            False
        """
        if not isinstance(other, FiniteField_flint_fq):
            return cmp(type(self), type(other))
        return cmp((self.cardinality(), self.variable_name(), self._modulus),
                   (other.cardinality(), other.variable_name(), other._modulus))

    def __richcmp__(left, right, op):
        """
        Compare ``left`` with ``right``.

        EXAMPLE::

            sage: k = FiniteField(2^17, 'a', impl='flint_fq')
            sage: j = FiniteField(2^18, 'b', impl='flint_fq')
            sage: k == j
            False

            sage: k == copy(k)
            True
        """
        return left._richcmp_helper(right, op)

    def gen(self, n=0):
        """
        Return a generator of the finite field.

        INPUT:

        - ``n`` -- ignored

        OUTPUT:

        A generator of the finite field.

        This generator is a root of the defining polynomial of the
        finite field.

        .. WARNING::

            This generator is not guaranteed to be a generator
            for the multiplicative group.  To obtain the latter, use
            :meth:`~sage.rings.finite_rings.finite_field_base.FiniteFields.multiplicative_generator()`.

        EXAMPLE::

            sage: R.<x> = PolynomialRing(GF(2))
            sage: FiniteField(2^4, 'b', impl='flint_fq').gen()
            b
            sage: k = FiniteField(3^4, 'alpha', impl='flint_fq')
            sage: a = k.gen()
            sage: a
            alpha
            sage: a^4
            alpha^3 + 1
        """
        return self._gen

    def characteristic(self):
        """
        Return the characteristic of ``self``.

        EXAMPLE::

            sage: F = FiniteField(3^4, 'a', impl='flint_fq')
            sage: F.characteristic()
            3
        """
        # This works since self is not its own prime field.
        return self.base_ring().characteristic()

    def degree(self):
        """
        Returns the degree of ``self`` over its prime field.

        EXAMPLE::

            sage: F = FiniteField(3^20, 'a', impl='flint_fq')
            sage: F.degree()
            20
        """
        return self._degree

    def polynomial(self):
        """
        Return the minimal polynomial of the generator of ``self`` in
        ``self.polynomial_ring()``.

        EXAMPLES::

            sage: F = FiniteField(3^2, 'a', impl='flint_fq')
            sage: F.polynomial()
            a^2 + 2*a + 2

            sage: F = FiniteField(7^20, 'a', impl='flint_fq')
            sage: f = F.polynomial(); f
            a^20 + a^12 + 6*a^11 + 2*a^10 + 5*a^9 + 2*a^8 + 3*a^7 + a^6 + 3*a^5 + 3*a^3 + a + 3
            sage: f(F.gen())
            0
        """
        return self._modulus

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

            sage: k = FiniteField(3^4, 'a', impl='flint_fq')
            sage: b = k(5) # indirect doctest
            sage: b.parent()
            Finite Field in a of size 3^4
            sage: a = k.gen()
            sage: k(a + 2)
            a + 2

        Univariate polynomials coerce into finite fields by evaluating
        the polynomial at the field's generator::

            sage: R.<x> = QQ[]
            sage: k, a = FiniteField(5^2, 'a', impl='flint_fq').objgen()
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

            sage: Q, q = FiniteField(5^7, 'q', impl='flint_fq').objgen()
            sage: L = GF(5)
            sage: LL.<xx> = L[]
            sage: Q(xx^2 + 2*xx + 4)
            q^2 + 2*q + 4

            sage: k = FiniteField(3^11, 't', impl='flint_fq')
            sage: k.polynomial()
            t^11 + 2*t^2 + 1
            sage: P = k.polynomial_ring()
            sage: k(P.0^11)
            t^2 + 2

        An element can be specified by its vector of coordinates with
        respect to the basis consisting of powers of the generator:

            sage: k = FiniteField(3^11, 't', impl='flint_fq')
            sage: V = k.vector_space()
            sage: V
            Vector space of dimension 11 over Finite Field of size 3
            sage: v = V([0,1,2,0,1,2,0,1,2,0,1])
            sage: k(v)
            t^10 + 2*t^8 + t^7 + 2*t^5 + t^4 + 2*t^2 + t

        Multivariate polynomials only coerce if constant::

            sage: k = FiniteField(5^2, 'a', impl='flint_fq')
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

            sage: F = FiniteField(2^3, 'a', impl='flint_fq')
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

            sage: F = FiniteField(2^4, 'a', impl='flint_fq')
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

            sage: k = FiniteField(3^2, 'a', impl='flint_fq')
            sage: a = k(11); a
            2
            sage: a.parent()
            Finite Field in a of size 3^2
            sage: V = k.vector_space(); v = V((1,2))
            sage: k(v)
            2*a + 1

        We create elements using a list and verify that :trac:`10486` has
        been fixed::

            sage: k = FiniteField(3^11, 't', impl='flint_fq')
            sage: x = k([1,0,2,1]); x
            t^3 + 2*t^2 + 1
            sage: x + x + x
            0
            sage: pari(x)
            t^3 + 2*t^2 + 1

        If the list is longer than the degree, we just get the result
        modulo the modulus::

            sage: from sage.rings.finite_rings.finite_field_flint_fq import FiniteField_flint_fq
            sage: R.<a> = PolynomialRing(GF(5))
            sage: k = FiniteField_flint_fq(5, a^2 - 2, 't')
            sage: x = k([0,0,0,1]); x
            2*t
            sage: pari(x)
            2*t

        When initializing from a list, the elements are first coerced
        to the prime field (:trac:`11685`)::

            sage: k = FiniteField(3^11, 't', impl='flint_fq')
            sage: k([ 0, 1/2 ])
            2*t
            sage: k([ k(0), k(1) ])
            t
            sage: k([ GF(3)(2), GF(3^5,'u')(1) ])
            t + 2
            sage: R.<x> = PolynomialRing(k)
            sage: k([ R(-1), x/x ])
            t + 2
        """
        if isinstance(x, self.element_class) and x.parent() is self:
            return x
        else:
            return self.element_class(self, x)
