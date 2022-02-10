# -*- coding: utf-8 -*-
r"""
Fields
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#                2012-2014 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.category_singleton import Category_contains_method_by_parent_class
from sage.categories.euclidean_domains import EuclideanDomains
from sage.categories.division_rings import DivisionRings

from sage.structure.element import coerce_binop

class Fields(CategoryWithAxiom):
    """
    The category of (commutative) fields, i.e. commutative rings where
    all non-zero elements have multiplicative inverses

    EXAMPLES::

        sage: K = Fields()
        sage: K
        Category of fields
        sage: Fields().super_categories()
        [Category of euclidean domains, Category of division rings]

        sage: K(IntegerRing())
        Rational Field
        sage: K(PolynomialRing(GF(3), 'x'))
        Fraction Field of Univariate Polynomial Ring in x over
        Finite Field of size 3
        sage: K(RealField())
        Real Field with 53 bits of precision

    TESTS::

        sage: TestSuite(Fields()).run()
    """
    _base_category_class_and_axiom = (DivisionRings, "Commutative")

    def extra_super_categories(self):
        """
        EXAMPLES::

            sage: Fields().extra_super_categories()
            [Category of euclidean domains]

        """
        return [EuclideanDomains()]

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: GF(4, "a") in Fields()
            True
            sage: QQ in Fields()
            True
            sage: ZZ in Fields()
            False
            sage: IntegerModRing(4) in Fields()
            False
            sage: InfinityRing in Fields()
            False

        This implementation will not be needed anymore once every
        field in Sage will be properly declared in the category
        :class:`Fields`().

        Caveat: this should eventually be fixed::

            sage: gap.Rationals in Fields()
            False

        typically by implementing the method :meth:`category`
        appropriately for Gap objects::

            sage: GR = gap.Rationals
            sage: GR.category = lambda : Fields()
            sage: GR in Fields()
            True

        The following tests against a memory leak fixed in :trac:`13370`. In order
        to prevent non-deterministic deallocation of fields that have been created
        in other doctests, we introduced a strong reference to all previously created
        uncollected objects in :trac:`19244`. ::

            sage: import gc
            sage: _ = gc.collect()
            sage: permstore = [X for X in gc.get_objects() if isinstance(X, sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic)]
            sage: n = len(permstore)
            sage: for i in prime_range(100):
            ....:     R = ZZ.quotient(i)
            ....:     t = R in Fields()

        First, we show that there are now more quotient rings in cache than before::

            sage: len([X for X in gc.get_objects() if isinstance(X, sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic)]) > n
            True

        When we delete the last quotient ring created in the loop and then do a garbage
        collection, all newly created rings vanish::

            sage: del R
            sage: _ = gc.collect()
            sage: len([X for X in gc.get_objects() if isinstance(X, sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic)]) - n
            0

        """
        import sage.rings.ring
        try:
            return self._contains_helper(x) or sage.rings.ring._is_Field(x)
        except Exception:
            return False

    @lazy_class_attribute
    def _contains_helper(cls):
        """
        Helper for containment tests in the category of fields.

        This helper just tests whether the given object's category
        is already known to be a sub-category of the category of
        fields. There are, however, rings that are initialised
        as plain commutative rings and found out to be fields
        only afterwards. Hence, this helper alone is not enough
        for a proper containment test.

        TESTS::

            sage: P.<x> = QQ[]
            sage: Q = P.quotient(x^2+2)
            sage: Q.category()
            Category of commutative no zero divisors quotients of algebras
             over (number fields and quotient fields and metric spaces)
            sage: F = Fields()
            sage: F._contains_helper(Q)
            False
            sage: Q in F  # This changes the category!
            True
            sage: F._contains_helper(Q)
            True

        """
        return Category_contains_method_by_parent_class(cls())

    def _call_(self, x):
        """
        Construct a field from the data in ``x``

        EXAMPLES::

            sage: K = Fields()
            sage: K
            Category of fields
            sage: Fields().super_categories()
            [Category of euclidean domains, Category of division rings]

            sage: K(IntegerRing()) # indirect doctest
            Rational Field
            sage: K(PolynomialRing(GF(3), 'x')) # indirect doctest
            Fraction Field of Univariate Polynomial Ring in x over
            Finite Field of size 3
            sage: K(RealField())
            Real Field with 53 bits of precision
        """
        try:
            return x.fraction_field()
        except AttributeError:
            raise TypeError("unable to associate a field to %s"%x)

    Finite = LazyImport('sage.categories.finite_fields', 'FiniteFields', at_startup=True)

    class ParentMethods:

        def is_field( self, proof=True ):
            r"""
            Returns True as ``self`` is a field.

            EXAMPLES::

                sage: QQ.is_field()
                True
                sage: Parent(QQ,category=Fields()).is_field()
                True
            """
            return True

        def is_integrally_closed(self):
            r"""
            Return ``True``, as per :meth:`IntegralDomain.is_integrally_closed`:
            for every field `F`, `F` is its own field of fractions,
            hence every element of `F` is integral over `F`.

            EXAMPLES::

                sage: QQ.is_integrally_closed()
                True
                sage: QQbar.is_integrally_closed()
                True
                sage: Z5 = GF(5); Z5
                Finite Field of size 5
                sage: Z5.is_integrally_closed()
                True
            """
            return True

        def _gcd_univariate_polynomial(self, a, b):
            """
            Return the gcd of ``a`` and ``b`` as a monic polynomial.

            INPUT:

            - ``a``, ``b`` -- two univariate polynomials defined over ``self``

            .. WARNING::

                If the base ring is inexact, the results may not be
                entirely stable.

            EXAMPLES::

                sage: R.<x> = QQbar[]
                sage: QQbar._gcd_univariate_polynomial(2*x, 2*x^2)
                x

            TESTS::

                sage: for A in (RR, CC, QQbar):
                ....:     g = A._gcd_univariate_polynomial
                ....:     R.<x> = A[]
                ....:     z = R.zero()
                ....:     assert(g(2*x, 2*x^2) == x and
                ....:            g(z, 2*x) == x and
                ....:            g(2*x, z) == x and
                ....:            g(z, z) == z)

                sage: R.<x> = RR[]
                sage: (x^3).gcd(x^5+1)
                1.00000000000000
                sage: (x^3).gcd(x^5+x^2)
                x^2
                sage: f = (x+3)^2 * (x-1)
                sage: g = (x+3)^5
                sage: f.gcd(g)
                x^2 + 6.00000000000000*x + 9.00000000000000

            The following example illustrates the fact that for inexact
            base rings, the returned gcd is often 1 due to rounding::

                sage: f = (x+RR.pi())^2 * (x-1)
                sage: g = (x+RR.pi())^5
                sage: f.gcd(g)
                1.00000000000000

            Check :trac:`23012`::

                sage: R.<x> = QQ[]
                sage: Q = R.quotient(x^2-1)   # Not a field
                sage: P.<x> = Q[]
                sage: def always_True(*args, **kwds): return True
                sage: Q.is_field = always_True
                sage: Q in Fields()
                True
                sage: Q._gcd_univariate_polynomial(x, x)
                x
            """
            while b:
                q, r = a.quo_rem(b)
                a, b = b, r
            if a:
                a = a.monic()
            return a

        def _xgcd_univariate_polynomial(self, a, b):
            r"""
            Return an extended gcd of ``a`` and ``b``.

            INPUT:

            - ``a``, ``b`` -- two univariate polynomials defined over ``self``

            OUTPUT:

            A tuple ``(d, u, v)`` of polynomials such that ``d`` is the
            greatest common divisor (monic or zero) of ``a`` and ``b``,
            and ``u``, ``v`` satisfy ``d = u*a + v*b``.

            .. WARNING::

                If the base ring is inexact, the results may not be
                entirely stable.

            ALGORITHM:

            This uses the extended Euclidean algorithm; see for example
            [Coh1993]_, Algorithm 3.2.2.

            EXAMPLES::

                sage: P.<x> = QQ[]
                sage: F = (x^2 + 2)*x^3; G = (x^2+2)*(x-3)
                sage: g, u, v = QQ._xgcd_univariate_polynomial(F,G)
                sage: g, u, v
                (x^2 + 2, 1/27, -1/27*x^2 - 1/9*x - 1/3)
                sage: u*F + v*G
                x^2 + 2

            ::

                sage: g, u, v = QQ._xgcd_univariate_polynomial(x,P(0)); g, u, v
                (x, 1, 0)
                sage: g == u*x + v*P(0)
                True
                sage: g, u, v = QQ._xgcd_univariate_polynomial(P(0),x); g, u, v
                (x, 0, 1)
                sage: g == u*P(0) + v*x
                True

            TESTS::

                sage: for A in (RR, CC, QQbar):
                ....:     g = A._xgcd_univariate_polynomial
                ....:     R.<x> = A[]
                ....:     z, h = R(0), R(1/2)
                ....:     assert(g(2*x, 2*x^2) == (x, h, z) and
                ....:            g(z, 2*x) == (x, z, h) and
                ....:            g(2*x, z) == (x, h, z) and
                ....:            g(z, z) == (z, z, z))

                sage: P.<x> = QQ[]
                sage: F = (x^2 + 2)*x^3; G = (x^2+2)*(x-3)
                sage: g, u, v = QQ._xgcd_univariate_polynomial(F,G)
                sage: g, u, v
                (x^2 + 2, 1/27, -1/27*x^2 - 1/9*x - 1/3)
                sage: u*F + v*G
                x^2 + 2

            We check that the behavior of xgcd with zero elements is
            compatible with gcd (:trac:`17671`)::

                sage: R.<x> = QQbar[]
                sage: zero = R.zero()
                sage: zero.xgcd(2*x)
                (x, 0, 1/2)
                sage: (2*x).xgcd(zero)
                (x, 1/2, 0)
                sage: zero.xgcd(zero)
                (0, 0, 0)
            """
            R = a.parent()
            zero = R.zero()
            if not b:
                if not a:
                    return (zero, zero, zero)
                c = ~a.leading_coefficient()
                return (c*a, R(c), zero)
            elif not a:
                c = ~b.leading_coefficient()
                return (c*b, zero, R(c))
            (u, d, v1, v3) = (R.one(), a, zero, b)
            while v3:
                q, r = d.quo_rem(v3)
                (u, d, v1, v3) = (v1, v3, u - v1*q, r)
            v = (d - a*u) // b
            if d:
                c = ~d.leading_coefficient()
                d, u, v = c*d, c*u, c*v
            return d, u, v

        def is_perfect(self):
            r"""
            Return whether this field is perfect, i.e., its characteristic is
            `p=0` or every element has a `p`-th root.

            EXAMPLES::

                sage: QQ.is_perfect()
                True
                sage: GF(2).is_perfect()
                True
                sage: FunctionField(GF(2), 'x').is_perfect()
                False

            """
            if self.characteristic() == 0:
                return True
            else:
                raise NotImplementedError

        def _test_characteristic_fields(self, **options):
            """
            Run generic tests on the method :meth:`.characteristic`.

            EXAMPLES::

                sage: QQ._test_characteristic_fields()

            .. NOTE::

                We cannot call this method ``_test_characteristic`` since that
                would overwrite the method in the super category, and for
                cython classes just calling
                ``super(sage.categories.fields.Fields().parent_class,
                self)._test_characteristic`` doesn't have the desired effect.

            .. SEEALSO::

                :meth:`sage.categories.rings.Rings.ParentMethods._test_characteristic`
            """
            tester = self._tester(**options)
            try:
                char = self.characteristic()
                tester.assertTrue(char.is_zero() or char.is_prime())
            except AttributeError:
                return
                # raised when self.one() does not have a additive_order() [or when char is an int and not an Integer which is already checked by _test_characteristic for rings]
            except NotImplementedError:
                return

        def fraction_field(self):
            r"""
            Returns the *fraction field* of ``self``, which is ``self``.

            EXAMPLES::

                sage: QQ.fraction_field() is QQ
                True
            """
            return self

        def _squarefree_decomposition_univariate_polynomial(self, f):
            r"""
            Return the square-free decomposition of ``f`` over this field.

            This is a helper method for
            :meth:`sage.rings.polynomial.squarefree_decomposition`.

            INPUT:

            - ``f`` -- a univariate non-zero polynomial over this field

            ALGORITHM: For rings of characteristic zero, we use the algorithm
            described in [Yun1976]_. Other fields may provide their own
            implementation by overriding this method.

            EXAMPLES::

                sage: x = polygen(QQ)
                sage: p = 37 * (x-1)^3 * (x-2)^3 * (x-1/3)^7 * (x-3/7)
                sage: p.squarefree_decomposition()
                (37*x - 111/7) * (x^2 - 3*x + 2)^3 * (x - 1/3)^7
                sage: p = 37 * (x-2/3)^2
                sage: p.squarefree_decomposition()
                (37) * (x - 2/3)^2
                sage: x = polygen(GF(3))
                sage: x.squarefree_decomposition()
                x
                sage: f = QQbar['x'](1)
                sage: f.squarefree_decomposition()
                1
            """
            from sage.structure.factorization import Factorization
            if f.degree() == 0:
                return Factorization([], unit=f[0])
            if self.characteristic() != 0:
                raise NotImplementedError("square-free decomposition not implemented for this polynomial.")

            factors = []
            cur = f
            f = [f]
            while cur.degree() > 0:
                cur = cur.gcd(cur.derivative())
                f.append(cur)

            g = []
            for i in range(len(f) - 1):
                g.append(f[i] // f[i+1])

            a = []
            for i in range(len(g) - 1):
                a.append(g[i] // g[i+1])
            a.append(g[-1])

            unit = f[-1]
            for i in range(len(a)):
                if a[i].degree() > 0:
                    factors.append((a[i], i+1))
                else:
                    unit = unit * a[i].constant_coefficient() ** (i + 1)

            return Factorization(factors, unit=unit, sort=False)

        def vector_space(self, *args, **kwds):
            r"""
            Gives an isomorphism of this field with a vector space over a subfield.

            This method is an alias for ``free_module``, which may have more documentation.

            INPUT:

            - ``base`` -- a subfield or morphism into this field (defaults to the base field)

            - ``basis`` -- a basis of the field as a vector space
              over the subfield; if not given, one is chosen automatically

            - ``map`` -- whether to return maps from and to the vector space

            OUTPUT:

            - ``V`` -- a vector space over ``base``
            - ``from_V`` -- an isomorphism from ``V`` to this field
            - ``to_V`` -- the inverse isomorphism from this field to ``V``

            EXAMPLES::

                sage: K.<a> = Qq(125)                                   # optional - sage.rings.padics
                sage: V, fr, to = K.vector_space()                      # optional - sage.rings.padics
                sage: v = V([1, 2, 3])                                  # optional - sage.rings.padics
                sage: fr(v, 7)                                          # optional - sage.rings.padics
                (3*a^2 + 2*a + 1) + O(5^7)
            """
            return self.free_module(*args, **kwds)

    class ElementMethods:
        def euclidean_degree(self):
            r"""
            Return the degree of this element as an element of an Euclidean
            domain.

            In a field, this returns 0 for all but the zero element (for
            which it is undefined).

            EXAMPLES::

                sage: QQ.one().euclidean_degree()
                0
            """
            if self.is_zero():
                raise ValueError("euclidean degree not defined for the zero element")
            from sage.rings.integer_ring import ZZ
            return ZZ.zero()

        def quo_rem(self, other):
            r"""
            Return the quotient with remainder of the division of this element
            by ``other``.

            INPUT:

            - ``other`` -- an element of the field

            EXAMPLES::

                sage: f,g = QQ(1), QQ(2)
                sage: f.quo_rem(g)
                (1/2, 0)
            """
            if other.is_zero():
                raise ZeroDivisionError
            return (self/other, self.parent().zero())

        def is_unit( self ):
            r"""
            Returns True if ``self`` has a multiplicative inverse.

            EXAMPLES::

                sage: QQ(2).is_unit()
                True
                sage: QQ(0).is_unit()
                False
            """
            return not self.is_zero()

        # Fields are unique factorization domains, so, there is gcd and lcm
        # Of course, in general gcd and lcm in a field are not very interesting.
        # However, they should be implemented!
        @coerce_binop
        def gcd(self,other):
            """
            Greatest common divisor.

            .. NOTE::

                Since we are in a field and the greatest common divisor is only
                determined up to a unit, it is correct to either return zero or
                one. Note that fraction fields of unique factorization domains
                provide a more sophisticated gcd.

            EXAMPLES::

                sage: K = GF(5)
                sage: K(2).gcd(K(1))
                1
                sage: K(0).gcd(K(0))
                0
                sage: all(x.gcd(y) == (0 if x == 0 and y == 0 else 1) for x in K for y in K)
                True

            For field of characteristic zero, the gcd of integers is considered
            as if they were elements of the integer ring::

                sage: gcd(15.0,12.0)
                3.00000000000000

            But for other floating point numbers, the gcd is just `0.0` or `1.0`::

                sage: gcd(3.2, 2.18)
                1.00000000000000

                sage: gcd(0.0, 0.0)
                0.000000000000000

            AUTHOR:

            - Simon King (2011-02) -- :trac:`10771`
            - Vincent Delecroix (2015) -- :trac:`17671`
            """
            P = self.parent()
            try:
                has_zero_char = P.characteristic() == 0
            except (AttributeError, NotImplementedError):
                has_zero_char = False
            if has_zero_char:
                from sage.rings.integer_ring import ZZ
                try:
                    return P(ZZ(self).gcd(ZZ(other)))
                except TypeError:
                    pass

            if self == P.zero() and other == P.zero():
                return P.zero()
            return P.one()

        @coerce_binop
        def lcm(self, other):
            """
            Least common multiple.

            .. NOTE::

                Since we are in a field and the least common multiple is only
                determined up to a unit, it is correct to either return zero or
                one. Note that fraction fields of unique factorization domains
                provide a more sophisticated lcm.

            EXAMPLES::

                sage: GF(2)(1).lcm(GF(2)(0))
                0
                sage: GF(2)(1).lcm(GF(2)(1))
                1

            For field of characteristic zero, the lcm of integers is considered
            as if they were elements of the integer ring::

                sage: lcm(15.0,12.0)
                60.0000000000000

            But for others floating point numbers, it is just `0.0` or `1.0`::

                sage: lcm(3.2, 2.18)
                1.00000000000000

                sage: lcm(0.0, 0.0)
                0.000000000000000

            AUTHOR:

            - Simon King (2011-02) -- :trac:`10771`
            - Vincent Delecroix (2015) -- :trac:`17671`
            """
            P = self.parent()
            try:
                has_zero_char = P.characteristic() == 0
            except (AttributeError, NotImplementedError):
                has_zero_char = False
            if has_zero_char:
                from sage.rings.integer_ring import ZZ
                try:
                    return P(ZZ(self).lcm(ZZ(other)))
                except TypeError:
                    pass

            if self.is_zero() or other.is_zero():
                return P.zero()
            return P.one()

        @coerce_binop
        def xgcd(self, other):
            r"""
            Compute the extended gcd of ``self`` and ``other``.

            INPUT:

            - ``other`` -- an element with the same parent as ``self``

            OUTPUT:

            A tuple ``(r, s, t)`` of elements in the parent of ``self`` such
            that ``r = s * self + t * other``. Since the computations are done
            over a field, ``r`` is zero if ``self`` and ``other`` are zero,
            and one otherwise.

            AUTHORS:

            - Julian Rueth (2012-10-19): moved here from
              :class:`sage.structure.element.FieldElement`

            EXAMPLES::

                sage: K = GF(5)
                sage: K(2).xgcd(K(1))
                (1, 3, 0)
                sage: K(0).xgcd(K(4))
                (1, 0, 4)
                sage: K(1).xgcd(K(1))
                (1, 1, 0)
                sage: GF(5)(0).xgcd(GF(5)(0))
                (0, 0, 0)

            The xgcd of non-zero floating point numbers will be a triple of
            floating points. But if the input are two integral floating points
            the result is a floating point version of the standard gcd on
            `\ZZ`::

                sage: xgcd(12.0, 8.0)
                (4.00000000000000, 1.00000000000000, -1.00000000000000)

                sage: xgcd(3.1, 2.98714)
                (1.00000000000000, 0.322580645161290, 0.000000000000000)

                sage: xgcd(0.0, 1.1)
                (1.00000000000000, 0.000000000000000, 0.909090909090909)
            """
            P = self.parent()
            try:
                has_zero_char = P.characteristic() == 0
            except (AttributeError, NotImplementedError):
                has_zero_char = False
            if has_zero_char:
                from sage.rings.integer_ring import ZZ
                try:
                    return tuple(P(x) for x in ZZ(self).xgcd(ZZ(other)))
                except TypeError:
                    pass

            if not self.is_zero():
                return (P.one(), ~self, P.zero())
            if not other.is_zero():
                return (P.one(), P.zero(), ~other)
            # else both are 0
            return (P.zero(), P.zero(), P.zero())

        def factor(self):
            """
            Return a factorization of ``self``.

            Since ``self`` is either a unit or zero, this function is trivial.

            EXAMPLES::

                sage: x = GF(7)(5)
                sage: x.factor()
                5
                sage: RR(0).factor()
                Traceback (most recent call last):
                ...
                ArithmeticError: factorization of 0.000000000000000 is not defined
            """
            if not self:
                raise ArithmeticError("factorization of {!r} is not defined".format(self))
            from sage.structure.factorization import Factorization
            return Factorization([], self)  # No factor; "self" as unit

        def inverse_of_unit(self):
            r"""
            Return the inverse of this element.

            EXAMPLES::

                sage: NumberField(x^7+2,'a')(2).inverse_of_unit()
                1/2

            Trying to invert the zero element typically raises a
            ``ZeroDivisionError``::

                sage: QQ(0).inverse_of_unit()
                Traceback (most recent call last):
                ...
                ZeroDivisionError: rational division by zero

            To catch that exception in a way that also works for non-units
            in more general rings, use something like::

                sage: try:
                ....:    QQ(0).inverse_of_unit()
                ....: except ArithmeticError:
                ....:    pass

            Also note that some “fields” allow one to invert the zero element::

                sage: RR(0).inverse_of_unit()
                +infinity
            """
            return ~self
