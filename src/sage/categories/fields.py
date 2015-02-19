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

import sage.rings.ring
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

        The following tests against a memory leak fixed in :trac:`13370`::

            sage: import gc
            sage: _ = gc.collect()
            sage: n = len([X for X in gc.get_objects() if isinstance(X, sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic)])
            sage: for i in prime_range(100):
            ...     R = ZZ.quotient(i)
            ...     t = R in Fields()
            sage: _ = gc.collect()
            sage: len([X for X in gc.get_objects() if isinstance(X, sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic)]) - n
            1

        """
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
            Join of Category of integral domains
             and Category of commutative algebras over Rational Field
             and Category of subquotients of monoids
             and Category of quotients of semigroups
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

        def _gcd_univariate_polynomial(self, f, g):
            """
            Return the greatest common divisor of ``f`` and ``g``, as a
            monic polynomial.

            INPUT:

                - ``f``, ``g`` -- two polynomials defined over ``self``

            .. NOTE::

                This is a helper method for
                :meth:`sage.rings.polynomial.polynomial_element.Polynomial.gcd`.

            EXAMPLES::

                sage: R.<x> = QQbar[]
                sage: QQbar._gcd_univariate_polynomial(2*x,2*x^2)
                x

            """
            ret = EuclideanDomains().ElementMethods().gcd(f,g)
            c = ret.leading_coefficient()
            if c.is_unit():
                return (1/c)*ret
            return ret

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
            else: raise NotImplementedError

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
            descriped in [Yun]_. Other fields may provide their own
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

            REFERENCES:

            .. [Yun] Yun, David YY. On square-free decomposition algorithms.
               In Proceedings of the third ACM symposium on Symbolic and algebraic
               computation, pp. 26-35. ACM, 1976.

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

        def __pow__(self, n):
            r"""
            Returns the vector space of dimension `n` over ``self``.

            EXAMPLES::

                sage: QQ^4
                Vector space of dimension 4 over Rational Field
            """
            from sage.modules.all import FreeModule
            return FreeModule(self, n)

    class ElementMethods:
        def euclidean_degree(self):
            r"""
            Return the degree of this element as an element of a euclidean
            domain.

            In a field, this returns 0 for all but the zero element (for
            which it is undefined).

            EXAMPLES::

                sage: QQ.one().euclidean_degree()
                0
            """
            if self.is_zero():
                raise ValueError("euclidean degree not defined for the zero element")
            from sage.rings.all import ZZ
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
        def gcd(self,other):
            """
            Greatest common divisor.

            NOTE:

            Since we are in a field and the greatest common divisor is
            only determined up to a unit, it is correct to either return
            zero or one. Note that fraction fields of unique factorization
            domains provide a more sophisticated gcd.

            EXAMPLES::

                sage: GF(5)(1).gcd(GF(5)(1))
                1
                sage: GF(5)(1).gcd(GF(5)(0))
                1
                sage: GF(5)(0).gcd(GF(5)(0))
                0

            For fields of characteristic zero (i.e., containing the
            integers as a sub-ring), evaluation in the integer ring is
            attempted. This is for backwards compatibility::

                sage: gcd(6.0,8); gcd(6.0,8).parent()
                2
                Integer Ring

            If this fails, we resort to the default we see above::

                sage: gcd(6.0*CC.0,8*CC.0); gcd(6.0*CC.0,8*CC.0).parent()
                1.00000000000000
                Complex Field with 53 bits of precision

            AUTHOR:

            - Simon King (2011-02): Trac ticket #10771

            """
            P = self.parent()
            try:
                other = P(other)
            except (TypeError, ValueError):
                raise ArithmeticError("The second argument can not be interpreted in the parent of the first argument. Can't compute the gcd")
            from sage.rings.integer_ring import ZZ
            if ZZ.is_subring(P):
                try:
                    return ZZ(self).gcd(ZZ(other))
                except TypeError:
                    pass
            # there is no custom gcd, so, we resort to something that always exists
            # (that's new behaviour)
            if self==0 and other==0:
                return P.zero()
            return P.one()

        def lcm(self,other):
            """
            Least common multiple.

            NOTE:

            Since we are in a field and the least common multiple is
            only determined up to a unit, it is correct to either return
            zero or one. Note that fraction fields of unique factorization
            domains provide a more sophisticated lcm.

            EXAMPLES::

                sage: GF(2)(1).lcm(GF(2)(0))
                0
                sage: GF(2)(1).lcm(GF(2)(1))
                1

            If the field contains the integer ring, it is first
            attempted to compute the gcd there::

                sage: lcm(15.0,12.0); lcm(15.0,12.0).parent()
                60
                Integer Ring

            If this fails, we resort to the default we see above::

                sage: lcm(6.0*CC.0,8*CC.0); lcm(6.0*CC.0,8*CC.0).parent()
                1.00000000000000
                Complex Field with 53 bits of precision
                sage: lcm(15.2,12.0)
                1.00000000000000

            AUTHOR:

            - Simon King (2011-02): Trac ticket #10771

            """
            P = self.parent()
            try:
                other = P(other)
            except (TypeError, ValueError):
                raise ArithmeticError("The second argument can not be interpreted in the parent of the first argument. Can't compute the lcm")
            from sage.rings.integer_ring import ZZ
            if ZZ.is_subring(P):
                try:
                    return ZZ(self).lcm(ZZ(other))
                except TypeError:
                    pass
            # there is no custom lcm, so, we resort to something that always exists
            if self==0 or other==0:
                return P.zero()
            return P.one()

        @coerce_binop
        def xgcd(self, other):
            """
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

                sage: (1/2).xgcd(2)
                (1, 2, 0)
                sage: (0/2).xgcd(2)
                (1, 0, 1/2)
                sage: (0/2).xgcd(0)
                (0, 0, 0)
            """
            R = self.parent()
            if not self.is_zero():
                return (R.one(), ~self, R.zero())
            if not other.is_zero():
                return (R.one(), R.zero(), ~other)
            # else both are 0
            return (R.zero(), R.zero(), R.zero())

