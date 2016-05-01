r"""
Unique factorization domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_attribute import lazy_class_attribute
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_singleton import Category_contains_method_by_parent_class
from sage.categories.gcd_domains import GcdDomains

class UniqueFactorizationDomains(Category_singleton):
    """
    The category of unique factorization domains
    constructive unique factorization domains, i.e. where one can constructively
    factor members into a product of a finite number of irreducible elements

    EXAMPLES::

        sage: UniqueFactorizationDomains()
        Category of unique factorization domains
        sage: UniqueFactorizationDomains().super_categories()
        [Category of gcd domains]

    TESTS::

        sage: TestSuite(UniqueFactorizationDomains()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: UniqueFactorizationDomains().super_categories()
            [Category of gcd domains]
        """
        return [GcdDomains()]

    def additional_structure(self):
        """
        Return whether ``self`` is a structure category.

        .. SEEALSO:: :meth:`Category.additional_structure`

        The category of unique factorization domains does not define
        additional structure: a ring morphism between unique factorization
        domains is a unique factorization domain morphism.

        EXAMPLES::

            sage: UniqueFactorizationDomains().additional_structure()
        """
        return None

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: GF(4, "a") in UniqueFactorizationDomains()
            True
            sage: QQ in UniqueFactorizationDomains()
            True
            sage: ZZ in UniqueFactorizationDomains()
            True
            sage: IntegerModRing(4) in UniqueFactorizationDomains()
            False
            sage: IntegerModRing(5) in UniqueFactorizationDomains()
            True

        This implementation will not be needed anymore once every
        field in Sage will be properly declared in the category
        :class:`UniqueFactorizationDomains`().
        """
        try:
            return self._contains_helper(x) or x.is_unique_factorization_domain()
        except Exception:
            return False

    @lazy_class_attribute
    def _contains_helper(cls):
        """
        Helper for containment tests in the category of unique
        factorization domains.

        This helper just tests whether the given object's category
        is already known to be a sub-category of the category of
        unique factorization domains. There are, however, rings that
        are initialised as plain commutative rings and found out to be
        unique factorization domains only afterwards. Hence, this helper
        alone is not enough for a proper containment test.

        TESTS::

            sage: R = Zmod(7)
            sage: R.category()
            Join of Category of finite commutative rings
                and Category of subquotients of monoids
                and Category of quotients of semigroups
                and Category of finite enumerated sets
            sage: ID = UniqueFactorizationDomains()
            sage: ID._contains_helper(R)
            False
            sage: R in ID  # This changes the category!
            True
            sage: ID._contains_helper(R)
            True
        """
        return Category_contains_method_by_parent_class(cls())

    class ParentMethods:
        def is_unique_factorization_domain(self, proof=True):
            """
            Return True, since this in an object of the category of unique factorization domains.

            EXAMPLES::

                sage: Parent(QQ,category=UniqueFactorizationDomains()).is_unique_factorization_domain()
                True

            """
            return True

        def _gcd_univariate_polynomial(self, f, g):
            """
            Return the greatest common divisor of ``f`` and ``g``.

            INPUT:

            - ``f``, ``g`` -- two polynomials defined over this UFD.

            .. NOTE::

                This is a helper method for
                :meth:`sage.rings.polynomial.polynomial_element.Polynomial.gcd`.

            ALGORITHM:

            Algorithm 3.3.1 in [GTM138]_, based on pseudo-division.

            EXAMPLES::

                sage: R.<x> = PolynomialRing(ZZ, sparse=True)
                sage: S.<T> = R[]
                sage: p = (-3*x^2 - x)*T^3 - 3*x*T^2 + (x^2 - x)*T + 2*x^2 + 3*x - 2
                sage: q = (-x^2 - 4*x - 5)*T^2 + (6*x^2 + x + 1)*T + 2*x^2 - x
                sage: quo,rem=p.pseudo_quo_rem(q); quo,rem
                ((3*x^4 + 13*x^3 + 19*x^2 + 5*x)*T + 18*x^4 + 12*x^3 + 16*x^2 + 16*x,
                 (-113*x^6 - 106*x^5 - 133*x^4 - 101*x^3 - 42*x^2 - 41*x)*T - 34*x^6 + 13*x^5 + 54*x^4 + 126*x^3 + 134*x^2 - 5*x - 50)
                sage: (-x^2 - 4*x - 5)^(3-2+1) * p == quo*q + rem
                True

            REFERENCES:

            .. [GTM138] Henri Cohen. A Course in Computational Number Theory.
               Graduate Texts in Mathematics, vol. 138. Springer, 1993.
            """
            if f.degree() < g.degree():
                A,B = g, f
            else:
                A,B = f, g

            if B.is_zero():
                return A

            a = b = self.zero()
            for c in A.coefficients():
                a = a.gcd(c)
                if a.is_one():
                    break
            for c in B.coefficients():
                b = b.gcd(c)
                if b.is_one():
                    break

            d = a.gcd(b)
            A = A // a
            B = B // b
            g = h = 1

            delta = A.degree()-B.degree()
            _,R = A.pseudo_quo_rem(B)

            while R.degree() > 0:
                A = B
                B = R // (g*h**delta)
                g = A.leading_coefficient()
                h = h*g**delta // h**delta
                delta = A.degree() - B.degree()
                _, R = A.pseudo_quo_rem(B)

            if R.is_zero():
                b = self.zero()
                for c in B.coefficients():
                    b = b.gcd(c)
                    if b.is_one():
                        break

                return d*B // b

            return d

    class ElementMethods:
        # prime?
        # squareFree
        # factor

        def nth_root(self, n):
            r"""
            Return a `n`-th root of this element.

            This generic method relies on factorization.

            EXAMPLES::

                sage: R.<x> = ZZ[]
                sage: a = 27 * (x+3)**6 * (x+5)**3
                sage: a.nth_root(3)
                3*x^3 + 33*x^2 + 117*x + 135

                sage: b = 25 * (x^2 + x + 1)
                sage: b.nth_root(2)
                Traceback (most recent call last):
                ...
                ValueError: 25*x^2 + 25*x + 25 is not a 2nd power
                sage: R(0).nth_root(3)
                0
                sage: R.<x> = QQ[]
                sage: a = 1/4 * (x/7 + 3/2)^2 * (x/2 + 5/3)^4
                sage: a.nth_root(2)
                1/56*x^3 + 103/336*x^2 + 365/252*x + 25/12

                sage: R.<x,y,z> = QQ[]
                sage: a = 32 * (x*y + 1)^5 * (x+y+z)^5
                sage: a.nth_root(5)
                2*x^2*y + 2*x*y^2 + 2*x*y*z + 2*x + 2*y + 2*z
                sage: b = x + 2*y + 3*z
                sage: b.nth_root(42)
                Traceback (most recent call last):
                ...
                ValueError: x + 2*y + 3*z is not a 42th power

                sage: R.<x> = NumberField(QQ['x'].gen()^2 - 2)
                sage: a = (1 + x)^3 * (2*x - 5)^6
                sage: b = a.nth_root(3); b
                13*x - 7
                sage: b^3 == a
                True
            """
            from sage.rings.integer_ring import ZZ

            n = ZZ.coerce(n)

            if n <= 0:
                raise ValueError("n (={}) must be positive".format(n))
            elif n.is_one() or self.is_zero():
                return self
            else:
                f = self.factor()
                u = f.unit()

                if u.is_one():
                    ans = u.base_ring().one()
                else:
                    # here we need to compute a n-th root of the unit ``u``. In
                    # most case, it is not doable. But in some cases (e.g.
                    # polynomial ring) we can hope that the base ring can handle
                    # the computation.
                    R = self.parent()

                    if u.parent() == R:
                        try:
                            u = R.base_ring()(u)
                        except (ValueError, TypeError):
                            pass

                    if u.parent() == R or not hasattr(u, 'nth_root'):
                        # we raise an error as otherwise calling nth_root will
                        # raise an AttributeError or lead to an infinite recursion
                        raise NotImplementedError("nth root not implemented for {}".format(u.parent()))

                    ans = R(u.nth_root(n))

                for (v, exp) in f:
                    if exp % n:
                        if n == 2:
                            postfix = 'nd'
                        elif n == 3:
                            postfix = 'rd'
                        else:
                            postfix = 'th'
                        raise ValueError("{} is not a {}{} power".format(self, n, postfix))
                    ans *= v ** (exp // n)

                return ans
