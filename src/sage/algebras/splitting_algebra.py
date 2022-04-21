# -*- coding: utf-8 -*-
r"""
Splitting Algebras

*Splitting algebras* have been considered by Dan Laksov, Anders Thorup,
Torsten Ekedahl and others (see references below) in order to study
intersection theory of Grassmann and other flag schemes. Similarly as
*splitting fields* they can be considered as extensions of rings containing
all the roots of a given monic polynomial over that ring under the
assumption that its Galois group is the symmetric group of order equal
to the polynomial's degree.

Thus they can be used as a tool to express elements of a ring generated by
`n` indeterminates in terms of symmetric functions in these indeterminates.

This realization of splitting algebras follows the approach of a recursive
quotient ring construction splitting off some linear factor of the
polynomial in each recursive step. Accordingly it is inherited from
:class:`PolynomialQuotientRing_domain`.

AUTHORS:

- Sebastian Oehms (April 2020): initial version
"""

# ****************************************************************************
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from warnings import warn

from sage.misc.verbose import verbose
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_domain
from sage.rings.polynomial.polynomial_quotient_ring_element import PolynomialQuotientRingElement


# ------------------------------------------------------------------------------------------------------------------
# Element class for the splitting algebra
# --------------------------------------------------------------------------------------------------------
class SplittingAlgebraElement(PolynomialQuotientRingElement):
    r"""
    Element class for :class:`SplittingAlgebra`.

    EXAMPLES::

        sage: from sage.algebras.splitting_algebra import SplittingAlgebra
        sage: cp6 = cyclotomic_polynomial(6)
        sage: CR6.<e6> = SplittingAlgebra(cp6)
        sage: type(e6)
        <class 'sage.algebras.splitting_algebra.SplittingAlgebra_with_category.element_class'>

        sage: type(CR6(5))
        <class 'sage.algebras.splitting_algebra.SplittingAlgebra_with_category.element_class'>
    """
    def __invert__(self):
        r"""
        Return the inverse of ``self``.

        Support inversion of special elements attached to the construction
        of the parent and which are recorded in the list
        ``self.parent()._invertible_elements``.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: CR3.<e3> = SplittingAlgebra(cyclotomic_polynomial(3))
            sage: ~e3
            -e3 - 1
            sage: ~(e3 + 5)
            Traceback (most recent call last):
            ...
            NotImplementedError: The base ring (=Integer Ring) is not a field
        """
        inv_elements = self.parent()._invertible_elements
        if self in inv_elements:
            return inv_elements[self]

        return super(SplittingAlgebraElement, self).__invert__()


    def is_unit(self):
        r"""
        Return ``True`` if ``self`` is invertible.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: CR3.<e3> = SplittingAlgebra(cyclotomic_polynomial(3))
            sage: e3.is_unit()
            True
        """
        inv_elements = self.parent()._invertible_elements
        if self in inv_elements:
            return True

        return super(SplittingAlgebraElement, self).is_unit()

    def dict(self):
        r"""
        Return the dictionary of ``self`` according to its lift to the cover.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: CR3.<e3> = SplittingAlgebra(cyclotomic_polynomial(3))
            sage: (e3 + 42).dict()
            {0: 42, 1: 1}
        """
        return self.lift().dict()


# ------------------------------------------------------------------------------------------------------------------
# Parent class of the splitting algebra
# --------------------------------------------------------------------------------------------------------
class SplittingAlgebra(PolynomialQuotientRing_domain):
    r"""
    For a given monic polynomial `p(t)` of degree `n` over a commutative
    ring `R`, the splitting algebra is the universal `R`-algebra in which
    `p(t)` has `n` roots, or, more precisely, over which `p(t)` factors,

    .. MATH::

        p(t) = (t - \xi_1) \cdots (t - \xi_n).

    This class creates an algebra as extension over the base ring of a
    given polynomial `p` such that `p` splits into linear factors over
    that extension. It is assumed (and not checked in general) that the
    Galois group of `p` is the symmetric Group `S(n)`. The construction
    is recursive (following [LT2012]_, 1.3).

    INPUT:

    - ``monic_polynomial`` -- the monic polynomial which should be split
    - ``names``  -- names for the indeterminates to be adjoined to the
      base ring of ``monic_polynomial``
    - ``warning`` -- (default: ``True``) can be used (by setting to ``False``)
      to suppress a warning which will be thrown whenever it cannot be
      checked that the Galois group of ``monic_polynomial`` is maximal

    EXAMPLES::

        sage: from sage.algebras.splitting_algebra import SplittingAlgebra
        sage: Lc.<w> = LaurentPolynomialRing(ZZ)
        sage: PabLc.<u,v> = Lc[]; t = polygen(PabLc)
        sage: S.<x, y> = SplittingAlgebra(t^3 - u*t^2 + v*t - w)
        doctest:...: UserWarning: Assuming x^3 - u*x^2 + v*x - w to have maximal
                                  Galois group!

        sage: roots = S.splitting_roots(); roots
        [x, y, -y - x + u]
        sage: all(t^3 -u*t^2 +v*t -w == 0 for t in roots)
        True
        sage: xi = ~x; xi
        (w^-1)*x^2 + ((-w^-1)*u)*x + (w^-1)*v
        sage: ~xi == x
        True
        sage: ~y
        ((-w^-1)*x)*y + (-w^-1)*x^2 + ((w^-1)*u)*x
        sage: zi = ((w^-1)*x)*y; ~zi
        -y - x + u

        sage: cp3 = cyclotomic_polynomial(3).change_ring(GF(5))
        sage: CR3.<e3> = SplittingAlgebra(cp3)
        sage: CR3.is_field()
        True
        sage: CR3.cardinality()
        25
        sage: F.<a> = cp3.splitting_field()
        sage: F.cardinality()
        25
        sage: E3 = cp3.change_ring(F).roots()[0][0]; E3
        3*a + 3
        sage: f = CR3.hom([E3]); f
        Ring morphism:
          From: Splitting Algebra of x^2 + x + 1
                with roots [e3, 4*e3 + 4]
                over Finite Field of size 5
          To:   Finite Field in a of size 5^2
          Defn: e3 |--> 3*a + 3

    REFERENCES:

    - [EL2002]_
    - [Lak2010]_
    - [Tho2011]_
    - [LT2012]_
    """
    Element = SplittingAlgebraElement

    def __init__(self, monic_polynomial, names='X', iterate=True, warning=True):
        r"""
        Python constructor.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: Lw.<w> = LaurentPolynomialRing(ZZ)
            sage: PuvLw.<u,v> = Lw[]; t = polygen(PuvLw)
            sage: S.<x, y> = SplittingAlgebra(t^3 - u*t^2 + v*t - w, warning=False)
            sage: TestSuite(S).run()
        """

        # ---------------------------------------------------------------------------------
        # checking input parameters
        # ---------------------------------------------------------------------------------

        base_ring = monic_polynomial.base_ring()
        if not monic_polynomial.is_monic():
            raise ValueError("given polynomial must be monic")
        deg = monic_polynomial.degree()

        from sage.structure.category_object import normalize_names
        self._root_names  = normalize_names(deg-1, names)
        root_names = list(self._root_names)
        verbose("Create splitting algebra to base ring %s and polynomial %s (%s %s)"
                % (base_ring, monic_polynomial, iterate, warning))

        self._defining_polynomial = monic_polynomial
        self._iterate = iterate

        try:
            if not base_ring.is_integral_domain():
                raise TypeError("base_ring must be an integral domain")
        except NotImplementedError:
            from sage.rings.ring import Ring
            if not isinstance(base_ring, Ring):
                raise TypeError("base_ring must be an instance of ring")
            if warning:
                warn('Assuming %s to be an integral domain!' % (base_ring))

        if deg < 1:
            raise ValueError("the degree of the polynomial must positive")

        self._splitting_roots     = []
        self._coefficients_list   = []
        self._invertible_elements = {}

        if isinstance(base_ring, SplittingAlgebra):
            self._invertible_elements = base_ring._invertible_elements

        # ------------------------------------------------------------------------------------
        # taking next root_name
        # ------------------------------------------------------------------------------------
        root_name = root_names[0]
        p = monic_polynomial.change_variable_name(root_name)
        P = p.parent()

        self._set_modulus_irreducible_ = False
        try:
            if not p.is_irreducible():
                raise ValueError("monic_polynomial must be irreducible")
        except (NotImplementedError, AttributeError):
            # assuming this has been checked mathematically before
            self._set_modulus_irreducible_ = True
            if warning:
                warn('Assuming %s to have maximal Galois group!' % (monic_polynomial))
                warning = False # one warning must be enough

        verbose("P %s defined:" % (P))

        if deg > 2 and iterate:
            # ------------------------------------------------------------------------------------
            # successive solution via recursion (on base_ring_step)
            # ------------------------------------------------------------------------------------
            base_ring_step = SplittingAlgebra(monic_polynomial, tuple(root_names), iterate=False, warning=False)
            first_root = base_ring_step.gen()

            verbose("base_ring_step %s defined:" % (base_ring_step))

            # ------------------------------------------------------------------------------------
            # splitting first root off
            # ------------------------------------------------------------------------------------
            from copy import copy
            root_names_reduces = copy(root_names)
            root_names_reduces.remove(root_name)

            P = base_ring_step[root_names_reduces[0]]
            p  =  P(monic_polynomial.dict())
            q, r = p.quo_rem( (P.gen()-first_root) )

            verbose("Invoking recursion with: %s" % (q,))

            SplittingAlgebra.__init__(self, q, root_names_reduces, warning=False)

            splitting_roots   = base_ring_step._splitting_roots   + self._splitting_roots
            coefficients_list = base_ring_step._coefficients_list + self._coefficients_list

            verbose("Adding roots: %s" % (splitting_roots))

            self._splitting_roots   = splitting_roots
            self._coefficients_list = coefficients_list
        else:
            PolynomialQuotientRing_domain.__init__(self, P, p, root_name)

            first_root = self.gen()
            self._splitting_roots.append(first_root)
            self._coefficients_list = [monic_polynomial.coefficients(sparse=False)]

            if not iterate:
                verbose("pre ring defined splitting_roots: %s" % (self._splitting_roots))
                return

            verbose("final ring defined splitting_roots: %s" % (self._splitting_roots))

        if deg == 2:
            coefficients = monic_polynomial.coefficients(sparse=False)
            lin_coeff = coefficients[1]
            self._splitting_roots.append(-lin_coeff - first_root)

        self._root_names = names
        self._splitting_roots = [self(root) for root in self._splitting_roots]
        verbose("splitting_roots: %s embedded" % (self._splitting_roots))


        # -------------------------------------------------------------------------------------------
        # try to calculate inverses of the roots. This is possible if the original polynomial
        # has an invertible constant term. For example let cf = [-w, v,-u, 1] that is
        # p = h^3 -u*h^2 + v*h -w, than u = x + y + z, v = x*y + x*z + y*z, w = x*y*z. If
        # w is invertible then 1/x = (v -(u-x)*x)/w, 1/y = (v -(u-y)*y)/w, 1/z = (v -(u-z)*z)/w
        # -------------------------------------------------------------------------------------------
        # first find the polynomial with invertible constant coefficient
        # -------------------------------------------------------------------------------------------
        cf0_inv = None
        for cf in self._coefficients_list:
            cf0 = cf[0]
            try:
                cf0_inv = ~(cf[0])
                cf0_inv = self(cf0_inv)
                verbose("invertible coefficient: %s found" %(cf0_inv))
                break
            except NotImplementedError:
                verbose("constant coefficient: %s not invertibe" %(cf0))


        # ----------------------------------------------------------------------------------
        # assuming that cf splits into linear factors over self and the _splitting_roots
        # are its roots we can calculate inverses
        # ----------------------------------------------------------------------------------
        if cf0_inv is not None:
            deg_cf = len(cf)-1
            pf  =  P(cf)
            for root in self._splitting_roots:
                check = self(pf)
                if not check.is_zero():
                    continue
                root_inv = self.one()
                for pos in range(deg_cf-1 ):
                    root_inv =  (-1 )**(pos+1 ) * cf[deg_cf-pos-1 ] - root_inv * root
                verbose("inverse %s of root %s" % (root_inv, root))
                root_inv = (-1 )**(deg_cf) * cf0_inv * root_inv
                self._invertible_elements.update({root:root_inv})
                verbose("adding inverse %s of root %s" % (root_inv, root))
            invert_items = [(k,v) for k, v in self._invertible_elements.items()]
            for k, v in invert_items:
                self._invertible_elements.update({v: k})
        return


    #######################################################################################################################
    # ---------------------------------------------------------------------------------------------------------------------
    # overloaded inherited methods
    # ---------------------------------------------------------------------------------------------------------------------
    #######################################################################################################################
    def __reduce__(self):
        r"""
        Used in pickling.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: L.<t, u, v, w > = LaurentPolynomialRing(ZZ); x = polygen(L)
            sage: S = SplittingAlgebra(x^4 -t*x^3 - u*x^2 - v*x + w, ('X', 'Y', 'Z'), warning=False)
            sage: S.__reduce__()
            (<class 'sage.algebras.splitting_algebra.SplittingAlgebra_with_category'>,
            (x^4 - t*x^3 - u*x^2 - v*x + w, ('X', 'Y', 'Z'), True, False))
            sage: S.base_ring().__reduce__()
            (<class 'sage.algebras.splitting_algebra.SplittingAlgebra_with_category'>,
            (Y^3 + (X - t)*Y^2 + (X^2 - t*X - u)*Y + X^3 - t*X^2 - u*X - v,
            ('Y', 'Z'),
            False,
            False))

            sage: TestSuite(S).run()
        """
        defining_polynomial = self.defining_polynomial()
        definig_coefficients = self._coefficients_list[0]
        if defining_polynomial.coefficients(sparse=False) != definig_coefficients:
            # case of factorization algebra (intermediate construction step)
            par_pol = self.cover_ring()
            defining_polynomial = par_pol(definig_coefficients)
        return self.__class__, (defining_polynomial, self._root_names, self._iterate, False)


    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: L.<u, v> = PolynomialRing(ZZ)
            sage: t = polygen(L)
            sage: Spl.<S, T> = SplittingAlgebra(t^3 - (u^2-v)*t^2 + (v+u)*t - 1)
            sage: Spl._repr_()
            'Splitting Algebra of x^3 + (-u^2 + v)*x^2 + (u + v)*x - 1
             with roots [S, T, -T - S + u^2 - v]
             over Multivariate Polynomial Ring in u, v over Integer Ring'
            sage: Spl.base_ring()    # indirect doctest
            Factorization Algebra of x^3 + (-u^2 + v)*x^2 + (u + v)*x - 1
             with roots [S] over Multivariate Polynomial Ring in u, v over Integer Ring
        """
        if self.is_completely_split():
            return ('Splitting Algebra of %s with roots %s over %s'
                    % (self.defining_polynomial(), self.splitting_roots(), self.scalar_base_ring()))
        else:
            return ('Factorization Algebra of %s with roots %s over %s'
                    % (self.defining_polynomial(), self.splitting_roots(), self.scalar_base_ring()))

    def _first_ngens(self, n):
        r"""
        Used by the preparser for ``R.<x> = ...``.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: L.<u, v> = PolynomialRing(ZZ)
            sage: t = polygen(L)
            sage: S.<X, Y> = SplittingAlgebra(t^3 - (u^2-v)*t^2 + (v+u)*t - 1)  # indirect doctest
            sage: X.parent()
            Splitting Algebra of x^3 + (-u^2 + v)*x^2 + (u + v)*x - 1
             with roots [X, Y, -Y - X + u^2 - v]
             over Multivariate Polynomial Ring in u, v over Integer Ring
            sage: S._first_ngens(4)
            (X, Y, u, v)
        """
        srts = self.splitting_roots()
        k = len(srts)-1
        gens = srts[:k] + list(self.scalar_base_ring().gens())
        return tuple(gens[:n])

    def _element_constructor_(self, x):
        r"""
        Make sure ``x`` is a valid member of ``self``, and return the constructed element.

        TESTS::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: L.<u, v, w> = LaurentPolynomialRing(ZZ); x = polygen(L)
            sage: S.<X, Y> = SplittingAlgebra(x^3 - u*x^2 + v*x - w)
            sage: S(u + v)
            u + v
            sage: S(X*Y + X)
            X*Y + X
            sage: TestSuite(S).run()                   # indirect doctest
        """
        if isinstance(x, SplittingAlgebraElement):
            # coercion from covering fixes pickling problems
            return self(x.lift())
        return super(SplittingAlgebra, self)._element_constructor_(x)

    def hom(self, im_gens, codomain=None, check=True, base_map=None):
        r"""
        This version keeps track with the special recursive structure
        of :class:`SplittingAlgebra`

        Type ``Ring.hom?`` to see the general documentation of this method.
        Here you see just special examples for the current class.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: L.<u, v, w> = LaurentPolynomialRing(ZZ); x = polygen(L)
            sage: S = SplittingAlgebra(x^3 - u*x^2 + v*x - w, ('X', 'Y'))
            sage: P.<x, y, z> = PolynomialRing(ZZ)
            sage: F = FractionField(P)
            sage: im_gens = [F(g) for g in [y, x, x + y + z, x*y+x*z+y*z, x*y*z]]
            sage: f = S.hom(im_gens)
            sage: f(u), f(v), f(w)
            (x + y + z, x*y + x*z + y*z, x*y*z)
            sage: roots = S.splitting_roots(); roots
            [X, Y, -Y - X + u]
            sage: [f(r) for r in roots]
            [x, y, z]
        """
        base_ring = self.base_ring()

        if not isinstance(im_gens, (list,tuple)):
            im_gens = [im_gens]

        all_gens = self.gens_dict_recursive()
        if len(im_gens) != len(all_gens):
            return super(SplittingAlgebra, self).hom(im_gens, codomain=codomain, check=check, base_map=base_map)

        num_gens = len(self.gens())
        im_gens_start = [img for img in im_gens if im_gens.index(img) <  num_gens]
        im_gens_end   = [img for img in im_gens if im_gens.index(img) >= num_gens]

        if not im_gens_end:
            return super(SplittingAlgebra, self).hom(im_gens, codomain=codomain, check=check, base_map=base_map)

        verbose('base %s im_gens_end %s codomain %s check %s base_map %s' % (base_ring, im_gens_end, codomain, check, base_map))
        hom_on_base_recurs = base_ring.hom(im_gens_end, codomain=codomain, check=check, base_map=base_map)
        verbose('hom_on_base_recurs %s' % (hom_on_base_recurs))

        cover_ring = self.cover_ring()
        hom_from_cover = cover_ring.hom(im_gens_start, codomain=codomain, check=check, base_map=hom_on_base_recurs)
        lift = self.lifting_map()
        return hom_from_cover*lift



    #######################################################################################################################
    # ---------------------------------------------------------------------------------------------------------------------
    # local methods
    # ---------------------------------------------------------------------------------------------------------------------
    #######################################################################################################################


    #######################################################################################################################
    # ---------------------------------------------------------------------------------------------------------------------
    # global methods
    # ---------------------------------------------------------------------------------------------------------------------
    #######################################################################################################################

    def is_completely_split(self):
        r"""
        Return True if the defining polynomial of ``self`` splits into linear factors over ``self``.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: L.<u, v, w > = LaurentPolynomialRing(ZZ); x = polygen(L)
            sage: S.<a,b> = SplittingAlgebra(x^3 - u*x^2 + v*x - w)
            sage: S.is_completely_split()
            True
            sage: S.base_ring().is_completely_split()
            False
        """
        return len(self.splitting_roots()) >= self.defining_polynomial().degree()

    @cached_method
    def lifting_map(self):
        r"""
        Return a section map from ``self`` to the cover ring. It is implemented according
        to the same named method of :class:`~sage.rings.quotient_ring.QuotientRing_nc`.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: x = polygen(ZZ)
            sage: S = SplittingAlgebra(x^2+1, ('I',))
            sage: lift = S.lifting_map()
            sage: lift(5)
            5
            sage: r1, r2 =S.splitting_roots()
            sage: lift(r1)
            I
        """
        from sage.rings.morphism import RingMap_lift
        return RingMap_lift(self, self.cover_ring())

    def splitting_roots(self):
        r"""
        Return the roots of the split equation.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: x = polygen(ZZ)
            sage: S = SplittingAlgebra(x^2+1, ('I',))
            sage: S.splitting_roots()
            [I, -I]
        """
        return self._splitting_roots

    @cached_method
    def scalar_base_ring(self):
        r"""
        Return the ring of scalars of ``self`` (considered as an algebra)

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: L.<u, v, w > = LaurentPolynomialRing(ZZ)
            sage: x = polygen(L)
            sage: S = SplittingAlgebra(x^3 - u*x^2 + v*x - w, ('X', 'Y'))
            sage: S.base_ring()
            Factorization Algebra of x^3 - u*x^2 + v*x - w with roots [X]
             over Multivariate Laurent Polynomial Ring in u, v, w over Integer Ring
            sage: S.scalar_base_ring()
            Multivariate Laurent Polynomial Ring in u, v, w over Integer Ring
        """
        base_ring = self.base_ring()
        if isinstance(base_ring, SplittingAlgebra):
            if base_ring.is_completely_split():
                # another splitting algebra independent of self
                return base_ring
            else:
                return base_ring.scalar_base_ring()
        return base_ring

    @cached_method
    def defining_polynomial(self):
        r"""
        Return the defining polynomial of ``self``.

        EXAMPLES::

            sage: from sage.algebras.splitting_algebra import SplittingAlgebra
            sage: L.<u, v, w > = LaurentPolynomialRing(ZZ)
            sage: x = polygen(L)
            sage: S = SplittingAlgebra(x^3 - u*x^2 + v*x - w, ('X', 'Y'))
            sage: S.defining_polynomial()
            x^3 - u*x^2 + v*x - w
        """
        base_ring = self.base_ring()
        if isinstance(base_ring, SplittingAlgebra):
            if base_ring.is_completely_split():
                # another splitting algebra independent of self
                return self._defining_polynomial
            else:
                return base_ring.defining_polynomial()
        return self._defining_polynomial


# --------------------------------------------------------------------------------------------
# ============================================================================================
# Utility function to create the roots of a polynomial in an appropriate extension ring
# ============================================================================================
# --------------------------------------------------------------------------------------------

def solve_with_extension(monic_polynomial, root_names=None, var='x', flatten=False, warning=True):
    r"""
    Return all roots of a monic polynomial in its base ring or in an appropriate
    extension ring, as far as possible.

    INPUT:

    - ``monic_polynomial`` -- the monic polynomial whose roots should be created
    - ``root_names``  -- names for the indeterminates needed to define the
      splitting algebra of the ``monic_polynomial`` (if necessary and possible)
    - ``var``  -- (default: ``'x'``) for the indeterminate needed to define the
      splitting field of the ``monic_polynomial`` (if necessary and possible)
    - ``flatten`` -- (default: ``True``) if ``True`` the roots will not be
      given as a list of pairs ``(root, multiplicity)`` but as a list of
      roots repeated according to their multiplicity
    - ``warning`` -- (default: ``True``) can be used (by setting to ``False``)
      to suppress a warning which will be thrown whenever it cannot be checked
      that the Galois group of ``monic_polynomial`` is maximal

    OUTPUT:

    List of tuples ``(root, multiplicity)`` respectively list of roots repeated
    according to their multiplicity if option ``flatten`` is ``True``.

    EXAMPLES::

        sage: from sage.algebras.splitting_algebra import solve_with_extension
        sage: t = polygen(ZZ)
        sage: p = t^2 -2*t +1
        sage: solve_with_extension(p, flatten=True )
        [1, 1]
        sage: solve_with_extension(p)
        [(1, 2)]

        sage: cp5 = cyclotomic_polynomial(5, var='T').change_ring(UniversalCyclotomicField())
        sage: solve_with_extension(cp5)
        [(E(5), 1), (E(5)^4, 1), (E(5)^2, 1), (E(5)^3, 1)]
        sage: _[0][0].parent()
        Universal Cyclotomic Field
    """
    def create_roots(monic_polynomial, warning=True):
        r"""
        This internal function creates all roots of a polynomial in an
        appropriate extension ring assuming that none of the roots is
        contained its base ring.

        It first tries to create the splitting field of the given polynomial.
        If this is not faithful the splitting algebra will be created.

        INPUT:

        - ``monic_polynomial`` -- the monic polynomial whose roots should
          be created
        - ``warning`` -- (default: ``True``) can be used (by setting to ``False``)
          to suppress a warning which will be thrown whenever it cannot be
          checked that the Galois group of ``monic_polynomial`` is maximal
        """
        parent = monic_polynomial.parent()
        base_ring = parent.base_ring()

        try:
            ext_field, embed = monic_polynomial.splitting_field(var, map=True)

            if embed.domain() != base_ring:
                # in this case the SplittingAlgebra is preferred
                raise NotImplementedError

            # -------------------------------------------------------------------------------------
            # in some cases the embedding of the base_ring in ext_field can not be obtained
            # as coercion
            # -------------------------------------------------------------------------------------
            reset_coercion = False
            from sage.rings.number_field.number_field import NumberField_generic
            if isinstance(base_ring, NumberField_generic):
                reset_coercion = True
            elif base_ring.is_finite() and not base_ring.is_prime_field():
                reset_coercion = True
            if reset_coercion:
                ext_field._unset_coercions_used()
                ext_field.register_coercion(embed)
                ext_field.register_conversion(embed)

            verbose("splitting field %s defined" % (ext_field))
            pol_emb = monic_polynomial.change_ring(ext_field)
            roots = pol_emb.roots()
        except NotImplementedError:
            ext_ring = SplittingAlgebra(monic_polynomial, name_list, warning=warning)
            verbose("splitting algebra %s defined" % (ext_ring))
            roots = [(r, 1) for r in ext_ring.splitting_roots()]
        return roots


    deg_pol = monic_polynomial.degree()
    if not root_names:
        from sage.structure.category_object import normalize_names
        root_names = normalize_names(deg_pol-1, 'r')
    name_list = list(root_names)
    root_list = []
    try:
        root_list = monic_polynomial.roots()
    except (TypeError, ValueError, NotImplementedError):
        pass

    if not root_list:
        # ------------------------------------------------------------------------------
        # no roots found: find roots in an appropriate extension ring
        # ------------------------------------------------------------------------------
        verbose("no roots in base_ring")
        if len(name_list) > deg_pol -1:
            name_list = [name_list[i] for i in range(deg_pol-1)]
        roots = create_roots(monic_polynomial, warning=warning)

    else:
        # ------------------------------------------------------------------
        # root calculation was possible but maybe some more roots in
        # an appropriate extension ring can be constructed.
        # ------------------------------------------------------------------
        num_roots = sum(m for r, m in root_list)
        if num_roots < deg_pol:
            h = monic_polynomial.variables()[0]
            divisor = monic_polynomial.base_ring().one()
            for r, m in root_list:
                divisor *= (h - r)**m
            q, r = monic_polynomial.quo_rem(divisor)
            if len(name_list) > deg_pol - num_roots - 1:
                name_list = [name_list[i] for i in range(deg_pol - num_roots - 1)]
            verbose("%d root found in base ring, now solving %s" % (num_roots,q))
            missing_roots = create_roots(q, warning=True)
            roots = root_list + missing_roots
        else:
            roots = root_list
            verbose("all roots in base ring")

    if flatten:
        from sage.misc.flatten import flatten
        return flatten([[rt]*m for rt, m in roots])
    return roots
