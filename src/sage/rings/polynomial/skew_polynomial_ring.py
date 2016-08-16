r"""
Skew Univariate Polynomial Rings

This module provides the :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general`
which constructs a general dense skew univariate polynomials over commutative base rings with
automorphisms over the base rings. This is the set of formal polynomials where the coefficients
are written on the left of the variable of the skew polynomial ring. The modified multiplication
operation over elements of the base ring is extended to all elements of the skew poynomial ring
by associativity and distributivity.

This module also provides :class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_finite_field`
which is a specialized class for skew polynomial rings over finite fields. It inherits from
:class:`~sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general` but allows for
the more efficient computations in the case of finite fields.

AUTHOR:

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests and refactored classes and methods

- Johan Rosenkilde (2016-08-03): changes for bug fixes, docstring and doctest errors

"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.functions.all import log
from sage.misc.cachefunc import cached_method, cached_function
import sage.algebras.algebra
import sage.categories.basic as categories
from sage.rings.integer import Integer
from sage.structure.category_object import normalize_names
from sage.misc.prandom import randint
from sage.rings.ring import Field
from sage.categories.morphism import Morphism
from sage.categories.morphism import IdentityMorphism
import sage.misc.latex as latex
from sage.rings.polynomial.skew_polynomial_element import SkewPolynomial

#########################################################################################


def _base_ring_to_fraction_field(S):
    """
    TODO
    Return the indeterminate generator, one and twist map of
    ``self``, or of the skew polynomial ring `S` associated
    with the fraction field of the base ring of ``self`` if
    the base ring of ``self`` is not a field (internal method).

    OUTPUT:

    - ``x`` -- indeterminate generator of ``self`` or `S`

    - ``one`` -- one of ``self`` or `S`

    - ``sigma`` -- twist map of ``self`` or `S`

    EXAMPLES:

        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: T.<x> = R['x', sigma]
        sage: points = [ (t, 3*t^2 + 4*t + 4) , (t^2, 4*t) ]
        sage: p = T.lagrange_polynomial(points); p #indirect doctest
        ((3*t^3 + 4*t^2)/(-t - 1))*x + (-3*t^3 - 7*t^2 - 4*t - 4)/-t
    """
    R = S.base_ring()
    if isinstance(R, Field):
        return S
    else:
        Q = R.fraction_field()
        gens = R.gens()
        sigmaS = S.twist_map()
        # try:
        sigmaQ = Q.hom([ Q(sigmaS(g)) for g in gens ])
        return Q[S.variable_name(), sigmaQ]
        # except Exception, e:
        #     raise ValueError("unable to lift the twist map to a twist map over %s (error was: %s)" % (Q, e))


def _minimal_vanishing_polynomial(R, eval_pts):
    """
    Return the minimal vanishing polynomial (internal method).

    See the documentation for
    :meth:`SkewPolynomialRing.minimal_vanishing_polynomial` for a description.

    INPUT:

    - ``R`` -- A skew polynomial ring over a field.

    - ``eval_pts`` -- a frozen set TODO of evaluation points which are linearly
        independent over the fixed field of the twist map of the associated
        skew polynomial ring

    OUTPUT:

    The minimal vanishing polynomial.

    EXAMPLES:

        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]
        sage: eval_pts = [1, t, t^2]
        sage: b = S.minimal_vanishing_polynomial(eval_pts); b #indirect doctest
        x^3 + 4
    """
    l = len(eval_pts)
    if l == 0:
        return R.one()
    elif l == 1:
        if eval_pts[0].is_zero():
            return R.one()
        else:
            return R.gen() - (R.twist_map()(eval_pts[0]) / eval_pts[0])
    else:
        t = l//2
        A = eval_pts[:t]
        B = eval_pts[t:]
        M_A = _minimal_vanishing_polynomial(R, A)
        M_A_B = M_A.multi_point_evaluation(B)
        M_M_A_B = _minimal_vanishing_polynomial(R, M_A_B)
        return M_M_A_B * M_A


def _lagrange_interpolation(R, eval_pts, values):
    """
    Return the lagrange polynomial (internal method).

    TODO

    TODO: args should be like lagrange_polynomial. Name as well?

    INPUT:

    - ``x`` -- the generator of the associated skew polynomial ring

    - ``one`` -- the one of the associated skew polynomial ring

    - ``sigma`` -- the twist map of the associated skew polynomial ring

    - ``eval_pts`` -- vector of evaluation points which are linearly
        independent over the fixed field of the twist map of the associated
        skew polynomial ring

    - ``values`` -- list of values that the lagrange polynomial takes
        at the respective `eval_pts`

    OUTPUT:

    The lagrange polynomial.

    EXAMPLES:

        sage: k.<t> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x',Frob]
        sage: points = [ (t, 3*t^2 + 4*t + 4) , (t^2, 4*t) ]
        sage: d = S.lagrange_polynomial(points); d #indirect doctest
        x + t
    """
    l = len(eval_pts)
    if l == 1:
        return (values[0]/eval_pts[0])*R.one()
    else:
        t = l//2
        A = eval_pts[:t]
        B = eval_pts[t:]
        M_A = _minimal_vanishing_polynomial(R, A)
        M_B = _minimal_vanishing_polynomial(R, B)
        A_ = M_B.multi_point_evaluation(A)
        B_ = M_A.multi_point_evaluation(B)
        I_1 = _lagrange_interpolation(R, A_, values[:t])
        I_2 = _lagrange_interpolation(R, B_, values[t:])
        return I_1 * M_B + I_2 * M_A


#########################################################################################


class SkewPolynomialRing_general(sage.algebras.algebra.Algebra,UniqueRepresentation):
    """
    General Skew Univariate polynomial ring over a ring.

    DEFINITION:

    Let `R` be a commutative ring and let `\sigma` be an automorphism over `R`. An
    automorphism (also called twist map) is a structure preserving map from a
    mathematical object (in this case, R) onto itself that also admits an inverse.
    The ring of skew polynomials over an indeterminate variable `X` is defined then,
    as the ring structure on the set `R` as:
    `R[X, \sigma] = { a_{n-1}X^{n-1} + ... + a_{1}X + a_{0} | a_{i} \in R and n \in N }`
    where the addition operation on `R[X, \sigma]` is given by the usual abelian
    group polynomial addition rule and the multiplication operation is defined
    by the modified rule `X*a = \sigma(a)X`.

    This ring is non-commutative and its elements are all such skew polynomials whose
    coefficients come from `R`.

    Reference: "Theory of Non-Commutative Polynomials" - Oystein Ore

    .. TODO::

        Add derivations.

    EXAMPLES::

        sage: R.<t> = ZZ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = SkewPolynomialRing(R,sigma); S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

    One can also use a shorter syntax::

        sage: S.<x> = R['x',sigma]; S
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1

    Be careful, with the latter syntax, one cannot omit the name of the
    variable neither in LHS nor in RHS. If we omit it in LHS, the variable
    is not created::

        sage: Sy = R['y',sigma]; Sy
        Skew Polynomial Ring in y over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        sage: y.parent()
        Traceback (most recent call last):
        ...
        NameError: name 'y' is not defined

    If we omit it in RHS, sage tries to create a polynomial ring and fails::

        sage: Sz.<z> = R[sigma]
        Traceback (most recent call last):
        ...
        ValueError: variable name 'Ring endomorphism of Univariate Polynomial Ring in t over Integer Ring\n
            Defn: t |--> t + 1' is not alphanumeric

    As for polynomials, skew polynomial rings with different variable names
    are not equal::

        sage: R['x',sigma] == R['y',sigma]
        False

    Of course, skew polynomial rings with different twist maps are not
    equal as well::

        sage: R['x',sigma] == R['x',sigma^2]
        False

    Saving and loading of polynomial rings works::

        sage: loads(dumps(R['x',sigma])) == R['x',sigma]
        True

    There is a coercion map from the base ring of the skew polynomial rings::

        sage: S.has_coerce_map_from(R)
        True
        sage: x.parent()
        Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        sage: t.parent()
        Univariate Polynomial Ring in t over Integer Ring
        sage: y = x+t; y
        x + t
        sage: y.parent() is S
        True

    .. SEE ALSO::

        :meth:`sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing`
        :mod:`sage.rings.polynomial.skew_polynomial_element`
    """
    @staticmethod
    def __classcall__(cls, base_ring, map, name=None, sparse=False, element_class=None):
        if not element_class:
            if sparse:
                raise NotImplementedError("sparse skew polynomials are not implemented")
            else:
                from sage.rings.polynomial import skew_polynomial_element
                element_class = skew_polynomial_element.SkewPolynomial_generic_dense
        return super(SkewPolynomialRing_general,cls).__classcall__(cls,base_ring,map,name,sparse,element_class)

    def __init__(self, base_ring, map, name, sparse, element_class):
        """
        This method is a constructor for a general, dense univariate skew polynomial ring.

        INPUT::

        - ``base_ring`` -- a commutative ring

        - ``map`` -- an automorphism of the base ring

        - ``name`` -- string or list of strings representing the name of the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``element_class`` -- class representing the type of element to be used in ring

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma); S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: S([1]) + S([-1])
            0
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<x> = k['x', Frob]; T
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        category = categories.Rings()
        self.__is_sparse = sparse
        self._polynomial_class = element_class
        if map is None:
            self._map = IdentityMorphism(base_ring)
        else:
            if isinstance (map, Morphism):
                if map.domain () == base_ring and map.codomain () == base_ring:
                    self._map = map
                else:
                    raise TypeError("given map is not an automorphism of %s" % base_ring)
            else:
                raise TypeError("given map is not a ring homomorphism")
        self._maps = { 0:IdentityMorphism(base_ring), 1:self._map }
        self._center = { }
        self._center_variable = None
        self._no_generic_basering_coercion = True
        sage.algebras.algebra.Algebra.__init__(self, base_ring, names=name, normalize=True, category=category)
        self.__generator = self._polynomial_class(self, [0,1], is_gen=True)
        base_inject = sage.rings.polynomial.skew_polynomial_element.SkewPolynomialBaseringInjection(base_ring,self)
        self._populate_coercion_lists_(
                coerce_list = [base_inject],
                convert_list = [list, base_inject])

    def __reduce__(self):
        """
        Return the globally unique skew polynomial ring based on
        given arguments
        """
        import sage.rings.polynomial.skew_polynomial_ring_constructor
        return (sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing,
                (self.base_ring(), self.twist_map(), self.variable_name(), self.is_sparse()))

    def _element_constructor_(self, x=None, check=True, is_gen = False, construct=False, **kwds):
        """
        Convert ``x`` into an element of this univariate skew polynomial ring,
        possibly non-canonically.

        INPUT:

        - ``x`` -- an element of the base ring of ``self`` or a ring that
          has a coerce map from ``self`` (default: ``None``).

        - ``check`` -- boolean (default: ``True``)

        - ``is_gen`` -- boolean (default: ``False``)

        - ``construct`` -- boolean (default: ``False``)

        OUTPUT:

        An element of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S(1 + x + x^2 + x^3)
            x^3 + x^2 + x + 1
            sage: S(1 + t)
            t + 1
            sage: S(1 + t).degree()
            0
            sage: S(0).list()
            []
        """
        C = self._polynomial_class
        if isinstance(x, list):
            return C(self, x, check=check, is_gen=False,construct=construct)
        if isinstance(x, Element):
            P = x.parent()
            def build(check):
                if x.is_zero():
                    return P.zero()
                else:
                    return C(self, [x], check=check, is_gen=False, construct=construct)
            if P is self:
                return x
            elif P is self.base_ring():
                build(False)
            elif P == self.base_ring() or self.base_ring().has_coerce_map_from(P):
                build(True)
        try:
            return x._polynomial_(self)
        except AttributeError:
            pass
        if isinstance(x,str):
            try:
                from sage.misc.parser import Parser, LookupNameMaker
                R = self.base_ring()
                p = Parser(Integer, R, LookupNameMaker({self.variable_name(): self.gen()}, R))
                return self(p.parse(x))
            except NameError:
                raise TypeError("unable to coerce string")
        return C(self, x, check, is_gen, construct=construct, **kwds)

    def _coerce_map_from_(self, P):
        """
        Check whether ``self`` has a coerce map from ``P``.

        The rings that canonically coerce into this ring are:

        - this ring itself

        - any ring that canonically coerces to the base ring of this ring

        - skew polynomial rings in the same variable and automorphism over
          any base ring that canonically coerces to the base ring of this ring

        INPUT:

        - ``P`` -- a ring.

        OUTPUT:

        Return ``True`` or ``False``.

        .. NOTE::

            Sparse skew polynomials are not implemented.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.has_coerce_map_from(S)
            True
            sage: S.has_coerce_map_from(R)
            True
            sage: S.has_coerce_map_from(ZZ)
            True
            sage: S.has_coerce_map_from(GF(5^3))
            False

            sage: S.coerce_map_from(ZZ)
            Composite map:
                From: Integer Ring
                To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
                Defn:   Polynomial base injection morphism:
                        From: Integer Ring
                        To:   Univariate Polynomial Ring in t over Integer Ring
                    then
                        Skew Polynomial base injection morphism:
                        From: Univariate Polynomial Ring in t over Integer Ring
                        To:   Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
            sage: S.coerce_map_from(S)
            Identity endomorphism of Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        try:
            connecting = self.base_ring().coerce_map_from(P)
            if connecting is not None:
                return self.coerce_map_from(self.base_ring()) * connecting
        except TypeError:
            pass
        try:
            if isinstance(P, SkewPolynomialRing_general):
                if self.__is_sparse and not P.is_sparse():
                    return False
                if P.variable_name() == self.variable_name():
                    if P.base_ring() is self.base_ring() and \
                            self.base_ring() is ZZ_sage:
                       if self._implementation_names == ('NTL',):
                            return False
                    return self.base_ring().has_coerce_map_from(P.base_ring())
        except AttributeError:
            pass

    def _repr_(self):
        """
        Return string representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Integer Ring twisted by t |--> t + 1
        """
        s = "Skew Polynomial Ring in %s over %s twisted by %s"%(self.variable_name(), self.base_ring(), self._map._repr_short())
        if self.is_sparse():
            s = "Sparse " + s
        return s

    def _latex_(self):
        """
        Return latex representation of ``self``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: latex(S)
            \Bold{Z}[t][x,\begin{array}{l}
            \text{\texttt{Ring{ }endomorphism{ }of{ }Univariate{ }Polynomial{ }Ring{ }in{ }t{ }over{ }Integer{ }Ring}}\\
            \text{\texttt{{ }{ }Defn:{ }t{ }|{-}{-}>{ }t{ }+{ }1}}
            \end{array}]
        """
        return "%s[%s,%s]"%(latex.latex(self.base_ring()), self.latex_variable_names()[0], latex.latex(self._map))

    def change_var(self, var):
        r"""
        Return the skew polynomial ring in variable ``var`` over the same base
        ring.

        INPUT:

        - ``var`` -- a string representing the name of the new variable of ``self``

        OUTPUT:

        ``self`` with variable name name changed to ``var``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: R.<x> = SkewPolynomialRing(k,Frob); R
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: R.change_var('y')
            Skew Polynomial Ring in y over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        from sage.rings.polynomial.skew_polynomial_ring_constructor import SkewPolynomialRing
        return SkewPolynomialRing(self.base_ring(), self.twist_map(), names = var, sparse=self.is_sparse())

    def characteristic(self):
        """
        Return the characteristic of the base ring of this skew polynomial ring.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: R['x',sigma].characteristic()
            0

            sage: k.<u> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: k['y',Frob].characteristic()
            5
        """
        return self.base_ring().characteristic()

    def twist_map(self, n=1):
        """
        Return the twist map, otherwise known as the automorphism over the base ring of
        ``self``, iterated `n` times.

        INPUT:

        -  ``n`` - a relative integer (default: 1)

        OUTPUT:

        -  The `n`-th iterative of the twist map of this skew polynomial ring.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.twist_map()
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 1
            sage: S.twist_map() == sigma
            True
            sage: S.twist_map(10)
            Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 10

        If `n` in negative, Sage tries to compute the inverse of the twist map.
        Sometimes it fails (even if the twist map is actually invertible)::

            sage: S.twist_map(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twist map Ring endomorphism of Univariate Polynomial Ring in t over Rational Field
                  Defn: t |--> t + 1

        Sometimes it succeeds::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<y> = k['y',Frob]
            sage: T.twist_map(-1)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3
        """
        try:
            return self._maps[n]
        except KeyError:
            if n >= 0:
                map = self._map**n
                self._maps[n] = map
                return map
            else:
                try:
                    map = self._map**n
                except TypeError:
                    raise NotImplementedError("inversion of the twist map %s" % self._map)
                self._maps[n] = map
                return map

    def gen(self, n=0):
        """
        Return the indeterminate generator of this skew polynomial ring.

        INPUT:

        - ``n`` -- an integer (default: 0)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]; S
            Skew Polynomial Ring in x over Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
            sage: y = S.gen(); y
            x
            sage: y == x
            True
            sage: y is x
            True
            sage: y.is_gen()
            True
        """
        if n != 0:
            raise IndexError("generator n not defined")
        return self.__generator

    def gens_dict(self):
        """
        Return a dictionary whose entries are ``{name:variable,...}``,
        where ``name`` stands for the variable names of this
        object (as strings) and ``variable`` stands for the corresponding
        generators (as elements of this object).

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.gens_dict()
            {'x': x}
        """
        return dict(zip(self.variable_names(), self.gens()))

    def parameter(self):
        """
        Return the generator of this skew polynomial ring.

        This is the same as ``self.gen()``.
        """
        return self.gen()

    def is_finite(self):
        """
        Return ``False`` since skew polynomial rings are not finite (unless the
        base ring is 0.)

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: k.is_finite()
            True
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_finite()
            False
        """
        R = self.base_ring()
        if R.is_finite() and R.order() == 1:
            return True
        return False

    def is_exact(self):
        """
        Return ``True`` if elements of this skew polynomial ring are exact.
        It happens if and only if elements of the base ring are exact.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_exact()
            True
            sage: S.base_ring().is_exact()
            True

            sage: R.<u> = k[[]]
            sage: sigma = R.hom([u+u^2])
            sage: T.<y> = R['y',sigma]
            sage: T.is_exact()
            False
            sage: T.base_ring().is_exact()
            False
        """
        return self.base_ring().is_exact()

    def is_sparse(self):
        """
        Return ``True`` if elements of this polynomial ring have a sparse
        representation.

        Since sparse skew polynomials are not yet implemented, this
        function always returns ``False``.

        EXAMPLES:

            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.is_sparse()
            False
        """
        return self.__is_sparse

    def ngens(self):
        """
        Return the number of generators of this skew polynomial ring, which is 1.

        EXAMPLES::

            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: S.ngens()
            1
        """
        return 1

    def random_element(self, degree=2, monic=False, *args, **kwds):
        r"""
        Return a random skew polynomial.

        INPUT:

        -  ``degree`` - Integer with degree (default: 2)
           or a tuple of integers with minimum and maximum degrees

        -  ``monic`` - if True, returns a monic skew polynomial
           (default: ``False``)

        -  ``*args, **kwds`` - Passed on to the ``random_element`` method for
           the base ring

        OUTPUT:

        -  Skew polynomial such that the coefficients of `x^i`, for `i` up
           to ``degree``, are random elements from the base ring, randomized
           subject to the arguments ``*args`` and ``**kwds``

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.random_element()  # random
            (2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2
            sage: S.random_element(monic=True)  # random
            x^2 + (2*t^2 + t + 1)*x + 3*t^2 + 3*t + 2

        If a tuple of two integers is given for the degree argument, a random
        integer will be chosen between the first and second element of the
        tuple as the degree::

            sage: S.random_element(degree=(2,7))  # random
            (3*t^2 + 1)*x^4 + (4*t + 2)*x^3 + (4*t + 1)*x^2 + (t^2 + 3*t + 3)*x + 3*t^2 + 2*t + 2

        If the minimal degree is greater than the maximal degree, sage raises
        a ValueError::

            sage: S.random_element(degree=(5,4))
            Traceback (most recent call last):
            ...
            ValueError: minimum degree must be less or equal than maximum degree

        When ``monic`` is false, the returned skew polynomial may have a degree
        less than ``degree`` (it happens when the random ``leading coefficient``
        is zero)::

            sage: S.random_element(degree=4) #random
            (3*t^2 + t)*x^3 + (2*t + 2)*x^2 + (3*t^2 + 2*t + 2)*x + t

        However, if ``monic`` is true, this can't happen::

            sage: S.random_element(degree=4,monic=True)  # random
            x^4 + (t^2 + 3*t + 3)*x^3 + (t^2 + 4*t + 3)*x^2 + (2*t^2 + 4*t + 4)*x + 2*t^2 + 2*t
        """
        R = self.base_ring()
        if isinstance(degree, (list, tuple)):
            if len(degree) != 2:
                raise ValueError("degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)")
            if degree[0] > degree[1]:
                raise ValueError("minimum degree must be less or equal than maximum degree")
            degree = randint(*degree)
        if monic:
            return self ([ R.random_element (*args, **kwds) for _ in range (degree) ] + [ R.one() ])
        else:
            return self ([ R.random_element (*args, **kwds) for _ in range (degree+1) ])

    def is_commutative(self):
        """
        Return ``True`` if this skew polynomial ring is commutative
        (i.e. if the twist map is the identity).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.is_commutative()
            False

            sage: T.<y> = k['y',Frob^3]
            sage: T.is_commutative()
            True
        """
        return self.twist_map().is_identity()

    def minimal_vanishing_polynomial(self, eval_pts):
        """
        Return the minimal-degree, monic skew polynomial which vanishes at all
        the given evaluation points.

        The degree of the vanishing polynomial is at most the length of
        ``eval_pts``. Equality holds if and only if the elements of ``eval_pts``
        are linearly independent over the fixed field of ``self.twist_map()``.

        INPUT:

        - ``eval_pts`` -- list of evaluation points which are linearly
          independent over the fixed field of the twist map of the associated
          skew polynomial ring

        OUTPUT:

        The minimal vanishing polynomial.

        EXAMPLES:

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: eval_pts = [1, t, t^2]
            sage: b = S.minimal_vanishing_polynomial(eval_pts); b
            x^3 + 4

        The minimal vanishing polynomial evaluates to 0 at each of the evaluation points::

            sage: eval = b.multi_point_evaluation(eval_pts); eval
            [0, 0, 0]

        If the evaluation points are linearly dependent over the fixed field of
        the twist map, then the returned polynomial has lower degree than the
        number of evaluation points:

            sage: S.minimal_vanishing_polynomial([t])
            x + 3*t^2 + 3*t
            sage: S.minimal_vanishing_polynomial([t, 3*t])
            x + 3*t^2 + 3*t
        """
        return _minimal_vanishing_polynomial(_base_ring_to_fraction_field(self), eval_pts)

    def lagrange_polynomial(self, points, check=True):
        """
        Return the minimal-degree polynomial which interpolates the given
        points, if one exists.

        More precisely, given `n` pairs `(x_1, y_1), ..., (x_n, y_n) \in R^2`,
        where `R` is ``self.base_ring()``, determine whether a skew polynomial
        `p(x)` such that `p(x_i) = y_i` for each `i` exists, and return it if
        it does.

        If the `x_i` are linearly independent over the fixed field of
        ``self.twist_map()`` then such a polynomial is guaranteed to exist.

        INPUT:

        - ``points`` -- a list of pairs ``(x_1, y_1),..., (x_n, y_n)`` of
          elements of the base ring of ``self``.

        # - ``check`` -- boolean (default: ``True``): if ``True`` then the computed polynomial is 
        # that verifies whether the
        #   `eval_pts` are linearly independent in the fixed field of twist map of
        #   the associated skew polynomial ring

        TODO Better doc of ``check`` etc.

        TODO Example showing linear dependence

        OUTPUT:

        The lagrange polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: points = [(t, 3*t^2 + 4*t + 4), (t^2, 4*t)]
            sage: d = S.lagrange_polynomial(points); d
            x + t

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: T.<x> = R['x', sigma]
            sage: points = [ (1, t^2 + 3*t + 4), (t, 2*t^2 + 3*t + 1), (t^2, t^2 + 3*t + 4) ]
            sage: p = T.lagrange_polynomial(points); p
            ((-t^4 - 2*t - 3)/-2)*x^2 + (-t^4 - t^3 - t^2 - 3*t - 2)*x + (-t^4 - 2*t^3 - 4*t^2 - 10*t - 9)/-2
            sage: p.multi_point_evaluation([1, t, t^2]) == [ t^2 + 3*t + 4, 2*t^2 + 3*t + 1, t^2 + 3*t + 4 ]
            True
        """
        l = len(points)
        if not all( len(pair) == 2 for pair in points ):
            raise TypeError("supplied points must be pairs of elements of base ring")
        eval_pts, values = zip(*points) #unzip
        eval_pts = list(eval_pts)
        values = list(values)
        
        if l > len(set(eval_pts)):
            raise TypeError("the evaluation points must be distinct")
        zero_i = [ i for i in range(l) if eval_pts[i].is_zero() ]
        if zero_i and not values[zero_i[0]].is_zero():
            raise TypeError("a skew polynomial always evaluates to 0 at 0, but a non-zero value was requested.")

        p = _lagrange_interpolation(_base_ring_to_fraction_field(self), eval_pts, values)

        if check and p.multi_point_evaluation(eval_pts) != values:
            raise ValueError("the requested Lagrange polynomial does not exist (the evaluation points are linearly dependent over the fixed field of the twist map and the sought values do not match.)")

        return p

class SkewPolynomialRing_finite_field(SkewPolynomialRing_general):
    """
    A specialized class for skew polynomial rings over finite fields.

    .. SEEALSO::

        :meth:`sage.rings.polynomial.skew_polynomial_ring_constructor.SkewPolynomialRing`
        :class:`sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general`
        :mod:`sage.rings.polynomial.skew_polynomial_finite_field`

    .. TODO::

        Add methods related to center of skew polynomial ring, irreducibility, karatsuba
        multiplication and factorization.
    """
    @staticmethod
    def __classcall__(cls, base_ring, map, name=None, sparse=False, element_class=None):
        if not element_class:
            if sparse:
                raise NotImplementedError("sparse skew polynomials are not implemented")
            else:
                from sage.rings.polynomial import skew_polynomial_finite_field
                element_class = skew_polynomial_finite_field.SkewPolynomial_finite_field_dense
                return super(SkewPolynomialRing_general,cls).__classcall__(cls,base_ring,map,name,sparse,element_class)

    def __init__(self, base_ring, map, name, sparse, element_class):
        """
        This method is a constructor for a general, dense univariate skew polynomial ring
        over a finite field.

        INPUT::

        - ``base_ring`` -- a commutative ring

        - ``map`` -- an automorphism of the base ring

        - ``name`` -- string or list of strings representing the name of the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``element_class`` -- class representing the type of element to be used in ring

        ..NOTE::

            Multivariate and Sparse rings are not implemented.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<x> = k['x', Frob]; T
            Skew Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        self._order = -1
        try:
            self._order = map.order()
        except (AttributeError,NotImplementedError):
            pass
        if self._order < 0:
            try:
                if map.is_identity():
                    self._order = 1
            except (AttributeError,NotImplementedError):
                pass
        if self._order < 0:
            raise NotImplementedError("unable to determine the order of %s" % map)
        SkewPolynomialRing_general.__init__ (self, base_ring, map, name, sparse, element_class)
        self._maps = [ map**i for i in range(self._order) ]

    def twist_map(self, n=1):
        """
        Return the twist map, otherwise known as the automorphism over the base ring of
        ``self``, iterated `n` times.

        INPUT:

        -  ``n`` - a relative integer (default: 1)

        OUTPUT:

        -  The `n`-th iterative of the twist map of this skew polynomial ring.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: S.twist_map()
            Frobenius endomorphism t |--> t^5 on Finite Field in t of size 5^3
            sage: S.twist_map(11)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3
            sage: S.twist_map(3)
            Identity endomorphism of Finite Field in t of size 5^3

        It also works if `n` is negative::

            sage: S.twist_map(-1)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3
        """
        return self._maps[n%self._order]
