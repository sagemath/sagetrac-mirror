r"""
Ring of Puiseux Polynomials

If `R` is a commutative ring, then the ring of Puiseux polynomials in `n`
variables over `R` is `R[x_1^{\pm 1}, x_2^{\pm 1}, \ldots, x_n^{\pm 1}]`.
We implement it as a quotient ring

.. MATH::

    R[x_1, y_1, x_2, y_2, \ldots, x_n, y_n] / (x_1 y_1 - 1, x_2 y_2 - 1, \ldots, x_n y_n - 1).

TESTS::

    sage: P.<q> = PuiseuxPolynomialRing(QQ)
    sage: qi = q^(1/4)
    sage: qi in P
    True
    sage: P(qi)
    q^{1/4}

    sage: A.<Y> = QQ[]
    sage: R.<X> = PuiseuxPolynomialRing(A)
    sage: matrix(R,2,2,[X,0,0,1])
    [X 0]
    [0 1]

"""
#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed@math.harvard.edu>,
#                          William Stein <wstein@gmail.com>,
#                          Mike Hansen <mhansen@gmail.com>
#                          Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import
from six import iteritems, iterkeys
from six.moves import range

from sage.structure.element import parent
from sage.structure.parent import Parent
from sage.rings.infinity import infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.latex import latex
from sage.rings.polynomial.puiseux_polynomial import PuiseuxPolynomial_univariate
from sage.rings.ring import CommutativeRing

def is_PuiseuxPolynomialRing(R):
    """
    Returns True if and only if R is a Puiseux polynomial ring.

    EXAMPLES::

        sage: from sage.rings.polynomial.laurent_polynomial_ring import is_PuiseuxPolynomialRing
        sage: P = PolynomialRing(QQ,2,'x')
        sage: is_PuiseuxPolynomialRing(P)
        False

        sage: R = PuiseuxPolynomialRing(QQ,3,'x')
        sage: is_PuiseuxPolynomialRing(R)
        True
    """
    return isinstance(R, PuiseuxPolynomialRing_generic)


def PuiseuxPolynomialRing(base_ring, *args, **kwds):
    r"""
    Return the globally unique univariate or multivariate Puiseux polynomial
    ring with given properties and variable name or names.

    There are four ways to call the Puiseux polynomial ring constructor:

    1. ``PuiseuxPolynomialRing(base_ring, name,    sparse=False)``
    2. ``PuiseuxPolynomialRing(base_ring, names,   order='degrevlex')``
    3. ``PuiseuxPolynomialRing(base_ring, name, n, order='degrevlex')``
    4. ``PuiseuxPolynomialRing(base_ring, n, name, order='degrevlex')``

    The optional arguments sparse and order *must* be explicitly
    named, and the other arguments must be given positionally.

    INPUT:

    - ``base_ring`` -- a commutative ring
    - ``name`` -- a string
    - ``names`` -- a list or tuple of names, or a comma separated string
    - ``n`` -- a positive integer
    - ``sparse`` -- bool (default: False), whether or not elements are sparse
    - ``order`` -- string or
      :class:`~sage.rings.polynomial.term_order.TermOrder`, e.g.,

        - ``'degrevlex'`` (default) -- degree reverse lexicographic
        - ``'lex'`` -- lexicographic
        - ``'deglex'`` -- degree lexicographic
        - ``TermOrder('deglex',3) + TermOrder('deglex',3)`` -- block ordering

    OUTPUT:

    ``PuiseuxPolynomialRing(base_ring, name, sparse=False)`` returns a
    univariate Puiseux polynomial ring; all other input formats return a
    multivariate Puiseux polynomial ring.

    UNIQUENESS and IMMUTABILITY: In Sage there is exactly one
    single-variate Puiseux polynomial ring over each base ring in each choice
    of variable and sparseness.  There is also exactly one multivariate
    Puiseux polynomial ring over each base ring for each choice of names of
    variables and term order.

    ::

        sage: R.<x,y> = PuiseuxPolynomialRing(QQ,2); R
        Multivariate Puiseux Polynomial Ring in x, y over Rational Field
        sage: f = x^2 - 2*y^-2

    You can't just globally change the names of those variables.
    This is because objects all over Sage could have pointers to
    that polynomial ring.

    ::

        sage: R._assign_names(['z','w'])
        Traceback (most recent call last):
        ...
        ValueError: variable names cannot be changed after object creation.


    EXAMPLES:

    1. ``PuiseuxPolynomialRing(base_ring, name, sparse=False)``

       ::

           sage: PuiseuxPolynomialRing(QQ, 'w')
           Univariate Puiseux Polynomial Ring in w over Rational Field

       Use the diamond brackets notation to make the variable
       ready for use after you define the ring::

           sage: R.<w> = PuiseuxPolynomialRing(QQ)
           sage: (1 + w)^3
           1 + 3*w + 3*w^2 + w^3

       You must specify a name::

           sage: PuiseuxPolynomialRing(QQ)
           Traceback (most recent call last):
           ...
           TypeError: you must specify the names of the variables

           sage: R.<abc> = PuiseuxPolynomialRing(QQ, sparse=True); R
           Univariate Puiseux Polynomial Ring in abc over Rational Field

           sage: R.<w> = PuiseuxPolynomialRing(PolynomialRing(GF(7),'k')); R
           Univariate Puiseux Polynomial Ring in w over Univariate Polynomial Ring in k over Finite Field of size 7

       Rings with different variables are different::

           sage: PuiseuxPolynomialRing(QQ, 'x') == PuiseuxPolynomialRing(QQ, 'y')
           False

    2. ``PuiseuxPolynomialRing(base_ring, names,   order='degrevlex')``

       ::

           sage: R = PuiseuxPolynomialRing(QQ, 'a,b,c'); R
           Multivariate Puiseux Polynomial Ring in a, b, c over Rational Field

           sage: S = PuiseuxPolynomialRing(QQ, ['a','b','c']); S
           Multivariate Puiseux Polynomial Ring in a, b, c over Rational Field

           sage: T = PuiseuxPolynomialRing(QQ, ('a','b','c')); T
           Multivariate Puiseux Polynomial Ring in a, b, c over Rational Field

       All three rings are identical.

       ::

           sage: (R is S) and  (S is T)
           True

       There is a unique Puiseux polynomial ring with each term order::

           sage: R = PuiseuxPolynomialRing(QQ, 'x,y,z', order='degrevlex'); R
           Multivariate Puiseux Polynomial Ring in x, y, z over Rational Field
           sage: S = PuiseuxPolynomialRing(QQ, 'x,y,z', order='invlex'); S
           Multivariate Puiseux Polynomial Ring in x, y, z over Rational Field
           sage: S is PuiseuxPolynomialRing(QQ, 'x,y,z', order='invlex')
           True
           sage: R == S
           False


    3. ``PuiseuxPolynomialRing(base_ring, name, n, order='degrevlex')``

       If you specify a single name as a string and a number of
       variables, then variables labeled with numbers are created.

       ::

           sage: PuiseuxPolynomialRing(QQ, 'x', 10)
           Multivariate Puiseux Polynomial Ring in x0, x1, x2, x3, x4, x5, x6, x7, x8, x9 over Rational Field

           sage: PuiseuxPolynomialRing(GF(7), 'y', 5)
           Multivariate Puiseux Polynomial Ring in y0, y1, y2, y3, y4 over Finite Field of size 7

           sage: PuiseuxPolynomialRing(QQ, 'y', 3, sparse=True)
           Multivariate Puiseux Polynomial Ring in y0, y1, y2 over Rational Field

       By calling the
       :meth:`~sage.structure.category_object.CategoryObject.inject_variables`
       method, all those variable names are available for interactive use::

           sage: R = PuiseuxPolynomialRing(GF(7),15,'w'); R
           Multivariate Puiseux Polynomial Ring in w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14 over Finite Field of size 7
           sage: R.inject_variables()
           Defining w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14
           sage: (w0 + 2*w8 + w13)^2
           w0^2 + 4*w0*w8 + 4*w8^2 + 2*w0*w13 + 4*w8*w13 + w13^2
    """
    from sage.rings.polynomial.polynomial_ring import is_PolynomialRing

    R = PolynomialRing(base_ring, *args, **kwds)

    if is_PolynomialRing(R):
        # univariate case
        P = PuiseuxPolynomialRing_univariate(R)
    else:
        raise NotImplementedError

    return P


class PuiseuxPolynomialRing_generic(CommutativeRing, Parent):
    """
    Puiseux polynomial ring (base class).

    EXAMPLES:

    This base class inherits from :class:`~sage.rings.ring.CommutativeRing`.
    Since :trac:`11900`, it is also initialised as such::

        sage: R.<x1,x2> = PuiseuxPolynomialRing(QQ)
        sage: R.category()
        Category of commutative rings
        sage: TestSuite(R).run()

    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: R = PuiseuxPolynomialRing(QQ,2,'x')
            sage: R == loads(dumps(R))
            True
        """
        self._n = R.ngens()
        self._R = R
        names = R.variable_names()
        CommutativeRing.__init__(self, R.base_ring(), names=names)
        self._populate_coercion_lists_(init_no_parent=True)

    def ngens(self):
        """
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').ngens()
            2
            sage: PuiseuxPolynomialRing(QQ,1,'x').ngens()
            1
        """
        return self._n

    def gen(self, i=0):
        r"""
        Returns the `i^{th}` generator of self.  If i is not specified, then
        the first generator will be returned.

        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').gen()
            x0
            sage: PuiseuxPolynomialRing(QQ,2,'x').gen(0)
            x0
            sage: PuiseuxPolynomialRing(QQ,2,'x').gen(1)
            x1

        TESTS::

            sage: PuiseuxPolynomialRing(QQ,2,'x').gen(3)
            Traceback (most recent call last):
            ...
            ValueError: generator not defined
        """
        if i < 0 or i >= self._n:
            raise ValueError("generator not defined")
        try:
            return self.__generators[i]
        except AttributeError:
            self.__generators = tuple(self(x) for x in self._R.gens())
            return self.__generators[i]

    def variable_names_recursive(self, depth=infinity):
        r"""
        Return the list of variable names of this ring and its base rings,
        as if it were a single multi-variate Puiseux polynomial.

        INPUT:

        - ``depth`` -- an integer or :mod:`Infinity <sage.rings.infinity>`.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: T = PuiseuxPolynomialRing(QQ, 'x')
            sage: S = PuiseuxPolynomialRing(T, 'y')
            sage: R = PuiseuxPolynomialRing(S, 'z')
            sage: R.variable_names_recursive()
            ('x', 'y', 'z')
            sage: R.variable_names_recursive(2)
            ('y', 'z')
        """
        if depth <= 0:
            return ()
        elif depth == 1:
            return self.variable_names()
        else:
            my_vars = self.variable_names()
            try:
               return self.base_ring().variable_names_recursive(depth - len(my_vars)) + my_vars
            except AttributeError:
                return my_vars

    def is_integral_domain(self, proof = True):
        """
        Returns True if self is an integral domain.

        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').is_integral_domain()
            True

        The following used to fail; see :trac:`7530`::

            sage: L = PuiseuxPolynomialRing(ZZ, 'X')
            sage: L['Y']
            Univariate Polynomial Ring in Y over Univariate Puiseux Polynomial Ring in X over Integer Ring
        """
        return self.base_ring().is_integral_domain(proof)

    def is_noetherian(self):
        """
        Returns True if self is Noetherian.

        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').is_noetherian()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def remove_var(self, var):
        """
        EXAMPLES::

            sage: R = PuiseuxPolynomialRing(QQ,'x,y,z')
            sage: R.remove_var('x')
            Multivariate Puiseux Polynomial Ring in y, z over Rational Field
            sage: R.remove_var('x').remove_var('y')
            Univariate Puiseux Polynomial Ring in z over Rational Field
        """
        vars = list(self.variable_names())
        vars.remove(str(var))
        return PuiseuxPolynomialRing(self.base_ring(), vars)

    def _coerce_map_from_(self, R):
        """
        EXAMPLES::

            sage: L.<x,y> = PuiseuxPolynomialRing(QQ)
            sage: L.coerce_map_from(QQ)
            Composite map:
              From: Rational Field
              To:   Multivariate Puiseux Polynomial Ring in x, y over Rational Field
              Defn:   Polynomial base injection morphism:
                      From: Rational Field
                      To:   Multivariate Polynomial Ring in x, y over Rational Field
                    then
                      Call morphism:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Multivariate Puiseux Polynomial Ring in x, y over Rational Field

        Let us check that coercion between Puiseux Polynomials over
        different base rings works (:trac:`15345`)::

            sage: R = PuiseuxPolynomialRing(ZZ, 'x')
            sage: T = PuiseuxPolynomialRing(QQ, 'x')
            sage: R.gen() + 3*T.gen()
            4*x
        """
        if R is self._R or (isinstance(R, PuiseuxPolynomialRing_generic)
            and self._R.has_coerce_map_from(R._R)):
            from sage.structure.coerce_maps import CallableConvertMap
            return CallableConvertMap(R, self, self._element_constructor_,
                                      parent_as_first_arg=False)
        elif isinstance(R, PuiseuxPolynomialRing_generic) and \
             R.variable_names() == self.variable_names() and \
             self.base_ring().has_coerce_map_from(R.base_ring()):
            return True

        f = self._R.coerce_map_from(R)
        if f is not None:
            from sage.categories.homset import Hom
            from sage.categories.morphism import CallMorphism
            return CallMorphism(Hom(self._R, self)) * f

    def __eq__(self, right):
        """
        Check whether ``self`` is equal to ``right``.

        EXAMPLES::

            sage: R = PuiseuxPolynomialRing(QQ,'x,y,z')
            sage: P = PuiseuxPolynomialRing(ZZ,'x,y,z')
            sage: Q = PuiseuxPolynomialRing(QQ,'x,y')

            sage: R == R
            True
            sage: R == Q
            False
            sage: Q == P
            False
            sage: P == R
            False
        """
        if type(self) != type(right):
            return False
        return self._R == right._R

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: R = PuiseuxPolynomialRing(QQ,'x,y,z')
            sage: P = PuiseuxPolynomialRing(ZZ,'x,y,z')
            sage: Q = PuiseuxPolynomialRing(QQ,'x,y')

            sage: R != R
            False
            sage: R != Q
            True
            sage: Q != P
            True
            sage: P != R
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: h1 = hash(PuiseuxPolynomialRing(ZZ,'x,y,z'))
            sage: h2 = hash(PuiseuxPolynomialRing(ZZ,'x,y,z'))
            sage: h3 = hash(PuiseuxPolynomialRing(QQ,'x,y,z'))
            sage: h4 = hash(PuiseuxPolynomialRing(ZZ,'x,y'))
            sage: h1 == h2 and h1 != h3 and h1 != h4
            True
        """
        return hash(self._R) ^ 12059065606945654693

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(PuiseuxPolynomialRing(QQ,2,'x'))
            \Bold{Q}[x_{0}^{\QQ}, x_{1}^{\QQ}]
        """
        vars = ', '.join([a + r'^{\QQ}' for a in self.latex_variable_names()])
        return "%s[%s]"%(latex(self.base_ring()), vars)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        EXAMPLES::

            sage: L.<x,y> = PuiseuxPolynomialRing(QQ)
            sage: L._is_valid_homomorphism_(QQ, (1/2, 3/2))
            True
        """
        if not codomain.has_coerce_map_from(self.base_ring()):
            # we need that elements of the base ring
            # canonically coerce into codomain.
            return False
        return True

    def term_order(self):
        """
        Returns the term order of self.

        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').term_order()
            Degree reverse lexicographic term order
        """
        return self._R.term_order()

    def is_finite(self):
        """
        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').is_finite()
            False

        """
        return False

    def is_field(self, proof = True):
        """
        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').is_field()
            False
        """
        return False

    def polynomial_ring(self):
        """
        Return the polynomial ring associated with self.

        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').polynomial_ring()
            Multivariate Polynomial Ring in x0, x1 over Rational Field
            sage: PuiseuxPolynomialRing(QQ,1,'x').polynomial_ring()
            Multivariate Polynomial Ring in x over Rational Field
        """
        return self._R

    def characteristic(self):
        """
        Returns the characteristic of the base ring.

        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').characteristic()
            0
            sage: PuiseuxPolynomialRing(GF(3),2,'x').characteristic()
            3
        """
        return self.base_ring().characteristic()

    def krull_dimension(self):
        """
        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').krull_dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def random_element(self, low_degree = -2, high_degree = 2, terms = 5, choose_degree=False,*args, **kwds):
        """
        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').random_element()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_exact(self):
        """
        Returns True if the base ring is exact.

        EXAMPLES::

            sage: PuiseuxPolynomialRing(QQ,2,'x').is_exact()
            True
            sage: PuiseuxPolynomialRing(RDF,2,'x').is_exact()
            False
        """
        return self.base_ring().is_exact()

    def change_ring(self, base_ring=None, names=None, sparse=False, order=None):
        """
        EXAMPLES::

            sage: R = PuiseuxPolynomialRing(QQ,2,'x')
            sage: R.change_ring(ZZ)
            Multivariate Puiseux Polynomial Ring in x0, x1 over Integer Ring
        """
        if base_ring is None:
            base_ring = self.base_ring()
        if names is None:
            names = self.variable_names()
        if self._n == 1:
            return PuiseuxPolynomialRing(base_ring, names[0], sparse = sparse)

        if order is None:
            order = self.polynomial_ring().term_order()
        return PuiseuxPolynomialRing(base_ring, self._n, names, order = order)


class PuiseuxPolynomialRing_univariate(PuiseuxPolynomialRing_generic):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: L = PuiseuxPolynomialRing(QQ,'x')
            sage: type(L)
            <class 'sage.rings.polynomial.laurent_polynomial_ring.PuiseuxPolynomialRing_univariate_with_category'>
            sage: L == loads(dumps(L))
            True


        TESTS::

            sage: TestSuite(PuiseuxPolynomialRing(Zmod(4), 'y')).run()
            sage: TestSuite(PuiseuxPolynomialRing(ZZ, 'u')).run()
            sage: TestSuite(PuiseuxPolynomialRing(Zmod(4)['T'], 'u')).run()
        """
        if R.ngens() != 1:
            raise ValueError("must be 1 generator")
        PuiseuxPolynomialRing_generic.__init__(self, R)

    Element = PuiseuxPolynomial_univariate

    def _repr_(self):
        """
        TESTS::

            sage: PuiseuxPolynomialRing(QQ,'x')  # indirect doctest
            Univariate Puiseux Polynomial Ring in x over Rational Field
        """
        return "Univariate Puiseux Polynomial Ring in %s over %s" % (self._R.variable_name(), self._R.base_ring())

    def _element_constructor_(self, x):
        """
        EXAMPLES::

            sage: L = PuiseuxPolynomialRing(QQ, 'x')
            sage: L(1/2)
            1/2

            sage: L(x + 3/x)
            3*x^-1 + x

        ::

            sage: L(exp(x))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert e^x to a rational

        ::

            sage: U = PuiseuxPolynomialRing(QQ, 'a')
            sage: V = PuiseuxPolynomialRing(QQ, 'c')
            sage: L.<a, b, c, d> = PuiseuxPolynomialRing(QQ)
            sage: M = PuiseuxPolynomialRing(QQ, 'c, d')
            sage: Mc, Md = M.gens()
            sage: N = PuiseuxPolynomialRing(M, 'a, b')
            sage: Na, Nb = N.gens()
            sage: U(Na)
            a
            sage: V(Mc)
            c

            sage: M(L(0))
            0
            sage: N(L(0))
            0
            sage: L(M(0))
            0
            sage: L(N(0))
            0

        ::

            sage: A.<a> = PuiseuxPolynomialRing(QQ)
            sage: B.<b> = PuiseuxPolynomialRing(A)
            sage: B(a)
            a
            sage: C.<c> = PuiseuxPolynomialRing(B)
            sage: B(C(b))
            b
            sage: D.<d, e> = PuiseuxPolynomialRing(B)
            sage: B(D(b))
            b
        """
        from sage.symbolic.expression import Expression
        if isinstance(x, Expression):
            raise NotImplementedError
            # return x.puiseux_polynomial(ring=self)

        elif isinstance(x, (PuiseuxPolynomial_univariate)):
            P = x.parent()
            if len(self.variable_names()) == len(P.variable_names()):
                x = x.dict()

        return self.element_class(self, x)

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: L = PuiseuxPolynomialRing(QQ, 'x')
            sage: loads(dumps(L)) == L
            True
        """
        return PuiseuxPolynomialRing_univariate, (self._R,)


