r"""
`q`-Commuting Polynomials

AUTHORS:

- Travis Scrimshaw (2022-08-23): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 Travis Scrimshaw  <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.categories.algebras import Algebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.monoids.free_abelian_monoid import FreeAbelianMonoid
from sage.matrix.constructor import matrix
from sage.structure.element import Matrix

class qCommutingPolynomials(CombinatorialFreeModule):
    r"""
    The algebra of `q`-commuting polynomials.

    Let `R` be a commutative ring, and fix an element `q \in R`. Let
    B = (B_{xy})_{x,y \in I}`  be a skew-symmetric bilinear form with
    index set `I`. Let `R[I]_{q,B}` denote the polynomial ring in the variables
    `I` such that we have the `q`-*commuting* relation for `x, y \in I`:

    .. MATH::

        y x = q^{B_{xy}} \cdot x y.

    This is a graded `R`-algebra with a natural basis given by monomials
    written in increasing order with respect to some total order on `I`.

    When `B_{xy} = 1` and `B_{yx} = -1` for all `x < y`, then we have
    a `q`-analog of the classical binomial coefficient theorem:

    .. MATH::

        (x + y)^n = \sum_{k=0}^n \binom{n}{k}_q x^k y^{n-k}.

    EXAMPLES::

        sage: q = ZZ['q'].fraction_field().gen()
        sage: R.<x,y> = algebras.qCommutingPolynomials(q)

    We verify a case of the `q`-binomial theorem::

        sage: f = (x + y)^10
        sage: all(f[b] == q_binomial(10, b.list()[0]) for b in f.support())
        True

    We now do a computation with a non-standard `B` matrix::

        sage: B = matrix([[0,1,2],[-1,0,3],[-2,-3,0]])
        sage: B
        [ 0  1  2]
        [-1  0  3]
        [-2 -3  0]
        sage: q = ZZ['q'].gen()
        sage: R.<x,y,z> = algebras.qCommutingPolynomials(q, B)
        sage: y * x
        q*x*y
        sage: z * x
        q^2*x*z
        sage: z * y
        q^3*y*z

        sage: f = (x + z)^10
        sage: all(f[b] == q_binomial(10, b.list()[0], q^2) for b in f.support())
        True

        sage: f = (y + z)^10
        sage: all(f[b] == q_binomial(10, b.list()[1], q^3) for b in f.support())
        True
    """
    @staticmethod
    def __classcall_private__(cls, q, n=None, B=None, base_ring=None, names=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R1.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R2 = algebras.qCommutingPolynomials(q, base_ring=q.parent(), names='x,y,z')
            sage: R3 = algebras.qCommutingPolynomials(q, names=['x', 'y', 'z'])
            sage: R1 is R2 is R3
            True
        """
        if base_ring is not None:
            q = base_ring(q)

        if B is None and isinstance(n, Matrix):
            n, B = B, n

        if names is None:
            raise ValueError("the names of the variables must be given")
        from sage.structure.category_object import normalize_names
        if n is None:
            if isinstance(names, str):
                n = names.count(',') + 1
            else:
                n = len(names)
        names = normalize_names(n, names)
        n = len(names)
        if B is None:
            B = matrix.zero(ZZ, n)
            for i in range(n):
                for j in range(i+1, n):
                    B[i,j] = 1
                    B[j,i] = -1
            B.set_immutable()
        else:
            if not B.is_skew_symmetric():
                raise ValueError("the matrix must be skew symmetric")
            B = B.change_ring(ZZ)
            B.set_immutable()
        return super().__classcall__(cls, q=q, B=B, names=names)

    def __init__(self, q, B, names):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: TestSuite(R).run()
        """
        self._q = q
        self._B = B
        base_ring = q.parent()
        indices = FreeAbelianMonoid(len(names), names)
        category = Algebras(base_ring).WithBasis().Graded()
        CombinatorialFreeModule.__init__(self, base_ring, indices,
                                         bracket=False, prefix='',
                                         sorting_key=qCommutingPolynomials._term_key,
                                         names=indices.variable_names(), category=category)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R
            q-commuting polynomial ring in x, y, z over Fraction Field of
             Univariate Polynomial Ring in q over Integer Ring with matrix:
            [ 0  1  1]
            [-1  0  1]
            [-1 -1  0]
        """
        names = ", ".join(self.variable_names())
        return "{}-commuting polynomial ring in {} over {} with matrix:\n{}".format(self._q, names, self.base_ring(), self._B)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: latex(R)
            \mathrm{Frac}(\Bold{Z}[q])[x, y, z]_{q}
        """
        from sage.misc.latex import latex
        names = ", ".join(self.variable_names())
        return "{}[{}]_{{{}}}".format(latex(self.base_ring()), names, self._q)

    @staticmethod
    def _term_key(x):
        r"""
        Compute a key for ``x`` for comparisons.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: elt = (x*y^3*z^2).leading_support()
            sage: R._term_key(elt)
            (6, [2, 3, 1])
        """
        L = x.list()
        L.reverse()
        return (sum(L), L)

    def gen(self, i):
        r"""
        Return the ``i``-generator of ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.gen(0)
            x
            sage: R.gen(2)
            z
        """
        return self.monomial(self._indices.gen(i))

    @cached_method
    def gens(self):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.gens()
            (x, y, z)
        """
        return tuple([self.monomial(g) for g in self._indices.gens()])

    @cached_method
    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.algebra_generators()
            Finite family {'x': x, 'y': y, 'z': z}
        """
        d = {v: self.gen(i) for i,v in enumerate(self.variable_names())}
        return Family(self.variable_names(), d.__getitem__, name="generator")

    @cached_method
    def one_basis(self):
        r"""
        Return the basis index of the element `1`.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.one_basis()
            1
        """
        return self._indices.one()

    def degree_on_basis(self, m):
        r"""
        Return the degree of the monomial index by ``m``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.degree_on_basis(R.one_basis())
            0
            sage: f = (x + y)^3 + z^3
            sage: f.degree()
            3
        """
        return sum(m.list())

    def dimension(self):
        r"""
        Return the dimension of ``self``, which is `\infty`.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.dimension()
            +Infinity
        """
        return infinity

    @cached_method
    def product_on_basis(self, x, y):
        r"""
        Return the product of two monomials given by ``x`` and ``y``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y> = algebras.qCommutingPolynomials(q)
            sage: R.product_on_basis(x.leading_support(), y.leading_support())
            x*y
            sage: R.product_on_basis(y.leading_support(), x.leading_support())
            q*x*y

            sage: x * y
            x*y
            sage: y * x
            q*x*y
            sage: y^2 * x
            q^2*x*y^2
            sage: y * x^2
            q^2*x^2*y
            sage: x * y * x
            q*x^2*y
            sage: y^2 * x^2
            q^4*x^2*y^2
            sage: (x + y)^2
            x^2 + (q+1)*x*y + y^2
            sage: (x + y)^3
            x^3 + (q^2+q+1)*x^2*y + (q^2+q+1)*x*y^2 + y^3
            sage: (x + y)^4
            x^4 + (q^3+q^2+q+1)*x^3*y + (q^4+q^3+2*q^2+q+1)*x^2*y^2 + (q^3+q^2+q+1)*x*y^3 + y^4

        With a non-standard `B` matrix::

            sage: B = matrix([[0,1,2],[-1,0,3],[-2,-3,0]])
            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q, B=B)
            sage: x * y
            x*y
            sage: y * x^2
            q^2*x^2*y
            sage: z^2 * x
            q^4*x*z^2
            sage: z^2 * x^3
            q^12*x^3*z^2
            sage: z^2 * y
            q^6*y*z^2
            sage: z^2 * y^3
            q^18*y^3*z^2
        """
        # Special case for multiplying by 1
        if x == self.one_basis():
            return self.monomial(y)
        if y == self.one_basis():
            return self.monomial(x)

        Lx = x.list()
        Ly = y.list()

        # This could be made more efficient
        qpow = sum(exp * sum(self._B[j,i] * val for j, val in enumerate(Ly[:i])) for i,exp in enumerate(Lx))
        return self.term(x * y, self._q ** qpow)

