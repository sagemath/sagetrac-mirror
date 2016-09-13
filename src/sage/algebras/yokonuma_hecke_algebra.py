"""
Yokonuma-Hecke Algebras

AUTHORS:

- Travis Scrimshaw (2015-11): initial version
"""

#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.all import QQ
from sage.categories.algebras import Algebras
from sage.categories.rings import Rings
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutations
from sage.sets.family import Family

class YokonumaHeckeAlgebra(CombinatorialFreeModule):
    r"""
    The Yokonuma-Hecke algebra `Y_{d,n}(q)`.

    Let `R` be a commutative ring and `q` be a unit in `R`. The
    *Yokonuma-Hecke algebra* `Y_{d,n}(q)` is the associative, unital
    `R`-algebra generated by `t_1, t_2, \ldots, t_n, g_1, g_2, \ldots,
    g_{n-1}` and subject to the relations:

    - `g_i g_j = g_j g_i` for all `|i - j| > 1`,
    - `g_i g_{i+1} g_i = g_{i+1} g_i g_{i+1}`,
    - `t_i t_j = t_j t_i`,
    - `t_j g_i = g_i t_{j s_i}`, and
    - `t_j^d = 1`,

    where `s_i` is the simple transposition `(i, i+1)`, along with
    the quadratic relation

    .. MATH::

        g_i^2 = 1 + \frac{(q - q^{-1})}{d} \left( \sum_{s=0}^{d-1}
        t_i^s t_{i+1}^{-s} \right) g_i.

    Thus the Yokonuma-Hecke algebra can be considered a quotient of
    the framed braid group `(\ZZ / d\ZZ) \wr B_n`, where `B_n` is the
    classical braid group on `n` strands, by the quadratic relations.
    Moreover, all of the algebra generators are invertible. In
    particular, we have

    .. MATH::

        g_i^{-1} = g_i - (q - q^{-1}) e_i.

    When we specialize `q = \pm 1`, we obtain the group algebra of
    the complex reflection group `G(d, 1, n) = (\ZZ / d\ZZ) \wr S_n`.
    Moreover for `d = 1`, the Yokonuma-Hecke algebra is equal to the
    :class`Iwahori-Hecke <IwahoriHeckeAlgebra>` of type `A_{n-1}`.

    INPUT:

    - ``d`` -- the maximum power of `t`
    - ``n`` -- the number of generators
    - ``q`` -- (optional) an invertible element in a commutative ring;
      the default is `q \in \QQ[q,q^{-1}]`
    - ``R`` -- (optional) a commutative ring containing ``q``; the
      default is the parent of `q`

    EXAMPLES:

    We construct `Y_{4,3}` and do some computations::

        sage: Y = algebras.YokonumaHecke(4, 3)
        sage: g1, g2, t1, t2, t3 = Y.algebra_generators()
        sage: g1 * g2
        g[1,2]
        sage: t1 * g1
        t1*g[1]
        sage: g2 * t2
        t3*g[2]
        sage: g2 * t3
        t2*g[2]
        sage: (g2 + t1) * (g1 + t2*t3)
        g[2,1] + t2*t3*g[2] + t1*g[1] + t1*t2*t3
        sage: g1 * g1
        1 - (1/4*q^-1-1/4*q)*g[1] - (1/4*q^-1-1/4*q)*t1*t2^3*g[1]
         - (1/4*q^-1-1/4*q)*t1^2*t2^2*g[1] - (1/4*q^-1-1/4*q)*t1^3*t2*g[1]
        sage: g2 * g1 * t1
        t3*g[2,1]

    We construct the elements `e_i` and show that they are idempotents::

        sage: e1 = Y.e(1); e1
        1/4 + 1/4*t1*t2^3 + 1/4*t1^2*t2^2 + 1/4*t1^3*t2
        sage: e1 * e1 == e1
        True
        sage: e2 = Y.e(2); e2
        1/4 + 1/4*t2*t3^3 + 1/4*t2^2*t3^2 + 1/4*t2^3*t3
        sage: e2 * e2 == e2
        True

    REFERENCES:

    .. [CL13] Maria Chlouveraki and Sofia Lambropoulou.
       *The Yokonuma-Hecke algebras and the HOMFLYPT polynomial*.
       (2015) :arxiv:`1204.1871v4`.

    .. [CPdA14] Maria Chlouveraki and Loic Poulain d'Andecy.
       *Representation theory of the Yokonuma-Hecke algebra*.
       (2014) :arxiv:`1302.6225v2`.

    .. [ERH15] Jorge Espanoza and Steen Ryom-Hansen.
       *Cell structures for the Yokonuma-Hecke algebra and the algebra
       of braids and ties*. (2015) :arxiv:`1506.00715`.

    .. [JPdA15] \N. Jacon and L. Poulain d'Andecy.
       *An isomorphism theorem for Yokonuma-Hecke algebras and
       applications to link invariants*. (2015) :arxiv:`1501.06389v3`.
    """
    @staticmethod
    def __classcall_private__(cls, d, n, q=None, R=None):
        """
        Standardize input to ensure a unique representation.

        TESTS::

            sage: Y1 = algebras.YokonumaHecke(5, 3)
            sage: q = LaurentPolynomialRing(QQ, 'q').gen()
            sage: Y2 = algebras.YokonumaHecke(5, 3, q)
            sage: Y3 = algebras.YokonumaHecke(5, 3, q, q.parent())
            sage: Y1 is Y2 and Y2 is Y3
            True
        """
        if q is None:
            q = LaurentPolynomialRing(QQ, 'q').gen()
        if R is None:
            R = q.parent()
        q = R(q)
        if R not in Rings().Commutative():
            raise TypeError("base ring must be a commutative ring")
        return super(YokonumaHeckeAlgebra, cls).__classcall__(cls, d, n, q, R)

    def __init__(self, d, n, q, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 3)
            sage: elts = Y.some_elements() + list(Y.algebra_generators())
            sage: TestSuite(Y).run(elements=elts)
        """
        self._d = d
        self._n = n
        self._q = q
        self._Pn = Permutations(n)
        import itertools
        C = itertools.product(*([range(d)]*n))
        indices = list( itertools.product(C, self._Pn))
        cat = Algebras(R).WithBasis()
        CombinatorialFreeModule.__init__(self, R, indices, prefix='Y',
                                         category=cat)
        self._assign_names(self.algebra_generators().keys())

    def _repr_(self):
        """ 
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.YokonumaHecke(5, 2)
            Yokonuma-Hecke algebra of rank 5 and order 2 with q=q
             over Univariate Laurent Polynomial Ring in q over Rational Field
        """
        return "Yokonuma-Hecke algebra of rank {} and order {} with q={} over {}".format(
            self._d, self._n, self._q, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 2)
            sage: latex(Y)
            \mathcal{Y}_{5,2}(q)
        """
        return "\\mathcal{Y}_{%s,%s}(%s)"%(self._d, self._n, self._q)

    def _repr_term(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: Y._repr_term( ((1, 0, 2), Permutation([3,2,1])) )
            't1*t3^2*g[2,1,2]'
        """
        gen_str = lambda e: '' if e == 1 else '^%s'%e
        lhs = '*'.join('t%s'%(j+1) + gen_str(i) for j,i in enumerate(m[0]) if i > 0)
        redword = m[1].reduced_word()
        if not redword:
            if not lhs:
                return '1'
            return lhs
        rhs = 'g[{}]'.format(','.join(str(i) for i in redword))
        if not lhs:
            return rhs
        return lhs + '*' + rhs

    def _latex_term(self, m):
        r"""
        Return a latex representation for the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: Y._latex_term( ((1, 0, 2), Permutation([3,2,1])) )
            't_{1} t_{3}^2 g_{2} g_{1} g_{2}'
        """
        gen_str = lambda e: '' if e == 1 else '^%s'%e
        lhs = ' '.join('t_{%s}'%(j+1) + gen_str(i) for j,i in enumerate(m[0]) if i > 0)
        redword = m[1].reduced_word()
        if not redword:
            if not lhs:
                return '1'
            return lhs
        return lhs + ' ' + ' '.join("g_{%d}"%i for i in redword)

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 3)
            sage: dict(Y.algebra_generators())
            {'g1': g[1], 'g2': g[2], 't1': t1, 't2': t2, 't3': t3}
        """
        one = self._Pn.one()
        zero = [0]*self._n
        d = {}
        for i in range(self._n):
            r = list(zero) # Make a copy
            r[i] = 1
            d['t%s'%(i+1)] = self.monomial( (tuple(r), one) )
        G = self._Pn.group_generators()
        for i in range(1, self._n):
            d['g%s'%i] = self.monomial( (tuple(zero), G[i]) )
        return Family(sorted(d), lambda i: d[i])

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 3)
            sage: Y.gens()
            (g[1], g[2], t1, t2, t3)
        """
        return tuple(self.algebra_generators())

    @cached_method
    def one_basis(self):
        """
        Return the index of the basis element of `1`.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(5, 3)
            sage: Y.one_basis()
            ((0, 0, 0), [1, 2, 3])
        """
        one = self._Pn.one()
        zero = [0]*self._n
        return (tuple(zero), one)

    @cached_method
    def e(self, i):
        """
        Return the element `e_i`.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: Y.e(1)
            1/4 + 1/4*t1*t2^3 + 1/4*t1^2*t2^2 + 1/4*t1^3*t2
            sage: Y.e(2)
            1/4 + 1/4*t2*t3^3 + 1/4*t2^2*t3^2 + 1/4*t2^3*t3
        """
        if i < 1 or i >= self._n:
            raise ValueError("invalid index")
        c = ~self.base_ring()(self._d)
        zero = [0]*self._n
        one = self._Pn.one()
        d = {}
        for s in range(self._d):
            r = list(zero) # Make a copy
            r[i-1] = s
            if s != 0:
                r[i] = self._d - s
            d[(tuple(r), one)] = c
        return self._from_dict(d, remove_zeros=False)

    def g(self, i=None):
        """
        Return the generator(s) `g_i`.

        INPUT:

        - ``i`` -- (default: ``None``) the generator `g_i` or if ``None``,
          then the list of all generators `g_i`

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(8, 3)
            sage: Y.g(1)
            g[1]
            sage: Y.g()
            [g[1], g[2]]
        """
        G = self.algebra_generators()
        if i is None:
            return [G['g%s'%i] for i in range(1, self._n)]
        return G['g%s'%i]

    def t(self, i=None):
        """
        Return the generator(s) `t_i`.

        INPUT:

        - ``i`` -- (default: ``None``) the generator `t_i` or if ``None``,
          then the list of all generators `t_i`

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(8, 3)
            sage: Y.t(2)
            t2
            sage: Y.t()
            [t1, t2, t3]
        """
        G = self.algebra_generators()
        if i is None:
            return [G['t%s'%i] for i in range(1, self._n+1)]
        return G['t%s'%i]

    def product_on_basis(self, m1, m2):
        """
        Return the product of the basis elements indexed by ``m1`` and ``m2``.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: m = ((1, 0, 2), Permutations(3)([2,1,3]))
            sage: 4 * Y.product_on_basis(m, m)
            -(q^-1-q)*t2^2*g[1] + 4*t1*t2 - (q^-1-q)*t1*t2*g[1]
             - (q^-1-q)*t1^2*g[1] - (q^-1-q)*t1^3*t2^3*g[1]

        Check that we apply the permutation correctly on `t_i`::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: g1, g2, t1, t2, t3 = Y.algebra_generators()
            sage: g21 = g2 * g1
            sage: g21 * t1
            t3*g[2,1]
        """
        t1,g1 = m1
        t2,g2 = m2
        # Commmute g1 and t2, then multiply t1 and t2
        #ig1 = g1
        t = [(t1[i] + t2[g1.index(i+1)]) % self._d for i in range(self._n)]
        one = self._Pn.one()
        if g1 == one:
            return self.monomial((tuple(t), g2))
        ret = self.monomial((tuple(t), g1))
        # We have to reverse the reduced word due to Sage's convention
        #   for permutation multiplication
        for i in g2.reduced_word():
            ret = self.linear_combination((self._product_by_basis_gen(m, i), c)
                                          for m,c in ret)
        return ret

    def _product_by_basis_gen(self, m, i):
        r"""
        Return the product `t g_w g_i`.

        If the quadratic relation is `g_i^2 = 1 + (q + q^{-1})e_i g_i`,
        then we have

        .. MATH::

            g_w g_i = \begin{cases}
            g_{ws_i} & \text{if } \ell(ws_i) = \ell(w) + 1, \\
            g_{ws_i} - (q - q^{-1}) g_w e_i & \text{if }
            \ell(w s_i) = \ell(w) - 1.
            \end{cases}

        INPUT:

        - ``m`` -- a pair ``[t, w]``, where ``t`` encodes the monomial
          and ``w``  is an element of the permutation group
        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(4, 3)
            sage: m = ((1, 0, 2), Permutations(3)([2,1,3]))
            sage: 4 * Y._product_by_basis_gen(m, 1)
            -(q^-1-q)*t2*t3^2*g[1] + 4*t1*t3^2 - (q^-1-q)*t1*t3^2*g[1]
             - (q^-1-q)*t1^2*t2^3*t3^2*g[1] - (q^-1-q)*t1^3*t2^2*t3^2*g[1]
        """
        t, w = m
        # We have to flip the side due to Sage's multiplication
        #   convention for permutations
        wi = w.apply_simple_reflection(i, side="left")
        if not w.has_descent(i, side="left"):
            return self.monomial((t, wi))

        R = self.base_ring()
        c = (self._q - ~self._q) * ~R(self._d)
        d = {(t, wi): R.one()}
        # We commute g_w and e_i and then multiply by t
        for s in range(self._d):
            r = list(t)
            r[w[i-1]-1] = (r[w[i-1]-1] + s) % self._d
            if s != 0:
                r[w[i]-1] = (r[w[i]-1] + self._d - s) % self._d
            d[(tuple(r), w)] = c
        return self._from_dict(d, remove_zeros=False)

    @cached_method
    def inverse_g(self, i):
        r"""
        Return the inverse of the generator `g_i`.

        From the quadratic relation, we have

        .. MATH::

            g_i^{-1} = g_i - (q - q^{-1}) e_i.

        EXAMPLES::

            sage: Y = algebras.YokonumaHecke(2, 4)
            sage: [2*Y.inverse_g(i) for i in range(1, 4)]
            [(q^-1+q) + 2*g[1] + (q^-1+q)*t1*t2,
             (q^-1+q) + 2*g[2] + (q^-1+q)*t2*t3,
             (q^-1+q) + 2*g[3] + (q^-1+q)*t3*t4]
        """
        if i < 1 or i >= self._n:
            raise ValueError("invalid index")
        return self.g(i) + (~self._q + self._q) * self.e(i)

    class Element(CombinatorialFreeModule.Element):
        def inverse(self):
            r"""
            Return the inverse if ``self`` is a basis element.

            EXAMPLES::

                sage: Y = algebras.YokonumaHecke(3, 3)
                sage: t = prod(Y.t()); t
                t1*t2*t3
                sage: ~t
                t1^2*t2^2*t3^2
                sage: [3*~(t*g) for g in Y.g()]
                [(q^-1+q)*t2*t3^2 + (q^-1+q)*t1*t3^2
                   + (q^-1+q)*t1^2*t2^2*t3^2 + 3*t1^2*t2^2*t3^2*g[1],
                 (q^-1+q)*t1^2*t3 + (q^-1+q)*t1^2*t2
                   + (q^-1+q)*t1^2*t2^2*t3^2 + 3*t1^2*t2^2*t3^2*g[2]]
            """
            if len(self) != 1:
                raise NotImplementedError("inverse only implemented for basis elements (monomials in the generators)"%self)
            H = self.parent()
            t,w = self.support_of_term()
            telt = H.monomial( (tuple((H._d - e) % H._d for e in t), H._Pn.one()) )
            return telt * H.prod(H.inverse_g(i) for i in reversed(w.reduced_word()))

        __invert__ = inverse

