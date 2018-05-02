"""
Ariki-Koike Algebras

The **Ariki-Koike algebras** were introduced by Ariki and Koike [AK94]_ as
a natural generalisation of the Iwahori-Hecke algebras of types `A` and `B`
(see class:`~sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra`).
Soon afterwards,  Brou\'e and Malle defined analogues of the Hecke
algebras for all complex reflection groups

Fix non-negative integers `r` an `n`. The Ariki-Koike algebras are
deformations of the group algebra of the complex reflection group
`G(r, 1, n) = \ZZ / r\ZZ \wr \mathfrak{S}_n`. If `R` is a ring containing a
*Hecke parameter* `q` and *cyclotomic parameters* `u_1, \ldots, i_r` then
the Ariki-Koike algebra `H_n(q, u_1, \ldots, u_r)` is the unital associative
`r`-algebra with generators `L_1, T_1, \ldots, T_{n-1}` an relations:

.. MATH::

    \begin{aligned}
        \prod_{i=0}^{r-1} (T_0 - u_i) & = 0, \\
        T_i^2 & = (q - 1) T_i + q \qquad \text{for } 1 \leq i < n, \\
        T_0 T_1 T_0 T_1 & = T_1 T_0 T_1 T_0, \\
        T_i T_j & = T_j T_i \qquad \text{if } |i - j| \geq 2, \\
        T_i T_{i+1} T_i = T_{i+1} T_i T_{i+1} \qquad \text{for } 1 \leq i < n.
    \end{aligned}

AUTHORS:

- Travis Scrimshaw (2016-04): initial version

REFERENCES:

.. [AK94] Susumu Ariki and Kazuhiko Koike.
   *A Hecke algebra of* `(\ZZ / r\ZZ) \wr \mathfrak{S}_n`
   *and construction of its irreducible representations*.
   Advances in Mathematics **106**, (1994) pp. 216-243.

.. [BM93] \M. Brou\'e and G. Malle, *Zyklotomische Heckealgebren*,
   Asterisque, **212** (1993), 119-89.

.. [MM98] Gunter Malle and Andrew Mathas.
   *Symmetric cyclotomic Hecke algebras*
   J. Algebra. **205** (1998) pp. 275-293.
"""

#*****************************************************************************
#  Copyright (C) 2016-2018 Travis Scrimshaw <tcsrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function

from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.all import ZZ
from sage.categories.algebras import Algebras
from sage.categories.rings import Rings
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutations
from sage.sets.family import Family

class ArikiKoikeAlgebra(CombinatorialFreeModule):
    r"""
    The Ariki-Koike algebra `H_{r,n}(q, u)`.

    Let `R` be an unital integral domain.
    Let `q, u_0, \ldots, u_{r-1} \in R` such that `q^{-1} \in R`.
    The *Ariki-Koike algebra* is the unital associative algebra
    `H_{r,n}(q, u)` generated by `T_0, \ldots, T_{n-1}` that satisfies
    the following relations:

    .. MATH::

        \begin{aligned}
            \prod_{i=0}^{r-1} (T_0 - u_i) & = 0, \\
            T_i^2 & = (q - 1) T_i + q \qquad \text{for } 1 \leq i < n, \\
            T_0 T_1 T_0 T_1 & = T_1 T_0 T_1 T_0, \\
            T_i T_j & = T_j T_i \qquad \text{if } |i - j| \geq 2, \\
            T_i T_{i+1} T_i = T_{i+1} T_i T_{i+1} \qquad \text{for } 1 \leq i < n.
        \end{aligned}

    Thus the Ariki-Koike algebra is a deformation of the group algebra of
    the complex reflection group `G(r, 1, n) = \ZZ / r\ZZ \wr \mathfrak{S}_n`.

    Next, we define **Jucys-Murphy elements**

    .. MATH::

        L_i = q^{-1} T_{i-1} \cdots T_1 T_0 T_1 \cdots T_{i-1}

    for `1 \leq i \leq n`. (Note: these element different by a power of `q` from
    the corresponding elements in [AK94]_, however, these elements are more
    commonly used because they lead to nicer representation theoretic
    formulas). Ariki and Koike [AK94]_ showed that `H_{r,n}(q, u)` is a free
    `R`-module with a basis given by

    .. MATH::

        \{L_1^c_i \cdots L_n^c_n T_w \mid w \in S_n, 0 \leq c_i < r\}.

    In particular, we have `\dim H_{r,n}(q,u) = r^n n! = |G(r, 1, n)|`.
    Moreover, we have `L_i L_j = L_i L_j` for all `1 \leq i, j \leq n`.

    The Ariki-Koike algebra `H_{r,n}(q, u)` can be considered as a quotient
    of the group algebra of the braid group for `G(r, 1, n)` by the ideal
    generated by `\prod_{i=0}^{r-1} (T_0 - u_i)` and `(T_i - q)(T_i + 1)`.
    Furthermore, `H_{r,n}(q, u)` can be constructed as a quotient of the
    extended affine Hecke algebra of type `A_{n-1}^{(1)}` by
    `\prod_{i=0}^{r-1} (X_1 - u_i)`.

    Since the Ariki-Koike algebra is a quotient of the group algebra
    of the braid group of `G(r, 1, n)`, we can recover the
    group algebra of `G(r, 1, n)` as follows. Consider
    `u = (1, \zeta_r, \ldots, \zeta_r^{r-1})`, where `\zeta_r` is
    a primitive `r`-th root of unity, then we have

    .. MATH::

        R G(r, 1, n) = H_{r,n}(1, u).

    INPUT:

    - ``r`` -- the maximum power of `L_i`
    - ``n`` -- the rank `S_n`
    - ``q`` -- (optional) an invertible element in a commutative ring;
      the default is `q \in R[q,q^{-1}]`, where `R` is the ring containing
      the variables ``u``
    - ``u`` -- (optional) the variables `u_1, \ldots, u_r`; the
      default is the generators of `\ZZ[u_1, \ldots, u_r]`
    - ``R`` -- (optional) a commutative ring containing ``q`` and ``u``; the
      default is the parent of `q` and `u_1, \ldots, u_r`

    EXAMPLES:

    We start by constructing an Ariki-Koike algebra where the
    values `q, u` are generic and do some computations::

        sage: H = algebras.ArikiKoike(3, 4)
        sage: H.inject_variables()
        Defining L1, L2, L3, L4, T1, T2, T3
        sage: T1 * T2 * T1 * T2
        q*T[2,1] - (1-q)*T[2,1,2]
        sage: T1 * L1 * T2 * L3 * T1 * T2
        -(q-q^2)*L2*L3*T[2] + q*L1*L2*T[2,1] - (1-q)*L1*L2*T[2,1,2]
        sage: L1^3
        u0*u1*u2 + ((-u0*u1-u0*u2-u1*u2))*L1 + ((u0+u1+u2))*L1^2
        sage: L3 * L2 * L1 
        L1*L2*L3
        sage: u = H.u()
        sage: q = H.q()
        sage: (q + 2*u[0]) * (T1 * T2) * L3
        (-2*u0+(2*u0-1)*q+q^2)*L3*T[1] + (-2*u0+(2*u0-1)*q+q^2)*L2*T[2]
         + (2*u0+q)*L1*T[1,2]

    We check the defining relations::

        sage: prod(L1 - val for val in u) == H.zero()
        True
        sage: L1 * T1 * L1 * T1 == T1 * L1 * T1 * L1
        True
        sage: T1 * T2 * T1 == T2 * T1 * T2
        True
        sage: T2 * T3 * T2 == T3 * T2 * T3
        True
        sage: L2 == q**-1 * T1 * L1 * T1
        True
        sage: L3 == q**-2 * T2 * T1 * L1 * T1 * T2
        True

    We construct an Ariki-Koike algebra with `u = (1, \zeta_3, \zeta_3^2)`,
    where `\zeta_3` is a primitive third root of unity::

        sage: F = CyclotomicField(3)
        sage: zeta3 = F.gen()
        sage: R.<q> = LaurentPolynomialRing(F)
        sage: H = algebras.ArikiKoike(3, 4, q=q, u=[1, zeta3, zeta3^2], R=R)
        sage: H.inject_variables()
        Defining L1, L2, L3, L4, T1, T2, T3
        sage: L1^3
        1
        sage: L2^3
        1 - (q^-1-1)*T[1] - (q^-1-1)*L1*L2^2*T[1] - (q^-1-1)*L1^2*L2*T[1]

    Next, we take `q = 1` to obtain the group algebra of `G(r, 1, n)`::

        sage: F = CyclotomicField(3)
        sage: zeta3 = F.gen()
        sage: H = algebras.ArikiKoike(3, 4, q=1, u=[1, zeta3, zeta3^2], R=F)
        sage: H.inject_variables()
        Defining L1, L2, L3, L4, T1, T2, T3
        sage: A = ColoredPermutations(3, 4).algebra(F)
        sage: s1, s2, s3, s0 = list(A.algebra_generators())
        sage: all(L^3 == H.one() for L in H.L())
        True
        sage: J = [s0, s3*s0*s3, s2*s3*s0*s3*s2, s1*s2*s3*s0*s3*s2*s1]
        sage: all(Ji^3 == A.one() for Ji in J)
        True

    """
    @staticmethod
    def __classcall_private__(cls, r, n, q=None, u=None, R=None):
        """
        Standardize input to ensure a unique representation.

        TESTS::

            sage: H1 = algebras.ArikiKoike(4, 3)
            sage: S = PolynomialRing(ZZ, 'u', 4)
            sage: R.<q> = LaurentPolynomialRing(S)
            sage: H2 = algebras.ArikiKoike(4, 3, q=q)
            sage: H3 = algebras.ArikiKoike(4, 3, q, S.gens(), R)
            sage: H1 is H2
            True
            sage: H2 is H3
            True
        """
        if u is None:
            if q is not None:
                R = q.parent()
            if R is None:
                R = PolynomialRing(ZZ, 'u', r)
                u = R.gens()
                if q is None:
                    R = LaurentPolynomialRing(R, 'q')
                    q = R.gen()
            else:
                u = PolynomialRing(ZZ, 'u', r).gens()
                if q is None:
                    q = 'q'
        else:
            if not isinstance(u, (list,tuple)):
                u = [u]*r
            if R is None:
                from sage.structure.element import get_coercion_model
                cm = get_coercion_model()
                if q is None:
                    R = cm.common_parent(*[val.parent() for val in u])
                    R = LaurentPolynomialRing(R, 'q')
                    q = R.gen()
                else:
                    R = cm.common_parent(q.parent(), *[val.parent() for val in u])
            elif q is None:
                q = 'q'
            u = [R(val) for val in u]
        if R not in Rings().Commutative():
            raise TypeError("base ring must be a commutative ring")
        q = R(q)
        u = tuple(u)
        return super(ArikiKoikeAlgebra, cls).__classcall__(cls, r, n, q, u, R)

    def __init__(self, r, n, q, u, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: TestSuite(H).run()
            sage: H = algebras.ArikiKoike(1, 4)
            sage: TestSuite(H).run()
            sage: H = algebras.ArikiKoike(2, 3)
            sage: TestSuite(H).run()
            sage: H = algebras.ArikiKoike(3, 4)
            sage: TestSuite(H).run() # long time
        """
        self._r = r
        self._n = n
        self._q = q
        self._u = u
        self._zero_tuple = tuple([0] * n) # it seems more efficient to copy this as we need it a lot
        self._Pn = Permutations(n)
        self._one_perm = self._Pn.one()
        import itertools
        C = itertools.product(*([range(r)]*n))
        indices = list(itertools.product(C, self._Pn))
        CombinatorialFreeModule.__init__(self, R, indices, prefix='T',
                                         category = Algebras(R).WithBasis())
        self._assign_names(self.algebra_generators().keys())

    def _repr_(self):
        """ 
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.ArikiKoike(5, 2)
            Ariki-Koike algebra of rank 5 and order 2
             with q=q and u=(u0, u1, u2, u3, u4)
             over Univariate Laurent Polynomial Ring in q
             over Multivariate Polynomial Ring in u0, u1, u2, u3, u4
             over Integer Ring
        """
        return "Ariki-Koike algebra of rank {} and order {} with q={} and u={} over {}".format(
            self._r, self._n, self._q, self._u, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 2)
            sage: latex(H)
            \mathcal{H}_{5,2}(q)
        """
        return "\\mathcal{H}_{%s,%s}(%s)"%(self._r, self._n, self._q)

    def _repr_term(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(4, 3)
            sage: H._repr_term( ((1, 0, 2), Permutation([3,2,1])) )
            'L1*L3^2*T[2,1,2]'
        """
        gen_str = lambda e: '' if e == 1 else '^%s'%e
        lhs = '*'.join('L%s'%(j+1) + gen_str(i) for j,i in enumerate(m[0]) if i > 0)
        redword = m[1].reduced_word()
        if not redword:
            if not lhs:
                return '1'
            return lhs
        rhs = 'T[{}]'.format(','.join(str(i) for i in redword))
        if not lhs:
            return rhs
        return lhs + '*' + rhs

    def _latex_term(self, m):
        r"""
        Return a latex representation for the basis element indexed by ``m``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(4, 3)
            sage: H._latex_term( ((1, 0, 2), Permutation([3,2,1])) )
            'L_{1} L_{3}^{2} T_{2} T_{1} T_{2}'
        """
        gen_str = lambda e: '' if e == 1 else '^{%s}'%e
        lhs = ' '.join('L_{%s}'%(j+1) + gen_str(i) for j,i in enumerate(m[0]) if i > 0)
        redword = m[1].reduced_word()
        if not redword:
            if not lhs:
                return '1'
            return lhs
        return lhs + ' ' + ' '.join("T_{%d}"%i for i in redword)

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: dict(H.algebra_generators())
            {'L1': L1, 'L2': L2, 'L3': L3, 'T1': T[1], 'T2': T[2]}

            sage: H = algebras.ArikiKoike(1, 4)
            sage: dict(H.algebra_generators())
            {'T1': T[1], 'T2': T[2], 'T3': T[3]}
        """
        d = {}
        if self._r != 1:
            for i in range(self._n):
                r = list(self._zero_tuple) # Make a copy
                r[i] = 1
                d['L%s'%(i+1)] = self.monomial( (tuple(r), self._one_perm) )
        G = self._Pn.group_generators()
        for i in range(1, self._n):
            d['T%s'%i] = self.monomial( (self._zero_tuple, G[i]) )
        return Family(sorted(d), lambda i: d[i])

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: H.gens()
            (L1, L2, L3, T[1], T[2])
        """
        return tuple(self.algebra_generators())

    @cached_method
    def one_basis(self):
        """
        Return the index of the basis element of `1`.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: H.one_basis()
            ((0, 0, 0), [1, 2, 3])
        """
        return (self._zero_tuple, self._one_perm)

    def q(self):
        """
        Return the variable `q`.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: H.q()
            q
        """
        return self._q

    def u(self):
        """
        Return the variables `u`.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: H.u()
            (u0, u1, u2, u3, u4)
        """
        return self._u

    def T(self, i=None):
        """
        Return the generator(s) `T_i`.

        INPUT:

        - ``i`` -- (default: ``None``) the generator `T_i` or if ``None``,
          then the list of all generators `T_i`

        EXAMPLES::

            sage: H = algebras.ArikiKoike(8, 3)
            sage: H.T(1)
            T[1]
            sage: H.T()
            [T[1], T[2]]
        """
        G = self.algebra_generators()
        if i is None:
            return [G['T%s'%j] for j in range(1, self._n)]
        return G['T%s'%i]

    def L(self, i=None):
        """
        Return the generator(s) `L_i`.

        INPUT:

        - ``i`` -- (default: ``None``) the generator `L_i` or if ``None``,
          then the list of all generators `L_i`

        EXAMPLES::

            sage: H = algebras.ArikiKoike(8, 3)
            sage: H.L(2)
            L2
            sage: H.L()
            [L1, L2, L3]

            sage: H = algebras.ArikiKoike(1, 3)
            sage: H.L(2)
            u + (-u*q^-1+u)*T[1]
            sage: H.L()
            [u,
             u + (-u*q^-1+u)*T[1],
             u + (-u*q^-1+u)*T[2] + (-u*q^-2+u*q^-1)*T[2,1,2]]
        """
        G = self.algebra_generators()
        if i is None:
            if self._r == 1:
                return [self._Li_power(j, 1) for j in range(1, self._n+1)]
            return [G['L%s'%j] for j in range(1, self._n+1)]
        if self._r == 1:
            return self._Li_power(i, 1)
        return G['L%s'%i]

    def dimension(self):
        """
        Return the dimension of ``self``.

        The dimension of `H_{r,n}(q, u)` is `r^n n!`.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(8, 3)
            sage: H.dimension()
            3072
            sage: H = algebras.ArikiKoike(6, 3)
            sage: H.dimension()
            1296
            sage: H = algebras.ArikiKoike(3, 5)
            sage: H.dimension()
            29160
        """
        from sage.functions.other import factorial
        return self._r**self._n * factorial(self._n)

    def some_elements(self):
        """
        Return a list of elements of ``self``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(4, 3)
            sage: H.some_elements()
            [2 + 2*T[2] + 3*T[1], T[2] + T[1] + L3 + L2 + L1,
             L1, L2, L3, T[1], T[2]]
        """
        G = self.algebra_generators()
        return [self.an_element(), self.sum(G)] + list(G)

    @cached_method
    def product_on_basis(self, m1, m2):
        """
        Return the product of the basis elements indexed by ``m1`` and ``m2``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(6, 3)
            sage: m = ((1, 0, 2), Permutations(3)([2,1,3]))
            sage: H.product_on_basis(m, m)
            q*L1*L2*L3^4

            sage: H = algebras.ArikiKoike(4, 3)
            sage: L1,L2,L3,T1,T2 = H.algebra_generators()
            sage: L1 * T1 * L1^2 * T1
            q*L1*L2^2 + (1-q)*L1^2*L2*T[1]
            sage: L1^2 * T1 * L1^2 * T1
            q*L1^2*L2^2 + (1-q)*L1^3*L2*T[1]
            sage: L1^3 * T1 * L1^2 * T1
            (-u0*u1*u2*u3+u0*u1*u2*u3*q)*L2*T[1]
             + ((u0*u1*u2+u0*u1*u3+u0*u2*u3+u1*u2*u3)+(-u0*u1*u2-u0*u1*u3-u0*u2*u3-u1*u2*u3)*q)*L1*L2*T[1]
             + ((-u0*u1-u0*u2-u1*u2-u0*u3-u1*u3-u2*u3)+(u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q)*L1^2*L2*T[1]
             + ((u0+u1+u2+u3)+(-u0-u1-u2-u3)*q)*L1^3*L2*T[1] + q*L1^3*L2^2

            sage: H = algebras.ArikiKoike(4, 3)
            sage: L1^2 * T1 * L1^3 * T1
            (-u0*u1*u2*u3+u0*u1*u2*u3*q)*L2*T[1]
             + ((u0*u1*u2+u0*u1*u3+u0*u2*u3+u1*u2*u3)+(-u0*u1*u2-u0*u1*u3-u0*u2*u3-u1*u2*u3)*q)*L1*L2*T[1]
             + ((-u0*u1-u0*u2-u1*u2-u0*u3-u1*u3-u2*u3)+(u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q)*L1^2*L2*T[1]
             + q*L1^2*L2^3
             + ((u0+u1+u2+u3)+(-u0-u1-u2-u3)*q)*L1^3*L2*T[1]
             + (1-q)*L1^3*L2^2*T[1]

            sage: L1^2 * T1*T2*T1 * L2 * L3 * T2
            (q-2*q^2+q^3)*L1^2*L2*L3 - (1-2*q+2*q^2-q^3)*L1^2*L2*L3*T[2]
             - (q-q^2)*L1^3*L3*T[1] + (1-2*q+q^2)*L1^3*L3*T[1,2]
             + q*L1^3*L2*T[2,1] - (1-q)*L1^3*L2*T[2,1,2]

            sage: H = algebras.ArikiKoike(2, 3)
            sage: L3 = H.L(3)
            sage: x = H.an_element()
            sage: (x * L3) * L3 == x * (L3 * L3)
            True
        """
        # Although it is tempting to make this recursive some care must be taken
        # here to ensure that the various "helper" methods" return linear
        # combinations of "standard" basis elements of the form (L,w), where
        # L is an n-tuple and w is a permutation because otherwise we may end up
        # in an infinite loop...

        # Product is of the form L1*T1*L2*T2: separate the L's and permutations
        L1,T1 = m1
        L2,T2 = m2

        if sum(L2) == 0:
            # Compute and return the product of T1 and T2, whilst fixing L
            return self._product_LTwTv(L1, T1, T2)

        # If T1 is trivial then we just have L1*L2*T2 we only need to rewrite
        # all of the "large" powers that appear in L1*L2. Unfortunately, this
        # will almost certainly introduce more T_w's and it will be recursive
        # because L_n^r, for example, will introduce many powers of L_k for k<n.
        if T1 == self._one_perm:
            Lbig = list(self._zero_tuple)   # separate the "big" and small
            Lsmall = list(self._zero_tuple) # powers of the Lk's
            for i in range(self._n):
                s = L1[i] + L2[i]
                if s < self._r:
                    Lsmall[i] = s
                else:
                    Lbig[i] = s
            if Lbig == list(self._zero_tuple):
                # if no big powers we only need to combine Lsmall and T2
                return self.monomial((tuple(Lsmall), T2))

            # The l variables all commute, so we can multiply them in any order
            # that we like. For improved efficiency, however, we move the Ls to
            # the left as soon as we can. For efficiency, we multiply the
            # "big" powers in the order L_n^N L_{n-1}^N...L_1^N as this
            # way we have to expand few powers the of the Lk's later.
            return (self.monomial((tuple(Lsmall), self._one_perm)) 
                    * prod(self._Li_power(i+1, Lbig[i])
                           for i in reversed(range(self._n)) if Lbig[i] > 0)
                    * self.monomial((self._zero_tuple, T2))
                    )

        # If we are still here then both T1 and L2 are non-trivial. Using the
        # method _product_Tw_L we expand the product T1*L2 as a linear
        # combination of standard basis elements using the method and then,
        # recursively, multiply on the left and right by L1 and T2,
        # respectively. In other words, we multiply as L1*(T1*L2)*T2.
        return ( self.monomial((L1, self._one_perm))
                 * self._product_Tw_L(T1, L2)
                 * self.monomial((self._zero_tuple, T2)) )

    def _product_LTwTv(self, L, w, v):
        r"""
        Return the product `L * T_w * Tv` as a linear combinations of
        terms of the form `L*T_x`.

        INPUT:

        - ``L`` -- an `n`-tuple
        - ``w`` -- the permutation ``w``
        - ``v`` -- the permutation ``v``

        The main point of this method is that it computes the product
        `L T_w T_v` and returns it as a linear combination of standard
        basis elements. That is, terms of the form `L T_x`. The monomial
        ``L`` does not play a role in this calculation and, instead, it
        is kept as a place holder for this "L-component" of the product.

        For this calculation the most important point is that

        .. MATH::

            T_i T_v = \begin{cases}
                T_{s_i v},              & \text{if } \ell(s_iv) > \ell(v),\\
                q T_{s_i v} + (q-1)T_v, & \text{if } \ell(s_iv) < \ell(v).
            \end{cases}

        This observation is used to rewrite the product `L T_w T_v`
        as a linear combination of standard basis elements. 

        This method is not intended to be called directly and, instead,
        is used by :meth:`product_on_basis`.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 4)
            sage: P4 = Permutations(4)
            sage: H._product_LTwTv((0, 3, 2, 4), P4([1,3,2,4]), P4([1,3,2,4]))
            q*L2^3*L3^2*L4^4 - (1-q)*L2^3*L3^2*L4^4*T[2]
            sage: H._product_LTwTv((0, 3, 2, 4), P4([1,3,2,4]), P4([1,3,4,2]))
            q*L2^3*L3^2*L4^4*T[3] - (1-q)*L2^3*L3^2*L4^4*T[2,3]
            sage: H._product_LTwTv((0, 3, 2, 4), P4([1,4,3,2]), P4([1,4,3,2]))
            q^3*L2^3*L3^2*L4^4 - (q^2-q^3)*L2^3*L3^2*L4^4*T[3]
             - (q^2-q^3)*L2^3*L3^2*L4^4*T[2]
             + (q-2*q^2+q^3)*L2^3*L3^2*L4^4*T[2,3]
             + (q-2*q^2+q^3)*L2^3*L3^2*L4^4*T[3,2]
             - (1-2*q+2*q^2-q^3)*L2^3*L3^2*L4^4*T[3,2,3]
        """
        Lwv = self.monomial( (tuple(L), v) )
        qm1 = self._q - self.base_ring().one()
        # TODO: Do not use elements of self but just a dictionary
        #   given by the permutations and the blas_dict methods.
        #   We can then add the L portion at the end to each term.
        for i in w.reduced_word()[::-1]:
            Lwvi = self.zero() # start from 0
            for (v,c) in Lwv:
                # We have to flip the side due to Sage's 
                # convention for multiplying permutations
                vi = v[1].apply_simple_reflection(i, side="right")
                if v[1].has_descent(i, side="right"):
                    Lwvi += self._from_dict({v: c * qm1, (v[0], vi): c * self._q},
                                            remove_zeros=False)
                else:
                    Lwvi += self._from_dict({(v[0], vi): c},
                                            remove_zeros=False)
            Lwv = Lwvi
        return Lwv

    def _product_Tw_L(self, w, L):
        r"""
        Given a permutation ``w`` and a monomial ``L`` return the product
        `T_w L` as a linear combination of terms of the form `L_v T_v`.

        To do this we write `w = s_{i_1} \cdots s_{i_k}` and then push each
        `T_{i_a}` past `L` using Lemma 3.2 of [MM98]_ (cf. Lemma 3.3 and
        Proposition 3.4 of [AK94]_), which says

        .. MATH::

            T_i L_i^a L_{i+1}^b = L_i^b L_{i+1}^a T_i + \begin{cases}
              (1-q) sum_{k=0}^{a-1} L_i^{a+k} L_{i+1}^{b-k}, &\text{if } a \leq b,\\
              (q-1) sum_{k=0}^{b-1} L_i^{b+k} L_{i+1}^{a-k}, &\text{if } a \geq b.
            \end{cases}

        Of course, `T_i` commutes with `L_k`, for `k \neq i,i+1`.

        This method is not intended to be called directly and, instead,
        is used by :meth:`product_on_basis`.

        INPUT:

        - ``w`` -- a permutation
        - ``L`` -- an tuple `(a_1, \ldots, a_n)`

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 4)
            sage: P4 = Permutations(4)
            sage: H._product_Tw_L(P4([1,3,2,4]), (0,2,2,0))
            L2^2*L3^2*T[2]
            sage: H._product_Tw_L(P4([1,3,2,4]), (0,1,3,0))
            -(1-q)*L2*L3^3 - (1-q)*L2^2*L3^2 + L2^3*L3*T[2]
            sage: H._product_Tw_L(P4([1,3,2,4]), (0,3,1,0))
            (1-q)*L2*L3^3 + L2*L3^3*T[2] + (1-q)*L2^2*L3^2
            sage: H._product_Tw_L(P4([1,3,2,4]), (2,3,1,3))
            (1-q)*L1^2*L2*L3^3*L4^3 + L1^2*L2*L3^3*L4^3*T[2] + (1-q)*L1^2*L2^2*L3^2*L4^3
        """
        wL = self.monomial( (L, self._one_perm) ) # initialise to L: this is what we will eventually return
        q = self._q
        one = q.parent().one()
        for i in w.reduced_word()[::-1]:
            iL = self.zero() # this will become T_i * L, written in standard form
            for (lv, c) in wL:
                L = list(lv[0]) # make a copy
                v = lv[1]
                a, b = L[i-1], L[i]
                L[i-1], L[i] = L[i], L[i-1] # swap L_i=L[i-1] and L_{i+1}=L[i]
                # the term L_1^{a_1}...L_i^{a_{i+1}}L_{i+1}^{a_i}...L_n^{a_n} T_i T_v
                # always appears
                iL += c * self._product_LTwTv(L, self._Pn.simple_reflections()[i], v) # need T_i*T_v

                if a < b:
                    Ls = [ list(L) for k in range(b-a) ] # make copies of L
                    for k in range(b-a):
                        Ls[k][i-1] = a + k
                        Ls[k][i] = b - k
                    c *= (q - one)
                    iL += self.sum_of_terms( ((tuple(l), v), c) for l in Ls )
                elif a > b:
                    Ls = [ list(L) for k in range(a-b) ] # make copies of L
                    for k in range(a-b):
                        Ls[k][i-1] = b + k
                        Ls[k][i] = a - k
                    c *= (one - q)
                    iL += self.sum_of_terms( ((tuple(l), v), c) for l in Ls )
            wL = iL # replace wL with iL and repeat
        return wL

    @cached_method
    def _Li_power(self, i, m):
        r"""
        Return `L_i^m`, where `m \geq 0`.

        To compute `L_i^m` we use Corollary 3.4 of [MM98]_ which says that

        .. MATH::

            L_i^m = q^{-1} T_{i-1} L_{i-1}^m T_{i-1}
              + (1 - q^{-1}) \sum_{c=1}^{m-1} L_i^c L_{i-1}^{m-c} T_{i-1}

        EXAMPLES::

            sage: H = algebras.ArikiKoike(3, 3)
            sage: for i in range(1,4): 
            ....:     for m in range(4): 
            ....:         print('L_{}^{} = {}'.format(i,m,H._Li_power(i,m)))
            L_1^0 = 1
            L_1^1 = L1
            L_1^2 = L1^2
            L_1^3 = u0*u1*u2 + ((-u0*u1-u0*u2-u1*u2))*L1 + ((u0+u1+u2))*L1^2
            L_2^0 = 1
            L_2^1 = L2
            L_2^2 = L2^2
            L_2^3 = u0*u1*u2 + (-u0*u1*u2*q^-1+u0*u1*u2)*T[1]
             + ((-u0*u1-u0*u2-u1*u2))*L2 + ((u0+u1+u2))*L2^2
             + ((u0+u1+u2)*q^-1+(-u0-u1-u2))*L1*L2*T[1]
             - (q^-1-1)*L1*L2^2*T[1] - (q^-1-1)*L1^2*L2*T[1]
            L_3^0 = 1
            L_3^1 = L3
            L_3^2 = L3^2
            L_3^3 = u0*u1*u2 + (-u0*u1*u2*q^-1+u0*u1*u2)*T[2]
            + (-u0*u1*u2*q^-2+u0*u1*u2*q^-1)*T[2,1,2]
            + ((-u0*u1-u0*u2-u1*u2))*L3 + ((u0+u1+u2))*L3^2
            + ((u0+u1+u2)*q^-1+(-u0-u1-u2))*L2*L3*T[2]
            - (q^-1-1)*L2*L3^2*T[2] - (q^-1-1)*L2^2*L3*T[2]
            + ((u0+u1+u2)*q^-2+(-2*u0-2*u1-2*u2)*q^-1+(u0+u1+u2))*L1*L3*T[1,2]
            + ((u0+u1+u2)*q^-2+(-u0-u1-u2)*q^-1)*L1*L3*T[2,1,2]
            - (q^-2-2*q^-1+1)*L1*L3^2*T[1,2] - (q^-2-q^-1)*L1*L3^2*T[2,1,2]
            - (q^-2-2*q^-1+1)*L1*L2*L3*T[1,2] - (q^-2-2*q^-1+1)*L1^2*L3*T[1,2]
            - (q^-2-q^-1)*L1^2*L3*T[2,1,2]
        """
        # shorthand for returning a tuple of the form (0,...,a,b,...,0) with a,b
        # in the (i-1)th and i-th positions, respectively
        def Ltuple(a, b):
            return tuple([b if j == i else a if j == i-1 else 0
                          for j in range(1,self._n+1)])

        # return "small" powers of the generators without change
        if m < self._r:
            return self.monomial( (Ltuple(0, m), self._one_perm) )

        def L(exp):
            return tuple([exp if j == i else 0 for j in range(1, self._n+1)])

        if i > 1:
            si = self._Pn.simple_reflections()[i-1]
            qsum = self.base_ring().one() - self._q**-1
            # by calling _Li_power we avoid infinite recursion here
            return ( self.sum_of_terms( ((Ltuple(c, m-c), si), qsum) for c in range(1, m) )
                     + self._q**-1 * self.T(i-1) * self._Li_power(i-1, m) * self.T(i-1) )

        # now left with the case i = 1 and m >= r
        if m > self._r:
            return self.monomial((Ltuple(0,1), self._one_perm)) * self._Li_power(i,m-1)

        z = PolynomialRing(self.base_ring(), 'DUMMY').gen()
        p = list(prod(z - val for val in self._u))[:-1]
        zero = self.base_ring().zero()
        return self._from_dict({(Ltuple(0,exp), self._one_perm): -coeff
                                for exp,coeff in enumerate(p) if coeff != zero},
                               remove_zeros=False)

    @cached_method
    def inverse_T(self, i):
        r"""
        Return the inverse of the generator `T_i`.

        From the quadratic relation, we have

        .. MATH::

            T_i^{-1} = q^{-1} T_i + (q^{-1} - 1).

        EXAMPLES::

            sage: H = algebras.ArikiKoike(3, 4)
            sage: [H.inverse_T(i) for i in range(1, 4)]
            [(q^-1-1) + (q^-1)*T[1],
             (q^-1-1) + (q^-1)*T[2],
             (q^-1-1) + (q^-1)*T[3]]

        TESTS::

            sage: H = algebras.ArikiKoike(4, 4)
            sage: all(H.inverse_T(i) * H.T(i) == H.one() for i in range(1, 4))
            True
            sage: all(H.T(i) * H.inverse_T(i) == H.one() for i in range(1, 4))
            True
        """
        c = ~self._q - self.base_ring().one()
        m = self.T(i).leading_support()
        return self._from_dict({m: ~self._q, self.one_basis(): c})

    class Element(CombinatorialFreeModule.Element):
        def inverse(self):
            r"""
            Return the inverse if ``self`` is a basis element.

            EXAMPLES::

                sage: H = algebras.ArikiKoike(3, 4)
                sage: t = prod(H.T()); t
                T[1,2,3]
                sage: t.inverse()
                (q^-3-3*q^-2+3*q^-1-1) + (q^-3-2*q^-2+q^-1)*T[3]
                 + (q^-3-2*q^-2+q^-1)*T[2] + (q^-3-q^-2)*T[3,2]
                 + (q^-3-2*q^-2+q^-1)*T[1] + (q^-3-q^-2)*T[1,3]
                 + (q^-3-q^-2)*T[2,1] + (q^-3)*T[3,2,1]
            """
            if len(self) != 1:
                raise NotImplementedError("inverse only implemented for monomials")
            l,w = self.support_of_term()
            if sum(l) != 0:
                raise NotImplementedError("inverse only implemented for monomials in T variables")
            H = self.parent()
            return ~self[l,w] * H.prod(H.inverse_T(i) for i in reversed(w.reduced_word()))

        __invert__ = inverse
