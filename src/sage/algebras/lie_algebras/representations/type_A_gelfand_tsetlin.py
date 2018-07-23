r"""
Irreducible representations of the Lie algebra $\mathfrak{gl}_n$ implemented
using Gelfand-Tsetlin bases

AUTHORS:

- Michael Walter (2018-07-22): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Michael Walter <m.walter@uva.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import

from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.gelfand_tsetlin_patterns import (
    GelfandTsetlinPatterns,
    GelfandTsetlinPattern,
)
from sage.categories.action import Action
from sage.categories.lie_algebra_representations import LieAlgebraRepresentations
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.root_system import RootSystem
from sage.combinat.partition import Partition
from sage.functions.other import factorial
from sage.matrix.constructor import matrix, diagonal_matrix
from sage.misc.all import prod
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.all import QQ, ZZ
from sage.structure.element import Element
from sage.matrix.matrix_space import is_MatrixSpace
import sage.combinat.root_system.type_A as type_A
import operator


class IrreducibleRepresentation(CombinatorialFreeModule):
    r"""
    An irreducible finite-dimensional complex representation of $\mathfrak{gl}_n$
    with a distinguished orthogonal (but non-normalized!) basis, labeled by
    Gelfand-Tsetlin patterns. See [Mol2002]_ for more mathematical details.

    Let $E_{ij}$ denote the elementary matrix with one in the $i$-th row and $j$-th
    column. We will work with the following generators of $\mathfrak{gl}_n$:

    - $h_k = E_{kk}$ for $k=1,\dots,n$
    - $e_k = E_{k,k+1}$ for $k=1,\dots,n-1$
    - $f_k = E_{k+1,k}$ for $k=1,\dots,n-1$

    The representation matrices of all generators are defined over $\QQ$.

    TODO: Naming is unfortunate since sometimes $h_k$ is used for the Cartan-Weyl
    generator of the k-th simple root.

    INPUT:

    - ``highest_weight`` -- the highest weight (an element in the weight lattice of
      some root system of type A)
    - ``K`` -- the base field

    EXAMPLES::

        sage: V = lie_algebra_representation([3, 1], "A1")
        sage: V
        Lie algebra representation of ['A', 1] with highest weight (3, 1)

        sage: Lambda = RootSystem(['A', 2]).weight_lattice().fundamental_weights()
        sage: V = lie_algebra_representation(Lambda[1])
        sage: V
        Lie algebra representation of ['A', 2] with highest weight Lambda[1]

        sage: L = RootSystem(['A', 2]).ambient_lattice()
        sage: V = lie_algebra_representation(L([7, 3, 3]))
        sage: V
        Lie algebra representation of ['A', 2] with highest weight (7, 3, 3)

    TESTS::

        sage: for n in [2, 3]:
        ....:     for p in PartitionsInBox(n, 5):
        ....:         V = lie_algebra_representation(p, ["A", n-1])
        ....:         TestSuite(V).run()
    """

    def __init__(self, highest_weight, K):
        self._highest_weight = highest_weight
        self._root_system = highest_weight.parent().root_system
        self._cartan_type = self._root_system.cartan_type()
        self._n = self._cartan_type.rank() + 1
        if not isinstance(self._cartan_type, type_A.CartanType):
            raise ValueError("highest weight should be for Cartan type A")

        # Gelfand-Tsetlin patterns
        L = self._root_system.ambient_lattice()
        assert L.dimension() == self._n
        self._patterns = GelfandTsetlinPatterns(top_row=vector(L(highest_weight)))
        self._patterns = [
            GelfandTsetlinPattern(pat) for pat in self._patterns
        ]  # FIXME #25919

        # set up category
        cat = LieAlgebraRepresentations(K).WithBasis()
        cat = cat.FiniteDimensional()

        CombinatorialFreeModule.__init__(self, K, self._patterns, category=cat)

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: V = lie_algebra_representation([4, 3, 1], "A2")
            sage: V.cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def root_system(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: V = lie_algebra_representation([4, 3, 1], "A2")
            sage: V.root_system()
            Root system of type ['A', 2]
        """
        return self._root_system

    def highest_weight(self):
        r"""
        Return the highest weight of ``self``.

        EXAMPLES::

            sage: V = lie_algebra_representation([4, 3, 1], "A2")
            sage: V.highest_weight()
            (4, 3, 1)
        """
        return self._highest_weight

    def highest_weight_pattern(self):
        r"""
        Return Gelfand-Tsetlin pattern of highest weight.

        EXAMPLES::

            sage: V = lie_algebra_representation([4, 3, 1], "A2")
            sage: V.highest_weight_pattern()
            [[4, 3, 1], [4, 3], [4]]
        """
        L = self._root_system.ambient_lattice()
        top_row = list(vector(L(self._highest_weight)))
        return GelfandTsetlinPattern([top_row[: self._n - k] for k in range(self._n)])

    def highest_weight_vector(self):
        r"""
        Return the highest weight vector.

        EXAMPLES::

            sage: V = lie_algebra_representation([4, 3, 1], "A2")
            sage: V.highest_weight_vector()
            B[[[4, 3, 1], [4, 3], [4]]]
        """
        return self.term(self.highest_weight_pattern())

    def _repr_(self):
        return "Lie algebra representation of %s with highest weight %s" % (
            self._cartan_type,
            self._highest_weight,
        )

    @cached_method
    def h_on_basis(self, k, p):
        r"""
        Return the result of acting by the generator $h_k$ on a basis vector.

        INPUT:

        - ``k`` -- index of the generator
        - ``p`` -- Gelfand-Tsetlin pattern labeling the basis vector

        EXAMPLES::

            sage: V = lie_algebra_representation([3, 1], "A1")
            sage: mu0, mu1, mu2 = V.basis().keys(); mu0, mu1, mu2
            ([[3, 1], [1]], [[3, 1], [2]], [[3, 1], [3]])
            sage: V.h_on_basis(1, mu0), V.h_on_basis(2, mu0)
            (B[[[3, 1], [1]]], 3*B[[[3, 1], [1]]])
            sage: V.h_on_basis(1, mu1), V.h_on_basis(2, mu1)
            (2*B[[[3, 1], [2]]], 2*B[[[3, 1], [2]]])
            sage: V.h_on_basis(1, mu2), V.h_on_basis(2, mu2)
            (3*B[[[3, 1], [3]]], B[[[3, 1], [3]]])
        """
        w = p.weight()[k - 1]
        return self.term(p, w)

    @cached_method
    def e_on_basis(self, k, p):
        r"""
        Return the result of acting by the generator $e_k$ on a basis vector.

        INPUT:

        - ``k`` -- index of the generator
        - ``p`` -- Gelfand-Tsetlin pattern labeling the basis vector

        EXAMPLES::

            sage: V = lie_algebra_representation([3, 1], "A1")
            sage: mu0, mu1, mu2 = V.basis().keys(); mu0, mu1, mu2
            ([[3, 1], [1]], [[3, 1], [2]], [[3, 1], [3]])

            sage: V.e_on_basis(1, mu0)
            2*B[[[3, 1], [2]]]
            sage: V.e_on_basis(1, mu1)
            2*B[[[3, 1], [3]]]
            sage: V.e_on_basis(1, mu2)
            0
        """
        n = self._n

        def GT(k, i):
            return p[n - k][i - 1]

        def l(k, i):
            return p[n - k][i - 1] - i + 1

        result = self.zero()
        for i in range(1, k + 1):
            # check that replacing GT(k, i) by GT(k, i) + 1 results in a valid GT pattern
            if k > 1 and i > 1 and GT(k - 1, i - 1) == GT(k, i):
                continue
            if k < n and GT(k + 1, i) == GT(k, i):
                continue
            c = [list(row) for row in p]
            c[n - k][i - 1] += 1
            c = GelfandTsetlinPattern(c)

            # compute action
            numer = prod([l(k, i) - l(k + 1, j) for j in range(1, k + 2)])
            denom = prod([l(k, i) - l(k, j) for j in range(1, k + 1) if j != i])
            result += self.term(c, -QQ((numer, denom)))

        return result

    @cached_method
    def f_on_basis(self, k, p):
        r"""
        Return the result of acting by the generator $f_k$ on a basis vector.

        INPUT:

        - ``k`` -- index of the generator
        - ``p`` -- Gelfand-Tsetlin pattern labeling the basis vector

        EXAMPLES::

            sage: V = lie_algebra_representation([3, 1], "A1")
            sage: mu0, mu1, mu2 = V.basis().keys(); mu0, mu1, mu2
            ([[3, 1], [1]], [[3, 1], [2]], [[3, 1], [3]])
            sage: V.f_on_basis(1, mu0)
            0
            sage: V.f_on_basis(1, mu1)
            B[[[3, 1], [1]]]
            sage: V.f_on_basis(1, mu2)
            B[[[3, 1], [2]]]
        """
        n = self._n

        def GT(k, i):
            return p[n - k][i - 1]

        def l(k, i):
            return p[n - k][i - 1] - i + 1

        result = self.zero()
        for i in range(1, k + 1):
            # check that replacing GT(k, i) by GT(k, i) - 1 results in a valid GT pattern
            if k > 1 and i < k and GT(k, i) == GT(k - 1, i):
                continue
            if k < n and GT(k, i) == GT(k + 1, i + 1):
                continue
            c = [list(row) for row in p]
            c[n - k][i - 1] -= 1
            c = GelfandTsetlinPattern(c)

            # compute action
            numer = prod([l(k, i) - l(k - 1, j) for j in range(1, k)])
            denom = prod([l(k, i) - l(k, j) for j in range(1, k + 1) if j != i])
            result += self.term(c, QQ((numer, denom)))

        return result

    @cached_method
    def root_action(self, root):
        r"""
        Return the action of the Cartan-Weyl generator labeled by the given root.

        INPUT:

        - ``root`` -- the root labeling the Cartan-Weyl generator

        OUTPUT:

        An endomorphism of ``self``.

        EXAMPLES::

            sage: V = lie_algebra_representation([1, 0, 0], "A2")
            sage: alpha = V.root_system().root_lattice().simple_roots()
            sage: V.root_action(alpha[1]).matrix()
            [0 0 0]
            [0 0 0]
            [0 1 0]
            sage: V.root_action(alpha[2]).matrix()
            [0 0 0]
            [2 0 0]
            [0 0 0]
            sage: V.root_action(alpha[1] + alpha[2]).matrix()
            [0 0 0]
            [0 0 0]
            [2 0 0]
            sage: V.root_action(-alpha[1]).matrix()
            [0 0 0]
            [0 0 1]
            [0 0 0]
            sage: V.root_action(-alpha[2]).matrix()
            [  0 1/2   0]
            [  0   0   0]
            [  0   0   0]
            sage: V.root_action(-alpha[1] - alpha[2]).matrix()
            [  0   0 1/2]
            [  0   0   0]
            [  0   0   0]
        """

        def fn(elem):
            return self.act_by_root(root, elem)

        return self.module_morphism(function=fn, codomain=self)

    def act_by_root(self, root, elem):
        r"""
        Return result of acting by the Cartan-Weyl generator labeled by the
        given root on the given element.

        INPUT:

        - ``root`` -- the root labeling the Cartan-Weyl generator
        - ``elem`` -- the element to act on
        
        EXAMPLES::

            sage: V = lie_algebra_representation([1, 0, 0], "A2")
            sage: alpha = V.root_system().root_lattice().simple_roots()
            sage: mu0, mu1, mu2 = V.basis().keys()

            sage: V.act_by_root(alpha[1], V(mu0) + 3*V(mu1) + 7*V(mu2)) == 3*V(mu2)
            True
            sage: V.act_by_root(alpha[2], V(mu0) + 3*V(mu1) + 7*V(mu2)) == 2*V(mu1)
            True
            sage: V.act_by_root(alpha[1] + alpha[2], V(mu0) + 3*V(mu1) + 7*V(mu2)) == 2*V(mu2)
            True
            sage: V.act_by_root(-alpha[1], V(mu0) + 3*V(mu1) + 7*V(mu2)) == 7*V(mu1)
            True
            sage: V.act_by_root(-alpha[2], V(mu0) + 3*V(mu1) + 7*V(mu2)) == 3/2*V(mu0)
            True
            sage: V.act_by_root(-alpha[1] - alpha[2], V(mu0) + 3*V(mu1) + 7*V(mu2)) == 7/2*V(mu0)
            True
        """
        L = self._root_system.ambient_lattice()
        assert root in L.positive_roots() or root in L.negative_roots()
        x = tuple(vector(L(root)))
        i, j = x.index(1) + 1, x.index(-1) + 1
        if i < j:
            if j == i + 1:
                return elem.e(i)
            most = L.root(i - 1, j - 2)
            return self.act_by_root(most, elem.e(j - 1)) - self.act_by_root(
                most, elem
            ).e(j - 1)
        else:
            if i == j + 1:
                return elem.f(j)
            rest = L.root(i - 1, j)
            return self.act_by_root(rest, elem.f(j)) - self.act_by_root(rest, elem).f(j)

    def action(self, X):
        r"""
        Return the action of a Lie algebra element.

        INPUT:

        - ``X`` -- a Lie algebra element 

            TODO: Right now, ``X`` is just a matrix.

        OUTPUT:

        An endomorphism of ``self``.

        EXAMPLES::
        
            sage: V = lie_algebra_representation([1, 0, 0], "A2")
            sage: V.action(matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])).matrix()
            [  9   4 7/2]
            [ 12   5   4]
            [  6   2   1]
        """
        L = self._root_system.ambient_lattice()
        m = diagonal_matrix(
            [
                vector(p.weight()).dot_product(vector(X.diagonal()))
                for p in self._patterns
            ]
        )
        for i in range(self._n):
            for j in range(i + 1, self._n):
                m += X[i, j] * self.root_action(L.root(i, j)).matrix()
                m += X[j, i] * self.root_action(L.root(j, i)).matrix()
        return self.module_morphism(matrix=m, codomain=self)

    def _get_action_(self, G, op, self_on_left):
        if is_MatrixSpace(G) and op == operator.mul and not self_on_left:
            return LieAlgebraAction(G, self)

    @cached_method
    def basis_norm_squared(self, p):
        r"""
        Return norm squared of a basis vector.

        INPUT:

        - ``p`` -- Gelfand-Tsetlin pattern labeling the basis vector

        EXAMPLES::

            sage: V = lie_algebra_representation([3, 1], "A1")
            sage: mu0, mu1, mu2 = V.basis().keys(); mu0, mu1, mu2
            ([[3, 1], [1]], [[3, 1], [2]], [[3, 1], [3]])
            sage: V.basis_norm_squared(mu0)
            4
            sage: V.basis_norm_squared(mu1)
            2
            sage: V.basis_norm_squared(mu2)
            1

            sage: V.basis_norm_squared(V.highest_weight_pattern())
            1
        """
        n = self._n

        def l(k, i):
            return p[n - k][i - 1] - i + 1

        return prod(
            prod(
                QQ(
                    (
                        factorial(l(k, i) - l(k - 1, j)),
                        factorial(l(k - 1, i) - l(k - 1, j)),
                    )
                )
                for i in range(1, k)
                for j in range(i, k)
            )
            * prod(
                QQ(
                    (
                        factorial(l(k, i) - l(k, j) - 1),
                        factorial(l(k - 1, i) - l(k, j) - 1),
                    )
                )
                for i in range(1, k + 1)
                for j in range(i + 1, k + 1)
            )
            for k in range(2, n + 1)
        )

    def _test_commutators(self, **options):
        """
        Verify the defining property of a Lie algebra representation for the
        action of two random elements and their commutant.
        """
        tester = self._tester(**options)
        n = self._n

        A = matrix.random(ZZ, n, n)
        B = matrix.random(ZZ, n, n)
        a = self.action(A)
        b = self.action(B)
        tester.assertEqual(
            self.action(A.commutator(B)).matrix(), a.matrix().commutator(b.matrix())
        )

    def _test_adjoint(self, **options):
        """
        Test that the actions of $E_alpha$ and $E_{-alpha}$ are each other's
        adjoints.
        """
        tester = self._tester(**options)
        L = self._root_system.ambient_lattice()
        for alpha in L.positive_roots():
            D = diagonal_matrix([self.basis_norm_squared(b) for b in self._patterns])
            D_inv = diagonal_matrix(
                [1 / self.basis_norm_squared(b) for b in self._patterns]
            )
            E_alpha = self.root_action(alpha).matrix()
            F_alpha = self.root_action(-alpha).matrix()
            tester.assertEqual(D_inv * F_alpha.transpose() * D, E_alpha)

    def _test_highest_weight_vector_normalized(self, **options):
        """
        test that the highest weight vector has unit norm.
        """
        tester = self._tester(**options)
        tester.assertEqual(self.basis_norm_squared(self.highest_weight_pattern()), 1)

    def _test_irrep_with_right_highest_weight(self, **options):
        """
        Test that there is a unique highest weight vector with the appropriate
        highest weight.
        """
        tester = self._tester(**options)
        L = self._root_system.ambient_lattice()
        roots = L.positive_roots()
        V = self.root_action(roots[0]).matrix().right_kernel()
        for alpha in roots[1:]:
            K = self.root_action(alpha).matrix().right_kernel()
            V = V.intersection(K)

        # single highest weight vector
        tester.assertEqual(V.dimension(), 1)

        # correct highest weight
        v = self.from_vector(V.basis()[0])
        lam = vector(L(self.highest_weight()))
        for k in range(self._n):
            tester.assertEqual(v.h(k + 1), lam[k] * v)


class IrreducibleRepresentationElement(CombinatorialFreeModule.Element):
    def inner_product(self, other):
        r"""
        Return the inner product between ``self`` and ``other``.

        EXAMPLES::

            sage: V = lie_algebra_representation([3, 1], "A1")
            sage: mu0, mu1, mu2 = V.basis().keys(); mu0, mu1, mu2
            ([[3, 1], [1]], [[3, 1], [2]], [[3, 1], [3]])
            sage: v = V(mu0) - V(mu1) + 3 * V(mu2)
            sage: w = V(mu0) + V(mu1) + 3 * V(mu2)
            sage: v.inner_product(w)
            11
        """
        return sum(self.parent().basis_norm_squared(p) * c * other[p] for p, c in self)

    def h(self, k):
        r"""
        Return the result of acting by the generator $h_k$.

        INPUT:

        - ``k`` -- index of the generator

        EXAMPLES::

            sage: V = lie_algebra_representation([3, 1], "A1")
            sage: mu0, mu1, mu2 = V.basis().keys(); mu0, mu1, mu2
            ([[3, 1], [1]], [[3, 1], [2]], [[3, 1], [3]])
            sage: (V(mu0) + 3*V(mu1) + 7*V(mu2)).h(1) == V(mu0) + 6*V(mu1) + 21*V(mu2)
            True
            sage: (V(mu0) + 3*V(mu1) + 7*V(mu2)).h(2) == 3*V(mu0) + 6*V(mu1) + 7*V(mu2)
            True
        """
        V = self.parent()
        return V.linear_combination((V.h_on_basis(k, p), c) for p, c in self)

    def e(self, k):
        r"""
        Return the result of acting by the generator $e_k$.

        INPUT:

        - ``k`` -- index of the generator

        EXAMPLES::

            sage: V = lie_algebra_representation([3, 1], "A1")
            sage: mu0, mu1, mu2 = V.basis().keys(); mu0, mu1, mu2
            ([[3, 1], [1]], [[3, 1], [2]], [[3, 1], [3]])
            sage: (V(mu0) + 3*V(mu1) + 7*V(mu2)).e(1)
            2*B[[[3, 1], [2]]] + 6*B[[[3, 1], [3]]]
        """
        V = self.parent()
        return V.linear_combination((V.e_on_basis(k, p), c) for p, c in self)

    def f(self, k):
        r"""
        Return the result of acting by the generator $f_k$.

        INPUT:

        - ``k`` -- index of the generator

        EXAMPLES::

            sage: V = lie_algebra_representation([3, 1], "A1")
            sage: mu0, mu1, mu2 = V.basis().keys(); mu0, mu1, mu2
            ([[3, 1], [1]], [[3, 1], [2]], [[3, 1], [3]])

            sage: (V(mu0) + 3*V(mu1) + 7*V(mu2)).f(1)
            3*B[[[3, 1], [1]]] + 7*B[[[3, 1], [2]]]
        """
        V = self.parent()
        return V.linear_combination((V.f_on_basis(k, p), c) for p, c in self)

    dot_product = inner_product


IrreducibleRepresentation.Element = IrreducibleRepresentationElement


class LieAlgebraAction(Action):
    r"""
    The action of the Lie algebra $G = \mathfrak{gl}_n$ on a representation $V$.

    EXAMPLES::

        sage: V = lie_algebra_representation([1, 0, 0], "A2")
        sage: mu0, mu1, mu2 = V.basis().keys()
        sage: matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) * V(mu0)
        9*B[[[1, 0, 0], [0, 0], [0]]] + 12*B[[[1, 0, 0], [1, 0], [0]]] + 6*B[[[1, 0, 0], [1, 0], [1]]]
    """

    def __init__(self, G, V):
        Action.__init__(self, G, V, True, operator.mul)

    def _call_(self, X, elem):
        return self.codomain().action(X)(elem)
