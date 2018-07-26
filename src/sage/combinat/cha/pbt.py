# -*- coding: utf-8 -*-
"""
Planar Binary Tree Hopf algebra

AUTHOR:

- Jean-Baptiste Priez
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.graded_hopf_algebras_with_basis import \
    GradedHopfAlgebrasWithBasis
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.category import Category
from sage.combinat.cha.tools.generic_basis import GenericBasis
from sage.structure.parent import Parent
from sage.rings.integer import Integer
from sage.misc.lazy_import import LazyImport
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.combinat.permutation import Permutation, Permutation_class, \
    Permutations
from sage.combinat.binary_tree import BinaryTrees
from sage.sets.set import Set

from categories.bidendriform_bialgebras import BidendriformBialgebras


class PlanarBinaryTreeFunctions(UniqueRepresentation, Parent):
    '''

    '''
    def __init__(self, R, left_to_right=False):
        '''
        EXAMPLES::

            sage: FQS = FreeQuasisymmetricFunctions(ZZ);FQS
            The Hopf algebra of free quasisymmetric functions over the Integer Ring
        '''
        from sage.categories.all import Rings
        assert(R in Rings()), "first argument must be a ring"
        self._base = R
        Parent.__init__(
            self, base=R,
            category=Category.join((
                GradedHopfAlgebras(R).WithRealizations(),
                BidendriformBialgebras(R).WithRealizations()
            ))
        )

        assert(isinstance(left_to_right, bool)), \
            "second argument must be `right` or `left`"
        self._left_to_right = left_to_right

    def BST(self, sigma):
        '''
        TESTS::

            sage: pbt = PBT(ZZ)
            sage: pbt.BST([3,1,2])
            [[., .], [., .]]
            sage: pbt = PBT(ZZ, True)
            sage: pbt.BST(Permutation([3,1,2]))
            [[., [., .]], .]
        '''
        if not isinstance(sigma, Permutation_class):
            sigma = Permutation(sigma)
        return sigma.binary_search_tree(self._left_to_right).shape()

    def canonical_permutation(self, bt):
        ''' Make the canonical standard word of bt

        TESTS::

            sage: P = PBT(QQ).P()
            sage: P.canonical_permutation(BinaryTree([[],[]]))
            [1, 3, 2]
            sage: P(BinaryTree([[],[]])) == P[1,3,2]
            True
            sage: P = PBT(QQ, True).P()
            sage: P(BinaryTree([[],[]])) == P[1,3,2]
            False
            sage: P(BinaryTree([[],[]])) == P[2,1,3]
            True
        '''
        return bt.canonical_permutation(self._left_to_right)

    def sylvester_class(self, bt):
        ''' Compute all permutations 'sigma' in the sylvester class
        associate to 'bt'.

        TESTS::

            sage: B = BinaryTree
            sage: b = B([B([]),B([])]);b # the shape tree 312
            [[., .], [., .]]
            sage: pbt = PBT(QQ, True)
            sage: pbt.sylvester_class(b)
            [[2, 1, 3], [2, 3, 1]]
            sage: pbt = PBT(QQ, False)
            sage: pbt.sylvester_class(b)
            [[3, 1, 2], [1, 3, 2]]
        '''
        def __classeSylv(forest, w=[]):
            '''
                Compute all permutations which forms the shape
                of bt by insertion in binary search tree.
            '''
            if len(forest) == 0:
                return [w]
            res = []
            for t in forest:
                tmp = list(forest)
                tmp.remove(t)
                for ti in t:
                    if not ti.is_empty():
                        tmp += [ti]
                if self._left_to_right:
                    res += __classeSylv(tmp, w + [int(t.label())])
                else:
                    res += __classeSylv(tmp, [int(t.label())] + w)
            return res
        # compute the canonical permutation : sigma
        # w = self.canonical_permutation( bt )
        if bt.node_number() == 0:
            return [self.canonical_permutation(bt)]
        return map(
            lambda li: Permutation(li),
            __classeSylv([bt.canonical_labelling()])
        )

    def _repr_(self):
        r"""
        EXAMPLES

        """
        return "The combinatorial Hopf algebra of Planar Binary Trees " + \
            "Functions over the %s" % self.base_ring()

    def a_realization(self):
        return self.P()

    def dual(self):
        return self

    P = Fundamental = LazyImport(
        "sage.combinat.cha._pbt.fundamental_basis", "Fundamental")
    Q = FundamentalDual = LazyImport(
        "sage.combinat.cha._pbt.fundamental_dual_basis", "FundamentalDual")
    H = Complete = LazyImport(
        "sage.combinat.cha._pbt.complete_basis", "Complete")
    E = Elementary = LazyImport(
        "sage.combinat.cha._pbt.elementary_basis", "Elementary")

    def indices(self):
        return BinaryTrees()

    class Bases(Category_realization_of_parent):

        def super_categories(self):
            R = self.base().base_ring()
            return [self.base().Realizations(),
                    GradedHopfAlgebrasWithBasis(R).Realizations(),
                    BidendriformBialgebras.WithBasis(R).Realizations()]

        class ParentMethods:

            def build_morphisms(self):
                '''
                Define morphisms associated to the current basis
                '''

            def __init_extra__(self):
                self.build_morphisms()

            def canonical_permutation(self, bt):
                return self.realization_of().canonical_permutation(bt)

            def BST(self, sigma):
                return self.realization_of().BST(sigma)

            def counit_on_basis(self, tree):
                if tree.node_number() == 0:
                    return self.base().one()
                else:
                    return self.base().zero()

        class Base(GenericBasis):
            _prefix = "** TO DEFINE **"
            _basis_indices = BinaryTrees()

            def one_basis(self):
                return self.basis().keys()()

            def __getitem__(self, c, *rest):
                if isinstance(c, (int, Integer)):
                    assert(c == 1), "it must be like a permutation or a tree"
                    res = self.basis().keys()([])
                elif len(c) > 1 and all(
                        [isinstance(i, (int, Integer)) for i in c]):
                    s = Set(c)
                    assert(min(s) == 1 and max(s) == s.cardinality()), \
                        "it must be like a permutation or a tree"
                    res = self.basis().keys()(self.BST(list(c) + list(rest)))
                elif len(c) == 0 and len(rest) == 0:
                    res = self.BST([])
                else:
                    res = {
                        BinaryTrees(): lambda x: x,
                        Permutations(): lambda x: self.BST(x)
                    }[parent(c)](c)
                return self.monomial(res)

            def _repr_term(self, c):
                s = self.canonical_permutation(c)
                return self.prefix() + str(s)

            def _latex_term(self, m):
                from sage.misc.latex import latex
                prefix = self.print_options()['latex_prefix']
                if len(m) == 0:
                    return "1"
                return prefix + \
                    "_{\\vcenter{\\hbox{\\scalebox{.3}\n{" + latex(m) + "}}}}"

            def indices(self):
                return self.realization_of().indices()
