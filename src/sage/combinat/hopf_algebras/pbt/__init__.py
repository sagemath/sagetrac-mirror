# -*- coding: utf-8 -*-
"""
The combinatorial Hopf algebra of Planar Binary Tree

This module implements methods related to the Hopf algebra of planar binary
trees also called the Loday-Ronco Hopf algebra (see [HNT05]_, [LR02]_ and
[LR98]_).

AUTHOR:

- Jean-Baptiste Priez

References
----------

.. [HNT05] The algebra of binary search trees,
    Florent Hivert,
    Jean-Christophe Novelli and
    Jean-Yves Thibon

.. [LR98] Hopf algebra of the planar binary trees,
    Jean-Louis Loday and
    María O. Ronco

.. [LR02] Order structure on the algebra of permutations and of planar
    binary trees,
    Jean-Louis Loday and
    María O. Ronco

Description
-----------

The *Loday-Ronco Hopf algebra* is the planar binary tree Hopf algebra.
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.bidendriform_bialgebras import BidendriformBialgebras
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.category import Category
from sage.combinat.hopf_algebras import GenericGradedConnexeHopfAlgebras, GenericBasisOfGCHopfAlgebra
from sage.combinat.hopf_algebras.categories.diese_product import DieseProductAlgebras
from sage.combinat.ncsf_qsym.generic_basis_code import GradedModulesWithInternalProduct
from sage.misc.functional import parent
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.integer import Integer
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.binary_tree import BinaryTrees
from sage.sets.set import Set


class PlanarBinaryTreeFunctions(GenericGradedConnexeHopfAlgebras):

    def the_category(self, R):
        return Category.join((
            HopfAlgebras(R).Graded().Connected().WithRealizations(),
            BidendriformBialgebras(R).WithRealizations()
        ))

    def __init__(self, R, left_to_right=False):
        """
        INPUT:

            - *left_to_right*: Set if the canonical labelling of the
        binary tree is associated to a binary search tree insertion
        from the left to the right of a permutation.

        EXAMPLES::

            sage: P = PBT(QQ, False).P() # = PBT(QQ).P()
            sage: ascii_art(P[3,1,2])
            P
               o
              / \
             o   o
            sage: Pr = PBT(QQ, True).P()
            sage: ascii_art(Pr[3,1,2])
            P
               o
              /
             o
              \
               o
        """
        GenericGradedConnexeHopfAlgebras.__init__(self, R)
        self._left_to_right = left_to_right

    def _get_permutation(self, bt):
        """
        TESTS::

            sage: bt = BinaryTree([[],[]]); ascii_art(bt)
              o
             / \
            o   o
            sage: PBT(QQ, False)._get_permutation(bt)
            [1, 3, 2]
            sage: PBT(QQ, True)._get_permutation(bt)
            [2, 1, 3]
        """
        if self._left_to_right:
            sigma = []
            bt.canonical_labelling().iterative_pre_order_traversal(
                lambda node: sigma.append(node.label())
            )
            return Permutation(sigma)
        else:
            return bt.to_312_avoiding_permutation()

    def _get_sylvester_class(self, bt):
        """
        TESTS::

            sage: bt = BinaryTree([[],[]]); ascii_art(bt)
              o
             / \
            o   o
            sage: PBT(QQ, False)._get_sylvester_class(bt)
            [[1, 3, 2], [3, 1, 2]]
            sage: PBT(QQ, True)._get_sylvester_class(bt)
            [[2, 1, 3], [2, 3, 1]]

        """
        return map(lambda li: Permutation(li), bt.sylvester_class(self._left_to_right))

    def _get_tree(self, sigma):
        """
        TESTS::

            sage: ascii_art(PBT(QQ, False)._get_tree(Permutation([3,1,2])))
              o
             / \
            o   o
            sage: ascii_art(PBT(QQ, True)._get_tree(Permutation([3,1,2])))
              o
             /
            o
             \
              o
        """
        return sigma.binary_search_tree_shape(self._left_to_right)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: PlanarBinaryTreeFunctions(QQ)
            The combinatorial Hopf algebra of Planar Binary Trees Functions over the Rational Field
            sage: PBT(QQ, True)
            The combinatorial Hopf algebra of Planar Binary Trees Functions over the Rational Field

        """
        return "The combinatorial Hopf algebra of Planar Binary Trees " + \
            "Functions over the %s" % self.base_ring()

    def a_realization(self):
        return self.P()

    def dual(self):
        return self

    def indices(self):
        return BinaryTrees()

    class Bases(Category_realization_of_parent):

        def super_categories(self):
            R = self.base().base_ring()
            return [self.base().Realizations(),
                    HopfAlgebrasWithBasis(R).Graded().Connected().Realizations(),
                    BidendriformBialgebras(R).Realizations().WithBasis(),
                    DieseProductAlgebras(R).WithBasis().Realizations(),
                    GradedModulesWithInternalProduct(R).WithBasis().Realizations()
            ]

        class ParentMethods:

            def _get_permutation(self, bt):
                return self.realization_of()._get_permutation(bt)

            def _get_sylvester_class(self, bt):
                return self.realization_of()._get_sylvester_class(bt)

            def _get_tree(self, sigma):
                return self.realization_of()._get_tree(sigma)

            def counit_on_basis(self, tree):
                if tree.node_number() == 0:
                    return self.base().one()
                else:
                    return self.base().zero()

        class Base(GenericBasisOfGCHopfAlgebra):
            _prefix = "** TO DEFINE **"

            @lazy_attribute
            def _basis_indices(self):
                return BinaryTrees()

            def one_basis(self):
                return self.basis().keys()()

            def __getitem__(self, c, *rest):
                """
                TESTS::

                    sage: P = PBT(QQ).P()
                    sage: P[3,1,2]
                    P[1, 3, 2]
                    sage: P[[]]
                    P[]
                    sage: P[1]
                    P[1]
                    sage: P[BinaryTree([])]
                    P[1]
                    sage: P[BinaryTree([[],[]])]
                    P[1, 3, 2]
                    sage: P[3,3,2]
                    Traceback (most recent call last):
                    ...
                    AssertionError: it must be like a permutation or a tree

                """
                # case: tree with one node
                if isinstance(c, (int, Integer)):
                    assert(c == 1), "it must be like a permutation or a tree"
                    res = self.basis().keys()([])
                # case: if c is a permutation
                elif len(c) > 1 and all([isinstance(i, (int, Integer)) for i in c]):
                    assert(list(c) + list(rest) in Permutations()), \
                        "it must be like a permutation or a tree"
                    res = self._get_tree(Permutation(list(c) + list(rest)))
                # case: empty list -> 1
                elif len(c) == 0 and len(rest) == 0:
                    res = self.basis().keys()()
                # may be the others case...
                else:
                    res = {
                        BinaryTrees(): lambda x: x,
                        Permutations(): lambda x: self._get_tree(x)
                    }[parent(c)](c)
                return self.monomial(res)

            def _repr_term(self, bt):
                """
                TESTS::

                    sage: P = PBT(QQ).P()
                    sage: P[3,1,2]
                    P[1, 3, 2]
                """
                sigma = self._get_permutation(bt)
                return self.prefix() + str(sigma)

            def _latex_term(self, m):
                """
                TESTS::

                    sage: latex(PBT(QQ).P()[3,1,2])
                    \mathsf{P}_{\vcenter{\hbox{\scalebox{.3}
                    { { \newcommand{\nodea}{\node[draw,circle] (a) {$$}
                    ;}\newcommand{\nodeb}{\node[draw,circle] (b) {$$}
                    ;}\newcommand{\nodec}{\node[draw,circle] (c) {$$}
                    ;}\begin{tikzpicture}[auto]
                    \matrix[column sep=.3cm, row sep=.3cm,ampersand replacement=\&]{
                             \& \nodea  \&         \\
                     \nodeb  \&         \& \nodec  \\
                    };
                    <BLANKLINE>
                    \path[ultra thick, red] (a) edge (b) edge (c);
                    \end{tikzpicture}} }}}}
                """
                from sage.misc.latex import latex
                prefix = self.print_options()['prefix']
                if len(m) == 0:
                    return "1"
                return "\\mathsf{" + prefix + "}" + \
                    "_{\\vcenter{\\hbox{\\scalebox{.3}\n{" + latex(m) + "}}}}"

