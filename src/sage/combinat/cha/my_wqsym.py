# -*- coding: utf-8 -*-
r"""
The combinatorial Hopf algebra of Word Quasi-Symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>,
#                          Rémi Maurice <maurice@univ-mlv.fr>.
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
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.lazy_import import LazyImport
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.graded_hopf_algebras_with_basis import \
    GradedHopfAlgebrasWithBasis
from sage.combinat.ncsf_qsym.generic_basis_code import \
    GradedModulesWithInternalProduct
from sage.combinat.packed_words import PackedWords

from tools.generic_basis import GenericBasis


class WordQuasiSymmetricFunctions(UniqueRepresentation, Parent):
    r'''

    '''

    def __init__(self, R):
        '''
        TESTS::

            sage: WordQuasiSymmetricFunctions(QQ)
            The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field
        '''
        from sage.categories.rings import Rings
        assert(R in Rings()), '%s must be a ring' % R
        Parent.__init__(
            self, base=R,
            category=GradedHopfAlgebras(R).WithRealizations()
        )

    S = Fundamental = LazyImport(
        "sage.combinat.cha._my_wqsym.fundamental_basis", "Fundamental")
    M = FundamentalDual = LazyImport(
        "sage.combinat.cha._my_wqsym.fundamental_dual_basis", "FundamentalDual")
    Me = Elementary = LazyImport(
        "sage.combinat.cha._my_wqsym.elementary_basis", "Elementary")
    Ms = Homogene = LazyImport(
        "sage.combinat.cha._my_wqsym.homogene_basis", "Homogene")
    L = LeftWeakOrder = LazyImport(
        "sage.combinat.cha._my_wqsym.left_weak_order_basis", "LeftWeakOrder")
    R = Homogene = LazyImport(
        "sage.combinat.cha._my_wqsym.right_weak_order_basis", "RightWeakOrder")
    GL = GreaterLeft = LazyImport(
        "sage.combinat.cha._my_wqsym.greater_left_basis", "GreaterLeft")
    SR = SmallerRight = LazyImport(
        "sage.combinat.cha._my_wqsym.smaller_right_basis", "SmallerRight")
    GR = GreaterRight = LazyImport(
        "sage.combinat.cha._my_wqsym.greater_right_basis", "GreaterRight")
    SL = SmallerLeft = LazyImport(
        "sage.combinat.cha._my_wqsym.smaller_left_basis", "SmallerLeft")
    GLD = GreaterLeftDual = LazyImport(
        "sage.combinat.cha._my_wqsym.greater_left_dual_basis", "GreaterLeftDual")
    SRD = SmallerRightDual = LazyImport(
        "sage.combinat.cha._my_wqsym.smaller_right_dual_basis", "SmallerRightDual")
        

    def _repr_(self):
        return "The combinatorial Hopf algebra of Word Quasi-Symmetric" + \
            " Functions over the %s" % self.base_ring()

    def dual(self):
        return self

    def a_realization(self):
        return self.Fundamental()

    class Bases(Category_realization_of_parent):

        def super_categories(self):
            R = self.base().base_ring()
            return [GradedHopfAlgebras(R).Realizations(),
                    GradedModulesWithInternalProduct(R).Realizations(),
                    GradedHopfAlgebrasWithBasis(R)]

        class ParentMethods:

            def build_morphisms(self):
                '''
                Define morphisms associated to the current basis
                '''

            def __init_extra__(self):
                self.build_morphisms()

        class Base(GenericBasis):
            _prefix = "** TO DEFINE **"
            _basis_indices = PackedWords()
