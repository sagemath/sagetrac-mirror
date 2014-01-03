# -*- coding: utf-8 -*-
r"""
The combinatorial Hopf algebra of Word Quasi-Symmetric functions

References:
-----------

.. [NovThi06] Polynomial realizations of some trialgebras,
    J.-C. Novelli, and J.-Y. Thibon,
    2006.
}
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.bidendriform_bialgebras import BidendriformBialgebras
from sage.categories.category import Category
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.combinat.hopf_algebras import GenericGradedConnexeHopfAlgebras, words_like_getitem, \
    GenericBasisOfGCHopfAlgebra
from sage.combinat.hopf_algebras.categories.diese_product import DieseProductAlgebras
from sage.combinat.ncsf_qsym.generic_basis_code import GradedModulesWithInternalProduct
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.packed_word import PackedWords
from sage.misc.lazy_attribute import lazy_attribute


class WordQuasiSymmetricFunctions(GenericGradedConnexeHopfAlgebras):
    r"""

    """

    def the_category(self, R):
        return Category.join((
            HopfAlgebras(R).Graded().Connected().WithRealizations(),
            BidendriformBialgebras(R).WithRealizations()
        ))

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
            return [self.base().Realizations(),
                HopfAlgebrasWithBasis(R).Graded().Connected().Realizations(),
                BidendriformBialgebras(R).WithBasis().Realizations(),
                DieseProductAlgebras(R).WithBasis().Realizations(),
                GradedModulesWithInternalProduct(R).WithBasis().Realizations()]

        class Base(GenericBasisOfGCHopfAlgebra):
            _prefix = "** TO DEFINE **"

            @lazy_attribute
            def _basis_indices(self):
                return PackedWords()

            __getitem__ = words_like_getitem
