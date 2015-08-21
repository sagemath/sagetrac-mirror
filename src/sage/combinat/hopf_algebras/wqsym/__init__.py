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
from sage.combinat.hopf_algebras import GenericGradedConnexeHopfAlgebra, \
    words_like_getitem
from sage.combinat.hopf_algebras.categories.diese_product import \
    DieseProductAlgebras
from sage.combinat.hopf_algebras.categories.polynomial_realization import \
    PolynomialRealizationAlgebras
from sage.combinat.ncsf_qsym.generic_basis_code import \
    GradedModulesWithInternalProduct
from sage.combinat.packed_word import PackedWords


class WordQuasiSymmetricFunctions(GenericGradedConnexeHopfAlgebra):
    r"""

    """
    _default_basis_indices_ = PackedWords()

    def _extra_categories_(self, R):
        return [BidendriformBialgebras(R).WithRealizations()]

    def a_realization(self):
        return self.Fundamental()

    def dual(self):
        return self

    def _repr_(self):
        return "The combinatorial Hopf algebra of Word Quasi-Symmetric" + \
            " Functions over the %s" % self.base_ring()

    class _Basis(GenericGradedConnexeHopfAlgebra._Basis):

        def _extra_categories_(self):
            R = self.realization_of().base_ring()
            return [
                BidendriformBialgebras(R).WithBasis().Realizations(),
                DieseProductAlgebras(R).WithBasis().Realizations(),
                GradedModulesWithInternalProduct(R).WithBasis().Realizations(),
                PolynomialRealizationAlgebras(R).WithBasis().Realizations()
            ]

        __getitem__ = words_like_getitem
