# -*- coding: utf-8 -*-
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.wqsym import WordQuasiSymmetricFunctions


class Micron(WordQuasiSymmetricFunctions._Basis):

    _prefix_ = "Mi"

    def _morphisms_(self):
        """
        TESTS::

            sage: Mi = WQSym(QQ).Mi()
            sage: S = WQSym(QQ).S()
            sage: Mi(S[1,1,2,3])
            Mi[1, 1, 2, 3] + Mi[1, 1, 3, 2] + Mi[2, 2, 1, 3] + Mi[2, 2, 3, 1] + Mi[3, 3, 1, 2] + Mi[3, 3, 2, 1]
        """
        # FIXME: ça semble faux
        S = self.realization_of().S()

        S_to_Mi = S.module_morphism(
            on_basis=lambda pw: self.sum_of_monomials(
                        pw.bruhat_greater(side="left")
            ), codomain=self,
            triangular="lower"
        )
        S_to_Mi.register_as_coercion()
        (~S_to_Mi).register_as_coercion()