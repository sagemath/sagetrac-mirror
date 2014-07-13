# -*- coding: utf-8 -*-
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.wqsym import WordQuasiSymmetricFunctions


class Nabla(WordQuasiSymmetricFunctions._Basis):
    """
    TESTS::

        sage: N = WQSym(QQ).N()
        sage: N[1,2,1] * N[1]
        N[1, 2, 1, 3] + N[1, 2, 3, 1] + N[1, 3, 1, 2] + 2*N[1, 3, 2, 1] + N[2, 3, 2, 1] + N[3, 1, 2, 1]

    """

    _prefix_ = "N"

    def _morphisms_(self):
        """
        TESTS::

            sage: G = WQSym(QQ).G()
            sage: N = WQSym(QQ).N()
            sage: N(G[1,1,2,3])
            N[1, 1, 2, 3] + N[1, 1, 3, 2] + N[1, 2, 1, 3] + N[1, 2, 3, 1] + N[1, 3, 1, 2] + N[1, 3, 2, 1] + N[2, 1, 1, 3] + N[2, 1, 3, 1] + N[2, 3, 1, 1] + N[3, 1, 1, 2] + N[3, 1, 2, 1] + N[3, 2, 1, 1]
            sage: N(G[1,3,2,1])
            N[1, 3, 2, 1] + N[3, 1, 2, 1] + N[3, 2, 1, 1]

        """
        # ça semble correcte.
        G = self.realization_of().G()

        G_to_N = G.module_morphism(
            on_basis=lambda pw: self.sum_of_monomials(
                        pw.bruhat_greater(side="right")
            ), codomain=self,
            triangular="lower"
        )
        G_to_N.register_as_coercion()
        (~G_to_N).register_as_coercion()

    
