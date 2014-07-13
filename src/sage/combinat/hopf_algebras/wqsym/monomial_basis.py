#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.wqsym import WordQuasiSymmetricFunctions
from sage.combinat.hopf_algebras.wqsym.fundamental_dual_basis import \
    _restriction


class Monomial(WordQuasiSymmetricFunctions._Basis):

    _prefix_ = "G"

    def _morphisms_(self):
        """
        TESTS::

            sage: G =WQSym(QQ).G()
            sage: M = WQSym(QQ).M()
            sage: G(M[1,2,3])
            G[1, 1, 1] + G[1, 1, 2] + G[1, 2, 2] + G[1, 2, 3]
            sage: G(M[1,3,2])
            G[1, 2, 1] + G[1, 3, 2]
            sage: G(M[2,1,3])
            G[2, 1, 2] + G[2, 1, 3]
            sage: G(M[2,3,1])
            G[2, 2, 1] + G[2, 3, 1]
            sage: G(M[3,1,2])
            G[2, 1, 1] + G[3, 1, 2]

        ::

            sage: M(G[1,1,1])
            M[1, 1, 1] - M[1, 1, 2] - M[1, 2, 2] + M[1, 2, 3]
            sage: M(G[2,2,3,2,1,4])
            M[2, 2, 3, 2, 1, 3] + M[2, 2, 3, 2, 1, 4]
            sage: M(G[2,4,5,1,2,3])
            M[2, 3, 3, 1, 2, 2] + M[2, 3, 4, 1, 2, 2] + M[2, 4, 4, 1, 2, 3] + M[2, 4, 5, 1, 2, 3]
        """

        # FIXME: the coercion fails somewhere...
        M = self.realization_of().M()

        G_to_M = self.module_morphism(
            on_basis=lambda pw: M.sum_of_monomials(
                        pw.greater_zabrocki_bergeron()
            ), codomain=M,
            triangular="upper"
        )
        G_to_M.register_as_coercion()
        (~G_to_M).register_as_coercion()


    def product_on_basis(self, p1, p2):
        """
        TESTS::

            sage: G = WQSym(QQ).G()
            sage: G[1] * G[1,1]
            G[1, 2, 2] + G[2, 1, 1]
        """
        return self.sum_of_monomials(map(lambda osp: osp.to_packed_word(),
                   p1.to_ordered_set_partition().shifted_shuffle(
                    p2.to_ordered_set_partition())
        ))

    def coproduct_on_basis(self, e):
        """
        TESTS::

            sage: G = WQSym(QQ).G()
            sage: M[1,1].coproduct()
            G[] # G[1, 1] + G[1, 1] # G[]
            sage: M[[]].coproduct()
            G[] # G[]
            sage: G[2,1,2,1,4,3].coproduct()
            G[] # G[2, 1, 2, 1, 4, 3] + G[1, 1] # G[1, 1, 3, 2] + G[2, 1, 2, 1] # G[2, 1] + G[2, 1, 2, 1, 3] # G[1] + G[2, 1, 2, 1, 4, 3] # G[]
        """
        return self.tensor_square().sum_of_monomials(
            _restriction(self, e, i) for i in set([-1] + list(e))
        )

