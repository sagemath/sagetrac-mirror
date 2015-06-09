# -*- coding: utf-8 -*-
r"""
The homogeneous basis of WQSym Hopf algebra.

Ms-basis of WQSym
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.wqsym import WordQuasiSymmetricFunctions


class Homogene(WordQuasiSymmetricFunctions.Bases.Base):
        """

        """
        _prefix = "Ms"

        def build_morphisms(self):
            """
            TESTS::

                sage: M = WQSym(QQ).M(); Ms = WQSym(QQ).Ms()
                sage: M(Ms[2,1,1])
                M[1, 1, 1] + M[1, 1, 2] + M[1, 2, 2] + M[1, 2, 3] + M[2, 1, 1] + M[2, 1, 2] + M[2, 1, 3] + M[3, 1, 2]
                sage: Ms(M[2,1,1])
                -Ms[1, 1, 1] + Ms[1, 1, 2] + Ms[2, 1, 1] - Ms[3, 1, 2]
            """
            M = self.realization_of().M()
            self.module_morphism(
                on_basis=lambda pw: \
                    M.sum_of_monomials(pw.pseudo_permutohedron_smaller()),
                codomain=M,
                triangular="upper",
                cmp=lambda p1, p2: p2.is_smaller_than(p1)
            ).register_as_coercion()

        def product_on_basis(self, e1, e2):
            """
            TESTS::

                sage: Ms = WQSym(QQ).Ms()
                sage: Ms[1,1]*Ms[1]
                Ms[2, 2, 1]
            """
            return self.monomial(e2.shifted_concatenation(e1, side="left"))
