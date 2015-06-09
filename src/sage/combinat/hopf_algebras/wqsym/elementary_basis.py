# -*- coding: utf-8 -*-
r"""
The elementary basis of WQSym Hopf algebra.

Me-basis of WQSym
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.wqsym import WordQuasiSymmetricFunctions


class Elementary(WordQuasiSymmetricFunctions.Bases.Base):
        """

        """
        _prefix = "Me"

        def build_morphisms(self):
            """
            TESTS::

                sage: M = WQSym(QQ).M(); Me = WQSym(QQ).Me()
                sage: M(Me[2,1,1])
                M[2, 1, 1] + M[3, 2, 1]
                sage: Me(M[2,1,1])
                Me[2, 1, 1] - Me[3, 2, 1]
            """
            M = self.realization_of().M()
            self.module_morphism(
                on_basis=lambda pw: M.sum_of_monomials(pw.pseudo_permutohedron_greater()),
                codomain=M,
                triangular="lower",
                cmp=lambda p1, p2: p1.is_smaller_than(p2)
            ).register_as_coercion()

        def product_on_basis(self, e1, e2):
            """
            TESTS::

                sage: Me = WQSym(QQ).Me()
                sage: Me[1,1]*Me[1]
                Me[1, 1, 2]
            """
            return self.monomial(e1.shifted_concatenation(e2))
