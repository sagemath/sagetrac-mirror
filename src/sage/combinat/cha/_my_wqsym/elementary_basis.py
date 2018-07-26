# -*- coding: utf-8 -*-
r"""
The fundamental basis of WQSym Hopf algebra.
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
from sage.combinat.cha.wqsym import WordQuasiSymmetricFunctions


class Elementary(WordQuasiSymmetricFunctions.Bases.Base):
        '''

        '''
        _prefix = "Me"

#        def build_morphisms(self):
#            '''
#            TESTS::
#
#                sage: M = WQSym(QQ).M(); Me = WQSym(QQ).Me()
#                sage: M(Me[2,1,1])
#                M[2, 1, 1] + M[3, 2, 1]
#                sage: Me(M[2,1,1])
#                Me[2, 1, 1] - Me[3, 2, 1]
#            '''
#            import sage.combinat.packed_words
#            M = self.realization_of().M()
#            self.module_morphism(
#                on_basis=lambda pw: \
#                    M.sum_of_monomials(pw.quasi_permutohedron_greater()),
#                codomain=M,
#                triangular="lower",
#                cmp=sage.combinat.packed_words.quasi_cmp
#            ).register_as_coercion()

        def product_on_basis(self, e1, e2):
            '''
            TESTS::

                sage: Me = WQSym(QQ).Me()
                sage: Me[1,1]*Me[1]
                Me[1, 1, 2]
            '''
            if len(e1) == 0:
                return self.monomial(e2)
            return self.monomial(
                self.basis().keys()(list(e1) + [i + max(e1) for i in e2])
            )
