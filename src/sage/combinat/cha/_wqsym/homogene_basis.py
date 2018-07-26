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


class Homogene(WordQuasiSymmetricFunctions.Bases.Base):
        '''

        '''
        _prefix = "Ms"

#        def build_morphisms(self):
#            '''
#            TESTS::
#
#                sage: M = WQSym(QQ).M(); Ms = WQSym(QQ).Ms()
#                sage: M(Ms[2,1,1])
#                M[1, 1, 1] + M[1, 1, 2] + M[1, 2, 2] + M[1, 2, 3] + M[2, 1, 1] + M[2, 1, 2] + M[2, 1, 3] + M[3, 1, 2]
#                sage: Ms(M[2,1,1])
#                -Ms[1, 1, 1] + Ms[1, 1, 2] + Ms[2, 1, 1] - Ms[3, 1, 2]
#            '''
#            import sage.combinat.packed_words
#            M = self.realization_of().M()
#            self.module_morphism(
#                on_basis=lambda pw: \
#                    M.sum_of_monomials(pw.quasi_permutohedron_smaller()),
#                codomain=M,
#                triangular="upper",
#                cmp=sage.combinat.packed_words.quasi_cmp
#            ).register_as_coercion()

        def product_on_basis(self, e1, e2):
            '''
            TESTS::

                sage: Ms = WQSym(QQ).Ms()
                sage: Ms[1,1]*Ms[1]
                Ms[2, 2, 1]
            '''
            if len(e2) == 0:
                return self.monomial(e1)
            return self.monomial(
                self.basis().keys()([i + max(e2) for i in e1] + list(e2)))
