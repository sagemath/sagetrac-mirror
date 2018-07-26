# -*- coding: utf-8 -*-
r"""
The elementary of PBT Hopf algebra.
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>.
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
from sage.combinat.cha.pbt import PlanarBinaryTreeFunctions


class Elementary(PlanarBinaryTreeFunctions.Bases.Base):
        r'''
        TESTS::
        '''
        _prefix = "E"

        def build_morphisms(self):
            # P <-> E
            P = self.realization_of().P()
            E_to_P = self.module_morphism(
                on_basis=lambda bt: \
                    P.sum_of_monomials(bt.sylvestrohedron_greater()),
                triangular="upper",
                codomain=P
            )
            E_to_P.register_as_coercion()
            (~E_to_P).register_as_coercion()

        def product_on_basis(self, t1, t2):
            r"""
            The product of the E basis

            TESTS::

                sage: E = PBT(QQ).Elementary()
                sage: E[3,1,2] * E[4,5,2,1,3]
                E[1, 3, 2, 5, 4, 7, 8, 6]
                sage: E[4,3,5,1,2] * E[4,3,1,2]
                E[1, 4, 3, 5, 2, 6, 9, 8, 7]
            """
            return self.monomial(t1._backslash_(t2))
