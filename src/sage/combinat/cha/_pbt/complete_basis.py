# -*- coding: utf-8 -*-
r"""
The complete basis of PBT Hopf algebra.
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


class Complete(PlanarBinaryTreeFunctions.Bases.Base):
        '''

        '''
        _prefix = "H"

        def build_morphisms(self):
            # P <-> H
            P = self.realization_of().P()
            H_to_P = self.module_morphism(
                on_basis=lambda bt: \
                    P.sum_of_monomials(bt.sylvestrohedron_smaller()),
                triangular="lower",
                codomain=P
            )
            H_to_P.register_as_coercion()
            (~H_to_P).register_as_coercion()

        def product_on_basis(self, t1, t2):
            r"""
            The product of the H basis

            TESTS::

                sage: H = PBT(QQ).Complete()
                sage: H[3,1,2] * H[4,5,2,1,3]
                H[1, 5, 4, 7, 8, 6, 3, 2]
            """
            return self.monomial(t1 / t2)
