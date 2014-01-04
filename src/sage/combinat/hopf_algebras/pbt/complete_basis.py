# -*- coding: utf-8 -*-
r"""
The complete basis of PBT Hopf algebra.
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.pbt import PlanarBinaryTreeFunctions


class Complete(PlanarBinaryTreeFunctions.Bases.Base):
    """
    The complete basis of PBT.

    .. see: [HNT05]_ §4.5 Multiplicative bases.

    EXAMPLES::

        sage: H = PBT(QQ).H()
        sage: a = H[3,1,2]; ascii_art(a)
        H
           o
          / \
         o   o
        sage: b = H[1,2,3,5,4]; ascii_art(b)
        H
               o
              / \
             o   o
            /
           o
          /
         o
        sage: ascii_art(a * b)
        H
           __o__
          /     \
         o       o
                  \
                   o
                  / \
                 o   o
                /
               o
              /
             o

    EXAMPLES::

        sage: H = PBT(QQ).H(); P = PBT(QQ).P()
        sage: a = P(H[2,1,3]); ascii_art(a)
        P    + P
           o        o
          /        /
         o        o
          \      /
           o    o
        sage: ascii_art(H(a))
        H
           o
          /
         o
          \
           o

    """
    _prefix = "H"

    def build_morphisms(self):
        """
        TESTS::

            sage: H = PBT(QQ).H(); P = PBT(QQ).P()
            sage: a = P(H[2,1,3]); a
            P[2, 1, 3] + P[1, 2, 3]
            sage: H(a)
            H[2, 1, 3]
            sage: b = H.monomial(BinaryTrees(5).random_element())
            sage: c = P(b)
            sage: b == H(c)
            True
        """
        # P <-> H
        P = self.realization_of().P()
        H_to_P = self.module_morphism(
            on_basis=lambda bt: \
                P.sum_of_monomials(bt.tamari_smaller()),
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
