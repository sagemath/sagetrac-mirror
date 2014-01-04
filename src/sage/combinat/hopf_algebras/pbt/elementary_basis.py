# -*- coding: utf-8 -*-
r"""
The elementary of PBT Hopf algebra.
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.pbt import PlanarBinaryTreeFunctions


class Elementary(PlanarBinaryTreeFunctions.Bases.Base):
    r"""
    The elementary basis of PBT.

    .. see: [HNT05]_ §4.5 Multiplicative bases.

    EXAMPLES::

        sage: E = PBT(QQ).E()
        sage: a = E[4,1,2,3]; ascii_art(a)
        E
             o
            / \
           o   o
          /
         o
        sage: b = E[2,1,4,3]; ascii_art(b)
        E
           _o_
          /   \
         o     o
          \
           o
        sage: ascii_art(a * b)
        E
                  _o__
                 /    \
               _o_     o
              /   \
             o     o
            / \
           o   o
          /
         o

    EXAMPLES::

        sage: P = PBT(QQ).P()
        sage: b = E[2,1,4,3]; ascii_art(b)
        E
           _o_
          /   \
         o     o
          \
           o
        sage: c = P(b); ascii_art(c)
        P        + P      + P
         o          o          _o_
          \          \        /   \
           o          o      o     o
            \        / \      \
             o      o   o      o
              \
               o
        sage: ascii_art(E(c))
        E
           _o_
          /   \
         o     o
          \
           o
    """
    _prefix = "E"

    def build_morphisms(self):
        """
        TESTS::

            sage: P = PBT(QQ).P(); E = PBT(QQ).E()
            sage: a = P(E[2,1,3]); a
            P[3, 2, 1] + P[2, 3, 1] + P[2, 1, 3]
            sage: E(a)
            E[2, 1, 3]
            sage: b = E.monomial(BinaryTrees(5).random_element())
            sage: c = P(b)
            sage: b == E(c)
            True
        """
        # P <-> E
        P = self.realization_of().P()
        E_to_P = self.module_morphism(
            on_basis=lambda bt: \
                P.sum_of_monomials(bt.tamari_greater()),
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
