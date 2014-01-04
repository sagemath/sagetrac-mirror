# -*- coding: utf-8 -*-
r"""
The fundamental dual basis of PBT Hopf algebra.

"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.pbt import PlanarBinaryTreeFunctions


class FundamentalDual(PlanarBinaryTreeFunctions.Bases.Base):
    """
    This is the fundamental dual basis of ``PBT``::

    .. see: [HNT05]_ §4.1 Duality.

        sage: Q = PBT(QQ).Q(); Q
        The combinatorial Hopf algebra of Planar Binary Trees Functions over the Rational Field on the FundamentalDual basis

    .. MATH::

        \mathbb{F}_{\sigma} =
        \sum_{w\in \mathfrak{A}^*; std(w) = \sigma} w

    where `\mathfrak A` is an infinite and totally ordered alphabet such
    that `\mathfrak A^*` is free monoid. (``FQSym`` defined as a
    sub-algebra of the free algebra `\mathbb{K} \\langle \mathfrak{A}
    \\rangle`.)

    The product of `(\mathbb{F}_\sigma)_{\sigma\in\mathfrak{G}}` is
    described by

    .. MATH::

        \mathbb{F}_{\sigma} \\times \mathbb{F}_{\mu} =
        \sum_{\gamma \in \sigma \Cup \mu} \mathbb{F}_{\gamma}

    where `\sigma \Cup \mu` is the shifted shuffle of `\sigma` and `\mu`
    with. We denote by `std` the standardization map which  associate to
    a word `w` a permutation `\sigma` (see ``Permutations``).

    EXAMPLES::

        sage: ascii_art(Q[1,2] * Q[2,1])
        Q      + Q      + Q      + Q        + Q        + Q
         o        o          o        o          _o_          o
          \        \          \      / \        /   \        / \
           o        o          o    o   o      o     o      o   o
            \      / \        /          \          /      /
             o    o   o      o            o        o      o
            /               /
           o               o

    And the coproduct is described by

    .. MATH::

        \Delta(\mathbb{F}_{\sigma}) =
        \sum_{i\in [n]} \mathbb{F}_{\sigma_{\mid i[}} \otimes 
        \mathbb{F}_{\sigma_{\mid [i}}

    where `[n] := [1,..,n]`, `i[` is the sub-interval of `[n]` defined by
    `[1,..,i-1]` and `[i` the sub-interval defined by `[i,..,n]`.

    EXAMPLES::

        sage: ascii_art(Q[3,1,2,4].coproduct())
        1 # Q      + Q  # Q    + Q    # Q    + Q      # Q  + Q      # 1
                 o    o      o      o      o      o      o        o
                /           /      /      /      / \             /
               o           o      o      o      o   o           o
              / \           \                                  / \
             o   o           o                                o   o

    (See [MalReut]_ and [NCSF-VI]_.)


    TESTS::

        sage: Q = PBT(QQ).Q()
        sage: #TestSuite(Q).run()
    """
    _prefix = "Q"

    def build_morphisms(self):
        self._morph_G_to_Q_FQSym()
        self._morph_Q_P()

    def _morph_G_to_Q_FQSym(self):
        """
        TESTS::

            sage: G = FQSym(QQ).G(); Q = PBT(QQ).Q()
            sage: Q(G[3,1,2])
            Q[1, 3, 2]
            sage: Q(G[1,2,3,6,5,4])
            Q[1, 2, 3, 6, 5, 4]

        """
        from sage.combinat.hopf_algebras.all import FQSym
        G = FQSym(self.base()).G()
        G.module_morphism(
            on_basis=lambda sigma: self(self._get_tree(sigma)),
            codomain=self
        ).register_as_coercion()

    def _morph_Q_P(self):
        """
        TESTS::
            sage: Q = PBT(QQ).Q(); P = PBT(QQ).P()
            sage: Q(P[2,3,1])
            Q[1, 3, 2]
            sage: Q(P[3,1,2])
            Q[2, 3, 1] + Q[1, 3, 2]
            sage: Q(P[4,2,1,3])
            Q[3, 2, 4, 1] + Q[1, 3, 4, 2] + Q[2, 1, 4, 3]
            sage: Q(P[5,4,2,1,3])
            Q[4, 3, 5, 2, 1] + Q[2, 4, 5, 3, 1] + Q[3, 2, 5, 4, 1] + Q[1, 4, 5, 3, 2] + Q[1, 3, 5, 4, 2] + Q[2, 1, 5, 4, 3]
            sage: Q(P[5,3,4,1,2])
            Q[2, 4, 5, 3, 1] + Q[3, 2, 5, 4, 1] + Q[2, 3, 5, 4, 1] + Q[1, 4, 5, 3, 2] + 2*Q[1, 3, 5, 4, 2] + Q[2, 1, 5, 4, 3] + Q[1, 2, 5, 4, 3]

        """
        P = self.realization_of().P()
        P.module_morphism(
            on_basis=lambda bt: self.sum_of_monomials(map(
                lambda sigma: self._get_tree(sigma.inverse()),
                self._get_sylvester_class(bt)
            )), codomain=self
            #### FIXME:: must be invertible ####
            #triangular="upper",
            #inverse_on_support=lambda bt: self._get_tree(self._get_permutation(bt).inverse()),
            #cmp=lambda t1, t2: self._get_permutation(t1) >= self._get_permutation(t2)
        ).register_as_coercion()

    def dual_basis(self):
        return self.realization_of().P()

    def product_on_basis(self, c1, c2):
        r"""
        The product of the F basis : concatenation of tree

        EXAMPLES::

            sage: Q = PBT(QQ).Q()
            sage: ascii_art(Q[2,1]*Q[3,1,2])
            Q          + Q          + Q        + Q          + Q          + Q          +
               o            _o_          _o_        _o_          __o__        __o__
              / \          /   \        /   \      /   \        /     \      /     \
             o   o        o     o      o     o    o     o      o       o    o       o
                  \            / \          /      \     \      \     /      \
                   o          o   o        o        o     o      o   o        o
                    \                       \                                  \
                     o                       o                                  o
            <BLANKLINE>
            Q          + Q          + Q          + Q
                 o            _o_          _o_          o
                / \          /   \        /   \        / \
               o   o        o     o      o     o      o   o
              /     \      /     /      / \          /
             o       o    o     o      o   o        o
                                                     \
                                                      o

        .. FIXME::

            improve the product : don't use permutation, compute directly
            on tree
            ... may be false ...
        """
        sigma = self._get_permutation(c1).inverse()
        mu = self._get_permutation(c2).inverse()
        return self.sum_of_monomials(map(
            lambda sigma: self._get_tree(sigma.inverse()),
            sigma.shifted_shuffle(mu)
        ))

    def coproduct_on_basis(self, p):
        r"""
        The coproduct of the Q basis : de-shuffle of tree

        TESTS::

            sage: Q = PBT(QQ).Q()
            sage: Q[[]].coproduct()
            Q[] # Q[]
            sage: Q[1].coproduct()
            Q[] # Q[1] + Q[1] # Q[]
            sage: Q[1,2].coproduct()
            Q[] # Q[1, 2] + Q[1] # Q[1] + Q[1, 2] # Q[]
            sage: Q[1,2,3].coproduct()
            Q[] # Q[1, 2, 3] + Q[1] # Q[1, 2] + Q[1, 2] # Q[1] + Q[1, 2, 3] # Q[]

        EXAMPLES::

            sage: Q = PBT(QQ).Q()
            sage: ascii_art(Q[6,4,5,2,1,3].coproduct())
            1 # Q            + Q  # Q          + Q    # Q      + Q    # Q      + Q
                   __o__        o      _o_        o      o          o      o        _o_
                  /     \             /   \        \      \        /      / \      /   \
                 o       o           o     o        o      o      o      o   o    o     o
                  \     / \               / \             / \      \               \
                   o   o   o             o   o           o   o      o               o
            <BLANKLINE>
             # Q    + Q          # Q  + Q            # 1
                o        __o__      o      __o__
                 \      /     \           /     \
                  o    o       o         o       o
                        \     /           \     / \
                         o   o             o   o   o
        """
        def restrict(bt):
            if bt.is_empty():
                yield (bt, bt)
                return

            for bm, bM in restrict(bt[1]):
                yield (self.basis().keys()([bt[0], bm]), bM)
            for bm, bM in restrict(bt[0]):
                yield (bm, self.basis().keys()([bM, bt[1]]))
        return self.tensor_square().sum_of_monomials(restrict(p))
