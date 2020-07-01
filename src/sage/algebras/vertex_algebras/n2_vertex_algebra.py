r"""
N=2 Super Vertex Algebra

The `N=2` super vertex algebra of central charge `c` is an extension of
the Virasoro vertex algebra of central charge `c` (with generator `L`)
by an even generator `J` which is primary of conformal weight `1` and
two odd generators `G_1,G_2` which are primary of conformal weight
`3/2`. The remaining `\lambda`-brackets are given by:

.. MATH::

    [J_\lambda J] &= \frac{\lambda}{3} c |0\rangle, \\
    [J_\lambda G_1] &= G_1, \\
    [J_\lambda G_2] &= -G_2, \\
    [{G_1}_\lambda G_1] &= [{G_2}_\lambda G_2 ] = 0, \\
    [{G_1}_\lambda G_2] &= L + \frac{1}{2} TJ + \lambda J +
    \frac{\lambda^2}{6}c |0\rangle.

AUTHORS:

- Reimundo Heluani (2020-06-09): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from .universal_enveloping_vertex_algebra import \
                                                UniversalEnvelopingVertexAlgebra
from sage.sets.family import Family

class N2VertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self, R, c=0):
        """
        The N=2 super vertex algebra.

        INPUT:

        - ``R`` -- a commutative ring; the base ring.
        - ``c`` -- an element of ``R`` (default: ``0``); the central
          charge.

        EXAMPLES::

            sage: V = vertex_algebras.N2(QQ, 1)
            sage: V.inject_variables()
            Defining L, J, G1, G2
            sage: G1*G2 + G2*G1
            L_-3|0>

        The topological twist is a Virasoro vector with central
        charge 0::

            sage: L2 = L - 1/2*J.T()
            sage: L2.bracket(L2) == {0: L2.T(), 1: 2*L2}
            True

        A singular vector in conformal weight 2::

            sage: V.find_singular(2)
            (L_-2|0> - 3/2*J_-1J_-1|0>,)
        """
        from sage.algebras.lie_conformal_algebras.\
             n2_lie_conformal_algebra import N2LieConformalAlgebra

        ML = lie_conformal_algebras.N2(R)
        cp = Family({ML.gen(-1):c})
        super(N2VertexAlgebra,self).__init__(R, ML, central_parameters=cp)
        self._c = c

    def _repr_(self):
        """
        The name of this vertex algebra.

        EXAMPLES::

            sage: V = vertex_algebras.N2(QQ); V
            The N=2 super vertex algebra of central charge 0 over Rational Field
        """
        return "The N=2 super vertex algebra of central charge {} over {}"\
                .format(self.central_charge(),self.base_ring())



