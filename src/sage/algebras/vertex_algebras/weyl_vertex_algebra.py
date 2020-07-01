r"""
Weyl Vertex Algebra

Given a commutative ring `R`, a free `R`-module `M` and a
non-degenerate, skew-symmetric, bilinear pairing
`\langle \cdot,\cdot\rangle: M \otimes_R M \rightarrow R`. The *Weyl*
vertex algebra associated to this datum is the vertex algebra generated
by `M` with `\lambda`-brackets given by:

.. MATH::

    [v_\lambda w] = \langle v, w\rangle |0\rangle.

This is not an H-graded vertex algebra. The choice of a
Lagrangian decomposition `M = L \oplus L^*` determines an H-graded
structure, however half the generators have conformal weight ``0`` and
therefore the graded pieces of this vertex algebra are not finite
dimensional.

.. SEEALSO::

    :mod:`Bosonic Ghosts Lie conformal algebra<sage.algebras.lie_conformal_algebras.bosonic_ghosts_lie_conformal_algebra>`

AUTHORS:

- Reimundo Heluani (2019-08-09): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
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

class WeylVertexAlgebra(UniversalEnvelopingVertexAlgebra):

    def __init__(self, R, ngens=None, gram_matrix=None, names=None,
                 index_set=None):
        r"""
        The Weyl vertex algebra.

        INPUT:

        - ``R`` -- a commutative ring; the base ring of this vertex
          algebra

        - ``ngens``: an even positive Integer (Default `2`); The number
          of generators of this vertex algebra.

        - ``gram_matrix`` -- a non-singular, skew-symmetric square
          matrix with values in ``R`` (default: `None`); the Gram
          matrix of the inner product of the generators

        - ``names`` -- a list or tuple of ``str``; alternative names
          for the generators

        - ``index_set`` -- an enumerated set; alternative indexing set
          for the generators

        OUTPUT:

        The Weyl Lie conformal algebra with generators
         `\alpha_i`, `i=1,...,ngens` and `\lambda`-brackets

        .. MATH::

            [{\alpha_i}_{\lambda} \alpha_j] = M_{ij} |0\rangle,

        where `M` is the ``gram_matrix`` above.

        .. NOTE::

            The returned vertex algebra is not `H`-graded.

        EXAMPLES::

            sage: V = vertex_algebras.Weyl(QQ, 4)
            sage: V.inject_variables()
            Defining alpha0, alpha1, alpha2, alpha3
            sage: alpha0.bracket(alpha2*alpha2)
            {0: 2*alpha2_(-1)|0>}
        """
        from sage.algebras.lie_conformal_algebras.weyl_lie_conformal_algebra import\
                                                    WeylLieConformalAlgebra
        ML = lie_conformal_algebras.Weyl(R, ngens=ngens, gram_matrix=gram_matrix,
                                      names=names, index_set=index_set)
        cp = Family({ML.gen(-1):R.one()})
        super(WeylVertexAlgebra,self).__init__(R, ML,
                                            central_parameters=cp)

    def _repr_(self):
        """
        The name of this vertex algebra.

        EXAMPLES::

            sage: V = vertex_algebras.Weyl(QQ, 4); V
            The Weyl vertex algebra with generators (alpha0_(-1)|0>, alpha1_(-1)|0>, alpha2_(-1)|0>, alpha3_(-1)|0>) over Rational Field
        """
        return "The Weyl vertex algebra with generators {} over {}".format(
                self.gens(),self.base_ring())

    def gram_matrix(self):
        """
        The Gram matrix of the inner product of the generators.

        EXAMPLES::

            sage: V = vertex_algebras.Weyl(QQ, 4)
            sage: V.gram_matrix()
            [ 0  0| 1  0]
            [ 0  0| 0  1]
            [-----+-----]
            [-1  0| 0  0]
            [ 0 -1| 0  0]
        """
        return self._lca.gram_matrix()
