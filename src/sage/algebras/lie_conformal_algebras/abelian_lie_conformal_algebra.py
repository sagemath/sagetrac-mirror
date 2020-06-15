"""
Abelian Lie Conformal Algebra.

AUTHORS:

- Reimundo Heluani (06-15-2020): Initial implementation.
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

from .graded_lie_conformal_algebra import GradedLieConformalAlgebra
from sage.structure.indexed_generators import standardize_names_index_set

class AbelianLieConformalAlgebra(GradedLieConformalAlgebra):
    def __init__(self, R, ngens=1, weights=None,
                 parity=None, names=None, index_set=None):
        r"""
        The Abelian Lie conformal algebra.

        INPUT:

        - ``R`` -- a commutative ring.
        - ``ngens`` -- a positive integer (default: ``1``); the number
          of generators of this Lie conformal algebra
        - ``weights`` -- a list of numbers (default: ``1`` for each
          generator); the weights of the generators. The resulting 
          Lie conformal algebra is `H`-graded.
        - ``parity`` -- ``None`` or a list of ``0`` or ``1`` (default:
           ``None``); The parity of the generators. If not ``None`` the
           resulting Lie Conformal algebra is a Super Lie conformal
           algebra

        OUTPUT:

        The Abelian Lie conformal algebra with generators `a_i`,
        `i=1,...,ngens` and vanishing `\lambda`-brackets.

        EXAMPLES::

            sage: R = AbelianLieConformalAlgebra(QQ,2)
            sage: R.inject_variables()
            Defining a0, a1
            sage: a0.bracket(a1.T(2))
            {}

        .. TODO::

            implement its own class to speed up arithmetics in this
            case.
        """
        if (names is None) and (index_set is None):
            names = 'a'
            self._latex_names = tuple(r'a_{%d}' % i for i in range(ngens))

        names,index_set = standardize_names_index_set(names=names,
                                                      index_set=index_set,
                                                      ngens=ngens)
        abeliandict = {} 
        
        GradedLieConformalAlgebra.__init__(self, R, abeliandict, names=names,
                                           index_set=index_set, weights=weights,
                                           parity=parity)

    def _repr_(self):
        return "The Abelian Lie conformal algebra with generators {}"\
               " over {}".format(self.gens(), self.base_ring())

