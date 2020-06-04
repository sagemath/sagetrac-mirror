"""
Bosonic Ghosts Lie Conformal Algebra

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation.
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

from sage.matrix.special import identity_matrix
from sage.structure.indexed_generators import standardize_names_index_set
from .graded_lie_conformal_algebra import GradedLieConformalAlgebra
class BosonicGhostsLieConformalAlgebra(GradedLieConformalAlgebra):

    def __init__(self,R,ngens=2,names=None,index_set=None):
        r"""
        The Bosonic ghosts or `\beta-\gamma`-system Lie conformal
        algebra.

        INPUT:

        - ``R``: a commutative ring.
        - ``ngens``: an even positive Integer (default: ``2``); the
          number of non-central generators of this Lie conformal
          algebra.

        EXAMPLES::

            sage: R = BosonicGhostsLieConformalAlgebra(QQ); R
            The Bosonic ghosts Lie conformal algebra with generators (beta, gamma, K) over Rational Field.
            sage: R.inject_variables(); beta.bracket(gamma)
            Defining beta, gamma, K
            {0: K}
            sage: beta.degree()
            1
            sage: gamma.degree()
            0

            sage: R = BosonicGhostsLieConformalAlgebra(QQbar, ngens = 4, names = 'abcd'); R
            The Bosonic ghosts Lie conformal algebra with generators (a, b, c, d, K) over Algebraic Field.
            sage: R.structure_coefficients()
            Finite family {('a', 'c'): ((0, K),),  ('b', 'd'): ((0, K),),  ('c', 'a'): ((0, -K),),  ('d', 'b'): ((0, -K),)}

        TESTS::

            sage: BosonicGhostsLieConformalAlgebra(AA).category()
            Category of finitely generated H-graded Lie conformal algebras with basis over Algebraic Real Field
        """
        try:
            assert (ngens > 0 and ngens % 2 == 0)
        except AssertionError:
            raise ValueError("ngens should be an even positive integer, " +
                             "got {}".format(ngens))
        if (names is None) and (index_set is None):
            from sage.misc.defaults import variable_names as varnames
            from sage.misc.defaults import latex_variable_names as laxnames
            names = varnames(ngens/2,'beta') + varnames(ngens/2,'gamma')
            self._latex_variable_names =  tuple(laxnames(ngens/2,r'\beta') +\
                                          laxnames(ngens/2,r'\gamma')) + ('K',)

        names,index_set = standardize_names_index_set(names=names,
                                                      index_set=index_set,
                                                      ngens=ngens)
        A = identity_matrix(R,ngens/2)
        from sage.matrix.special import block_matrix
        gram_matrix = block_matrix([[R.zero(),A],[-A,R.zero()]])
        ghostsdict = { (i,j): {0: {('K',0): gram_matrix[index_set.rank(i),
                    index_set.rank(j)]}} for i in index_set for j in index_set}
        weights = (1,)*(ngens//2) + (0,)*(ngens//2)
        super(BosonicGhostsLieConformalAlgebra,self).__init__(R,
                                           ghostsdict,names=names,
                                           index_set=index_set,
                                           weights=weights,
                                           central_elements=('K',))

    def _repr_(self):
        return "The Bosonic ghosts Lie conformal algebra with generators {} "\
               "over {}.".format(self.gens(),self.base_ring())

