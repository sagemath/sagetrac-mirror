r"""
Convex Cone of Positive-Semidefinite Matrices
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.geometry.convex_set import ConvexSet_closed, ConvexSet_relatively_open

class PositiveSemidefiniteMatrices_base(UniqueRepresentation):


    def __init__(self, matrix_space):
        r"""
        TESTS::

            sage: from sage.geometry.semialgebraic.semidefinite import PositiveSemidefiniteMatrices, PositiveDefiniteMatrices
            sage: M = MatrixSpace(QQ, 2)
            sage: M_psd = PositiveSemidefiniteMatrices(M)
            sage: TestSuite(M_psd).run()
            sage: M_pd = PositiveDefiniteMatrices(M)
            sage: TestSuite(M_pd).run()
        """
        self._matrix_space = matrix_space

    def ambient_vector_space(self, base_field=None):
        r"""
        Return the ambient vector space.
        """
        return self._matrix_space

    ambient = ambient_vector_space

    def ambient_dim(self):
        return self.ambient_vector_space().dimension()

    def dim(self):
        # This is currently assuming a full matrix space
        n = self._matrix_space.ncols()
        return n * (n + 1) / 2

    def contains(self, point):
        if point not in self._matrix_space:
            return False
        return self._predicate(point)

    __contains__ = contains

class PositiveSemidefiniteMatrices(PositiveSemidefiniteMatrices_base, ConvexSet_closed):
    r"""
    The convex cone of positive semidefinite matrices

    EXAMPLES::

        sage: from sage.geometry.semialgebraic.semidefinite import PositiveSemidefiniteMatrices

    A positive-definite matrix::

        sage: A = matrix(QQ, [ [2,1],
        ....:                  [1,2] ] )
        sage: M_psd = PositiveSemidefiniteMatrices(A.parent()); M_psd
        <sage.geometry.semialgebraic.semidefinite.PositiveSemidefiniteMatrices object at ...>
        sage: A in M_psd
        True

    A positive-semidefinite (but not positive-definite) matrix::

        sage: A = matrix(QQ, [ [1,1],
        ....:                  [1,1] ] )
        sage: A in M_psd
        True

    And finally, an indefinite matrix::

        sage: A = matrix(QQ, [ [0,1],
        ....:                  [1,0] ] )
        sage: A in M_psd
        False
    """
    def relative_interior(self):
        return PositiveDefiniteMatrices(self._matrix_space)

    def _predicate(self, point):
        return point.is_positive_semidefinite()

    def is_relatively_open(self):
        return self.dimension() == 0

    def is_compact(self):
        return self.dimension() == 0

class PositiveDefiniteMatrices(PositiveSemidefiniteMatrices_base, ConvexSet_relatively_open):
    r"""
    The convex set of positive definite matrices

    EXAMPLES::

        sage: from sage.geometry.semialgebraic.semidefinite import PositiveDefiniteMatrices

    A positive-definite matrix::

        sage: A = matrix(QQ, [ [2,1],
        ....:                  [1,2] ] )
        sage: M_pd = PositiveDefiniteMatrices(A.parent()); M_pd
        <sage.geometry.semialgebraic.semidefinite.PositiveDefiniteMatrices object at ...>
        sage: A in M_pd
        True

    A positive-semidefinite (but not positive-definite) matrix::

        sage: A = matrix(QQ, [ [1,1],
        ....:                  [1,1] ] )
        sage: A in M_pd
        False

    And finally, an indefinite matrix::

        sage: A = matrix(QQ, [ [0,1],
        ....:                  [1,0] ] )
        sage: A in M_pd
        False
    """
    def closure(self):
        return PositiveSemidefiniteMatrices(self._matrix_space)

    def _predicate(self, point):
        return point.is_positive_definite()

    def is_closed(self):
        return self.dimension() == 0
