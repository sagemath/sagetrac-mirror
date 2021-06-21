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

import itertools
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.fields import Fields
from sage.matrix.matrix_space import is_MatrixSpace, MatrixSpace
from sage.geometry.convex_set import ConvexSet_closed_conic, ConvexSet_relatively_open

class SemidefiniteMatrices_base(UniqueRepresentation):
    r"""
    Base class for convex sets of semidefinite matrices.
    """

    @staticmethod
    def __classcall__(cls, *args, **kwds):
        """
        Normalize init arguments for :class:`UniqueRepresentation`.

        TESTS::

            sage: from sage.geometry.semialgebraic.semidefinite import (
            ....:    PositiveSemidefiniteMatrices, PositiveDefiniteMatrices)
            sage: PositiveDefiniteMatrices(QQ, 2) is PositiveDefiniteMatrices(MatrixSpace(QQ, 2))
            True
        """
        if len(args) == 1 and is_MatrixSpace(args[0]):
            matrix_space = args[0]
        else:
            matrix_space = MatrixSpace(*args, **kwds)
        if matrix_space.nrows() != matrix_space.ncols():
            raise TypeError('matrix space must be square')
        base_ring = matrix_space.base_ring()
        if base_ring not in Fields():
            raise NotImplementedError('semigroups of (semi)definite matrices over non-field rings are not implemented')
        if not matrix_space.nrows():
            # Trivial case of 0x0 matrices. The unique element has all properties.
            # We could just "return matrix_space" but then the result
            # would not have all methods that we provide.
            cls = PositiveSemidefiniteMatrices
        return super().__classcall__(cls, matrix_space)

    def __init__(self, matrix_space):
        r"""
        TESTS::

            sage: from sage.geometry.semialgebraic.semidefinite import (
            ....:    PositiveSemidefiniteMatrices, PositiveDefiniteMatrices)
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

        This is the full matrix space.

        EXAMPLES::

            sage: from sage.geometry.semialgebraic.semidefinite import (
            ....:     PositiveSemidefiniteMatrices, PositiveDefiniteMatrices)
            sage: PositiveSemidefiniteMatrices(QQ, 0).ambient_vector_space()
            Full MatrixSpace of 0 by 0 dense matrices over Rational Field
            sage: PositiveSemidefiniteMatrices(AA, 1).ambient_vector_space()
            Full MatrixSpace of 1 by 1 dense matrices over Algebraic Real Field
        """
        return self._matrix_space

    @cached_method
    def an_affine_basis(self):
        r"""
        Return a tuple of matrices that form an affine basis of ``self``.

        EXAMPLES::

            sage: from sage.geometry.semialgebraic.semidefinite import (
            ....:     PositiveSemidefiniteMatrices, PositiveDefiniteMatrices)
            sage: M_psd = PositiveDefiniteMatrices(QQ, 3)
            sage: B = M_psd.an_affine_basis(); B
            (
            [0 0 0]  [1 0 0]  [0 0 0]  [0 0 0]  [0 1 0]  [0 0 1]  [0 0 0]
            [0 0 0]  [0 0 0]  [0 1 0]  [0 0 0]  [1 0 0]  [0 0 0]  [0 0 1]
            [0 0 0], [0 0 0], [0 0 0], [0 0 1], [0 0 0], [1 0 0], [0 1 0]
            )
            sage: len(B) == 1 + M_psd.dimension()
            True

        """
        B = [self._matrix_space.zero()]
        n = self._matrix_space.nrows()
        def immutable_matrix(dict):
            M = self._matrix_space.matrix(dict)
            M.set_immutable()
            return M
        B.extend(immutable_matrix({(i, i): 1})
                 for i in range(n))
        B.extend(immutable_matrix({(i, j): 1, (j, i): 1})
                 for i, j in itertools.combinations(range(n), 2))
        return tuple(B)

    def affine_hull(self):
        r"""
        Return the affine hull of ``self``.

        It is the linear space of symmetric matrices.

        EXAMPLES::

            sage: from sage.geometry.semialgebraic.semidefinite import (
            ....:     PositiveSemidefiniteMatrices, PositiveDefiniteMatrices)
            sage: M_psd = PositiveDefiniteMatrices(QQ, 3); M_psd
            Cone of positive-definite matrices
             of Full MatrixSpace of 3 by 3 dense matrices over Rational Field
            sage: M_sym = M_psd.affine_hull(); M_sym
            Free module generated by {0, 1, 2, 3, 4, 5} over Rational Field
            sage: M_sym.dimension() == M_psd.dimension()
            True

        """
        return self._matrix_space.submodule(self.an_affine_basis()[1:])

    def dim(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.geometry.semialgebraic.semidefinite import (
            ....:     PositiveSemidefiniteMatrices, PositiveDefiniteMatrices)
            sage: PositiveSemidefiniteMatrices(QQ, 0).dimension()
            0
            sage: PositiveDefiniteMatrices(AA, 1).dimension()
            1
            sage: PositiveSemidefiniteMatrices(RR, 2).dimension()
            3
        """
        # This is currently assuming a full matrix space
        n = self._matrix_space.ncols()
        return n * (n + 1) // 2

    @abstract_method
    def _predicate(self, point):
        r"""
        Containment test for elements of the ambient space

        Subclasses must implement this method.

        TESTS::

            sage: from sage.geometry.semialgebraic.semidefinite import SemidefiniteMatrices_base
            sage: class NegativeDefiniteMatrices(SemidefiniteMatrices_base):
            ....:     def _predicate(self, point):
            ....:         return (-point).is_positive_definite()
            sage: M_nd = NegativeDefiniteMatrices(QQ, 2)
            sage: matrix.identity(2) in M_nd
            False
            sage: -matrix.identity(2) in M_nd
            True
        """

    def contains(self, point):
        r"""
        Test whether ``self`` contains the given ``point`` (matrix).

        EXAMPLES::

            sage: from sage.geometry.semialgebraic.semidefinite import PositiveDefiniteMatrices
            sage: M_pd_2 = PositiveDefiniteMatrices(QQ, 2); M_pd_2
            Cone of positive-definite matrices of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: matrix.identity(2) in M_pd_2
            True
            sage: -matrix.identity(2) in M_pd_2
            False
            sage: [[2, 1], [1, 2]] in M_pd_2    # via the element constructor
            True

        TESTS:

        No errors, just ``False``, if the dimensions of the matrix are wrong::

            sage: matrix.identity(0) in M_pd_2
            False

        The `0 \times 0` matrix is positive-definite because all of its
        eigenvalues are positive::

            sage: M_pd_0 = PositiveDefiniteMatrices(QQ, 0); M_pd_0
            Cone of positive-semidefinite matrices of Full MatrixSpace of 0 by 0 dense matrices over Rational Field
            sage: matrix(QQ, 0, 0, []) in M_pd_0
            True

        No errors, just ``False``, if the ``point`` does not make sense::

            sage: "positive" in M_pd_0
            False
        """
        try:
            point = self._matrix_space(point)
        except (TypeError, ValueError):
            return False
        return self._predicate(point)

    __contains__ = contains

    def _some_elements_(self):
        r"""
        Generate some elements in ``self``.

        EXAMPLES::

            sage: from sage.geometry.semialgebraic.semidefinite import PositiveSemidefiniteMatrices
            sage: M_psd = PositiveSemidefiniteMatrices(QQ, 2); M_psd
            Cone of positive-semidefinite matrices
             of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: M_psd.some_elements()  # indirect doctest
            [
            [13/2  3/4]  [1 0]  [  1 1/2]  [  1 1/2]  [0 0]
            [ 3/4    4], [0 0], [1/2   1], [1/2   1], [0 1]
            ]
        """
        for A in self._matrix_space.some_elements():
            t = 0
            A = (A + A.T) / 2
            while A not in self:
                t = 4 * t + 1
                A += t * self._matrix_space.identity_matrix()
            yield A


class PositiveSemidefiniteMatrices(SemidefiniteMatrices_base, ConvexSet_closed_conic):
    r"""
    The convex cone of positive-semidefinite symmetric matrices

    INPUT:

    - either an instance of :class:`MatrixSpace` or parameters for the constructor
      of :class:`MatrixSpace`

    EXAMPLES::

        sage: from sage.geometry.semialgebraic.semidefinite import PositiveSemidefiniteMatrices
        sage: PositiveSemidefiniteMatrices(QQ, 2)
        Cone of positive-semidefinite matrices
         of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: PositiveSemidefiniteMatrices(MatrixSpace(AA, 2, sparse=True))
        Cone of positive-semidefinite matrices
         of Full MatrixSpace of 2 by 2 sparse matrices over Algebraic Real Field
        sage: PositiveSemidefiniteMatrices(SR, 3)
        Cone of positive-semidefinite matrices
         of Full MatrixSpace of 3 by 3 dense matrices over Symbolic Ring

    We can test the containment of matrices in the cone.  First, a positive-definite matrix::

        sage: A = matrix(QQ, [ [2,1],
        ....:                  [1,2] ] )
        sage: M_psd = PositiveSemidefiniteMatrices(A.parent()); M_psd
        Cone of positive-semidefinite matrices
         of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: A in M_psd
        True
        sage: 1/17 * A in M_psd
        True

    A positive-semidefinite (but not positive-definite) matrix::

        sage: A = matrix(QQ, [ [1,1],
        ....:                  [1,1] ] )
        sage: A in M_psd
        True

    An indefinite matrix::

        sage: A = matrix(QQ, [ [0,1],
        ....:                  [1,0] ] )
        sage: A in M_psd
        False

    And, finally, a non-symmetric matrix::

        sage: A = matrix(QQ, [ [2,2],
        ....:                  [1,2] ] )
        sage: A in M_psd
        False

    Matrices with symbolic constants::

        sage: SRc = SR.subring(no_variables=True); SRc
        Symbolic Constants Subring
        sage: M_psd_SRc = PositiveSemidefiniteMatrices(SRc, 2); M_psd_SRc
        Cone of positive-semidefinite matrices
         of Full MatrixSpace of 2 by 2 dense matrices over Symbolic Constants Subring
        sage: A = matrix(SRc, [[e, sqrt(2)], [sqrt(2), e]]); A
        [      e sqrt(2)]
        [sqrt(2)       e]
        sage: A in M_psd_SRc
        True
    """
    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.geometry.semialgebraic.semidefinite import PositiveSemidefiniteMatrices
            sage: M_psd = PositiveSemidefiniteMatrices(AA, 2)
            sage: repr(M_psd)  # indirect doctest
            'Cone of positive-semidefinite matrices of Full MatrixSpace of 2 by 2 dense matrices over Algebraic Real Field'
        """
        return r"Cone of positive-semidefinite matrices of " + repr(self._matrix_space)

    def relative_interior(self):
        r"""
        Return the relative interior of ``self``.

        It is the convex set of positive-definite matrices.

        EXAMPLES::

            sage: from sage.geometry.semialgebraic.semidefinite import PositiveSemidefiniteMatrices
            sage: M_psd = PositiveSemidefiniteMatrices(AA, 2)
            sage: M_psd.relative_interior()
            Cone of positive-definite matrices
             of Full MatrixSpace of 2 by 2 dense matrices over Algebraic Real Field

            sage: M_psd = PositiveSemidefiniteMatrices(AA, 0)
            sage: M_psd.relative_interior() is M_psd
            True
        """
        return PositiveDefiniteMatrices(self._matrix_space)

    _predicate = attrcall('is_positive_definite')

    def is_relatively_open(self):
        return self.dimension() == 0


class PositiveDefiniteMatrices(SemidefiniteMatrices_base, ConvexSet_relatively_open):
    r"""
    The convex set of positive-definite symmetric matrices

    EXAMPLES::

        sage: from sage.geometry.semialgebraic.semidefinite import PositiveDefiniteMatrices

    A positive-definite matrix::

        sage: A = matrix(QQ, [ [2,1],
        ....:                  [1,2] ] )
        sage: M_pd = PositiveDefiniteMatrices(A.parent()); M_pd
        Cone of positive-definite matrices
         of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: A in M_pd
        True

    A positive-semidefinite (but not positive-definite) matrix::

        sage: A = matrix(QQ, [ [1,1],
        ....:                  [1,1] ] )
        sage: A in M_pd
        False

    An indefinite matrix::

        sage: A = matrix(QQ, [ [0,1],
        ....:                  [1,0] ] )
        sage: A in M_pd
        False

    And, finally, a non-symmetric matrix::

        sage: A = matrix(QQ, [ [2,2],
        ....:                  [1,2] ] )
        sage: A in M_pd
        False
    """
    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.geometry.semialgebraic.semidefinite import PositiveDefiniteMatrices
            sage: M_pd = PositiveDefiniteMatrices(AA, 2)
            sage: repr(M_pd)  # indirect doctest
            'Cone of positive-definite matrices of Full MatrixSpace of 2 by 2 dense matrices over Algebraic Real Field'
        """
        return r"Cone of positive-definite matrices of " + repr(self._matrix_space)

    def closure(self):
        r"""
        Return the topological closure of ``self``.

        It is the convex cone of positive-semidefinite matrices.

        EXAMPLES::

            sage: from sage.geometry.semialgebraic.semidefinite import PositiveDefiniteMatrices
            sage: M_pd = PositiveDefiniteMatrices(AA, 2)
            sage: M_pd.closure()
        """
        return PositiveSemidefiniteMatrices(self._matrix_space)

    _predicate = attrcall('is_positive_definite')

    def is_closed(self):
        return self.dimension() == 0
