r"""
Morphisms of chain complexes

AUTHORS:

- Benjamin Antieau <d.ben.antieau@gmail.com> (2009.06)

- Travis Scrimshaw (2012-08-18): Made all simplicial complexes immutable to
  work with the homset cache.

This module implements morphisms of chain complexes. The input is a dictionary
whose keys are in the grading group of the chain complex and whose values are
matrix morphisms.

EXAMPLES::

    from sage.matrix.constructor import zero_matrix
    sage: S = simplicial_complexes.Sphere(1)
    sage: S
    Simplicial complex with vertex set (0, 1, 2) and facets {(1, 2), (0, 2), (0, 1)}
    sage: C = S.chain_complex()
    sage: C.differential()
    {0: [], 1: [ 1  1  0]
    [ 0 -1 -1]
    [-1  0  1], 2: []}
    sage: f = {0:zero_matrix(ZZ,3,3),1:zero_matrix(ZZ,3,3)}
    sage: G = Hom(C,C)
    sage: x = G(f)
    sage: x
    Chain complex morphism from Chain complex with at most 2 nonzero terms over Integer Ring to Chain complex with at most 2 nonzero terms over Integer Ring
    sage: x._matrix_dictionary
    {0: [0 0 0]
    [0 0 0]
    [0 0 0], 1: [0 0 0]
    [0 0 0]
    [0 0 0]}
"""

#*****************************************************************************
# Copyright (C) 2009 D. Benjamin Antieau <d.ben.antieau@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#
#*****************************************************************************

import sage.matrix.all as matrix
from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ

def is_ChainComplexMorphism(x):
    """
    Returns ``True`` if and only if ``x`` is a chain complex morphism.

    EXAMPLES::

        sage: from sage.homology.chain_complex_morphism import is_ChainComplexMorphism
        sage: S = simplicial_complexes.Sphere(14)
        sage: H = Hom(S,S)
        sage: i = H.identity()  # long time (8s on sage.math, 2011)
        sage: S = simplicial_complexes.Sphere(6)
        sage: H = Hom(S,S)
        sage: i = H.identity()
        sage: x = i.associated_chain_complex_morphism()
        sage: x # indirect doctest
        Chain complex morphism from Chain complex with at most 7 nonzero terms over
        Integer Ring to Chain complex with at most 7 nonzero terms over Integer Ring
        sage: is_ChainComplexMorphism(x)
        True
    """
    return isinstance(x,ChainComplexMorphism)

class ChainComplexMorphism(SageObject):
    """
    An element of this class is a morphism of chain complexes.
    """
    def __init__(self, matrices, C, D, check=True):
        """
        Create a morphism from a dictionary of matrices.

        EXAMPLES::

            from sage.matrix.constructor import zero_matrix
            sage: S = simplicial_complexes.Sphere(1)
            sage: S
            Simplicial complex with vertex set (0, 1, 2) and facets {(1, 2), (0, 2), (0, 1)}
            sage: C = S.chain_complex()
            sage: C.differential()
            {0: [], 1: [ 1  1  0]
            [ 0 -1 -1]
            [-1  0  1], 2: []}
            sage: f = {0:zero_matrix(ZZ,3,3),1:zero_matrix(ZZ,3,3)}
            sage: G = Hom(C,C)
            sage: x = G(f)
            sage: x
            Chain complex morphism from Chain complex with at most 2 nonzero terms
            over Integer Ring to Chain complex with at most 2 nonzero terms over 
            Integer Ring
            sage: x._matrix_dictionary
            {0: [0 0 0]
            [0 0 0]
            [0 0 0], 1: [0 0 0]
            [0 0 0]
            [0 0 0]}

        Check that the bug in :trac:`13220` has been fixed::

            sage: X = simplicial_complexes.Simplex(1)
            sage: Y = simplicial_complexes.Simplex(0)
            sage: g = Hom(X,Y)({0:0, 1:0})
            sage: g.associated_chain_complex_morphism()
            Chain complex morphism from Chain complex with at most 2 nonzero 
            terms over Integer Ring to Chain complex with at most 1 nonzero terms 
            over Integer Ring

        Check that an error is raised if the matrices are the wrong size::

            sage: C = ChainComplex({0: zero_matrix(ZZ, 0, 1)})
            sage: D = ChainComplex({0: zero_matrix(ZZ, 0, 2)})
            sage: Hom(C,D)({0: matrix(1, 2, [1, 1])})  # 1x2 is the wrong size.
            Traceback (most recent call last):
            ...
            ValueError: matrix in degree 0 is not the right size
            sage: Hom(C,D)({0: matrix(2, 1, [1, 1])})  # 2x1 is right.
            Chain complex morphism from Chain complex with at most 1 nonzero terms over Integer Ring to Chain complex with at most 1 nonzero terms over Integer Ring
        """
        if not C.base_ring() == D.base_ring():
            raise NotImplementedError('morphisms between chain complexes of different'
                                      ' base rings are not implemented')
        d = C.degree_of_differential()
        if d != D.degree_of_differential():
            raise ValueError('degree of differential does not match')
            
        from sage.misc.misc import uniq
        degrees = uniq(C.differential().keys() + D.differential().keys())
        initial_matrices = dict(matrices)
        matrices = dict()
        for i in degrees:
            if i - d not in degrees:
                if not (C.free_module_rank(i) == D.free_module_rank(i) == 0):
                    raise ValueError('{} and {} are not rank 0 in degree {}'.format(C, D, i))
                continue
            try:
                matrices[i] = initial_matrices.pop(i)
            except KeyError:
                matrices[i] = matrix.zero_matrix(C.base_ring(),
                                                 D.differential(i).ncols(),
                                                 C.differential(i).ncols(), sparse=True)
        if check:
            # All remaining matrices given must be 0x0.
            if not all(m.ncols() == m.nrows() == 0 for m in initial_matrices.values()):
                raise ValueError('the remaining matrices are not empty')
            # Check sizes of matrices.
            for i in matrices:
                if (matrices[i].nrows() != D.free_module_rank(i) or
                    matrices[i].ncols() != C.free_module_rank(i)):
                    raise ValueError('matrix in degree {} is not the right size'.format(i))
            # Check commutativity.
            for i in degrees:
                if i - d not in degrees:
                    if not (C.free_module_rank(i) == D.free_module_rank(i) == 0):
                        raise ValueError('{} and {} are not rank 0 in degree {}'.format(C, D, i))
                    continue
                if i + d not in degrees:
                    if not (C.free_module_rank(i+d) == D.free_module_rank(i+d) == 0):
                        raise ValueError('{} and {} are not rank 0 in degree {}'.format(C, D, i+d))
                    continue
                Dm = D.differential(i) * matrices[i]
                mC = matrices[i+d] * C.differential(i)
                if mC != Dm:
                    raise ValueError('matrices must define a chain complex morphism')
        self._matrix_dictionary = matrices
        self._domain = C
        self._codomain = D

    def in_degree(self, n):
        """
        The matrix representing this morphism in degree n

        INPUT:

        - ``n`` -- degree

        EXAMPLES::

            sage: C = ChainComplex({0: identity_matrix(ZZ, 1)})
            sage: D = ChainComplex({0: zero_matrix(ZZ, 1), 1: zero_matrix(ZZ, 1)})
            sage: f = Hom(C,D)({0: identity_matrix(ZZ, 1), 1: zero_matrix(ZZ, 1)})
            sage: f.in_degree(0)
            [1]

        Note that if the matrix is not specified in the definition of
        the map, it is assumed to be zero::

            sage: f.in_degree(2)
            []
            sage: f.in_degree(2).nrows(), f.in_degree(2).ncols()
            (1, 0)
            sage: C.free_module(2)
            Ambient free module of rank 0 over the principal ideal domain Integer Ring
            sage: D.free_module(2)
            Ambient free module of rank 1 over the principal ideal domain Integer Ring
        """
        try:
            return self._matrix_dictionary[n]
        except KeyError:
            from sage.matrix.constructor import zero_matrix
            rows = self._codomain.free_module_rank(n)
            cols = self._domain.free_module_rank(n)
            return zero_matrix(self._domain.base_ring(), rows, cols)

    def __neg__(self):
        """
        Returns ``-x``.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: x = i.associated_chain_complex_morphism()
            sage: w = -x
            sage: w._matrix_dictionary
            {0: [-1  0  0  0]
            [ 0 -1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1],
             1: [-1  0  0  0  0  0]
            [ 0 -1  0  0  0  0]
            [ 0  0 -1  0  0  0]
            [ 0  0  0 -1  0  0]
            [ 0  0  0  0 -1  0]
            [ 0  0  0  0  0 -1],
             2: [-1  0  0  0]
            [ 0 -1  0  0]
            [ 0  0 -1  0]
            [ 0  0  0 -1]}

        """
        f = dict()
        for i in self._matrix_dictionary.keys():
            f[i] = -self._matrix_dictionary[i]
        return ChainComplexMorphism(f, self._domain, self._codomain)

    def __add__(self,x):
        """
        Returns ``self + x``.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: x = i.associated_chain_complex_morphism()
            sage: z = x+x
            sage: z._matrix_dictionary
            {0: [2 0 0 0]
            [0 2 0 0]
            [0 0 2 0]
            [0 0 0 2],
             1: [2 0 0 0 0 0]
            [0 2 0 0 0 0]
            [0 0 2 0 0 0]
            [0 0 0 2 0 0]
            [0 0 0 0 2 0]
            [0 0 0 0 0 2],
             2: [2 0 0 0]
            [0 2 0 0]
            [0 0 2 0]
            [0 0 0 2]}

        """
        if not isinstance(x,ChainComplexMorphism) or self._codomain != x._codomain or self._domain != x._domain or self._matrix_dictionary.keys() != x._matrix_dictionary.keys():
            raise TypeError("Unsupported operation.")
        f = dict()
        for i in self._matrix_dictionary.keys():
            f[i] = self._matrix_dictionary[i] + x._matrix_dictionary[i]
        return ChainComplexMorphism(f, self._domain, self._codomain)

    def __mul__(self,x):
        """
        Return ``self * x`` if ``self`` and ``x`` are composable morphisms
        or if ``x`` is an element of the base ring.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: x = i.associated_chain_complex_morphism()
            sage: y = x*2
            sage: y._matrix_dictionary
            {0: [2 0 0 0]
            [0 2 0 0]
            [0 0 2 0]
            [0 0 0 2],
             1: [2 0 0 0 0 0]
            [0 2 0 0 0 0]
            [0 0 2 0 0 0]
            [0 0 0 2 0 0]
            [0 0 0 0 2 0]
            [0 0 0 0 0 2],
             2: [2 0 0 0]
            [0 2 0 0]
            [0 0 2 0]
            [0 0 0 2]}
            sage: z = y*y
            sage: z._matrix_dictionary
            {0: [4 0 0 0]
            [0 4 0 0]
            [0 0 4 0]
            [0 0 0 4],
             1: [4 0 0 0 0 0]
            [0 4 0 0 0 0]
            [0 0 4 0 0 0]
            [0 0 0 4 0 0]
            [0 0 0 0 4 0]
            [0 0 0 0 0 4],
             2: [4 0 0 0]
            [0 4 0 0]
            [0 0 4 0]
            [0 0 0 4]}

        TESTS:

        Make sure that the product is taken in the correct order
        (``self * x``, not ``x * self`` -- see :trac:`19065`)::

            sage: C = ChainComplex({0: zero_matrix(ZZ, 0, 2)})
            sage: D = ChainComplex({0: zero_matrix(ZZ, 0, 1)})
            sage: f = Hom(C,D)({0: matrix(1, 2, [1, 1])})
            sage: g = Hom(D,C)({0: matrix(2, 1, [1, 1])})
            sage: (f*g).in_degree(0)
            [2]

        Before :trac:`19065`, the following multiplication produced a
        ``KeyError`` because `f` was not explicitly defined in degree 2::

            sage: C0 = ChainComplex({0: zero_matrix(ZZ, 0, 1)})
            sage: C1 = ChainComplex({1: zero_matrix(ZZ, 0, 1)})
            sage: C2 = ChainComplex({2: zero_matrix(ZZ, 0, 1)})
            sage: f = ChainComplexMorphism({}, C0, C1)
            sage: g = ChainComplexMorphism({}, C1, C2)
            sage: g * f
            Chain complex morphism from Chain complex with at most 1 nonzero terms over Integer Ring to Chain complex with at most 1 nonzero terms over Integer Ring
            sage: f._matrix_dictionary
            {0: [], 1: []}
            sage: g._matrix_dictionary
            {1: [], 2: []}
        """
        if not isinstance(x,ChainComplexMorphism) or self._domain != x._codomain:
            try:
                y = self._domain.base_ring()(x)
            except TypeError:
                raise TypeError("multiplication is not defined")
            f = dict()
            for i in self._matrix_dictionary:
                f[i] = self._matrix_dictionary[i] * y
            return ChainComplexMorphism(f,self._domain,self._codomain)
        f = dict()
        for i in self._matrix_dictionary:
            f[i] = self._matrix_dictionary[i]*x.in_degree(i)
        return ChainComplexMorphism(f,x._domain,self._codomain)

    def __rmul__(self,x):
        """
        Returns ``x * self`` if ``x`` is an element of the base ring.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: x = i.associated_chain_complex_morphism()
            sage: 2*x == x*2
            True
            sage: 3*x == x*2
            False
        """
        try:
            y = self._domain.base_ring()(x)
        except TypeError:
            raise TypeError("multiplication is not defined")
        f = dict()
        for i in self._matrix_dictionary.keys():
            f[i] = y * self._matrix_dictionary[i]
        return ChainComplexMorphism(f,self._domain,self._codomain)

    def __sub__(self,x):
        """
        Return ``self - x``.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: x = i.associated_chain_complex_morphism()
            sage: y = x-x
            sage: y._matrix_dictionary
            {0: [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0],
             1: [0 0 0 0 0 0]
            [0 0 0 0 0 0]
            [0 0 0 0 0 0]
            [0 0 0 0 0 0]
            [0 0 0 0 0 0]
            [0 0 0 0 0 0],
             2: [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]}

        """
        return self + (-x)

    def __eq__(self,x):
        """
        Return ``True`` if and only if ``self == x``.

        EXAMPLES::

            sage: S = SimplicialComplex(is_mutable=False)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: x = i.associated_chain_complex_morphism()
            sage: x
            Chain complex morphism from Trivial chain complex over Integer Ring
            to Trivial chain complex over Integer Ring
            sage: f = x._matrix_dictionary
            sage: C = S.chain_complex()
            sage: G = Hom(C,C)
            sage: y = G(f)
            sage: x == y
            True
        """
        return isinstance(x,ChainComplexMorphism) \
                and self._codomain == x._codomain \
                and self._domain == x._domain \
                and self._matrix_dictionary == x._matrix_dictionary

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: S = SimplicialComplex(is_mutable=False)
            sage: H = Hom(S,S)
            sage: i = H.identity()
            sage: x = i.associated_chain_complex_morphism()
            sage: x
            Chain complex morphism from Trivial chain complex over Integer Ring
            to Trivial chain complex over Integer Ring
            sage: x._repr_()
            'Chain complex morphism from Trivial chain complex over Integer Ring
            to Trivial chain complex over Integer Ring'
        """
        return "Chain complex morphism from {} to {}".format(self._domain, self._codomain)

