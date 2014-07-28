# cython: profile = True
r"""
Implements helper routines for checking existance of a represented minor

A represented minor of a `GF(q)-` representable matroid ``M`` is also a
`GF(q)-` representable matroid. Its representation matrix can be obtained
from the representation matrix of ``M`` corresponding to some basis by
deleting rows (contraction) and deleting columns (deletion) and by permuting
the leftover rows and columns. We are given a `GF(q)-` representable matroid
``M``, another `GF(q)` representable matroid ``N`` , a basis ``B_M`` of
``M`` and basis ``B_N`` of ``N``. We are looking for injective but not
necessarily surjective maps from the ``B_N`` to ``B_M`` and from ``E(N)-B_N``
to ``E(M)-B_M``. This extension is designed to find such maps (provided that
they exist) by reducing the problem to that of induced sungraph isomorphism
testing. Essentially we would test whether the fundamental graph (see [Oxley]
pg. 185) of ``N`` corresponding to basis ``B_N`` can be obtained as a vertex
induced subgraph of the fundamental graph of ``M`` corresponding to basis
``B_M``. The algorithm used is the one in [Ullman] which employs a depth
first search like procedure with pruning. While the functions in this
extension are meant to be used by various childern of LinearMatroid class, one
can directly access them using:

    sage: from sage.matroids.advanced import *

AUTHORS:

- Jayant Apte (2005-01-03): initial version

EXAMPLES::


REFERENCES
==========
..  [Ullman] J. R. Ullmann. 1976. An Algorithm for Subgraph Isomorphism.
    J. ACM 23, 1 (January 1976), 31-42.

"""

# *****************************************************************************
#       Copyright (C) 2013 Jayant Apte <jayant91089@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
include 'sage/misc/bitset.pxi'
from lean_matrix cimport BinaryMatrix
from sage.rings.all import GF
cdef bint npruned_not_defined = True
cdef npruned

cpdef init_iso_matrices(BinaryMatrix M_rmat, BinaryMatrix N_rmat):
    """
    Return a list of two matrices [M1_0, M2_0] s.t. [M1_0 0; 0 M2_0] will be
    the bitwise logical OR of all degree preserving isomorphisms

    INPUT:

    - ``M_rmat`` -- reduced representation matrix of type ``BinaryMatrix``
       of a binary matroid.
    - ``N_rmat`` -- reduced representation matrix of type ``BinaryMatrix``
       of a binary matroid.

    OUTPUT:

    A list [M1_0, M2_0] of 2 ``BinaryMatrix`` objects.

    ..NOTE::

        The dimensions of ``N_rmat`` are assumed to be less that equal to
        those of ``M_rmat``.

    .. WARNING::

            Intended for internal use. This method does no checks of any kind.
    """
    cdef long r1, r2
    cdef BinaryMatrix M1_0, M2_0
    # initialize M1_0 and M2_0
    M1_0 = BinaryMatrix(N_rmat.nrows(), M_rmat.nrows())
    M2_0 = BinaryMatrix(N_rmat.ncols(), M_rmat.ncols())
    for r1 in xrange(N_rmat.nrows()):
        d1 = N_rmat.row_len(r1)
        for r2 in xrange(M_rmat.nrows()):
            d2 = M_rmat.row_len(r2)
            if d1 <= d2:
                M1_0.set(r1, r2)
    for c1 in xrange(N_rmat.ncols()):
        d1 = col_len(N_rmat, c1)
        for c2 in xrange(M_rmat.ncols()):
            d2 = col_len(M_rmat, c2)
            if d1 <= d2:
                M2_0.set(c1, c2)
    return [M1_0, M2_0]

cpdef col_len(BinaryMatrix mat, long c):
    """
    Return hamming weight of ``c`` th column of ``mat``.

    INPUT:

    - ``mat`` -- A ``BinaryMatrix`` object.
    - ``c`` -- Index of a column column of ``mat``.

    OUTPUT:

    An integer specfying hamming weight of ``c`` th column of ``mat``.

    .. WARNING::

            Intended for internal use. This method does no checks of any kind.
    """
    cdef int r, d
    d = 0
    for r in xrange(mat.nrows()):
        if mat.is_nonzero(r, c):
            d += 1
    return d

cpdef _check_bin_minor(BinaryMatrix M_rmat, BinaryMatrix N_rmat, indices, M,
                       Npcl, long nloops, bint debug=False):
    """
    Return ``True`` if the bipartite graph corresponding to ``N_rmat`` is
    vertex induced subgraph of bipartite graph corresponding to ``M_rmat``

    INPUT:

    - ``M_rmat`` -- reduced representation matrix of type ``BinaryMatrix``
       of a simple binary matroid.
    - ``N_rmat`` -- reduced representation matrix of type ``BinaryMatrix``
       of a simple binary matroid.
    - ``indices`` -- A list of 5 lists that specify labels of rows and columns
       of ``M_rmat``,``N_rmat`` and rows of reduced representation matrix of
       ``M`` respectively.
    - ``M`` -- A ``BinaryMatroid`` such that ``M_rmat`` is reduced
       representation matrix of the associated simple matroid.
    - ``nloops`` -- Number of loops in matroid ``M``.
    - ``debug`` -- (default: ``False``) Placeholder for a debugging flag.

    OUTPUT:

    ``True`` if the bipartite graph corresponding to ``N_rmat`` is vertex
    induced subgraph of bipartite graph corresponding to ``M_rmat``.

    .. WARNING::

            Intended for internal use. This method does no checks of any kind.
    """
    cdef long i, x
    cdef bitset_t unused_cols
    cdef BinaryMatrix M_rmatT, N_rmatT, M1, M2
    M_rmatT = mat_transpose(M_rmat)
    N_rmatT = mat_transpose(N_rmat)
    bitset_init(unused_cols, M_rmat.nrows() + M_rmat.ncols())
    [M1_0, M2_0] = init_iso_matrices(M_rmat, N_rmat)
    for i in xrange(M_rmat.nrows()+M_rmat.ncols()):
        bitset_add(unused_cols, i)
    M1 = copy_mat(M1_0)
    M2 = BinaryMatrix(M2_0.nrows(), M1_0.ncols()).augment(copy_mat(M2_0))
    if not degrees_are_sane(M1, M2):
        return False
    else:
        x = recurse(unused_cols, 0, M_rmat, N_rmat, M_rmatT, N_rmatT, M1, M2,
                    indices, M, Npcl, nloops)
        return x == 1

cpdef degrees_are_sane(BinaryMatrix M1, BinaryMatrix M2):
    """
    Return ``True`` if there is an all-zero row in M1 or M2 and return
    ``False`` otherwise.

    INPUT:

    - ``M1`` -- A ``BinaryMatrix`` object
    - ``M2`` -- A ``BinaryMatrix`` object

    OUTPUT:

    A ``bool`` that is ``True`` if neither ``M1`` nor ``M2`` have an all-zero
    row and ``False`` otherwise.

    .. WARNING::

            Intended for internal use; does no input checking.

    """
    cdef long i
    for i in xrange(M1.nrows()):
        if bitset_isempty(M1._M[i]):
            return False
    for i in xrange(M2.nrows()):
        if bitset_isempty(M2._M[i]):
            return False
    return True


cdef _pcl_have_maps(M, Npcl, iso_map, indices, long nloops):
    """
    Return ``True`` if parallel elememts and loops of ``M`` have consistent
    mapping under ``iso_map`` and ``False`` otherwise.

    INPUT:

    - ``M`` -- A ``BinaryMatroid``.
    - ``Npcl`` -- A list of lists specifying parallel classes of a candidate
      minor of ``M``.
    - ``iso_map`` -- A dictionary specifying an injective but not necessarily
      surjective map from ground set of candidate minor ``N`` to that of
      simple matroid associated with ``M``.
    - ``nloops`` -- Number of loops in candidate minor ``N``.
    - ``indices`` --A list of 5 lists such that ``indices[4]`` specifies the
      row labels of reduced representation matrix of ``M``

    OUTPUT:

    A ``bool`` that is ``True`` if parallel elements and loops of ``M`` have
    consistent mapping under ``iso_map`` and ``False`` otherwise.

    .. WARNING::

            Intended for internal use; does no input checking.

    """
    # contract the rows that are unmapped
    M2 = M.contract(set(indices[4]) - set(iso_map.values()))
    if len(M2.loops()) < nloops:
        return False
    # for each repmap check if it has big enough parallel class
    for pcl in Npcl:
        rep = list(set(pcl) & set(iso_map.keys()))[0]  # representative
        repmap = iso_map[rep]
        repmap_pcl_M = set(M2.closure({repmap})) - M2.loops()
        if len(repmap_pcl_M) < len(pcl):
            # pcl doesn't have consistant map
            return False
    return True

cdef iso_mats_to_dict(BinaryMatrix M1, BinaryMatrix M2, indices):
    """
    Return a dictionary corresponding to maps specified by ``M1`` and ``M2``.

    INPUT:

    - ``M1`` -- A ``BinaryMatrix`` object.
    - ``M2`` -- A ``BinaryMatrix`` object.
    - ``indices`` -- A list of lists specifying labels of rows and columns of
      ``M1`` and ``M2`` respectively. Lists indices[4] onwards are ignored.

    OUTPUT:

    A dictionary mapping labels of rows and columns of ``M1`` to those of `M2``
    respectively.

    .. WARNING::

            Intended for internal use; does no input checking.

    """
    M_R = indices[0]
    M_C = indices[1]
    N_R = indices[2]
    N_C = indices[3]
    iso_map = {}
    for i in xrange(M1.nrows()):
        iso_map[N_R[i]] = M_R[bitset_list(M1._M[i])[0]]
    for i in xrange(M2.nrows()):
        iso_map[N_C[i]] = M_C[bitset_list(M2._M[i])[0] - M1.ncols()]
    return iso_map


cdef recurse(bitset_t unused_cols, long cur_row, BinaryMatrix M_rmat,
             BinaryMatrix N_rmat, BinaryMatrix M_rmatT, BinaryMatrix N_rmatT,
             BinaryMatrix M1, BinaryMatrix M2, indices, M, Npcl, long nloops):
    """
    Return ``True`` if there exists an induced subgraph isomorphism and
    ``False`` otherwise.

    INPUT:

    - ``unused_cols`` -- A ``bitset_t`` object specifying unmapped vertices of
      bipartite graph corresponding to ``M_rmat``.
    - ``cur_row`` -- Depth of recursion.
    - ``M_rmat`` -- A ``BinaryMatrix`` that is adjacency matrix of a bipartite
      graph.
    - ``N_rmat`` -- A ``BinaryMatrix`` that is adjacency matrix of a bipartite
      graph.
    - ``M_rmatT`` -- Transpose of ``M_rmat``.
    - ``N_rmatT`` -- Transpose of ``N_rmat``.
    - ``M1,M2`` -- Non-zero submatrices of isomorphism matrix.
      ``[[M1 0],[0 M2]]`` between two bipartite graphs.
    - ``indices`` -- ``indices`` -- A list of 5 lists that specify labels of
      rows and columns of ``M_rmat``,``N_rmat`` and rows of reduced
      representation matrix of ``M`` respectively.
    - ``M`` -- A ``BinaryMatroid``.
    - ``Npcl`` -- A list of lists specifying parallel classes of a candidate
      minor of ``M``.
    - ``nloops`` -- Number of loops in candidate minor ``N``

    OUTPUT:

    A ``bool`` that is ``True`` if there is an induced subgrph isomorphism
    between fundamental graphs ``N_rmat`` and ``M_rmat`` of simple matroids
    such that the ``nloops`` loops and ``len(Npcl)`` parallel elements of the
    non-simple candidate minor have consistent mappings.

    .. WARNING::

            Intended for internal use; does no input checking.

    """
    cdef long i, j
    cpdef bitset_t valid_cols
    cpdef bitset_t curr_cols
    cpdef bitset_t saved_row
    cdef BinaryMatrix iso, saved_M1, saved_M2
    cdef long found = 0
    if cur_row == M1.nrows() + M2.nrows():
        found = is_weak_induced_isomorphism(M1, M2, M_rmat, N_rmat)
        if found == 1:
            iso_map = iso_mats_to_dict(M1, M2, indices)
            if _pcl_have_maps(M, Npcl, iso_map, indices, nloops):
                return found
            else:
                return 0
        else:
            return 0
    else:
        # save M1, M2 before pruning
        saved_M1 = copy_mat(M1)
        saved_M2 = copy_mat(M2)
        if prune(cur_row, M1, M2, M_rmat, N_rmat, M_rmatT, N_rmatT) is False:
            if cur_row < M1.nrows():
                restore_mat(M1, saved_M1, cur_row)
                restore_mat(M2, saved_M2, 0)
            else:
                restore_mat(M2, saved_M2, cur_row - M1.nrows())
            return found
        if cur_row < M1.nrows():
            # cur_row in in M1
            bitset_init(curr_cols, M_rmat.nrows() + M_rmat.ncols())
            bitset_init(valid_cols, M_rmat.nrows() + M_rmat.ncols())
            for i in xrange(M1.ncols()):
                bitset_add(curr_cols, i)
            bitset_intersection(valid_cols, curr_cols, M1._M[cur_row])
            bitset_intersection(valid_cols, valid_cols, unused_cols)
            bitset_init(saved_row, M1.ncols())
            bitset_copy(saved_row, M1._M[cur_row])
        else:
            # cur_row is in M2
            bitset_init(curr_cols, M_rmat.nrows()+M_rmat.ncols())
            bitset_init(valid_cols, M_rmat.nrows() + M_rmat.ncols())
            for i in xrange(M1.ncols(), M1.ncols() + M2.ncols()):
                bitset_add(curr_cols, i)
            bitset_intersection(valid_cols, curr_cols,
                                M2._M[cur_row - M1.nrows()])
            bitset_intersection(valid_cols, valid_cols, unused_cols)
            bitset_init(saved_row, M2.ncols())
            bitset_copy(saved_row, M2._M[cur_row - M1.nrows()])
    while not bitset_isempty(valid_cols):
        i = bitset_pop(valid_cols)
        bitset_remove(unused_cols, i)
        if cur_row < M1.nrows():
            # set up M1 and recurse
            bitset_clear(M1._M[cur_row])
            bitset_add(M1._M[cur_row], i)
            found = recurse(unused_cols, cur_row + 1, M_rmat, N_rmat, M_rmatT,
                            N_rmatT, M1, M2, indices, M, Npcl, nloops)
            bitset_add(unused_cols, i)
            if found == 1:
                break
        else:
            # set up M2 and recurse
            bitset_clear(M2._M[cur_row - M1.nrows()])
            bitset_add(M2._M[cur_row - M1.nrows()], i)
            found = recurse(unused_cols, cur_row + 1, M_rmat, N_rmat, M_rmatT,
                            N_rmatT, M1, M2, indices, M, Npcl, nloops)
            bitset_add(unused_cols, i)
            if found == 1:
                break
    if cur_row < M1.nrows():
        restore_mat(M1, saved_M1, cur_row)
        restore_mat(M2, saved_M2, 0)
    else:
        restore_mat(M2, saved_M2, cur_row - M1.nrows())
    return found

cdef restore_mat(BinaryMatrix mat, BinaryMatrix saved_mat, long cur_row):
    """
    Copy ``cur_row`` onwards rows of ``saved mat`` into ``mat``

    INPUT:

    - ``mat`` -- A ``BinaryMatrix``.
    - ``saved_mat`` -- A ``BinaryMatrix``.
    - ``cur_row`` -- Row index beyond which the ``mat`` will be restored.

    .. NOTE::

        Dimensions of ``mat`` are assumed to be same as that of ``saved_mat``
        and `0\leq` ``cur_row`` `\leq` ``mat.nrows()`` `-1`

    .. WARNING::

            Intended for internal use; does no input checking.
    """
    cdef long r
    for r in xrange(cur_row, saved_mat.nrows()):
        bitset_copy(mat._M[r], saved_mat._M[r])
    return

cpdef prune(long cur_row, BinaryMatrix M1, BinaryMatrix M2,
            BinaryMatrix M_rmat, BinaryMatrix N_rmat, BinaryMatrix M_rmatT,
            BinaryMatrix N_rmatT):
    """
    Prune the search tree based on consistency of mapping of neighbours

    INPUT:

    - ``cur_row`` -- depth of recursion.
    - ``M1`` -- A ``BinaryMatrix``.
    - ``M2`` -- A ``BinaryMatrix``.
    - ``M_rmat`` -- A ``BinaryMatrix`` that is adjacency matrix of a bipartite
      graph.
    - ``N_rmat`` -- A ``BinaryMatrix`` that is adjacency matrix of a bipartite
      graph.
    - ``M_rmatT`` -- Transpose of ``M_rmat``.
    - ``N_rmatT`` -- Transpose of ``N_rmat``.

    OUTPUT:

    A ``bool`` that is ``True`` if pruning doesn't create any all-zero rows
    and ``False`` otherwise.

    .. WARNING::

            Intended for internal use; does no input checking.
    """
    global npruned, npruned_not_defined
    if npruned_not_defined:
        npruned = 0
        npruned_not_defined = False
    cdef long r, c, flag
    for r in xrange(max(0, cur_row), M1.nrows()):
        for c in xrange(M1.ncols()):
            if M1.is_nonzero(r, c):
                for x in _neighbours(N_rmat, N_rmatT, r, 1):
                    flag = 0
                    for y in _neighbours(M_rmat, M_rmatT, c, 1):
                        if M2.is_nonzero(x, y + M1.ncols()):
                            flag = 1
                            break
                    if flag == 0:
                        npruned += 1
                        M1.set_unsafe(r, c, 0)
    for r in xrange(max(0, cur_row), M2.nrows()):
        for c in xrange(M1.ncols(), M2.ncols()):
            if M2.is_nonzero(r, c):
                for x in _neighbours(N_rmat, N_rmatT, r, 0):
                    flag = 0
                    for y in _neighbours(M_rmat, M_rmatT, c - M1.ncols(), 0):
                        if M1.is_nonzero(x, y):
                            flag = 1
                            break
                    if flag == 0:
                        npruned += 1
                        M2.set_unsafe(r, c, 0)
    # check degree sanity
    if degrees_are_sane(M1, M2):
        return True
    else:
        return False

cpdef BinaryMatrix copy_mat(BinaryMatrix mat):
    """
    Return a copy of ``mat``

    INPUT:

    - ``mat`` -- A ``BinaryMatrix``.

    OUTPUT:

    A ``BinaryMatrix`` that is copy of ``mat``

    .. WARNING::

            Intended for internal use; does no input checking.
    """
    cdef long i
    cdef BinaryMatrix mat_copy
    mat_copy = BinaryMatrix(mat.nrows(), mat.ncols())
    for i in xrange(mat.nrows()):
        bitset_copy(mat_copy._M[i], mat._M[i])
    return mat_copy

cpdef _neighbours(BinaryMatrix rmat, BinaryMatrix rmatT, long i, long rc):
    """
    Return a list of neighbours of vertex with index ``i``

    INPUT:

    - ``rmat`` -- A ``BinaryMatrix`` that is adjacency matrix of a bipartite
      graph.
    - ``rmatT`` -- Transpose of ``rmat``.
    - ``i`` -- Index of row or column of ``rmat``
    - ``rc`` -- A binary number specifying whether ``i`` is index of a row or
      a column of ``rmat``.

    OUTPUT:

    A list of indices of opsitions that are non-zero in ``i`` th row or
    column.

    .. WARNING::

            Intended for internal use; does no input checking.
    """
    cdef long c
    neighbours = []
    if rc == 1:  # neighbors of lhs , look through ith row
        for c in xrange(rmat.ncols()):
            if rmat.is_nonzero(i, c):
                neighbours.append(c)
    if rc == 0:
        for c in xrange(rmatT.ncols()):
            if rmatT.is_nonzero(i, c):
                neighbours.append(c)
    return neighbours

cpdef is_weak_induced_isomorphism(BinaryMatrix M1_1, BinaryMatrix M2_1,
                                  BinaryMatrix  M_rmat, BinaryMatrix N_rmat):
    """
    Return ``True`` if ``M1_1`` and ``M2_1`` specify an induced graph
    isomorphism

    INPUT:

    - ``M1_1, M2_1`` -- ``BinaryMatrix`` objects. ``[[M1_1, 0],[0, M2_1]]`` is
     a candidate vertex induced bipartite isomorphism
    - ``M_rmat`` -- Adjacency matrix corresponding to a bipartite graph
    - ``N_rmat`` -- Adjacency matrix corresponding to a bipartite graph

    OUTPUT:

    A ``bool`` that is ``True`` if ``[[M1_1, 0],[0, M2_1]]`` is a vertex
    induced isomorphism.

    ALGORITHM:

    Adapted from equation (1) in [Ullman] with double implication to test for
    induced subgraph isomorphism.

    .. WARNING::

            Intended for internal use; does no input checking.
    """
    cdef long i, k
    cdef int j
    cdef int * roworders
    cdef int * colorders
    cdef BinaryMatrix rmat1
    cdef BinaryMatrix M1, M2
    M1 = M1_1.copy()
    M2 = M2_1.copy()
    # rmat corresponding to given isomorphism
    rmat1 = BinaryMatrix(N_rmat.nrows(), N_rmat.ncols())
    roworders = <int * > sage_malloc(N_rmat.nrows() * sizeof(int))
    for i in xrange(M1.nrows()):
        j = bitset_pop(M1._M[i])
        # row i maps to row j
        roworders[i] = j
    colorders = <int * > sage_malloc(N_rmat.ncols() * sizeof(int))
    for i in xrange(M2.nrows()):
        j = bitset_pop(M2._M[i])
        # col i maps to row j
        colorders[i] = j - M1.ncols()
    for i in xrange(M1.nrows()):
        for k in xrange(M2.nrows()):
            rmat1.set_unsafe(i, k, M_rmat.get_unsafe(roworders[i],
                                                     colorders[k]))
    sage_free(roworders)
    sage_free(colorders)
    for i in xrange(N_rmat.nrows()):
        if not bitset_eq(N_rmat._M[i], rmat1._M[i]):
            return 0
    return 1

cdef mat_transpose(BinaryMatrix mat):
    """
    Return transpose of matrix ``mat``

    INPUT:

    - ``mat`` -- A ``BinaryMatrix``.

    OUTPUT:

    A ``BinaryMatrix`` that is transpose of ``mat``

    .. WARNING::

            Intended for internal use; does no input checking.
    """
    cdef long i, j
    matT = BinaryMatrix(mat.ncols(), mat.nrows())
    for i in xrange(mat.nrows()):
        for j in xrange(mat.ncols()):
            matT.set_unsafe(j, i, mat.get_unsafe(i, j))
    return matT

cpdef bint is_new_rmat(BinaryMatrix cand, used_rmats):
    """
    Return ``True`` if matrix ``cand`` is present in list ``used_rmats``,
    ``False`` otherwise.

    INPUT:

    - ``cand`` -- A ``BinaryMatrix``.
    -``used_rmats`` -- A list of ``BinaryMatrix`` objects

    OUTPUT:

    A ``bint`` that is ``True`` if cand is not present in the list
    ``used_rmats``.

    .. WARNING::

            Intended for internal use; does no input checking.
    """
    cdef long i
    for i in xrange(len(used_rmats)):
        if mats_equal(cand, used_rmats[i]):
            return False
    return True

cpdef bint mats_equal(BinaryMatrix M1, BinaryMatrix M2):
    """
    Return ``True`` if ``M1`` is elementwise equal to ``M2``

    INPUT:

    - ``M1`` -- A ``BinaryMatrix``.
    - ``M2`` -- A ``BinaryMatrix``.

    OUTPUT:

    A ``bint`` that is ``True`` if ``M1=M2`` and ``False`` otherwise.

    .. WARNING::

            Intended for internal use; does no input checking.
    """
    cdef long i
    for i in xrange(M1.nrows()):
        if not bitset_eq(M1._M[i], M2._M[i]):
            return False
    return True
