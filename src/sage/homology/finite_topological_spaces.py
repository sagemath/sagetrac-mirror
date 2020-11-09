r"""
Finite topological spaces

This module implements finite topological spaces and related concepts.

A *finite topological space* is a topological space with finitely many points and
a *finite preordered set* is a finite set with a transitive and reflexive relation.
Finite spaces and finite preordered sets are basically the same objects considered
from different perspectives. Given a finite topological space `X`, for every point
`x\in X`, define the *minimal open set* `U_x` as the intersection of all the open
sets which contain `x` (it is an open set since arbitrary intersections of open
sets in finite spaces are open). The minimal open sets constitute a basis for the
topology of `X`. Indeed, any open set `U` of `X` is the union of the sets `U_x`
with `x\in U`. This basis is called the *minimal basis of* `X`. A preorder on `X`
is given  by `x\leqslant y` if `x\in U_y`.

If `X` is now a finite preordered set, one can define a topology on `X` given by
the basis `\lbrace y\in X\vert y\leqslant x\rbrace_{x\in X}`. Note that if `y\leqslant x`,
then `y` is contained in every basic set containing `x`, and therefore `y\in U_x`.
Conversely, if `y\in U_x`, then `y\in\lbrace z\in X\vert z\leqslant x\rbrace`.
Therefore `y\leqslant x` if and only if `y\in U_x`. This shows that these two
applications, relating topologies and preorders on a finite set, are mutually
inverse. This simple remark, made in first place by Alexandroff [Ale1937]_, allows
us to study finite spaces by combining Algebraic Topology with the combinatorics
arising from their intrinsic preorder structures. The antisymmetry of a finite
preorder corresponds exactly to the `T_0` separation axiom. Recall that a topological
space `X` is said to be `T_0` if for any pair of points in `X` there exists an
open set containing one and only one of them. Therefore finite `T_0`-spaces are
in correspondence with finite partially ordered sets (posets) [Bar2011]_.

Now, if `X = \lbrace x_1, x_2, \ldots , x_n\rbrace` is a finite space and for
each `i` the unique minimal open set containing `x_i` is denoted by `U_i`, a
*topogenous matrix* of the space is the `n \times n` matrix `A = \left[a_{ij}\right]`
defined by `a_{ij} = 1` if `x_i \in U_j` and `a_{ij} = 0` otherwise (this is the
transposed matrix of the Definition 1 in [Shi1968]_). A finite space `X` is `T_0`
if and only if the topogenous matrix `A` defined above is similar (via a permutation
matrix) to a certain upper triangular matrix [Shi1968]_. This is the reason one
can assume that the topogenous matrix of a finite `T_0`-space is upper triangular.


AUTHOR::

- Julian Cuevas-Rozo (2020): Initial version

REFERENCES:

- [Ale1937]_
- [Bar2011]_
- [Shi1968]_

"""
# ****************************************************************************
#       Copyright (C) 2020 Julian Cuevas-Rozo <jlcrozo@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.parent import Parent
from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_sparse import Matrix_integer_sparse
from sage.combinat.posets.posets import Poset

from sage.rings.integer_ring import ZZ
from sage.homology.homology_group import HomologyGroup

from sage.libs.ecl import EclObject, ecl_eval, EclListIterator
from sage.interfaces import kenzo
from sage.features.kenzo import Kenzo

###############################################################
# This section will be included to src/sage/interfaces/kenzo.py

kenzonames = ['2h-regularization',
               'copier-matrice',
               'creer-matrice',
               'convertarray',
               'dvfield-aux',
               'edges-to-matrice',
               'h-regular-dif',
               'h-regular-dif-dvf-aux',
               'matrice-to-lmtrx',
               'mtrx-prdc',
               'newsmith-equal-matrix',
               'newsmith-mtrx-prdc',
               'random-top-2space',
               'randomtop',
               'vector-to-list']
               
if Kenzo().is_present():
    ecl_eval("(require :kenzo)")
    ecl_eval("(in-package :cat)")
    ecl_eval("(setf *HOMOLOGY-VERBOSE* nil)")
    for s in kenzonames:
        name = '__{}__'.format(s.replace('-', '_'))
        exec('{} = EclObject("{}")'.format(name, s))


def quotient_group_matrices(*matrices, left_null=False, right_null=False, check=True):
    r"""
    Return a presentation of the homology group `\ker M1/ \im M2`.

    INPUT:

    - ``matrices`` -- A tuple of ECL matrices. The length `L` of this parameter
      can take the value 0, 1 or 2.

    - ``left_null`` -- (default ``False``) A boolean.

    - ``right_null`` -- (default ``False``) A boolean.

    - ``check`` -- (default ``True``) A boolean. If it is ``True`` and `L=2`, it
      checks that the product of the ``matrices`` is the zero matrix.

    OUTPUT:

    - If `L=0`, it returns the trivial group.

    - If `L=1` (``matrices`` = M), then one of the parameters ``left_null`` or
      ``right_null`` must be ``True``: in case ``left_null`` == ``True``, it
      returns the homology group `\ker 0/ \im M` and in case ``right_null`` == ``True``,
      it returns the homology group `\ker M/ \im 0`.

    - If `L=2` (``matrices`` = (M1, M2)), it returns the homology group `\ker M1/ \im M2`.

    EXAMPLES::

        sage: from sage.homology.finite_topological_spaces import quotient_group_matrices, __convertarray__
        sage: from sage.interfaces.kenzo import s2k_matrix
        sage: quotient_group_matrices()
        0
        sage: s_M1 = matrix(2, 3, [1, 2, 3, 4, 5, 6])
        sage: M1 = __convertarray__(s2k_matrix(s_M1))
        sage: quotient_group_matrices(M1, left_null=True)
        C3
        sage: quotient_group_matrices(M1, right_null=True)
        Z
        sage: s_M2 = matrix(2, 2, [1, -1, 1, -1])
        sage: M2 = __convertarray__(s2k_matrix(s_M2))
        sage: s_M3 = matrix(2, 2, [1, 0, 1, 0])
        sage: M3 = __convertarray__(s2k_matrix(s_M3))
        sage: quotient_group_matrices(M2, M3)
        0
        sage: s_M4 = matrix(2, 2, [0, 0, 1, 0])
        sage: M4 = __convertarray__(s2k_matrix(s_M4))
        sage: quotient_group_matrices(M2, M4)
        Traceback (most recent call last):
        ...
        AssertionError: m1*m2 must be zero
    """
    assert not (left_null and right_null), "left_null and right_null must not be both True"
    if len(matrices)==0:
        return HomologyGroup(0, ZZ)
    elif len(matrices)==1:
        if left_null==True:
            m2 = matrices[0]
            m1 = __creer_matrice__(0, kenzo.__nlig__(m2))
        elif right_null==True:
            m1 = matrices[0]
            m2 = __creer_matrice__(kenzo.__ncol__(m1), 0)
        else:
            raise AssertionError("left_null or right_null must be True")
    elif len(matrices)==2:
        m1, m2 = matrices
        if check==True:
            rowsm1 = kenzo.__nlig__(m1)
            colsm1 = kenzo.__ncol__(m1)
            rowsm2 = kenzo.__nlig__(m2)
            colsm2 = kenzo.__ncol__(m2)
            assert colsm1==rowsm2, "Number of columns of m1 must be equal to the number of rows of m2"
            assert __newsmith_equal_matrix__(__newsmith_mtrx_prdc__(m1, m2), \
                                                   __creer_matrice__(rowsm1, colsm2)).python(), \
                                               "m1*m2 must be zero"
    homology = kenzo.__homologie__(__copier_matrice__(m1), __copier_matrice__(m2))
    lhomomology = [i for i in EclListIterator(homology)]
    res = []
    for component in lhomomology:
        pair = [i for i in EclListIterator(component)]
        res.append(pair[0].python())
    return HomologyGroup(len(res), ZZ, res)

def k2s_binary_matrix_sparse(kmatrix):
    r"""
    Converts a Kenzo binary sparse matrice (type `matrice`) to a matrix in SageMath.

    INPUT:

    - ``kmatrix`` -- A Kenzo binary sparse matrice (type `matrice`).

    EXAMPLES::

        sage: from sage.homology.finite_topological_spaces import k2s_binary_matrix_sparse, \
              s2k_binary_matrix_sparse, __randomtop__
        sage: KM2 = __randomtop__(6,1)
        sage: k2s_binary_matrix_sparse(KM2)
        [1 1 1 1 1 1]
        [0 1 1 1 1 1]
        [0 0 1 1 1 1]
        [0 0 0 1 1 1]
        [0 0 0 0 1 1]
        [0 0 0 0 0 1]
        sage: KM = __randomtop__(100, float(0.8))
        sage: SM = k2s_binary_matrix_sparse(KM)
        sage: SM == k2s_binary_matrix_sparse(s2k_binary_matrix_sparse(SM))
        True
    """
    data = __vector_to_list__(__matrice_to_lmtrx__(kmatrix)).python()
    dim = len(data)
    mat_dict = {}
    for j in range(dim):
        colj = data[j]
        for entry in colj:
            mat_dict[(entry[0], j)] = 1
    return matrix(dim, mat_dict)

def s2k_binary_matrix_sparse(smatrix):
    r"""
    Converts a binary matrix in SageMath to a Kenzo binary sparse matrice (type `matrice`).

    INPUT:

    - ``smatrix`` -- A binary matrix.

    EXAMPLES::

        sage: from sage.homology.finite_topological_spaces import k2s_binary_matrix_sparse, \
              s2k_binary_matrix_sparse
        sage: SM2 = matrix.ones(5)
        sage: s2k_binary_matrix_sparse(SM2)
        <ECL: 
        ========== MATRIX 5 lines + 5 columns =====
        L1=[C1=1][C2=1][C3=1][C4=1][C5=1]
        L2=[C1=1][C2=1][C3=1][C4=1][C5=1]
        L3=[C1=1][C2=1][C3=1][C4=1][C5=1]
        L4=[C1=1][C2=1][C3=1][C4=1][C5=1]
        L5=[C1=1][C2=1][C3=1][C4=1][C5=1]
        ========== END-MATRIX>
    """
    dim = smatrix.nrows()
    entries = []
    for entry in smatrix.dict().keys():
        entries.append([entry[0]+1, entry[1]+1])
    kentries = EclObject(entries)
    return __edges_to_matrice__(kentries, dim)

###############################################################


def FiniteSpace(data, elements=None, is_T0=False):
    r"""
    Construct a finite topological space from various forms of input data.

    INPUT:

    - ``data`` -- different input are accepted by this constructor:

      1. A dictionary representing the minimal basis of the space.

      2. A list or tuple of minimal open sets (in this case the elements of the 
         space are assumed to be ``range(n)`` where ``n`` is the length of ``data``).

      3. A topogenous matrix (assumed sparse). If ``elements=None``, the elements
         of the space are assumed to be ``range(n)`` where ``n`` is the dimension
         of the matrix.

      4. A finite poset (by now if ``poset._is_facade = False``, the methods are
         not completely tested).

    - ``elements`` -- (default ``None``) it is ignored when data is of type 1, 2
      or 4. When ``data`` is a topogenous matrix, this parameter gives the 
      underlying set of the space.
      
    - ``is_T0`` -- (default ``False``) it is a boolean that indicates, when it is 
      previously known, if the finite space is `T_0.

    EXAMPLES:

    A dictionary as ``data``::

        sage: from sage.homology.finite_topological_spaces import FiniteSpace
        sage: T = FiniteSpace({'a': {'a', 'c'}, 'b': {'b'}, 'c':{'a', 'c'}}) ; T
        Finite topological space of 3 points with minimal basis
         {'a': {'a', 'c'}, 'b': {'b'}, 'c': {'a', 'c'}}
        sage: type(T)
        <class 'sage.homology.finite_topological_spaces.FiniteTopologicalSpace'>
        sage: FiniteSpace({'a': {'a', 'b'}})
        Traceback (most recent call last):
        ...
        ValueError: The data does not correspond to a valid dictionary
        sage: FiniteSpace({'a': {'a', 'b'}, 'b': {'a', 'b'}, 'c': {'a', 'c'}})
        Traceback (most recent call last):
        ...
        ValueError: The introduced data does not define a topology

    When ``data`` is a tuple or a list, the elements are in ``range(n)`` where
    ``n`` is the length of ``data``::

        sage: from sage.homology.finite_topological_spaces import FiniteSpace
        sage: T = FiniteSpace([{0, 3}, {1, 3}, {2, 3}, {3}]) ; T
        Finite T0 topological space of 4 points with minimal basis
         {0: {3, 0}, 1: {3, 1}, 2: {3, 2}, 3: {3}}
        sage: type(T)
        <class 'sage.homology.finite_topological_spaces.FiniteTopologicalSpace_T0'>
        sage: T.elements()
        [3, 0, 1, 2]
        sage: FiniteSpace(({0, 2}, {0, 2}))
        Traceback (most recent call last):
        ...
        ValueError: This kind of data assume the elements are in range(2)

    If ``data`` is a topogenous matrix, the parameter ``elements``, when it is not
    ``None``, determines the list of elements of the space::

        sage: from sage.homology.finite_topological_spaces import FiniteSpace
        sage: mat_dict = {(0, 0): 1, (0, 3): 1, (0, 4): 1, (1, 1): 1, (1, 2): 1, (2, 1): 1, \
        ....:             (2, 2): 1, (3, 3): 1, (3, 4): 1, (4, 3): 1, (4, 4): 1}
        sage: mat = matrix(mat_dict) ; mat
        [1 0 0 1 1]
        [0 1 1 0 0]
        [0 1 1 0 0]
        [0 0 0 1 1]
        [0 0 0 1 1]
        sage: T = FiniteSpace(mat) ; T
        Finite topological space of 5 points with minimal basis
         {0: {0}, 1: {1, 2}, 2: {1, 2}, 3: {0, 3, 4}, 4: {0, 3, 4}}
        sage: T.elements()
        [0, 1, 2, 3, 4]
        sage: M = FiniteSpace(mat, elements=(5, 'e', 'h', 0, 'c')) ; M
        Finite topological space of 5 points with minimal basis
         {5: {5}, 'e': {'e', 'h'}, 'h': {'e', 'h'}, 0: {5, 0, 'c'}, 'c': {5, 0, 'c'}}
        sage: M.elements()
        [5, 'e', 'h', 0, 'c']
        sage: FiniteSpace(mat, elements=[5, 'e', 'h', 0, 0])
        Traceback (most recent call last):
        ...
        AssertionError: Not valid list of elements

    Finally, when ``data`` is a finite poset, the corresponding finite T0 space
    is constructed::

        sage: from sage.homology.finite_topological_spaces import FiniteSpace
        sage: P = Poset([[1, 2], [4], [3], [4], []])
        sage: T = FiniteSpace(P) ; T
        Finite T0 topological space of 5 points with minimal basis
         {0: {0}, 1: {0, 1}, 2: {0, 2}, 3: {0, 2, 3}, 4: {0, 1, 2, 3, 4}}
        sage: type(T)
        <class 'sage.homology.finite_topological_spaces.FiniteTopologicalSpace_T0'>
        sage: T.poset() == P
        True
    """
    if hasattr(data, '_hasse_diagram'): # isinstance(data, FinitePosets): # type 4
        minimal_basis = {x: set(data.order_ideal([x])) for x in data.list()}
        topogenous = data.lequal_matrix()
        return FiniteTopologicalSpace_T0(elements=data.list(), minimal_basis=minimal_basis,
                                         topogenous=topogenous, poset=data)

    topogenous = None

    if isinstance(data, dict): # type 1
        n = len(data)
        eltos = set()
        for B in data.values():
            eltos = eltos.union(B)
        if not eltos==set(data):
            raise ValueError("The data does not correspond to a valid dictionary")
        basis = data

    if isinstance(data, (list, tuple)): # type 2
        n = len(data)
        eltos = set()
        # In this case, the elements are assumed being range(n)
        for B in data:
            eltos = eltos.union(B)
        if not eltos==set(range(n)):
            raise ValueError("This kind of data assume the elements are in range({})" \
                             .format(n))
        basis = dict(zip(range(n), data))

    if isinstance(data, Matrix_integer_sparse): # type 3
        n = data.dimensions()[0]
        assert n==data.dimensions()[1], \
               "Topogenous matrices are square"
        assert set(data.dict().values())=={1}, \
               "Topogenous matrices must have entries in {0,1}"
        basis = {}
        # Extracting a minimal basis from the topogenous matrix info
        if elements:
            if not isinstance(elements, (list, tuple)):
                raise ValueError("Parameter 'elements' must be a list or a tuple")
            assert len(set(elements))==n, \
                   "Not valid list of elements"
            for j in range(n):
                Uj = set([elements[i] for i in data.nonzero_positions_in_column(j)])
                basis[elements[j]] = Uj
            eltos = elements
        else:
            for j in range(n):
                Uj = set(data.nonzero_positions_in_column(j))
                basis[j] = Uj
            eltos = range(n)

    # This fixes a topological sort (it guarantees an upper triangular topogenous matrix)
    eltos = list(eltos)
    sorted_str_eltos = sorted([str(x) for x in eltos])
    eltos.sort(key = lambda x: (len(basis[x]), sorted_str_eltos.index(str(x))))

    # Now, check that 'basis' effectively defines a minimal basis for a topology
    nonzero = {(eltos.index(x), j):1 for j in range(n) \
               for x in basis[eltos[j]]}
    topogenous = matrix(n, nonzero)
    squared = topogenous*topogenous
    if not topogenous.nonzero_positions() == squared.nonzero_positions():
        raise ValueError("The introduced data does not define a topology")

    if is_T0:
        return FiniteTopologicalSpace_T0(elements=eltos, minimal_basis=basis,
                                         topogenous=topogenous)
    # Determine if the finite space is T0
    partition = []
    eltos2 = eltos.copy()
    while eltos2:
        x = eltos2.pop(0)
        Ux = basis[x] - set([x])
        equiv_class = set([x])
        for y in Ux:
            if x in basis[y]:
                equiv_class = equiv_class.union(set([y]))
                eltos2.remove(y)
        partition.append(equiv_class)

    if len(partition)==n:
        return FiniteTopologicalSpace_T0(elements=eltos, minimal_basis=basis,
                                         topogenous=topogenous)
    result = FiniteTopologicalSpace(elements=eltos, minimal_basis=basis,
                                    topogenous=topogenous)
    setattr(result, '_T0', partition)
    return result

def RandomFiniteT0Space(*args):
    r"""
    Returns a random finite `T_0` space.

    INPUT:

    - ``args`` -- A tuple of two arguments. The first argument must be an integer
     number, while the second argument must be either a number between 0 and 1, or
     ``True``.

    OUTPUT:

    - If ``args[1]``=``True``, a random finite `T_0` space of cardinality ``args[0]``
      of height 3 without beat points is returned.

    - If ``args[1]`` is a number, a random finite `T_0` space of cardinality ``args[0]``
      and density ``args[1]`` of ones in its topogenous matrix is returned.

    EXAMPLES::

        sage: from sage.homology.finite_topological_spaces import RandomFiniteT0Space
        sage: RandomFiniteT0Space(5, 0)
        Finite T0 topological space of 5 points with minimal basis 
         {0: {0}, 1: {1}, 2: {2}, 3: {3}, 4: {4}}
        sage: RandomFiniteT0Space(5, 2)
        Finite T0 topological space of 5 points with minimal basis
         {0: {0}, 1: {0, 1}, 2: {0, 1, 2}, 3: {0, 1, 2, 3}, 4: {0, 1, 2, 3, 4}}
        sage: RandomFiniteT0Space(6, True)
        Finite T0 topological space of 6 points with minimal basis
         {0: {0}, 1: {1}, 2: {0, 1, 2}, 3: {0, 1, 3}, 4: {0, 1, 2, 3, 4}, 5: {0, 1, 2, 3, 5}}
        sage: RandomFiniteT0Space(150, 0.2)
        Finite T0 topological space of 150 points
        sage: RandomFiniteT0Space(5, True)
        Traceback (most recent call last):
        ...
        AssertionError: The first argument must be an integer number greater than 5
    """
    assert len(args)==2, "Two arguments must be given"
    assert args[0].is_integer(), "The first argument must be an integer number"
    if args[1]==True:
        assert args[0]>5, "The first argument must be an integer number greater than 5"
        kenzo_top = __random_top_2space__(args[0])
    else:
        kenzo_top = __randomtop__(args[0], EclObject(float(args[1])))
    topogenous = k2s_binary_matrix_sparse(kenzo_top)
    basis = {j:set(topogenous.nonzero_positions_in_column(j)) for j in range(args[0])}
    return FiniteTopologicalSpace_T0(elements=list(range(args[0])), minimal_basis=basis,
                                     topogenous=topogenous)


class FiniteTopologicalSpace(Parent):
    r"""
    Finite topological spaces.

    Users should not call this directly, but instead use :func:`FiniteSpace`.
    See that function for more documentation.
    """
    def __init__(self, elements, minimal_basis, topogenous):
        r"""
        Define a finite topological space.

        INPUT:

        - ``elements`` -- list of the elements of the space. 

        - ``minimal_basis`` -- a dictionary where the values are sets representing
          the minimal open sets containing the respective key.

        - ``topogenous`` -- a topogenous matrix of the finite space corresponding
          to the order given by ``elements``.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteTopologicalSpace
            sage: elements = [1, 2, 'a', 3]
            sage: minimal_basis = {'a': {3, 'a'}, 3: {3, 'a'}, 2: {2, 1}, 1: {1}}
            sage: mat_dict = {(0, 0): 1, (0, 1): 1, (1, 1): 1, (2, 2): 1, \
            ....:             (2, 3): 1, (3, 2): 1, (3, 3): 1}
            sage: T = FiniteTopologicalSpace(elements, minimal_basis, matrix(mat_dict)) ; T
            Finite topological space of 4 points with minimal basis
             {'a': {'a', 3}, 3: {'a', 3}, 2: {1, 2}, 1: {1}}
            sage: T.topogenous_matrix() == matrix(mat_dict)
            True
        """
        # Assign attributes
        self._cardinality = len(elements)
        self._elements = elements
        self._minimal_basis = minimal_basis
        self._topogenous = topogenous

    def space_sorting(self, element):
        r"""
        Return a pair formed by the index of `element` in `self._elements` and 
        the index of `str(element)` in the sorted list consisting of the strings of
        elements in `self._elements`.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace({0: {3, 0}, 3: {3, 0}, 2: {2, 1}, 1: {1}})
            sage: T._elements
            [1, 0, 2, 3]
            sage: T.space_sorting(1)
            (0, 1)
            sage: T.space_sorting(2)
            (2, 2)
        """
        eltos = self._elements
        sorted_str_eltos = sorted([str(x) for x in eltos])
        return (eltos.index(element), sorted_str_eltos.index(str(element)))

    def _repr_(self):
        r"""
        Print representation.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: FiniteSpace({0: {0, 1}, 1: {0, 1}})
            Finite topological space of 2 points with minimal basis
             {0: {0, 1}, 1: {0, 1}}
            sage: Q = Poset((divisors(120), attrcall("divides")), linear_extension=True)
            sage: FiniteSpace(Q)
            Finite T0 topological space of 16 points
        """
        n = self._cardinality
        if n < 10:
            sorted_minimal_basis = {x: sorted(self._minimal_basis[x], key=self.space_sorting)
                                    for x in self._minimal_basis}
            return "Finite topological space of {} points with minimal basis \n {}" \
                   .format(n, sorted_minimal_basis).replace('[', '{').replace(']', '}')
        else:
            return "Finite topological space of {} points".format(n)

    def __contains__(self, x):
        r"""
        Return ``True`` if ``x`` is an element of the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: P = Poset((divisors(6), attrcall("divides")), linear_extension=True)
            sage: T = FiniteSpace(P)
            sage: 3 in T
            True
            sage: 4 in T
            False
        """
        return x in self._elements

    def elements(self):
        r"""
        Return the list of elements in the underlying set of the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace(({0}, {1}, {2, 3}, {3}))
            sage: T.elements()
            [0, 1, 3, 2]
        """
        return self._elements
        
    def underlying_set(self):
        r"""
        Return the underlying set of the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace(({0}, {1}, {2, 3}, {3}))
            sage: T.underlying_set()
            {0, 1, 2, 3}
        """
        return set(self._elements)

    def subspace(self, points=None, is_T0=False):
        r"""
        Return the subspace whose elements are in ``points``.

        INPUT:

        - ``points`` -- (default ``None``) A tuple, list or set contained in ``self.elements()``.

        - ``is_T0`` -- if it is known that the resulting subspace is T0, fix ``True``
          (default ``False``).

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace(({0}, {1, 3, 4}, {0, 2, 5}, {1, 3, 4}, {1, 3, 4}, {0, 2, 5}))
            sage: T.subspace((0, 3, 5))
            Finite T0 topological space of 3 points with minimal basis
             {0: {0}, 3: {3}, 5: {0, 5}}
            sage: T.subspace([4])
            Finite T0 topological space of 1 points with minimal basis
             {4: {4}}
            sage: T.subspace() == T
            True
        """
        if points is None:
            return self
        assert isinstance(points, (tuple, list, set)), \
               "Parameter must be of type tuple, list or set"
        points = set(points)
        assert points <= set(self._elements), \
               "There are points that are not in the space"
        if points==set(self._elements):
            return self
        minimal_basis = {x: self._minimal_basis[x] & points for x in points}
        return FiniteSpace(minimal_basis, is_T0=is_T0)

    def cardinality(self):
        r"""
        Return the number of elements in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: P = Poset((divisors(360), attrcall("divides")), linear_extension=True)
            sage: T = FiniteSpace(P)
            sage: T.cardinality() == P.cardinality()
            True
        """
        return self._cardinality

    def minimal_basis(self):
        r"""
        Return the minimal basis that generates the topology of the finite space.

        OUTPUT:

        - A dictionary whose keys are the elements of the space and the values
          are the minimal open sets containing the respective element.

        EXAMPLES::
        
            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace(({0}, {0, 1, 2}, {0, 1, 2}, {3, 4}, {3, 4}))
            sage: T.minimal_basis()
            {0: {0}, 1: {0, 1, 2}, 2: {0, 1, 2}, 3: {3, 4}, 4: {3, 4}}
            sage: M = T.equivalent_T0()
            sage: M.minimal_basis()
            {0: {0}, 1: {0, 1}, 3: {3}}
        """
        return self._minimal_basis

    def minimal_open_set(self, x):
        r"""
        Return the minimal open set containing ``x``.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace(({0}, {0, 1, 2}, {0, 1, 2}, {3, 4}, {3, 4}))
            sage: T.minimal_open_set(1)
            {0, 1, 2}
        """
        if not x in self:
            raise ValueError("The point {} is not an element of the space".format(x))
        else:
            return self._minimal_basis[x]

    def topogenous_matrix(self):
        r"""
        Return the topogenous matrix of the finite space.

        OUTPUT:

        - A binary matrix whose `(i,j)` entry is equal to 1 if and only if ``self._elements[i]``
          is in ``self._minimal_basis[self._elements[j]]``.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace(({0}, {1, 3, 4}, {0, 2, 5}, {1, 3, 4}, {1, 3, 4}, {0, 2, 5}))
            sage: T.topogenous_matrix()
            [1 0 1 0 0 1]
            [0 1 0 1 1 0]
            [0 0 1 0 0 1]
            [0 1 0 1 1 0]
            [0 1 0 1 1 0]
            [0 0 1 0 0 1]
            sage: T0 = T.equivalent_T0()
            sage: T0.topogenous_matrix()
            [1 0 1]
            [0 1 0]
            [0 0 1]
        """
        return self._topogenous

    def is_T0(self):
        r"""
        Return ``True`` if the finite space satisfies the T0 separation axiom.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0}, {1}, {2, 3}, {2, 3}])
            sage: T.is_T0()
            False
            sage: T.equivalent_T0().is_T0()
            True
        """
        return isinstance(self, FiniteTopologicalSpace_T0)

    def equivalent_T0(self, points=None, check=True):
        r"""
        Return a finite T0 space homotopy equivalent to ``self``.

        INPUT:

        - ``points`` -- (default ``None``) a tuple, list or set of representatives
          elements of the equivalent classes induced by the partition ``self._T0``.

        - ``check`` -- if ``True`` (default), it is checked that ``points`` effectively
          defines a set of representatives of the partition ``self._T0``.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace(({0}, {1, 3, 4}, {0, 2, 5}, {1, 3, 4}, {1, 3, 4}, {0, 2, 5}))
            sage: T.is_T0()
            False
            sage: T._T0
            [{0}, {1, 3, 4}, {2, 5}]
            sage: M1 = T.equivalent_T0()
            sage: M1.is_T0()
            True
            sage: M1.elements()
            [0, 1, 2]
            sage: M2 = T.equivalent_T0(points={0,4,5}, check=False)
            sage: M2.elements()
            [0, 4, 5]
            sage: T.equivalent_T0(points={0,3,4})
            Traceback (most recent call last):
            ...
            ValueError: Parameter 'points' is not a valid set of representatives
        """
        if self._T0 is True:
            return self
        else:
            if points is None:
                points = [list(A)[0] for A in self._T0]
            elif check:
                assert isinstance(points, (tuple, list, set)), \
                       "Parameter 'points' must be of type tuple, list or set"
                assert len(points)==len(self._T0), \
                       "Parameter 'points' does not have a valid length"
                points2 = set(points.copy())
                partition = self._T0.copy()
                while points2:
                    x = points2.pop()
                    class_x = None
                    for k in range(len(partition)):
                        if x in partition[k]:
                            class_x = k
                            partition.pop(k)
                            break
                    if class_x is None:
                            raise ValueError("Parameter 'points' is not a valid set of representatives")
            return self.subspace(points, is_T0=True)

    def Ux(self, x):
        r"""
        Return the list of the elements in the minimal open set containing ``x``.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {5}, 6: {5, 6}, 3: {3, 5}, 2: {2, 5, 6}, \
                                   4: {2, 4, 5, 6}, 1: {1, 5}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.Ux(2)
            [5, 6, 2]

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.2))
            sage: x = random.choice(T._elements)
            sage: T.is_contractible(T.Ux(x))
            True
        """
        return sorted(self._minimal_basis[x], key=self.space_sorting)

    def Fx(self, x):
        r"""
        Return the list of the elements in the closure of `\lbrace x\rbrace`.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {5}, 6: {5, 6}, 3: {3, 5}, 2: {2, 5, 6}, \
                                   4: {2, 4, 5, 6}, 1: {1, 5}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.Fx(2)
            [2, 4]

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.2))
            sage: x = random.choice(T._elements)
            sage: T.is_contractible(T.Fx(x))
            True
        """
        result = [y for y in self._elements if x in self._minimal_basis[y]]
        if result==[]:
            raise ValueError("The point {} is not an element of the space".format(x))
        return result

    def Cx(self, x):
        r"""
        Return the list of the elements in the star of ``x``.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {5}, 6: {5, 6}, 3: {3, 5}, 2: {2, 5, 6}, \
                                   4: {2, 4, 5, 6}, 1: {1, 5}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.Cx(2)
            [5, 6, 2, 4]

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.2))
            sage: x = random.choice(T._elements)
            sage: T.is_contractible(T.Cx(x))
            True
        """
        return self.Ux(x) + self.Fx(x)[1:]

    def Ux_tilded(self, x):
        r"""
        Return the list of the elements in `\widehat{U}_x = U_x \minus \lbrace x\rbrace`.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {5}, 6: {5, 6}, 3: {3, 5}, 2: {2, 5, 6}, \
                                   4: {2, 4, 5, 6}, 1: {1, 5}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.Ux_tilded(2)
            [5, 6]
        """
        return self.Ux(x)[:-1]

    def Fx_tilded(self, x):
        r"""
        Return the list of the elements in `\widehat{F}_x = F_x \minus \lbrace x\rbrace`.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {5}, 6: {5, 6}, 3: {3, 5}, 2: {2, 5, 6}, \
                                   4: {2, 4, 5, 6}, 1: {1, 5}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.Fx_tilded(2)
            [4]
        """
        return self.Fx(x)[1:]

    def Cx_tilded(self, x):
        r"""
        Return the list of the elements in `\widehat{C}_x = C_x \minus \lbrace x\rbrace`.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {5}, 6: {5, 6}, 3: {3, 5}, 2: {2, 5, 6}, \
                                   4: {2, 4, 5, 6}, 1: {1, 5}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.Cx_tilded(2)
            [5, 6, 4]
        """
        return self.Ux(x)[:-1] + self.Fx(x)[1:]

    def opposite(self):
        r"""
        Return the opposite space of ``self``.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: mat_dict = {(0, 0): 1, (0, 3): 1, (0, 4): 1, (1, 1): 1, (1, 2): 1, (2, 1): 1, \
            ....:             (2, 2): 1, (3, 3): 1, (3, 4): 1, (4, 3): 1, (4, 4): 1}
            sage: T = FiniteSpace(matrix(mat_dict))
            sage: T
            Finite topological space of 5 points with minimal basis 
             {0: {0}, 1: {1, 2}, 2: {1, 2}, 3: {0, 3, 4}, 4: {0, 3, 4}}
            sage: T.opposite()
            Finite topological space of 5 points with minimal basis 
             {0: {3, 4, 0}, 1: {1, 2}, 2: {1, 2}, 3: {3, 4}, 4: {3, 4}}
            sage: T.topogenous_matrix()
            [1 0 0 1 1]
            [0 1 1 0 0]
            [0 1 1 0 0]
            [0 0 0 1 1]
            [0 0 0 1 1]
            sage: T.opposite().topogenous_matrix()
            [1 1 0 0 0]
            [1 1 0 0 0]
            [0 0 1 1 1]
            [0 0 1 1 1]
            [0 0 0 0 1]
        """
        minimal_basis_op = {x:set(self.Fx(x)) for x in self._elements}
        T0 = isinstance(self, FiniteTopologicalSpace_T0)
        return FiniteSpace(minimal_basis_op, is_T0=T0)

    def is_interior_point(self, x, E):
        r"""
        Return ``True`` if ``x`` is an interior point of ``E`` in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.is_interior_point(1, {1, 2, 3})
            True
            sage: T.is_interior_point(2, {1, 2, 3})
            False
            sage: T.is_interior_point(1, set())
            False
            sage: T.is_interior_point(3, T.underlying_set())
            True
        """
        assert x in self.underlying_set() , "Parameter 'x' must be an element of the space"
        assert E <= self.underlying_set() , "Parameter 'E' must be a subset of the underlying set"
        if not x in E:
            return False
        return self._minimal_basis[x] <= E

    def interior(self, E):
        r"""
        Return the interior of a subset in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.interior({1, 2, 3})
            {1}
            sage: T.interior({1, 2, 3, 4})
            {1, 2, 3, 4}
            sage: T.interior({2, 3})
            set()

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.5))
            sage: X = T.underlying_set()
            sage: k = randint(0,len(X))
            sage: E = set(random.sample(X, k))
            sage: Int = T.interior(E)
            sage: T.is_open(Int)
            True
            sage: T.interior(Int) == Int
            True
            sage: Int == X - T.closure(X - E)
            True
            sage: m = randint(0,len(X))
            sage: M = set(random.sample(X, m))
            sage: T.interior(E & M) == Int & T.interior(M)
            True
        """
        X = self.underlying_set()
        if E == X or E == set():
            return E
        assert E < X , "The parameter must be a subset of the underlying set"
        return set([x for x in E if self.is_interior_point(x, E)])

    def is_open(self, E):
        r"""
        Return ``True`` if ``E`` is an open subset of the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.is_open({0})
            False
            sage: T.is_open({0, 1})
            True
            sage: T.is_open({0, 1, 4})
            True
            sage: T.is_open(set())
            True
        """
        return E == self.interior(E)

    def is_closed(self, E):
        r"""
        Return ``True`` if ``E`` is a closed subset of the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace({'a':{'a','b'},'b':{'a','b'},'c':{'c','d'},'d':{'d'}})
            sage: T.is_closed({'a','b','c'})
            True
            sage: T.is_closed({'b'})
            False
        """
        X = self.underlying_set() 
        return self.is_open(X - E)

    def is_exterior_point(self, x, E):
        r"""
        Return ``True`` if ``x`` is an exterior point of ``E`` in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.is_exterior_point(1, {2, 3})
            True
            sage: T.is_exterior_point(3, {0, 1, 2})
            False
        """
        return self._minimal_basis[x].isdisjoint(E)

    def exterior(self, E):
        r"""
        Return the exterior of a subset in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.exterior({2})
            {0, 1, 4}
            sage: T.exterior({2, 4})
            {0, 1}

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.5))
            sage: X = T.underlying_set()
            sage: k = randint(0,len(X))
            sage: E = set(random.sample(X, k))
            sage: Ext = T.exterior(E)
            sage: Ext.isdisjoint(E)
            True
            sage: Ext == T.interior(X - E)
            True
            sage: Ext == X - T.closure(E)
            True
            sage: T.interior(E) <= T.exterior(Ext)
            True
        """
        X = self.underlying_set()
        if E == X:
            return set()
        if E == set():
            return X
        assert E < X , "The parameter must be a subset of the underlying set"
        return set([x for x in X - E if self.is_exterior_point(x, E)])

    def is_boundary_point(self, x, E):
        r"""
        Return ``True`` if ``x`` is a boundary point of ``E`` in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.is_boundary_point(0, {1, 2, 3})
            True
            sage: T.is_boundary_point(1, {2, 3, 4})
            False
        """
        Ux = self._minimal_basis[x]
        return bool(Ux & E) and not bool(Ux <= E)

    def boundary(self, E):
        r"""
        Return the boundary of a subset in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.boundary({1})
            {0}
            sage: T.boundary({2, 3})
            {2, 3}

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.5))
            sage: X = T.underlying_set()
            sage: k = randint(0,len(X))
            sage: E = set(random.sample(X, k))
            sage: Fr = T.boundary(E)
            sage: T.is_closed(Fr)
            True
            sage: Fr == T.boundary(X - E)
            True
            sage: Fr == T.closure(E) - T.interior(E)
            True
            sage: Fr == T.closure(E) & T.closure(X - E)
            True
            sage: T.interior(E) == E - Fr
            True
            sage: T.boundary(Fr) <= Fr
            True
            sage: T.boundary(T.boundary(Fr)) == T.boundary(Fr)
            True
            sage: X == Fr.union(T.interior(E), T.exterior(E))
            True
        """
        X = self.underlying_set()
        if E == X or E == set():
            return set()
        assert E < X , "The parameter must be a subset of the underlying set"
        return set([x for x in X if self.is_boundary_point(x, E)])

    def is_limit_point(self, x, E):
        r"""
        Return ``True`` if ``x`` is a limit point of ``E`` in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.is_limit_point(0, {1})
            True
            sage: T.is_limit_point(1, {0, 1})
            False
        """
        Ux_minus_x = self._minimal_basis[x] - {x}
        return not Ux_minus_x.isdisjoint(E)

    def derived(self, E):
        r"""
        Return the derived set of a subset in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.derived({0, 1, 2})
            {0, 3}
            sage: T.derived({3, 4})
            {2, 3}

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.5))
            sage: X = T.underlying_set()
            sage: k = randint(0,len(X))
            sage: E = set(random.sample(X, k))
            sage: Der = T.derived(E)
            sage: T.derived(Der) <= E.union(Der)
            True
            sage: T.closure(E) == E.union(Der)
            True
        """
        X = self.underlying_set()
        if E == X or E == set():
            return E
        assert E < X , "The parameter must be a subset of the underlying set"
        return set([x for x in X if self.is_limit_point(x, E)])

    def is_closure_point(self, x, E):
        r"""
        Return ``True`` if ``x`` is a point of closure of ``E`` in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.is_closure_point(3, {1})
            False
            sage: T.is_closure_point(3, {1,2})
            True
        """
        return not self._minimal_basis[x].isdisjoint(E)

    def closure(self, E):
        r"""
        Return the closure of a subset in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.closure({0, 2})
            {0, 2, 3}
            sage: T.closure({0})
            {0}

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.5))
            sage: X = T.underlying_set()
            sage: k = randint(0,len(X))
            sage: E = set(random.sample(X, k))
            sage: Cl = T.closure(E)
            sage: T.is_closed(Cl)
            True
            sage: T.closure(Cl) == Cl
            True
            sage: Cl == X - T.interior(X - E)
            True
            sage: T.interior(T.boundary(Cl)) == set()
            True
            sage: Cl == E.union(T.boundary(E))
            True
            sage: m = randint(0,len(X))
            sage: M = set(random.sample(X, m))
            sage: T.closure(E.union(M)) == Cl.union(T.closure(M))
            True
        """
        X = self.underlying_set()
        if E == X or E == set():
            return E
        assert E < X , "The parameter must be a subset of the underlying set"
        return E.union(set([x for x in X - E if self.is_closure_point(x, E)]))

    def is_dense(self, E):
        r"""
        Return ``True`` if ``E`` is dense in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1, 2}, {0, 1, 2}, {2}])
            sage: T.is_dense({2})
            True
            sage: T.is_dense({0, 1})
            False
        """
        return self.closure(E) == self.underlying_set()

    def is_isolated_point(self, x, E=None):
        r"""
        Return ``True`` if ``x`` is an isolated point of ``E`` in the finite space.
        If ``E`` is ``None``, return ``True`` if ``x`` is an isolated point of 
        the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.is_isolated_point(0)
            False
            sage: T.is_isolated_point(0, {0, 2, 3, 4})
            True
        """
        if E:
            return (self._minimal_basis[x] & E) == set([x])
        else:
            return self._minimal_basis[x] == set([x])

    def isolated_set(self, E=None):
        r"""
        Return the set of isolated points of a subset in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace([{0, 1}, {1}, {2, 3, 4}, {2, 3, 4}, {4}])
            sage: T.isolated_set()
            {1, 4}
            sage: T.isolated_set({0, 2, 3, 4})
            {0, 4}

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.5))
            sage: X = T.underlying_set()
            sage: k = randint(0,len(X))
            sage: E = set(random.sample(X, k))
            sage: Iso = T.isolated_set(E)
            sage: T.closure(E) == Iso.union(T.derived(E))
            True
        """
        if E is None: 
            E = self.underlying_set()
        return set([x for x in E if self.is_isolated_point(x, E)])


class FiniteTopologicalSpace_T0(FiniteTopologicalSpace):
    r"""
    Finite topological spaces satisfying the T0 separation axiom (Kolmogorov spaces).

    Users should not call this directly, but instead use :func:`FiniteSpace`.
    See that function for more documentation.
    """
    def __init__(self, elements, minimal_basis, topogenous, poset=None):
        r"""
        Define a finite T0 topological space.

        INPUT:

        - ``elements`` -- list of the elements of the space. 

        - ``minimal_basis`` -- a dictionary where the values are sets representing
          the minimal open sets containing the respective key.

        - ``topogenous`` -- a topogenous matrix of the finite space corresponding
          to the order given by ``elements`` (it is assumed upper triangular).

        - ``poset`` -- a poset corresponding to the finite space (Alexandroff
          correspondence) (default ``None``).

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteTopologicalSpace_T0
            sage: elements = [0, 1, 2, 3]
            sage: minimal_basis = {0: {0}, 1: {0, 1}, 2: {0, 1, 2}, 3: {0, 3}}
            sage: mat_dict = {(0, 0): 1, (0, 1): 1, (0, 2): 1, (0, 3): 1, \
            ....:             (1, 1): 1, (1, 3): 1, (2, 2): 1, (3, 3): 1}
            sage: T = FiniteTopologicalSpace_T0(elements, minimal_basis, matrix(mat_dict)); T
            Finite T0 topological space of 4 points with minimal basis
             {0: {0}, 1: {0, 1}, 2: {0, 1, 2}, 3: {0, 3}}
        """
        FiniteTopologicalSpace.__init__(self, elements, minimal_basis, topogenous)
        if poset:
            # isinstance(poset, FinitePosets)
            assert hasattr(poset, '_hasse_diagram'), \
                   "Parameter 'poset' must be a real poset!"
            # Verify the coherence of the parameters
            assert set(self._elements)==set(poset.list()), \
                   "Elements of poset and minimal_basis do not coincide"
            self._elements = poset.list()
        else:
            # Construct the associated poset
            elmts = self._elements
            f = lambda x, y: self._topogenous[elmts.index(x), elmts.index(y)]==1
            poset = Poset((elmts, f), linear_extension=True) 
        self._poset = poset
        self._T0 = True

    def _repr_(self):
        r"""
        Print representation.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: P = Poset((divisors(6), attrcall("divides")), linear_extension=True)
            sage: FiniteSpace(P)
            Finite T0 topological space of 4 points with minimal basis
             {1: {1}, 2: {1, 2}, 3: {1, 3}, 6: {1, 2, 3, 6}}
            sage: Q = Poset((divisors(120), attrcall("divides")), linear_extension=True)
            sage: FiniteSpace(Q)
            Finite T0 topological space of 16 points
        """
        n = self._cardinality
        if n < 10:
            sorted_minimal_basis = {x: sorted(self._minimal_basis[x], key=self.space_sorting)
                                    for x in self._minimal_basis}
            return "Finite T0 topological space of {} points with minimal basis \n {}" \
                   .format(n, sorted_minimal_basis).replace('[', '{').replace(']', '}')
        else:
            return "Finite T0 topological space of {} points".format(n)

    def poset(self):
        r"""
        Return the corresponding poset of the finite space (Alexandroff correspondence).

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = ({0}, {0, 1}, {0, 1, 2}, {0, 3})
            sage: T = FiniteSpace(minimal_basis) ; T
            Finite T0 topological space of 4 points with minimal basis
             {0: {0}, 1: {0, 1}, 2: {0, 1, 2}, 3: {0, 3}}
            sage: T.poset()
            Finite poset containing 4 elements with distinguished linear extension
            sage: P = Poset((divisors(12), attrcall("divides")), linear_extension=True)
            sage: T = FiniteSpace(P)
            sage: T.poset() == P
            True
        """
        return self._poset

    def show(self, highlighted_edges=None):
        r"""
        Displays the Hasse diagram of the poset ``self._poset``.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: T = FiniteSpace(posets.RandomPoset(15, 0.2))
            sage: T.show()
            Graphics object consisting of 31 graphics primitives
        """
        if highlighted_edges:
            return self._poset.plot(cover_colors = {'blue': highlighted_edges})
        return self._poset.plot()

    def stong_matrix(self):
        r"""
        Returns the Stong matrix of the finite `T_0` space i.e. the adjacency matrix
        of the Hasse diagram of its associated poset, with ones in its diagonal.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: covers = [[9, 13], [7, 13], [4, 13], [8, 12], [7, 12], [5, 12],
            ....:           [9, 11], [6, 11], [5, 11], [8, 10], [6, 10], [4, 10],
            ....:           [3, 9], [2, 9], [3, 8], [2, 8], [3, 7], [1, 7], [3, 6],
            ....:           [1, 6], [2, 5], [1, 5], [2, 4], [1, 4]]
            sage: P = Poset((list(range(1,14)), covers), cover_relations=True)
            sage: X = FiniteSpace(P)
            sage: X.topogenous_matrix()
            [1 0 1 1 0 0 0 1 1 1 1 1 1]
            [0 1 1 1 0 1 1 0 1 1 0 1 1]
            [0 0 1 0 0 0 0 0 1 0 0 1 0]
            [0 0 0 1 0 0 0 0 0 1 0 0 1]
            [0 0 0 0 1 1 1 1 1 1 1 1 1]
            [0 0 0 0 0 1 0 0 1 0 0 0 1]
            [0 0 0 0 0 0 1 0 0 1 0 1 0]
            [0 0 0 0 0 0 0 1 1 1 0 0 0]
            [0 0 0 0 0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 1 1 1]
            [0 0 0 0 0 0 0 0 0 0 0 1 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 1]
            sage: X.stong_matrix()
            [1 0 1 1 0 0 0 1 0 0 1 0 0]
            [0 1 1 1 0 1 1 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 1 0 0 1 0]
            [0 0 0 1 0 0 0 0 0 1 0 0 1]
            [0 0 0 0 1 1 1 1 0 0 1 0 0]
            [0 0 0 0 0 1 0 0 1 0 0 0 1]
            [0 0 0 0 0 0 1 0 0 1 0 1 0]
            [0 0 0 0 0 0 0 1 1 1 0 0 0]
            [0 0 0 0 0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 1 1 1]
            [0 0 0 0 0 0 0 0 0 0 0 1 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 1]
        """
        return self._poset._hasse_diagram.adjacency_matrix(sparse=True) + matrix.identity(self._cardinality)

    def order_complex(self):
        r"""
        Return the order complex of the finite space i.e. the simplicial complex
        whose simplices are the nonempty chains of ``self.poset()``.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = ({0}, {0, 1}, {0, 1, 2}, {0, 3})
            sage: T = FiniteSpace(minimal_basis) ; T
            Finite T0 topological space of 4 points with minimal basis
             {0: {0}, 1: {0, 1}, 2: {0, 1, 2}, 3: {0, 3}}
            sage: T.order_complex()
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 3), (0, 1, 2)}
        """
        return self._poset.order_complex()

    def barycentric_subdivision(self):
        r"""
        Return the barycentric subdivision of the finite space i.e. the face poset
        of its order complex.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = ({0}, {0, 1}, {0, 1, 2}, {0, 3})
            sage: T = FiniteSpace(minimal_basis) ; T
            Finite T0 topological space of 4 points with minimal basis
             {0: {0}, 1: {0, 1}, 2: {0, 1, 2}, 3: {0, 3}}
            sage: T.barycentric_subdivision()
            Finite T0 topological space of 9 points with minimal basis 
             {(3,): {(3,)}, (2,): {(2,)}, (1,): {(1,)}, (1, 2): {(2,), (1,), (1, 2)},
             (0,): {(0,)}, (0, 1): {(1,), (0,), (0, 1)}, (0, 2): {(2,), (0,), (0, 2)},
             (0, 1, 2): {(2,), (1,), (1, 2), (0,), (0, 1), (0, 2), (0, 1, 2)},
             (0, 3): {(3,), (0,), (0, 3)}}
        """
        return FiniteSpace(self._poset.order_complex().face_poset(), is_T0=True)

    def is_down_beat_point(self, x, subspace=None):
        r"""
        Return ``True`` if ``x`` is a down beat point of the subspace of ``self``
        determined by ``subspace``.

        INPUT:

        - ``x`` - an element of the finite space. In case ``subspace`` is not
          ``None``, `x`` must be one of its elements.

        - ``subspace`` -- (default ``None``) a list of elements in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {5}, 4: {4}, 2: {2}, 6: {2, 4, 6}, \
                                   1: {1, 4}, 3: {1, 3, 4}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.is_down_beat_point(6)
            False
            sage: T.is_down_beat_point(6, [3, 4, 5, 6])
            True
        """
        xindex = self._elements.index(x)
        if subspace is None:
            subspaceindex = [i for i in range(xindex - 1,-1,-1) \
                             if self._topogenous[i, xindex]==1]
        else:
            sortsubspace = sorted(subspace, key=self._elements.index, reverse=True)
            subspaceindex = [self._elements.index(i) for i in sortsubspace \
                             if self._topogenous[self._elements.index(i), xindex]==1 \
                                and self._elements.index(i)!=xindex]
        if subspaceindex==[]:
            return False
        maximal = subspaceindex[0]
        for i in subspaceindex:
            if not self._topogenous[i, maximal]==self._topogenous[i, xindex]:
                return False
        return True

    def is_up_beat_point(self, x, subspace=None):
        r"""
        Return ``True`` if ``x`` is an up beat point of the subspace of ``self``
        determined by ``subspace``.

        INPUT:

        - ``x`` - an element of the finite space. In case ``subspace`` is not
          ``None``, `x`` must be one of its elements.

        - ``subspace`` -- (default ``None``) a list of elements in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {5}, 4: {4}, 2: {2}, 6: {2, 4, 6}, \
                                   1: {1, 4}, 3: {1, 3, 4}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.is_up_beat_point(4)
            False
            sage: T.is_up_beat_point(4, [1, 2, 3, 4, 5])
            True
        """
        xindex = self._elements.index(x)
        if subspace is None:
            subspaceindex = [j for j in range(xindex + 1, self._cardinality) \
                             if self._topogenous[xindex, j]==1]
        else:
            sortsubspace = sorted(subspace, key=self._elements.index)
            subspaceindex = [self._elements.index(i) for i in sortsubspace \
                             if self._topogenous[xindex, self._elements.index(i)]==1 \
                                and self._elements.index(i)!=xindex]
        if subspaceindex==[]:
            return False
        minimal = subspaceindex[0]
        for j in subspaceindex:
            if not self._topogenous[minimal, j]==self._topogenous[xindex, j]:
                return False
        return True

    def is_beat_point(self, x, subspace=None):
        r"""
        Return ``True`` if ``x`` is a beat point of the subspace of ``self``
        determined by ``subspace``.

        INPUT:

        - ``x`` - an element of the finite space. In case ``subspace`` is not
          ``None``, `x`` must be one of its elements.

        - ``subspace`` -- (default ``None``) a list of elements in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {5}, 4: {4}, 2: {2}, 6: {2, 4, 6}, \
                                   1: {1, 4}, 3: {1, 3, 4}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.is_beat_point(2)
            True
            sage: T.is_beat_point(2, [2, 3, 4, 5])
            False
        """
        if self._elements.index(x) < self._cardinality / 2:
            return self.is_down_beat_point(x, subspace) or self.is_up_beat_point(x, subspace)
        else:
            return self.is_up_beat_point(x, subspace) or self.is_down_beat_point(x, subspace)

    def core_list(self, subspace=None):
        r"""
        Return a list of elements in a core of the subspace of ``self`` determined
        by ``subspace``.

        INPUT:

        - ``subspace`` -- (default ``None``) a list of elements in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {4, 5}, 4: {4}, 2: {2}, 6: {2, 4, 6}, \
                                   1: {1, 2, 4}, 3: {1, 2, 4, 3}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.core_list()
            [2, 4, 6, 3]
            sage: T.core_list([3, 2, 1, 4, 5, 6])
            [2, 1, 4, 6]
            sage: T.core_list([1, 2, 3, 4, 5])
            [5]

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.2))
            sage: X = T._elements
            sage: k = randint(0,len(X))
            sage: E1 = random.sample(X, k)
            sage: E2 = random.sample(E1, k)
            sage: len(T.core_list(E1)) == len(T.core_list(E2)) # cores are homeomorphic
            True
        """
        if subspace is None:
            subspace = self._elements
        beatpoint = None
        for x in subspace:
            if self.is_beat_point(x, subspace):
                beatpoint = x
                break
        if beatpoint is None:
            return subspace
        else:
            return self.core_list([y for y in subspace if y != beatpoint])

    def core(self, subspace=None):
        r"""
        Return a core of the subspace of ``self`` determined by ``subspace``.

        INPUT:

        - ``subspace`` -- (default ``None``) a list of elements in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {4, 5}, 4: {4}, 2: {2}, 6: {2, 4, 6}, \
                                   1: {1, 2, 4}, 3: {1, 2, 4, 3}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.core()
            Finite T0 topological space of 4 points with minimal basis 
             {2: {2}, 3: {2, 4, 3}, 4: {4}, 6: {2, 4, 6}}
            sage: T.core([3,2,1,4,5,6])
            Finite T0 topological space of 4 points with minimal basis 
             {1: {2, 4, 1}, 2: {2}, 4: {4}, 6: {2, 4, 6}}
            sage: T.core([1,2,3,4,5])
            Finite T0 topological space of 1 points with minimal basis 
             {5: {5}}
        """
        return self.subspace(self.core_list(subspace), is_T0=True)

    def is_contractible(self, subspace=None):
        r"""
        Return ``True`` if the finite space is contractible (in the setting of finite spaces,
        this is equivalent to say that its cores are singletons).

        INPUT:

        - ``subspace`` -- (default ``None``) a list of elements in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {5: {4, 5}, 4: {4}, 2: {2}, 6: {2, 4, 6}, \
                                   1: {1, 2, 4}, 3: {1, 2, 4, 3}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.is_contractible()
            False
            sage: T.is_contractible([1,2,3,4,5])
            True

        TESTS::

            sage: import random
            sage: P = posets.RandomPoset(20, 0.5)
            sage: X = P.list()
            sage: k = randint(0,len(X))
            sage: E = random.sample(X, k)
            sage: S = P.subposet(E)
            sage: F = FiniteSpace(S)
            sage: S.has_top()==False or F.is_contractible()
            True
            sage: S.has_bottom()==False or F.is_contractible()
            True
        """
        return len(self.core_list(subspace))==1

    def is_weak_point(self, x, subspace=None):
        r"""
        Return ``True`` if ``x`` is a weak beat point of the subspace of ``self``
        determined by ``subspace``.

        INPUT:

        - ``x`` - an element of the finite space. In case ``subspace`` is not
          ``None``, `x`` must be one of its elements.

        - ``subspace`` -- (default ``None``) a list of elements in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {2: {2}, 5: {2, 5}, 1: {1}, 3: {1, 2, 3}, \
                                   4: {1, 2, 4}, 7: {1, 2, 3, 4, 7}, \
                                   6: {1, 2, 3, 4, 5, 6, 7}, \
                                   8: {1, 2, 3, 4, 5, 6, 7, 8}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.is_beat_point(1)
            False
            sage: T.is_weak_point(1)
            True

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.2))
            sage: X = T._elements
            sage: k = randint(0,len(X))
            sage: E = random.sample(X, k)
            sage: x = random.choice(E)
            sage: T.is_beat_point(x, E)==False or T.is_beat_point(x, E)==T.is_weak_point(x, E)
            True
        """
        subspaceU = self.Ux_tilded(x)
        subspaceF = self.Fx_tilded(x)
        if subspace is not None:
            subspaceU = list(set(subspaceU) & set(subspace))
            subspaceF = list(set(subspaceF) & set(subspace))
        if self._elements.index(x) < self._cardinality / 2:
            return self.is_contractible(subspaceU) or self.is_contractible(subspaceF)
        else:
            return self.is_contractible(subspaceF) or self.is_contractible(subspaceU)

    def weak_core_list(self, subspace=None):
        r"""
        Return a list of elements in a weak core (finite space with no weak points)
        of the subspace of ``self`` determined by ``subspace``.

        INPUT:

        - ``subspace`` -- (default ``None``) a list of elements in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {2: {2}, 5: {2, 5}, 1: {1}, 3: {1, 2, 3}, \
                                   4: {1, 2, 4}, 7: {1, 2, 3, 4, 7}, \
                                   6: {1, 2, 3, 4, 5, 6, 7}, \
                                   8: {1, 2, 3, 4, 5, 6, 7, 8}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.weak_core_list()
            [8]
            sage: T.weak_core_list([1,2,3,4,5])
            [1, 2, 3, 4]

        TESTS::

            sage: import random
            sage: T = FiniteSpace(posets.RandomPoset(30, 0.5))
            sage: X = T._elements
            sage: k = randint(0,len(X))
            sage: E = random.sample(X, k)
            sage: len(T.weak_core_list(E)) <= len(T.core_list(E))
            True
        """
        realsubspace = subspace or self._elements
        weakpoint = None
        for x in realsubspace:
            if self.is_beat_point(x, subspace) or self.is_weak_point(x, subspace):
                weakpoint = x
                break
        if weakpoint is None:
            return realsubspace
        else:
            return self.weak_core_list([y for y in realsubspace if y != weakpoint])

    def weak_core(self, subspace=None):
        r"""
        Return a weak core (finite space with no weak points) of the subspace of
        ``self`` determined by ``subspace``.

        INPUT:

        - ``subspace`` -- (default ``None``) a list of elements in the finite space.

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: minimal_basis = {2: {2}, 5: {2, 5}, 1: {1}, 3: {1, 2, 3}, \
                                   4: {1, 2, 4}, 7: {1, 2, 3, 4, 7}, \
                                   6: {1, 2, 3, 4, 5, 6, 7}, \
                                   8: {1, 2, 3, 4, 5, 6, 7, 8}}
            sage: T = FiniteSpace(minimal_basis)
            sage: T.weak_core()
            Finite T0 topological space of 1 points with minimal basis 
             {8: {8}}
            sage: T.weak_core([1,2,3,4,5])
            Finite T0 topological space of 4 points with minimal basis 
             {1: {1}, 2: {2}, 3: {1, 2, 3}, 4: {1, 2, 4}}             
        """
        return self.subspace(self.weak_core_list(subspace), is_T0=True)

    def discrete_vector_field(self, h_admissible=None):
        r"""
        Return a discrete vector field on the finite `T_0` space i.e. a homologically
        admissible Morse matching on the Hasse diagram of the associated poset.

        INPUT:

        - ``h_admissible`` -- (default ``None``) If it is ``True``, all the edges
          `(x, y)` of the Hasse diagram are assumed to be homologically admissible
          i.e. tha subspace `\widehat{U}_y - \lbrace x\rbrace` is homotopically 
          trivial (this can be assumed when the finite space is a barycentric
          subdivision).

        EXAMPLES::

            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: Pcovers = [[1, 2], [2, 3], [3, 4], [3, 5], [4, 6], [5, 6], [6, 7],
            ....:            [6, 8], [8, 9], [9, 10], [1, 11], [7, 12], [9, 12],
            ....:            [7, 13], [10, 13], [11, 13], [8, 15], [7, 16], [8, 16],
            ....:            [11, 16], [15, 17], [2, 19], [6, 20], [18, 20]]
            sage: P = Poset((list(range(1,21)), Pcovers), cover_relations=True)
            sage: X = FiniteSpace(P)
            sage: dvf = X.discrete_vector_field(); dvf
            [(2, 3), (4, 6), (8, 9), (7, 12), (15, 17), (18, 20), (10, 13), (11, 16)]
            sage: X.show(dvf)
            Graphics object consisting of 45 graphics primitives
            sage: Qcovers = [[1, 2], [2, 3], [3, 4], [3, 5]]
            sage: Q = Poset((list(range(1,6)), Qcovers), cover_relations=True)
            sage: Y = FiniteSpace(Q)
            sage: Z = Y.barycentric_subdivision()
            sage: dvf = Z.discrete_vector_field(h_admissible=True)
            sage: Z.show(dvf)
            Graphics object consisting of 71 graphics primitives
        """
        kenzo_top = s2k_binary_matrix_sparse(self._topogenous)
        kenzo_dvfield = EclListIterator(__dvfield_aux__(kenzo_top, None, h_admissible))
        result = []
        for vector in kenzo_dvfield:
            vectorpy = vector.python()
            result.append((self._elements[vectorpy[0]-1], self._elements[vectorpy[1]-1]))
        return result

    def hregular_homology(self, deg=None, dvfield=None):
        r"""
        The homology of an h-regular finite space.

        INPUT:

        - ``deg`` -- an element of the grading group for the chain
          complex (default ``None``); the degree in which
          to compute homology -- if this is ``None``, return the
          homology in every degree in which the chain complex is
          possibly nonzero.

        - ``dvfield`` -- (default ``None``) a list of edges representing a discrete
          vector field on the finite space.

        EXAMPLES::
        
            sage: from sage.homology.finite_topological_spaces import FiniteSpace
            sage: covers = [[9, 13], [7, 13], [4, 13], [8, 12], [7, 12], [5, 12], [9, 11],
            ....:           [6, 11], [5, 11], [8, 10], [6, 10], [4, 10], [3, 9], [2, 9],
            ....:           [3, 8], [2, 8], [3, 7], [1, 7], [3, 6], [1, 6], [2, 5], [1, 5],
            ....:           [2, 4], [1, 4]]
            sage: P = Poset((list(range(1,14)), covers), cover_relations = True)
            sage: X = FiniteSpace(P)
            sage: X.hregular_homology()
            {0: Z, 1: C2, 2: 0}
            sage: dvf = X.discrete_vector_field()
            sage: X.show(dvf)
            Graphics object consisting of 38 graphics primitives
            sage: X.hregular_homology(dvfield = dvf)
            {0: Z, 1: C2, 2: 0}
        """
        assert deg==None or deg.is_integer(), "The degree must be an integer number or None"
        height = self._poset.height()
        if deg and (deg < 0 or deg >= height):
            return HomologyGroup(0, ZZ)
        kenzo_stong = s2k_binary_matrix_sparse(self.stong_matrix())
        if dvfield:
            kenzo_targets = EclObject([self._elements.index(edge[1])+1 for edge in dvfield])
            kenzo_sources = EclObject([self._elements.index(edge[0])+1 for edge in dvfield])
            matrices = __h_regular_dif_dvf_aux__(kenzo_stong, kenzo_targets, kenzo_sources)
        else:
            matrices = __h_regular_dif__(kenzo_stong)
        if deg is not None:
            if deg == height - 1:
                M1 = __copier_matrice__(kenzo.__nth__(height-1, matrices))
                return quotient_group_matrices(M1, right_null=True)
            else:
                M1 = __copier_matrice__(kenzo.__nth__(deg, matrices))
                M2 = __copier_matrice__(kenzo.__nth__(deg+1, matrices))
                return quotient_group_matrices(M1, M2, check=False)
        else:
            result = {}
            for dim in range(0, height - 1):
                M1 = __copier_matrice__(kenzo.__nth__(dim, matrices))
                M2 = __copier_matrice__(kenzo.__nth__(dim+1, matrices))
                result[dim] = quotient_group_matrices(M1, M2, check=False)
            Mh = __copier_matrice__(kenzo.__nth__(height-1, matrices))
            result[height-1] = quotient_group_matrices(Mh, right_null=True)
            return result
                   
