"""
Unpickling methods

Python saves objects by providing a pair ``(f, data)`` such that ``f(data)``
reconstructs the object. This module collects the loading (_unpickling_ in
Python terminology) functions for Sage's matroids.

.. NOTE::

    The reason this code was separated out from the classes was to make it
    play nice with lazy importing of the ``Matroid()`` and ``matroids``
    keywords.

AUTHORS:

- Rudi Pendavingh, Stefan van Zwam (2013-07-01): initial version
"""
#*****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

include 'sage/data_structures/bitset.pxi'
import sage.matroids.matroid
import sage.matroids.basis_exchange_matroid
from .minor_matroid import MinorMatroid
from .dual_matroid import DualMatroid
from .circuit_closures_matroid cimport CircuitClosuresMatroid
from .basis_matroid cimport BasisMatroid
from .linear_matroid cimport LinearMatroid, RegularMatroid, BinaryMatroid, TernaryMatroid, QuaternaryMatroid
from .lean_matrix cimport GenericMatrix, BinaryMatrix, TernaryMatrix, QuaternaryMatrix, IntegerMatrix


#############################################################################
# BasisMatroid
#############################################################################

def unpickle_basis_matroid(version, data):
    """
    Unpickle a BasisMatroid.

    *Pickling* is Python's term for the loading and saving of objects.
    Functions like these serve to reconstruct a saved object. This all happens
    transparently through the ``load`` and ``save`` commands, and you should
    never have to call this function directly.

    INPUT:

    - ``version`` -- an integer, expected to be 0
    - ``data`` -- a tuple ``(E, R, name, BB)`` in which ``E`` is the groundset
      of the matroid, ``R`` is the rank, ``name`` is a custom name, and ``BB``
      is the bitpacked list of bases, as pickled by Sage's ``bitset_pickle``.

    OUTPUT:

    A matroid.

    .. WARNING::

        Users should never call this function directly.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = BasisMatroid(matroids.named_matroids.Vamos())
        sage: M == loads(dumps(M))  # indirect doctest
        True

    """
    cdef BasisMatroid M
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    E, R, name, BB = data
    M = BasisMatroid(groundset=E, rank=R)
    bitset_unpickle(M._bb, BB)
    M._reset_invariants()
    M.reset_current_basis()
    if name is not None:
        M.rename(name)
    return M


#############################################################################
# CircuitClosuresMatroid
#############################################################################

def unpickle_circuit_closures_matroid(version, data):
    """
    Unpickle a CircuitClosuresMatroid.

    *Pickling* is Python's term for the loading and saving of objects.
    Functions like these serve to reconstruct a saved object. This all happens
    transparently through the ``load`` and ``save`` commands, and you should
    never have to call this function directly.

    INPUT:

    - ``version`` -- an integer, expected to be 0
    - ``data`` -- a tuple ``(E, CC, name)`` in which ``E`` is the groundset
      of the matroid, ``CC`` is the dictionary of circuit closures, and
      ``name`` is a custom name.

    OUTPUT:

    A matroid.

    .. WARNING::

        Users should never call this function directly.

    EXAMPLES::

        sage: M = matroids.named_matroids.Vamos()
        sage: M == loads(dumps(M))  # indirect doctest
        True
    """
    cdef CircuitClosuresMatroid M
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    M = CircuitClosuresMatroid(groundset=data[0], circuit_closures=data[1])
    if data[2] is not None:
        M.rename(data[2])
    return M


#############################################################################
# DualMatroid
#############################################################################

def unpickle_dual_matroid(version, data):
    """
    Unpickle a DualMatroid.

    *Pickling* is Python's term for the loading and saving of objects.
    Functions like these serve to reconstruct a saved object. This all happens
    transparently through the ``load`` and ``save`` commands, and you should
    never have to call this function directly.

    INPUT:

    - ``version`` -- an integer, expected to be 0
    - ``data`` -- a tuple ``(M, name)`` in which ``M`` is
      the internal matroid, and ``name`` is a custom name.

    OUTPUT:

    A matroid.

    .. WARNING::

        Users should not call this function directly. Instead, use load/save.

    EXAMPLES::

        sage: M = matroids.named_matroids.Vamos().dual()
        sage: M == loads(dumps(M))  # indirect doctest
        True
    """
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    M = DualMatroid(data[0])
    if data[1] is not None:
        M.rename(data[1])
    return M


#############################################################################
# LeanMatrix subclasses
#############################################################################

def unpickle_generic_matrix(version, data):
    """
    Reconstruct a ``GenericMatrix`` object (internal Sage data structure).

    .. WARNING::

        Users should not call this method directly.

    EXAMPLES::

        sage: from sage.matroids.lean_matrix import *
        sage: A = GenericMatrix(2, 5, ring=QQ)
        sage: A == loads(dumps(A))  # indirect doctest
        True
    """
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    cdef GenericMatrix A = GenericMatrix(0, 0, ring=data[2])
    A._entries = data[3][:]
    A._nrows = data[0]
    A._ncols = data[1]
    return A


def unpickle_binary_matrix(version, data):
    """
    Reconstruct a ``BinaryMatrix`` object (internal Sage data structure).

    .. WARNING::

        Users should not call this method directly.

    EXAMPLES::

        sage: from sage.matroids.lean_matrix import *
        sage: A = BinaryMatrix(2, 5)
        sage: A == loads(dumps(A))  # indirect doctest
        True
        sage: C = BinaryMatrix(2, 2, Matrix(GF(2), [[1, 1], [0, 1]]))
        sage: C == loads(dumps(C))
        True
    """
    cdef BinaryMatrix A
    cdef long i
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    nrows, ncols, versionB, size, limbs, longsize, M = data
    A = BinaryMatrix(nrows, ncols)
    for i from 0 <= i < nrows:
        bitset_unpickle(A._M[i], (versionB, size, limbs, longsize, M[i]))
    return A


def unpickle_ternary_matrix(version, data):
    """
    Reconstruct a ``TernaryMatrix`` object (internal Sage data structure).

    .. WARNING::

        Users should not call this method directly.

    EXAMPLES::

        sage: from sage.matroids.lean_matrix import *
        sage: A = TernaryMatrix(2, 5)
        sage: A == loads(dumps(A))  # indirect doctest
        True
        sage: C = TernaryMatrix(2, 2, Matrix(GF(3), [[1, 1], [0, 1]]))
        sage: C == loads(dumps(C))
        True
    """
    cdef TernaryMatrix A
    cdef long i
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    nrows, ncols, versionB, size, limbs, longsize, M0, M1 = data
    A = TernaryMatrix(nrows, ncols)
    for i from 0 <= i < nrows:
        bitset_unpickle(A._M0[i], (versionB, size, limbs, longsize, M0[i]))
        bitset_unpickle(A._M1[i], (versionB, size, limbs, longsize, M1[i]))
    return A


def unpickle_quaternary_matrix(version, data):
    """
    Reconstruct a ``QuaternaryMatrix`` object (internal Sage data structure).

    .. WARNING::

        Users should not call this method directly.

    EXAMPLES::

        sage: from sage.matroids.lean_matrix import *
        sage: A = QuaternaryMatrix(2, 5, ring=GF(4, 'x'))
        sage: A == loads(dumps(A))  # indirect doctest
        True
        sage: C = QuaternaryMatrix(2, 2, Matrix(GF(4, 'x'), [[1, 1], [0, 1]]))
        sage: C == loads(dumps(C))
        True
    """
    cdef QuaternaryMatrix A
    cdef long i
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    nrows, ncols, ring, versionB, size, limbs, longsize, M0, M1 = data
    A = QuaternaryMatrix(nrows, ncols, ring=ring)
    for i from 0 <= i < nrows:
        bitset_unpickle(A._M0[i], (versionB, size, limbs, longsize, M0[i]))
        bitset_unpickle(A._M1[i], (versionB, size, limbs, longsize, M1[i]))
    return A


def unpickle_integer_matrix(version, data):
    """
    Reconstruct an ``IntegerMatrix`` object (internal Sage data structure).

    .. WARNING::

        Users should not call this method directly.

    EXAMPLES::

        sage: from sage.matroids.lean_matrix import *
        sage: A = IntegerMatrix(2, 5)
        sage: A == loads(dumps(A))  # indirect doctest
        True
    """
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    cdef IntegerMatrix A = IntegerMatrix(data[0], data[1])
    cdef long i
    for i from 0 <= i < A._nrows * A._ncols:
        A._entries[i] = data[2][i]
    return A


#############################################################################
# LinearMatroid and subclasses
#############################################################################

def unpickle_linear_matroid(version, data):
    """
    Unpickle a LinearMatroid.

    *Pickling* is Python's term for the loading and saving of objects.
    Functions like these serve to reconstruct a saved object. This all happens
    transparently through the ``load`` and ``save`` commands, and you should
    never have to call this function directly.

    INPUT:

    - ``version`` -- an integer (currently 0).
    - ``data`` -- a tuple ``(A, E, reduced, name)`` where ``A`` is the
      representation matrix, ``E`` is the groundset of the matroid,
      ``reduced`` is a boolean indicating whether ``A`` is a reduced matrix,
      and ``name`` is a custom name.

    OUTPUT:

    A :class:`LinearMatroid` instance.

    .. WARNING::

        Users should never call this function directly.

    EXAMPLES::

        sage: M = Matroid(Matrix(GF(7), [[1, 0, 0, 1, 1], [0, 1, 0, 1, 2],
        ....:                                               [0, 1, 1, 1, 3]]))
        sage: M == loads(dumps(M))  # indirect doctest
        True
        sage: M.rename("U35")
        sage: loads(dumps(M))
        U35
    """
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    A, gs, reduced, name = data
    if not reduced:
        M = LinearMatroid(groundset=gs, matrix=A, keep_initial_representation=True)
    else:
        M = LinearMatroid(groundset=gs, reduced_matrix=A)
    if name is not None:
        M.rename(name)
    return M


def unpickle_binary_matroid(version, data):
    """
    Unpickle a BinaryMatroid.

    *Pickling* is Python's term for the loading and saving of objects.
    Functions like these serve to reconstruct a saved object. This all happens
    transparently through the ``load`` and ``save`` commands, and you should
    never have to call this function directly.

    INPUT:

    - ``version`` -- an integer (currently 0).
    - ``data`` -- a tuple ``(A, E, B, name)`` where ``A`` is the
      representation matrix, ``E`` is the groundset of the matroid, ``B`` is
      the currently displayed basis, and ``name`` is a custom name.

      OUTPUT:

      A :class:`BinaryMatroid` instance.

    .. WARNING::

        Users should never call this function directly.

    EXAMPLES::

        sage: M = Matroid(Matrix(GF(2), [[1, 0, 0, 1], [0, 1, 0, 1],
        ....:                            [0, 0, 1, 1]]))
        sage: M == loads(dumps(M))  # indirect doctest
        True
        sage: M.rename("U34")
        sage: loads(dumps(M))
        U34
    """
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    A, gs, basis, name = data
    if basis is None:
        M = BinaryMatroid(groundset=gs, matrix=A, keep_initial_representation=True)
    else:
        M = BinaryMatroid(groundset=gs, matrix=A, basis=basis)
    if name is not None:
        M.rename(name)
    return M


def unpickle_ternary_matroid(version, data):
    """
    Unpickle a TernaryMatroid.

    *Pickling* is Python's term for the loading and saving of objects.
    Functions like these serve to reconstruct a saved object. This all happens
    transparently through the ``load`` and ``save`` commands, and you should
    never have to call this function directly.

    INPUT:

    - ``version`` -- an integer (currently 0).
    - ``data`` -- a tuple ``(A, E, B, name)`` where ``A`` is the
      representation matrix, ``E`` is the groundset of the matroid, ``B`` is
      the currently displayed basis, and ``name`` is a custom name.

    OUTPUT:

    A :class:`TernaryMatroid` instance.

    .. WARNING::

        Users should never call this function directly.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = TernaryMatroid(Matrix(GF(3), [[1, 0, 0, 1], [0, 1, 0, 1],
        ....:           [0, 0, 1, 1]]))
        sage: M == loads(dumps(M))  # indirect doctest
        True
        sage: M.rename("U34")
        sage: loads(dumps(M))
        U34
    """
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    A, gs, basis, name = data
    if basis is None:
        M = TernaryMatroid(groundset=gs, matrix=A, keep_initial_representation=True)
    else:
        M = TernaryMatroid(groundset=gs, matrix=A, basis=basis)
    if name is not None:
        M.rename(name)
    return M


def unpickle_quaternary_matroid(version, data):
    """
    Unpickle a QuaternaryMatroid.

    *Pickling* is Python's term for the loading and saving of objects.
    Functions like these serve to reconstruct a saved object. This all happens
    transparently through the ``load`` and ``save`` commands, and you should
    never have to call this function directly.

    INPUT:

    - ``version`` -- an integer (currently 0).
    - ``data`` -- a tuple ``(A, E, B, name)`` where ``A`` is the
      representation matrix, ``E`` is the groundset of the matroid, ``B`` is
      the currently displayed basis, and ``name`` is a custom name.

    OUTPUT:

    A :class:`TernaryMatroid` instance.

    .. WARNING::

        Users should never call this function directly.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = QuaternaryMatroid(Matrix(GF(3), [[1, 0, 0, 1], [0, 1, 0, 1],
        ....:          [0, 0, 1, 1]]))
        sage: M == loads(dumps(M))  # indirect doctest
        True
        sage: M.rename("U34")
        sage: loads(dumps(M))
        U34
        sage: M = QuaternaryMatroid(Matrix(GF(4, 'x'), [[1, 0, 1],
        ....:                                           [1, 0, 1]]))
        sage: loads(dumps(M)).representation()
        [1 0 1]
        [1 0 1]
    """
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    A, gs, basis, name = data
    if basis is None:
        M = QuaternaryMatroid(groundset=gs, matrix=A, keep_initial_representation=True)
    else:
        M = QuaternaryMatroid(groundset=gs, matrix=A, basis=basis)
    if name is not None:
        M.rename(name)
    return M


def unpickle_regular_matroid(version, data):
    """
    Unpickle a RegularMatroid.

    *Pickling* is Python's term for the loading and saving of objects.
    Functions like these serve to reconstruct a saved object. This all happens
    transparently through the ``load`` and ``save`` commands, and you should
    never have to call this function directly.

    INPUT:

    - ``version`` -- an integer (currently 0).
    - ``data`` -- a tuple ``(A, E, reduced, name)`` where ``A`` is the
      representation matrix, ``E`` is the groundset of the matroid,
      ``reduced`` is a boolean indicating whether ``A`` is a reduced matrix,
      and ``name`` is a custom name.

    OUTPUT:

    A :class:`RegularMatroid` instance.

    .. WARNING::

        Users should never call this function directly.

    EXAMPLES::

        sage: M = matroids.named_matroids.R10()
        sage: M == loads(dumps(M))  # indirect doctest
        True
        sage: M.rename("R_{10}")
        sage: loads(dumps(M))
        R_{10}
    """
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    A, gs, reduced, name = data
    if not reduced:
        M = RegularMatroid(groundset=gs, matrix=A, keep_initial_representation=True)
    else:
        M = RegularMatroid(groundset=gs, reduced_matrix=A)
    if name is not None:
        M.rename(name)
    return M


#############################################################################
# Minor matroids
#############################################################################

def unpickle_minor_matroid(version, data):
    """
    Unpickle a MinorMatroid.

    *Pickling* is Python's term for the loading and saving of objects.
    Functions like these serve to reconstruct a saved object. This all happens
    transparently through the ``load`` and ``save`` commands, and you should
    never have to call this function directly.

    INPUT:

    - ``version`` -- an integer, currently `0`.
    - ``data`` -- a tuple ``(M, C, D, name)``, where ``M`` is the original
      matroid of which the output is a minor, ``C`` is the set of
      contractions, ``D`` is the set of deletions, and ``name`` is a custom
      name.

    OUTPUT:

    A :class:`MinorMatroid` instance.

    .. WARNING::

        Users should never call this function directly.

    EXAMPLES::

        sage: M = matroids.named_matroids.Vamos().minor('abc', 'g')
        sage: M == loads(dumps(M))  # indirect doctest
        True
    """
    if version != 0:
        raise TypeError("object was created with newer version of Sage. Please upgrade.")
    M = MinorMatroid(matroid=data[0], contractions=data[1], deletions=data[2])
    if data[3] is not None:
        M.rename(data[3])
    return M
