
from .matroid cimport Matroid


cdef class MatroidUnion(Matroid):
    r"""
    Matroid Union.

    The matroid union of a set of matroids `\{(E_1,I_1),\ldots,(E_n,I_n)\}` is
    a matroid `(E,I)` where `E= \bigcup_{i=1}^n E_i` and

        `I= \{\bigcup_{i=1}^n J_i | J_i \in I_i \}`.

    INPUT:

    - ``matroids`` -- a iterator of matroids.

    OUTPUT:

    A ``MatroidUnion`` instance, it's a matroid union of all matroids in ``matroids``.
    """
    def __init__(self, matroids):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: MatroidUnion([matroids.Uniform(2,4),matroids.Uniform(5,8)])
            Matroid of rank 7 on 8 elements as matroid union of
            Matroid of rank 2 on 4 elements with circuit-closures
            {2: {{0, 1, 2, 3}}}
            Matroid of rank 5 on 8 elements with circuit-closures
            {5: {{0, 1, 2, 3, 4, 5, 6, 7}}}
        """
        self.matroids = list(matroids)
        E = set()
        for M in self.matroids:
            E.update(M.groundset())
        self._groundset = frozenset(E)

    cpdef groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT:

        A set.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: M = MatroidUnion([matroids.Uniform(2,4),matroids.Uniform(5,8)])
            sage: sorted(M.groundset())
            [0, 1, 2, 3, 4, 5, 6, 7]
        """
        return self._groundset

    cpdef _rank(self, X):
        r"""
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and ``X`` may be assumed to
        have the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: M = MatroidSum([matroids.Uniform(2,4),matroids.Uniform(2,4)])
            sage: M._rank([(0,0),(1,0)])
            2
            sage: M._rank([(0,0),(0,1),(0,2),(1,0),(1,1)])
            4

        ALGORITHM:

            Matroid intersection of a matroid sum and partition matroid.

        """
        summands = []
        for e in self.matroids:
            summands.append(e.delete(e.groundset()-X))
        sum_matroid = MatroidSum(summands)
        d = {}
        for (i,x) in sum_matroid.groundset():
            if x not in d:
                d[x]=set()
            d[x].add(i)
        part_matroid = PartitionMatroid([[(i,x) for i in d[x]] for x in d])
        return len(sum_matroid._intersection_unweighted(part_matroid))

    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: MatroidUnion([matroids.Uniform(2,4),matroids.Uniform(5,8)])
            Matroid of rank 7 on 8 elements as matroid union of
            Matroid of rank 2 on 4 elements with circuit-closures
            {2: {{0, 1, 2, 3}}}
            Matroid of rank 5 on 8 elements with circuit-closures
            {5: {{0, 1, 2, 3, 4, 5, 6, 7}}}
        """
        S = Matroid._repr_(self) + " as matroid union of \n"
        for M in self.matroids:
            S = S + M._repr_() +"\n"
        return S[:-1]

cdef class MatroidSum(Matroid):
    r"""
    Matroid Sum.

    The matroid sum of a list of matroids `(E_1,I_1),\ldots,(E_n,I_n)` is a matroid
    `(E,I)` where `E= \bigsqcup_{i=1}^n E_i` and `I= \bigsqcup_{i=1}^n I_i`.

    INPUT:

    - ``matroids`` -- a iterator of matroids.

    OUTPUT:

    A ``MatroidSum`` instance, it's a matroid sum of all matroids in ``matroids``.
    """
    def __init__(self, summands):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: MatroidSum([matroids.Uniform(2,4),matroids.Uniform(5,8)])
            Matroid of rank 7 on 12 elements as matroid sum of
            Matroid of rank 2 on 4 elements with circuit-closures
            {2: {{0, 1, 2, 3}}}
            Matroid of rank 5 on 8 elements with circuit-closures
            {5: {{0, 1, 2, 3, 4, 5, 6, 7}}}
        """
        self.summands = list(summands)
        E = set()
        for i in range(len(self.summands)):
            g = self.summands[i].groundset()
            E.update(zip([i]*len(g),g))
        self._groundset = frozenset(E)

    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: MatroidSum([matroids.Uniform(2,4),matroids.Uniform(2,4)])
            Matroid of rank 4 on 8 elements as matroid sum of
            Matroid of rank 2 on 4 elements with circuit-closures
            {2: {{0, 1, 2, 3}}}
            Matroid of rank 2 on 4 elements with circuit-closures
            {2: {{0, 1, 2, 3}}}
        """
        S = Matroid._repr_(self) + " as matroid sum of \n"
        for M in self.summands:
            S = S + M._repr_() +"\n"
        return S[:-1]

    cpdef groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT:

        A set.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: M = MatroidSum([matroids.Uniform(2,4),matroids.Uniform(2,4)])
            sage: sorted(M.groundset())
            [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3)]
        """
        return self._groundset

    cpdef _rank(self, X):
        r"""
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and ``X`` may be assumed to
        have the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: M = MatroidSum([matroids.Uniform(2,4),matroids.Uniform(2,4)])
            sage: M._rank([(0,0),(1,0)])
            2
            sage: M._rank([(0,0),(0,1),(0,2),(1,0),(1,1)])
            4
        """
        partition = {}
        for (i,x) in X:
            if i not in partition:
                partition[i] = set()
            partition[i].add(x)
        rk = 0
        for i, Xi in partition.iteritems():
            rk+= self.summands[i]._rank(Xi)
        return rk

cdef class PartitionMatroid(Matroid):
    r"""
    Partition Matroid.

    Given a set of disjoint sets `S=\{S_1,\ldots,S_n\}`, the partition matroid
    on `S` is `(E,I)` where `E=\bigcup_{i=1}^n S_i` and

        `I= \{X| |X\cap S_i|\leq 1,X\subset E \}`.

    INPUT:

    - ``partition`` -- an iterator of disjoint sets.

    OUTPUT:

    A ``PartitionMatroid`` instance, it's partition matroid of the ``partition``.
    """
    def __init__(self, partition):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: PartitionMatroid([[1,2,3],[4,5,6]])
            Partition Matroid of rank 2 on 6 elements
            sage: PartitionMatroid([[1,2],[2,3]])
            Traceback (most recent call last):
            ...
            ValueError: not an iterator of disjoint sets
            sage: PartitionMatroid([])
            Partition Matroid of rank 0 on 0 elements
        """
        P2 = [set(p) for p in partition]
        if P2:
            if len(set.union(*P2)) != sum(len(pi) for pi in P2):
                raise ValueError("not an iterator of disjoint sets")
        self.p = {}
        for i, pi in enumerate(P2):
            for j in pi:
                self.p[j] = i
        E = set()
        for P in partition:
            E.update(P)
        self._groundset = frozenset(E)

    cpdef groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT:

        A set.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: M = PartitionMatroid([[1,2,3],[4,5,6]])
            sage: sorted(M.groundset())
            [1, 2, 3, 4, 5, 6]
        """
        return self._groundset

    cpdef _rank(self, X):
        r"""
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and ``X`` may be assumed to
        have the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: M = PartitionMatroid([[1,2,3],[4,5,6]])
            sage: M._rank([1,5])
            2
            sage: M._rank([1,2])
            1
        """
        return len(set(map(self.p.get, X)))

    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: from sage.matroids.union_matroid import *
            sage: PartitionMatroid([[1,2,3],[4,5,6]])
            Partition Matroid of rank 2 on 6 elements
        """
        return "Partition " + Matroid._repr_(self)
