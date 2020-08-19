r"""
Vector Partitions

AUTHORS:

- Amritanshu Prasad (2013): Initial version
- D\. K\. Sunko (2020): Added option 'algorithm'

"""
#*****************************************************************************
#       Copyright (C) 2013 Amritanshu Prasad <amri@imsc.res.in>
#       Copyright (C) 2020 Denis Sunko <dks@phy.hr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import warnings

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.combinat import CombinatorialElement
from sage.combinat.fast_vector_partitions import fast_vector_partitions
from sage.combinat.partition import Partition


def find_min(vect):
    """
    Return a string of ``0``'s with one ``1`` at the location where the list
    vect has its last entry which is not equal to ``0``.

    INPUT:

    - ``vec`` -- A list of integers

    OUTPUT:

    A list of the same length with ``0``'s everywhere, except for a ``1``
    at the last position where ``vec`` has an entry not equal to ``0``.

    EXAMPLES::

        sage: from sage.combinat.vector_partition import find_min
        sage: find_min([2, 1])
        [0, 1]
        sage: find_min([2, 1, 0])
        [0, 1, 0]
    """
    i = len(vect)
    while vect[i - 1] == 0 and i > 0:
        i = i - 1
    min = [0] * len(vect)
    if i > 0:
        min[i - 1] = 1
    return min


def IntegerVectorsIterator(vect, min=None):
    """
    Return an iterator over the list of integer vectors which are componentwise
    less than or equal to ``vect``, and lexicographically greater than or equal
    to ``min``.

    INPUT:

    - ``vect`` -- A list of non-negative integers
    - ``min`` -- A list of non-negative integers dominated elementwise by ``vect``

    OUTPUT:

    A list in lexicographic order of all integer vectors (as lists) which are
    dominated elementwise by ``vect`` and are greater than or equal to ``min`` in
    lexicographic order.

    EXAMPLES::

        sage: from sage.combinat.vector_partition import IntegerVectorsIterator
        sage: list(IntegerVectorsIterator([1, 1]))
        [[0, 0], [0, 1], [1, 0], [1, 1]]

        sage: list(IntegerVectorsIterator([1, 1], min = [1, 0]))
        [[1, 0], [1, 1]]
    """
    if not vect:
        yield []
    else:
        if min is None:
            min = [0] * len(vect)
        if vect < min:
            return
        else:
            for vec in IntegerVectorsIterator(vect[1:], min=min[1:]):
                yield [min[0]] + vec
            for j in range(min[0] + 1, vect[0] + 1):
                for vec in IntegerVectorsIterator(vect[1:]):
                    yield [j] + vec


class VectorPartition(CombinatorialElement):
    r"""
    A vector partition is a multiset of integer vectors.
    """

    @staticmethod
    def __classcall_private__(cls, vecpar):
        """
        Create a vector partition.

        EXAMPLES::

            sage: VectorPartition([[3, 2, 1], [2, 2, 1]])
            [[2, 2, 1], [3, 2, 1]]

        The parent class is the class of vector partitions of the sum of the
        vectors in ``vecpar``::

            sage: V = VectorPartition([[3, 2, 1], [2, 2, 1]])
            sage: V.parent()._vec
            (5, 4, 2)
        """
        vec = [sum([vec[i] for vec in vecpar]) for i in range(len(vecpar[0]))]
        P = VectorPartitions(vec)
        return P(vecpar)

    def __init__(self, parent, vecpar):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: elt =  VectorPartition([[3, 2, 1], [2, 2, 1]])
            sage: TestSuite(elt).run()
        """
        CombinatorialElement.__init__(self, parent, sorted(vecpar))

    def sum(self):
        """
        Return the sum vector as a list.

        EXAMPLES::

            sage: V = VectorPartition([[3, 2, 1], [2, 2, 1]])
            sage: V.sum()
            [5, 4, 2]
        """
        return list(self.parent()._vec)

    def partition_at_vertex(self, i):
        """
        Return the partition obtained by sorting the ``i``-th elements of
        the vectors in the vector partition.

        EXAMPLES::

            sage: V = VectorPartition([[1, 2, 1], [2, 4, 1]])
            sage: V.partition_at_vertex(1)
            [4, 2]
        """
        return Partition(sorted([vec[i] for vec in self._list], reverse=True))


class VectorPartitions(UniqueRepresentation, Parent):
    r"""
    Class of all vector partitions of ``vec`` with all parts greater than
    or equal to ``min`` in lexicographic order.

    A vector partition of ``vec`` is a list of vectors with non-negative
    integer entries whose sum is ``vec``.

    INPUT:

    - ``vec`` -- a list of non-negative integers.
    - ``min`` -- (optional) list of non-negative integers.
    - ``algorithm`` -- (optional) string specifying the algorithm used.
    
    OUTPUT:
    
    Vector partitions of the vector ``vec``, with parts lexicographically
    not smaller than ``min``.
    
    If ``min`` is not provided, all partitions will be generated.
    
    If ``algorithm = 'lex'``, partitions will be in lexicographic order. If
    ``algorithm = 'neglex'`` (the default), they will be in negative
    lexicographic order. If ``algorithm`` is anything else, a warning is raised
    and the default is used.

    EXAMPLES:

    If ``min`` is not specified, then the class of all vector partitions of
    ``vec`` is created::

        sage: VP = VectorPartitions([2, 2])
        sage: for vecpar in VP:
        ....:     print(vecpar)
        [[2, 2]]
        [[1, 1], [1, 1]]
        [[1, 0], [1, 2]]
        [[0, 2], [2, 0]]
        [[0, 2], [1, 0], [1, 0]]
        [[0, 1], [2, 1]]
        [[0, 1], [1, 0], [1, 1]]
        [[0, 1], [0, 1], [2, 0]]
        [[0, 1], [0, 1], [1, 0], [1, 0]]

    If ``min`` is specified, then the class consists of only those vector
    partitions whose parts are all greater than or equal to ``min`` in
    lexicographic order::

        sage: VP = VectorPartitions([2, 2], min = [1, 0], algorithm = 'lex')
        sage: for vecpar in VP:
        ....:     print(vecpar)
        [[1, 0], [1, 2]]
        [[1, 1], [1, 1]]
        [[2, 2]]
        
    Illustration of the ``algorithm`` option::
    
        sage: forw = VectorPartitions([3,3,3], algorithm = 'lex')                       
        sage: back = VectorPartitions([3,3,3])                                          
        sage: forward = list(forw)                                                      
        sage: backward = list(back)                                                     
        sage: forward == backward[::-1]                                                 
        True

    .. WARNING::

        The default algorithm has changed from 'lex' in Sage 9.1 to 'neglex' in
        Sage 9.2, which uses :func:`fast_vector_partitions` internally.
    """

    @staticmethod
    def __classcall_private__(cls, vec, min=None, algorithm='neglex'):
        r"""
        Create the class of vector partitions of ``vec`` where all parts
        are greater than or equal to the vector ``min``.

        EXAMPLES::

            sage: VP1 = VectorPartitions([2, 1])
            sage: VP2 = VectorPartitions((2, 1), min = [0, 1])
            sage: VP1 is VP2
            True
        """
        if min is None:
            min = find_min(vec)  #tuple([0 for v in vec[:-1]]+[1])
        min = tuple(min)
        vec = tuple(vec)
        return super(VectorPartitions, cls).__classcall__(
            cls, vec, min, algorithm)

    def __init__(self, vec, min, algorithm):
        r"""
        Initialize ``self``.

        TESTS::

            sage: VP = VectorPartitions([2, 2])
            sage: TestSuite(VP).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._vec = vec
        self._min = min
        self._algorithm = algorithm

    def _element_constructor_(self, vecpar):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: VP = VectorPartitions([2, 2])
            sage: elt = VP([[1, 0], [1, 2]]); elt
            [[1, 0], [1, 2]]
            sage: elt.parent() is VP
            True
        """
        return self.element_class(self, vecpar)

    Element = VectorPartition

    def _iter_neglex_(self):
        for p in fast_vector_partitions(self._vec, self._min):
            yield self.element_class(self, p)

    def _iter_lex_(self):
        for vec in IntegerVectorsIterator(
                list(self._vec), min=list(self._min)):  # choose the first part
            if tuple(vec) == self._vec:
                yield self.element_class(self, [vec])
            else:  # recursively find all possibilities for the rest of the vector partition
                for smaller_partition in VectorPartitions(
                    [x - vec[i] for i, x in enumerate(self._vec)],
                        min=vec,
                        algorithm=self._algorithm):
                    yield self.element_class(self,
                                             [vec] + list(smaller_partition))

    def __iter__(self):
        r"""
        Iterator for vector partitions.

        EXAMPLES::

            sage: VP = VectorPartitions([2, 2])
            sage: VP.cardinality()
            9
        """
        if self._algorithm == 'neglex':
            iterator = self._iter_neglex_
        elif self._algorithm == 'lex':
            iterator = self._iter_lex_
        else:
            warnings.warn("Unrecognized option algorithm = '" +
                          self._algorithm + "', using algorithm = 'neglex'")
            iterator = self.iter_neglex

        if all(coord == 0 for coord in self._vec):
            # the zero vector has only the empty partition
            yield self.element_class(self, [])
        else:
            yield from iterator()
