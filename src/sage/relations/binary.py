"""
Finite Binary Relations
"""
#*****************************************************************************
#  Copyright (C) 2017 Victor Porton <porton@narod.ru>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.sets.set import Set, Set_object_enumerated
from sage.graphs.graph import DiGraph

def BinaryRelation(points=[]):
    r"""
    Create a binary relation.

    :param points:  a set, list, or tuple of points (like `\{(1,2),(3,4)\}` or `\{[1,2],[3,4]\}`)
    :return: A binary relation

    Create a binary relation. Currently only finite sets of points are supported.

    EXAMPLES::

        sage: from sage.relations.binary import BinaryRelation
        sage: id = BinaryRelation({(0,0),(1,1),(2,2)})
        sage: id
        relation {(0, 0), (1, 1), (2, 2)}
        sage: (id.x_range, id.y_range)
        ({0, 1, 2}, {0, 1, 2})
        sage: r1 = BinaryRelation({(0,1)})
        sage: r2 = BinaryRelation({(1,2)})
        sage: r1.reverse
        relation {(1, 0)}
        sage: r2.reverse
        relation {(2, 1)}
        sage: r2.compose(r1)
        relation {(0, 2)}
        sage: id.compose(r1) == r1
        True
        sage: r1.compose(id) == r1
        True
        sage: id.is_symmetric
        True
        sage: id.is_antisymmetric
        True
        sage: id.is_transitive
        True
        sage: (r1 | r2).is_transitive
        False
        sage: id.is_reflexive(range(3))
        True
        sage: (r1 | r2).transitive_closure
        relation {(0, 1), (1, 2), (0, 2)}
        sage: id.to_digraph()
        Looped digraph on 3 vertices

        sage: from sage.relations.binary import FiniteBinaryRelation
        sage: from sage.graphs.graph import Graph
        sage: FiniteBinaryRelation.from_graph(Graph([[0,1,'x'], [1,1,'y']], loops=True))
        relation {(0, 1), (1, 0), (1, 1)}

        Set operations also work:

        sage: r1.union(r2)
        relation {(0, 1), (1, 2)}

        Product relations:

        sage: FiniteBinaryRelation.square(xrange(2)) - FiniteBinaryRelation.identity(xrange(2))
        relation {(0, 1), (1, 0)}
    """
    return FiniteBinaryRelation(points)

#################################################################
class FiniteBinaryRelation(Set_object_enumerated):
    def __init__(self, points=[]):
        points = [(p[0],p[1]) for p in points]
        super(Set_object_enumerated, self).__init__(points)
        if not self.is_finite():
            raise ValueError("FiniteBinaryRelation must be finite.")

    def _repr_(self):
        return "relation {" + ', '.join(["(%s, %s)" % (repr(p[0]), repr(p[1])) for p in self]) + "}"

    def contains(self, x, y):
        r"""
        Test if the relation contains the given point.
        """
        return (x, y) in self

    @property
    def x_range(self):
        r"""
        :return: The x range for the relation.
        """
        return Set([p[0] for p in self])

    @property
    def y_range(self):
        r"""
        :return: The y range for the relation.
        """
        return Set([p[1] for p in self])

    @property
    def reverse(self):
        r"""
        :return: The reverse relation (x and y ranges are interchanged).
        """
        return FiniteBinaryRelation(points=[(p[1],p[0]) for p in self])

    def compose(self, rel2):
        r"""
        :return: Composition of the binary relation `rel2` with `self`.

        Composition of the binary relation `rel2` with `self`. See
        https://en.wikipedia.org/wiki/Composition_of_relations
        """
        result = set()

        yz_map = {}
        for p in self:
            if p[0] in self:
                yz_map[p[0]].add(p[1])
            else:
                yz_map[p[0]] = [p[1]]

        for p in rel2:
            if p[1] in yz_map:
                result.update([(p[0], z) for z in yz_map[p[1]]])

        return FiniteBinaryRelation(result)

    # Below we convert to list() to avoid "multilevel" linked datastructures

    def union(self, other):
        return FiniteBinaryRelation(self.set().union(other.set()))

    def intersection(self, other):
        return FiniteBinaryRelation(self.set().intersection(other.set()))

    def difference(self, other):
        return FiniteBinaryRelation(self.set().difference(other.set()))

    def symmetric_difference(self, other):
        return FiniteBinaryRelation(self.set().symmetric_difference(other.set()))

    @property
    def is_symmetric(self):
        r"""
        Check if the relation is symmetric.
        """
        for p in self:
            if (p[1], p[0]) not in self:
                return False
        return True

    @property
    def is_antisymmetric(self):
        r"""
        Check if the relation is antisymmetric.
        """
        for p in self:
            if p[0] != p[1] and (p[1], p[0]) in self:
                return False
        return True

    @property
    def is_transitive(self):
        r"""
        Check if the relation is transitive.
        """
        return self.compose(self).issubset(self)  # TODO: Can be made more efficient?

    def is_reflexive(self, on=None):
        r"""
        Check if the relation is reflexive on a set.

        Do not forget to specify the range, otherwise the union of x range and y range will be used,
        what may be not what you need!
        """
        if on is None:
            on = self.x_range().union(self.y_range)
        for i in on:
            if not (i, i) in self:
                return False
        return True

    @property
    def is_preorder(self):
        r"""
        Check if the relation is a preorder.
        """
        return self.reflexive and self.transitive

    @property
    def is_strict_order(self):
        r"""
        Check if the relation is a strict partial order.
        """
        return self.is_preorder and self.is_symmetric

    @property
    def is_nonstrict_order(self):
        r"""
        Check if the relation is a nonstrict partial order.
        """
        return self.is_preorder and self.is_antisymmetric

    @property
    def transitive_closure(self):
        result = self
        while True:
            next = result | result.compose(result)
            if next == result:
                return result
            result = next

    def to_digraph(self):
        return DiGraph(self.list(), loops=True)

    @staticmethod
    def from_graph(graph):
        r"""
        Returns the adjacency relation for the (di)graph.

        TODO: Should instead be a method in GenericGraph class.
        """
        points = set()
        directed = graph.is_directed()
        for i,j,l in graph.edge_iterator():
            points.add((i,j))
            if not directed:
                points.add((j, i))
        return FiniteBinaryRelation(points)

    @staticmethod
    def identity(A):
        return FiniteBinaryRelation([(x, x) for x in A])

    @staticmethod
    def product(A, B):
        r"""
        Cartesian product relation.
        """
        return FiniteBinaryRelation([(x,y) for x in A for y in B])

    @staticmethod
    def square(A):
        r"""
        Cartesian product of a relation with itself.
        """
        return FiniteBinaryRelation.product(A, A)

    # TODO: subsets() method?