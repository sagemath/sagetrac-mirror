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

from sage.sets.set import Set_object_enumerated

def BinaryRelation(x_range, y_range, points=[]):
    r"""
    :param x_range: a set, list, or tuple
    :param y_range: a set, list, or tuple
    :param points:  a set, list, or tuple of points (like `{(1,2),(3,4)}` or `{[1,2],[3,4]}`)
    :return: A binary relation

    Create a binary relation. Currently only finite ranges and finite set of points are supported.

    That points are in the specified range, is currently not tested.

    EXAMPLES::

        sage: from sage.relations.binary import BinaryRelation
        sage: id = BinaryRelation(xrange(3), xrange(3), {(0,0),(1,1),(2,2)})
        sage: id
        {(0, 0), (1, 1), (2, 2)}
        sage: r1 = BinaryRelation(xrange(3), xrange(3), {(0,1)})
        sage: r2 = BinaryRelation(xrange(3), xrange(3), {(1,2)})
        sage: r1.reverse
        {(1, 0)}
        sage: r2.reverse
        {(2, 1)}
        sage: r2.compose(r1)
        {(0, 2)}

    Set operations also work:

        sage: r1.union(r2)
        {(0, 1), (1, 2)}

    """
    return FiniteBinaryRelation(x_range, y_range, points)

#################################################################
class FiniteBinaryRelation(Set_object_enumerated):
    def __init__(self, x_range, y_range, points=[]):
        self._x_range = Set_object_enumerated(x_range)
        self._y_range = Set_object_enumerated(y_range)
        if not self._x_range.is_finite() or not self._y_range.is_finite():
            raise ValueError("FiniteBinaryRelation ranges must be finite.")
        points = [ (p[0],p[1]) for p in points ]
        super(Set_object_enumerated, self).__init__(points)
        if not self.is_finite():
            raise ValueError("FiniteBinaryRelation must be finite.")

    def contains(self, point_x, point_y):
        r"""
        Test if the relation contains the given point.
        """
        return (point_x,point_y) in self

    @property
    def x_range(self):
        r"""
        :return: The x range for the relation.
        """
        return self._x_range

    @property
    def y_range(self):
        r"""
        :return: The y range for the relation.
        """
        return self._y_range

    @property
    def reverse(self):
        r"""
        :return: The reverse relation (x and y ranges are interchanged).
        """
        return FiniteBinaryRelation(self.y_range, self.x_range, points=[ (p[1],p[0]) for p in self ])

    def compose(self, rel2):
        r"""
        :return: Composition of the binary relation rel2 with self.
        """
        if self.x_range != rel2.y_range:
            raise ValueError("y_range of composed binary relation is not equal to x_range of self.")
        result = set()
        for x in rel2.x_range:
            for y in rel2.y_range:
                for z in self.y_range:
                    if (x,y) in rel2 and (y,z) in self:
                        result.add((x,z))
        return FiniteBinaryRelation(rel2.x_range, self.y_range, result)
