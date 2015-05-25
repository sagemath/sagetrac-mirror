r"""
SkewTableau Element class and subclasses.

See SkewTableaux and its subclasses for the corresponding Parent
classes.

AUTHORS:

- Mike Hansen: Initial version
- Travis Scrimshaw, Arthur Lubovsky (2013-02-11): Factored out
  ``CombinatorialClass``
- Josh Swanson (and others) (2015): tableaux refactoring/cleanup

CLASSES:

.. TODO:: List classes.
"""
#*****************************************************************************
#       Copyright (C) 2015 Josh Swanson,
#                     2015 Jan Keitel
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
#*****************************************************************************

import six

from sage.combinat.tableaux.bad_shape_tableau import BadShapeTableau

class SkewTableau(BadShapeTableau):
    r"""
    A tableau of skew shape.

    If `\lambda` and `\mu` are two partitions such that
    `\mu \subseteq \lambda`, then a tableau of skew shape
    `\lambda / \mu` means a map that sends the cells of this
    skew shape to some objects. This map "knows" the
    partitions `\lambda` and `\mu`.

    Cell locations are pairs of non-negative integers, though values are
    unrestricted. Cells must form a skew shape.
    """
    _generic_parent = parent_class('SkewTableaux')

    def __init__(self, parent, skp):
        """
        TESTS::

            sage: st = SkewTableau([[None, 1],[2,3]])
            sage: st = SkewTableau([[None,1,1],[None,2],[4]])
            sage: TestSuite(st).run()

        A skew tableau is immutable, see :trac:`15862`::

            sage: T = SkewTableau([[None,2],[2]])
            sage: t0 = T[0]
            sage: t0[1] = 3
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
            sage: T[0][1] = 5
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
        """
        try:
            skp = tuple(tuple(row) for row in skp)
        except TypeError:
            raise TypeError("each element of the skew tableau must be an iterable")
        self._skp = skp

        AbstractTableau.__init__(self, parent)

    @cached_method
    def _dict_unsafe(self):
        d = {(i, j): k for i, row in enumerate(self._tuple())
                       for j, k in enumerate(row)}
        return d

    def _repr_(self):
        r"""
        Return the representation string.

        OUTPUT:

        A string.
        """
        return repr(self.to_list())

    def to_tuple(self):
        r"""
        Return the tuple of tuples representing the entries of
        ``self``.
        """
        return self._skp

    def to_list(self):
        r"""
        Return the list of lists representing the entries of
        ``self``.

        This is a copy of the defining data.

        OUTPUT:

        A list of lists.
        """
        return [list(i) for i in self._skp]
