r"""
Elements of `\Q/n\Z`.

EXAMPLES::

    sage: A = QQ / (3*ZZ)
    sage: x = A(11/3); x
    2/3
    sage: x*14
    1/3
    sage: x.additive_order()
    9
    sage: x / 3
    2/9
"""

#*****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import AdditiveGroupElement
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.structure.richcmp import richcmp, op_EQ, op_NE

class QmodnZ_Element(AdditiveGroupElement):
    r"""
    The ``QmodnZ_Element`` class represents an element of the abelian group `\Q/n\Z`.

    INPUT:

    - ``q`` -- a rational number.

    - ``parent`` -- the parent abelian group `\Q/n\Z`.

    OUTPUT:

    The element `q` of abelian group `\Q/n\Z`, in standard form.

    EXAMPLES::

        sage: G = QQ/(19*ZZ)
        sage: G(400/19)
        39/19
    """
    def __init__(self, parent, x, construct=False):
        r"""
        Create an element of `\Q/n\Z`.

        EXAMPLES::

            sage: G = QQ/(3*ZZ)
            sage: G.random_element()
            47/16
        """

        AdditiveGroupElement.__init__(self, parent)
        # x = (a/b) = q(n/m) + r/mb
        # am = q(nb) + r
        # r < nb so r/mb < n/m
        if construct:
            self._x = x
            return
        n = parent.n.numerator()
        if n == 0:
            self._x = x
        else:
            m = parent.n.denominator()
            a = x.numerator()
            b = x.denominator()
            q, r = (a*m).quo_rem(n*b)
            self._x = r/(m*b)

    def lift(self):
        r"""
        Return the smallest non-negative rational number reducing to this element.

        EXAMPLES::

            sage: G = QQ/(5*ZZ)
            sage: g = G(2/4); g
            1/2
            sage: q = lift(g); q
            1/2

        TESTS::

            sage: q.parent() is QQ
            True
        """
        return self._x

    def __neg__(self):
        r"""
        Return the additive inverse of this element in `\Q/n\Z`.

        EXAMPLES::

            sage: from sage.groups.additive_abelian.qmodnz import QmodnZ
            sage: G = QmodnZ(5/7)
            sage: g = G(13/21)
            sage: -g
            2/21

        TESTS::

            sage: G = QmodnZ(19/23)
            sage: g = G(15/23)
            sage: -g
            4/23
            sage: g + -g == G(0)
            True
        """
        if self._x == 0:
            return self
        else:
            QZ = self.parent()
            return QZ.element_class(QZ, QZ.n - self._x, True)

    def _add_(self, other):
        r"""
        Return the sum of two elements in `\Q/n\Z`.

        EXAMPLES::

            sage: from sage.groups.additive_abelian.qmodnz import QmodnZ
            sage: G = QmodnZ(9/10)
            sage: g = G(5)
            sage: h = G(1/2)
            sage: g + h
            1/10
            sage: g + h == G(1/10)
            True

        TESTS::

            sage: h + g == G(1/10)
            True
        """
        QZ = self.parent()
        ans = self._x + other._x
        if ans >= QZ.n:
            ans -= QZ.n
        return QZ.element_class(QZ, ans, True)

    def _sub_(self, other):
        r"""
        Returns the difference of two elements in `\Q/n\Z`.

        EXAMPLES::

            sage: from sage.groups.additive_abelian.qmodnz import QmodnZ
            sage: G = QmodnZ(9/10)
            sage: g = G(4)
            sage: h = G(1/2)
            sage: g - h
            4/5
            sage: h - g
            1/10
            sage: g - h == G(4/5)
            True
            sage: h - g == G(1/10)
            True
        """
        QZ = self.parent()
        ans = self._x - other._x
        if ans < 0:
            ans += QZ.n
        return QZ.element_class(QZ, ans, True)

    def _rmul_(self, c):
        r"""
        Returns the (right) scalar product of this element by ``c`` in `\Q/n\Z`.

        EXAMPLES::

            sage: from sage.groups.additive_abelian.qmodnz import QmodnZ
            sage: G = QmodnZ(5/7)
            sage: g = G(13/21)
            sage: g*6
            1/7
        """
        QZ = self.parent()
        return QZ.element_class(QZ, self._x * c)

    def _lmul_(self, c):
        r"""
        Returns the (left) scalar product of this element by ``c`` in `\Q/n\Z`.

        EXAMPLES::

            sage: from sage.groups.additive_abelian.qmodnz import QmodnZ
            sage: G = QmodnZ(5/7)
            sage: g = G(13/21)
            sage: 6*g
            1/7

        TESTS::

            sage: 6*g == g*6
            True
            sage: 6*g == 5*g
            False
        """
        return self._rmul_(c)

    def __div__(self, other):
        r"""
        Division.

        .. WARNING::

            Division of `x` by `m` does not yield a well defined
            result, since there are `m` elements `y` of `\Q/n\Z`
            with the property that `x = my`.  We return the one
            with the smallest non-negative lift.

        EXAMPLES::

            sage: G = QQ/(4*ZZ)
            sage: x = G(3/8)
            sage: x / 4
            3/32
        """
        #TODO: This needs to be implemented.
        QZ = self.parent()
        other = ZZ(other)
        return QZ.element_class(QZ, self._x / other, True)

    def _repr_(self):
        r"""
        Display the element.

        EXAMPLES::

            sage: G = QQ/(8*ZZ)
            sage: g = G(25/7); g;
            25/7
        """
        return repr(self._x)

    def __hash__(self):
        r"""
        Hashing.

        TESTS::

            sage: G = QQ/(4*ZZ)
            sage: g = G(4/5)
            sage: hash(g)
            -7046029254386353128
            sage: hash(G(3/4))
            3938850096065010962
            sage: hash(G(1))
            1
        """
        return hash(self._x)

    def _richcmp_(self, right, op):
        r"""
        Compare two elements.

        EXAMPLES::

            sage: G = QQ/(4*ZZ)
            sage: g = G(4/5)
            sage: h = G(6/7)
            sage: g == h
            False
            sage: g == g
            True
        """
        if op == op_EQ or op == op_NE:
            return richcmp(self._x, right._x, op)
        else:
            return NotImplemented

    def additive_order(self):
        r"""
        Returns the order of this element in the abelian group `\Q/n\Z`.

        EXAMPLES::

            sage: G = QQ/(12*ZZ)
            sage: g = G(5/3)
            sage: g.additive_order()
            36
            sage: (-g).additive_order()
            36
        """
        # a/b * k = n/m * r
        QZ = self.parent()
        if QZ.n == 0:
            if self._x == 0:
                return ZZ(1)
            return infinity
        return (self._x / QZ.n).denominator()
