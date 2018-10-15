# -*- coding: utf-8 -*-
r"""
AUTHORS:

- Matthew Lancellotti (2018): Initial version
"""
#*****************************************************************************
#  Copyright (C) 2018 Matthew Lancellotti <mvlancellotti@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.sets.non_negative_integers import NonNegativeIntegers

from sage.combinat.partition import Partition, Partitions, k_rectangle_dimension_list

def is_k_shape(self, k):
    r"""
    A partition is a `k`-*shape* if its `k`-boundary has row-shape and col-shape that are partitions themselves. (Definition 2.1 of .next_within_bounds_)

    Given a :class:`Partition` ``self`` and a natural number ``k``, returns ``True`` if and only if ``self`` is a `k`-shape.

    Given a :class:`Partition` ``self`` *only*, returns ``True`` if and only if there exists some `k \in [1, n-1]` such that ``self`` is a `k`-shape.

    EXAMPLES::

        sage: from sage.combinat.k_shape import is_k_shape
        sage: is_k_shape(Partition([3, 1]), 1)
        False
        sage: is_k_shape(Partition([3, 1]), 2)
        True
    """
    ptn = _Partitions(self)
    if k is None:
        # see if it's a k-shape for any k in [1, n-1].
        # (note that every partition is a 0-shape and an n-shape)
        n = ptn.size()
        lis = [is_k_shape(ptn, kk) for kk in range(1, n)]
        return any(lis)
    else:
        k_bdy = ptn.k_boundary(k)
        return k_bdy.is_linked()


class kShape(Partition):


    def __init__(self, p, k):
        # ideally we move k = NatNum(k) and p = Ptn(p) below into classcall_private
        assert k in NonNegativeIntegers()
        self.k = k
        assert is_k_shape(p, self.k)
        p = _Partitions(p)
        # Partition.__init__(self, Partitions, p)
        super(Partition, self)

    def h_bounds(self, width):
        r"""
        Recall the `H_i` as defined in Definition 3.3 of [HM2011]_.

        Given a natural number ``k`` (used for the `k`-shape or `k`-boundary) and a width ``width``, returns `(y_\text{min}, y_\text{max})`, the two vertical coordinates which define the horizontal strip.

        EXAMPLES:

        The 4-boundary of partition (10, 7, 4, 2, 2, 2, 1, 1, 1, 1) is shown below on a cartesian plane with the vertical lines corresponding to the vertical bounds shown in blue.

        .. image:: _static/h-v-bounds.png
            :height: 240px
            :align: center
            :alt: The skew-shape (10, 7, 4, 2, 2, 2, 1, 1, 1, 1) / (7, 4, 2, 1, 1, 1) is shown on a cartesian plane with the vertical lines y = 0, y = 1, y = 2, and y = 10 labelled.

        ::

            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).h_bounds(width=3)
            (0, 2)
            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).h_bounds(width=2)
            (2, 3)
            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).h_bounds(width=1)
            (3, 10)

        ::

            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).h_bounds(width=99)
            (0, 0)
            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).h_bounds(width=0)
            Traceback (most recent call last):
            ...
            ValueError: min() arg is an empty sequence

        ..  SEEALSO::

            :meth:`v_bounds`
        """
        k = self.k
        r = self.k_row_lengths(k)
        # pad with a row of infinite length and a row of length 0
        r = [float('inf')] + r + [0]
        y_min = max([j for j in range(0, len(r)) if r[j] > width])
        y_max = min([j for j in range(0, len(r)) if r[j] < width]) - 1
        return (y_min, y_max)

    def v_bounds(self, height):
        r"""
        This is `V_i`, the vertical analog of :meth:`h_bounds`.

        EXAMPLES:

        The 4-boundary of partition (10, 7, 4, 2, 2, 2, 1, 1, 1, 1) is shown below on a cartesian plane with the vertical lines corresponding to the vertical bounds shown in blue.

        .. image:: _static/h-v-bounds.png
            :height: 240px
            :align: center
            :alt: The skew-shape (10, 7, 4, 2, 2, 2, 1, 1, 1, 1) / (7, 4, 2, 1, 1, 1) is shown on a cartesian plane with the vertical lines y = 0, y = 1, y = 2, and y = 10 labelled.

        ::

            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).v_bounds(width=4)
            (0, 1)
            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).v_bounds(width=3)
            (1, 2)
            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).v_bounds(width=2)
            (2, 2)
            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).v_bounds(width=1)
            (2, 10)

        ::

            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).v_bounds(width=99)
            (0, 0)
            sage: kShape([10, 7, 4, 2, 2, 2, 1, 1, 1, 1], k=4).v_bounds(width=0)
            Traceback (most recent call last):
            ...
            ValueError: min() arg is an empty sequence

        ..  SEEALSO::

            :meth:`h_bounds`
        """
        return self.h_bounds(self.conjugate(), height)

    def is_k_reducible_by_rectangle(self, height, width):
        r"""
        Checks if this `k`-shape is `k`-reducible for a `k`-rectangle of dimensions ``height`` x ``width``.

        See Proposition 3.8 in Combinatorics of k-shapes and Genocchi numbers.

        INPUTS:

        - ``self`` -- a kShape

        - ``height`` -- the height of the rectangle.

        - ``width`` -- the width of the rectangle.

        EXAMPLES:

            sage: Partition([1]).is_k_reducible_by_rectangle(1, (1,1))
            sage: True
            sage: Partition([2, 1]).is_k_reducible_by_rectangle(1, (1,1))
            sage: True
            sage: Partition([1, 1]).is_k_reducible_by_rectangle(2, (1,1))
            sage: False
            sage: Partition([1, 1]).is_k_reducible_by_rectangle(2, (1,2))
            sage: True
            sage: Partition([1, 1]).is_k_reducible_by_rectangle(2, (2,1))
            sage: False

        TESTS::

			sage: p = kShape([3, 2, 1], k=3)
			sage: p.is_k_reducible_by_rectangle(3, 1)
			False

			sage: p = kShape([3, 2, 1], k=3)
			sage: p.is_k_reducible_by_rectangle(2, 2)
			True

			sage: p = kShape([3, 2, 1], k=3)
			sage: p.is_k_reducible_by_rectangle(1, 3)
			False

			sage: p = kShape([3, 2, 1], k=3)
			sage: p.is_k_reducible_by_rectangle(2, 1)
			False

			sage: p = kShape([3, 2, 1], k=3)
			sage: p.is_k_reducible_by_rectangle(1, 2)
			False

        ..  SEEALSO::

            :meth:`is_reducible`
        """
        k = self.k
        h = height
        w = width
        assert h + w - 1 == k or h + w - 1 == k - 1
        # get intersection H_a \cap V_b \cap k_rim
        rim = self.k_rim(k)
        (y_min, y_max) = self.h_bounds(h)
        (x_min, x_max) = self.v_bounds(w)
        intersection_rim = [(x, y) for (x, y) in rim
                            if x_min <= x <= x_max and y_min <= y <= y_max]
        # check condition (iii) of Proposition 3.8
        if not intersection_rim:
            return False
        else:
            # min_y is DIFFERENT than y_min
            min_y = intersection_rim[0][1]
            max_y = intersection_rim[-1][1]
            return max_y - min_y >= w

    def is_reducible(self):
        r"""
        Return ``True`` if and only if ``self`` is reducible.

        A `k`-shape ``self`` is called *reducible* if there exists a `k`- or `k-1`-rectangle corresponding to both the `k`-row-shape and `k`-column-shape of `self`.

        For a more rigorous definition, see Definition 3.7 of [HM2011]_.

        Note that this is different than the definition of a reducible partition!

        (Also, a `k`-shape is reducible if and only if it is not irreducible.)

        EXAMPLES:

        The partition [3, 2, 1] has 3-row-shape [2, 2, 1] and 3-column-shape [2, 2, 1].  It is 3-reducible because there exists a 2x2-rectangle R in the 3-row-shape and the cells that make up R when viewed in the 3-column-shape form a 2x2-rectangle (you can't see it, but the 2's are switched here)::

            sage: kShape([3, 2, 1]).is_reducible(k=3)
            True

        In this example, no good rectangle can be found::

            sage: kShape([5, 3, 2, 1, 1]).is_reducible(k=4)
            False

        TESTS:

			sage: p = kShape([1], k=1)
			sage: p.is_reducible()
			True

			sage: p = kShape([2, 1], k=1)
			sage: p.is_reducible()
			True

			sage: p = kShape([1, 1], k=2)
			sage: p.is_reducible()
			True

			sage: p = kShape([2, 1, 1], k=2)
			sage: p.is_reducible()
			True

			sage: p = kShape([2, 1, 1], k=3)
			sage: p.is_reducible()
			True

        ..  SEEALSO::

            :meth:`is_irreducible`, :meth:`k_to_irreducible_k_shapes`
        """
        k = self.k
        rect_dim_list = k_rectangle_dimension_list(
            k) + k_rectangle_dimension_list(k-1)
        for (a, b) in rect_dim_list:
            if self.is_k_reducible_by_rectangle(a, b):
                return True
        return False

    def is_irreducible(self):
        r"""
        Return ``True`` if and only if ``self`` is irreducible.

        A `k`-shape ``self`` is called *irreducible* if there does *not* exist a `k`- or `k-1`-rectangle corresponding to both the `k`-row-shape and `k`-column-shape of ``self``.

        For a more rigorous definition, see Definition 3.7 of [HM2011]_.

        (Also, a `k`-shape is irreducible if and only if it is not reducible.)

        EXAMPLES:

        The partition [3, 2, 1] has 3-row-shape [2, 2, 1] and 3-column-shape [2, 2, 1].  It is not 3-irreducible because there exists a 2x2-rectangle R in the 3-row-shape and the cells that make up R when viewed in the 3-column-shape form a 2x2-rectangle (you can't see it, but the 2's are switched here)::

            sage: kShape([3, 2, 1]).is_irreducible(k=3)
            False

        In this example, no good rectangle can be found, making it irreducible::

            sage: kShape([5, 3, 2, 1, 1]).is_irreducible(k=4)
            True

        ..  SEEALSO::

            :meth:`is_reducible`, :meth:`k_to_irreducible_k_shapes`
        """
        return not self.is_reducible()

############# GETTER FUNCS ############
def k_to_irreducible_k_shapes(k):
    r"""
    Given a natural number ``k``, return a list of all irreducible ``k``-shapes.

    Note that the algorithm runs very slowly after `k=4` :(.

    EXAMPLES::

        sage: from sage.combinat.k_shape import k_to_irreducible_k_shapes
        sage: k_to_irreducible_k_shapes(3)
        [[], [1], [2, 1]]

    TESTS::

        sage: k_to_irreducible_k_shapes(1)
        [[]]

        sage: k_to_irreducible_k_shapes(2)
        [[]]

    ..  SEEALSO::

        :meth:`is_reducible`, :meth:`is_irreducible`
    """
    bound = (k-1)*k/2
    n_bound = bound**2
    ptns = []
    for n in range(0, n_bound+1):
        ptns += Partitions(n, max_length=bound, max_part=bound)
        k_irr_k_shapes = [kShape(p, k=k) for p in ptns
                          if is_k_shape(p, k) and kShape(p, k=k).is_irreducible()]
    return k_irr_k_shapes
