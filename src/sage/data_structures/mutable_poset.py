r"""
Asymptotic Ring
"""
#*****************************************************************************
# Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************

import sage


# *****************************************************************************


class MutablePosetElement(sage.structure.sage_object.SageObject):
    r"""
    An element of a mutable poset.
    """
    def __init__(self, poset, value):
        r"""
        See :class:`MutablePosetElement` for details.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: MutablePosetElement(P, (1, 2))
            (1, 2)
        """
        self._poset_ = poset
        self._value_ = value
        self.predecessors = set()
        self.successors = set()


    @property
    def poset(self):
        r"""
        The poset to which the element belongs.

        INPUT:

        Nothing.

        OUTPUT:

        A poset.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: e.poset is P
            True
        """
        return self._poset_


    @property
    def value(self):
        r"""
        The value of the element.

        INPUT:

        Nothing.

        OUTPUT:

        An object.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: e.value
            (1, 2)
        """
        return self._value_


    def __repr__(self):
        r"""
        Return the representation of the element.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        This methods usually returns the representation string of its
        :meth:`value`. The only exception is if this value is
        ``None``. In this case either ``'zero'`` or ``'oo'`` is
        returned depending in the nonexistence of predecessors and
        sucessors respectively.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: repr(MutablePosetElement(P, (1, 2)))  # indirect doctest
            '(1, 2)'
            sage: repr(P._zero_)  # indirect doctest
            'zero'
            sage: repr(P._oo_)  # indirect doctest
            'oo'
        """
        if self.value is None:
            if not self.predecessors:
                return 'zero'
            if not self.successors:
                return 'oo'
        return repr(self.value)


    def __hash__(self):
        r"""
        Return the hash of the element.

        INPUT:

        Nothing.

        OUTPUT:

        A hash value.

        This returns the hash value of the element.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: hash(MutablePosetElement(P, (1, 2))) == hash((1, 2))
            True
        """
        return hash(self.value)


    def __le__(left, right):
        r"""
        Return if ``left`` is less or equal to ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        ``True`` or ``False``.

        This methods usually returns if the values of the given
        elements are less or equal. The only exception is if this
        value is ``None``. In this case the elements are considered as
        special elements: If it has no predecessors, then it is
        interpreted as an element smaller than any other, if it has no
        successors, then as larger than any other.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: z = P._zero_
            sage: oo = P._oo_
            sage: z <= e
            True
            sage: e <= oo
            True
            sage: z <= oo
            True
            sage: oo <= z
            False
            sage: oo <= e
            False
            sage: e <= z
            False
            sage: z <= z
            True
            sage: oo <= oo
            True
            sage: e <= e
            True
        """
        if left.value is None:
            if not left.predecessors:
                # zero on the left
                return True
            else:
                # oo on the left
                if right.value is None:
                    # zero or oo on the right
                    return not right.successors
                else:
                    # not zero, not oo on the right
                    return False
        if right.value is None:
            if not right.successors:
                # oo on the right
                return True
            else:
                # zero on the right
                if left.value is None:
                    # zero or oo on the left
                    return not left.predecessors
                else:
                    # not zero, not oo on the right
                    return False
        return left.value <= right.value


# *****************************************************************************


def _sort_set_by_tuple_iter_(S, T):
    r"""
    Return an iterator over ``S`` respecting the order given by ``T``.

    INPUT:

    - ``S`` -- a set (or something which supports containment test).

    - ``T`` -- a tuple (or other iterable).

    OUTPUT:

    An iterator.

    In the iterator all elements of ``T``, which are also in ``S``
    appear. The order given by ``T`` is kept.

    EXAMPLES::

        sage: from sage.data_structures.mutable_poset import _sort_set_by_tuple_iter_
        sage: tuple(_sort_set_by_tuple_iter_({3, 4, 6}, (5, 4, 1, 2, 3, 6)))
        (4, 3, 6)
    """
    return iter(ell for ell in T if ell in S)


# *****************************************************************************


class MutablePoset(sage.structure.sage_object.SageObject):
    r"""
    A mutable poset.
    """
    def __init__(self, data=None):
        r"""
        See :class:`MutablePoset` for details.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: MP()
            poset()
        """

        if data is not None:
            raise NotImplementedError

        self._zero_ = MutablePosetElement(self, None)
        self._oo_ = MutablePosetElement(self, None)
        self._zero_.successors.add(self._oo_)
        self._oo_.predecessors.add(self._zero_)
        self._elements_ = {}


    def iter_elements(self):
        r"""
        Return an iterator over all elements.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: tuple(P.iter_elements())
            ()
        """
        return self._elements_.itervalues()


    def iter_all(self):
        r"""
        Return an iterator over all elements including a smallest
        element (`0`) and a largest element (`\infty`).

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: tuple(P.iter_all())
            (oo, zero)
        """
        yield self._oo_
        for e in self._elements_.itervalues():
            yield e
        yield self._zero_


    def repr(self):
        r"""
        Return a representation of the poset.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: print MP().repr()
            poset()
        """
        s = 'poset('
        s += ', '.join(repr(element) for element in self.iter_elements())
        s += ')'
        return s


    def repr_full(self):
        r"""
        Return a representation with ordering details of the poset.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: print MP().repr_full()
            poset()
            +-- oo
            |   +-- no successors
            |   +-- predecessors: zero
            +-- zero
            |   +-- successors:   oo
            |   +-- no predecessors
        """
        sortedelements = tuple(self.iter_all())
        strings = [self.repr()]
        for element in sortedelements:
            s = '+-- ' + repr(element) + '\n'
            if element.successors:
                s += '|   +-- successors:   '
                s += ', '.join(repr(e) for e in
                               _sort_set_by_tuple_iter_(element.successors,
                                                        sortedelements))
            else:
                s += '|   +-- no successors'
            s += '\n'
            if element.predecessors:
                s += '|   +-- predecessors: '
                s += ', '.join(repr(e) for e in
                               _sort_set_by_tuple_iter_(element.predecessors,
                                                        sortedelements))
            else:
                s += '|   +-- no predecessors'
            strings.append(s)
        return '\n'.join(strings)


    __repr__ = repr


# *****************************************************************************
