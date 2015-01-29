r"""
Mutable Poset

This module provides a class representing a finite partially ordered
set (poset) for the purpose of being used as a data structure. Thus
the here introduced posets are mutable, i.e., elements can be added and
removed from a poset at any time.

To get in touch with Sage's "usual" posets, start with the page
:mod:`Posets <sage.combinat.posets.__init__>` in the reference manual.


.. _mutable_poset_intro:

Introduction
============


.. _mutable_poset_examples:

Examples
========

First Steps
-----------

We start by creating an empty poset. This is simply done by

::

    sage: from sage.data_structures.mutable_poset import MutablePoset as MP
    sage: P = MP()
    sage: P
    poset()

A poset should contain elements, thus let us add them with

::

    sage: P.add(42)
    sage: P.add(7)
    sage: P.add(13)
    sage: P.add(3)

Let us look at the poset again::

    sage: P
    poset(3, 7, 13, 42)

We see that they elements are sorted using `\leq` which exists on the
integers `\ZZ`. Since this is even a total order, we could have used
the more efficient :class:`MutableToset`. (Note that also other data
structures are suitable for totally ordered sets.)


A less boring Example
---------------------

Let us continue with a less boring example. We define the class

::

    sage: class T(tuple):
    ....:     def __le__(left, right):
    ....:         return all(l <= r for l, r in zip(left, right))

It is equipped with a `\leq`-operation which makes `a \leq b` if all
entries of `a` are at most `b`. For example, we have

::

    sage: a = T((1,1))
    sage: b = T((2,1))
    sage: c = T((1,2))
    sage: a <= b, a <= c, b <= c
    (True, True, False)

The last comparison gives ``False``, since the first entries give `2 \leq 1`.

Now, let us add such elements to a poset::

    sage: Q = MP()
    sage: Q.add(T((1, 1)))
    sage: Q.add(T((3, 3)))
    sage: Q.add(T((4, 1)))
    sage: Q.add(T((3, 2)))
    sage: Q.add(T((2, 3)))
    sage: Q.add(T((2, 2)))
    sage: Q
    poset((1, 1), (2, 2), (2, 3), (3, 2), (3, 3), (4, 1))

In the representation above, the elements are sorted topologically,
smallest first. This does not show (directly) more structural
information. We can overcome this and display a "wiring layout" by
typing::

    sage: print Q.repr_full(reverse=True)
    poset((3, 3), (2, 3), (3, 2), (2, 2), (4, 1), (1, 1))
    +-- oo
    |   +-- no successors
    |   +-- predecessors:   (3, 3), (4, 1)
    +-- (3, 3)
    |   +-- successors:   oo
    |   +-- predecessors:   (2, 3), (3, 2)
    +-- (2, 3)
    |   +-- successors:   (3, 3)
    |   +-- predecessors:   (2, 2)
    +-- (3, 2)
    |   +-- successors:   (3, 3)
    |   +-- predecessors:   (2, 2)
    +-- (2, 2)
    |   +-- successors:   (2, 3), (3, 2)
    |   +-- predecessors:   (1, 1)
    +-- (4, 1)
    |   +-- successors:   oo
    |   +-- predecessors:   (1, 1)
    +-- (1, 1)
    |   +-- successors:   (2, 2), (4, 1)
    |   +-- predecessors:   null
    +-- null
    |   +-- successors:   (1, 1)
    |   +-- no predecessors

Note that we use ``reverse=True`` to let the elements appear from
largest (on the top) to smallest (on the bottom).

If you look at the output above, you'll see two additional elements,
namely ``oo`` (`\infty`) and ``null`` (`\emptyset`). So what are these
strange animals? The answer is simple and maybe you can guess it
already. The `\infty`-element is larger than every other element,
therefore a successor of the maximal elements in the poset. Similarly,
the `\emptyset`-element is smaller than any other element, therefore a
predecessor of the poset's minimal elements. Both do not have to scare
us; they are just there and sometimes useful.


AUTHORS:

- Daniel Krenn (2015-01-21): initial version

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the Austrian Science Fund (FWF): P 24644-N26.

Classes and their Methods
=========================
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


class MutablePosetShell(sage.structure.sage_object.SageObject):
    r"""
    A shell for an element of a :class:`mutable poset <MutablePoset>`.
    """
    def __init__(self, poset, element):
        r"""
        See :class:`MutablePosetShell` for details.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: MutablePosetShell(P, (1, 2))
            (1, 2)
        """
        self._poset_ = poset
        self._element_ = element
        self._predecessors_ = set()
        self._successors_ = set()


    @property
    def poset(self):
        r"""
        The poset to which this shell belongs.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: e = MutablePosetShell(P, (1, 2))
            sage: e.poset is P
            True
        """
        return self._poset_


    @property
    def element(self):
        r"""
        The element contained in this shell.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: e = MutablePosetShell(P, (1, 2))
            sage: e.element
            (1, 2)
        """
        return self._element_


    @property
    def key(self):
        r"""
        The key of the element contained in this shell.

        The element is converted by the poset to the key.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: P = MP()
            sage: e = MutablePosetShell(P, (1, 2))
            sage: e.key
            (1, 2)
            sage: Q = MP(key=lambda k: k[0])
            sage: f = MutablePosetShell(Q, (1, 2))
            sage: f.key
            1
        """
        return self.poset.get_key(self._element_)


    def predecessors(self, reverse=False):
        r"""
        Return the predecessors of this shell.

        INPUT:

        - ``reverse`` -- (default: ``False``) if set, then returns
          successors instead.

        OUTPUT:

        A set.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: e = MutablePosetShell(P, (1, 2))
            sage: e.predecessors()
            set()
        """
        if reverse:
            return self._successors_
        return self._predecessors_


    def successors(self, reverse=False):
        r"""
        Return the successors of this shell.

        INPUT:

        - ``reverse`` -- (default: ``False``) if set, then returns
          predecessors instead.

        OUTPUT:

        A set.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: e = MutablePosetShell(P, (1, 2))
            sage: e.successors()
            set()
        """
        if reverse:
            return self._predecessors_
        return self._successors_


    def is_special(self):
        r"""
        Return if this shell contains either the null-element, i.e., the
        element smaller than any possible other element or the
        infinity-element, i.e., the element larger than any possible
        other element.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.null.is_special()
            True
            sage: P.oo.is_special()
            True
        """
        return self.element is None


    def is_null(self, reverse=False):
        r"""
        Return if this shell contains the null-element, i.e., the element
        smaller than any possible other element.

        INPUT:

        - ``reverse`` -- (default: ``False``) if set, then returns
          if the element is the largest possible.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.null.is_null()
            True
            sage: P.oo.is_null()
            False
        """
        return self.element is None and not self.predecessors(reverse)


    def is_oo(self, reverse=False):
        r"""
        Return if this shell contains the infinity-element, i.e., the element
        larger than any possible other element.

        INPUT:

        - ``reverse`` -- (default: ``False``) if set, then returns
          if the element is the smallest possible.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.null.is_oo()
            False
            sage: P.oo.is_oo()
            True
        """
        return self.element is None and not self.successors(reverse)


    def __repr__(self):
        r"""
        Return the representation of this shell.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        This methods usually returns the representation string of its
        :meth:`element`. The only exception is if this element is
        ``None``. In this case either ``'null'`` or ``'oo'`` is
        returned depending in the nonexistence of predecessors and
        sucessors respectively.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: repr(MutablePosetShell(P, (1, 2)))  # indirect doctest
            '(1, 2)'
            sage: repr(P.null)  # indirect doctest
            'null'
            sage: repr(P.oo)  # indirect doctest
            'oo'
        """
        if self.element is None:
            if not self.predecessors():
                return 'null'
            if not self.successors():
                return 'oo'
        return repr(self.element)


    def __hash__(self):
        r"""
        Return the hash of this shell.

        INPUT:

        Nothing.

        OUTPUT:

        A hash value.

        This returns the hash value of the key of the element
        contained in this shell.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: hash(MutablePosetShell(P, (1, 2))) == hash((1, 2))
            True
        """
        return hash(self.key)


    def le(left, right, reverse=False):
        r"""
        Return if ``left`` is less or equal to ``right``.

        INPUT:

        - ``left`` -- a shell.

        - ``right`` -- a shell.

        - ``reverse`` -- (default: ``False``) if set, then return
          ``right <= left`` instead.

        OUTPUT:

        ``True`` or ``False``.

        This methods usually returns if the keys of the given elements
        (in the shells) are less or equal. The only exception is if
        the element is ``None``. In this case the shell contains special
        elements: If it has no predecessors, then it is interpreted as
        an element smaller than any other, if it has no successors,
        then as larger than any other.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: e = MutablePosetShell(P, (1, 2))
            sage: z = P.null
            sage: oo = P.oo
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

        ::

            sage: z.le(e, reverse=True)
            False
            sage: e.le(oo, reverse=True)
            False
            sage: z.le(oo, reverse=True)
            False
            sage: oo.le(z, reverse=True)
            True
            sage: oo.le(e, reverse=True)
            True
            sage: e.le(z, reverse=True)
            True
            sage: z.le(z, reverse=True)
            True
            sage: oo.le(oo, reverse=True)
            True
            sage: e.le(e, reverse=True)
            True
        """
        if reverse:
            left, right = (right, left)

        if left.element is None:
            if not left.predecessors():
                # null on the left
                return True
            else:
                # oo on the left
                if right.element is None:
                    # null or oo on the right
                    return not right.successors()
                else:
                    # not null, not oo on the right
                    return False
        if right.element is None:
            if not right.successors():
                # oo on the right
                return True
            else:
                # null on the right
                if left.element is None:
                    # null or oo on the left
                    return not left.predecessors()
                else:
                    # not null, not oo on the right
                    return False
        return left.key <= right.key


    __le__ = le


    def eq(left, right):
        r"""
        Return if ``left`` is equal to ``right``.

        INPUT:

        - ``left`` -- a shell.

        - ``right`` -- a shell.

        OUTPUT:

        ``True`` or ``False``.

        This method compares the keys of the elements contained in the shells.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: e = MutablePosetShell(P, (1, 2))
            sage: z = P.null
            sage: oo = P.oo
            sage: z == z
            True
            sage: oo == oo
            True
            sage: e == e
            True
            sage: z == e
            False
            sage: e == oo
            False
            sage: oo == z
            False
        """
        if left.element is None and right.element is None:
            return left.is_null() == right.is_null()
        return left.key == right.key


    __eq__ = eq


    def _copy_all_linked_(self, memo, poset):
        r"""
        Helper function for :meth:`MutablePoset.copy`.

        INPUT:

        - ``memo`` -- a dictionary which assigns to the id of the
          calling shell to a copy of it.

        - ``poset`` -- the poset to which the newly created shells belongs.

        OUTPUT:

        A new shell.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: Q = MP()
            sage: memo = {}
            sage: z = P.null._copy_all_linked_(memo, Q)
            sage: z.poset is Q
            True
            sage: oo = z.successors().pop()
            sage: oo == P.oo
            True
        """
        try:
            return memo[id(self)]
        except KeyError:
            pass

        new = self.__class__(poset, self.element)
        memo[id(self)] = new

        for reverse in (False, True):
            for e in self.successors(reverse):
                new.successors(reverse).add(e._copy_all_linked_(memo, poset))

        return new


    def _search_covers_(self, covers, shell, reverse=False):
        r"""
        Helper function for :meth:`covers`.

        INPUT:

        - ``covers`` -- a set which finally contains all covers.

        - ``shell`` -- the shell for which to find the covering shells.

        - ``reverse`` -- (default: ``False``) if not set, then find
          the lower covers, otherwise find the upper covers.

        OUTPUT:

        ``True`` or ``False``.

        Note that ``False`` is returned if we do not have
        ``self <= shell``.

       TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1, 1)))
            sage: P.add(T((1, 3, 1)))
            sage: P.add(T((2, 1, 2)))
            sage: P.add(T((4, 4, 2)))
            sage: P.add(T((1, 2, 2)))
            sage: P.add(T((2, 2, 2)))
            sage: e = P.shell(T((2, 2, 2))); e
            (2, 2, 2)
            sage: covers = set()
            sage: P.null._search_covers_(covers, e)
            True
            sage: sorted(covers, key=lambda c: repr(c.element))
            [(1, 2, 2), (2, 1, 2)]
        """
        if not self.le(shell, reverse) or self == shell:
            return False
        if not any([e._search_covers_(covers, shell, reverse)
                    for e in self.successors(reverse)]):
            covers.add(self)
        return True


    def covers(self, shell, reverse=False):
        r"""
        Return the covers of the given shell (considering only
        shells which originate from itself).

        INPUT:

        - ``shell`` -- the shell for which to find the covering shells.

        - ``reverse`` -- (default: ``False``) if not set, then find
          the lower covers, otherwise find the upper covers.

        OUTPUT:

        A set of the covers.

        Suppose ``reverse`` is ``False``. This method returns all the
        lower covers of the given ``shell``, i.e., shells in the
        poset, which are at most the given shell and maximal with
        this property. Only shells which are (not necessarily
        direct) successors of the calling shell are considered.

        If ``reverse`` is ``True``, then the reverse direction is
        taken, i.e., in the text above replace lower covers by upper
        covers, maximal by minimal, and successors by predecessors.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: P.add(T((2, 2)))
            sage: e = P.shell(T((2, 2))); e
            (2, 2)
            sage: sorted(P.null.covers(e),
            ....:        key=lambda c: repr(c.element))
            [(1, 2), (2, 1)]
            sage: sorted(P.oo.covers(e, reverse=True),
            ....:        key=lambda c: repr(c.element))
            [(4, 4)]
        """
        covers = set()
        self._search_covers_(covers, shell, reverse)
        return covers


    def _iter_depth_first_visit_(self, marked, reverse=False, key=None):
        r"""
        Helper function for :meth:`iter_depth_first`.

        INPUT:

        - ``marked`` -- a set in which marked shells are stored.

        - ``reverse`` -- (default: ``False``) -- if set, reverses the order.

        - ``key`` -- (default: ``None``) a function used for sorting
          the successors. If this is ``None``, no sorting occurs.

        OUTPUT:

        An iterator.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(42)
            sage: P.add(5)
            sage: marked = set()
            sage: list(P.oo._iter_depth_first_visit_(marked, True))
            [oo, 42, 5, null]
        """
        if self in marked:
            return
        marked.add(self)
        yield self
        S = self.successors(reverse)
        if key is not None:
            S = sorted(S, key=key)
        for shell in S:
            for e in shell._iter_depth_first_visit_(marked, reverse, key):
                yield e


    def iter_depth_first(self, reverse=False, key=None):
        r"""
        Iterates over all shells in depth first order.

        INPUT:

        - ``reverse`` -- (default: ``False``) -- if set, reverses the
          order, i.e., ``False`` starts at bottom (`\emptyset`),
          ``True`` starts at top (`\infty`).

        - ``key`` -- (default: ``None``) a function used for sorting
          the direct successors of a shell (used in case of a
          tie). If this is ``None``, no sorting occurs.

        OUTPUT:

        An iterator.

        ALGORITHM:

        See :wikipedia:`Depth-first_search`.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: P.add(T((2, 2)))
            sage: list(P.null.iter_depth_first(reverse=False, key=repr))
            [null, (1, 1), (1, 2), (1, 3), (4, 4), oo, (2, 2), (2, 1)]
            sage: list(P.oo.iter_depth_first(reverse=True, key=repr))
            [oo, (4, 4), (1, 3), (1, 2), (1, 1), null, (2, 2), (2, 1)]
        """
        marked = set()
        return self._iter_depth_first_visit_(marked, reverse, key)


    def _iter_topological_visit_(self, marked, reverse=False, key=None):
        r"""
        Helper function for :meth:`iter_topological`.

        INPUT:

        - ``marked`` -- a set in which marked shells are stored.

        - ``reverse`` -- (default: ``False``) -- if set, reverses the order.

        - ``key`` -- (default: ``None``) a function used for sorting
          the direct successors of a shell (used in case of a
          tie). If this is ``None``, no sorting occurs.

        OUTPUT:

        An iterator.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(42)
            sage: P.add(5)
            sage: marked = set()
            sage: list(P.null._iter_topological_visit_(marked, True))
            [oo, 42, 5, null]
        """
        if self in marked:
            return
        marked.add(self)
        S = self.predecessors(reverse)
        if key is not None:
            S = sorted(S, key=key)
        for shell in S:
            for e in shell._iter_topological_visit_(marked, reverse, key):
                yield e
        yield self


    def iter_topological(self, reverse=False, key=None):
        r"""
        Iterates over all shells in topological order.

        INPUT:

        - ``reverse`` -- (default: ``False``) -- if set, reverses the
          order, i.e., ``False`` gives smallest shells first,
          ``True`` gives largest first.

        - ``key`` -- (default: ``None``) a function used for sorting
          the direct successors of a shell (used in case of a
          tie). If this is ``None``, no sorting occurs.

        OUTPUT:

        An iterator.

        ALGORITHM:

        Here a simplified version of the algorithm found in [T1976]_
        and [CLRS2001]_ is used. See also
        :wikipedia:`Topological_sorting`.

        .. [T1976] Robert E. Tarjan, *Edge-disjoint spanning trees and
           depth-first search*, Acta Informatica 6 (2), 1976, 171-185,
           :doi:`10.1007/BF00268499`.

        .. [CLRS2001] Thomas H. Cormen, Charles E. Leiserson, Ronald
           L. Rivest and Clifford Stein, *Section 22.4: Topological
           sort*, Introduction to Algorithms (2nd ed.), MIT Press and
           McGraw-Hill, 2001, 549-552, ISBN 0-262-03293-7.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: P.add(T((2, 2)))

        ::

            sage: for e in P.shells_topological(include_special=True,
            ....:                                 reverse=True):
            ....:     print e
            ....:     print list(e.iter_topological(reverse=True, key=repr))
            oo
            [oo]
            (4, 4)
            [oo, (4, 4)]
            (1, 3)
            [oo, (4, 4), (1, 3)]
            (2, 2)
            [oo, (4, 4), (2, 2)]
            (1, 2)
            [oo, (4, 4), (1, 3), (2, 2), (1, 2)]
            (2, 1)
            [oo, (4, 4), (2, 2), (2, 1)]
            (1, 1)
            [oo, (4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1)]
            null
            [oo, (4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1), null]

        ::

            sage: for e in P.shells_topological(include_special=True,
            ....:                                 reverse=True):
            ....:     print e
            ....:     print list(e.iter_topological(reverse=False, key=repr))
            oo
            [null, (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4), oo]
            (4, 4)
            [null, (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4)]
            (1, 3)
            [null, (1, 1), (1, 2), (1, 3)]
            (2, 2)
            [null, (1, 1), (1, 2), (2, 1), (2, 2)]
            (1, 2)
            [null, (1, 1), (1, 2)]
            (2, 1)
            [null, (1, 1), (2, 1)]
            (1, 1)
            [null, (1, 1)]
            null
            [null]
        """
        marked = set()
        return self._iter_topological_visit_(marked, reverse, key)


# *****************************************************************************


def sorted_set_by_tuple(S, T):
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

        sage: from sage.data_structures.mutable_poset import sorted_set_by_tuple
        sage: tuple(sorted_set_by_tuple({3, 4, 6}, (5, 4, 1, 2, 3, 6)))
        (4, 3, 6)
    """
    return iter(ell for ell in T if ell in S)


# *****************************************************************************


def is_MutablePoset(P):
    r"""
    Tests if ``P`` inherits from :class:`MutablePoset`.

    TESTS::

        sage: from sage.data_structures.mutable_poset import MutablePoset as MP
        sage: from sage.data_structures.mutable_poset import is_MutablePoset
        sage: P = MP()
        sage: is_MutablePoset(P)
        True
    """
    return isinstance(P, MutablePoset)


class MutablePoset(sage.structure.sage_object.SageObject):
    r"""
    A mutable poset (partially ordered set) as data structure.

    INPUT:

    - ``element_exists_hook`` -- a function. It is called during
      :meth:`add` when an element (more precisely its key) is already
      in this poset. This function has the following properties:

      - It gets as a first parameter the existing element mentioned above.

      - As a second parameter it gets the element from :meth:`add`.

      - It should return the element which to put in the poset
      instead of the existing one. If this return value is ``None``,
      the existing element is removed out of the poset.

      If ``element_exists_hook`` is ``None`` (default) the
      (existing element) is not changed, i.e., this is equivalent
      to ``element_exists_hook`` returns the first parameter. Note
      that it is not allowed that the key of this new element
      differs from the key of existing.


    .. TODO::

        Implement the following methods of
        :class:`~sage.combinat.posets.posets.FinitePoset`:

        - :meth:`~sage.combinat.posets.posets.FinitePoset.bottom`: Returns the bottom element of the poset, if it exists.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.cardinality`: Returns the number of elements in the poset.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.compare_elements`: Compares x and y in the poset.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.covers`: Returns True if y covers x and False otherwise.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.has_bottom`: Returns True if the poset has a unique minimal element.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.has_top`: Returns True if the poset contains a unique maximal element, and False otherwise.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.is_bounded`: Returns True if the poset contains a unique maximal element and a unique minimal element, and False otherwise.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.is_chain`: Returns True if the poset is totally ordered, and False otherwise.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.is_gequal`: Returns True if x is greater than or equal to y in the poset, and False otherwise.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.is_greater_than`: Returns True if x is greater than but not equal to y in the poset, and False otherwise.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.is_lequal`: Returns True if x is less than or equal to y in the poset, and False otherwise.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.is_less_than`: Returns True if x is less than but not equal to y in the poset, and False otherwise.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.linear_extension`: Returns a linear extension of this poset.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.list`: List the elements of the poset. This just returns the result of linear_extension().
        - :meth:`~sage.combinat.posets.posets.FinitePoset.lower_covers_iterator`: Returns an iterator for the lower covers of the element y. An lower cover of y is an element x such that y x is a cover relation.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.lower_covers`: Returns a list of lower covers of the element y. An lower cover of y is an element x such that y x is a cover relation.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.maximal_elements`: Returns a list of the maximal elements of the poset.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.minimal_elements`: Returns a list of the minimal elements of the poset.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.relations_iterator`: Returns an iterator for all the relations of the poset.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.relations`: Returns a list of all relations of the poset.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.top`: Returns the top element of the poset, if it exists.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.upper_covers_iterator`: Returns an iterator for the upper covers of the element y. An upper cover of y is an element x such that y x is a cover relation.
        - :meth:`~sage.combinat.posets.posets.FinitePoset.upper_covers`: Returns a list of upper covers of the element y. An upper cover of y is an element x such that y x is a cover relation.
    """
    def __init__(self, data=None, key=None, element_exists_hook=None):
        r"""
        See :class:`MutablePoset` for details.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: MP()
            poset()

        ::

            sage: P = MP()
            sage: P.add(42)
            sage: MP(P)
            poset(42)

        ::

            sage: MP([3, 5, 7])
            poset(3, 5, 7)

        ::

            sage: MP(33)
            Traceback (most recent call last):
            ...
            TypeError: 33 is not iterable; do not know what to do with it.
        """
        if is_MutablePoset(data):
            if key is not None:
                raise TypeError('Cannot use key when data is a poset.')
            self._copy_shells_(data)

        else:
            self.clear()

            if key is None:
                self._key_ = lambda k: k
            else:
                self._key_ = key

            self._element_exists_hook_ = element_exists_hook

            if data is not None:
                try:
                    it = iter(data)
                except TypeError:
                    raise TypeError('%s is not iterable; do not know what to '
                                    'do with it.' % (data,))
                self.union_update(it)


    def clear(self):
        r"""
        Remove all elements from this poset.

        INPUT:

        Nothing.

        OUTPUT:

        Nothing.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(42); P
            poset(42)
            sage: P.clear()
            sage: print P.repr_full()
            poset()
            +-- null
            |   +-- no predecessors
            |   +-- successors:   oo
            +-- oo
            |   +-- predecessors:   null
            |   +-- no successors
        """
        self._null_ = MutablePosetShell(self, None)
        self._oo_ = MutablePosetShell(self, None)
        self._null_.successors().add(self._oo_)
        self._oo_.predecessors().add(self._null_)
        self._shells_ = {}


    @property
    def null(self):
        r"""
        The shell `\emptyset` whose element is smaller than any other element.

        EXAMPLES:

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: z = P.null; z
            null
            sage: z.is_null()
            True
        """
        return self._null_


    @property
    def oo(self):
        r"""
        The shell `\infty` whose element is larger than any other element.

        EXAMPLES:

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: oo = P.oo; oo
            oo
            sage: oo.is_oo()
            True
        """
        return self._oo_


    def shell(self, key):
        r"""
        Return the shell of the element corresponding to ``key``.

        INPUT:

        ``key`` -- the key of an object.

        OUTPUT:

        An instance of :class:`MutablePosetShell`.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(42)
            sage: e = P.shell(42); e
            42
            sage: type(e)
            <class 'sage.data_structures.mutable_poset.MutablePosetShell'>
        """
        return self._shells_[key]


    def get_key(self, element):
        r"""
        Return the key corresponding to the given element.

        INPUT:

        - ``element`` -- an object.

        OUTPUT:

        An object (the key of ``element``).

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: from sage.data_structures.mutable_poset import MutablePosetShell
            sage: P = MP()
            sage: P.get_key(None) is None
            True
            sage: P.get_key((1, 2))
            (1, 2)
            sage: Q = MP(key=lambda k: k[0])
            sage: Q.get_key((1, 2))
            1
        """
        if element is None:
            return None
        return self._key_(element)


    def _copy_shells_(self, other):
        r"""
        Helper function for copying shells.

        INPUT:

        - ``other`` -- the mutable poset from which the shells
          should be copied this poset.

        OUTPUT:

        Nothing.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: Q = MP()
            sage: Q._copy_shells_(P)
            sage: P.repr_full() == Q.repr_full()
            True
        """
        from copy import copy
        self._key_ = copy(other._key_)
        memo = {}
        self._null_ = other._null_._copy_all_linked_(memo, self)
        self._oo_ = memo[id(other._oo_)]
        self._shells_ = dict((f.key, f) for f in
                              iter(memo[id(e)]
                                   for e in other._shells_.itervalues()))


    def copy(self):
        r"""
        Creates a shallow copy.

        INPUT:

        Nothing.

        OUTPUT:

        A poset with the same content as ``self``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: Q = copy(P)  # indirect doctest
            sage: P.repr_full() == Q.repr_full()
            True
        """
        new = self.__class__()
        new._copy_shells_(self)
        return new


    __copy__ = copy


    def shells(self, include_special=False, reverse=False):
        r"""
        Return an iterator over all shells.

        INPUT:

        - ``include_special`` -- (default: ``False``) if set, then
          including shells containing a smallest element (`\emptyset`)
          and a largest element (`\infty`).

        - ``reverse`` -- (default: ``False``) if set, the order is
          reversed. This only affects the shells `\emptyset` and `\infty`.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: tuple(P.shells())
            ()
            sage: tuple(P.shells(include_special=True))
            (null, oo)
            sage: tuple(P.shells(include_special=True, reverse=True))
            (oo, null)
        """
        if include_special:
            yield self.null if not reverse else self.oo
        for e in self._shells_.itervalues():
            yield e
        if include_special:
            yield self.oo if not reverse else self.null


    def shells_topological(self, include_special=False,
                             reverse=False, key=None):
        r"""
        Return an iterator over all shells in topological order.

        INPUT:

        - ``include_special`` -- (default: ``False``) if set, then
          including shells containing a smallest element (`\emptyset`)
          and a largest element (`\infty`).
 
        - ``reverse`` -- (default: ``False``) -- if set, reverses the
          order, i.e., ``False`` gives smallest elements first,
          ``True`` gives largest first.

        - ``key`` -- (default: ``None``) a function used for sorting
          the direct successors of a shell (used in case of a
          tie). If this is ``None``, no sorting according to the reprsentation
          string occurs.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: P.add(T((2, 2)))
            sage: list(P.shells_topological())
            [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4)]
            sage: list(P.shells_topological(reverse=True))
            [(4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1)]
            sage: list(P.shells_topological(include_special=True))
            [null, (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4), oo]
            sage: list(P.shells_topological(
            ....:     include_special=True, reverse=True))
            [oo, (4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1), null]
        """
        if key is None:
            key = repr
        shell = self.oo if not reverse else self.null
        return iter(e for e in shell.iter_topological(reverse, key)
                    if include_special or not e.is_special())


    def elements(self, **kwargs):
        r"""
        Return an iterator over all elements.

        INPUT:

        - ``kwargs`` -- arguments are passed to :meth:`shells`.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3)
            sage: P.add(42)
            sage: P.add(7)
            sage: [(v, type(v)) for v in sorted(P.elements())]
            [(3, <type 'sage.rings.integer.Integer'>),
             (7, <type 'sage.rings.integer.Integer'>),
             (42, <type 'sage.rings.integer.Integer'>)]

        Note that

        ::

            sage: it = iter(P)
            sage: sorted(it)
            [3, 7, 42]

        returns all elements as well.
        """
        for shell in self.shells(**kwargs):
            yield shell.element


    __iter__ = elements


    def elements_topological(self, **kwargs):
        r"""
        Return an iterator over all elements in topological order.

        INPUT:

        - ``kwargs`` -- arguments are passed to :meth:`shells_topological`.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: P.add(T((2, 2)))
            sage: [(v, type(v)) for v in P.elements_topological()]
            [((1, 1), <class '__main__.T'>),
             ((1, 2), <class '__main__.T'>),
             ((1, 3), <class '__main__.T'>),
             ((2, 1), <class '__main__.T'>),
             ((2, 2), <class '__main__.T'>),
             ((4, 4), <class '__main__.T'>)]
        """
        for shell in self.shells_topological(**kwargs):
            yield shell.element


    def keys(self, **kwargs):
        r"""
        Return an iterator over all keys of the elements.

        INPUT:

        - ``kwargs`` -- arguments are passed to :meth:`shells`.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP(key=lambda c: -c)
            sage: P.add(3)
            sage: P.add(42)
            sage: P.add(7)
            sage: [(v, type(v)) for v in sorted(P.keys())]
            [(-42, <type 'sage.rings.integer.Integer'>),
             (-7, <type 'sage.rings.integer.Integer'>),
             (-3, <type 'sage.rings.integer.Integer'>)]

            sage: [(v, type(v)) for v in sorted(P.elements())]
            [(3, <type 'sage.rings.integer.Integer'>),
             (7, <type 'sage.rings.integer.Integer'>),
             (42, <type 'sage.rings.integer.Integer'>)]

            sage: [(v, type(v)) for v in sorted(P.shells(),
            ....:                               key=lambda c: c.element)]
            [(3, <class 'sage.data_structures.mutable_poset.MutablePosetShell'>),
             (7, <class 'sage.data_structures.mutable_poset.MutablePosetShell'>),
             (42, <class 'sage.data_structures.mutable_poset.MutablePosetShell'>)]
        """
        for shell in self.shells(**kwargs):
            yield shell.key


    def keys_topological(self, **kwargs):
        r"""
        Return an iterator over all keys of the elements in
        topological order.

        INPUT:

        - ``kwargs`` -- arguments are passed to :meth:`shells_topological`.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP(key=lambda c: c[0])
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: P.add(T((2, 2)))
            sage: [(v, type(v)) for v in P.keys_topological()]
            [(1, <type 'sage.rings.integer.Integer'>),
             (2, <type 'sage.rings.integer.Integer'>),
             (4, <type 'sage.rings.integer.Integer'>)]
            sage: [(v, type(v)) for v in P.elements_topological()]
            [((1, 1), <class '__main__.T'>),
             ((2, 1), <class '__main__.T'>),
             ((4, 4), <class '__main__.T'>)]
            sage: [(v, type(v)) for v in P.shells_topological()]
            [((1, 1), <class 'sage.data_structures.mutable_poset.MutablePosetShell'>),
             ((2, 1), <class 'sage.data_structures.mutable_poset.MutablePosetShell'>),
             ((4, 4), <class 'sage.data_structures.mutable_poset.MutablePosetShell'>)]
        """
        for shell in self.shells_topological(**kwargs):
            yield shell.key


    def repr(self, include_special=False, reverse=False):
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
        s += ', '.join(repr(shell) for shell in
                       self.shells_topological(include_special, reverse))
        s += ')'
        return s


    def repr_full(self, reverse=False):
        r"""
        Return a representation with ordering details of the poset.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: print MP().repr_full(reverse=True)
            poset()
            +-- oo
            |   +-- no successors
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   oo
            |   +-- no predecessors
        """
        sortedshells = tuple(
            self.shells_topological(include_special=True, reverse=reverse))
        strings = [self.repr(include_special=False, reverse=reverse)]
        for shell in sortedshells:
            strings.append('+-- ' + repr(shell))
            for rev in (not reverse, reverse):
                what = 'successors' if not rev else 'predecessors'
                if shell.successors(rev):
                    s = '|   +-- ' + what + ':   '
                    s += ', '.join(repr(e) for e in
                                   sorted_set_by_tuple(shell.successors(rev),
                                                       sortedshells))
                else:
                    s = '|   +-- no ' + what
                strings.append(s)
        return '\n'.join(strings)


    __repr__ = repr


    def contains(self, key):
        r"""
        Tests if ``key`` is encapsulated by one of the poset's elements.

        INPUT:

        - ``key`` -- an object.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: T((1, 1)) in P  # indirect doctest
            True
            sage: T((1, 2)) in P  # indirect doctest
            False
        """
        return key in self._shells_


    __contains__ = contains


    def add(self, element):
        r"""
        Add the given object as element to the poset.

        INPUT:

        - ``element`` -- an object (hashable and supporting comparison
          with the operator ``<=``.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: print P.repr_full(reverse=True)
            poset((4, 4), (1, 3), (1, 2), (2, 1), (1, 1))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4)
            +-- (4, 4)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3), (2, 1)
            +-- (1, 3)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2)
            +-- (1, 2)
            |   +-- successors:   (1, 3)
            |   +-- predecessors: (1, 1)
            +-- (2, 1)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 1)
            +-- (1, 1)
            |   +-- successors:   (1, 2), (2, 1)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 1)
            |   +-- no predecessors
            sage: P.add(T((2, 2)))
            sage: reprP = P.repr_full(reverse=True); print reprP
            poset((4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4)
            +-- (4, 4)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3), (2, 2)
            +-- (1, 3)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2)
            +-- (2, 2)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2), (2, 1)
            +-- (1, 2)
            |   +-- successors:   (1, 3), (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (2, 1)
            |   +-- successors:   (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (1, 1)
            |   +-- successors:   (1, 2), (2, 1)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 1)
            |   +-- no predecessors

        When adding an element which is already in the poset, nothing happens::

            sage: e = T((2, 2))
            sage: P.add(e)
            sage: P.repr_full(reverse=True) == reprP
            True

        We can influence the behavior when an element with existing key
        is to be inserted in the poset. For example, we can perform an
        addition on some argument of the elements::

            sage: def add_existing(existing, other):
            ....:     return (existing[0], existing[1] + other[1])
            sage: A = MP(key=lambda k: k[0], element_exists_hook=add_existing)
            sage: A.add((3, 'a'))
            sage: A
            poset((3, 'a'))
            sage: A.add((3, 'b'))
            sage: A
            poset((3, 'ab'))

        We can also deal with cancellations. If the return value of
        our hook-function is ``None``, then the element is removed out of
        the poset::

            sage: def add_existing_None(existing, other):
            ....:     s = existing[1] + other[1]
            ....:     if s == 0:
            ....:         return None
            ....:     return (existing[0], s)
            sage: B = MP(key=lambda k: k[0],
            ....:        element_exists_hook=add_existing_None)
            sage: B.add((7, 42))
            sage: B.add((7, -42))
            sage: B
            poset()

        TESTS::

            sage: R = MP(key=lambda k: T(k[2:3]))
            sage: R.add(T((1, 1, 42)))
            sage: R.add(T((1, 3, 42)))
            sage: R.add(T((2, 1, 7)))
            sage: R.add(T((4, 4, 42)))
            sage: R.add(T((1, 2, 7)))
            sage: R.add(T((2, 2, 7)))
            sage: print R.repr_full(reverse=True)
            poset((1, 1, 42), (2, 1, 7))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (1, 1, 42)
            +-- (1, 1, 42)
            |   +-- successors:   oo
            |   +-- predecessors: (2, 1, 7)
            +-- (2, 1, 7)
            |   +-- successors:   (1, 1, 42)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (2, 1, 7)
            |   +-- no predecessors
        """
        if element is None:
            raise ValueError('None is not allowed as element.')
        key = self.get_key(element)

        if key in self._shells_:
            if self._element_exists_hook_ is not None:
                existing = self.shell(key)
                new = self._element_exists_hook_(existing.element, element)
                if new is None:
                    self.remove(key)
                else:
                    existing._element_ = new
            return

        new = MutablePosetShell(self, element)
        smaller = self.null.covers(new, reverse=False)
        larger = self.oo.covers(new, reverse=True)

        for reverse in (False, True):
            sm = smaller if not reverse else larger
            la = larger if not reverse else smaller
            for shell in sm:
                for e in shell.successors(reverse).intersection(la):
                    e.predecessors(reverse).remove(shell)
                    shell.successors(reverse).remove(e)
                new.predecessors(reverse).add(shell)
                shell.successors(reverse).add(new)

        self._shells_[key] = new


    def remove(self, key, raise_key_error=True):
        r"""
        Remove the given object from the poset.

        INPUT:

        - ``key`` -- the key of an object.

        - ``raise_key_error`` -- (default: ``True``) switch raising
          ``KeyError`` on and off.

        OUTPUT:

        Nothing.

        If the element is not a member and ``raise_key_error`` is set
        (default), raise a ``KeyError``.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: P.add(T((2, 2)))
            sage: print P.repr_full(reverse=True)
            poset((4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4)
            +-- (4, 4)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3), (2, 2)
            +-- (1, 3)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2)
            +-- (2, 2)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2), (2, 1)
            +-- (1, 2)
            |   +-- successors:   (1, 3), (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (2, 1)
            |   +-- successors:   (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (1, 1)
            |   +-- successors:   (1, 2), (2, 1)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 1)
            |   +-- no predecessors
            sage: P.remove(T((1, 2)))
            sage: print P.repr_full(reverse=True)
            poset((4, 4), (1, 3), (2, 2), (2, 1), (1, 1))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4)
            +-- (4, 4)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3), (2, 2)
            +-- (1, 3)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 1)
            +-- (2, 2)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (2, 1)
            +-- (2, 1)
            |   +-- successors:   (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (1, 1)
            |   +-- successors:   (1, 3), (2, 1)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 1)
            |   +-- no predecessors

        TESTS::

            sage: Q = MP(key=lambda k: T(k[0:2]))
            sage: Q.add(T((1, 1, 42)))
            sage: Q.add(T((1, 3, 42)))
            sage: Q.add(T((2, 1, 7)))
            sage: Q.add(T((4, 4, 42)))
            sage: Q.add(T((1, 2, 7)))
            sage: Q.add(T((2, 2, 7)))
            sage: print Q.repr_full(reverse=True)
            poset((4, 4, 42), (1, 3, 42), (2, 2, 7),
                  (1, 2, 7), (2, 1, 7), (1, 1, 42))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4, 42)
            +-- (4, 4, 42)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3, 42), (2, 2, 7)
            +-- (1, 3, 42)
            |   +-- successors:   (4, 4, 42)
            |   +-- predecessors: (1, 2, 7)
            +-- (2, 2, 7)
            |   +-- successors:   (4, 4, 42)
            |   +-- predecessors: (1, 2, 7), (2, 1, 7)
            +-- (1, 2, 7)
            |   +-- successors:   (1, 3, 42), (2, 2, 7)
            |   +-- predecessors: (1, 1, 42)
            +-- (2, 1, 7)
            |   +-- successors:   (2, 2, 7)
            |   +-- predecessors: (1, 1, 42)
            +-- (1, 1, 42)
            |   +-- successors:   (1, 2, 7), (2, 1, 7)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 1, 42)
            |   +-- no predecessors
            sage: Q.remove((1,1))
            sage: print Q.repr_full(reverse=True)
            poset((4, 4, 42), (1, 3, 42), (2, 2, 7), (1, 2, 7), (2, 1, 7))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4, 42)
            +-- (4, 4, 42)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3, 42), (2, 2, 7)
            +-- (1, 3, 42)
            |   +-- successors:   (4, 4, 42)
            |   +-- predecessors: (1, 2, 7)
            +-- (2, 2, 7)
            |   +-- successors:   (4, 4, 42)
            |   +-- predecessors: (1, 2, 7), (2, 1, 7)
            +-- (1, 2, 7)
            |   +-- successors:   (1, 3, 42), (2, 2, 7)
            |   +-- predecessors: null
            +-- (2, 1, 7)
            |   +-- successors:   (2, 2, 7)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 2, 7), (2, 1, 7)
            |   +-- no predecessors
        """
        if key is None:
            raise ValueError('None is not allowed as key.')

        try:
            shell = self._shells_[key]
        except KeyError:
            if not raise_key_error:
                return
            raise KeyError('Key %s is not contained in this poset.' % (key,))

        for reverse in (False, True):
            for p in shell.predecessors(reverse):
                S = p.successors(reverse)
                S.remove(shell)
                D = set(s for s in p.iter_depth_first(reverse)
                        if s in shell.successors(reverse))
                S.update(shell.successors(reverse))
                S.difference_update(D)
        del self._shells_[key]


    def discard(self, key, raise_key_error=False):
        r"""
        Remove the given object from the poset.

        INPUT:

        - ``key`` -- the key of an object.

        - ``raise_key_error`` -- (default: ``False``) switch raising
          ``KeyError`` on and off.

        OUTPUT:

        Nothing.

        If the element is not a member and ``raise_key_error`` is set
        (not default), raise a ``KeyError``.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            sage: P.add(T((1, 3)))
            sage: P.add(T((2, 1)))
            sage: P.add(T((4, 4)))
            sage: P.add(T((1, 2)))
            sage: P.add(T((2, 2)))
            sage: P.discard(T((1, 2)))
            sage: P.remove(T((1, 2)))
            Traceback (most recent call last):
            ...
            KeyError: 'Key (1, 2) is not contained in this poset.'
            sage: P.discard(T((1, 2)))
        """
        return self.remove(key, raise_key_error)


    def pop(self, **kwargs):
        r"""
        Remove and return an arbitrary poset element.

        INPUT:

        - ``kwargs`` -- arguments are passed to :meth:`shells_topological`.

        OUTPUT:

        An object.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP(key=lambda c: -c)
            sage: P.add(3)
            sage: P
            poset(3)
            sage: P.pop()
            3
            sage: P
            poset()
            sage: P.pop()
            Traceback (most recent call last):
            ...
            KeyError: 'pop from an empty poset'
        """
        try:
            del kwargs['include_special']
        except KeyError:
            pass
        kwargs['include_special'] = False

        try:
            shell = next(self.shells_topological(**kwargs))
        except StopIteration:
            raise KeyError('pop from an empty poset')
        self.remove(shell.key)
        return shell.element


    def union(left, *right):
        r"""
        Return the union of the given posets as a new poset

        INPUT:

        - ``left`` -- a poset.

        - ``right`` -- a poset or an iterable. In the latter case the
          iterated objects are seen as elements of a poset.

        OUTPUT:

        A poset.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.union(Q)
            poset(3, 4, 7, 8, 42)

        TESTS::

            sage: P.union(P, Q, Q, P)
            poset(3, 4, 7, 8, 42)
       """
        new = left.copy()
        for r in right:
            new.update(r)
        return new


    def union_update(self, other):
        r"""
        Update the poset with the union of itself and another poset.

        INPUT:

        - ``other`` -- a poset or an iterable. In the latter case the
          iterated objects are seen as elements of a poset.

        OUTPUT:

        Nothing.

        .. TODO::

            Use the already existing information in the other poset to speed
            up this function. (At the moment each element of the other poset
            is inserted one by one and without using this information.)

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.union_update(Q)
            sage: P
            poset(3, 4, 7, 8, 42)

        TESTS::

            sage: Q.update(P)
            sage: Q
            poset(3, 4, 7, 8, 42)
        """
        try:
            it = other.elements()
        except AttributeError:
            it = iter(other)
        for element in it:
            self.add(element)


    update = union_update  # as in a Python set
    r"""
    Alias of :meth:`union_update`.
    """


    def difference(left, *right):
        r"""
        Return a new poset where all elements of this poset, which are
        contained in one of the other given posets, are removed.

        INPUT:

        - ``left`` -- a poset.

        - ``right`` -- a poset or an iterable. In the latter case the
          iterated objects are seen as elements of a poset.

        OUTPUT:

        A poset.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.difference(Q)
            poset(3, 7)

        TESTS::

            sage: P.difference(Q, Q)
            poset(3, 7)
            sage: P.difference(P)
            poset()
            sage: P.difference(Q, P)
            poset()
        """
        new = left.copy()
        for r in right:
            new.difference_update(r)
        return new


    def difference_update(self, other):
        r"""
        Remove all elements of another poset from this poset.

        INPUT:

        - ``other`` -- a poset or an iterable. In the latter case the
          iterated objects are seen as elements of a poset.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.difference_update(Q)
            sage: P
            poset(3, 7)
        """
        try:
            it = other.keys()
        except AttributeError:
            it = iter(other)
        for key in it:
            self.discard(key)


    def intersection(left, *right):
        r"""
        Return the intersection of the given posets as a new poset

        INPUT:

        - ``left`` -- a poset.

        - ``right`` -- a poset or an iterable. In the latter case the
          iterated objects are seen as elements of a poset.

        OUTPUT:

        A poset.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.intersection(Q)
            poset(42)

        TESTS::

            sage: P.intersection(P, Q, Q, P)
            poset(42)
        """
        new = left.copy()
        for r in right:
            new.intersection_update(r)
        return new


    def intersection_update(self, other):
        r"""
        Update the poset with the intersection of itself and another poset.

        INPUT:

        - ``other`` -- a poset or an iterable. In the latter case the
          iterated objects are seen as elements of a poset.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.intersection_update(Q)
            sage: P
            poset(42)
        """
        try:
            keys = tuple(self.keys())
        except AttributeError:
            keys = tuple(iter(self))
        for key in keys:
            if key not in other:
                self.discard(key)


    def symmetric_difference(left, right):
        r"""
        Return the symmetric difference of two posets as a new poset.

        INPUT:

        - ``left`` -- a poset.

        - ``right`` -- a poset.

        OUTPUT:

        A poset.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.symmetric_difference(Q)
            poset(3, 4, 7, 8)
        """
        new = left.copy()
        new.symmetric_difference_update(right)
        return new


    def symmetric_difference_update(self, other):
        r"""
        Update the poset with the symmetric difference of itself and
        another poset.

        INPUT:

        - ``other`` -- a poset.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.symmetric_difference_update(Q)
            sage: P
            poset(3, 4, 7, 8)
        """
        T = other.difference(self)
        self.difference_update(other)
        self.union_update(T)


    def is_disjoint(self, other):
        r"""
        Return if another poset is disjoint to this poset.

        INPUT:

        - ``other`` -- a poset or an iterable. In the latter case the
          iterated objects are seen as elements of a poset.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.is_disjoint(Q)
            False
            sage: P.is_disjoint(Q.difference(P))
            True
        """
        return all(key not in other for key in self.keys())


    isdisjoint = is_disjoint  # as in a Python set
    r"""
    Alias of :meth:`is_disjoint`.
    """


    def is_subset(self, other):
        r"""
        Return if another poset contains this poset, i.e., if this poset
        is a subset of the other poset.

        INPUT:

        - ``other`` -- a poset or an iterable. In the latter case the
          iterated objects are seen as elements of a poset.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.is_subset(Q)
            False
            sage: Q.is_subset(P)
            False
            sage: P.is_subset(P)
            True
            sage: P.is_subset(P.union(Q))
            True
        """
        return all(key in other for key in self.keys())


    issubset = is_subset  # as in a Python set
    r"""
    Alias of :meth:`is_subset`.
    """


    def is_superset(self, other):
        r"""
        Return if this poset contains another poset, i.e., if this poset
        is a superset of the other poset.

        INPUT:

        - ``other`` -- a poset or an iterable. In the latter case the
          iterated objects are seen as elements of a poset.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(3); P.add(42); P.add(7); P
            poset(3, 7, 42)
            sage: Q = MP()
            sage: Q.add(4); Q.add(8); Q.add(42); Q
            poset(4, 8, 42)
            sage: P.is_superset(Q)
            False
            sage: Q.is_superset(P)
            False
            sage: P.is_superset(P)
            True
            sage: P.union(Q).is_superset(P)
            True
        """
        try:
            it = other.keys()
        except AttributeError:
            it = iter(other)
        return all(key in self for key in it)


    issuperset = is_superset  # as in a Python set
    r"""
    Alias of :meth:`is_superset`.
    """


# *****************************************************************************


class MutableTosetShell(MutablePosetShell):
    r"""
    A shell containing an element of a mutable toset (totally ordered set).

    .. TODO::

        Implement this class.
    """
    pass


# *****************************************************************************


class MutableToset(MutablePoset):
    r"""
    A mutable toset (totally ordered set) as data structure.

    .. TODO::

        Implement this class.
    """
    pass


# *****************************************************************************
