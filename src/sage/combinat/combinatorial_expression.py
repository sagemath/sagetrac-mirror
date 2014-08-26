# -*- coding: utf-8 -*-
r"""

Examples
========

::

    sage: R = CombinatorialExpressionRing(); R
    Combinatorial Expression Ring


Non-adjacent Forms
------------------

::

    sage: # NAF = (0 + P0 + M0)* (P + M + e)  with (P = 1, M = -1)


Plane/Ordered Binary Trees
--------------------------

::

    sage: T = R(SR('T'), flavor='unlabeled', function=True); T
    T = None
    sage: z = R(SR('z')); z
    z
    sage: e = R(SR(1), size=0); e
    1
    sage: T.assign(e + z * T * T); T
    T = 1 + z*T*T

.. NOTE::

    be aware: there is sage.combinat.binary_tree.BinaryTrees




AUTHORS:

- Daniel Krenn (2014-04-01): initial version
- Daniel Krenn (2014-08-24): rewritten to parent/element structure

"""

#*****************************************************************************
#       Copyright (C) 2014 Daniel Krenn <devel@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import sage
from sage.rings.integer import Integer
from sage.misc.lazy_attribute import lazy_attribute
from itertools import izip

#*****************************************************************************
# Flavor
#*****************************************************************************

class _GenericFlavor_(sage.structure.sage_object.SageObject):

    def is_unlabeled(self):
        """
        Returns whether combinatorial expression is unlabeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``True`` if unlabeled, ``False`` if labeled, ``None`` otherwise.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _GenericFlavor_)
            sage: F = _GenericFlavor_()
            sage: F.is_unlabeled() is None
            True
        """
        return None


    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``True`` if labeled, ``False`` if unlabeled, ``None`` otherwise.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _GenericFlavor_)
            sage: F = _GenericFlavor_()
            sage: F.is_labeled() is None
            True
        """
        return None


    # ------------------------------------------------------------------------


    @staticmethod
    def _class_with_prefix_(prefix, classname):
        """
        Returns a combinatorial expression class specified by input
        parameters.

        INPUT:

        - ``prefix`` -- a string.

        - ``classname`` -- a string.

        OUTPUT:

        A class.

        EXAMPLES::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _GenericFlavor_)
            sage: _GenericFlavor_._class_with_prefix_('Generic', 'Atom')
            <class 'sage.combinat.combinatorial_expression.GenericAtom'>
        """
        return globals()[prefix + classname]


    # ------------------------------------------------------------------------


    @classmethod
    def class_generic(cls, classname):
        """
        Returns a combinatorial expression class with generic flavor.

        INPUT:

        - ``classname`` -- a string.

        OUTPUT:

        A class

        EXAMPLES::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _GenericFlavor_)
            sage: _GenericFlavor_.class_generic('Atom')
            <class 'sage.combinat.combinatorial_expression.GenericAtom'>
         """
        return cls._class_with_prefix_('Generic', classname)


    @classmethod
    def class_unlabeled(cls, classname):
        """
        Returns a combinatorial expression class with generic flavor.

        INPUT:

        - ``classname`` -- a string.

        OUTPUT:

        A class

        EXAMPLES::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _GenericFlavor_)
            sage: _GenericFlavor_.class_unlabeled('Atom')
            <class 'sage.combinat.combinatorial_expression.UnlabeledAtom'>
         """
        return cls._class_with_prefix_('Unlabeled', classname)


    @classmethod
    def class_labeled(cls, classname):
        """
        Returns a combinatorial expression class with generic flavor.

        INPUT:

        - ``classname`` -- a string.

        OUTPUT:

        A class

        EXAMPLES::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _GenericFlavor_)
            sage: _GenericFlavor_.class_labeled('Atom')
            <class 'sage.combinat.combinatorial_expression.LabeledAtom'>
            """
        return cls._class_with_prefix_('Labeled', classname)


# ----------------------------------------------------------------------------


class _UnlabeledFlavor_(_GenericFlavor_):

    def is_unlabeled(self):
        """
        Returns whether combinatorial expression is unlabeled or not.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` if unlabeled, ``False`` if labeled, ``None`` otherwise.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _UnlabeledFlavor_)
            sage: F = _UnlabeledFlavor_()
            sage: F.is_unlabeled()
            True
        """
        return True


    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` if labeled, ``False`` if unlabeled, ``None`` otherwise.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _UnlabeledFlavor_)
            sage: F = _UnlabeledFlavor_()
            sage: F.is_labeled()
            False
        """
        return False


# ----------------------------------------------------------------------------


class _LabeledFlavor_(_GenericFlavor_):

    def is_unlabeled(self):
        """
        Returns whether combinatorial expression is unlabeled or not.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` if unlabeled, ``False`` if labeled, ``None`` otherwise.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _LabeledFlavor_)
            sage: F = _LabeledFlavor_()
            sage: F.is_unlabeled()
            False
        """
        return False


    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` if labeled, ``False`` if unlabeled, ``None`` otherwise.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _LabeledFlavor_)
            sage: F = _LabeledFlavor_()
            sage: F.is_labeled()
            True
        """
        return True


# ----------------------------------------------------------------------------


class _EmptyFlavor_(_GenericFlavor_):

    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` since this instance is labeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _EmptyFlavor_)
            sage: F = _EmptyFlavor_()
            sage: F.is_labeled()
            True
        """
        return True


    def is_unlabeled(self):
        """
        Returns whether combinatorial expression is unlabeled or not.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` since this instance is unlabeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _EmptyFlavor_)
            sage: F = _EmptyFlavor_()
            sage: F.is_unlabeled()
            True
        """
        return True


#*****************************************************************************
# Element: Base
#*****************************************************************************


class GenericExpression(
    sage.structure.element.RingElement,
    _GenericFlavor_):
    """
    Abstract base class for all combinatorial expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A new (abstract) combinatorial expression.

    .. NOTE::

        This should not be used directly. Use a derived class.

    .. SEEALSO::

        :class:`UnlabeledExpression`,
        :class:`LabeledExpression'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     GenericFunction,
        ....:     GenericDisjointUnion,
        ....:     GenericCartesianProduct,
        ....:     GenericAtom)
        sage: R = CombinatorialExpressionRing()
        sage: z = GenericAtom(R, SR('z'))
        sage: z
        z
        sage: eps = GenericAtom(R, SR(1))
        sage: eps
        1
        sage: T = GenericFunction(R, SR('T')); T
        T = None
        sage: A = GenericCartesianProduct(R, z, T, T); A
        z*T*T
        sage: B = GenericDisjointUnion(R, eps, A); B
        1 + z*T*T
        sage: T.assign(B)
        sage: T
        T = 1 + z*T*T

    .. automethod:: _add_

    .. automethod:: _mul_
    """
    def __init__(self, parent, *operands):
        """
        See :class:`GenericExpression` for more information.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: z = GenericExpression(R); z  # indirect doctest
            GenericExpression()
            sage: z.parent()
            Combinatorial Expression Ring
        """
        self.assign(*operands)
        super(GenericExpression, self).__init__(parent)


    # ------------------------------------------------------------------------


    def assign(self, *operands):
        """
        Assign operands.

        INPUT:

        - ``*operands`` -- the operands of the operation defined by this class.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: T = R(SR('T'), function=True); T
            T = None
            sage: T.assign(R(SR(1), size=0))
            sage: T
            T = 1
        """
        self._operands_ = tuple(operands)


    # ------------------------------------------------------------------------


    def operands(self):
        """
        Return the operands.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple of operands.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: c = R(SR('a')) + R(SR('b')); c
            a + b
            sage: c.operands()
            (a, b)
        """
        return self._operands_


    def iter_operands(self):
        """
        Return an iterator of the operands.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: c = R(SR('a')) + R(SR('b')); c
            a + b
            sage: for o in c.iter_operands():
            ....:     print o, type(o)
            a <class 'sage.combinat.combinatorial_expression.GenericAtom_with_category'>
            b <class 'sage.combinat.combinatorial_expression.GenericAtom_with_category'>
        """
        return iter(self._operands_)


    #------------------------------------------------------------------------


    def iter_all_expressions(self, memo=None):
        """
        Return an iterator of all expressions contained in ``self``

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        .. NOTE::

            Each element is returned exactly once.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: T = R(SR('T'), function=True)
            sage: z = R(SR('z'))
            sage: e = R(SR(1), size=0)
            sage: T.assign(e + z * T * T); T
            T = 1 + z*T*T
            sage: for expr in T.iter_all_expressions():
            ....:     print expr
            T = 1 + z*T*T
            1 + z*T*T
            1
            z*T*T
            z*T
            z

        TESTS::

            sage: [expr for expr in (z*z*z).iter_all_expressions()]
            [z*z*z, z*z, z]
            sage: [expr for expr in
            ....:  R.Operators.cartesian_product(z, z, z).iter_all_expressions()]
            [z*z*z, z]
        """
        if memo is None:
            memo = {}

        if not self._update_memo_(memo):
            return

        yield self

        for o in self.iter_operands():
            for e in o.iter_all_expressions(memo):
                yield e


    #------------------------------------------------------------------------


    def _name_(self):
        """
        Return the name (the name of the class) of ``self``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        This function strips of a trailing ``_with_category``.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: R(SR('z'))._name_()
            'GenericAtom'
            sage: y = GenericExpression(R); y  # indirect doctest
            GenericExpression()
        """
        s = self.__class__.__name__
        i = s.rfind('_with_category')
        if i >= 0:
            return s[:i]
        else:
            return s


    #------------------------------------------------------------------------


    def _update_memo_(self, memo):
        """
        Update the memo during recursive evaluation.

        INPUT:

        - ``memo`` -- a dictionary.

        OUTPUT:

        ``True`` if ``self`` was not already in memo, otherwise ``False``.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: y = GenericExpression(R)
            sage: memo = {}; memo
            {}
            sage: y._update_memo_(memo)
            True
            sage: memo
            {...: True}
            sage: memo.popitem() == (id(y), True)
            True
        """
        key = id(self)
        m = memo.get(key)
        memo[key] = True
        return m is None


    #------------------------------------------------------------------------


    def _repr_(self):
        """
        Return a representation string.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: GenericExpression(R)._repr_()
            'GenericExpression()'
        """
        memo = {}
        return self._repr_main_(memo)


    def _repr_main_(self, memo):
        """
        Returns a representation string.

        INPUT:

        - ``memo`` -- a dictionary.

        OUTPUT:

        A string.

        This function is usually overridden in a derived class.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: GenericExpression(R)._repr_main_({})
            'GenericExpression()'
        """
        return self._repr_recursive_(memo)


    def _repr_recursive_(self, memo):
        """
        Recursively build the representation string.

        INPUT:

        - ``memo`` -- a dictionary.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: GenericExpression(R)._repr_recursive_({})
            'GenericExpression()'
            sage: z = GenericExpression(R)
            sage: z.assign(z)
            sage: z
            GenericExpression(REC(GenericExpression))
        """
        if not self._update_memo_(memo):
            return self._repr_recursion_()
        return self._repr_join_(
            self._repr_make_parenthesis_(
                o._repr_main_(memo), self._repr_need_parenthesis_(o))
            for o in self.iter_operands())


    @staticmethod
    def _repr_make_parenthesis_(s, make=True):
        """
        Make parenthesis around a string.

        INPUT:

        - ``s`` -- a string.

        - ``make`` -- a boolean indicating whether to make parenthesis
          or not.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: GenericExpression._repr_make_parenthesis_(
            ....:     '42', False)
            '42'
            sage: GenericExpression._repr_make_parenthesis_(
            ....:     '42', True)
            '(42)'
            """
        if make:
            return '(' + s + ')'
        else:
            return s


    def _repr_need_parenthesis_(self, operand):
        """
        Return if parenthesis are needed around the representation string
        of the given operand.

        INPUT:

        - ``operand`` -- the (inner) operand.

        OUTPUT:

        ``False``.

        This function is usually overridden in a derived class.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: GenericExpression(R)._repr_need_parenthesis_(
            ....:     GenericExpression(R))
            False
        """
        return False


    def _repr_join_(self, reprs_of_operands):
        """
        Join the given representation strings together.

        INPUT:

        - ``reprs_of_operands`` -- an iterable of strings.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: GenericExpression(R)._repr_join_(
            ....:     str(i) for i in srange(4))
            'GenericExpression(0, 1, 2, 3)'
        """
        return (self._name_() +
                self._repr_make_parenthesis_(
                ', '.join(o for o in reprs_of_operands)))


    def _repr_recursion_(self):
        """
        Return a string signaling that a recursion was encountered.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: z = GenericExpression(R)
            sage: z.assign(z)
            sage: z  # indirect doctest
            GenericExpression(REC(GenericExpression))
        """
        return 'REC' + self._repr_make_parenthesis_(self._name_())


    #------------------------------------------------------------------------


    def __cmp__(left, right, memo=None):
        """
        Return whether operators are equal or not.

        INPUT:

        - ``left`` -- left expression.

        - ``right`` -- right expression.

        OUTPUT:

        ``0`` if operands are equal, otherwise ``-1``.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: B = R(SR('B'), function=True)
            sage: T = R(SR('T'), function=True)
            sage: B == T
            True
            sage: y = R(SR('y'))
            sage: z = R(SR('z'))
            sage: y == z
            True
            sage: e = R(SR(1), size=0)
            sage: f = R(SR(1), size=0)
            sage: e == f
            True
            sage: B.assign(f + y * B * B); B
            B = 1 + y*B*B
            sage: B == T
            False
            sage: T.assign(e + z * T * T); T
            T = 1 + z*T*T
            sage: B == T
            True
            sage: T == B
            True

        TESTS::

            sage: a = R(SR('a'), size=3)
            sage: b = R(SR('b'), size=2)
            sage: a == b
            False
            sage: j = R(SR('j'))
            sage: k = R(SR('k'))
            sage: j * k == j * j
            False
            sage: j * j == j * k
            False
            sage: j * k == k * j
            True
        """
        memo = {}
        # _cmp_main_ is True if equal, otherwise False; convert it to
        # be compatible with cmp(...)
        return 0 if left._cmp_main_(right, memo) else -1


    def _cmp_main_(left, right, memo):
        """
        Return whether operators are equal or not.

        INPUT:

        - ``left`` -- left expression.

        - ``right`` -- right expression.

        - ``memo`` -- a dictionary.

        OUTPUT:

        ``True`` if operands are equal, otherwise ``False``.

        The keys of the dictionary ``memo`` are identifiers of the
        left combinatorial expressions, their values are pairs of
        identifiers of the right combinatorial expression and ``True``
        if they compare equal, otherwise ``False``.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: a = R(SR('a'), size=3)
            sage: b = R(SR('b'), size=2)
            sage: a == b  # indirect doctest
            False
            sage: j = R(SR('j'))
            sage: k = R(SR('k'))
            sage: j * k == j * j  # indirect doctest
            False
            sage: j * j == j * k  # indirect doctest
            False
            sage: j * k == k * j  # indirect doctest
            True
        """
        # keys for memo
        idleft = id(left)
        idright = id(right)

        # look if we already worked with one of the operands
        lookup = memo.get(idleft)
        if lookup is not None:
            idr, result = lookup
            if idr == idright:
                return result
            return False
        elif memo.get(idright) is not None:
            # we worked with right, but not with left, so switch operands
            return right._cmp_main_(left, memo)

        # add to memo to prevent cycling
        memo[idleft] = (idright, True)
        memo[idright] = (idleft, True)

        # check if equal
        if not left._cmp_without_operands_(right):
            result = False
        else:
            # check operands
            L = left.operands()
            R = right.operands()
            if len(L) != len(R):
                result = False
            else:
                result = all(l._cmp_main_(r, memo) for l, r in izip(L, R))

        # store result in memo
        memo[idleft] = (idright, result)
        memo[idright] = (idleft, result)

        return result


    def _cmp_without_operands_(left, right):
        """
        Return whether operators are equal or not (but without
        comparing the operands).

        INPUT:

        - ``left`` -- left expression.

        - ``right`` -- right expression.

        OUTPUT:

        ``True`` if operators are equal, otherwise ``False``.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: S = R(SR('S'), function=True)
            sage: T = R(SR('T'), function=True)
            sage: S._cmp_without_operands_(T)
            True
        """
        return left.__class__ == right.__class__


    #------------------------------------------------------------------------


    def is_singleton(self):
        """
        Return whether ``self`` is a singleton.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: B = R(SR('B'), function=True)
            sage: B.is_singleton()
            False
            sage: y = R(SR('y'))
            sage: y.is_singleton()
            True
            sage: e = R(SR(1), size=0)
            sage: e.is_singleton()
            True
        """
        return False


    def is_empty_singleton(self):
        """
        Return whether ``self`` is the empty singleton.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: B = R(SR('B'), function=True)
            sage: B.is_empty_singleton()
            False
            sage: y = R(SR('y'))
            sage: y.is_empty_singleton()
            False
            sage: e = R(SR(1), size=0)
            sage: e.is_empty_singleton()
            True
        """
        return False

    #------------------------------------------------------------------------


    def _apply_operator_(self, operatorclassname, *operands):
        """
        Return an instance of the specified class (as operator) with
        given operands.

        INPUT:

        - ``operatorclassname`` -- the name of the class (without flavor).

        - ``*operands`` -- the operands of the operation.

        OUTPUT:

        A new combinatorial expression.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: g = R(SR('g'))
            sage: u = R(SR('u'), unlabeled=True)
            sage: l = R(SR('l'), labeled=True)
            sage: a = g + g  # indirect doctest
            sage: a.is_unlabeled(), a.is_labeled()
            (None, None)
            sage: b = u + u  # indirect doctest
            sage: b.is_unlabeled(), b.is_labeled()
            (True, False)
            sage: c = l + l  # indirect doctest
            sage: c.is_unlabeled(), c.is_labeled()
            (False, True)
            sage: d = g + u  # indirect doctest
            sage: d.is_unlabeled(), d.is_labeled()
            (None, None)
            sage: e = u + l  # indirect doctest
            sage: e.is_unlabeled(), e.is_labeled()
            (None, None)
            sage: f = l + g  # indirect doctest
            sage: f.is_unlabeled(), f.is_labeled()
            (None, None)
        """
        if all(o.is_unlabeled() for o in operands):
            flavor = 'Unlabeled'
        elif all(o.is_labeled() for o in operands):
            flavor = 'Labeled'
        else:
            flavor = 'Generic'

        return globals()[flavor + operatorclassname](
            operands[0].parent(), *operands)


    #------------------------------------------------------------------------


    def disjoint_union(self, *others):
        """
        Returns the disjoint union of ``self`` and the given operands.

        INPUT:

        - ``*others`` -- operands.

        OUTPUT:

        A new combinatorial expression.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: y = R(SR('y'))
            sage: z = R(SR('z'))
            sage: a = y.disjoint_union(z); a
            y + z
            sage: b = R.Operators.disjoint_union(y, z); b
            y + z
            sage: c = y + z; c
            y + z
            sage: a == b == c
            True

        .. SEEALSO::

            :meth:`_add_`,
            :meth:`Operators.disjoint_union`.

            :meth:`cartesian_product`.
        """
        return Operators.disjoint_union(self, *others)


    def _disjoint_union_(self, *others):
        """
        Returns the disjoint union of ``self`` and the given operands.

        INPUT:

        - ``*others`` -- operands.

        OUTPUT:

        A new combinatorial expression.

        This helper function does the actual operation. It can be
        assumed that all operands have the same parent.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: y = R(SR('y'))
            sage: z = R(SR('z'))
            sage: y._disjoint_union_(z)
            y + z

        .. SEEALSO::

            :meth:`_add_`,
            :meth:`disjoint_union`,
            :meth:`Operators.disjoint_union`.
        """
        return self._apply_operator_('DisjointUnion', self, *others)


    def _add_(left, right):
        """
        Returns the disjoint union of the given operands.

        INPUT:

        - ``left`` -- left operand.

        - ``right`` -- right operand.

        OUTPUT:

        A new combinatorial expression.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: y = R(SR('y'))
            sage: z = R(SR('z'))
            sage: y + z
            y + z

        .. SEEALSO::

            :meth:`disjoint_union`,
            :meth:`Operators.disjoint_union`.

            :meth:`cartesian_product`.
        """
        return left.disjoint_union(right)


    #------------------------------------------------------------------------


    def cartesian_product(self, *others):
        """
        Returns the cartesian product of ``self`` and the given operands.

        INPUT:

        - ``*others`` -- operands.

        OUTPUT:

        A new combinatorial expression.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: y = R(SR('y'))
            sage: z = R(SR('z'))
            sage: a = y.cartesian_product(z); a
            y*z
            sage: b = R.Operators.cartesian_product(y, z); b
            y*z
            sage: c = y * z; c
            y*z
            sage: a == b == c
            True

        .. SEEALSO::

            :meth:`_mul_`,
            :meth:`Operators.cartesian_product`.

            :meth:`disjoint_union`.
        """
        return Operators.cartesian_product(self, *others)


    def _cartesian_product_(self, *others):
        """
        Returns the cartesian product of ``self`` and the given operands.

        INPUT:

        - ``*others`` -- operands.

        OUTPUT:

        A new combinatorial expression.

        This helper function does the actual operation. It can be
        assumed that all operands have the same parent.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: y = R(SR('y'))
            sage: z = R(SR('z'))
            sage: y._cartesian_product_(z)
            y*z

        .. SEEALSO::

            :meth:`_mul_`,
            :meth:`cartesian_product`,
            :meth:`Operators.cartesian_product`.
        """
        return self._apply_operator_('CartesianProduct', self, *others)


    def _mul_(left, right):
        """
        Returns the cartesian product of the given operands.

        INPUT:

        - ``left`` -- left operand.

        - ``right`` -- right operand.

        OUTPUT:

        A new combinatorial expression.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: y = R(SR('y'))
            sage: z = R(SR('z'))
            sage: y * z
            y*z

        .. SEEALSO::

            :meth:`cartesian_product`,
            :meth:`Operators.cartesian_product`.

            :meth:`disjoint_union`.
        """
        return left.cartesian_product(right)


    #------------------------------------------------------------------------


    def _preprocess_size_(self, size):
        """
        Determine the size and return it as dictionary.

        INPUT:

        - ``size`` -- either a dictionary (keys are singletons, values
          are their sizes) or a (single) size. If the latter, then it
          is assumed that ``self`` has a unique singleton of size at
          least `1`.

        OUTPUT:

        The size as dictionary.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: T = R(SR('T'), function=True)
            sage: z = R(SR('z'))
            sage: e = R(SR(1), size=0)
            sage: T.assign(e + z * T * T); T
            T = 1 + z*T*T
            sage: T._preprocess_size_(42)
            {z: 42}
            sage: T._preprocess_size_({z: 73})
            {z: 73}
            sage: T._preprocess_size_({T: 666})
            Traceback (most recent call last):
            ...
            TypeError: Parameter size is not well formed.
        """
        if isinstance(size, dict):
            if not all(k in self.parent() and
                       hasattr(k, 'size')
                       for k, _ in size.iteritems()):
                raise TypeError('Parameter size is not well formed.')
            return size

        # at this point we interpret size as a size

        singletons = tuple(s for s in self.iter_all_expressions()
                           if s.is_singleton() and not s.is_empty_singleton())
        if len(singletons) != 1:
            raise ValueError('Do not know to which singleton the '
                             'size %s should be associated (found %s '
                             'singletons).' % (size, len(singletons)))
        return {singletons[0]: size}


    #------------------------------------------------------------------------


    def iter_elements(self, size):
        """
        Return an iterator over all elements described by this
        combinatorial expression.

        INPUT:

        - ``size`` -- the desired size of each of the elements.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: TODO  # not tested -- add a lot of examples here
        """
        return self._iter_elements_(self._preprocess_size_(size))


    def _iter_elements_(self, size):
        """
        Return the actual iterator.

        INPUT:

        - ``size`` -- a dictionary.

        OUTPUT:

        An iterator.

        This method should be overridden in derived class.

        EXAMPLES::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: GenericExpression(R)._iter_elements_({})
            Traceback (most recent call last):
            ...
            NotImplementedError: Iteration for an instance of
            GenericExpression is not implemented.
        """
        raise NotImplementedError('Iteration for an instance of %s is '
                                  'not implemented.' % (self._name_(),))


    #------------------------------------------------------------------------


    def random_element(self, size):
        """
        Return a random element described by this
        combinatorial expression.

        INPUT:

        - ``size`` -- the desired size of the element.

        OUTPUT:

        An element.

        EXAMPLES::

            sage: TODO  # not tested -- add a lot of examples here
        """
        return self._random_element_(self._preprocess_size_(size))


    def _random_element_(self, size):
        """
        Return the actual random element.

        INPUT:

        - ``size`` -- a dictionary.

        OUTPUT:

        An element.

        This method should be overridden in derived class.

        EXAMPLES::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericExpression)
            sage: R = CombinatorialExpressionRing()
            sage: GenericExpression(R)._random_element_({})
            Traceback (most recent call last):
            ...
            NotImplementedError: Generation of an random element for
            an instance of GenericExpression not implemented.
        """
        raise NotImplementedError('Generation of an random element '
                                  'for an instance of %s '
                                  'not implemented.' % (self._name_(),))


# ----------------------------------------------------------------------------


class UnlabeledExpression(
    GenericExpression,
    _UnlabeledFlavor_):
    """
    Abstract base class for all unlabeled combinatorial expressions.

    For details see :class:`GenericExpression`.

    .. SEEALSO::

        :class:`GenericExpression`,
        :class:`LabeledExpression'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.
    """
    pass


# ----------------------------------------------------------------------------


class LabeledExpression(
    GenericExpression,
    _LabeledFlavor_):
    """
    Abstract base class for all labeled combinatorial expressions.

    For details see :class:`GenericExpression`.

    .. SEEALSO::

        :class:`GenericExpression`,
        :class:`UnlabeledExpression`.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.
    """
    pass


#*****************************************************************************
# Element: Function
#*****************************************************************************


class GenericFunction(GenericExpression):
    """
    A class representing a combinatorial function.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A new combinatorial function.

    .. SEEALSO::

        :class:`UnlabeledFunction`,
        :class:`LabeledFunction'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     GenericFunction)
        sage: R = CombinatorialExpressionRing()
        sage: GenericFunction(R, SR('T'))
        T = None
    """
    def __init__(self, parent, expression, *operands):
        """
        See :class:`GenericFunction` for more information.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericFunction)
            sage: R = CombinatorialExpressionRing()
            sage: z = GenericFunction(R, SR('z')); z  # indirect doctest
            z = None
        """
        super(GenericFunction, self).__init__(parent, *operands)

        # we allow any expressions, so the following is commented out
        #if expression not in self.parent().base_ring():
        #    raise ValueError('%s is not in %s' % (expression,
        #                                          self.parent().base_ring()))
        self._expression_ = expression


    #------------------------------------------------------------------------


    def _repr_(self):
        """
        Return a representation string.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: T = R(SR('T'), function=True)
            sage: T._repr_()
            'T = None'
            sage: T.assign(R(SR('z')))
            sage: T._repr_()
            'T = z'
        """
        memo = {}
        s = repr(self._expression_) + ' = '
        if self._operands_:
            return s + self._repr_recursive_(memo)
        else:
            return s + 'None'


    def _repr_main_(self, memo):
        """
        Returns a representation string.

        INPUT:

        - ``memo`` -- a dictionary.

        OUTPUT:

        A string.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: T = R(SR('T'), function=True)
            sage: T.assign(R(SR('z')) * T)
            sage: T._repr_()
            'T = z*T'
        """
        self._update_memo_(memo)
        return repr(self._expression_)


    def _repr_join_(self, reprs_of_operands):
        """
        Join the given representation strings together.

        INPUT:

        - ``reprs_of_operands`` -- an iterable of strings.

        OUTPUT:

        A string.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: T = R(SR('T'), function=True)
            sage: T.assign(R(SR('z')), T)
            sage: T._repr_()  # result looks weired...
            'T = z, T'
        """
        return ', '.join(o for o in reprs_of_operands)


# ----------------------------------------------------------------------------


class UnlabeledFunction(
    GenericFunction,
    UnlabeledExpression):
    """
    A class representing a unlabeled combinatorial function.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A new unlabeled combinatorial expression.

    .. SEEALSO::

        :class:`GenericFunction`,
        :class:`LabeledFunction',

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     UnlabeledFunction)
        sage: R = CombinatorialExpressionRing()
        sage: T = UnlabeledFunction(R, SR('T')); T
        T = None
        sage: T.is_unlabeled()
        True
        sage: T.is_labeled()
        False
    """
    pass


# ----------------------------------------------------------------------------

class LabeledFunction(
    GenericFunction,
    LabeledExpression):
    """
    A class representing a labeled combinatorial function.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A new labeled combinatorial expression.

    .. SEEALSO::

        :class:`GenericFunction`,
        :class:`UnlabeledFunction`.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     LabeledFunction)
        sage: R = CombinatorialExpressionRing()
        sage: LabeledFunction(R, SR('T'))
        T = None
    """
    pass


#*****************************************************************************
# Element: Singleton
#*****************************************************************************

class GenericSingleton(GenericExpression):
    """
    A class representing a singleton of given size.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``singleton`` -- the singleton.

    - ``size`` -- the size of the singleton.

    OUTPUT:

    A new singleton.

    .. SEEALSO::

        :class:`UnlabeledSingleton`,
        :class:`LabeledSingleton'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`.

    EXAMPLES::

        sage: R = CombinatorialExpressionRing()
        sage: y = R(SR('y'), size=2); y
        y
        sage: y.size()
        2
        sage: z = R(SR('z')); z
        z
        sage: z.size()
        1
        sage: e = R(SR(1), size=0); e
        1
        sage: e.size()
        0
    """

    def __init__(self, parent, singleton, size):
        """
        See :class:`GenericSingleton` for more information.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericSingleton)
            sage: R = CombinatorialExpressionRing()
            sage: z = GenericSingleton(R, SR('z'), size=2); z  # indirect doctest
            z
        """
        super(GenericSingleton, self).__init__(parent)
        self._set_singleton_(singleton, size)


    #------------------------------------------------------------------------


    def _set_singleton_(self, singleton, size):
        """
        Set data.

        INPUT:

        - ``singleton`` -- an element of the base ring.

        - ``size`` -- the size of the singleton.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: z = R(SR('z')); z  # indirect doctest
            z
            sage: z.size()
            1
            sage: A = R(Automaton()); A
            Automaton with 0 states
            sage: type(A)
            <class 'sage.combinat.combinatorial_expression.GenericAtom_with_category'>
        """
        # we allow any singletons, so the following is commented out
        #if singleton not in self.parent().base_ring():
        #    raise ValueError('%s is not in %s' % (singleton,
        #                                          self.parent().base_ring()))
        self._singleton_ = singleton
        self._size_ = size


    def singleton(self):
        """
        Return the singleton (element of the base ring).

        INPUT:

        Nothing.

        OUTPUT:

        An element of the base ring.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: z = R(SR('z'))
            sage: z.singleton()
            z
            sage: z.singleton().parent()
            Symbolic Ring
        """
        return self._singleton_


    def size(self):
        """
        Return the size of the singleton.

        INPUT:

        Nothing.

        OUTPUT:

        The size of ``self``.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: z = R(SR('z'))
            sage: z.size()
            1
        """
        return self._size_


    #------------------------------------------------------------------------


    def _repr_main_(self, memo):
        """
        Returns a representation string.

        INPUT:

        - ``memo`` -- a dictionary.

        OUTPUT:

        A string.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R(SR('z'))._repr_main_({})
            'z'
        """
        self._update_memo_(memo)
        return repr(self.singleton())


    #------------------------------------------------------------------------


    def _cmp_without_operands_(left, right):
        """
        Return whether operators are equal or not (but without
        comparing the operands).

        INPUT:

        - ``left`` -- left expression.

        - ``right`` -- right expression.

        OUTPUT:

        ``True`` if operators are equal, otherwise ``False``.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericSingleton)
            sage: R = CombinatorialExpressionRing()
            sage: a = R(SR('a'))
            sage: b = R(SR('b'))
            sage: c = R(SR('c'), size=2)
            sage: d = GenericSingleton(R, SR('d'), size=1)
            sage: e = GenericSingleton(R, SR('e'), size=0)
            sage: f = R(SR('1'), size=0)
            sage: a == b  # indirect doctest
            True
            sage: a == c  # indirect doctest
            False
            sage: a == d  # indirect doctest
            True
            sage: a == e  # indirect doctest
            False
            sage: a == f  # indirect doctest
            False
            sage: c == d  # indirect doctest
            False
            sage: c == e  # indirect doctest
            False
            sage: c == f  # indirect doctest
            False
            sage: d == e  # indirect doctest
            False
            sage: e == f  # indirect doctest
            True
        """
        return (left.is_singleton() == right.is_singleton() and
                left.size() == right.size())


    #------------------------------------------------------------------------


    def is_singleton(self):
        """
        Return whether ``self`` is a singleton.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: B = R(SR('B'), function=True)
            sage: B.is_singleton()
            False
            sage: y = R(SR('y'))
            sage: y.is_singleton()
            True
            sage: e = R(SR(1), size=0)
            sage: e.is_singleton()
            True
        """
        return True


    def is_empty_singleton(self):
        """
        Return whether ``self`` is the empty singleton.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: B = R(SR('B'), function=True)
            sage: B.is_empty_singleton()
            False
            sage: y = R(SR('y'))
            sage: y.is_empty_singleton()
            False
            sage: e = R(SR(1), size=0)
            sage: e.is_empty_singleton()
            True
        """
        return self.is_singleton() and self.size() == Integer(0)


# ----------------------------------------------------------------------------


class UnlabeledSingleton(
    GenericSingleton,
    UnlabeledExpression):
    """
    A class representing an unlabeled singleton.

    For details see :class:`GenericSingleton`.

    .. SEEALSO::

        :class:`GenericSingleton`,
        :class:`LabeledSingleton'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`.
    """
    pass


# ----------------------------------------------------------------------------


class LabeledSingleton(
    GenericSingleton,
    LabeledExpression):
    """
    A class representing a labeled singleton.

    For details see :class:`GenericSingleton`.

    .. SEEALSO::

        :class:`GenericSingleton`,
        :class:`UnlabeledSingleton`.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`.
    """
    pass


#*****************************************************************************
# Element: Atom
#*****************************************************************************


class GenericAtom(GenericSingleton):
    """
    A class representing an atom (a singleton of size 1).

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``atom`` -- the atom.

    OUTPUT:

    A new atom.

    .. SEEALSO::

        :class:`UnlabeledAtom`,
        :class:`LabeledAtom'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`,

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     GenericAtom)
        sage: R = CombinatorialExpressionRing()
        sage: GenericAtom(R, SR('z'))
        z
        sage: GenericAtom(R, SR(1))
        1
        sage: y = R(SR('y')); y
        y
        sage: type(y)
        <class 'sage.combinat.combinatorial_expression.GenericAtom_with_category'>
    """
    def __init__(self, parent, atom):
        """
        See :class:`GenericAtom` for more information.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericAtom)
            sage: R = CombinatorialExpressionRing()
            sage: z = GenericAtom(R, SR('z')); z  # indirect doctest
            z
        """
        super(GenericAtom, self).__init__(parent, atom, Integer(1))


# ----------------------------------------------------------------------------


class UnlabeledAtom(
    GenericAtom,
    UnlabeledSingleton):
    """
    A class representing an unlabeled singleton.

    For details see :class:`GenericAtom`.

    .. SEEALSO::

        :class:`GenericAtom`,
        :class:`LabeledAtom'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.
    """
    pass


# ----------------------------------------------------------------------------


class LabeledAtom(
    GenericAtom,
    LabeledSingleton):
    """
    A class representing a labeled singleton.

    For details see :class:`GenericAtom`.

    .. SEEALSO::

        :class:`GenericAtom`,
        :class:`UnlabeledAtom`.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.
    """
    pass


#*****************************************************************************
# Element: Empty
#*****************************************************************************


class GenericEmpty(_EmptyFlavor_, GenericSingleton):
    """
    A class representing an empty singleton (a singleton of size 0).

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``empty`` -- the empty singleton.

    OUTPUT:

    A new empty singleton.

    .. NOTE::

        The empty expression has all flavors, i.e., it is unlabeled,
        as well as labeled at the same time.

    .. SEEALSO::

        :class:`UnlabeledEmpty`,
        :class:`LabeledEmpty'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     GenericEmpty)
        sage: R = CombinatorialExpressionRing()
        sage: GenericEmpty(R, SR('z'))
        z
        sage: GenericEmpty(R, SR(1))
        1
        sage: type(R(SR(1), size=0))
        <class 'sage.combinat.combinatorial_expression.GenericEmpty_with_category'>
        sage: type(R(SR(1), size=0, unlabeled=True))
        <class 'sage.combinat.combinatorial_expression.GenericEmpty_with_category'>
        sage: type(R(SR(1), size=0, labeled=True))
        <class 'sage.combinat.combinatorial_expression.GenericEmpty_with_category'>
        sage: R(SR(1), size=0) == R(SR(1), size=0, unlabeled=True) == R(SR(1), size=0, labeled=True)
        True
    """
    def __init__(self, parent, empty):
        """
        See :class:`GenericEmpty` for more information.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericEmpty)
            sage: R = CombinatorialExpressionRing()
            sage: z = GenericEmpty(R, SR(1)); z  # indirect doctest
            1
        """
        super(GenericEmpty, self).__init__(parent, empty, Integer(0))


# ----------------------------------------------------------------------------


UnlabeledEmpty = GenericEmpty
LabeledEmpty = GenericEmpty


#*****************************************************************************
# Element: DisjointUnion
#*****************************************************************************


class GenericDisjointUnion(GenericExpression):
    """
    A class representing a disjoint union of combinatorial expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A disjoint union of the operands.

    .. SEEALSO::

        :class:`UnlabeledDisjointUnion`,
        :class:`LabeledDisjointUnion'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     GenericDisjointUnion,
        ....:     GenericAtom)
        sage: R = CombinatorialExpressionRing()
        sage: GenericDisjointUnion(
        ....:     R, GenericAtom(R, SR(1)), GenericAtom(R, SR('z')))
        1 + z
    """

    def _repr_join_(self, reprs_of_operands):
        """
        Joins the given representation strings together.

        INPUT:

        - ``reprs_of_operands`` -- an iterable of strings.

        OUTPUT:

        A string.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R(SR('a')) + R(SR('b'))  # indirect doctest
            a + b
        """
        return ' + '.join(reprs_of_operands)


    def _repr_need_parenthesis_(self, operand):
        """
        Returns if parenthesis are needed around the representation string
        of the given operand.

        INPUT:

        - ``operand`` -- the (inner) operand.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: (R(SR('a')) * R(SR('b'))) + R(SR('c'))  # indirect doctest
            a*b + c
            sage: (R(SR('a')) + R(SR('b'))) * R(SR('c'))  # indirect doctest
            (a + b)*c
        """
        return isinstance(operand, GenericDisjointUnion)


# ----------------------------------------------------------------------------


class UnlabeledDisjointUnion(
    GenericDisjointUnion,
    UnlabeledExpression):
    """
    A class representing an unlabeled disjoint union of combinatorial
    expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    An unlabeled disjoint union of the operands.

    .. SEEALSO::

        :class:`GenericDisjointUnion`,
        :class:`LabeledDisjointUnion'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     UnlabeledDisjointUnion,
        ....:     UnlabeledAtom)
        sage: R = CombinatorialExpressionRing()
        sage: UnlabeledDisjointUnion(
        ....:     R, UnlabeledAtom(R, SR(1)), UnlabeledAtom(R, SR('z')))
        1 + z
    """
    pass


# ----------------------------------------------------------------------------


class LabeledDisjointUnion(
    GenericDisjointUnion,
    LabeledExpression):
    """
    A class representing a labeled disjoint union of combinatorial
    expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A labeled disjoint union of the operands.

    .. SEEALSO::

        :class:`GenericDisjointUnion`,
        :class:`UnlabeledDisjointUnion`.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     LabeledDisjointUnion,
        ....:     LabeledAtom)
        sage: R = CombinatorialExpressionRing()
        sage: LabeledDisjointUnion(
        ....:     R, LabeledAtom(R, SR(1)), LabeledAtom(R, SR('z')))
        1 + z
    """
    pass


#*****************************************************************************
# Element: CartesianProduct
#*****************************************************************************


class GenericCartesianProduct(GenericExpression):
    """
    A class representing a cartesian product of combinatorial expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A cartesian product of the operands.

    .. SEEALSO::

        :class:`UnlabeledCartesianProduct`,
        :class:`LabeledCartesianProduct'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     GenericCartesianProduct,
        ....:     GenericDisjointUnion,
        ....:     GenericAtom)
        sage: R = CombinatorialExpressionRing()
        sage: GenericCartesianProduct(
        ....:     R, GenericAtom(R, SR(1)), GenericAtom(R, SR('z')))
        1*z
        sage: GenericCartesianProduct(
        ....:     R, GenericAtom(R, SR(1)), GenericDisjointUnion(
        ....:         R, GenericAtom(R, SR('y')), GenericAtom(R, SR('z'))))
        1*(y + z)
    """

    def _repr_join_(self, reprs_of_operands):
        """
        Join the given representation strings together.

        INPUT:

        - ``reprs_of_operands`` -- an iterable of strings.

        OUTPUT:

        A string.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R(SR('a')) * R(SR('b'))  # indirect doctest
            a*b
        """
        return '*'.join(reprs_of_operands)


    def _repr_need_parenthesis_(self, operand):
        """
        Return if parenthesis are needed around the representation string
        of the given operand.

        INPUT:

        - ``operand`` -- the (inner) operand.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: (R(SR('a')) * R(SR('b'))) + R(SR('c'))  # indirect doctest
            a*b + c
            sage: (R(SR('a')) + R(SR('b'))) * R(SR('c'))  # indirect doctest
            (a + b)*c
        """
        return isinstance(operand, GenericDisjointUnion)


# ----------------------------------------------------------------------------


class UnlabeledCartesianProduct(
    GenericCartesianProduct,
    UnlabeledExpression):
    """
    A class representing an unlabeled cartesian product of
    combinatorial expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    An unlabeled cartesian product of the operands.

    .. SEEALSO::

        :class:`GenericCartesianProduct`,
        :class:`LabeledCartesianProduct'.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     UnlabeledCartesianProduct,
        ....:     UnlabeledAtom)
        sage: R = CombinatorialExpressionRing()
        sage: UnlabeledCartesianProduct(
        ....:     R, UnlabeledAtom(R, SR(1)), UnlabeledAtom(R, SR('z')))
        1*z
    """
    pass


# ----------------------------------------------------------------------------


class LabeledCartesianProduct(
    GenericCartesianProduct,
    LabeledExpression):
    """
    A class representing a labeled cartesian product of
    combinatorial expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A labeled cartesian product of the operands.

    .. SEEALSO::

        :class:`GenericCartesianProduct`,
        :class:`UnlabeledCartesianProduct`.

        :class:`CombinatorialExpressionRing`,
        :class:`GenericAtom`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     LabeledCartesianProduct,
        ....:     LabeledAtom)
        sage: R = CombinatorialExpressionRing()
        sage: LabeledCartesianProduct(
        ....:     R, LabeledAtom(R, SR(1)), LabeledAtom(R, SR('z')))
        1*z
    """
    pass


#*****************************************************************************
# Element: TODO
#*****************************************************************************

# finite set ??
# multi set
# sequence


#*****************************************************************************
# Parent: CombinatorialExpressionRing
#*****************************************************************************


class CombinatorialExpressionRing(
    sage.structure.unique_representation.UniqueRepresentation,
    sage.rings.ring.Ring):
    """
    A ring of combinatorial expressions.

    INPUT:

    - ``base`` -- the base ring.

    OUTPUT:

    A ring.

    EXAMPLES::

        sage: R = CombinatorialExpressionRing(); R
        Combinatorial Expression Ring
        sage: R(SR('z'))
        z

    .. SEEALSO::

        :class:`GenericAtom`,
        :class:`GenericCartesianProduct`,
        :class:`GenericDisjointUnion`,
        :class:`GenericEmpty`,
        :class:`GenericExpression`,
        :class:`GenericFunction`,
        :class:`GenericSingleton`.
    """
    def __init__(self, base=None):
        """
        See :class:`CombinatorialExpressionRing` for more information.

        TESTS::

            sage: R = CombinatorialExpressionRing(); R
            Combinatorial Expression Ring
            sage: type(R)
            <class 'sage.combinat.combinatorial_expression.CombinatorialExpressionRing_with_category'>
            sage: R.category()
            Category of rings
            sage: R.element_class
            <class 'sage.combinat.combinatorial_expression.CombinatorialExpressionRing_with_category.element_class'>
        """
        if base is None:
            base = sage.symbolic.ring.SR
        if base != sage.symbolic.ring.SR:
            raise NotImplementedError("%s not allowed as base ring." % (base,))
        super(CombinatorialExpressionRing, self).__init__(base=base)

        self.Operators = Operators  # this should be static, no idea how...


    def _repr_(self):
        """
        Return a representation string.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: repr(CombinatorialExpressionRing())  # indirect doctest
            'Combinatorial Expression Ring'
        """
        return "Combinatorial Expression Ring"


    def base_ring(self):
        """
        Return the base ring.

        INPUT:

        Nothing.

        OUTPUT:

        A ring.

        EXAMPLES::

            sage: CombinatorialExpressionRing(base=SR).base_ring()
            Symbolic Ring
        """
        return self.base().base_ring()


    #------------------------------------------------------------------------


    Element = GenericExpression


    #------------------------------------------------------------------------


    @lazy_attribute
    def element_class_generic_function(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_generic_function.__name__
            'GenericFunction_with_category'
        """
        return self.__make_element_class__(GenericFunction)


    @lazy_attribute
    def element_class_unlabeled_function(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_unlabeled_function.__name__
            'UnlabeledFunction_with_category'
        """
        return self.__make_element_class__(UnlabeledFunction)


    @lazy_attribute
    def element_class_labeled_function(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_labeled_function.__name__
            'LabeledFunction_with_category'
        """
        return self.__make_element_class__(LabeledFunction)


    #------------------------------------------------------------------------


    @lazy_attribute
    def element_class_generic_empty(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_generic_empty.__name__
            'GenericEmpty_with_category'
        """
        return self.__make_element_class__(GenericEmpty)


    @lazy_attribute
    def element_class_unlabeled_empty(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_unlabeled_empty.__name__
            'GenericEmpty_with_category'
        """
        return self.__make_element_class__(UnlabeledEmpty)


    @lazy_attribute
    def element_class_labeled_empty(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_labeled_empty.__name__
            'GenericEmpty_with_category'
        """
        return self.__make_element_class__(LabeledEmpty)


    #------------------------------------------------------------------------


    @lazy_attribute
    def element_class_generic_atom(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_generic_atom.__name__
            'GenericAtom_with_category'
        """
        return self.__make_element_class__(GenericAtom)


    @lazy_attribute
    def element_class_unlabeled_atom(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_unlabeled_atom.__name__
            'UnlabeledAtom_with_category'
        """
        return self.__make_element_class__(UnlabeledAtom)


    @lazy_attribute
    def element_class_labeled_atom(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_labeled_atom.__name__
            'LabeledAtom_with_category'
        """
        return self.__make_element_class__(LabeledAtom)


    #------------------------------------------------------------------------


    @lazy_attribute
    def element_class_generic_singleton(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_generic_singleton.__name__
            'GenericSingleton_with_category'
        """
        return self.__make_element_class__(GenericSingleton)


    @lazy_attribute
    def element_class_unlabeled_singleton(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_unlabeled_singleton.__name__
            'UnlabeledSingleton_with_category'
        """
        return self.__make_element_class__(UnlabeledSingleton)


    @lazy_attribute
    def element_class_labeled_singleton(self):
        """
        The class of the specified element.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R.element_class_labeled_singleton.__name__
            'LabeledSingleton_with_category'
        """
        return self.__make_element_class__(LabeledSingleton)


    #------------------------------------------------------------------------


    def _get_element_class_(self, what, flavor):
        """
        Helper function for getting the element class.

        INPUT:

        - ``what`` -- classname.

        - ``flavor`` -- the desired flavor.

        OUTPUT:

        A class.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R._get_element_class_('singleton', 'labeled')
            <class 'sage.combinat.combinatorial_expression.LabeledSingleton_with_category'>
        """
        return getattr(self, 'element_class_' + flavor + '_' + what)


    #------------------------------------------------------------------------


    def _element_constructor_(self, *args, **kwargs):
        """
        Returns a combinatorial expression out of the given data.

        INPUT:

        - ``*args`` -- the data out of which we construct the
          combinatorial expression.

        - ``unlabeled`` -- build an unlabeled expression.

        - ``labeled`` -- build a labeled expression.

        - ``size`` -- a size (e.g. an integer) or ``None`` (default).

        - ``empty`` -- a boolean (default: ``False``) indicating
          whether result should be empty (i.e. a singleton of size `0`).

        - ``function`` -- a boolean (default: ``False``) indicating
          whether result should be a combinatorial function or not.

        OUTPUT:

        A combinatorial expression.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: R(SR('T'), function=True)  # indirect doctest
            T = None
            sage: z = R(SR('z'), labeled=True, size=1); z  # indirect doctest
            z
            sage: type(z)
            <class 'sage.combinat.combinatorial_expression.LabeledAtom_with_category'>
            sage: z.size()
            1
            sage: e = R(SR(1), size=0); e  # indirect doctest
            1
            sage: type(e)
            <class 'sage.combinat.combinatorial_expression.GenericEmpty_with_category'>
            sage: e.size()
            0
            sage: y = R(SR('y'), size=2); y  # indirect doctest
            y
            sage: type(y)
            <class 'sage.combinat.combinatorial_expression.GenericSingleton_with_category'>
            sage: y.size()
            2
            sage: R(SR('T'), function=True)
            T = None
        """
        empty = kwargs.pop('empty', not args)
        function = kwargs.pop('function', False)
        size = kwargs.pop('size', None)

        if size is not None and function:
            raise ValueError('A combinatorial function does not have a size, '
                             'but size=%s given' % (size,))
        if size is not None and size != Integer(0) and empty:
            raise ValueError('An empty expression has zero size, '
                             'but size=%s given' % (size,))
        if sum(1 for p in (empty, function) if p) > 1:
            raise ValueError('More than one of the keywords '
                             '"empty", "function" given.')

        labeled = kwargs.pop('labeled', None)
        unlabeled = kwargs.pop('unlabeled', None)

        if unlabeled and labeled and not empty and size != Integer(0):
            raise ValueError('An (non-empty) expression cannot be unlabeled '
                             'and labeled at the same time.')
        if unlabeled and labeled:
            flavor = 'generic'
        elif unlabeled:
            flavor = 'unlabeled'
        elif labeled:
            flavor = 'labeled'
        else:
            flavor = 'generic'

        if not args:
            if hasattr(self.base_ring(), 'one'):
                data = self.base_ring().one()
            elif hasattr(self.base_ring(), 'zero'):
                data = self.base_ring().zero()
            else:
                raise ValueError('No argument given '
                                 '(and cannot find something appropriate).')
        elif len(args) > 1:
            raise ValueError('Too many input arguments given.')
        else:
            data = args[0]

        if isinstance(data, GenericExpression):
            raise NotImplementedError('We should make a copy...')

        if empty or size == 0:
            return self._get_element_class_('empty', flavor)(self, data)

        if size == 1:
            return self._get_element_class_('atom', flavor)(self, data)

        if size:
            return self._get_element_class_('singleton', flavor)(
                self, data, size=size)

        if function:
            return self._get_element_class_('function', flavor)(self, data)

        # At this point, we know that we have processed what user
        # wanted, so from now on we try to find out, what the user
        # might wanted ;)

        try:
            P = data.parent()
        except AttributeError:
            P = None

        if P == sage.symbolic.ring.SR:
            if data.operator() is None:
                return self._get_element_class_('atom', flavor)(self, data)
            raise NotImplementedError('At the moment only symbolic '
                                      'expressions without '
                                      'operators can be converted to '
                                      'combinatorial expressions directly.')

        # we have something else, so we interpret it as an atom
        return self._get_element_class_('atom', flavor)(self, data)


    #------------------------------------------------------------------------


    def _an_element_(self):
        """
        Returns a combinatorial expression.

        INPUT:

        Nothing.

        OUTPUT:

        A combinatorial expression.

        TESTS::

            sage: R = CombinatorialExpressionRing(base=SR)
            sage: s = R._an_element_(); s
            some_variable
            sage: type(s)
            <class 'sage.combinat.combinatorial_expression.GenericAtom'>
            sage: s.parent() == R
            True
        """
        return GenericAtom(self, self.base_ring().an_element())


#*****************************************************************************
# Operators
#*****************************************************************************


class Operators(sage.structure.sage_object.SageObject):
    """
    This class contains a collection of operators used to create
    combinatorial expressions.
    """

    @staticmethod
    def _apply_operator_(operatorname, *operands, **kwargs):
        """
        Calls the operator of the first operand with other operands as
        arguments.

        INPUT:

        - ``operatorname`` -- a string.

        - ``*operands`` -- the arguments.

        - ``parent`` -- the desired parent of the result. If this is
          ``None``, then the parent is determined automatically.

        OUTPUT:

        A combinatorial expression.

        TESTS::

            sage: R = CombinatorialExpressionRing()
            sage: y = R(SR('y'))
            sage: z = R(SR('z'))
            sage: R.Operators._apply_operator_(
            ....:     '_disjoint_union_', y, z)
            y + z
        """
        parent = kwargs.pop('parent', None)  # Python 2.x workaround,
                                             # since parent=None
                                             # CANNOT be used AFTER
                                             # *operands.
        if not operands:
            raise NotImplementedError('Operands must be given.')
        # TODO: do coercion here
        return getattr(operands[0], operatorname)(*operands[1:])


    #------------------------------------------------------------------------


    @classmethod
    def disjoint_union(cls, *operands, **kwargs):
        """
        Returns the disjoint union of ``self`` and the given operands.

        INPUT:

        - ``*others`` -- operands.

        OUTPUT:

        A new combinatorial expression.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: y = R(SR('y'))
            sage: z = R(SR('z'))
            sage: a = y.disjoint_union(z); a
            y + z
            sage: b = R.Operators.disjoint_union(y, z); b
            y + z
            sage: c = y + z; c
            y + z
            sage: a == b == c
            True

        .. SEEALSO::

            :meth:`GenericExpression._mul_`,
            :meth:`GenericExpression.disjoint_union`.

            :meth:`disjoint_union`.
        """
        parent = kwargs.pop('parent', None)
        return cls._apply_operator_('_disjoint_union_',
                                    *operands,
                                    parent=parent)


    #------------------------------------------------------------------------


    @classmethod
    def cartesian_product(cls, *operands, **kwargs):
        """
        Returns the cartesian product of ``self`` and the given operands.

        INPUT:

        - ``*others`` -- operands.

        OUTPUT:

        A new combinatorial expression.

        EXAMPLES::

            sage: R = CombinatorialExpressionRing()
            sage: y = R(SR('y'))
            sage: z = R(SR('z'))
            sage: a = y.cartesian_product(z); a
            y*z
            sage: b = R.Operators.cartesian_product(y, z); b
            y*z
            sage: c = y * z; c
            y*z
            sage: a == b == c
            True

        .. SEEALSO::

            :meth:`GenericExpression._mul_`,
            :meth:`GenericExpression.cartesian_product`.

            :meth:`disjoint_union`.
        """
        parent = kwargs.pop('parent', None)
        return cls._apply_operator_('_cartesian_product_',
                                    *operands,
                                    parent=parent)


# ----------------------------------------------------------------------------

