# -*- coding: utf-8 -*-
r"""

AUTHORS:

- Daniel Krenn (2014-04-01): initial version

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
from sage.structure.sage_object import SageObject  # TODO löschen
from itertools import izip

#*****************************************************************************
# Flavor
#*****************************************************************************

class _GenericFlavor_(sage.structure.sage_object.SageObject):

q    def is_labeled(self):
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


    @staticmethod
    def _class_with_prefix_(prefix, classname):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return globals()[prefix + classname]


    @classmethod
    def class_generic(cls, classname):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return cls._class_with_prefix_('Generic', classname)


    @classmethod
    def class_unlabeled(cls, classname):
        return cls._class_with_prefix_('Unlabeled', classname)


    @classmethod
    def class_labeled(cls, classname):
        return cls._class_with_prefix_('Labeled', classname)


# ----------------------------------------------------------------------------


class _UnlabeledFlavor_(_GenericFlavor_):

    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``False`` since this instance is unlabeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _UnlabeledFlavor_)
            sage: F = _UnlabeledFlavor_()
            sage: F.is_labeled()
            False
        """
        return False


    def is_unlabeled(self):
        """
        Returns whether combinatorial expression is unlabeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``True`` since this instance is unlabeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _UnlabeledFlavor_)
            sage: F = _UnlabeledFlavor_()
            sage: F.is_unlabeled()
            True
        """
        return True


# ----------------------------------------------------------------------------


class _LabeledFlavor_(_GenericFlavor_):

    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``True`` since this instance is labeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _LabeledFlavor_)
            sage: F = _LabeledFlavor_()
            sage: F.is_labeled()
            True
        """
        return True


    def is_unlabeled(self):
        """
        Returns whether combinatorial expression is unlabeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``False`` since this instance is labeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _LabeledFlavor_)
            sage: F = _LabeledFlavor_()
            sage: F.is_unlabeled()
            False
        """
        return False


# ----------------------------------------------------------------------------


class _EmptyFlavor_(_GenericFlavor_):

    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing

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

        Nothing

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


class GenericBase(
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

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     GenericExpression,
        ....:     GenericDisjointUnion,
        ....:     GenericCartesianProduct,
        ....:     GenericAtom)
        sage: R = CombinatorialExpressionRing(SR)
        sage: z = GenericAtom(R, var('z'))
        sage: z
        z
        sage: eps = GenericAtom(R, SR(1))
        sage: eps
        1
        sage: T = GenericExpression(R, var('T')); T
        T == None
        sage: A = GenericCartesianProduct(R, z, T, T); A
        z*T*T
        sage: B = GenericDisjointUnion(R, eps, A); B
        1 + z*T*T
        sage: T.assign(B)
        sage: T
        T == 1 + z*T*T
    """
    def __init__(self, parent, *operands):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        self.assign(*operands)
        super(GenericBase, self).__init__(parent)


    # ------------------------------------------------------------------------


    def assign(self, *operands):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        self._operands_ = tuple(operands)


    def operands(self):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return self._operands_


    def iter_operands(self):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return iter(self._operands_)


    #------------------------------------------------------------------------


    def _name_(self):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return self.__class__.__name__


    #------------------------------------------------------------------------


    def _update_memo_(self, memo):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        #if not memo:
        #    memo = {}
        key = id(self)
        m = memo.get(key)
        memo[key] = True
        return m is None


    #------------------------------------------------------------------------


    def _repr_(self):
        """
        Returns a representation string.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: TODO  # not tested
        """
        memo = {}
        return self._repr_main_(memo)


    def _repr_main_(self, memo):
        """
        Returns a repesentation string.

        INPUT:

        - ``memo`` -- a dictionary.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: TODO  # not tested
        """
        return self._repr_recursive_(memo)


    def _repr_recursive_(self, memo):
        """
        Recursively builds the representation string

        INPUT:

        - ``memo`` -- a dictionary.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: TODO  # not tested
        """
        if not self._update_memo_(memo):
            return self._repr_recursion_()
        return self._repr_join_(
            self._repr_make_parenthesis_(
                o._repr_main_(memo), self._repr_need_parenthesis_(o))
            for o in self.operands())


    @staticmethod
    def _repr_make_parenthesis_(s, make=True):
        """
        Makes parenthesis around a string.

        INPUT:

        - ``s`` -- a string.

        - ``make`` -- a boolean indicating whether to make parenthesis
          or not.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     GenericBase)
            sage: GenericBase._repr_make_parenthesis_(
            ....:     '42', False)
            '42'
            sage: GenericBase._repr_make_parenthesis_(
            ....:     '42', True)
            '(42)'
            """
        if make:
            return '(' + s + ')'
        else:
            return s


    def _repr_need_parenthesis_(self, operand):
        """
        Returns if parenthesis are needed around the representation string
        of the given operand.

        INPUT:

        - ``operand`` -- the (inner) operand.

        OUTPUT:

        ``False``.
        """
        return False


    def _repr_join_(self, reprs_of_operands):
        """
        Joins the given representation strings together.

        INPUT:

        - ``reprs_of_operands`` -- an iterable of strings.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: TODO  # not tested
        """
        return (self._name_() +
                self._repr_make_parenthesis_(
                ', '.join(o for o in reprs_of_operands)))


    def _repr_recursion_(self):
        """
        Returns a string signalling that an recursion was encountered.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: TODO  # not tested
        """
        return 'REC'


    #------------------------------------------------------------------------


    def __cmp__(left, right, memo=None):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        if not memo:
            memo = {}
        key = (id(left), id(right))
        m = memo.get(key)
        if m is not None:
            return m
        memo[key] = True

        if left.__class__ != right.__class__:
            result = False
        else:
            S = left.operands()
            O = right.operands()
            if len(S) != len(O):
                result = False
            else:
                result = all(s.__cmp__(o, memo) for s, o in izip(S, O))

        memo[key] = result
        return result


    #------------------------------------------------------------------------


    def _apply_operator_(self, operatorclassname, *operands):
        """
        Returns an instance of the specified class (as operator) with
        given operands.

        INPUT:

        - ``operatorclassname`` -- the name of the class (without flavor).

        - ``*operands`` -- the operands of the operation.

        OUTPUT:

        A new combinatorial expression.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     CombinatorialExpressionRing, Operators)
            sage: R = CombinatorialExpressionRing(SR)
            sage: g = R(var('g'))
            sage: u = R(var('u'), unlabeled=True)
            sage: l = R(var('l'), labeled=True)
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
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return Operators.disjoint_union(self, *others)


    def _disjoint_union_(self, *others):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return self._apply_operator_('DisjointUnion', self, *others)


    def _add_(left, right):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return left.disjoint_union(right)


    #------------------------------------------------------------------------


    def cartesian_product(self, *others):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return Operators.cartesian_product(self, *others)


    def _cartesian_product_(self, *others):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return self._apply_operator_('CartesianProduct', self, *others)


    def _mul_(left, right):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return left.cartesian_product(right)


    #------------------------------------------------------------------------


    def iter_elements(self, size):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        # do preprocessing of size here
        return self._iter_elements_(size)


    def _iter_elements_(self, size):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return NotImplementedError('Iteration for an instance of %s is '
                                   'not implemented.' % (self._name_(),))


    #------------------------------------------------------------------------


    def random_element(self, size):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        # do preprocessing of size here
        return self._random_element_(size)


    def _random_element_(self, size):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return NotImplementedError('Generation of an random element '
                                   'for an instance of %s '
                                   'not implemented.' % (self._name_(),))


# ----------------------------------------------------------------------------


class UnlabeledBase(
    GenericBase,
    _UnlabeledFlavor_):
    pass


# ----------------------------------------------------------------------------


class LabeledBase(
    GenericBase,
    _LabeledFlavor_):
    pass


#*****************************************************************************
# Element: Expression
#*****************************************************************************


class GenericExpression(GenericBase):
    """
    A class representing a combinatorial expression.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A new combinatorial expression.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     GenericExpression)
        sage: R = CombinatorialExpressionRing(SR)
        sage: GenericExpression(R, var('T'))
        T == None
    """
    def __init__(self, parent, expression, *operands):
        super(GenericExpression, self).__init__(parent, *operands)
        self._expression_ = expression


    def _repr_(self):
        memo = {}
        s = repr(self._expression_) + ' == '
        if self._operands_:
            return s + self._repr_recursive_(memo)
        else:
            return s + 'None'


    def _repr_main_(self, memo):
        self._update_memo_(memo)
        return repr(self._expression_)


    def _repr_join_(self, reprs_of_operands):
        """
        Joins the given representation strings together.

        INPUT:

        - ``reprs_of_operands`` -- an iterable of strings.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: TODO  # not tested
        """
        return ', '.join(o for o in reprs_of_operands)

# ----------------------------------------------------------------------------


class UnlabeledExpression(
    GenericExpression,
    UnlabeledBase):
    """
    A class representing a unlabeled combinatorial expression.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A new unlabeled combinatorial expression.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     UnlabeledExpression)
        sage: R = CombinatorialExpressionRing(SR)
        sage: T = UnlabeledExpression(R, var('T')); T
        T == None
        sage: T.is_unlabeled()
        True
        sage: T.is_labeled()
        False
    """
    pass


# ----------------------------------------------------------------------------

class LabeledExpression(
    GenericExpression,
    LabeledBase):
    """
    A class representing a labeled combinatorial expression.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A new labeled combinatorial expression.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     LabeledExpression)
        sage: R = CombinatorialExpressionRing(SR)
        sage: LabeledExpression(R, var('T'))
        T == None
    """
    pass


#*****************************************************************************
# Element: Singleton
#*****************************************************************************

class GenericSingleton(GenericBase):
    """
    A class representing a singleton of given size.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``singleton`` -- the singleton.

    - ``size`` -- the size of the singleton.

    OUTPUT:

    A new singleton.

    EXAMPLES::

        sage: TODO  # not tested
    """

    def __init__(self, parent, singleton, size):
        """
        TODO

        TESTS::

            sage: TODO  # not tested
        """
        super(GenericSingleton, self).__init__(parent)
        self._set_singleton_(singleton, size)


    def _set_singleton_(self, singleton, size):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self._singleton_ = singleton
        self._size_ = size


    def size(self):
        return self._size_


    def _repr_main_(self, memo):
        self._update_memo_(memo)
        return repr(self._singleton_)


# ----------------------------------------------------------------------------


class UnlabeledSingleton(
    GenericSingleton,
    UnlabeledBase):
    pass


# ----------------------------------------------------------------------------


class LabeledSingleton(
    GenericSingleton,
    LabeledBase):
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

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     GenericAtom)
        sage: R = CombinatorialExpressionRing(SR)
        sage: GenericAtom(R, var('z'))
        z
        sage: GenericAtom(R, SR(1))
        1
        sage: y = R(var('y')); y
        y
        sage: type(y)
        <class 'sage.combinat.combinatorial_expression.GenericAtom_with_category'>
    """
    def __init__(self, parent, atom):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        super(GenericAtom, self).__init__(parent, atom, Integer(1))


# ----------------------------------------------------------------------------


class UnlabeledAtom(
    GenericAtom,
    UnlabeledSingleton):
    pass


# ----------------------------------------------------------------------------


class LabeledAtom(
    GenericAtom,
    LabeledSingleton):
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

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     GenericEmpty)
        sage: R = CombinatorialExpressionRing(SR)
        sage: GenericEmpty(R, var('z'))
        z
        sage: GenericEmpty(R, SR(1))
        1
        sage: type(R(SR(1)))
        <class 'sage.combinat.combinatorial_expression.GenericEmpty_with_category'>
        sage: type(R(SR(1), unlabeled=True))
        <class 'sage.combinat.combinatorial_expression.GenericEmpty_with_category'>
        sage: type(R(SR(1), labeled=True))
        <class 'sage.combinat.combinatorial_expression.GenericEmpty_with_category'>
        sage: R(SR(1)) == R(SR(1), unlabeled=True) == R(SR(1), labeled=True)
    """
    def __init__(self, parent, empty):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        super(GenericEmpty, self).__init__(parent, empty, Integer(0))


# ----------------------------------------------------------------------------


UnlabeledEmpty = GenericEmpty
LabeledEmpty = GenericEmpty


#*****************************************************************************
# Element: DisjointUnion
#*****************************************************************************


class GenericDisjointUnion(GenericBase):
    """
    A class representing a disjoint union of combinatorial expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A disjoint union of the operands.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     GenericDisjointUnion,
        ....:     GenericAtom)
        sage: R = CombinatorialExpressionRing(SR)
        sage: GenericDisjointUnion(
        ....:     R, GenericAtom(R, SR(1)), GenericAtom(R, var('z')))
        1 + z
    """
    def _repr_join_(self, reprs_of_operands):
        """
        Joins the given representation strings together.

        INPUT:

        - ``reprs_of_operands`` -- an iterable of strings.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: TODO  # not tested
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
        """
        return isinstance(operand, GenericDisjointUnion)


# ----------------------------------------------------------------------------


class UnlabeledDisjointUnion(
    GenericDisjointUnion,
    UnlabeledBase):
    """
    A class representing an unlabeled disjoint union of combinatorial
    expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    An unlabeled disjoint union of the operands.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     UnlabeledDisjointUnion,
        ....:     UnlabeledAtom)
        sage: R = CombinatorialExpressionRing(SR)
        sage: UnlabeledDisjointUnion(
        ....:     R, UnlabeledAtom(R, SR(1)), UnlabeledAtom(R, var('z')))
        1 + z
    """
    pass


# ----------------------------------------------------------------------------


class LabeledDisjointUnion(
    GenericDisjointUnion,
    LabeledBase):
    """
    A class representing a labeled disjoint union of combinatorial
    expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A labeled disjoint union of the operands.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     LabeledDisjointUnion,
        ....:     LabeledAtom)
        sage: R = CombinatorialExpressionRing(SR)
        sage: LabeledDisjointUnion(
        ....:     R, LabeledAtom(R, SR(1)), LabeledAtom(R, var('z')))
        1 + z
    """
    pass


#*****************************************************************************
# Element: CartesianProduct
#*****************************************************************************


class GenericCartesianProduct(GenericBase):
    """
    A class representing a cartesian product of combinatorial expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A cartesian product of the operands.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     GenericCartesianProduct,
        ....:     GenericDisjointUnion,
        ....:     GenericAtom)
        sage: R = CombinatorialExpressionRing(SR)
        sage: GenericCartesianProduct(
        ....:     R, GenericAtom(R, SR(1)), GenericAtom(R, var('z')))
        1*z
        sage: GenericCartesianProduct(
        ....:     R, GenericAtom(R, SR(1)), GenericDisjointUnion(
        ....:         R, GenericAtom(R, var('y')), GenericAtom(R, var('z'))))
        1*(y + z)
    """
    def _repr_join_(self, reprs_of_operands):
        """
        Joins the given representation strings together.

        INPUT:

        - ``reprs_of_operands`` -- an iterable of strings.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: TODO  # not tested
        """
        return '*'.join(reprs_of_operands)


    def _repr_need_parenthesis_(self, operand):
        """
        Returns if parenthesis are needed around the representation string
        of the given operand.

        INPUT:

        - ``operand`` -- the (inner) operand.

        OUTPUT:

        ``True`` or ``False``.
        """
        return isinstance(operand, GenericDisjointUnion)


# ----------------------------------------------------------------------------


class UnlabeledCartesianProduct(
    GenericCartesianProduct,
    UnlabeledBase):
    """
    A class representing an unlabeled cartesian product of
    combinatorial expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    An unlabeled cartesian product of the operands.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     UnlabeledCartesianProduct,
        ....:     UnlabeledAtom)
        sage: R = CombinatorialExpressionRing(SR)
        sage: UnlabeledCartesianProduct(
        ....:     R, UnlabeledAtom(R, SR(1)), UnlabeledAtom(R, var('z')))
        1*z
    """
    pass


# ----------------------------------------------------------------------------


class LabeledCartesianProduct(
    GenericCartesianProduct,
    LabeledBase):
    """
    A class representing a labeled cartesian product of
    combinatorial expressions.

    INPUT:

    - ``parent`` -- the parent of the new object.

    - ``*operands`` -- the operands of the operation defined by this class.

    OUTPUT:

    A labeled cartesian product of the operands.

    TESTS::

        sage: from sage.combinat.combinatorial_expression import (
        ....:     CombinatorialExpressionRing,
        ....:     LabeledCartesianProduct,
        ....:     LabeledAtom)
        sage: R = CombinatorialExpressionRing(SR)
        sage: LabeledCartesianProduct(
        ....:     R, LabeledAtom(R, SR(1)), LabeledAtom(R, var('z')))
        1*z
    """
    pass


#*****************************************************************************
# Element: TODO
#*****************************************************************************

# finite set ??
# multi set
# sequence

