r"""
Groups of elements on the unit circle (e.g. roots of unity).

AUTHORS:

- Daniel Krenn (2018)

Classes and Methods
===================
"""
#*****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#               2014--2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.structure.element import MultiplicativeGroupElement
from sage.structure.factory import UniqueFactory
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp_by_eq_and_lt
from sage.structure.unique_representation import UniqueRepresentation


class AbstractArgument(MultiplicativeGroupElement):
    r"""
    An element of :class:`AbstractArgumentGroup`. This abstract class
    encapsulates an element of the parent's base, i.e. it can be seen
    as a wrapper class.

    INPUT:

    - ``parent`` -- a SageMath parent

    - ``element`` -- an element of parent's base

    - ``normalize`` -- a boolean (default: ``True``)
    """

    def __init__(self, parent, element, normalize=True):
        r"""
        See :class:`AbstractArgument` for more information.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=1/2)  # indirect doctest
            e^(2*pi*0.500000000000000)
            sage: C(exponent=3/2)
            e^(2*pi*0.500000000000000)

            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0)
            1
            sage: U(exponent=1)
            1
            sage: U(exponent=2/3) == U(exponent=5/3)
            True
        """
        if parent is None:
            raise ValueError('parent must be provided')
        super(AbstractArgument, self).__init__(parent=parent)

        element = parent.base()(element)
        if normalize:
            element = self._normalize_(element)
        self._element_ = element

    @staticmethod
    def _normalize_(element):
        r"""
        Normalizes the given element.

        INPUT:

        - ``element`` -- an element of the parent's base

        OUTPUT:

        An element.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import AbstractArgument
            sage: AbstractArgument._normalize_(3/2)
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented in concrete realizations
        """
        raise NotImplementedError('only implemented in concrete realizations')

    def __hash__(self):
        r"""
        Return a hash value of this argument.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: hash(C(exponent=1/3))  # indirect doctest, random
            42
        """
        return hash((self.parent(), self._element_))

    _richcmp_ = richcmp_by_eq_and_lt("_eq_", "_lt_")

    def _eq_(self, other):
        r"""
        Return whether this argument equals ``other``.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0) == U(exponent=1)
            True
            sage: U(exponent=2/3) == U(exponent=5/3)
            True
            sage: U(exponent=2/3) == U(exponent=-2/3)
            False

        ::

            sage: from sage.groups.roots_of_unity_group import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(I) == C(I)
            True

        As we do not have normalization in :class:`ArgumentByElement`,
        then following, although equal, is not equal::

            sage: C(I) == C(2*I)
            False
        """
        return self._element_ == other._element_

    def _lt_(self, other):
        r"""
        Raise an error since points on the unit circle are not comparable.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0) < U(exponent=0)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for the roots of unity 1 and 1
            sage: U(exponent=0) < U(exponent=1/2)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for the roots of unity 1 and -1
            sage: U(exponent=0) > U(exponent=1/2)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for the roots of unity -1 and 1
        """
        raise RuntimeError("cannot decide '<' "
                           "for the roots of unity "
                           "{} and {}".format(self, other))

class AbstractArgumentGroup(UniqueRepresentation, Parent):
    r"""
    A group whose elements represent (complex) arguments.

    INPUT:

    - ``base`` -- a SageMath parent

    - ``category`` -- a category
    """

    Element = AbstractArgument

    @staticmethod
    def __classcall__(cls, base, category=None):
        r"""
        See :class:`AbstractArgumentGroup` for more information.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup(RR).category()  # indirect doctest
            Category of commutative groups
        """
        category = cls._determine_category_(category)
        return super(AbstractArgumentGroup, cls).__classcall__(
            cls, base, category)

    @staticmethod
    def _determine_category_(category):
        r"""
        Return the category of this argument group.

        INPUT:

        - ``category`` -- a category or ``None`` (in which case the output
          equals ``category``)

        OUTPUT:

        A category.

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup._determine_category_(None)
            Category of commutative groups
            sage: UnitCircleGroup._determine_category_(Groups())
            Category of groups
        """
        if category is None:
            from sage.categories.groups import Groups
            category = Groups().Commutative()
        return category

    def __init__(self, base, category):
        r"""
        See :class:`AbstractArgumentGroup` for more information.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup(RR).base()  # indirect doctest
            Real Field with 53 bits of precision
        """
        super(AbstractArgumentGroup, self).__init__(category=category,
                                              base=base)

    def __hash__(self):
        r"""
        Return a hash value of this argument group.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: hash(UnitCircleGroup(RR))  # indirect doctest, random
            42
        """
        return hash((self.__class__, self.base()))

    def _an_element_(self):
        r"""
        Return an element of this argument group.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup(RR).an_element()  # indirect doctest
            e^(2*pi*0.000000000000000)
        """
        return self.element_class(self, self.base().an_element())


class UnitCirclePoint(AbstractArgument):
    r"""
    An element of :class:`UnitCircleGroup`
    which is `e^{2\pi\mathit{exponent}}`.

    INPUT:

    - ``parent`` -- a SageMath parent

    - ``exponent`` -- a number (of a subset of the reals)

    - ``normalize`` -- a boolean (default: ``True``)
    """

    @staticmethod
    def _normalize_(exponent):
        r"""
        Normalizes the given exponent so that it is in `[0,1)`.

        INPUT:

        - ``exponent`` -- an element of the parent's base

        OUTPUT:

        An element.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCirclePoint
            sage: UnitCirclePoint._normalize_(3/2)
            1/2
        """
        return exponent - exponent.floor()

    @property
    def exponent(self):
        r"""
        The exponent of this point on the unit circle.

        EXAMPLES:

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=4/3).exponent
            0.333333333333333
        """
        return self._element_

    def _repr_(self):
        r"""
        Return a representation string of this point on the unit circle.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=1/3)
            e^(2*pi*0.333333333333333)
        """
        return 'e^(2*pi*{})'.format(self.exponent)

    def _mul_(self, other):
        r"""
        Return the product of this point on the unit circle and ``other``.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=0.3) * C(exponent=0.4)
            e^(2*pi*0.700000000000000)
        """
        P = self.parent()
        return P.element_class(P, self.exponent + other.exponent)

    def __pow__(self, exponent):
        r"""
        Return the power of this point on the unit circle
        to the given ``exponent``.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup

            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=0.1)^2
            e^(2*pi*0.200000000000000)
            sage: C(exponent=0.1)^QQ(2/1)
            e^(2*pi*0.200000000000000)

            sage: U = RootsOfUnityGroup()
            sage: a = U(exponent=1/7); a
            zeta7
            sage: a^(7/3)
            zeta3
        """
        P = self.parent()
        return P.element_class(P, self.exponent * exponent)

    def __invert__(self):
        r"""
        Return the inverse of this point on the unit circle.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: ~C(exponent=0.4)
            e^(2*pi*0.600000000000000)
            sage: C(1) / C(exponent=0.4)
            e^(2*pi*0.600000000000000)
        """
        P = self.parent()
        return P.element_class(P, -self.exponent)

    def is_one(self):
        r"""
        Return whether this point on the unit circle is `1`.

        EXAMPLES:

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: C = UnitCircleGroup(QQ)
            sage: C(exponent=0).is_one()
            True
            sage: C(exponent=1/2).is_one()
            False
            sage: C(exponent=2/3).is_one()
            False
        """
        return self.exponent == 0

    def is_minus_one(self):
        r"""
        Return whether this point on the unit circle is `-1`.

        EXAMPLES:

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: C = UnitCircleGroup(QQ)
            sage: C(exponent=0).is_minus_one()
            False
            sage: C(exponent=1/2).is_minus_one()
            True
            sage: C(exponent=2/3).is_minus_one()
            False
        """
        from sage.rings.rational_field import QQ
        return self.exponent == QQ(1)/QQ(2)


class UnitCircleGroup(AbstractArgumentGroup):
    r"""
    A group of points on the unit circle. These points are
    represented by `e^{2\pi\mathit{exponent}}`.

    INPUT:

    - ``base`` -- a SageMath parent representing a subset of the reals

    - ``category`` -- a category

    EXAMPLES::

        sage: from sage.groups.roots_of_unity_group import UnitCircleGroup

        sage: R = UnitCircleGroup(RR); R
        Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
        sage: R(exponent=2.42)
        e^(2*pi*0.420000000000000)

        sage: Q = UnitCircleGroup(QQ); Q
        Unit Circle Group with Exponents in Rational Field modulo ZZ
        sage: Q(exponent=6/5)
        e^(2*pi*1/5)
    """

    Element = UnitCirclePoint

    def _repr_(self):
        r"""
        Return a representation string of this unit circle group.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup(RR)  # indirect doctest
            Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
        """
        return 'Unit Circle Group with Exponents in {} modulo ZZ'.format(self.base())

    def _repr_short_(self):
        r"""
        Return a short representation string of this unit circle group.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup(RR)._repr_short_()
            'U_RR'
        """
        from sage.rings.asymptotic.misc import parent_to_repr_short
        s = parent_to_repr_short(self.base())
        if ' ' in s:
            s = '({})'.format(s)
        return 'U_{}'.format(s)

    def _element_constructor_(self, data, exponent=None):
        r"""
        Construct an element out of the given data.

        INPUT:

        - ``data`` -- an object

        - ``exponent`` -- a number (of a subset of the reals) or ``None``

        OUTPUT:

        A :class:`UnitCirclePoint`.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=1/2)
            e^(2*pi*0.500000000000000)

            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0)
            1
            sage: U(exponent=1)
            1
            sage: U(exponent=1/2)
            -1
            sage: U(exponent=1/4)
            I
            sage: U(exponent=1/3)
            zeta3

            sage: C.<z> = CyclotomicField(6)
            sage: z, U(z)
            (z, zeta6)
            sage: z^2, U(z^2)
            (z - 1, zeta3)

            sage: U(ZZ(-1))
            -1
            sage: U(QQ(-1))
            -1
            sage: U(int(-1))
            -1
        """
        from sage.groups.generic import discrete_log
        from sage.rings.asymptotic.misc import combine_exceptions
        from sage.rings.rational_field import QQ
        from sage.rings.number_field.number_field import NumberField_cyclotomic

        if exponent is None:
            if isinstance(data, int) and data == 0:
                raise ValueError('no input specified')

            elif isinstance(data, self.element_class):
                if data.parent() == self:
                    return data
                exponent = data.exponent

            elif data == 1 or data == '1':
                exponent = 0

            elif data == -1 or data == '-1':
                exponent = QQ(1)/QQ(2)

            else:
                try:
                    P = data.parent()
                except AttributeError:
                    raise TypeError('{} is not in {}'.format(data, self))

                if isinstance(P, NumberField_cyclotomic):
                    zeta = P.gen()
                    n = zeta.multiplicative_order()
                    try:
                        exponent = QQ(discrete_log(data, zeta)) / QQ(n)
                    except ValueError as e:
                        raise combine_exceptions(
                            ValueError('{} is not in {}'.format(data, self)), e)

            if exponent is None:
                raise ValueError('{} is not in {}'.format(data, self))

        elif not isinstance(data, int) or data != 0:
            raise ValueError('input is ambigous: '
                             '{} as well as exponent={} '
                             'specified'.format(data, exponent))

        return self.element_class(self, exponent)


class RootOfUnity(UnitCirclePoint):
    r"""
    A root of unity (i.e. an element of :class:`RootsOfUnityGroup`)
    which is `e^{2\pi\mathit{exponent}}` for a rational ``exponent``.
    """

    def exponent_numerator(self):
        r"""
        Return the numerator of the rational quotient in `[0,1)`
        representing the exponent of this root of unity.

        EXAMPLES::

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: a = U(exponent=2/3); a
            zeta3^2
            sage: a.exponent_numerator()
            2
        """
        return self.exponent.numerator()

    def exponent_denominator(self):
        r"""
        Return the denominator of the rational quotient in `[0,1)`
        representing the exponent of this root of unity.

        EXAMPLES::

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: a = U(exponent=2/3); a
            zeta3^2
            sage: a.exponent_denominator()
            3
        """
        return self.exponent.denominator()

    def _repr_(self):
        r"""
        Return a representation string of this root of unity.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0)
            1
            sage: U(exponent=1/2)
            -1
            sage: U(exponent=1/4)
            I
            sage: U(exponent=3/4)
            -I
            sage: U(exponent=1/3)
            zeta3
            sage: U(exponent=2/3)
            zeta3^2
        """
        from sage.rings.rational_field import QQ
        if self.exponent == 0:
            return '1'
        if self.exponent == QQ(1)/QQ(2):
            return '-1'
        if self.exponent == QQ(1)/QQ(4):
            return 'I'
        if self.exponent == QQ(3)/QQ(4):
            return '-I'
        num = self.exponent_numerator()
        den = self.exponent_denominator()
        zeta = 'zeta{}'.format(den)
        if num == 1:
            return zeta
        return '{}^{}'.format(zeta, num)


class RootsOfUnityGroup(UnitCircleGroup):
    r"""
    The group of all roots of unity.

    INPUT:

    - ``category`` -- a category

    This is a specialized :class:`UnitCircleGroup` with base `\QQ`.

    EXAMPLES::

        sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
        sage: U = RootsOfUnityGroup(); U
        Group of Roots of Unity
        sage: U(exponent=1/4)
        I
    """

    Element = RootOfUnity

    @staticmethod
    def __classcall__(cls, category=None):
        r"""
        See :class:`RootsOfUnityGroup` for more information.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: RootsOfUnityGroup().category()  # indirect doctest
            Category of commutative groups
        """
        category = cls._determine_category_(category)
        return super(AbstractArgumentGroup, cls).__classcall__(
            cls, category)

    def __init__(self, category):
        r"""
        See :class:`RootsOfUnityGroup` for more information.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: RootsOfUnityGroup().base()  # indirect doctest
            Rational Field
        """
        from sage.rings.rational_field import QQ
        return super(RootsOfUnityGroup, self).__init__(base=QQ,
                                                         category=category)
    def _repr_(self):
        r"""
        Return a representation string of this roots of unity group.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: RootsOfUnityGroup()  # indirect doctest
            Group of Roots of Unity
        """
        return 'Group of Roots of Unity'

    def _repr_short_(self):
        r"""
        Return a short representation string of this roots of unity group.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: RootsOfUnityGroup()._repr_short_()
            'U'
        """
        return 'U'


class ArgumentByElement(AbstractArgument):
    r"""
    An element of :class:`ArgumentByElementGroup`.

    INPUT:

    - ``parent`` -- a SageMath parent

    - ``element`` -- a nonzero element of the parent's base

    - ``normalize`` -- a boolean (default: ``True``)
    """

    def __init__(self, parent, element, normalize=True):
        r"""
        See :class:`ArgumentByElement` for more information.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(1+2*I)  # indirect doctest
            e^(I*arg(1.00000000000000 + 2.00000000000000*I))
        """
        super(ArgumentByElement, self).__init__(parent, element, normalize=normalize)
        if self._element_ == 0:
            raise ValueError('{} is not allowed'.format(element))

    @staticmethod
    def _normalize_(element):
        r"""
        Normalizes the given element.

        This is the identity for :class:`ArgumentByElement`.

        INPUT:

        - ``element`` -- an element of the parent's base

        OUTPUT:

        An element.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentByElement
            sage: ArgumentByElement._normalize_(3/2)
            3/2
        """
        return element

    def _repr_(self):
        r"""
        Return a representation string of this argument by element.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(2+3*I)  # indirect doctest
            e^(I*arg(2.00000000000000 + 3.00000000000000*I))
        """
        return 'e^(I*arg({}))'.format(self._element_)

    def _mul_(self, other):
        r"""
        Return the product of this argument by element with ``other``.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(I) * C(1 + I)  # indirect doctest
            e^(I*arg(-1.00000000000000 + 1.00000000000000*I))
        """
        P = self.parent()
        return P.element_class(P, self._element_ * other._element_)

    def __pow__(self, exponent):
        r"""
        Return the power of this argument by element
        to the given ``exponent``.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(I)^5  # indirect doctest
            e^(I*arg(1.00000000000000*I))
            sage: C(1+I)^3  # indirect doctest
            e^(I*arg(-2.00000000000000 + 2.00000000000000*I))
        """
        P = self.parent()
        return P.element_class(P, self._element_ ** exponent)

    def __invert__(self):
        r"""
        Return the inverse of this argument by element.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: ~C(I)  # indirect doctest
            e^(I*arg(-1.00000000000000*I))
        """
        P = self.parent()
        return P.element_class(P, ~self._element_)


class ArgumentByElementGroup(AbstractArgumentGroup):
    r"""
    A group of (complex) arguments. The arguments are represented
    by a the formal argument of an element, i.e.,
    by `\mathrm{arg}(\mathit{element})`.

    INPUT:

    - ``base`` -- a SageMath parent representing a subset of the complex plane

    - ``category`` -- a category

    EXAMPLES::

        sage: from sage.groups.roots_of_unity_group import ArgumentByElementGroup
        sage: C = ArgumentByElementGroup(CC); C
        Unit Circle Group with Argument of Elements in
        Complex Field with 53 bits of precision
        sage: C(1 + 2*I)
        e^(I*arg(1.00000000000000 + 2.00000000000000*I))
    """

    Element = ArgumentByElement

    def _repr_(self):
        r"""
        Return a representation string of this argument by element group.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentByElementGroup
            sage: ArgumentByElementGroup(CC)  # indirect doctest
            Unit Circle Group with Argument of Elements in
            Complex Field with 53 bits of precision
        """
        return 'Unit Circle Group with Argument of Elements in {}'.format(self.base())

    def _repr_short_(self):
        r"""
        Return a short representation string of this  argument by element group.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentByElementGroup
            sage: ArgumentByElementGroup(CC)._repr_short_()
            'Arg_CC'
        """
        from sage.rings.asymptotic.misc import parent_to_repr_short
        s = parent_to_repr_short(self.base())
        if ' ' in s:
            s = '({})'.format(s)
        return 'Arg_{}'.format(s)

    def _element_constructor_(self, data):
        r"""
        Construct an element out of the given data.

        INPUT:

        - ``data`` -- an object

        OUTPUT:

        A :class:`ArgumentByElement`.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(1 + 2*I)  # indirect doctest
            e^(I*arg(1.00000000000000 + 2.00000000000000*I))
            sage: C(1)
            e^(I*arg(1.00000000000000))
            sage: C(-1)
            e^(I*arg(-1.00000000000000))
            sage: C(ZZ(-1))
            e^(I*arg(-1.00000000000000))
            sage: C(QQ(-1))
            e^(I*arg(-1.00000000000000))
            sage: C(int(-1))
            e^(I*arg(-1.00000000000000))
            sage: C('-1')
            e^(I*arg(-1.00000000000000))
        """
        if isinstance(data, int) and data == 0:
            raise ValueError('no input specified')

        elif isinstance(data, self.element_class):
            if data.parent() == self:
                return data
            element = data._element_

        elif data == 1 or data == '1':
            element = 1

        elif data == -1 or data == '-1':
            element = -1

        else:
            element = data

        return self.element_class(self, element)


def exactly_one_is_true(iterable):
    r"""
    Return whether exactly one element of ``iterable`` evaluates ``True``.

    INPUT:

    - ``iterable`` -- an iterable object

    OUTPUT:

    A boolean.

    .. NOTE::

        The implementation is suggested by
        `stackoverflow entry <https://stackoverflow.com/a/16801605/1052778>`_.

    EXAMPLES::

        sage: from sage.groups.roots_of_unity_group import exactly_one_is_true
        sage: exactly_one_is_true([])
        False
        sage: exactly_one_is_true([True])
        True
        sage: exactly_one_is_true([False])
        False
        sage: exactly_one_is_true([True, True])
        False
        sage: exactly_one_is_true([False, True])
        True
        sage: exactly_one_is_true([True, False, True])
        False
        sage: exactly_one_is_true([False, True, False])
        True
    """
    it = iter(iterable)
    return any(it) and not any(it)


class ArgumentGroupFactory(UniqueFactory):
    r"""
    A factory for creating argument groups.

    INPUT:

    - ``data`` -- an object

      The factory will analyze ``data`` and interpret it as
      ``specification`` or ``domain``.

    - ``specification`` -- a string

      The following is possible:

      - ``'U'`` give the :class:`RootsOfUnityGroup`,

      - ``'U_P'``, where ``'P'`` is
        a string representing a SageMath parent which is interpreted as
        ``exponents``.

      - ``'Arg_P'``, where ``'P'`` is
        a string representing a SageMath parent which is interpreted as
        ``domain``

    - ``domain`` -- a SageMath parent representing a subset of the complex plane.
      An instance of :class:`ArgumentByElementGroup` will be created with the given
      ``domain``.

    - ``exponents`` -- a SageMath parent representing a subset of the reals.
      An instance of :class`UnitCircleGroup` will be created with the given
      ``exponents``

    Exactly one of ``data``, ``specification``, ``exponents`` has to be provided.

    Further keyword parameters will be carried on to the initialization of
    the group.

    EXAMPLES::

        sage: from sage.groups.roots_of_unity_group import ArgumentGroup

        sage: ArgumentGroup('U')
        Group of Roots of Unity

        sage: ArgumentGroup(ZZ)
        Group of Roots of Unity
        sage: ArgumentGroup(QQ)
        Group of Roots of Unity
        sage: ArgumentGroup('U_QQ')
        Group of Roots of Unity
        sage: ArgumentGroup(AA)
        Group of Roots of Unity

        sage: ArgumentGroup(RR)
        Unit Circle Group with Exponents in
        Real Field with 53 bits of precision modulo ZZ
        sage: ArgumentGroup('Arg_RR')
        Unit Circle Group with Exponents in
        Real Field with 53 bits of precision modulo ZZ
        sage: ArgumentGroup(RIF)
        Unit Circle Group with Exponents in
        Real Interval Field with 53 bits of precision modulo ZZ
        sage: ArgumentGroup(RBF)
        Unit Circle Group with Exponents in
        Real ball field with 53 bits of precision modulo ZZ

        sage: ArgumentGroup(CC)
        Unit Circle Group with Exponents in
        Real Field with 53 bits of precision modulo ZZ
        sage: ArgumentGroup('Arg_CC')
        Unit Circle Group with Exponents in
        Real Field with 53 bits of precision modulo ZZ
        sage: ArgumentGroup(CIF)
        Unit Circle Group with Exponents in
        Real Interval Field with 53 bits of precision modulo ZZ
        sage: ArgumentGroup(CBF)
        Unit Circle Group with Exponents in
        Real ball field with 53 bits of precision modulo ZZ

        sage: ArgumentGroup(CyclotomicField(3))
        Unit Circle Group with Argument of Elements in
        Cyclotomic Field of order 3 and degree 2
    """
    def create_key_and_extra_args(self,
                                  data=None,
                                  specification=None,
                                  domain=None,
                                  exponents=None,
                                  **kwds):
        r"""
        Normalize the input.

        See :class:`ArgumentGroupFactory` for a description and examples.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentGroup

            sage: ArgumentGroup('U') is ArgumentGroup(exponents=QQ)  # indirect doctest
            True
            sage: ArgumentGroup('Arg_RR') is ArgumentGroup(exponents=RR)  # indirect doctest
            True
            sage: ArgumentGroup('Arg_CC') is ArgumentGroup(domain=CC)  # indirect doctest
            True
        """
        from sage.rings.complex_arb import ComplexBallField
        from sage.rings.complex_field import ComplexField_class
        from sage.rings.complex_interval_field import ComplexIntervalField_class
        from sage.rings.integer_ring import ZZ
        from sage.rings.qqbar import AA
        from sage.rings.rational_field import QQ
        from sage.rings.real_arb import RealBallField
        from sage.rings.real_mpfr import RealField_class
        from sage.rings.real_mpfi import RealIntervalField_class

        if not exactly_one_is_true(
                (data is not None,
                 specification is not None,
                 domain is not None,
                 exponents is not None)):
            raise ValueError('input ambigous')

        if data is not None:
            if isinstance(data, str):
                specification = data
            else:
                domain = data

        if specification is not None:
            if specification == 'U':
                return (RootsOfUnityGroup, ()), kwds
            elif specification.startswith('U_'):
                from sage.rings.asymptotic.misc import repr_short_to_parent
                exponents = repr_short_to_parent(specification[2:])
            elif specification.startswith('Arg_') or specification.startswith('arg_'):
                from sage.rings.asymptotic.misc import repr_short_to_parent
                domain = repr_short_to_parent(specification[4:])
            else:
                raise ValueError('unknown specification {}'.format(specification))

        if domain is not None:
            if domain in (ZZ, QQ, AA):
                # we only need +1 and -1
                return (RootsOfUnityGroup, ()), kwds
            elif isinstance(domain, (RealField_class,
                                     RealIntervalField_class,
                                     RealBallField)):
                return (UnitCircleGroup, (domain,)), kwds
            elif isinstance(domain, (ComplexField_class,
                                     ComplexIntervalField_class,
                                     ComplexBallField)):
                return (UnitCircleGroup, (domain._real_field(),)), kwds
            else:
                return (ArgumentByElementGroup, (domain,)), kwds

        elif exponents is not None:
            if exponents in (ZZ, QQ):
                return (RootsOfUnityGroup, ()), kwds
            else:
                return (UnitCircleGroup, (exponents,)), kwds

    def create_object(self, version, key, **kwds):
        r"""
        Create an object from the given arguments.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import ArgumentGroup
            sage: ArgumentGroup('U')  # indirect doctest
            Group of Roots of Unity
        """
        cls, args = key
        return cls(*args, **kwds)


ArgumentGroup = ArgumentGroupFactory('sage.groups.roots_of_unity_group.ArgumentGroup')
r"""
A factory for argument groups.

This is an instance of :class:`ArgumentGroupFactory` whose documentation
provides more details.
"""
