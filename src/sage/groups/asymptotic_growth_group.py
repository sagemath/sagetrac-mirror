r"""
Asymptotic Growth Group

AUTHORS:

- Benjamin Hackl (2015-01): initial version

"""

#*****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import re

from sage.structure.element import MultiplicativeGroupElement
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR



class GenericGrowthElement(MultiplicativeGroupElement):
    r"""
    Class for a generic asymptotic growth group element. These
    elements hold exactly one asymptotic term and can be compared to
    each other, multiplied and divided, but possess no explicit
    coefficient. At this stage, just the order of magnitude shall be
    managed. In this class, only base structure is handled.
    For a concrete realization see :class:`GrowthElementPower`.

    INPUT:

    - ``parent`` -- a GenericGrowthGroup, the parent of the
      element.

    OUTPUT:

    A generic element of a generic growth group.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.GenericGrowthGroup()
        sage: e = agg.GenericGrowthElement(P); e
        Generic element of a Generic Asymptotic Growth Group
        sage: e.parent()
        Generic Asymptotic Growth Group
    """

    def __init__(self, parent, raw_element):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GenericGrowthGroup()
            sage: e = agg.GenericGrowthElement(P)
            sage: e.category()
            Category of elements of Generic Asymptotic Growth Group

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: e = agg.GenericGrowthElement(None)
            Traceback (most recent call last):
            ...
            ValueError: The parent must be provided
        """
        if parent is None:
            raise ValueError("The parent must be provided")
        super(GenericGrowthElement, self).__init__(parent=parent)

        self._raw_element_ = parent.base()(raw_element)


    def _mul_(self, other):
        r"""
        Abstract multiplication method for generic asymptotic growth
        group elements. For a concrete realization see
        :meth:`GrowthElementPower._mul_`.

        INPUT:

        - ``other`` -- a GenericGrowthElement from the same parent as
          ``self``.

        OUTPUT:

        A GenericGrowthElement representing the product of ``self``
        and ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e2 = P(x=None, exponent=3)
            sage: e1._mul_(e2)
            x^5
            sage: e1 * e2 * e1
            x^7
        """
        raise NotImplementedError("Only implemented in concrete realizations")


    def _repr_(self):
        r"""
        Represent the abstract element of an abstract asymptotic growth
        group as "Generic element of a Generic Asymptotic Growth
        Group".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GenericGrowthGroup()
            sage: e = agg.GenericGrowthElement(P); e._repr_()
            'Generic element of a Generic Asymptotic Growth Group'
        """
        return 'GenericGrowthElement(%s)' % (self._raw_element_,)


    def is_le_one(self):
        r"""
        Abstract method for comparison with one.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GrowthGroupPower('x')
            sage: (~G.gen()).is_le_one()
            True
        """
        raise NotImplementedError("Only implemented in concrete realizations")


    def __le__(self, other):
        r"""
        Return if this growth element is at most (less than or equal
        to) ``other``.

        INPUT:

        - ``other`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

            The comparison of two elements with the same parent is done in
            :meth:`_le_`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G.an_element() <= G.an_element()
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations

        ::

            sage: P_ZZ = agg.GrowthGroupPower('x', base=ZZ)
            sage: P_QQ = agg.GrowthGroupPower('x', base=QQ)
            sage: P_ZZ.gen() <= P_QQ.gen()^2
            True
            sage: ~P_ZZ.gen() <= P_ZZ.gen()
            True
        """
        from sage.structure.element import have_same_parent
        if have_same_parent(self, other):
            return self._le_(other)

        from sage.structure.element import get_coercion_model
        import operator
        try:
            return get_coercion_model().bin_op(self, other, operator.le)
        except TypeError:
            return False


    def _le_(self, other):
        r"""
        Return if this :class:`GenericGrowthElement` is at most (less
        than or equal to) ``other``.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`GenericGrowthElement`.

        TESTS::

            sage: TODO
        """
        raise NotImplementedError("Only implemented in concrete realizations")


    def __eq__(self, other):
        r"""
        Return if this growth element is equal to ``other``.

        INPUT:

        - ``other`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

            The comparison of two elements with the same parent is done in
            :meth:`_eq_`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G.an_element() == G.an_element()
            True
            sage: G(raw_element=42) == G(raw_element=7)
            False

        ::

            sage: G_ZZ = agg.GenericGrowthGroup(ZZ)
            sage: G_QQ = agg.GenericGrowthGroup(QQ)
            sage: G_ZZ(raw_element=1) == G_QQ(raw_element=1)
            True

        ::

            sage: P_ZZ = agg.GrowthGroupPower('x', base=ZZ)
            sage: P_QQ = agg.GrowthGroupPower('x', base=QQ)
            sage: P_ZZ.gen() == P_QQ.gen()
            True
            sage: ~P_ZZ.gen() == P_ZZ.gen()
            False
            sage: ~P_ZZ(1) == P_ZZ(1)
            True
        """
        from sage.structure.element import have_same_parent
        if have_same_parent(self, other):
            return self._eq_(other)

        from sage.structure.element import get_coercion_model
        import operator
        try:
            return get_coercion_model().bin_op(self, other, operator.eq)
        except TypeError:
            return False


    def _eq_(self, other):
        r"""
        Return if this :class:`GenericGrowthElement` is equal to ``other``.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`GenericGrowthElement`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower('x')
            sage: e1 = P(raw_element=1)
            sage: e1._eq_(P.gen())
            True
            sage: e2 = e1^4
            sage: e2 == e1^2 * e1 * e1
            True
            sage: e2 == e1
            False
        """
        return self._raw_element_ == other._raw_element_

    def __hash__(self):
        r"""
        Return the hash of this element

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        The hash uses the representation string of this element
        produced by :meth:`_repr_`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ);
            sage: hash(G(raw_element=42))  # random
            5656565656565656
        """
        return hash(repr(self))


class GenericGrowthGroup(Parent, UniqueRepresentation):
    r"""
    Class for the generic asymptotic growth group. Only base
    structure is handled here, for a concrete realization see
    :class:`GrowthGroupPower`.

    INPUT:

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a
      subcategory of ``Join of Category of groups and Category
      of posets``. This is also the default category if ``None``
      is specified.

    - ``base`` -- The base of the parent associated to concrete
      implementations of this abstract base class.

    OUTPUT:

    A generic asymptotic growth group.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.GenericGrowthGroup(); P
        Generic Asymptotic Growth Group
    """

    # TODO: implement some sort of "assume", where basic assumptions
    # for the variables can be stored. --> within the cartesian product

    # TODO: implement a cartesian product class for the asymptotic
    # growth group. the "standard" cart. prod. produces an element
    # of the cart. prod. class, which does not have a order structure.

    # enable the category framework for elements
    Element = GenericGrowthElement


    def __init__(self, base, category=None):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GenericGrowthGroup().category()
            Join of Category of groups and Category of posets

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GenericGrowthGroup()
            sage: P.is_parent_of(agg.GenericGrowthElement(P))
            True
            sage: P2 = agg.GenericGrowthGroup(category=FiniteGroups()\
            ....:                               & Posets())
            sage: P2.category()
            Join of Category of finite groups and Category of finite posets
            sage: P3 = agg.GenericGrowthGroup(category=Rings())
            Traceback (most recent call last):
            ...
            ValueError: (Category of rings,) is not a subcategory of Join of Category of groups and Category of posets
        """
        from sage.categories.groups import Groups
        from sage.categories.posets import Posets

        if category is None:
            category = Groups() & Posets()
        else:
            if not isinstance(category, tuple):
                category = (category, )
            if not any(cat.is_subcategory(Groups() & Posets()) for cat in
                       category):
                raise ValueError("%s is not a subcategory of %s"
                                 % (category, Groups() & Posets()))
        super(GenericGrowthGroup, self).__init__(category=category,
                                                 base=base)


    def le(self, left, right):
        r"""
        Return if the growth element ``left`` is at most (less than or
        equal to) the growth element ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        Return whether the asymptotic order of magnitude of `x` is less
        than or equal to the asymptotic order of magnitude of `y`.

        INPUT:

        - ``x``, ``y`` -- elements of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GrowthGroupPower('x')
            sage: x = G.gen()
            sage: G.le(x, x^2)
            True
            sage: G.le(x^2, x)
            False
            sage: G.le(x^0, 1)
            True
        """
        return (self(left) / self(right)).is_le_one()


    def _repr_(self):
        r"""
        A representation string for this generic growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ)  # indirect doctest
            Generic Growth Group over Integer Ring
        """
        return "Generic Growth Group over %s" % (self.base(),)


    def __hash__(self):
        r"""
        Return the hash of this group.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: hash(agg.GenericGrowthGroup(ZZ))  # random
            4242424242424242
        """
        return hash((self.__class__, self.base()))


    def _an_element_(self):
        r"""
        Return an element of ``self``.

        INPUT:

        Nothing.

        OUTPUT:

        An element of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GenericGrowthGroup();
            sage: P.an_element()  # indirect doctest
            Generic element of a Generic Asymptotic Growth Group
        """
        return self.element_class(self, self.base().an_element())


    def _element_constructor_(self, data, raw_element=None):
        r"""
        Converts given object to this growth group.

        INPUT:

        - ``data`` -- an object representing the element to be
          initialized.

        - ``raw_element`` -- (default: ``None``) if given, then this is
          directly passed to the element constructor (i.e., no conversion
          is performed).

        OUTPUT:

        An element of an asymptotic power growth group.

        .. NOTE::

            This method calls :meth:`_convert_`, which does the actual
            conversion from ``data``.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G_ZZ = agg.GenericGrowthGroup(ZZ)
            sage: z = G_ZZ(raw_element=42); z
            GenericGrowthElement(42)
            sage: z is G_ZZ(z)
            True

        ::

            sage: G_QQ = agg.GenericGrowthGroup(QQ)
            sage: q = G_QQ(raw_element=42)
            sage: q is z
            False
            sage: G_ZZ(q)
            GenericGrowthElement(42)
            sage: G_QQ(z)
            GenericGrowthElement(42)
            sage: q is G_ZZ(q)
            False

        ::

            sage: G_ZZ()
            Traceback (most recent call last):
            ...
            ValueError: No input specified. Cannot continue.
            sage: G_ZZ('blub')
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert blub.
            sage: G_ZZ('x', raw_element=42)
            Traceback (most recent call last):
            ...
            ValueError: Input is ambigous: x as well as raw_element=42 are specified
        """
        if raw_element is None:
            if type(data) == self.element_class and data.parent() == self:
                return data
            elif isinstance(data, GenericGrowthElement):
                raw_element = data._raw_element_
            elif type(data) == int and data == 0:
                raise ValueError('No input specified. Cannot continue.')
            else:
                raw_element = self._convert_(data)
            if raw_element is None:
                raise ValueError('Cannot convert %s.' % (data,))
        elif type(data) != int or data != 0:
            raise ValueError('Input is ambigous: '
                             '%s as well as raw_element=%s '
                             'are specified' % (data, raw_element))

        return self.element_class(self, raw_element)


    def _convert_(self, data):
        r"""
        Converts given ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        .. NOTE::

            This method always returns ``None`` in this abstract base
            class and should be overridden in inherited class.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G._convert_('icecream') is None
            True
        """
        pass


    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this growth group.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            Another growth group ``S`` coerces into this growth group
            if and only if the base of ``S`` coerces into the base of
            this growth group.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G1 = agg.GrowthGroupPower('x', base=ZZ)
            sage: G2 = agg.GrowthGroupPower('x', base=QQ)
            sage: bool(G1._coerce_map_from_(G2))
            False
            sage: bool(G2._coerce_map_from_(G1))
            True
        """
        if isinstance(S, GenericGrowthGroup):
            if self.base().has_coerce_map_from(S.base()):
                return True


class GrowthElementPower(GenericGrowthElement):
    r"""
    Class for the concrete realization of asymptotic growth group
    elements in the case of polynomial growth. These elements hold
    exactly one asymptotic term.

    A power growth element represents a polynomial term
    `\operatorname{variable}^{\operatorname{exponent}}`.
    More complex constructions including logarithmic or exponential
    terms can be constructed via a cartesian product of the related
    growth groups. Asymptotic growth elements can be multiplied,
    divided, inverted, and compared to each other. However, they
    possess no explicit coefficient.

    The elements can be specified by either an expression ``x`` being
    a string, an element from the symbolic or a polynomial ring or the
    integer `1`. On the other hand, elements can also be specified
    by their exponent.

    INPUT:

    - ``parent`` -- a :class:`GrowthGroupPower`, the
      parent of the element.
    - ``x`` -- an expression (string, polynomial ring element,
      symbolic ring element, or the integer `1`) representing
      the element to be initialized.
    - ``exponent`` -- the exponent of the power element.

    OUTPUT:

    An asymptotic power growth element with the specified
    parent and magnitude of growth, i.e. exponent.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.GrowthGroupPower("x")
        sage: e1 = P(x=1); e1
        1
        sage: e2 = P(x=None, exponent=2); e2
        x^2
        sage: e1 == e2
        False
        sage: P.le(e1, e2)
        True
        sage: P.le(e1, P.gen()) and P.le(P.gen(), e2)
        True
    """
    # TODO: implement comparison for the cartesian product structure.


#    def __init__(self, parent, exponent):
#        r"""
#        See :class:`GrowthElementPower` for more information.
#
#        EXAMPLES::
#
#            sage: import sage.groups.asymptotic_growth_group as agg
#            sage: G = agg.GrowthGroupPower('x')
#            sage: e1 = G.gen(); e1
#            x
#            sage: e1.is_idempotent()
#            False
#            sage: e1.is_one()
#            False
#            sage: e1.parent()
#            Asymptotic Power Growth Group in x over Integer Ring
#            sage: e2 = G.one(); e2
#            1
#            sage: e2.is_idempotent() and e2.is_one()
#            True
#        """
#        super(GrowthElementPower, self).__init__(parent=parent, raw_element=exponent)


    @property
    def exponent(self):
        return self._raw_element_


    def _repr_(self):
        r"""
        Represent the asymptotic power growth element as
        ``variable^{exponent}``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x", base=QQ)
            sage: P(x=1)._repr_()
            '1'
            sage: P(x=None, exponent=5)._repr_()
            'x^5'
            sage: P(x=None, exponent=1/2)._repr_()
            'x^(1/2)'
        """
        if self.exponent == 0:
            return "1"
        elif self.exponent == 1:
            return self.parent()._var_
        elif self.exponent in ZZ:
            return self.parent()._var_ + "^" + str(self.exponent)
        else:
            return self.parent()._var_ + "^(" + str(self.exponent) + ")"


    def _mul_(self, other):
        r"""
        Multiply two asymptotic power growth elements from the
        same parent by adding their exponents.

        INPUT:

        - ``other`` -- the asymptotic growth element to be
          multiplied with ``self``.

        OUTPUT:

        An asymptotic power growth element representing the product
        of ``self`` and ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e2 = P(x=None, exponent=3)
            sage: e3 = e1._mul_(e2); e3
            x^5
            sage: e3 == e1 * e2
            True
        """
        cls = self.__class__
        return cls(self.parent(), self.exponent + other.exponent)


    def __invert__(self):
        r"""
        Return the multiplicative inverse from a given asymptotic power
        growth element.

        INPUT:

        Nothing.

        OUTPUT:

        The multiplicative inverse asymptotic power growth element
        of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e2 = e1.__invert__(); e2
            x^-2
            sage: e2 == ~e1
            True
        """
        cls = self.__class__
        return cls(self.parent(), -self.exponent)


    def _div_(self, other):
        r"""
        Divide two asymptotic power growth elements from the same
        parent by subtracting their exponents.

        INPUT:

        - ``other`` -- the asymptotic growth element which ``self``
          is divided by.

        OUTPUT:

        The result of the division of ``self`` by ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e1._div_(P.gen())
            x
            sage: e1._div_(P.gen()) == e1 / P.gen()
            True
        """
        cls = self.__class__
        return cls(self.parent(), self.exponent - other.exponent)


    def __pow__(self, power):
        r"""
        Return a asymptotic power element to the power of
        ``power``.

        INPUT:

        - ``power`` -- a rational number.

        OUTPUT:

        The asymptotic power growth element ``self`` to the power of
        ``power``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: P.gen().__pow__(5)
            x^5
            sage: P.gen().__pow__(1/2)
            x^(1/2)
            sage: P.gen()^7
            x^7
        """
        new_exponent = self.exponent * power
        if new_exponent in self.parent().base():
            return self.parent()(raw_element=self.exponent * power)

        if new_exponent in RR:
            pnt = GrowthGroupPower(self.parent()._var_,
                                   base=new_exponent.parent())
            return pnt(raw_element=new_exponent)
        else:
            raise NotImplementedError("Only real exponents are implemented.")


    def __eq__(self, other):
        r"""
        Tests for equality of the given elements (with taking care of
        different parents by using the coercion model).

        INPUT:

        - ``other`` -- power element to be compared with ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P1 = agg.GrowthGroupPower("x", base=ZZ)
            sage: P2 = agg.GrowthGroupPower("x", base=QQ)
            sage: P1.gen() == P2.gen()
            True
        """
        from sage.structure.element import have_same_parent
        if have_same_parent(self, other):
            return self._eq_(other)

        from sage.structure.element import get_coercion_model
        import operator
        try:
            return get_coercion_model().bin_op(self, other, operator.eq)
        except TypeError:
            return False


    def _eq_(self, other):
        r"""
        Return whether the asymptotic power growth elements
        ``self`` and ``other`` are equal and have the same parent.

        INPUT:

        - ``other`` -- an asymptotic power growth element to
          be compared to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P(x=None, exponent=1)
            sage: e1._eq_(P.gen())
            True
            sage: e2 = e1^4
            sage: e2 == e1^2 * e1 * e1
            True
            sage: e2 == e1
            False
        """
        return self.exponent == other.exponent


    def _le_(self, other):
        r"""
        Return whether the exponent of ``self`` is less than or equal
        to the exponent of ``other`` by calling
        :meth:`GrowthGroupPower.le`.

        INPUT:

        - ``other`` -- a growth power element from the same parent to
          be compared to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P.gen()
            sage: e2 = P(None, exponent=2)
            sage: e1._le_(e2)
            True
            sage: e2._le_(e1)
            False
        """
        return self.parent().le(self, other)


    def is_le_one(self):
        r"""
        Return whether or not the growth of the asymptotic power
        growth element ``self`` is less than or equal to the
        (constant) growth of `1`.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P.gen()
            sage: e1.is_le_one()
            False
            sage: (P.one() / P.gen()).is_le_one()
            True
        """
        return self.exponent <= 0

    def __hash__(self):
        r"""
        Return the hash of the tuple containing the exponent and the
        parent.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: hash(P.gen())  # random
            -3234005094684624010
        """
        return hash((self.exponent, self.parent()))


class GrowthGroupPower(GenericGrowthGroup):
    r"""
    Class for the concrete realization of asymptotic power growth
    groups. These are the parents for the
    :class:`GrowthElementPower` elements. These parents
    contain single monomial terms in a specified variable with
    exponents from a specified base.

    More complex growth elements can be constructed by constructing the
    cartesian product of various asymptotic growth groups.

    INPUT:

    - ``var`` -- either a symbol from the symbolic ring, a
      generator from a polynomial ring, or a alphanumeric string
      starting with a letter and optionally containing underscores.

    - ``base`` -- the base structure containing the exponents.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a
      subcategory of ``Join of Category of groups and Category of
      posets``. This is also the default category if ``None`` is
      specified.

    OUTPUT:

    An asymptotic power growth group.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.GrowthGroupPower("x"); P
        Asymptotic Power Growth Group in x over Integer Ring
    """
    # TODO: implement the cartesian product structure

    # enable the category framework for elements
    Element = GrowthElementPower


    @staticmethod
    def __classcall__(cls, var=None, base=ZZ, category=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: x1, x2, x3 = "x", PolynomialRing(ZZ, "x").gen(), var("x")
            sage: P1 = agg.GrowthGroupPower(x1)
            sage: P2 = agg.GrowthGroupPower(x2)
            sage: P3 = agg.GrowthGroupPower(x3)
            sage: P1 is P2 and P2 is P3
            True
            sage: x4 = buffer("xylophone", 0, 1)
            sage: P4 = agg.GrowthGroupPower(x4)
            sage: P1 is P4
            True
        """
        if hasattr(var, "is_symbol") and var.is_symbol():
            var = repr(var)
        elif hasattr(var, "is_gen") and var.is_gen():
            var = repr(var)
        elif isinstance(var, buffer):
            var = str(var)
        return super(GrowthGroupPower, cls).\
            __classcall__(cls, var, base, category)


    def __init__(self, var, base=ZZ, category=None):
        r"""
        For more information see :class:`GrowthGroupPower`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P1 = agg.GrowthGroupPower("x"); P1
            Asymptotic Power Growth Group in x over Integer Ring
            sage: var('n')
            n
            sage: P2 = agg.GrowthGroupPower(n, base=QQ); P2
            Asymptotic Power Growth Group in n over Rational Field
            sage: y = PolynomialRing(ZZ, "y").gen()
            sage: P3 = agg.GrowthGroupPower(y); P3
            Asymptotic Power Growth Group in y over Integer Ring
            sage: P4 = agg.GrowthGroupPower("y"); P4
            Asymptotic Power Growth Group in y over Integer Ring
            sage: P3 is P4
            True
        """
        if var is None:
            raise ValueError("Variable for initialization required.")
        else:
            if re.match("^[A-Za-z][A-Za-z0-9_]*$", var):
                self._var_ = str(var)
            else:
                raise ValueError("Only alphanumeric strings starting with a "
                                 "letter may be variables")
        super(GrowthGroupPower, self).\
            __init__(category=category, base=base)
        self._populate_coercion_lists_()


    def _convert_(self, data):
        r"""

        Converts given ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        .. NOTE::

            TODO (mention possibilites here?)

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower('x', base=ZZ)
            sage: P._convert_('icecream') is None
            True
            sage: P(1)  # indirect doctest
            1
            sage: P('x')  # indirect doctest
            x

        ::

            sage: P(x)  # indirect doctest
            x
            sage: P(x^-333)  # indirect doctest
            x^-333
            sage: P(log(x)^2)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert log(x)^2.
        """
        if data == 1:
            return self.base().zero()
        if str(data) == self._var_:
            return self.base().one()  # not sure if this will work for
                                      # all possible bases

        try:
            P = data.parent()
        except AttributeError:
            return  # this has to end here

        from sage.symbolic.ring import SR
        import operator
        if P is SR:
            if data.operator() == operator.pow:
                base, exponent = data.operands()
                if str(base) == self._var_:
                    return exponent
        #elif ...
        #TODO: SR; PolynomialRing

        ###TODO
        #if exponent is None and x.parent() is parent:
        #    self.exponent = x.exponent
        #    super(GrowthElementPower, self).__init__(parent=parent)
        #if exponent not in RR:
        #    raise NotImplementedError("Non-real exponents are not supported.")
        #else:
        #    self.exponent = parent.base()(exponent)


    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this growth group.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P_x_ZZ = agg.GrowthGroupPower('x', base=ZZ)
            sage: P_x_QQ = agg.GrowthGroupPower('x', base=QQ)
            sage: bool(P_x_ZZ._coerce_map_from_(P_x_QQ))
            False
            sage: bool(P_x_QQ._coerce_map_from_(P_x_ZZ))
            True
            sage: P_y_ZZ = agg.GrowthGroupPower('y', base=ZZ)
            sage: bool(P_y_ZZ._coerce_map_from_(P_x_ZZ))
            False
            sage: bool(P_x_ZZ._coerce_map_from_(P_y_ZZ))
            False
            sage: bool(P_y_ZZ._coerce_map_from_(P_x_QQ))
            False
            sage: bool(P_x_QQ._coerce_map_from_(P_y_ZZ))
            False
        """
        if super(GrowthGroupPower, self)._coerce_map_from_(S):
            if self._var_ == S._var_:
                return True
        # TODO: SR, PolynomialRing



    def _repr_(self):
        r"""
        Represent the asymptotic power growth group as
        "Asymptotic Power Growth Group in ``variable``
        over ``base``".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GrowthGroupPower("x")._repr_()
            'Asymptotic Power Growth Group in x over Integer Ring'
            sage: agg.GrowthGroupPower("v_107", base=QQ)._repr_()
            'Asymptotic Power Growth Group in v_107 over Rational Field'
        """
        return "Asymptotic Power Growth Group in %s over %s" \
               % (self._var_, self.base())


    def _an_element_(self):
        r"""
        Return an element of ``self``.

        INPUT:

        Nothing.

        OUTPUT:

        An element of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: P.an_element()  # indirect doctest
            x
        """
        return self.gen()

    def gens(self):
        r"""
        Return a list with all generators of ``self``. For power growth
        groups, this is a list with only one element: the variable
        to the power `1`.

        INPUT:

        Nothing.

        OUTPUT:

        A list of generators.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: P.gens()
            (x,)
        """
        return (self.gen(), )


    def ngens(self):
        r"""
        Return the number of generators of ``self``. For power growth
        groups, this is exactly `1`.

        INPUT:

        Nothing.

        OUTPUT:

        The number of generators of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: P.ngens()
            1
        """
        return 1


    def one(self):
        r"""
        Return the neutral element of the asymptotic power growth
        group.

        INPUT:

        Nothing.

        OUTPUT:

        The neutral element of the asymptotic power growth group
        ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: e1 = agg.GrowthGroupPower("x").one(); e1
            1
            sage: e1.is_idempotent()
            True
        """
        return self(1)


    def gen(self):
        r"""
        Return the asymptotic power growth element with
        exponent 1.

        INPUT:

        Nothing.

        OUTPUT:

        The asymptotic power growth element with exponent 1.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: e1 = agg.GrowthGroupPower("x").gen(); e1
            x
            sage: e1.exponent == 1
            True
        """
        return self(raw_element=self.base().gen())


    def __hash__(self):
        r"""
        Return the hash of the tuple containing the variable and the
        base.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: hash(P)  # random
            -8144479309627091876
        """
        return hash((self._var_, self.base()))
