r"""
Asymptotic Term Monoid

This module implements asymptotic term monoids. This structures build
upon the (asymptotic) growth groups. While growth elements only model
the growth of a function as it tends towards infinity, an asymptotic
term is basically a wrapper for the growth which specifies the "type"
of the term. The multiplication is the associated group operation.

Besides an abstract base term, :class:`GenericTerm`, this module
implements the following types of terms:

- :class:`OTerm` -- O terms in infinity, see
  :wikipedia:`Big_O_notation`.
- :class:`TermWithCoefficient` -- abstract base class for
  asymptotic terms with coefficients.
- :class:`ExactTerm` -- this class represents a growth element
  multiplied with some non-zero coefficient from a base ring.
- :class:`LTermGeneric` -- abstract base class for "L terms", i.e.
  "concrete" O terms.

A characteristic property of asymptotic terms is that some terms are
able to "absorb" other terms (see
:meth:`sage.monoids.AsymptoticTermMonoid.GenericTerm.absorb`). For
instance, `O(x^2)` is able to absorb `O(x)` (with result
`O(x^2)`), and `3*x^5` is able to absorb `-2*x^5` (with result
`1*x^5`). Essentially, absorption can be interpreted as the
addition of "compatible" terms.

.. TODO::

    - The behavior of L terms depends on their underlying growth
      group. These concrete realizations of L terms have to be
      implemented.

AUTHORS:

- Benjamin Hackl (2015-01): initial version
- Benjamin Hackl (2015-06): refactoring caused by refactoring growth groups

"""

# *****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************

import sage


class GenericTerm(sage.structure.element.MonoidElement):
    r"""
    Base class for asymptotic terms. Mainly the structure and
    several properties asymptotic terms have are handled here.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an element of the parent's ``growth_group``,
      specifying the growth of the asymptotic term.

    OUTPUT:

    A generic asymptotic term.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
        sage: P = atm.GenericTermMonoid(MG)
        sage: t1, t2 = P(x), P(x^2); (t1, t2)
        (Generic Term with growth x, Generic Term with growth x^2)
        sage: t1 * t2
        Generic Term with growth x^3
        sage: t1.can_absorb(t2)
        False
        sage: t1.absorb(t2)
        Traceback (most recent call last):
        ...
        ArithmeticError: Generic Term with growth x cannot absorb Generic Term with growth x^2
        sage: t1.can_absorb(t1)
        False
    """

    def __init__(self, parent, growth):
        r"""
        See :class:`GenericTerm` for more information.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: P = atm.GenericTermMonoid(MG)
            sage: P(x^2)
            Generic Term with growth x^2
        """
        from sage.groups.asymptotic_growth_group import GenericGrowthElement

        if parent is None:
            raise ValueError('The parent must be provided')
        if growth is None or not isinstance(growth, GenericGrowthElement):
            raise ValueError('The growth must be provided and must inherit '
                             'from GenericGrowthElement')
        else:
            if growth not in parent.growth_group():
                raise ValueError("%s is not in the parent's "
                                 "specified growth group" % growth)
            else:
                self.growth = growth
        super(GenericTerm, self).__init__(parent=parent)


    def _mul_(self, other):
        r"""
        Abstract multiplication method for generic terms.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A :class:`GenericTerm` representing the product of ``self``
        and ``other``.

        .. NOTE::

            This method os called by the coercion framework, thus,
            it can be assumed that this element, as well as ``other``
            are from a common parent.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: P = atm.GenericTermMonoid(MG)
            sage: t1, t2 = P(x), P(x^2)
            sage: t1, t2
            (Generic Term with growth x, Generic Term with growth x^2)
            sage: t1 * t2
            Generic Term with growth x^3
        """
        return self.parent()(self.growth * other.growth)


    def can_absorb(self, other):
        r"""
        Check, whether this asymptotic term is able to absorb
        the asymptotic term ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            For example, an :class:`OTerm`
            is able to *absorb* another :class:`OTerm` or an
            :class:`ExactTerm` with weaker (or equal) growth. For more
            information see :class:`OTerm`, :class:`ExactTerm`, and
            :class:`LTermGeneric`.

        EXAMPLES:

        We want to show step by step which terms can be absorbed
        by which other terms. We start by defining the necessary
        term monoids and some terms::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: OT = atm.OTermMonoid(growth_group=MG)
            sage: ET = atm.ExactTermMonoid(growth_group=MG, base_ring=QQ)
            sage: LT = atm.LTermGenericMonoid(growth_group=MG, base_ring=QQ)
            sage: ot1, ot2 = OT(x), OT(x^2)
            sage: et1 = ET(x^2, 2)
            sage: lt1 = LT(x^2, 2, start=0)

        :class:`OTerm` is able to absorb other :class:`OTerm`,
        :class:`LTermGeneric` (and descendants thereof) as well
        as :class:`ExactTerm`, as long as the growth of the other
        term is less than or equal to the growth of this element::

            sage: ot1, ot2
            (O(x), O(x^2))
            sage: ot1.can_absorb(ot2), ot2.can_absorb(ot1)
            (False, True)
            sage: et1
            2 * x^2
            sage: ot1.can_absorb(et1)
            False
            sage: ot2.can_absorb(et1)
            True
            sage: lt1
            2 * L(x^2, 0)
            sage: ot1.can_absorb(lt1)
            False
            sage: ot2.can_absorb(lt1)
            True

        :class:`ExactTerm` can only absorb another
        :class:`ExactTerm` if the growth coincides with the
        growth of this element::

            sage: et1.can_absorb(ET(x^2, 5))
            True
            sage: any(et1.can_absorb(t) for t in [ot1, ot2, lt1])
            False

        :class:`LTermGeneric` can absorb arbitrary other
        :class:`LTermGeneric`, and :class:`ExactTerm` whose
        growth is less than or equal to the growth of this element::

            sage: any(lt1.can_absorb(t) for t in [ot1, ot2])
            False
            sage: lt1.can_absorb(LT(x^5, 1, start=0))
            True
            sage: lt1.can_absorb(et1)
            True
        """
        if not isinstance(other, GenericTerm):
            raise TypeError('%s is not an asymptotic term' % other)
        if isinstance(self, OTerm):
            if isinstance(other, (OTerm, LTermGeneric, ExactTerm)):
                return other <= self
            else:
                return False
        elif isinstance(self, LTermGeneric):
            if isinstance(other, ExactTerm):
                return other <= self
            elif isinstance(other, LTermGeneric):
                return True
            else:
                return False
        elif isinstance(self, ExactTerm):
            if isinstance(other, ExactTerm):
                return self.growth == other.growth
            else:
                return False
        else:
            return False


    def absorb(self, other):
        r"""
        Absorb the asymptotic term ``other``, yielding a new
        asymptotic term (or ``None``). For a more detailed
        explanation of the *absorption* of asymptotic terms see
        the introduction of this module, or the following examples.

        INPUT:

        - ``other`` -- an asymptotic term that can be absorbed
          by this element.

        OUTPUT:

        An asymptotic term or ``None``.

        EXAMPLES:

        We want to demonstrate in which cases an asymptotic term
        is able to absorb another term, as well as explain the output
        of this operation. We start by defining some parents and
        elements::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG_QQ = agg.MonomialGrowthGroup(QQ, 'x'); x = MG_QQ.gen()
            sage: OT = atm.OTermMonoid(MG_QQ)
            sage: ET = atm.ExactTermMonoid(growth_group=MG_QQ, base_ring=QQ)
            sage: LT = atm.LTermGenericMonoid(growth_group=MG_QQ, base_ring=QQ)
            sage: ot1, ot2 = OT(x), OT(x^2)
            sage: et1, et2 = ET(x, 100), ET(x^2, 2)
            sage: et3, et4 = ET(x^2, 1), ET(x^2, -2)
            sage: lt1 = LT(x, 5, start=0)

        Because of the definition of `O` terms (see
        :wikipedia:`Big_O_notation`), they are able to absorb all
        other asymptotic terms with weaker or equal growth. The
        result of the absorption is the original `O` Term::

            sage: ot1.absorb(ot1)
            O(x)
            sage: ot1.absorb(et1)
            O(x)
            sage: ot1.absorb(lt1)
            O(x)

        This corresponds to `O(x) + O(x) = O(x)`,
        `O(x) + 100x = O(x)`, and `O(x) + 5\cdot L(x,0) = O(x)`.
        If absorb is called on a term that cannot be absorbed, an
        ``ArithmeticError`` is raised::

            sage: ot1.absorb(ot2)
            Traceback (most recent call last):
            ...
            ArithmeticError: O(x) cannot absorb O(x^2)

        This would only work the other way around::

            sage: ot2.absorb(ot1)
            O(x^2)

        :class:`ExactTerm` is able to absorb another
        :class:`ExactTerm` if the terms have the same growth. In this
        case, *absorption* is nothing else than an addition of the
        respective coefficients::

            sage: et2.absorb(et3)
            3 * x^2
            sage: et3.absorb(et2)
            3 * x^2
            sage: et3.absorb(et4)
            -1 * x^2

        Note that, for technical reasons, the coefficient `0` is not
        allowed, and thus ``None`` is returned if two exact terms
        cancel each other out::

            sage: et2.absorb(et4)
            sage: repr(et4.absorb(et2))
            'None'

        .. TODO:

            The absorption of `L` terms is implemented at a later
            point.
        """
        from sage.structure.element import have_same_parent

        if not self.can_absorb(other):
            raise ArithmeticError('%s cannot absorb %s' % (self, other))

        if have_same_parent(self, other):
            return self._absorb_(other)

        from sage.structure.element import get_coercion_model

        try:
            return get_coercion_model().bin_op(self, other,
                                               lambda self, other:
                                               self._absorb_(other))
        except TypeError:
            return False


    def _absorb_(self, other):
        r"""
        Let this element absorb ``other``.

        INPUT:

        - ``other`` -- an asymptotic term from the same parent as
          this element.

        OUTPUT:

        An asymptotic term or ``None``.

        .. NOTE::

            This is not implemented for abstract base classes. For
            concrete realizations see :meth:`OTerm._absorb_`
            or :meth:`ExactTerm._absorb_`.



        EXAMPLES:

        First, we define some asymptotic terms::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: P = atm.GenericTermMonoid(MG)
            sage: t1, t2 = P(x), P(x^2)

        When it comes to absorption, note that the method
        :meth:`can_absorb` (which is called before absorption takes
        place) does not allow the absorption of generic terms. Thus,
        an ``ArithmeticError`` is raised::

            sage: t2.absorb(t1)
            Traceback (most recent call last):
            ...
            ArithmeticError: Generic Term with growth x^2 cannot absorb Generic Term with growth x
        """
        raise NotImplementedError('Not implemented in abstract base classes')


    def __le__(self, other):
        r"""
        Return whether the growth of this term is less than
        or equal to the growth of ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        EXAMPLES:

        First, we define some asymptotic terms (and their parents)::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: GT = atm.GenericTermMonoid(MG)
            sage: OT = atm.OTermMonoid(MG)
            sage: ET_ZZ = atm.ExactTermMonoid(MG, ZZ)
            sage: ET_QQ = atm.ExactTermMonoid(MG, QQ)
            sage: g1, g2 = GT(x), GT(x^2); g1, g2
            (Generic Term with growth x, Generic Term with growth x^2)
            sage: o1, o2 = OT(x^-1), OT(x^3); o1, o2
            (O(1/x), O(x^3))
            sage: t1, t2 = ET_ZZ(x^2, 5), ET_QQ(x^3, 2/7); t1, t2
            (5 * x^2, 2/7 * x^3)

        In order for the comparison to work, the terms have come from
        or coerce into the same parent. Concretely, comparing
        :class:`GenericTerm` to, for example, an :class:`OTerm`
        always yields ``False``::

            sage: g1 <= g2
            True
            sage: o1, g1
            (O(1/x), Generic Term with growth x)
            sage: o1 <= g1
            False

        If the elements of the common parent do not possess
        coefficients, then only the growth is compared::

            sage: o1 <= o1
            True
            sage: o1 <= o2
            True
            sage: o1 <= t1 and t1 <= o2
            True

        Otherwise, for asymptotic exact terms and `L` terms, in case
        of equal growth, also the coefficient is compared::

            sage: t1 <= t2
            True
            sage: t1 <= ET_ZZ(x^2, 3)
            False
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
        Return whether this generic term grows at most (i.e. less
        or equal) like ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method is called by thge coercion framework, thus,
            it can be assumed that this element, as well as ``other``
            are from the same parent.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: GT = atm.GenericTermMonoid(MG)
            sage: t1, t2 = GT(x^-2), GT(x^5); t1, t2
            (Generic Term with growth x^(-2), Generic Term with growth x^5)
            sage: t1._le_(t2)
            True
            sage: t2._le_(t1)
            False
        """
        return self.parent().le(self, other)


    def _repr_(self):
        r"""
        Represent the generic term as ``Generic Term with growth
        ...``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: P = atm.GenericTermMonoid(growth_group=MG)
            sage: P(x)._repr_()
            'Generic Term with growth x'
            sage: P(x^7)._repr_()
            'Generic Term with growth x^7'
        """
        return 'Generic Term with growth ' + repr(self.growth)


class GenericTermMonoid(sage.structure.parent.Parent,
                        sage.structure.unique_representation.UniqueRepresentation):
    r"""
    Parent for generic asymptotic terms. In this class the base
    structure for asymptotic term monoids will be handled. These
    monoids are the parents for asymptotic terms (for example, see
    :class:`GenericTerm` or :class:`OTerm`). Basically, asymptotic
    terms consist of a ``growth`` (which is an asymptotic growth
    group element, for example
    :class:`~sage.groups.asymptotic_growth_group.GrowthElementPower`).

    INPUT:

    - ``growth_group`` -- an asymptotic growth group. This is the parent of
      the elements specifying the growth of the asymptotic terms.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a subcategory
      of ``Join of Category of Monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    OUTPUT:

    A generic asymptotic term monoid.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: MG_x = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG_x.gen()
        sage: MG_y = agg.MonomialGrowthGroup(QQ, 'y'); y = MG_y.gen()
        sage: GT_x_ZZ, GT_y_QQ = atm.GenericTermMonoid(MG_x), atm.GenericTermMonoid(MG_y)
        sage: GT_x_ZZ
        Generic Term Monoid over Monomial Growth Group in x over Integer Ring
        sage: GT_y_QQ
        Generic Term Monoid over Monomial Growth Group in y over Rational Field
    """

    # enable the category framework for elements
    Element = GenericTerm

    def __init__(self, growth_group=None, category=None):
        r"""
        See :class:`GenericTermMonoid` for more information.
        
        EXAMPLES::
        
            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG_x = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: GT_x = atm.GenericTermMonoid(MG_x); GT_x
            Generic Term Monoid over Monomial Growth Group in x over Integer Ring
            sage: GT_x.growth_group()
            Monomial Growth Group in x over Integer Ring
            sage: MG_y = agg.MonomialGrowthGroup(QQ, 'y')
            sage: GT_y = atm.GenericTermMonoid(MG_y); GT_y
            Generic Term Monoid over Monomial Growth Group in y over Rational Field
            sage: GT_x is GT_y
            False

        ::

            sage: atm.GenericTermMonoid()
            Traceback (most recent call last):
            ...
            ValueError: Growth Group has to be specified
        """

        from sage.categories.monoids import Monoids
        from sage.categories.posets import Posets
        from sage.groups.asymptotic_growth_group import GenericGrowthGroup

        if category is None:
            category = Monoids() & Posets()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(Monoids() & Posets()) for cat in
                       category):
                raise ValueError('%s is not a subcategory of %s'
                                 % (category, Monoids() & Posets()))
        if growth_group is None:
            raise ValueError('Growth Group has to be specified')
        else:
            if not isinstance(growth_group, GenericGrowthGroup):
                raise ValueError('%s does not inherit from %s'
                                 % (growth_group, GenericGrowthGroup()))
        self._growth_group = growth_group
        super(GenericTermMonoid, self).__init__(category=category)

    def growth_group(self):
        r"""
        Return the growth group underlying this term monoid.

        INPUT:

        Nothing.

        OUTPUT:

        An asymptotic growth group.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: atm.ExactTermMonoid(MG, ZZ).growth_group()
            Monomial Growth Group in x over Integer Ring
        """
        return self._growth_group

    def _repr_(self):
        r"""
        Represent the generic term monoid as 'Generic Term Monoid
        over `growth_group`'.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: GG = agg.GenericGrowthGroup(ZZ)
            sage: atm.GenericTermMonoid(GG)._repr_()
            'Generic Term Monoid over Generic Growth Group over Integer Ring'
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: atm.GenericTermMonoid(growth_group=MG)._repr_()
            'Generic Term Monoid over Monomial Growth Group in x over Integer Ring'
        """
        return 'Generic Term Monoid over %s' % repr(self.growth_group())

    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this term monoid.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            Another GenericTermMonoid ``S`` coerces into this term
            monoid if and only if the growth group of ``S`` coerces
            into the growth group of this term monoid.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: T_ZZ = atm.GenericTermMonoid(growth_group=MG_ZZ); T_ZZ
            Generic Term Monoid over Monomial Growth Group in x over Integer Ring
            sage: MG_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: T_QQ = atm.GenericTermMonoid(growth_group=MG_QQ); T_QQ
            Generic Term Monoid over Monomial Growth Group in x over Rational Field
            sage: T_QQ._coerce_map_from_(T_ZZ)
            True
        """
        if isinstance(S, self.__class__):
            if self.growth_group().coerce_map_from(
                    S.growth_group()) is not None:
                return True

    def _element_constructor_(self, x):
        r"""
        Convert the given object to to this term monoid.

        INPUT:

        - ``x`` -- an object representing the element to be
          initialized.

        OUTPUT:

        An element of this term monoid.

        .. NOTE::

            The object ``x`` is either an asymptotic term that is
            to be coerced into this term monoid, or an asymptotic
            growth element that is used for creating an element
            of this term monoid.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: MG_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P_ZZ = atm.GenericTermMonoid(growth_group=MG_ZZ)
            sage: P_QQ = atm.GenericTermMonoid(growth_group=MG_QQ)
            sage: term1 = P_ZZ(MG_ZZ.gen())
            sage: term2 = P_QQ(MG_QQ.gen()^2)
            sage: term1 <= term2  # in order for two terms to be compared,
            ....:                 # a coercion into one of the parents
            ....:                 # has to be found.
            True
            sage: P_QQ.coerce(term1)  # coercion does not fail
            Generic Term with growth x
        """
        from sage.groups.asymptotic_growth_group import GenericGrowthElement

        if isinstance(x, GenericGrowthElement):
            return self.element_class(self, x)

        if x is None:
            raise ValueError('The growth of the term has to be specified!')
        elif x.parent() is self:
            return x
        elif isinstance(x, GenericTerm):
            return self.element_class(self,
                                      self.growth_group().coerce(x.growth))
        else:
            raise ValueError('Input is not an asymptotic growth element.')

    def _an_element_(self):
        r"""
        Return an element of this term monoid.

        INPUT:

        Nothing.

        OUTPUT:

        An element of this term monoid.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: atm.OTermMonoid(MG).an_element()  # indirect doctest
            O(x)
            sage: atm.GenericTermMonoid(MG).an_element()  # indirect doctest
            Generic Term with growth x
        """
        return self(self.growth_group().an_element())

    def le(self, left, right):
        r"""
        Return whether the growth of term ``left`` is less than or equal
        to the growth of term ``right``.

        INPUT:

        - ``left``, ``right`` -- elements of this term monoid.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: P = atm.GenericTermMonoid(growth_group=MG)
            sage: t1, t2 = P(x), P(x^2)
            sage: P.le(t1,t2)
            True
        """
        return left.growth <= right.growth


class OTerm(GenericTerm):
    r"""
    Class for an asymptotic term representing an `O` term with
    specified growth. For the mathematical properties of `O` terms
    see :wikipedia:`Big_O_Notation`.

    `O` terms may *absorb* terms of weaker or equal growth.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an element of the parent's ``growth_group``,
      specifying the growth of the asymptotic term.

    OUTPUT:

    An asymptotic `O` term.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
        sage: OT = atm.OTermMonoid(MG)
        sage: t1, t2, t3 = OT(x^-7), OT(x^5), OT(x^42)
        sage: t1, t2, t3
        (O(x^(-7)), O(x^5), O(x^42))
        sage: t1.can_absorb(t2)
        False
        sage: t2.can_absorb(t1)
        True
        sage: t2.absorb(t1)
        O(x^5)
        sage: t1 <= t2 and t2 <= t3
        True
        sage: t3 <= t1
        False
    """

    def _repr_(self):
        r"""
        Represent this O term as ``O(growth)``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: OT = atm.OTermMonoid(MG)
            sage: t1, t2, t3 = OT(x), OT(x^2), OT(x^3)
            sage: t1._repr_(), t2._repr_()
            ('O(x)', 'O(x^2)')
            sage: t3
            O(x^3)
        """
        return 'O(%s)' % self.growth


    def _absorb_(self, other):
        r"""
        Let this `O` term absorb another `O` term ``other``.

        INPUT:

        - ``other`` -- an asymptotic `O` term.

        OUTPUT:

        An asymptotic `O` term.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element and ``other``
            have the same parent.

            Also, observe that the result of a "dominant" `O` term
            absorbing another `O` term, always is the "dominant"
            `O` term again.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: OT = atm.OTermMonoid(growth_group=MG)
            sage: ot1, ot2 = OT(x), OT(x^2)
            sage: ot1.absorb(ot1)
            O(x)
            sage: ot2.absorb(ot1)
            O(x^2)
            sage: ot1.absorb(ot2)
            Traceback (most recent call last):
            ...
            ArithmeticError: O(x) cannot absorb O(x^2)
        """
        return self


class OTermMonoid(GenericTermMonoid):
    r"""
    Parent for asymptotic big `O` terms.

    INPUT:

    - ``growth_group`` -- the underlying asymptotic growth group,
      compare :class:`GenericTermMonoid`.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    OUTPUT:

    A monoid for asymptotic `O` terms.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: MG_x_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
        sage: MG_y_QQ = agg.MonomialGrowthGroup(QQ, 'y')
        sage: OT_x_ZZ = atm.OTermMonoid(MG_x_ZZ); OT_x_ZZ
        Asymptotic O Term Monoid over Monomial Growth Group in x over Integer Ring
        sage: OT_y_QQ = atm.OTermMonoid(MG_y_QQ); OT_y_QQ
        Asymptotic O Term Monoid over Monomial Growth Group in y over Rational Field
    """
    # enable the category framework for elements
    Element = OTerm


    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this term monoid.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        ``True`` or ``None``.

        .. NOTE::

            Another term monoid ``S`` coerces into this term monoid
            if ``S`` is an instance of one of the following classes:

            - :class:`OTermMonoid`

            - :class:`LTermGenericMonoid` or a descendant thereof

            - :class:`ExactTermMonoid`

            Additionally, the growth group underlying ``S`` has to
            coerce into the growth group of this term monoid.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG_ZZ = agg.MonomialGrowthGroup(ZZ, 'x'); x_ZZ = MG_ZZ.gen()
            sage: MG_QQ = agg.MonomialGrowthGroup(QQ, 'x'); x_QQ = MG_QQ.gen()
            sage: OT_ZZ = atm.OTermMonoid(MG_ZZ)
            sage: OT_QQ = atm.OTermMonoid(MG_QQ)
            sage: ET = atm.ExactTermMonoid(MG_ZZ, ZZ)

        Now, the ``OTermMonoid`` whose growth group is over the
        Integer Ring has to coerce into the ``OTermMonoid`` with
        the growth group over the rational field, and the
        ``ExactTermMonoid`` also has to coerce in both the
        ``OTermMonoid``'s::

            sage: OT_QQ._coerce_map_from_(OT_ZZ)
            True
            sage: OT_QQ._coerce_map_from_(ET)
            True
            sage: ET._coerce_map_from_(OT_ZZ) is None
            True
        """
        if isinstance(S, (ExactTermMonoid, LTermGenericMonoid)):
            if self.growth_group().coerce_map_from(
                    S.growth_group()) is not None:
                return True
        else:
            return super(OTermMonoid, self)._coerce_map_from_(S)


    def _repr_(self):
        r"""
        Represent this `O` term monoid as 'Asymptotic O term monoid
        over ``growth_group``'.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: atm.OTermMonoid(MG)._repr_()
            'Asymptotic O Term Monoid over Monomial Growth Group in x over Integer Ring'
        """
        return 'Asymptotic O Term Monoid over %s' % self.growth_group()


class TermWithCoefficient(GenericTerm):
    r"""
    Base class for asymptotic terms possessing a coefficient. For
    example, the classes :class:`ExactTerm` and :class:`LTermGeneric`
    both directly inherit from this class.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element from
      ``parent.growth_group()``.

    - ``coefficient`` -- an element from ``parent.base_ring``.

    OUTPUT:

    An asymptotic term with coefficient.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
        sage: CT_ZZ = atm.TermWithCoefficientMonoid(MG, ZZ)
        sage: CT_QQ = atm.TermWithCoefficientMonoid(MG, QQ)
        sage: CT_ZZ(x^2, 5)
        Asymptotic Term with coefficient 5 and growth x^2
        sage: CT_QQ(x^3, 3/8)
        Asymptotic Term with coefficient 3/8 and growth x^3
    """

    def __init__(self, parent, growth, coefficient=1):
        r"""
        See :class:`TermWithCoefficient` for more information.

        EXAMPLES:

        First, we define some parent monoids::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: CT_ZZ = atm.TermWithCoefficientMonoid(MG, ZZ)
            sage: CT_QQ = atm.TermWithCoefficientMonoid(MG, QQ)

        The coefficients have to be from the given base ring::

            sage: t = CT_ZZ(x, 1/2)
            Traceback (most recent call last):
            ...
            ValueError: 1/2 is not in Integer Ring
            sage: t = CT_QQ(x, 1/2); t
            Asymptotic Term with coefficient 1/2 and growth x

        For technical reasons, the coefficient 0 is not allowed::

            sage: t = CT_ZZ(x^42, 0)
            Traceback (most recent call last):
            ...
            ValueError: 0 is not a valid coefficient
        """
        if coefficient not in parent.base_ring():
            raise ValueError('%s is not in %s' % (coefficient,
                                                  parent.base_ring()))
        elif coefficient == 0:
            raise ValueError('0 is not a valid coefficient')
        else:
            self.coefficient = coefficient
        super(TermWithCoefficient, self).__init__(parent=parent, growth=growth)


    def _repr_(self):
        r"""
        Represent the :class:`TermWithCoefficient` as 'Asymptotic
        Term with coefficient ``coefficient`` and growth ``growth``'.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: P = atm.TermWithCoefficientMonoid(MG, ZZ)
            sage: P(x^2, 5)._repr_()
            'Asymptotic Term with coefficient 5 and growth x^2'
        """
        return 'Asymptotic Term with coefficient %s and growth %s' % \
               (self.coefficient, self.growth)


    def _mul_(self, other):
        r"""
        Multiplication method for asymptotic terms with coefficients.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        An asymptotic term representing the product of this element
        and ``other``.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element and ``other`` have
            the same parent.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: CT = atm.TermWithCoefficientMonoid(MG, ZZ)
            sage: ET = atm.ExactTermMonoid(MG, ZZ)

            This method handels the multiplication of abstract terms
            with coefficient (i.e. :class:`TermWithCoefficient`) and
            exact terms (i.e. :class:`ExactTerm`). First, an example
            for abstract terms::

            sage: t1, t2 = CT(x^2, 2), CT(x^3, 3)
            sage: t1 * t2
            Asymptotic Term with coefficient 6 and growth x^5

            And now, an example for exact terms::

            sage: t1, t2 = ET(x^2, 2), ET(x^3, 3)
            sage: t1 * t2
            6 * x^5
        """
        return self.parent()(self.growth * other.growth,
                             self.coefficient * other.coefficient)


    def _le_(self, other):
        r"""
        Return whether this asymptotic term with coefficient grows
        at most (less or equal) like ``other``.

        INPUT:

        - ``other`` -- an asymptotic term with coefficient.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element and ``other`` are
            from the same parent.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: ET = atm.ExactTermMonoid(MG, QQ)
            sage: t1, t2, t3 = ET(x, 5), ET(x^2, 3), ET(x^2, 42)
            sage: t1 <= t2
            True
            sage: t2 <= t1
            False
            sage: t2 <= t3
            True
            sage: t3 <= t2
            False
        """
        if self.growth == other.growth:
            return self.coefficient <= other.coefficient
        else:
            return super(TermWithCoefficient, self)._le_(other)


class TermWithCoefficientMonoid(GenericTermMonoid):
    r"""
    This class implements the base structure for parents of
    asymptotic terms possessing a coefficient from some coefficient
    ring. In particular, this is also the parent for
    :class:`TermWithCoefficient`.

    INPUT:

    - ``growth_group`` -- The underlying asymptotic growth group,
      compare :class:`GenericTermMonoid`.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    - ``base_ring`` -- the ring which contains the
      coefficients of the elements.

    OUTPUT:

    A monoid for asymptotic terms with coefficients.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: MG_ZZ = agg.MonomialGrowthGroup(ZZ, 'x'); x_ZZ = MG_ZZ.gen()
        sage: MG_QQ = agg.MonomialGrowthGroup(QQ, 'x'); x_QQ = MG_QQ.gen()
        sage: TC_ZZ = atm.TermWithCoefficientMonoid(MG_ZZ, QQ); TC_ZZ
        Monoid for asymptotic terms with coefficients from Rational Field over Monomial Growth Group in x over Integer Ring
        sage: TC_QQ = atm.TermWithCoefficientMonoid(MG_QQ, QQ); TC_QQ
        Monoid for asymptotic terms with coefficients from Rational Field over Monomial Growth Group in x over Rational Field
        sage: TC_ZZ == TC_QQ or TC_ZZ is TC_QQ
        False
        sage: TC_QQ.coerce_map_from(TC_ZZ)
        Conversion map:
          From: Monoid for asymptotic terms with coefficients from Rational Field over Monomial Growth Group in x over Integer Ring
          To:   Monoid for asymptotic terms with coefficients from Rational Field over Monomial Growth Group in x over Rational Field
    """

    # enable the category framework for elements
    Element = TermWithCoefficient

    def __init__(self, growth_group=None, base_ring=None,
                 category=None):
        r"""
        For more information see :class:`TermWithCoefficientMonoid`.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: P_ZZ = atm.TermWithCoefficientMonoid(MG, ZZ); P_ZZ
            Monoid for asymptotic terms with coefficients from Integer Ring over Monomial Growth Group in x over Integer Ring
            sage: P_QQ = atm.TermWithCoefficientMonoid(MG, QQ); P_QQ
            Monoid for asymptotic terms with coefficients from Rational Field over Monomial Growth Group in x over Integer Ring
            sage: P_QQ.category()
            Join of Category of monoids and Category of posets
        """
        if base_ring is None:
            raise ValueError('Base ring is not specified')
        self._base_ring = base_ring
        super(TermWithCoefficientMonoid,
              self).__init__(growth_group=growth_group, category=category)


    def base_ring(self):
        r"""
        Return the base ring of this term monoid, i.e. the ring where
        the coefficients are from.

        INPUT:

        Nothing.

        OUTPUT:

        A ring.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: atm.ExactTermMonoid(MG, ZZ).base_ring()  # indirect doctest
            Integer Ring
        """
        return self._base_ring

    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this term monoid.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            Another term monoid ``S`` coerces into this exact term
            monoid if both, the base ring as well as the growth
            group underlying ``S`` coerce into the base ring and the
            growth group underlying this term monoid.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: MG_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: TC_ZZ = atm.TermWithCoefficientMonoid(MG_ZZ, ZZ); TC_ZZ
            Monoid for asymptotic terms with coefficients from Integer Ring over Monomial Growth Group in x over Integer Ring
            sage: TC_QQ = atm.TermWithCoefficientMonoid(MG_QQ, QQ); TC_QQ
            Monoid for asymptotic terms with coefficients from Rational Field over Monomial Growth Group in x over Rational Field
            sage: TC_QQ._coerce_map_from_(TC_ZZ)
            True
            sage: TC_ZZ._coerce_map_from_(TC_QQ) is None
            True
        """
        if isinstance(S, TermWithCoefficientMonoid):
            return (super(TermWithCoefficientMonoid, self).
                    _coerce_map_from_(S) and self.base_ring().
                    coerce_map_from(S.base_ring()) is not None)


    def _element_constructor_(self, x, coefficient=None):
        r"""
        Construct an asymptotic term with coefficient or convert
        the given object ``x`` to this term monoid.

        INPUT:

        - ``x`` -- a growth element or an object representing the
          element to be initialized.

        - ``coefficient`` -- an element of the ``base_ring``.

        OUTPUT:

        An asymptotic term.

        .. NOTE::

            The object ``x`` is either an asymptotic term with
            coefficient that is to be coerced into this term monoid,
            or an asymptotic growth element that is used together
            with ``coefficient`` in order to create an element of
            this term monoid.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: P = atm.TermWithCoefficientMonoid(MG, ZZ)
            sage: t1 = P(x^2, 5); t1  # indirect doctest
            Asymptotic Term with coefficient 5 and growth x^2

        TESTS::

            sage: P(x)
            Traceback (most recent call last):
            ...
            ValueError: The coefficient is not specified
        """
        if isinstance(x, TermWithCoefficient):
            return self.element_class(self, x.growth, x.coefficient)
        elif coefficient is None:
            raise ValueError('The coefficient is not specified')
        else:
            return self.element_class(self, x, coefficient)


    def _repr_(self):
        r"""
        Represent this term monoid as 'Monoid for asymptotic terms with
        coefficients from ``base_ring`` over
        ``growth_group``'.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: atm.TermWithCoefficientMonoid(MG, ZZ)._repr_()
            'Monoid for asymptotic terms with coefficients from Integer Ring over Monomial Growth Group in x over Integer Ring'
        """
        return 'Monoid for asymptotic terms with coefficients from %s ' \
               'over %s' % (self.base_ring(), self.growth_group())


    def _an_element_(self):
        r"""
        Return an element of this term monoid with coefficient.

        INPUT:

        Nothing.

        OUTPUT:

        An element of this term monoid.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: atm.TermWithCoefficientMonoid(MG, ZZ).an_element()  # indirect doctest
            Asymptotic Term with coefficient 1 and growth x
            sage: atm.ExactTermMonoid(MG, ZZ).an_element()  # indirect doctest
            1 * x
            sage: atm.LTermGenericMonoid(MG, ZZ).an_element()  # indirect doctest
            1 * L(x, 0)
        """
        return self(self.growth_group().an_element(),
                    self.base_ring().an_element())


class LTermGeneric(TermWithCoefficient):
    r"""
    Base class for asymptotic `L` terms, i.e. big `O` terms with
    an explicitly specified constant and starting point.

    .. NOTE::

        When adding ('absorbing') various different `L` terms, some
        sort of "lifting" occurs in general. This lifting affects the
        coefficient of the resulting `L` term, and works
        differently for various growth classes. This makes different
        implementations for different growth classes necessary.

    .. TODO::

        Implement *absorption* (:meth:`_absorb_`) for `L` terms.

    INPUT:

        - ``parent`` -- the parent of the asymptotic term.

        - ``growth`` -- an asymptotic growth element from
          ``parent.growth_group()``.

        - ``coefficient`` -- an element from
          ``parent.base_ring``.

        - ``start`` -- a real number representing the starting point
          of the estimations in the definition of asymptotic `O`
          terms (see :wikipedia:`Big_O_Notation`).

    OUTPUT:

    An asymptotic `L` term with given coefficient and starting point.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
        sage: LT_ZZ = atm.LTermGenericMonoid(MG, ZZ)
        sage: LT_QQ = atm.LTermGenericMonoid(MG, QQ)
        sage: lt_ZZ, lt_QQ = LT_ZZ(x^2, 42, 42), LT_QQ(x, 12/7, 0); lt_ZZ, lt_QQ
        (42 * L(x^2, 42), 12/7 * L(x, 0))
        sage: lt_ZZ <= lt_QQ
        False
        sage: lt_QQ <= lt_ZZ
        True
        sage: lt_ZZ * lt_QQ
        72 * L(x^3, 42)
        sage: lt_QQ.can_absorb(lt_ZZ)
        True
        sage: lt_QQ.absorb(lt_ZZ)
        Traceback (most recent call last):
        ...
        NotImplementedError: Not yet implemented
    """

    def __init__(self, parent, growth, coefficient=1, start=0):
        r"""
        See :class:`LTermGeneric` for more information and examples.
        """
        from sage.rings.real_mpfr import RR

        if start not in RR:
            raise ValueError('%s is not a real number' % start)
        else:
            self.start = start
        super(LTermGeneric, self).__init__(parent=parent, growth=growth,
                                           coefficient=coefficient)


    def _repr_(self):
        r"""
        Represent the generic `L` term as "``coefficient * L(growth,
        start)``".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: LT = atm.LTermGenericMonoid(MG, ZZ)
            sage: LT(x^2, 5, start=0)._repr_()
            '5 * L(x^2, 0)'
        """
        return '%s * L(%s, %s)' % (self.coefficient, self.growth, self.start)


    def _absorb_(self, other):
        r"""
        Absorb the `L` term ``other`` by this `L` term.

        INPUT:

        - ``other`` -- an asymptotic `L` term.

        OUTPUT:

        An asymptotic `L` term.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element and ``other`` have
            the same parent.

        .. TODO::

            Implement this method for specific growth parents.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: LT = atm.LTermGenericMonoid(MG, QQ)
            sage: t1, t2 = LT(x^2, 5, 0), LT(x, 42, 5); t1, t2
            (5 * L(x^2, 0), 42 * L(x, 5))
            sage: t1.can_absorb(t2)
            True
            sage: t1.absorb(t2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Not yet implemented
        """
        raise NotImplementedError('Not yet implemented')


    def _mul_(self, other):
        r"""
        Multiply this `L` term with ``other``.

        INPUT:

        - ``other`` -- an asymptotic `L` term.

        OUTPUT:

        An asymptotic `L` term.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element and ``other`` have
            the same parent.

            Furthermore, When taking the product of two asymptotic
            `L` terms, the growth of the product is the product of
            the growth elements of the factors, and analogue for the
            respective coefficient. The starting point of the
            resulting `L` term is the maximum starting point of the
            factors.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: LT = atm.LTermGenericMonoid(MG, QQ)
            sage: lt1, lt2 = LT(x^2, 2, 0), LT(x^-1, 3, 5)
            sage: lt1._mul_(lt1)
            4 * L(x^4, 0)
            sage: lt1._mul_(lt2)
            6 * L(x, 5)
        """
        result = super(LTermGeneric, self)._mul_(other)
        result.start = max(self.start, other.start)
        return result


class LTermGenericMonoid(TermWithCoefficientMonoid):
    r"""
    Base class for parents of asymptotic `L` terms. Also, the parent
    for :class:`LTermGeneric`. `L` terms are asymptotic terms which
    behave like `O` terms, with the difference that for `L` terms,
    the constant and the starting point for the inequality occurring
    in the definition of `O` terms (see :wikipedia:`Big_O_notation`)
    are explicitly given.

    INPUT:

    - ``growth_group`` -- The underlying asymptotic growth group,
      compare :class:`GenericTermMonoid`.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    - ``base_ring`` -- the ring which contains the coefficients of
      the elements.

    OUTPUT:

    A generic asymptotic `L` term monoid.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: MG_ZZ = agg.MonomialGrowthGroup(ZZ, 'x'); x_ZZ = MG_ZZ.gen()
        sage: MG_QQ = agg.MonomialGrowthGroup(QQ, 'x'); x_QQ = MG_QQ.gen()
        sage: LT_ZZ = atm.LTermGenericMonoid(MG_ZZ, QQ); LT_ZZ
        Generic L Term Monoid with coefficients from Rational Field over Monomial Growth Group in x over Integer Ring
        sage: LT_QQ = atm.LTermGenericMonoid(MG_QQ, QQ); LT_QQ
        Generic L Term Monoid with coefficients from Rational Field over Monomial Growth Group in x over Rational Field
        sage: LT_ZZ == LT_QQ or LT_ZZ is LT_QQ
        False
        sage: LT_QQ.coerce_map_from(LT_ZZ)
        Conversion map:
          From: Generic L Term Monoid with coefficients from Rational Field over Monomial Growth Group in x over Integer Ring
          To:   Generic L Term Monoid with coefficients from Rational Field over Monomial Growth Group in x over Rational Field
    """
    # enable the category framework for elements
    Element = LTermGeneric


    def _element_constructor_(self, x, coefficient=None, start=0):
        r"""
        Construct a generic `L` term or convert the given
        object ``x`` to this term monoid.

        INPUT:

        - ``x`` -- a growth element or an object representing the
          element to be initialized.

        - ``coefficient`` -- an element of the ``base_ring``.

        - ``start`` -- a real number indicating the point where
          the `L` term is valid.

        OUTPUT:

        A generic `L` term.

        .. NOTE::

            The object ``x`` is either a `L` term that is to be
            coerced into this term monoid, or an asymptotic growth
            element that is used together with ``coefficient`` and
            ``start`` in order to create an element of this term
            monoid.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: LT = atm.LTermGenericMonoid(MG, QQ)
            sage: ET = atm.ExactTermMonoid(MG, QQ)
            sage: lt, et = LT(x^3, 42/5, 3), ET(x^7, 5/9)
            sage: LT(lt)
            42/5 * L(x^3, 3)
            sage: LT(lt) == lt
            True
            sage: LT(et)
            5/9 * L(x^7, 0)
        """
        if isinstance(x, LTermGeneric):
            return self.element_class(self, x.growth, x.coefficient, x.start)
        elif isinstance(x, ExactTerm):
            return self.element_class(self, x.growth, x.coefficient, 0)
        else:
            return self.element_class(self, x, coefficient, start)


    def _repr_(self):
        r"""
        Represent the generic `L` term monoid as 'Generic L Term
        Monoid with coefficients from ``base_ring`` over
        ``growth_group``'.

        INPUT:

        Noting.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(QQ, 'x'); x = MG.gen()
            sage: atm.LTermGenericMonoid(MG, ZZ)._repr_()
            'Generic L Term Monoid with coefficients from Integer Ring over Monomial Growth Group in x over Rational Field'
        """
        return 'Generic L Term Monoid with coefficients from %s over %s' % \
               (self.base_ring(), self.growth_group())


    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this term monoid.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            Another term monoid ``S`` coerces into this `L` term
            monoid if both, the base ring as well as the growth group
            underlying ``S`` coerce into the base ring and the growth
            group underlying this term monoid.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: LT_QQ = atm.LTermGenericMonoid(MG, QQ)
            sage: LT_ZZ = atm.LTermGenericMonoid(MG, ZZ)
            sage: ET = atm.ExactTermMonoid(MG, ZZ)
            sage: LT_QQ.coerce_map_from(LT_ZZ)
            Conversion map:
              From: Generic L Term Monoid with coefficients from Integer Ring over Monomial Growth Group in x over Integer Ring
              To:   Generic L Term Monoid with coefficients from Rational Field over Monomial Growth Group in x over Integer Ring
            sage: LT_QQ.coerce_map_from(ET)
            Conversion map:
              From: Exact Term Monoid with coefficients from Integer Ring over Monomial Growth Group in x over Integer Ring
              To:   Generic L Term Monoid with coefficients from Rational Field over Monomial Growth Group in x over Integer Ring
        """
        if isinstance(S, ExactTermMonoid):
            if self.growth_group().coerce_map_from(S.growth_group()) is not None \
                    and self.base_ring().coerce_map_from(S.base_ring()) \
                    is not None:
                return True
        return super(LTermGenericMonoid, self)._coerce_map_from_(S)


class ExactTerm(TermWithCoefficient):
    r"""
    Class for asymptotic exact terms. These terms primarily consist of
    an asymptotic growth element as well as a coefficient specifying
    the growth of the asymptotic term.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element from
      ``parent.growth_group()``.

    - ``coefficient`` -- an element from ``parent.base_ring``.

    OUTPUT:

    An asymptotic exact term.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
        sage: ET = atm.ExactTermMonoid(MG, QQ)

    Asymptotic exact terms may be multiplied (with the usual rules
    applying)::

        sage: ET(x^2, 3) * ET(x, -1)
        -3 * x^3
        sage: ET(x^0, 4) * ET(x^5, 2)
        8 * x^5

    They may also be multiplied with `L` or `O` terms::

        sage: OT = atm.OTermMonoid(MG)
        sage: LT = atm.LTermGenericMonoid(MG, QQ)
        sage: ET(x^2, 42) * OT(x)
        O(x^3)
        sage: ET(x^2, 42) * LT(x, 1, 5)
        42 * L(x^3, 5)

    Absorption for asymptotic exact terms relates to addition::

        sage: ET(x^2, 5).can_absorb(ET(x^5, 12))
        False
        sage: ET(x^2, 5).can_absorb(ET(x^2, 1))
        True
        sage: ET(x^2, 5).absorb(ET(x^2, 1))
        6 * x^2

    Note that, as for technical reasons, `0` is not allowed as a
    coefficient for an asymptotic term with coefficient, ``None``
    is returned if two asymptotic exact terms cancel out each other
    during absorption::

        sage: ET(x^2, 42).can_absorb(ET(x^2, -42))
        True
        sage: ET(x^2, 42).absorb(ET(x^2, -42)) is None
        True
    """

    def _repr_(self):
        r"""
        Represent this exact term as ``coefficient *
        growth``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: ET = atm.ExactTermMonoid(MG, ZZ)
            sage: et1 = ET(x^2, 2); et1
            2 * x^2
        """
        return '%s * %s' % (self.coefficient, self.growth)


    def _absorb_(self, other):
        r"""
        Let this exact term absorb another exact term ``other``.

        INPUT:

        - ``other`` -- an exact term.

        OUTPUT:

        An exact term or ``None``.

        .. NOTE::

            In the context of exact terms, absorption translates
            to addition. As the coefficient `0` is not allowed,
            ``None`` is returned if the terms cancel out.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: ET = atm.ExactTermMonoid(MG, QQ)

        Asymptotic exact terms can absorb other asymptotic exact
        terms with the same growth::

            sage: et1, et2 = ET(x^2, 2), ET(x^2, -2)
            sage: et1.absorb(et1)
            4 * x^2
            sage: repr(et1.absorb(et2))
            'None'

        If the growth differs, an ``ArithmeticException`` is raised::

            sage: ET(x^5, 1).absorb(et1)
            Traceback (most recent call last):
            ...
            ArithmeticError: 1 * x^5 cannot absorb 2 * x^2
        """
        cls = self.__class__
        coef_new = self.coefficient + other.coefficient
        if coef_new == 0:
            return None
        else:
            return cls(self.parent(), self.growth, coef_new)


class ExactTermMonoid(TermWithCoefficientMonoid):
    r"""
    Parent for asymptotic exact terms, implemented in
    :class:`ExactTerm`.

    INPUT:

    - ``growth_group`` -- The underlying asymptotic growth group, compare
      :class:`GenericTermMonoid`.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    - ``base_ring`` -- the ring which contains the coefficients of
      the elements.

    OUTPUT:

    A monoid for asymptotic exact terms.

    EXAMPLES::


        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: MG_ZZ = agg.MonomialGrowthGroup(ZZ, 'x'); x_ZZ = MG_ZZ.gen()
        sage: MG_QQ = agg.MonomialGrowthGroup(QQ, 'x'); x_QQ = MG_QQ.gen()
        sage: ET_ZZ = atm.ExactTermMonoid(MG_ZZ, ZZ); ET_ZZ
        Exact Term Monoid with coefficients from Integer Ring over Monomial Growth Group in x over Integer Ring
        sage: ET_QQ = atm.ExactTermMonoid(MG_QQ, QQ); ET_QQ
        Exact Term Monoid with coefficients from Rational Field over Monomial Growth Group in x over Rational Field
        sage: ET_QQ.coerce_map_from(ET_ZZ)
        Conversion map:
          From: Exact Term Monoid with coefficients from Integer Ring over Monomial Growth Group in x over Integer Ring
          To:   Exact Term Monoid with coefficients from Rational Field over Monomial Growth Group in x over Rational Field
    """
    # enable the category framework for elements
    Element = ExactTerm

    def _repr_(self):
        r"""
        Represent the asymptotic exact term monoid as 'Exact Term
        Monoid with coefficients from ``base_ring`` over
        ``growth_group``'.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x'); x = MG.gen()
            sage: atm.ExactTermMonoid(MG, QQ)._repr_()
            'Exact Term Monoid with coefficients from Rational Field over Monomial Growth Group in x over Integer Ring'
        """
        return 'Exact Term Monoid with coefficients from %s over %s' % \
               (self.base_ring(), self.growth_group())


class TermMonoidFactory(sage.structure.factory.UniqueFactory):
    r"""
    Factory for asymptotic term monoids. Generates an
    :class:`OTermMonoid`, an :class:`ExactTermMonoid`, or
    an `L` term monoid.

    INPUT:

    - ``term`` -- The kind of term that shall be created. Either
      'exact', 'O', or 'L'.

    - ``growth_group`` -- A growth group.

    - ``base_ring`` -- The base ring for coefficients.

    - ``start`` -- A real number representing the starting point of
      the estimations in the definition of asymptotic `O` terms (see
      :wikipedia:`Big_O_Notation`). Only for `L` terms.

    OUTPUT:

    An asymptotic term monoid.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
        sage: OT = atm.TermMonoid('O', MG); OT
        Asymptotic O Term Monoid over Monomial Growth Group in x over Integer Ring
        sage: ET = atm.TermMonoid('exact', MG, ZZ); ET
        Exact Term Monoid with coefficients from Integer Ring over Monomial Growth Group in x over Integer Ring
    """
    def create_key_and_extra_args(self, term, growth_group, base_ring=None,
                                  **kwds):
        r"""
        Given the arguments and keyword, create a key that uniquely
        determines this object.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: atm.TermMonoid.create_key_and_extra_args('O', MG)
            (('O', Monomial Growth Group in x over Integer Ring, None), {})
            sage: atm.TermMonoid.create_key_and_extra_args('exact', MG, ZZ)
            (('exact', Monomial Growth Group in x over Integer Ring, Integer Ring), {})
            sage: atm.TermMonoid.create_key_and_extra_args('L', MG, None)
            Traceback (most recent call last):
            ...
            ValueError: a base ring has to be specified
        """
        if term not in ['O', 'exact', 'L']:
            raise ValueError("%s has to be either 'exact', 'O', or 'L'" % term)

        from sage.groups.asymptotic_growth_group import GenericGrowthGroup
        if not isinstance(growth_group, GenericGrowthGroup):
            raise ValueError("%s has to be an asymptotic growth group"
                             % growth_group)

        if term == 'O':
            if base_ring is not None:
                raise ValueError("O term monoids do not require a base ring")
        else:
            if base_ring is None:
                raise ValueError("a base ring has to be specified")

        return (term, growth_group, base_ring), kwds

    def create_object(self, version, key, **kwds):
        r"""
        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: atm.TermMonoid('O', MG)
            Asymptotic O Term Monoid over Monomial Growth Group in x over Integer Ring
            sage: atm.TermMonoid('exact', MG, ZZ)
            Exact Term Monoid with coefficients from Integer Ring over Monomial Growth Group in x over Integer Ring
        """

        term, growth_group, base_ring = key
        if term == 'O':
            return OTermMonoid(growth_group, **kwds)
        elif term == 'exact':
            return ExactTermMonoid(growth_group, base_ring, **kwds)
        else:
            # in this case we are generating an L term monoid. For now,
            # just return a GenericLTermMonoid. as soon as L terms are
            # actually implemented, L term monoids relative to the
            # growth group have to be constructed.
            return LTermGenericMonoid(growth_group, base_ring, **kwds)


TermMonoid = TermMonoidFactory("TermMonoid")