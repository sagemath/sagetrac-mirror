r"""
Growth Groups as Cartesian Products

AUTHORS:

- Daniel Krenn (2015-06-02): cartesian products
- Benjamin Hackl (2015-07)

.. WARNING::

    As this code is experimental, warnings are thrown when a growth
    group is created for the first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G = agg.GenericGrowthGroup(ZZ); G
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        Growth Group Generic(ZZ)
        sage: G = agg.MonomialGrowthGroup(ZZ, 'x'); G
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        Growth Group x^ZZ
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

import sage


class CartesianProductFactory(sage.structure.factory.UniqueFactory):
    r"""
    Create various types of cartesian products depending on its input.

    INPUT:

    - ``growth_groups`` -- a tuple (or other iterable) of growth groups.

    - ``order`` -- (default: ``None``) if specified, then this order
      is taken for comparing two cartesian product elements. If ``order`` is
      ``None`` this is determined automatically.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: A = GrowthGroup('x^ZZ'); A
        Growth Group x^ZZ
        sage: B = GrowthGroup('log(x)^ZZ'); B
        Growth Group log(x)^ZZ
        sage: C = cartesian_product([A, B]); C  # indirect doctest
        Growth Group x^ZZ * log(x)^ZZ
        sage: C._le_ == C.le_lex
        True
        sage: D = GrowthGroup('y^ZZ'); D
        Growth Group y^ZZ
        sage: E = cartesian_product([A, D]); E  # indirect doctest
        Growth Group x^ZZ * y^ZZ
        sage: E._le_ == E.le_product
        True
        sage: F = cartesian_product([C, D]); F  # indirect doctest
        Growth Group x^ZZ * log(x)^ZZ * y^ZZ
        sage: F._le_ == F.le_product
        True
        sage: cartesian_product([A, E]); G  # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: Growth groups (Growth Group x^ZZ, Growth Group x^ZZ * y^ZZ)
        do not have pairwise disjoint variables.
        sage: cartesian_product([A, B, D])  # indirect doctest
        Growth Group x^ZZ * log(x)^ZZ * y^ZZ

    TESTS::

        sage: from sage.rings.asymptotic.growth_group_cartesian import CartesianProductFactory
        sage: CartesianProductFactory('factory')([A, B], category=Sets())
        Growth Group x^ZZ * log(x)^ZZ
        sage: CartesianProductFactory('factory')([], category=Sets())
        Traceback (most recent call last):
        ...
        TypeError: Cannot create cartesian product without factors.
    """
    def create_key_and_extra_args(self, growth_groups, category, **kwds):
        r"""
        Given the arguments and keywords, create a key that uniquely
        determines this object.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group_cartesian import CartesianProductFactory
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: A = GrowthGroup('x^ZZ')
            sage: CartesianProductFactory('factory').create_key_and_extra_args(
            ....:     [A], category=Sets(), order='blub')
            (((Growth Group x^ZZ,), Category of sets), {'order': 'blub'})
        """
        return (tuple(growth_groups), category), kwds


    def create_object(self, version, args, **kwds):
        r"""
        Create an object from the given arguments.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: cartesian_product([GrowthGroup('x^ZZ')])  # indirect doctest
            Growth Group x^ZZ
        """
        growth_groups, category = args
        if not growth_groups:
            raise TypeError('Cannot create cartesian product without factors.')
        order = kwds.pop('order', None)
        if order is not None:
            return GenericProduct(growth_groups, category, order=order, **kwds)

        vg = tuple((g.variable_names(), g) for g in growth_groups)

        # check if all groups have a variable
        if not all(v for v, _ in vg):
            raise NotImplementedError('Growth groups %s have no variable.' %
                                      tuple(g for g in growth_groups
                                            if not g.variable_names()))

        # sort by variables
        from itertools import groupby, product
        vgs = tuple((v, tuple(gs)) for v, gs in
                    groupby(sorted(vg, key=lambda k: k[0]), key=lambda k: k[0]))

        # check if variables are pairwise disjoint
        for u, w in product(iter(v for v, _ in vgs), repeat=2):
            if u != w and set(u).intersection(set(w)):
                raise ValueError('Growth groups %s do not have pairwise disjoint '
                                 'variables.' % (growth_groups,))

        # build cartesian products
        u_groups = list()
        for _, gs in vgs:
            gs = tuple(g for _, g in gs)
            if len(gs) > 1:
                u_groups.append(UnivariateProduct(gs, category, **kwds))
            else:
                u_groups.append(gs[0])

        if len(u_groups) > 1:
            m_group = MultivariateProduct(tuple(u_groups), category, **kwds)
        else:
            m_group = u_groups[0]
        return m_group


CartesianProductGrowthGroups = CartesianProductFactory('CartesianProductGrowthGroups')


from sage.combinat.posets.cartesian_product import CartesianProductPosets
from sage.rings.asymptotic.growth_group import GenericGrowthGroup
class GenericProduct(CartesianProductPosets, GenericGrowthGroup):
    r"""
    A cartesian product of growth groups.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: P = agg.MonomialGrowthGroup(QQ, 'x')
        sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
        sage: C = cartesian_product([P, L], order='lex'); C  # indirect doctest
        Growth Group x^QQ * log(x)^ZZ
        sage: C.an_element()
        x^(1/2) * log(x)

    ::

        sage: Px = agg.MonomialGrowthGroup(QQ, 'x')
        sage: Lx = agg.MonomialGrowthGroup(ZZ, 'log(x)')
        sage: Cx = cartesian_product([Px, Lx], order='lex')  # indirect doctest
        sage: Py = agg.MonomialGrowthGroup(QQ, 'y')
        sage: C = cartesian_product([Cx, Py], order='product'); C  # indirect doctest
        Growth Group x^QQ * log(x)^ZZ * y^QQ
        sage: C.an_element()
        x^(1/2) * log(x) * y^(1/2)

    .. SEEALSO::

        :class:`~sage.sets.cartesian_product.CartesianProduct`,
        :class:`~sage.combinat.posets.cartesian_product.CartesianProductPosets`.
    """

    def _element_constructor_(self, data):
        r"""
        Converts the given object to an element of this cartesian
        product.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * y^ZZ')
            sage: G_log = GrowthGroup('x^ZZ * log(x)^ZZ * y^ZZ')

        Conversion from the symbolic ring works::

            sage: x,y = var('x y')
            sage: G(x^-3 * y^2)
            x^(-3) * y^2
            sage: G(x^4), G(y^2)
            (x^4, y^2)
            sage: G(1)
            1

        Even more complex expressions can be parsed::

            sage: G_log(x^42 * log(x)^-42 * y^42)
            x^42 * log(x)^(-42) * y^42

        TESTS::

            sage: G = GrowthGroup('x^ZZ * y^ZZ')
            sage: G('x'), G('y')
            (x, y)

        ::

            sage: G_log(log(x))
            log(x)
        """
        if data == 1:
            return self.one()

        if data is None:
            raise ValueError('%s cannot be converted.' % (data,))

        if isinstance(data, list):
            try:
                obj = super(GenericProduct,
                            self)._element_constructor_(data)
                return obj
            except ValueError:
                return self.prod(self(elem) for elem in data)

        if hasattr(data, 'parent'):
            P = data.parent()
            if P is self:
                return data

            elif P is sage.symbolic.ring.SR:
                import operator
                from sage.symbolic.operators import mul_vararg
                op = data.operator()
                if op == operator.pow or data.is_symbol() \
                        or isinstance(op, sage.functions.log.Function_log):
                    return self(self._convert_to_factor_(data))
                elif op == mul_vararg:
                    return self(data.operands())
            # room for other parents (e.g. polynomial ring et al.)

        # try to convert the input to one of the factors
        data_conv = self._convert_to_factor_(data)
        if data_conv is not None:
            factors = self.cartesian_factors()
            return self([data_conv if factor == data_conv.parent() else 1 for
                         factor in factors])

        # final attempt: try parsing the representation string
        str_lst = str(data).replace(' ', '').split('*')
        return self(str_lst)


    def _repr_(self, condense=False):
        r"""
        A representation string of this cartesian product of growth groups.

        INPUT:

        - ``condense`` -- (default: ``False``) if set, then a shorter
          output is returned, e.g. the prefix-string ``Growth Group``
          is not show in this case.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(QQ, 'x')
            sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
            sage: cartesian_product([P, L], order='lex')._repr_()
            'Growth Group x^QQ * log(x)^ZZ'
        """
        return GenericGrowthGroup._repr_(self, condense)


    def _repr_short_(self):
        r"""
        A short (shorter than :meth:`._repr_`) representation string
        for this cartesian product of growth groups.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(QQ, 'x')
            sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
            sage: cartesian_product([P, L], order='lex')._repr_short_()
            'x^QQ * log(x)^ZZ'
        """
        return ' * '.join(S._repr_short_() for S in self.cartesian_factors())


    def _convert_to_factor_(self, data):
        r"""
        Helper method. Try to convert some input ``data`` to an
        element of one of the cartesian factors of this product.

        INPUT:

        - ``data`` -- some input to be converted.

        OUTPUT:

        An element of an cartesian factor of this product,
        or ``None``.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * log(x)^QQ * y^QQ')
            sage: e1 = G._convert_to_factor_(x^2)
            sage: (e1, e1.parent())
            (x^2, Growth Group x^ZZ * log(x)^QQ)
            sage: G._convert_to_factor_('asdf') is None
            True
        """
        for factor in self.cartesian_factors():
            try:
                if hasattr(factor, '_convert_to_factor_'):
                    return factor(factor._convert_to_factor_(data))
                return factor(data)
            except (ValueError, TypeError):
                continue


    def _coerce_map_from_(self, S):
        r"""
        Return ``True`` if there is a coercion from ``S``.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ * log(x)^ZZ').has_coerce_map_from(QQ)  # indirect doctest
            False
        """
        pass


    def gens_monomial(self):
        r"""
        Return a tuple containing monomial generators of this growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple containing elements of this growth group.

        .. NOTE::

            This method calls the ``gens_monomial()`` method on the
            individual factors of this cartesian product and
            concatenates the respective outputs.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ * log(x)^ZZ * y^QQ * log(z)^ZZ')
            sage: G.gens_monomial()
            (x, y)
        """
        return sum(iter(factor.gens_monomial()
                        for factor in self.cartesian_factors()),
                   tuple())


    def variable_names(self):
        r"""
        Return the names of the variables.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ * log(x)^ZZ * y^QQ * log(z)^ZZ').variable_names()
            ('x', 'y', 'z')
        """
        vars = sum(iter(factor.variable_names()
                        for factor in self.cartesian_factors()),
                   tuple())
        from itertools import groupby
        return tuple(v for v, _ in groupby(vars))


    class Element(CartesianProductPosets.Element):
        def _repr_(self):
            r"""
            A representation string for this cartesian product element.

            INPUT:

            Nothing.

            OUTPUT:

            A string.

            EXAMPLES::

                sage: import sage.rings.asymptotic.growth_group as agg
                sage: P = agg.MonomialGrowthGroup(QQ, 'x')
                sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
                sage: cartesian_product([P, L], order='lex').an_element()._repr_()
                'x^(1/2) * log(x)'
            """
            s = ' * '.join(repr(v) for v in self.value if not v.is_one())
            if s == '':
                return '1'
            return s


    CartesianProduct = CartesianProductGrowthGroups


class UnivariateProduct(GenericProduct):
    def __init__(self, sets, category, **kwargs):
        r"""
        A cartesian product of growth groups with the same variables.

        TEST::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: type(GrowthGroup('x^ZZ * log(x)^ZZ'))  # indirect doctest
            <class 'sage.rings.asymptotic.growth_group_cartesian.UnivariateProduct_with_category'>
        """
        super(UnivariateProduct, self).__init__(
            sets, category, order='lex', **kwargs)


    CartesianProduct = CartesianProductGrowthGroups


class MultivariateProduct(GenericProduct):
    def __init__(self, sets, category, **kwargs):
        r"""
        A cartesian product of growth groups with the pairwise different variables.

        TEST::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: type(GrowthGroup('x^ZZ * y^ZZ'))  # indirect doctest
            <class 'sage.rings.asymptotic.growth_group_cartesian.MultivariateProduct_with_category'>
        """
        super(MultivariateProduct, self).__init__(
            sets, category, order='product', **kwargs)


    CartesianProduct = CartesianProductGrowthGroups
