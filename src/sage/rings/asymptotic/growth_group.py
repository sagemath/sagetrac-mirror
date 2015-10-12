r"""
(Asymptotic) Growth Groups

This module provides support for (asymptotic) growth groups.

Such groups are equipped with a partial order: the elements can be
seen as functions, and the behavior as their argument (or arguments)
gets large (tend to `\infty`) is compared.

Growth groups are used for the calculations done in the
:doc:`asymptotic ring <asymptotic_ring>`.

A Formal Definition
===================

The elements of a :doc:`growth group <growth_group>` are equipped with a partial
ordering and usually contain a variable. Examples are (among many
other possibilities)

- elements of the form `z^q` for some integer or rational `q` (growth
  groups ``z^ZZ`` or ``z^QQ``),

- elements of the form `\log(z)^q` for some integer or rational `q` (growth
  groups ``log(z)^ZZ`` or ``log(z)^QQ``),

- elements of the form `a^z` for some
  rational `a` (growth group ``QQ^z``), or

- more sophisticated constructions like products `x^r \log(x)^s \cdot
  a^y \cdot y^q` (this corresponds to an element of the growth group
  ``x^QQ * log(x)^ZZ * QQ^y * y^QQ``).

The ordering in all these examples is the growth as `x`, `y`, or `z`
(independently) tend to `\infty`. For elements only using the
variable `z` this means, `g_1 \leq g_2` if

.. MATH::

    \lim_{z\to\infty} \frac{g_2}{g_1} \leq 1.

.. WARNING::

    As this code is experimental, warnings are thrown when a growth
    group is created for the first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup, GrowthGroup
        sage: GenericGrowthGroup(ZZ)
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        Growth Group Generic(ZZ)
        sage: GrowthGroup('x^ZZ * log(x)^ZZ')
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        Growth Group x^ZZ * log(x)^ZZ

Overview
========

For many purposes the factory ``GrowthGroup`` (see
:class:`GrowthGroupFactory`) is the most convenient way to generate a
growth group.

    sage: from sage.rings.asymptotic.growth_group import GrowthGroup

Here are some examples::

    sage: GrowthGroup('z^ZZ')
    Growth Group z^ZZ
    sage: GrowthGroup('z^QQ')
    Growth Group z^QQ

Each of these two generated groups is a :class:`MonomialGrowthGroup`,
whose elements are powers of a fixed symbol (above ``'z'``).

Similarly, we can construct logarithmic factors by::

    sage: GrowthGroup('log(z)^QQ')
    Growth Group log(z)^QQ

which again creates a
:class:`MonomialGrowthGroup`. An :class:`ExponentialGrowthGroup` is generated in the same way. Our factory gives
::

    sage: E = GrowthGroup('QQ^z'); E
    Growth Group QQ^z

and a typical element looks like this::

    sage: E.an_element()
    (1/2)^z

More complex groups are created in a similar fashion. For example
::

    sage: C = GrowthGroup('QQ^z * z^QQ * log(z)^QQ'); C
    Growth Group QQ^z * z^QQ * log(z)^QQ

This contains elements of the form
::

    sage: C.an_element()
    (1/2)^z*z^(1/2)*log(z)^(1/2)

The group `C` itself is a cartesian product; to be precise a
:class:`~sage.rings.asymptotic.growth_group_cartesian.UnivariateProduct`. We
can see its factors::

    sage: C.cartesian_factors()
    (Growth Group QQ^z, Growth Group z^QQ, Growth Group log(z)^QQ)

Multivariate constructions are also possible::

    sage: GrowthGroup('x^QQ * y^QQ')
    Growth Group x^QQ * y^QQ

This gives a
:class:`~sage.rings.asymptotic.growth_group_cartesian.MultivariateProduct`.

Both these cartesian products are derived from the class
:class:`~sage.rings.asymptotic.growth_group_cartesian.GenericProduct`. Moreover
all growth groups have the abstract base class
:class:`GenericGrowthGroup` in common.

Some Examples
^^^^^^^^^^^^^

EXAMPLES::

    sage: import sage.rings.asymptotic.growth_group as agg
    sage: G_x = agg.GrowthGroup('x^ZZ'); repr(G_x)
    'Growth Group x^ZZ'
    sage: G_xy = agg.GrowthGroup('x^ZZ * y^ZZ'); G_xy
    Growth Group x^ZZ * y^ZZ
    sage: G_xy.an_element()
    x*y
    sage: x = G_xy('x'); y = G_xy('y')
    sage: x^2
    x^2
    sage: elem = x^21*y^21; elem^2
    x^42*y^42

A monomial growth group itself is totally ordered, all elements
are comparable. However, this does **not** hold for cartesian
products::

    sage: e1 = x^2*y; e2 = x*y^2
    sage: e1 <= e2 or e2 <= e1
    False

In terms of uniqueness, we have the following behaviour::

    sage: agg.GrowthGroup('x^ZZ * y^ZZ') is agg.GrowthGroup('y^ZZ * x^ZZ')
    True

The above is ``True`` since the order of the factors does not play a role here; they use different variables. But when using the same variable, it plays a role::

    sage: agg.GrowthGroup('x^ZZ * log(x)^ZZ') is agg.GrowthGroup('log(x)^ZZ * x^ZZ')
    False

In this case the components are ordered lexicographically, which
means that in the second growth group, ``log(x)`` is assumed to
grow faster than ``x`` (which is nonsense, mathematically). See
:class:`CartesianProduct <sage.rings.asymptotic.growth_group_cartesian.CartesianProductFactory>`
for more details.

Short notation also allows the construction of more complicated 
growth groups::

    sage: G = agg.GrowthGroup('QQ^x * x^ZZ * log(x)^QQ * y^QQ')
    sage: G.an_element()
    (1/2)^x*x*log(x)^(1/2)*y^(1/2)
    sage: x, y = var('x y')
    sage: G(2^x * log(x) * y^(1/2)) * G(x^(-5) * 5^x * y^(1/3))
    10^x*x^(-5)*log(x)*y^(5/6)

AUTHORS:

- Benjamin Hackl (2015-01-01): initial version
- Daniel Krenn (2015-05-29): initial version and review
- Daniel Krenn (2015-06-02): cartesian products
- Benjamin Hackl (2015-07-00): growth group factory
- Benjamin Hackl (2015-08-00): exponential growth group, initial version
- Daniel Krenn (2015-08-31): various improvements, review; documentation

ACKNOWLEDGEMENT:

- Benjamin Hackl, Clemens Heuberger and Daniel Krenn are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

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

import sage
from sage.misc.lazy_import import lazy_import
lazy_import('sage.rings.asymptotic.growth_group_cartesian', 'CartesianProductGrowthGroups')


def is_lt_one(self):
    r"""
    Return if this element is less than `1`.

    INPUT:

    Nothing.

    OUTPUT:

    A boolean.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G = GrowthGroup('x^ZZ'); x = G(x)
        sage: (x^42).is_lt_one()
        False
        sage: (x^(-42)).is_lt_one()
        True
    """
    one = self.parent().one()
    return self <= one and self != one

class Variable(sage.structure.unique_representation.CachedRepresentation,
               sage.structure.sage_object.SageObject):
    r"""
    A class managing the variable of a growth group.

    INPUT:

    - ``var`` -- an object whose representation string is used as the
      variable. It has to be a valid Python identifier. ``var`` can
      also be a tuple (or other iterable) of such objects.

    - ``repr`` -- (default: ``None``) if specified, then this string
      will be displayed instead of ``var``. Use this to get
      e.g. ``log(x)^ZZ``: ``var`` is then used to specify the variable `x`.


    TESTS::

        sage: from sage.rings.asymptotic.growth_group import Variable
        sage: v = Variable('x'); repr(v), v.variable_names()
        ('x', ('x',))
        sage: v = Variable('x1'); repr(v), v.variable_names()
        ('x1', ('x1',))
        sage: v = Variable('x_42'); repr(v), v.variable_names()
        ('x_42', ('x_42',))
        sage: v = Variable(' x'); repr(v), v.variable_names()
        ('x', ('x',))
        sage: v = Variable('x '); repr(v), v.variable_names()
        ('x', ('x',))
        sage: v = Variable(''); repr(v), v.variable_names()
        ('', ())

    ::

        sage: v = Variable(('x', 'y')); repr(v), v.variable_names()
        ('x, y', ('x', 'y'))
        sage: v = Variable(('x', 'log(y)')); repr(v), v.variable_names()
        ('x, log(y)', ('x', 'y'))
        sage: v = Variable(('x', 'log(x)')); repr(v), v.variable_names()
        Traceback (most recent call last):
        ...
        ValueError: Variable names ('x', 'x') are not pairwise distinct.

    ::

        sage: v = Variable('log(x)'); repr(v), v.variable_names()
        ('log(x)', ('x',))
        sage: v = Variable('log(log(x))'); repr(v), v.variable_names()
        ('log(log(x))', ('x',))

    ::

        sage: v = Variable('x', repr='log(x)'); repr(v), v.variable_names()
        ('log(x)', ('x',))
    """
    def __init__(self, var, repr=None):
        r"""
        See :class:`Variable` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('blub')
            blub
            sage: Variable('blub') is Variable('blub')
            True

        ::

            sage: Variable('(:-)')
            Traceback (most recent call last):
            ...
            TypeError: Malformed expression: (: !!! -)
            sage: Variable('(:-)', repr='icecream')
            Traceback (most recent call last):
            ...
            ValueError: '(:-)' is not a valid name for a variable.
        """
        from sage.symbolic.ring import isidentifier

        if not isinstance(var, (list, tuple)):
            var = (var,)
        var = tuple(str(v).strip() for v in var)

        if repr is None:
            var_bases = sum(iter(
                self.extract_variable_names(v)
                if not isidentifier(v) else (v,)
                for v in var), tuple())
            var_repr = ', '.join(var)
        else:
            for v in var:
                if not isidentifier(v):
                    raise ValueError("'%s' is not a valid name for a variable." % (v,))
            var_bases = var
            var_repr = str(repr).strip()

        if len(var_bases) != len(set(var_bases)):
            raise ValueError('Variable names %s are not pairwise distinct.' %
                             (var_bases,))
        self.var_bases = var_bases
        self.var_repr = var_repr


    def __hash__(self):
        r"""
        Return the hash of this variable.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: hash(Variable('blub'))  # random
            -123456789
        """
        return hash((self.var_repr,) + self.var_bases)


    def __eq__(self, other):
        r"""
        Compare whether this variable equals ``other``.

        INPUT:

        - ``other`` -- another variable.

        OUTPUT:

        A boolean.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x') == Variable('x')
            True
            sage: Variable('x') == Variable('y')
            False
        """
        return self.var_repr == other.var_repr and self.var_bases == other.var_bases


    def __ne__(self, other):
        r"""
        Return whether this variable does not equal ``other``.

        INPUT:

        - ``other`` -- another variable.

        OUTPUT:

        A boolean.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x') != Variable('x')
            False
            sage: Variable('x') != Variable('y')
            True
        """
        return not self == other


    def _repr_(self):
        r"""
        Return a representation string of this variable.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('blub')  # indirect doctest
            blub
        """
        return self.var_repr


    def variable_names(self):
        r"""
        Return the names of the variables.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x').variable_names()
            ('x',)
            sage: Variable('log(x)').variable_names()
            ('x',)
        """
        return self.var_bases


    def is_monomial(self):
        r"""
        Return whether this is a monomial variable.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x').is_monomial()
            True
            sage: Variable('log(x)').is_monomial()
            False
        """
        return len(self.var_bases) == 1 and self.var_bases[0] == self.var_repr


    @staticmethod
    def extract_variable_names(s):
        r"""
        Determine the name of the variable for the given string.

        INPUT:

        - ``s`` -- a string.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable.extract_variable_names('')
            ()
            sage: Variable.extract_variable_names('x')
            ('x',)
            sage: Variable.extract_variable_names('exp(x)')
            ('x',)
            sage: Variable.extract_variable_names('sin(cos(ln(x)))')
            ('x',)

        ::

            sage: Variable.extract_variable_names('log(77w)')
            ('w',)
            sage: Variable.extract_variable_names('log(x')
            Traceback (most recent call last):
            ....
            TypeError: Bad function call: log(x !!!
            sage: Variable.extract_variable_names('x)')
            Traceback (most recent call last):
            ....
            TypeError: Malformed expression: x) !!!
            sage: Variable.extract_variable_names('log)x(')
            Traceback (most recent call last):
            ....
            TypeError: Malformed expression: log) !!! x(
            sage: Variable.extract_variable_names('log(x)+y')
            ('x', 'y')
            sage: Variable.extract_variable_names('icecream(summer)')
            ('summer',)

        ::

            sage: Variable.extract_variable_names('a + b')
            ('a', 'b')
            sage: Variable.extract_variable_names('a+b')
            ('a', 'b')
            sage: Variable.extract_variable_names('a +b')
            ('a', 'b')
            sage: Variable.extract_variable_names('+a')
            ('a',)
            sage: Variable.extract_variable_names('a+')
            Traceback (most recent call last):
            ...
            TypeError: Malformed expression: a+ !!!
            sage: Variable.extract_variable_names('b!')
            ('b',)
            sage: Variable.extract_variable_names('-a')
            ('a',)
            sage: Variable.extract_variable_names('a*b')
            ('a', 'b')
            sage: Variable.extract_variable_names('2^q')
            ('q',)
            sage: Variable.extract_variable_names('77')
            ()

        ::

            sage: Variable.extract_variable_names('a + (b + c) + d')
            ('a', 'b', 'c', 'd')
        """
        from sage.symbolic.ring import SR
        if s == '':
            return ()
        return tuple(str(s) for s in SR(s).variables())


def log(self, base=None):
    r"""
    Return the logarithm of this element.

    INPUT:

    - ``base`` -- the base of the logarithm. If ``None``
      (default value) is used, the natural logarithm is taken.

    OUTPUT:

    A growth element.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G = GrowthGroup('x^ZZ * log(x)^ZZ')
        sage: x, = G.gens_monomial()
        sage: log(x)  # indirect doctest
        log(x)
        sage: log(x^5)
        Traceback (most recent call last):
        ...
        ArithmeticError: When calculating log(x^5) a factor 5 != 1 appeared,
        which is not contained in Growth Group x^ZZ * log(x)^ZZ.

    ::

        sage: G = GrowthGroup('QQ^x * x^ZZ')
        sage: x, = G.gens_monomial()
        sage: el = x.rpow(2); el
        2^x
        sage: log(el)
        Traceback (most recent call last):
        ...
        ArithmeticError: When calculating log(2^x) a factor log(2) != 1
        appeared, which is not contained in Growth Group QQ^x * x^ZZ.
        sage: log(el, base=2)
        x

    ::

        sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
        sage: x = GenericGrowthGroup(ZZ).an_element()
        sage: log(x)  # indirect doctest
        Traceback (most recent call last):
        ...
        NotImplementedError: Cannot determine logarithmized factorization of
        GenericGrowthElement(1) in abstract base class.

    ::

        sage: x = GrowthGroup('x^ZZ').an_element()
        sage: log(x)
        Traceback (most recent call last):
        ...
        ArithmeticError: Cannot build log(x) since log(x) is not in
        Growth Group x^ZZ.

    TESTS::

        sage: G = GrowthGroup("QQ['e']^x * x^ZZ")
        sage: x, = G.gens_monomial()
        sage: log(exp(x))
        x

    ::

        sage: G.one().log()
        Traceback (most recent call last):
        ...
        ArithmeticError: log(1) is zero, which is not contained in
        Growth Group QQ[e]^x * x^ZZ.
    """
    log_factor = self.log_factor(base=base)
    if not log_factor:
        raise ArithmeticError('log(%s) is zero, '
                              'which is not contained in %s.' %
                              (self, self.parent()))

    if len(log_factor) != 1:
        raise ArithmeticError('Calculating log(%s) results in a sum, '
                              'which is not contained in %s.' %
                              (self, self.parent()))
    g, c = log_factor[0]
    if c != 1:
        raise ArithmeticError('When calculating log(%s) a factor %s != 1 '
                              'appeared, which is not contained in %s.' %
                              (self, c, self.parent()))
    return g


def log_factor(self, base=None):
    r"""
    Return the logarithm of the factorization of this
    element.

    INPUT:

    - ``base`` -- the base of the logarithm. If ``None``
      (default value) is used, the natural logarithm is taken.

    OUTPUT:

    A tuple of pairs, where the first entry is a growth
    element and the second a multiplicative coefficient.

    ALGORITHM:

        This function factors the given element and calculates
        the logarithm of each of these factors.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G = GrowthGroup('QQ^x * x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ')
        sage: x, y = G.gens_monomial()
        sage: (x * y).log_factor()
        ((log(x), 1), (log(y), 1))
        sage: (x^123).log_factor()
        ((log(x), 123),)
        sage: (G('2^x') * x^2).log_factor(base=2)
        ((x, 1), (log(x), 2/log(2)))

    ::

        sage: G(1).log_factor()
        ()

    ::

        sage: log(x).log_factor()
        Traceback (most recent call last):
        ...
        ArithmeticError: Cannot build log(log(x)) since log(log(x)) is
        not in Growth Group QQ^x * x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ.

    .. SEEALSO::

        :meth:`~GenericGrowthElement.factors`,
        :meth:`~GenericGrowthElement.log`.

    TESTS::

        sage: G = GrowthGroup("QQ['e']^x * x^ZZ * log(x)^ZZ")
        sage: x, = G.gens_monomial()
        sage: (exp(x) * x).log_factor()
        ((x, 1), (log(x), 1))
    """
    log_factor = self._log_factor_(base=base)

    from growth_group import GenericGrowthGroup
    for g, c in log_factor:
        if hasattr(g, 'parent') and \
           isinstance(g.parent(), GenericGrowthGroup):
            continue
        raise ArithmeticError('Cannot build log(%s) since %s '
                              'is not in %s.' % (self, g, self.parent()))

    return log_factor


def rpow(self, base):
    r"""
    Calculate the power of ``base`` to this element.

    INPUT:

    - ``base`` -- an element.

    OUTPUT:

    A growth element.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G = GrowthGroup('QQ^x * x^ZZ')
        sage: x = G('x')
        sage: x.rpow(2)
        2^x
        sage: x.rpow(1/2)
        (1/2)^x

    ::

        sage: x.rpow(0)
        Traceback (most recent call last):
        ...
        ValueError: 0 is not an allowed base for calculating the power to x.
        sage: (x^2).rpow(2)
        Traceback (most recent call last):
        ...
        ArithmeticError: Cannot construct 2^(x^2) in Growth Group QQ^x * x^ZZ
        > *previous* TypeError: unsupported operand parent(s) for '*':
        'Growth Group QQ^x * x^ZZ' and 'Growth Group ZZ^(x^2)'

    ::

        sage: G = GrowthGroup('QQ^(x*log(x)) * x^ZZ * log(x)^ZZ')
        sage: x = G('x')
        sage: (x * log(x)).rpow(2)
        2^(x*log(x))
    """
    if base == 0:
        raise ValueError('%s is not an allowed base for calculating the '
                         'power to %s.' % (base, self))

    element = None
    try:
        element = self._rpow_element_(base)
    except ValueError:
        pass

    var = str(self)
    if '*' in var or '^' in var:
        var = '(' + var + ')'

    if element is None:
        if base == 'e':
            from sage.rings.integer_ring import ZZ
            base = ZZ[base](base)
        E = ExponentialGrowthGroup(base.parent(), var)
        element = E(raw_element=base)

    try:
        return self.parent().one() * element
    except (TypeError, ValueError) as e:
        from misc import combine_exceptions
        raise combine_exceptions(
            ArithmeticError('Cannot construct %s^%s in %s' %
                            (base, var, self.parent())), e)


class GenericGrowthElement(sage.structure.element.MultiplicativeGroupElement):
    r"""
    A basic implementation of a generic growth element.

    Growth elements form a group by multiplication, and (some of) the
    elements can be compared to each other, i.e., all elements form a
    poset.

    INPUT:

    - ``parent`` -- a :class:`GenericGrowthGroup`.

    - ``raw_element`` -- an element from the base of the parent.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G = agg.GenericGrowthGroup(ZZ)
        sage: g = agg.GenericGrowthElement(G, 42); g
        GenericGrowthElement(42)
        sage: g.parent()
        Growth Group Generic(ZZ)
        sage: G(raw_element=42) == g
        True
    """

    def __init__(self, parent, raw_element):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)
            GenericGrowthElement(42)

        TESTS::

            sage: G(raw_element=42).category()
            Category of elements of Growth Group Generic(ZZ)

        ::

            sage: agg.GenericGrowthElement(None, 0)
            Traceback (most recent call last):
            ...
            ValueError: The parent must be provided
        """
        if parent is None:
            raise ValueError('The parent must be provided')
        super(GenericGrowthElement, self).__init__(parent=parent)

        self._raw_element_ = parent.base()(raw_element)


    def _repr_(self):
        r"""
        A representation string for this generic element.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)  # indirect doctest
            GenericGrowthElement(42)
            sage: H = GenericGrowthGroup(ZZ, 'h')
            sage: H(raw_element=42)  # indirect doctest
            GenericGrowthElement(42, h)
        """
        vars = ', '.join(self.parent()._var_.variable_names())
        if vars:
            vars = ', ' + vars
        return 'GenericGrowthElement(%s%s)' % (self._raw_element_, vars)


    def __hash__(self):
        r"""
        Return the hash of this element.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ);
            sage: hash(G(raw_element=42))  # random
            5656565656565656
        """
        return hash((self.parent(), self._raw_element_))


    def _mul_(self, other):
        r"""
        Abstract multiplication method for generic elements.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`.

        OUTPUT:

        A :class:`GenericGrowthElement` representing the product with
        ``other``.

        .. NOTE::

            Inherited classes must override this.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: g = G.an_element()
            sage: g*g
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
        """
        raise NotImplementedError('Only implemented in concrete realizations.')


    def __invert__(self):
        r"""
        Return the inverse of this growth element.

        OUTPUT:

        An instance of :class:`GenericGrowthElement`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: ~G.an_element()
            Traceback (most recent call last):
            ...
            NotImplementedError: Inversion of GenericGrowthElement(1) not implemented
            (in this abstract method).
            sage: G.an_element()^7
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
            sage: P = GrowthGroup('x^ZZ')
            sage: ~P.an_element()
            x^(-1)
        """
        raise NotImplementedError('Inversion of %s not implemented '
                                  '(in this abstract method).' % (self,))


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

            sage: import sage.rings.asymptotic.growth_group as agg
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

            sage: P_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_QQ = agg.MonomialGrowthGroup(QQ, 'x')
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

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P(raw_element=1)
            sage: e1._eq_(P.gen())
            True
            sage: e2 = e1^4
            sage: e2 == e1^2*e1*e1
            True
            sage: e2 == e1
            False
        """
        return self._raw_element_ == other._raw_element_


    def __ne__(self, other):
        r"""
        Return if this growth element is not equal to ``other``.

        INPUT:

        - ``other`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

            The comparison of two elements with the same parent is done in
            :meth:`_eq_`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: G.one() != G(1)
            False
            sage: G.one() != G.one()
            False
            sage: G(1) != G(1)
            False
        """
        return not self == other


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

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_QQ = agg.MonomialGrowthGroup(QQ, 'x')
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

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: e1 = G(raw_element=1); e2 = G(raw_element=2)
            sage: e1 <= e2  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
        """
        raise NotImplementedError('Only implemented in concrete realizations.')


    log = log
    log_factor = log_factor


    def _log_factor_(self, base=None):
        r"""
        Helper method for calculating the logarithm of the factorization
        of this element.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        OUTPUT:

        A tuple of pairs, where the first entry is either a growth
        element or something out of which we can construct a growth element
        and the second a multiplicative coefficient.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(QQ)
            sage: G.an_element().log_factor()  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot determine logarithmized factorization of
            GenericGrowthElement(1/2) in abstract base class.
        """
        raise NotImplementedError('Cannot determine logarithmized factorization '
                                  'of %s in abstract base class.' % (self,))



    rpow = rpow


    def _rpow_element_(self, base):
        r"""
        Return an element which is the power of ``base`` to this
        element; it lives (in contrast to :meth:`rpow`) in its own group.

        INPUT:

        - ``base`` -- an element.

        OUTPUT:

        A growth element or ``None``.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('QQ^x')
            sage: x = G(raw_element=3)
            sage: x._rpow_element_(2) is None
            True
        """
        pass


    def factors(self):
        r"""
        Return the atomic factors of this growth element. An atomic factor
        cannot be split further.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple of growth elements.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: G.an_element().factors()
            (x,)
        """
        return (self,)

    is_lt_one = is_lt_one


class GenericGrowthGroup(
        sage.structure.unique_representation.UniqueRepresentation,
        sage.structure.parent.Parent):
    r"""
    A basic implementation for growth groups.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``).

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    .. NOTE::

        This class should be derived for concrete implementations.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G = agg.GenericGrowthGroup(ZZ); G
        Growth Group Generic(ZZ)

    .. SEEALSO::

        :class:`MonomialGrowthGroup`
    """
    # TODO: implement some sort of 'assume', where basic assumptions
    # for the variables can be stored. --> within the cartesian product

    # enable the category framework for elements
    Element = GenericGrowthElement


    @staticmethod
    def __classcall__(cls, base, var=None, category=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        For more information see :class:`MonomialGrowthGroup`.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P1 = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P2 = agg.MonomialGrowthGroup(ZZ, ZZ['x'].gen())
            sage: P3 = agg.MonomialGrowthGroup(ZZ, SR.var('x'))
            sage: P1 is P2 and P2 is P3
            True
            sage: P4 = agg.MonomialGrowthGroup(ZZ, buffer('xylophone', 0, 1))
            sage: P1 is P4
            True
            sage: P5 = agg.MonomialGrowthGroup(ZZ, 'x ')
            sage: P1 is P5
            True

        ::

            sage: L1 = agg.MonomialGrowthGroup(QQ, log(x))
            sage: L2 = agg.MonomialGrowthGroup(QQ, 'log(x)')
            sage: L1 is L2
            True
        """
        if not isinstance(base, sage.structure.parent.Parent):
            raise TypeError('%s is not a valid base.' % (base,))

        if var is None:
            var = Variable('')
        elif not isinstance(var, Variable):
            var = Variable(var)

        if category is None:
            from sage.categories.monoids import Monoids
            from sage.categories.posets import Posets
            category = Monoids() & Posets()

        return super(GenericGrowthGroup, cls).__classcall__(
            cls, base, var, category)


    @sage.misc.superseded.experimental(trac_number=17601)
    def __init__(self, base, var, category):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ).category()
            Join of Category of monoids and Category of posets

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x')
            Growth Group x^ZZ
            sage: agg.MonomialGrowthGroup(QQ, SR.var('n'))
            Growth Group n^QQ
            sage: agg.MonomialGrowthGroup(ZZ, ZZ['y'].gen())
            Growth Group y^ZZ
            sage: agg.MonomialGrowthGroup(QQ, 'log(x)')
            Growth Group log(x)^QQ

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.ExponentialGrowthGroup(QQ, 'x')
            Growth Group QQ^x
            sage: agg.ExponentialGrowthGroup(SR, ZZ['y'].gen())
            Growth Group SR^y

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G.is_parent_of(G(raw_element=42))
            True
            sage: G2 = agg.GenericGrowthGroup(ZZ, category=FiniteGroups() & Posets())
            sage: G2.category()
            Join of Category of finite groups and Category of finite posets

        ::

            sage: G = agg.GenericGrowthGroup('42')
            Traceback (most recent call last):
            ...
            TypeError: 42 is not a valid base.

        ::

            sage: agg.MonomialGrowthGroup('x', ZZ)
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base.
            sage: agg.MonomialGrowthGroup('x', 'y')
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base.

        ::

            sage: agg.ExponentialGrowthGroup('x', ZZ)
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base.
            sage: agg.ExponentialGrowthGroup('x', 'y')
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base.

        """
        self._var_ = var
        super(GenericGrowthGroup, self).__init__(category=category,
                                                 base=base)


    def _repr_short_(self):
        r"""
        A short representation string of this abstract growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: GenericGrowthGroup(QQ)._repr_short_()
            'Generic(QQ)'
            sage: GenericGrowthGroup(QQ)
            Growth Group Generic(QQ)
            sage: GenericGrowthGroup(QQ, ('a', 'b'))
            Growth Group Generic(QQ, a, b)
        """
        from misc import parent_to_repr_short
        vars = ', '.join(self._var_.variable_names())
        if vars:
            vars = ', ' + vars
        return 'Generic(%s%s)' % (parent_to_repr_short(self.base()), vars)


    def _repr_(self, condense=False):
        r"""
        A representations string of this growth group.

        INPUT:

        - ``condense`` -- (default: ``False``) if set, then a shorter
          output is returned, e.g. the prefix-string ``Growth Group``
          is not show in this case.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x')  # indirect doctest
            Growth Group x^ZZ
            sage: agg.MonomialGrowthGroup(QQ, 'log(x)')  # indirect doctest
            Growth Group log(x)^QQ

        TESTS::

            sage: agg.MonomialGrowthGroup(QQ, 'log(x)')._repr_(condense=True)
            'log(x)^QQ'
        """
        pre = 'Growth Group ' if not condense else ''
        return '%s%s' % (pre, self._repr_short_())


    def __hash__(self):
        r"""
        Return the hash of this group.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: hash(agg.GenericGrowthGroup(ZZ))  # random
            4242424242424242

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: hash(P)  # random
            -1234567890123456789

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.ExponentialGrowthGroup(ZZ, 'x')
            sage: hash(P)  # random
            -1234567890123456789
        """
        return hash((self.__class__, self.base(), self._var_))


    def _an_element_(self):
        r"""
        Return an element of ``self``.

        INPUT:

        Nothing.

        OUTPUT:

        An element of ``self``.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ).an_element()  # indirect doctest
            GenericGrowthElement(1)
            sage: agg.MonomialGrowthGroup(ZZ, 'z').an_element()  # indirect doctest
            z
            sage: agg.MonomialGrowthGroup(QQ, 'log(z)').an_element()  # indirect doctest
            log(z)^(1/2)
        """
        return self.element_class(self, self.base().an_element())


    def some_elements(self):
        r"""
        Return some elements of this growth group.

        See :class:`TestSuite` for a typical use case.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: tuple(agg.MonomialGrowthGroup(ZZ, 'z').some_elements())
            (1, z, z^(-1), z^2, z^(-2), z^3, z^(-3),
             z^4, z^(-4), z^5, z^(-5), ...)
            sage: tuple(agg.MonomialGrowthGroup(QQ, 'z').some_elements())
            (z^(1/2), z^(-1/2), z^2, z^(-2),
             1, z, z^(-1), z^42,
             z^(2/3), z^(-2/3), z^(3/2), z^(-3/2),
             z^(4/5), z^(-4/5), z^(5/4), z^(-5/4), ...)
        """
        return iter(self.element_class(self, e)
                    for e in self.base().some_elements())


    def _create_element_via_parent_(self, raw_element):
        r"""
        Create an element with a possibly other parent.

        INPUT:

        - ``raw_element`` -- the element data.

        OUTPUT:

        An element.

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('z^ZZ')
            sage: G._create_element_via_parent_(3).parent()
            Growth Group z^ZZ
            sage: G._create_element_via_parent_(1/2).parent()
            Growth Group z^QQ
        """
        if raw_element.parent() is self.base():
            parent = self
        else:
            from misc import underlying_class
            parent = underlying_class(self)(raw_element.parent(), self._var_,
                                            category=self.category())
        return parent(raw_element=raw_element)


    def le(self, left, right):
        r"""
        Return if the growth of ``left`` is at most (less than or
        equal to) the growth of ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: x = G.gen()
            sage: G.le(x, x^2)
            True
            sage: G.le(x^2, x)
            False
            sage: G.le(x^0, 1)
            True
        """
        return self(left) <= self(right)


    def _element_constructor_(self, data, raw_element=None):
        r"""
        Convert a given object to this growth group.

        INPUT:

        - ``data`` -- an object representing the element to be
          initialized.

        - ``raw_element`` -- (default: ``None``) if given, then this is
          directly passed to the element constructor (i.e., no conversion
          is performed).

        OUTPUT:

        An element of this growth group.

        .. NOTE::

            Either ``data`` or ``raw_element`` has to be given. If
            ``raw_element`` is specified, then no positional argument
            may be passed.

            This method calls :meth:`_convert_`, which does the actual
            conversion from ``data``.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_ZZ = agg.GenericGrowthGroup(ZZ)
            sage: z = G_ZZ(raw_element=42); z  # indirect doctest
            GenericGrowthElement(42)
            sage: z is G_ZZ(z)  # indirect doctest
            True

        ::

            sage: G_QQ = agg.GenericGrowthGroup(QQ)
            sage: q = G_QQ(raw_element=42)  # indirect doctest
            sage: q is z
            False
            sage: G_ZZ(q)  # indirect doctest
            GenericGrowthElement(42)
            sage: G_QQ(z)  # indirect doctest
            GenericGrowthElement(42)
            sage: q is G_ZZ(q)  # indirect doctest
            False

        ::

            sage: G_ZZ()
            Traceback (most recent call last):
            ...
            ValueError: No input specified. Cannot continue.
            sage: G_ZZ('blub')  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: blub is not in Growth Group Generic(ZZ).
            sage: G_ZZ('x', raw_element=42)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Input is ambigous: x as well as raw_element=42 are specified.

        ::

            sage: G_x = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: x = G_x(raw_element=1)  # indirect doctest
            sage: G_y = agg.MonomialGrowthGroup(ZZ, 'y')
            sage: G_y(x)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: x is not in Growth Group y^ZZ.
        """
        if raw_element is None:
            if isinstance(data, self.element_class):
                if data.parent() == self:
                    return data
                try:
                    if self._var_ != data.parent()._var_:
                        raise ValueError('%s is not in %s.' % (data, self))
                except AttributeError:
                    pass
                raw_element = data._raw_element_
            elif isinstance(data, int) and data == 0:
                raise ValueError('No input specified. Cannot continue.')
            else:
                raw_element = self._convert_(data)
            if raw_element is None:
                raise ValueError('%s is not in %s.' % (data, self))
        elif not isinstance(data, int) or data != 0:
            raise ValueError('Input is ambigous: '
                             '%s as well as raw_element=%s '
                             'are specified.' % (data, raw_element))

        return self.element_class(self, raw_element)


    def _convert_(self, data):
        r"""
        Convert ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        .. NOTE::

            This method always returns ``None`` in this abstract base
            class, and should be overridden in inherited class.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
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

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: G_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: G_ZZ.has_coerce_map_from(G_QQ)  # indirect doctest
            False
            sage: G_QQ.has_coerce_map_from(G_ZZ)  # indirect doctest
            True

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P_x_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_x_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P_x_ZZ.has_coerce_map_from(P_x_QQ)  # indirect doctest
            False
            sage: P_x_QQ.has_coerce_map_from(P_x_ZZ)  # indirect doctest
            True
            sage: P_y_ZZ = agg.MonomialGrowthGroup(ZZ, 'y')
            sage: P_y_ZZ.has_coerce_map_from(P_x_ZZ)  # indirect doctest
            False
            sage: P_x_ZZ.has_coerce_map_from(P_y_ZZ)  # indirect doctest
            False
            sage: P_y_ZZ.has_coerce_map_from(P_x_QQ)  # indirect doctest
            False
            sage: P_x_QQ.has_coerce_map_from(P_y_ZZ)  # indirect doctest
            False

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P_x_ZZ = agg.GrowthGroup('ZZ^x')
            sage: P_x_QQ = agg.GrowthGroup('QQ^x')
            sage: P_x_ZZ.has_coerce_map_from(P_x_QQ)  # indirect doctest
            False
            sage: P_x_QQ.has_coerce_map_from(P_x_ZZ)  # indirect doctest
            True
            sage: P_y_ZZ = agg.GrowthGroup('ZZ^y')
            sage: P_y_ZZ.has_coerce_map_from(P_x_ZZ)  # indirect doctest
            False
            sage: P_x_ZZ.has_coerce_map_from(P_y_ZZ)  # indirect doctest
            False
            sage: P_y_ZZ.has_coerce_map_from(P_x_QQ)  # indirect doctest
            False
            sage: P_x_QQ.has_coerce_map_from(P_y_ZZ)  # indirect doctest
            False

        ::

            sage: agg.GrowthGroup('x^QQ').has_coerce_map_from(agg.GrowthGroup('QQ^x')) # indirect doctest
            False
        """
        if isinstance(S, type(self)) and self._var_ == S._var_:
            if self.base().has_coerce_map_from(S.base()):
                return True


    def _pushout_(self, other):
        r"""
        Construct the pushout of this and the other growth group. This is called by
        :func:`sage.categories.pushout.pushout`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.categories.pushout import pushout
            sage: cm = sage.structure.element.get_coercion_model()
            sage: A = GrowthGroup('QQ^x')
            sage: B = GrowthGroup('y^ZZ')

        When using growth groups with disjoint variable lists, then a
        pushout can be constructed::

            sage: A._pushout_(B)
            Growth Group QQ^x * y^ZZ
            sage: cm.common_parent(A, B)
            Growth Group QQ^x * y^ZZ

        In general, growth groups of the same variable cannot be
        combined automatically, since there is no order relation between the two factors::

            sage: C = GrowthGroup('x^QQ')
            sage: cm.common_parent(A, C)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Growth Group QQ^x' and 'Growth Group x^QQ'

        However, combining is possible if the factors with the same variable
        overlap::

            sage: cm.common_parent(GrowthGroup('x^ZZ * log(x)^ZZ'), GrowthGroup('exp(x)^ZZ * x^ZZ'))
            Growth Group exp(x)^ZZ * x^ZZ * log(x)^ZZ
            sage: cm.common_parent(GrowthGroup('x^ZZ * log(x)^ZZ'), GrowthGroup('y^ZZ * x^ZZ'))
            Growth Group x^ZZ * log(x)^ZZ * y^ZZ

        ::

            sage: cm.common_parent(GrowthGroup('x^ZZ'), GrowthGroup('y^ZZ'))
            Growth Group x^ZZ * y^ZZ

        ::

            sage: cm.record_exceptions()
            sage: cm.common_parent(GrowthGroup('x^ZZ'), GrowthGroup('y^ZZ'))
            Growth Group x^ZZ * y^ZZ
            sage: sage.structure.element.coercion_traceback()  # not tested
        """
        if not isinstance(other, GenericGrowthGroup) and \
           not (other.construction() is not None and
                isinstance(other.construction()[0], AbstractGrowthGroupFunctor)):
            return

        if set(self.variable_names()).isdisjoint(set(other.variable_names())):
            from sage.categories.cartesian_product import cartesian_product
            return cartesian_product([self, other])


    def gens_monomial(self):
        r"""
        Return a tuple containing monomial generators of this growth
        group.

        INPUT:

        Nothing.

        OUTPUT:

        An empty tuple.

        .. NOTE::

            A generator is called monomial generator if the variable
            of the underlying growth group is a valid identifier. For
            example, ``x^ZZ`` has ``x`` as a monomial generator,
            while ``log(x)^ZZ`` or ``icecream(x)^ZZ`` do not have
            monomial generators.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ).gens_monomial()
            ()
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GrowthGroup('ZZ^x').gens_monomial()
            ()
        """
        return tuple()


    def gens(self):
        r"""
        Return a tuple of all generators of this monomial growth
        group.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple whose entries are instances of
        :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.gens()
            (x,)
            sage: agg.MonomialGrowthGroup(ZZ, 'log(x)').gens()
            (log(x),)
        """
        return (self(raw_element=self.base().one()),)


    def gen(self, n=0):
        r"""
        Return the `n`-th generator (as a group) of this growth group.

        INPUT:

        - ``n`` -- default: `0`.

        OUTPUT:

        A :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.gen()
            x

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('QQ^x')
            sage: P.gen()
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        return self.gens()[n]


    def ngens(self):
        r"""
        Return the number of generators (as a group) of this growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A Python integer.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.ngens()
            1
            sage: agg.MonomialGrowthGroup(ZZ, 'log(x)').ngens()
            1

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('QQ^x')
            sage: P.ngens()
            0
        """
        return len(self.gens())


    def variable_names(self):
        r"""
        Return the names of the variables.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ).variable_names()
            ()

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ').variable_names()
            ('x',)
            sage: GrowthGroup('log(x)^ZZ').variable_names()
            ('x',)

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('QQ^x').variable_names()
            ('x',)
            sage: GrowthGroup('QQ^(x*log(x))').variable_names()
            ('x',)
        """
        return self._var_.variable_names()


    CartesianProduct = CartesianProductGrowthGroups


from sage.categories.pushout import ConstructionFunctor
class AbstractGrowthGroupFunctor(ConstructionFunctor):
    r"""
    A base class for the functors constructing growth groups.

    INPUT:

    - ``var`` -- a string or list of strings (or anything else
      :class:`Variable` accepts).

    - ``domain`` -- a category.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: GrowthGroup('z^QQ').construction()[0]  # indirect doctest
        MonomialGrowthGroup[z]

    .. SEEALSO::

        :doc:`asymptotic_ring`,
        :class:`ExponentialGrowthGroupFunctor`,
        :class:`MonomialGrowthGroupFunctor`,
        :class:`sage.rings.asymptotic.asymptotic_ring.AsymptoticRingFunctor`,
        :class:`sage.categories.pushout.ConstructionFunctor`.
    """

    _functor_name = 'AbstractGrowthGroup'

    rank = 13

    def __init__(self, var, domain):
        r"""
        See :class:`AbstractGrowthGroupFunctor` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import AbstractGrowthGroupFunctor
            sage: AbstractGrowthGroupFunctor('x', Groups())
            AbstractGrowthGroup[x]
        """
        if var is None:
            var = Variable('')
        elif not isinstance(var, Variable):
            var = Variable(var)
        self.var = var
        super(ConstructionFunctor, self).__init__(
            domain, sage.categories.monoids.Monoids() & sage.categories.posets.Posets())


    def _repr_(self):
        r"""
        Return a representation string of this functor.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('QQ^t').construction()[0]  # indirect doctest
            ExponentialGrowthGroup[t]
        """
        return '%s[%s]' % (self._functor_name, self.var)


    def merge(self, other):
        r"""
        Merge this functor with ``other`` of possible.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A functor or ``None``.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F = GrowthGroup('QQ^t').construction()[0]
            sage: G = GrowthGroup('t^QQ').construction()[0]
            sage: F.merge(F)
            ExponentialGrowthGroup[t]
            sage: F.merge(G) is None
            True
        """
        if self == other:
            return self


    def __eq__(self, other):
        r"""
        Return whether this functor is equal to ``other``.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F = GrowthGroup('QQ^t').construction()[0]
            sage: G = GrowthGroup('t^QQ').construction()[0]
            sage: F == F
            True
            sage: F == G
            False
        """
        return type(self) == type(other) and self.var == other.var


    def __ne__(self, other):
        r"""
        Return whether this functor is not equal to ``other``.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F = GrowthGroup('QQ^t').construction()[0]
            sage: G = GrowthGroup('t^QQ').construction()[0]
            sage: F != F
            False
            sage: F != G
            True
        """
        return not self == other


class MonomialGrowthElement(GenericGrowthElement):
    r"""
    An implementation of monomial growth elements.

    INPUT:

    - ``parent`` -- a :class:`MonomialGrowthGroup`.

    - ``raw_element`` -- an element from the base ring of the parent.

      This ``raw_element`` is the exponent of the created monomial
      growth element.

    A monomial growth element represents a term of the type
    `\operatorname{variable}^{\operatorname{exponent}}`. The multiplication
    corresponds to the addition of the exponents.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
        sage: e1 = P(1); e1
        1
        sage: e2 = P(raw_element=2); e2
        x^2
        sage: e1 == e2
        False
        sage: P.le(e1, e2)
        True
        sage: P.le(e1, P.gen()) and P.le(P.gen(), e2)
        True
    """

    @property
    def exponent(self):
        r"""
        The exponent of this growth element.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P(x^42).exponent
            42
        """
        return self._raw_element_


    def _repr_(self):
        r"""
        A representation string for this monomial growth element.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P(1)._repr_()
            '1'
            sage: P(x^5)  # indirect doctest
            x^5
            sage: P(x^(1/2))  # indirect doctest
            x^(1/2)

        TESTS::

            sage: P(x^-1)  # indirect doctest
            x^(-1)
            sage: P(x^-42)  # indirect doctest
            x^(-42)
        """
        from sage.rings.integer_ring import ZZ

        var = repr(self.parent()._var_)
        if self.exponent.is_zero():
            return '1'
        elif self.exponent.is_one():
            return var
        elif self.exponent in ZZ and self.exponent > 0:
            return var + '^' + str(self.exponent)
        else:
            return var + '^(' + str(self.exponent) + ')'


    def _mul_(self, other):
        r"""
        Multiply this monomial growth element with another.

        INPUT:

        - ``other`` -- a :class:`MonomialGrowthElement`

        OUTPUT:

        The product as a :class:`MonomialGrowthElement`.

        .. NOTE::

            Two monomial growth elements are multiplied by adding
            their exponents.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: a = P(x^2)
            sage: b = P(x^3)
            sage: c = a._mul_(b); c
            x^5
            sage: c == a*b
            True
            sage: a*b*a  # indirect doctest
            x^7
        """
        return self.parent()(raw_element=self.exponent + other.exponent)


    def __invert__(self):
        r"""
        Return the multiplicative inverse of this monomial growth element.

        INPUT:

        Nothing.

        OUTPUT:

        The multiplicative inverse as a :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P(raw_element=2)
            sage: e2 = e1.__invert__(); e2
            x^(-2)
            sage: e2 == ~e1
            True
        """
        return self.parent()._create_element_via_parent_(-self.exponent)


    def __pow__(self, exponent):
        r"""
        Calculate the power of this growth element to the given ``exponent``.

        INPUT:

        - ``exponent`` -- a number. This can be anything that is a
          valid right hand side of ``*`` with elements of the
          parent's base.

        OUTPUT:

        The result of this exponentiation, a :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('x^ZZ')
            sage: x = P.gen()
            sage: a = x^7; a
            x^7
            sage: a^(1/2)
            x^(7/2)
            sage: (a^(1/2)).parent()
            Growth Group x^QQ
            sage: P = GrowthGroup('x^QQ')
            sage: b = P.gen()^(7/2); b
            x^(7/2)
            sage: b^12
            x^42
        """
        return self.parent()._create_element_via_parent_(self.exponent * exponent)


    def _log_factor_(self, base=None):
        r"""
        Helper method for calculating the logarithm of the factorization
        of this element.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        OUTPUT:

        A tuple of pairs, where the first entry is either a growth
        element or something out of which we can construct a growth element
        and the second a multiplicative coefficient.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^QQ')
            sage: G('x').log_factor()  # indirect doctest
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot build log(x) since log(x) is not in
            Growth Group x^QQ.
        """
        if self.is_one():
            return tuple()
        coefficient = self.exponent
        if base is not None:
            from sage.functions.log import log
            coefficient = coefficient / log(base)
        return (('log(%s)' % (self.parent()._var_,), coefficient),)


    def _rpow_element_(self, base):
        r"""
        Return an element which is the power of ``base`` to this
        element; it lives (in contrast to :meth:`rpow`) in its own group.

        INPUT:

        - ``base`` -- an element.

        OUTPUT:

        A growth element or ``None``.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: x = G('x')
            sage: x._rpow_element_(2)
            Traceback (most recent call last):
            ...
            ValueError: Variable %s is not a log of something.
            sage: G = GrowthGroup('log(x)^ZZ')
            sage: lx = G(raw_element=1); lx
            log(x)
            sage: rp = lx._rpow_element_('e'); rp
            x
            sage: rp.parent()
            Growth Group x^ZZ
        """
        var = str(self.parent()._var_)
        if not(var.startswith('log(') and self.exponent.is_one()):
            raise ValueError('Variable %s is not a log of something.')
        new_var = var[4:-1]
        if base == 'e':
            from sage.rings.integer_ring import ZZ
            M = MonomialGrowthGroup(ZZ, new_var)
            return M(raw_element=ZZ(1))
        else:
            from sage.functions.log import log
            new_exponent = log(base)
            M = MonomialGrowthGroup(new_exponent.parent(), new_var)
            return M(raw_element=new_exponent)


    def _le_(self, other):
        r"""
        Return if this :class:`MonomialGrowthElement` is at most
        (less than or equal to) ``other``.

        INPUT:

        - ``other`` -- a :class:`MonomialGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`MonomialGrowthElement`.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P_ZZ.gen() <= P_QQ.gen()^2  # indirect doctest
            True
        """
        return self.exponent <= other.exponent


class MonomialGrowthGroup(GenericGrowthGroup):
    r"""
    A growth group dealing with powers of a fixed object/symbol.

    The elements :class:`MonomialGrowthElement` of this group represent powers
    of a fixed base; the group law is the multiplication, which corresponds
    to the addition of the exponents of the monomials.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``).

      As monomials are represented by this group, the elements in
      ``base`` are the exponents of these monomials.

    - ``var`` -- an object.

      The string representation of ``var`` acts as a base of the
      monomials represented by this group.

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: P = agg.MonomialGrowthGroup(ZZ, 'x'); P
        Growth Group x^ZZ
        sage: agg.MonomialGrowthGroup(ZZ, log(SR.var('y')))
        Growth Group log(y)^ZZ

    .. SEEALSO::

        :class:`GenericGrowthGroup`

    TESTS::

        sage: L1 = agg.MonomialGrowthGroup(QQ, log(x))
        sage: L2 = agg.MonomialGrowthGroup(QQ, 'log(x)')
        sage: L1 is L2
        True
    """

    # enable the category framework for elements
    Element = MonomialGrowthElement


    def _repr_short_(self):
        r"""
        A short representation string of this monomial growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'a')  # indirect doctest
            Growth Group a^ZZ


        TESTS::

            sage: agg.MonomialGrowthGroup(ZZ, 'a')._repr_short_()
            'a^ZZ'
            sage: agg.MonomialGrowthGroup(QQ, 'a')._repr_short_()
            'a^QQ'
            sage: agg.MonomialGrowthGroup(PolynomialRing(QQ, 'x'), 'a')._repr_short_()
            'a^QQ[x]'
        """
        from misc import parent_to_repr_short
        return '%s^%s' % (self._var_, parent_to_repr_short(self.base()))


    def _convert_(self, data):
        r"""
        Convert ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
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
            x^(-333)
            sage: P(log(x)^2)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: log(x)^2 is not in Growth Group x^ZZ.

        ::

            sage: PR.<x> = ZZ[]; x.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: P(x^2)  # indirect doctest
            x^2

        ::

            sage: PSR.<x> = ZZ[[]]
            sage: P(x^42)  # indirect doctest
            x^42
            sage: P(x^12 + O(x^17))
            Traceback (most recent call last):
            ...
            ValueError: x^12 + O(x^17) is not in Growth Group x^ZZ.

        ::

            sage: R.<w,x> = ZZ[]
            sage: P(x^4242)  # indirect doctest
            x^4242
            sage: P(w^4242)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: w^4242 is not in Growth Group x^ZZ.

        ::

            sage: PSR.<w,x> = ZZ[[]]
            sage: P(x^7)  # indirect doctest
            x^7
            sage: P(w^7)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: w^7 is not in Growth Group x^ZZ.

        ::

            sage: P('x^7')
            x^7
            sage: P('1/x')
            x^(-1)
            sage: P('x^(-2)')
            x^(-2)
            sage: P('x^-2')
            x^(-2)

        ::

            sage: P('1')
            1
        """
        if data == 1 or data == '1':
            return self.base().zero()
        var = repr(self._var_)
        if str(data) == var:
            return self.base().one()

        try:
            P = data.parent()
        except AttributeError:
            if var not in str(data):
                return  # this has to end here
            from sage.symbolic.ring import SR
            return self._convert_(SR(data))

        from sage.symbolic.ring import SR
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        from sage.rings.polynomial.multi_polynomial_ring_generic import \
            MPolynomialRing_generic
        from sage.rings.power_series_ring import PowerSeriesRing_generic
        import operator
        if P is SR:
            if data.operator() == operator.pow:
                base, exponent = data.operands()
                if str(base) == var:
                    return exponent
        elif isinstance(P, (PolynomialRing_general, MPolynomialRing_generic)):
            if data.is_monomial() and len(data.variables()) == 1:
                if var == str(data.variables()[0]):
                    return data.degree()
        elif isinstance(P, PowerSeriesRing_generic):
            if hasattr(data, 'variables') and len(data.variables()) == 1:
                from sage.rings.integer_ring import ZZ
                if data.is_monomial() and data.precision_absolute() not in ZZ:
                    if var == str(data.variables()[0]):
                        return data.degree()
            elif len(P.variable_names()) == 1 and \
                            var == str(data.variable()[0]):
                from sage.rings.integer_ring import ZZ
                if data.is_monomial() and data.precision_absolute() not in ZZ:
                    return data.degree()


    def gens_monomial(self):
        r"""
        Return a tuple containing monomial generators of this growth
        group.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple containing elements of this growth group.

        .. NOTE::

            A generator is called monomial generator if the variable
            of the underlying growth group is a valid identifier. For
            example, ``x^ZZ`` has ``x`` as a monomial generator,
            while ``log(x)^ZZ`` or ``icecream(x)^ZZ`` do not have
            monomial generators.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x').gens_monomial()
            (x,)
            sage: agg.MonomialGrowthGroup(QQ, 'log(x)').gens_monomial()
            ()
        """
        if not self._var_.is_monomial():
            return tuple()
        return (self(raw_element=self.base().one()),)


    def construction(self):
        r"""
        Return the construction of this growth group.

        OUTPUT:

        A pair whose first entry is a
        :class:`monomial construction functor <MonomialGrowthGroupFunctor>`
        and its second entry the base.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ').construction()
            (MonomialGrowthGroup[x], Integer Ring)
        """
        return MonomialGrowthGroupFunctor(self._var_), self.base()


class MonomialGrowthGroupFunctor(AbstractGrowthGroupFunctor):
    r"""
    A :class:`construction functor <sage.categories.pushout.ConstructionFunctor>`
    for :class:`monomial growth groups <MonomialGrowthGroup>`.

    INPUT:

    - ``var`` -- a string or list of strings (or anything else
      :class:`Variable` accepts).

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup, MonomialGrowthGroupFunctor
        sage: GrowthGroup('z^QQ').construction()[0]
        MonomialGrowthGroup[z]

    .. SEEALSO::

        :doc:`asymptotic_ring`,
        :class:`AbstractGrowthGroupFunctor`,
        :class:`ExponentialGrowthGroupFunctor`,
        :class:`sage.rings.asymptotic.asymptotic_ring.AsymptoticRingFunctor`,
        :class:`sage.categories.pushout.ConstructionFunctor`.

    TESTS::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup, MonomialGrowthGroupFunctor
        sage: cm = sage.structure.element.get_coercion_model()
        sage: A = GrowthGroup('x^QQ')
        sage: B = MonomialGrowthGroupFunctor('x')(ZZ['t'])
        sage: cm.common_parent(A, B)
        Growth Group x^QQ[t]
    """

    _functor_name = 'MonomialGrowthGroup'


    def __init__(self, var):
        r"""
        See :class:`MonomialGrowthGroupFunctor` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import MonomialGrowthGroupFunctor
            sage: MonomialGrowthGroupFunctor('x')
            MonomialGrowthGroup[x]
        """
        super(MonomialGrowthGroupFunctor, self).__init__(var,
            sage.categories.commutative_additive_monoids.CommutativeAdditiveMonoids())


    def _apply_functor(self, base):
        r"""
        Apply this functor to the given ``base``.

        INPUT:

        - ``base`` - anything :class:`MonomialGrowthGroup` accepts.

        OUTPUT:

        A monomial growth group.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F, R = GrowthGroup('z^QQ').construction()
            sage: F(R)  # indirect doctest
            Growth Group z^QQ
        """
        return MonomialGrowthGroup(base, self.var)


class ExponentialGrowthElement(GenericGrowthElement):
    r"""
    An implementation of exponential growth elements.

    INPUT:

    - ``parent`` -- an :class:`ExponentialGrowthGroup`.

    - ``raw_element`` -- an element from the base ring of the parent.

      This ``raw_element`` is the base of the created exponential
      growth element.

    An exponential growth element represents a term of the type
    `\operatorname{base}^{\operatorname{variable}}`. The multiplication
    corresponds to the multiplication of the bases.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: P = agg.GrowthGroup('ZZ^x')
        sage: e1 = P(1); e1
        1
        sage: e2 = P(raw_element=2); e2
        2^x
        sage: e1 == e2
        False
        sage: P.le(e1, e2)
        True
        sage: P.le(e1, P(1)) and P.le(P(1), e2)
        True
    """

    @property
    def base(self):
        r"""
        The base of this exponential growth element.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('ZZ^x')
            sage: P(42^x).base
            42
        """
        return self._raw_element_


    def _repr_(self):
        r"""
        A representation string for this exponential growth element.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('QQ^x')
            sage: P(1)._repr_()
            '1'
            sage: P(5^x)  # indirect doctest
            5^x
            sage: P((1/2)^x)  # indirect doctest
            (1/2)^x

        TESTS::

            sage: P((-1)^x)  # indirect doctest
            (-1)^x
        """
        from sage.rings.integer_ring import ZZ

        var = repr(self.parent()._var_)
        if self.base.is_one():
            return '1'
        elif not any(s in str(self.base) for s in '-/^'):
            return str(self.base) + '^' + var
        else:
            return '(' + str(self.base) + ')^' + var


    def _mul_(self, other):
        r"""
        Multiply this exponential growth element with another.

        INPUT:

        - ``other`` -- a :class:`ExponentialGrowthElement`

        OUTPUT:

        The product as a :class:`ExponentialGrowthElement`.

        .. NOTE::

            Two exponential growth elements are multiplied by
            multiplying their bases.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('ZZ^x')
            sage: a = P(2^x)
            sage: b = P(3^x)
            sage: c = a._mul_(b); c
            6^x
            sage: c == a*b
            True
            sage: a*b*a  # indirect doctest
            12^x
        """
        return self.parent()(raw_element=self.base * other.base)


    def __invert__(self):
        r"""
        Return the multiplicative inverse of this exponential growth element.

        INPUT:

        Nothing.

        OUTPUT:

        The multiplicative inverse as a :class:`ExponentialGrowthElement`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: P = GrowthGroup('ZZ^x')
            sage: e1 = P(raw_element=2)
            sage: e2 = e1.__invert__(); e2
            (1/2)^x
            sage: e2 == ~e1
            True
            sage: e2.parent()
            Growth Group QQ^x

        ::

            sage: (~P(raw_element=1)).parent()
            Growth Group QQ^x
        """
        return self.parent()._create_element_via_parent_(1 / self.base)


    def __pow__(self, exponent):
        r"""
        Calculate the power of this growth element to the given ``exponent``.

        INPUT:

        - ``exponent`` -- a number. This can anything that is valid to be
          on the right hand side of ``*`` with an elements of the
          parent's base.

        OUTPUT:

        The result of this exponentiation a :class:`ExponentialGrowthElement`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('ZZ^x')
            sage: a = P(7^x); a
            7^x
            sage: b = a^(1/2); b
            sqrt(7)^x
            sage: b.parent()
            Growth Group SR^x
            sage: b^12
            117649^x
        """
        return self.parent()._create_element_via_parent_(self.base ** exponent)


    def _log_factor_(self, base=None):
        r"""
        Helper method for calculating the logarithm of the factorization
        of this element.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        OUTPUT:

        A tuple of pairs, where the first entry is either a growth
        element or something out of which we can construct a growth element
        and the second a multiplicative coefficient.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('QQ^x')
            sage: G('4^x').log_factor(base=2)  # indirect doctest
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot build log(4^x) since x is not in
            Growth Group QQ^x.
        """
        if self.is_one():
            return tuple()
        b = self.base
        if base is None and hasattr(b, 'is_monomial') and b.is_monomial() and \
                        b.variable_name() == 'e':
            coefficient = b.valuation()
        elif base is None and str(b) == 'e':
            coefficient = self.parent().base().one()
        else:
            from sage.functions.log import log
            coefficient = log(b, base=base)

        return ((str(self.parent()._var_), coefficient),)


    def _le_(self, other):
        r"""
        Return if this :class:`ExponentialGrowthElement` is at most
        (less than or equal to) ``other``.

        INPUT:

        - ``other`` -- a :class:`ExponentialGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`ExponentialGrowthElement`.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P_ZZ = agg.GrowthGroup('ZZ^x')
            sage: P_SR = agg.GrowthGroup('SR^x')
            sage: P_ZZ(2^x) <= P_SR(sqrt(3)^x)^2  # indirect doctest
            True
        """
        return bool(abs(self.base) <= abs(other.base))


class ExponentialGrowthGroup(GenericGrowthGroup):
    r"""
    A growth group dealing with expressions involving a fixed
    variable/symbol as the exponent.

    The elements :class:`ExponentialGrowthElement` of this group
    represent exponential functions with bases from a fixed base
    ring; the group law is the multiplication.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``).

      As exponential expressions are represented by this group,
      the elements in ``base`` are the bases of these exponentials.

    - ``var`` -- an object.

      The string representation of ``var`` acts as an exponent of the
      elements represented by this group.

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: P = agg.ExponentialGrowthGroup(QQ, 'x'); P
        Growth Group QQ^x

    .. SEEALSO::

        :class:`GenericGrowthGroup`
    """

    # enable the category framework for elements
    Element = ExponentialGrowthElement


    def _repr_short_(self):
        r"""
        A short representation string of this exponential growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.ExponentialGrowthGroup(QQ, 'a')  # indirect doctest
            Growth Group QQ^a


        TESTS::

            sage: agg.ExponentialGrowthGroup(QQ, 'a')._repr_short_()
            'QQ^a'
            sage: agg.ExponentialGrowthGroup(PolynomialRing(QQ, 'x'), 'a')._repr_short_()
            'QQ[x]^a'
        """
        from misc import parent_to_repr_short
        return '%s^%s' % (parent_to_repr_short(self.base()), self._var_)


    def _convert_(self, data):
        r"""
        Converts given ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.ExponentialGrowthGroup(ZZ, 'x')
            sage: P._convert_('icecream') is None
            True
            sage: P(1)  # indirect doctest
            1
            sage: P('2^x')  # indirect doctest
            2^x

        ::

            sage: P(2^x)  # indirect doctest
            2^x
            sage: P((-333)^x)  # indirect doctest
            (-333)^x
            sage: P(0)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: 0 is not in Growth Group ZZ^x.

        ::

            sage: P('7^x')
            7^x
            sage: P('(-2)^x')
            (-2)^x

        ::

            sage: P = agg.GrowthGroup('SR^x')
            sage: P(sqrt(3)^x)
            sqrt(3)^x
            sage: P((3^(1/3))^x)
            (3^(1/3))^x
            sage: P(e^x)
            e^x
            sage: P(exp(2*x))
            (e^2)^x
        """
        if data == 1 or data == '1':
            return self.base().one()
        var = repr(self._var_)
        try:
            P = data.parent()
        except AttributeError:
            if var not in str(data):
                return  # this has to end here

            elif str(data).endswith('^' + var):
                return self.base()(str(data).replace('^' + var, '')
                                   .replace('(', '').replace(')', ''))
            else:
                return  # end of parsing

        from sage.symbolic.ring import SR
        import operator
        from sage.symbolic.operators import mul_vararg
        if P is SR:
            op = data.operator()
            if op == operator.pow:
                base, exponent = data.operands()
                if str(exponent) == var:
                    return base
                elif exponent.operator() == mul_vararg:
                    return base ** (exponent / SR(var))
            elif isinstance(op, sage.functions.log.Function_exp):
                from sage.functions.log import exp
                base = exp(1)
                exponent = data.operands()[0]
                if str(exponent) == var:
                    return base
                elif exponent.operator() == mul_vararg:
                    return base ** (exponent / SR(var))


    def gens(self):
        r"""
        Return a tuple of all generators of this exponential growth
        group.

        INPUT:

        Nothing.

        OUTPUT:

        An empty tuple.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: E = GrowthGroup('ZZ^x')
            sage: E.gens()
            ()
        """
        return tuple()


    def construction(self):
        r"""
        Return the construction of this growth group.

        OUTPUT:

        A pair whose first entry is an
        :class:`exponential construction functor <ExponentialGrowthGroupFunctor>`
        and its second entry the base.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('QQ^x').construction()
            (ExponentialGrowthGroup[x], Rational Field)
        """
        return ExponentialGrowthGroupFunctor(self._var_), self.base()


class ExponentialGrowthGroupFunctor(AbstractGrowthGroupFunctor):
    r"""
    A :class:`construction functor <sage.categories.pushout.ConstructionFunctor>`
    for :class:`exponential growth groups <ExponentialGrowthGroup>`.

    INPUT:

    - ``var`` -- a string or list of strings (or anything else
      :class:`Variable` accepts).

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup, ExponentialGrowthGroupFunctor
        sage: GrowthGroup('QQ^z').construction()[0]
        ExponentialGrowthGroup[z]

    .. SEEALSO::

        :doc:`asymptotic_ring`,
        :class:`AbstractGrowthGroupFunctor`,
        :class:`MonomialGrowthGroupFunctor`,
        :class:`sage.rings.asymptotic.asymptotic_ring.AsymptoticRingFunctor`,
        :class:`sage.categories.pushout.ConstructionFunctor`.

    TESTS::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup, ExponentialGrowthGroupFunctor
        sage: cm = sage.structure.element.get_coercion_model()
        sage: A = GrowthGroup('QQ^x')
        sage: B = ExponentialGrowthGroupFunctor('x')(ZZ['t'])
        sage: cm.common_parent(A, B)
        Growth Group QQ[t]^x
    """

    _functor_name = 'ExponentialGrowthGroup'


    def __init__(self, var):
        r"""
        See :class:`ExponentialGrowthGroupFunctor` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import ExponentialGrowthGroupFunctor
            sage: ExponentialGrowthGroupFunctor('x')
            ExponentialGrowthGroup[x]
        """
        super(ExponentialGrowthGroupFunctor, self).__init__(var,
            sage.categories.monoids.Monoids())


    def _apply_functor(self, base):
        r"""
        Apply this functor to the given ``base``.

        INPUT:

        - ``base`` - anything :class:`ExponentialGrowthGroup` accepts.

        OUTPUT:

        An exponential growth group.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: F, R = GrowthGroup('QQ^z').construction()
            sage: F(R)  # indirect doctest
            Growth Group QQ^z
        """
        return ExponentialGrowthGroup(base, self.var)


class GrowthGroupFactory(sage.structure.factory.UniqueFactory):
    r"""
    A factory creating asymptotic growth groups.

    INPUT:

    - ``specification`` -- a string.

    OUTPUT:

    An asymptotic growth group.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: GrowthGroup('x^ZZ')
        Growth Group x^ZZ
        sage: GrowthGroup('log(x)^QQ')
        Growth Group log(x)^QQ

    This factory can also be used to construct Cartesian products
    of growth groups::

        sage: GrowthGroup('x^ZZ * y^ZZ')
        Growth Group x^ZZ * y^ZZ
        sage: GrowthGroup('x^ZZ * log(x)^ZZ')
        Growth Group x^ZZ * log(x)^ZZ
        sage: GrowthGroup('x^ZZ * log(x)^ZZ * y^QQ')
        Growth Group x^ZZ * log(x)^ZZ * y^QQ
        sage: GrowthGroup('QQ^x * x^ZZ * y^QQ * QQ^z')
        Growth Group QQ^x * x^ZZ * y^QQ * QQ^z

    TESTS::

        sage: TestSuite(GrowthGroup('x^ZZ')).run(verbose=True)  # long time
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass

    ::

        sage: TestSuite(GrowthGroup('QQ^y')).run(verbose=True)  # long time
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass

    ::

        sage: TestSuite(GrowthGroup('x^QQ * log(x)^ZZ')).run(verbose=True)  # long time
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
    """
    def create_key_and_extra_args(self, specification, **kwds):
        r"""
        Given the arguments and keyword, create a key that uniquely
        determines this object.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GrowthGroup.create_key_and_extra_args('asdf')
            Traceback (most recent call last):
            ...
            ValueError: 'asdf' is not a valid string describing a growth group.
        """
        from misc import split_str_by_mul
        factors = split_str_by_mul(specification)
        factors = tuple(f.replace('**', '^') for f in factors)

        for f in factors:
            if '^' not in f and '**' not in f:
                raise ValueError("'%s' is not a valid string describing "
                                 "a growth group." % (f,))

        return factors, kwds


    def create_object(self, version, factors, **kwds):
        r"""
        Create an object from the given arguments.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GrowthGroup('as^df')  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: 'as^df' is not a valid string describing a growth group.
            > *previous* ValueError: Cannot create a parent out of 'as'.
            >> *previous* SyntaxError: unexpected EOF while parsing (<string>, line 1)
            > *and* ValueError: Cannot create a parent out of 'df'.
            >> *previous* NameError: name 'df' is not defined
            sage: agg.GrowthGroup('x^y^z')
            Traceback (most recent call last):
            ...
            ValueError: 'x^y^z' is not a valid string describing a growth group.
            > *previous* ValueError: Cannot create a parent out of 'x'.
            >> *previous* NameError: name 'x' is not defined
            > *and* ValueError: Cannot create a parent out of 'y^z'.
            >> *previous* NameError: name 'y' is not defined
        """
        from misc import repr_short_to_parent
        groups = []
        for factor in factors:
            b, _, e = factor.partition('^')

            try:
                B = repr_short_to_parent(b)
            except ValueError as exc_b:
                B = None
            try:
                E = repr_short_to_parent(e)
            except ValueError as exc_e:
                E = None

            if B is None and E is None:
                from misc import combine_exceptions
                raise combine_exceptions(
                    ValueError("'%s' is not a valid string describing "
                               "a growth group." % (factor,)), exc_b, exc_e)
            elif B is None and E is not None:
                groups.append(MonomialGrowthGroup(E, b, **kwds))
            elif B is not None and E is None:
                groups.append(ExponentialGrowthGroup(B, e, **kwds))
            else:
                raise ValueError("'%s' is an ambigous string of a growth group "
                                 "description." % (factor,))

        if len(groups) == 1:
            return groups[0]

        from sage.categories.cartesian_product import cartesian_product
        return cartesian_product(groups)


GrowthGroup = GrowthGroupFactory("GrowthGroup")
