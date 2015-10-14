r"""
Asymptotic Expansions --- Miscellaneous

AUTHORS:

- Daniel Krenn (2015)

ACKNOWLEDGEMENT:

- Benjamin Hackl, Clemens Heuberger and Daniel Krenn are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

- Benjamin Hackl is supported by the Google Summer of Code 2015.


Functions, Classes and Methods
==============================
"""

#*****************************************************************************
# Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage


def repr_short_to_parent(s):
    r"""
    Helper method for the growth group factory, which converts a short
    representation string to a parent.

    INPUT:

    - ``s`` -- a string, short representation of a parent.

    OUTPUT:

    A parent.

    The possible short representations are shown in the examples below.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import repr_short_to_parent
        sage: repr_short_to_parent('ZZ')
        Integer Ring
        sage: repr_short_to_parent('QQ')
        Rational Field
        sage: repr_short_to_parent('SR')
        Symbolic Ring
        sage: repr_short_to_parent('NN')
        Non negative integer semiring

    TESTS::

        sage: repr_short_to_parent('abcdef')
        Traceback (most recent call last):
        ...
        ValueError: Cannot create a parent out of 'abcdef'.
        > *previous* NameError: name 'abcdef' is not defined
    """
    from sage.misc.sage_eval import sage_eval
    try:
        P = sage_eval(s)
    except Exception as e:
        raise combine_exceptions(
            ValueError("Cannot create a parent out of '%s'." % (s,)), e)

    from sage.misc.lazy_import import LazyImport
    if type(P) is LazyImport:
        P = P._get_object()

    from sage.structure.parent import is_Parent
    if not is_Parent(P):
        raise ValueError("'%s' does not describe a parent." % (s,))
    return P


def parent_to_repr_short(P):
    r"""
    Helper method which generates a short(er) representation string
    out of a parent.

    INPUT:

    - ``P`` -- a parent.

    OUTPUT:

    A string.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import parent_to_repr_short
        sage: parent_to_repr_short(ZZ)
        'ZZ'
        sage: parent_to_repr_short(QQ)
        'QQ'
        sage: parent_to_repr_short(SR)
        'SR'
        sage: parent_to_repr_short(ZZ['x'])
        'ZZ[x]'
        sage: parent_to_repr_short(QQ['d, k'])
        '(QQ[d, k])'
        sage: parent_to_repr_short(QQ['e'])
        'QQ[e]'
        sage: parent_to_repr_short(SR[['a, r']])
        '(SR[[a, r]])'
        sage: parent_to_repr_short(Zmod(3))
        '(Ring of integers modulo 3)'
        sage: parent_to_repr_short(Zmod(3)['g'])
        '(Univariate Polynomial Ring in g over Ring of integers modulo 3)'
    """
    def abbreviate(P):
        if P is sage.rings.integer_ring.ZZ:
            return 'ZZ'
        elif P is sage.rings.rational_field.QQ:
            return 'QQ'
        elif P is sage.symbolic.ring.SR:
            return 'SR'
        raise ValueError('Cannot abbreviate %s.' % (P,))

    poly = sage.rings.polynomial.polynomial_ring.is_PolynomialRing(P) or \
           sage.rings.polynomial.multi_polynomial_ring_generic.is_MPolynomialRing(P)
    from sage.rings import multi_power_series_ring
    power = sage.rings.power_series_ring.is_PowerSeriesRing(P) or \
            multi_power_series_ring.is_MPowerSeriesRing(P)

    if poly or power:
        if poly:
            op, cl = ('[', ']')
        else:
            op, cl = ('[[', ']]')
        try:
            s = abbreviate(P.base_ring()) + op + ', '.join(P.variable_names()) + cl
        except ValueError:
            s = str(P)
    else:
        try:
            s = abbreviate(P)
        except ValueError:
            s = str(P)

    if ' ' in s:
        s = '(' + s + ')'
    return s


def split_str_by_op(string, op, strip_parentheses=True):
    r"""
    Split the given string into a tuple of substrings arising by
    splitting by '*' and taking care of parentheses.

    INPUT:

    - ``string`` -- a string.

    - ``op`` -- a string.

    - ``strip_parentheses`` -- (default: ``True``) a boolean.

    OUTPUT:

    A tuple of strings.

    TESTS::

        sage: from sage.rings.asymptotic.misc import split_str_by_op
        sage: split_str_by_op('x^ZZ', '*')
        ('x^ZZ',)
        sage: split_str_by_op('log(x)^ZZ * y^QQ', '*')
        ('log(x)^ZZ', 'y^QQ')
        sage: split_str_by_op('log(x)**ZZ * y**QQ', '*')
        ('log(x)**ZZ', 'y**QQ')
        sage: split_str_by_op('a^b * * c^d', '*')
        Traceback (most recent call last):
        ...
        ValueError: 'a^b * * c^d' is invalid since a '*' follows a '*'.
        sage: split_str_by_op('a^b * (c*d^e)', '*')
        ('a^b', 'c*d^e')

    ::

        sage: split_str_by_op('(a^b)^c', '^')
        ('a^b', 'c')
        sage: split_str_by_op('a^(b^c)', '^')
        ('a', 'b^c')
    """
    factors = list()
    balanced = True
    if string and op is not None and string.startswith(op):
        raise ValueError("'%s' is invalid since it starts with a '%s'." %
                         (string, op))
    for s in string.split(op):
        if not s:
            factors[-1] += op
            balanced = False
            continue
        if not s.strip():
            raise ValueError("'%s' is invalid since a '%s' follows a '%s'." %
                             (string, op, op))
        if not balanced:
            s = factors.pop() + op + s
        balanced = s.count('(') == s.count(')')
        factors.append(s)

    if not balanced:
        raise ValueError("Parentheses in '%s' are not balanced." % (string,))

    def strip(s):
        s = s.strip()
        if not s:
            return s
        if strip_parentheses and s[0] == '(' and s[-1] == ')':
            s = s[1:-1]
        return s.strip()

    return tuple(strip(f) for f in factors)


def repr_op(left, op, right=None):
    r"""
    Create a string ``left op right`` with
    taking care of parentheses in its operands.

    INPUT:

    - ``left`` -- an element.

    - ``op`` -- a string.

    - ``right`` -- an alement.

    OUTPUT:

    A string.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import repr_op
        sage: repr_op('a^b', '^', 'c')
        '(a^b)^c'
    """
    left = str(left)
    right = str(right) if right is not None else ''

    def add_parentheses(s, op):
        if op == '^':
            signals = ('^', '*', '+', ' ')
        else:
            return s
        if any(sig in s for sig in signals):
            return '(%s)' % (s,)
        else:
            return s

    return add_parentheses(left, op) + op +\
        add_parentheses(right, op)


def combine_exceptions(e, *f):
    r"""
    Helper function which combines the messages of the given exceptions.

    INPUT:

    - ``e`` -- an exception.

    - ``*f`` -- exceptions.

    OUTPUT:

    An exception.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import combine_exceptions
        sage: raise combine_exceptions(ValueError('Outer.'), TypeError('Inner.'))
        Traceback (most recent call last):
        ...
        ValueError: Outer.
        > *previous* TypeError: Inner.
        sage: raise combine_exceptions(ValueError('Outer.'),
        ....:                          TypeError('Inner1.'), TypeError('Inner2.'))
        Traceback (most recent call last):
        ...
        ValueError: Outer.
        > *previous* TypeError: Inner1.
        > *and* TypeError: Inner2.
        sage: raise combine_exceptions(ValueError('Outer.'),
        ....:                          combine_exceptions(TypeError('Middle.'),
        ....:                                             TypeError('Inner.')))
        Traceback (most recent call last):
        ...
        ValueError: Outer.
        > *previous* TypeError: Middle.
        >> *previous* TypeError: Inner.
    """
    import re
    msg = ('\n *previous* ' +
           '\n *and* '.join("%s: %s" % (ff.__class__.__name__, str(ff)) for ff in f))
    msg = re.sub(r'^([>]* \*previous\*)', r'>\1', msg, flags=re.MULTILINE)
    msg = re.sub(r'^([>]* \*and\*)', r'>\1', msg, flags=re.MULTILINE)
    msg = str(e.args if len(e.args) > 1 else e.args[0]) + msg
    e.args = (msg,)
    return e


def underlying_class(P):
    r"""
    Return the underlying class (class without the attached
    categories) of the given instance.

    OUTPUT:

    A class.

    EXAMPLES::

        sage: from sage.rings.asymptotic.misc import underlying_class
        sage: type(QQ)
        <class 'sage.rings.rational_field.RationalField_with_category'>
        sage: underlying_class(QQ)
        <class 'sage.rings.rational_field.RationalField'>
    """
    cls = type(P)
    if not P._is_category_initialized():
        return cls
    from sage.structure.misc import is_extension_type
    if is_extension_type(cls):
        return cls

    from sage.categories.sets_cat import Sets
    Sets_parent_class = Sets().parent_class
    while issubclass(cls, Sets_parent_class):
        cls = cls.__base__
    return cls


def merge_overlapping(A, B, key=None):
    r"""
    Merge the two overlapping tuples/lists.

    INPUT:

    - ``A`` -- a list or tuple (type has to coincide with type of ``B``).

    - ``B`` -- a list or tuple (type has to coincide with type of ``A``).

    - ``key`` -- (default: ``None``) a function. If ``None``, then the
      identity is used.  This ``key``-function applied on an element
      of the list/tuple is used for comparison. Thus elements with the
      same key are considered as equal.

    OUTPUT:

    A pair of lists or tuples (depending on the type of ``A`` and ``B``).

    .. NOTE::

        Suppose we can decompose the list `A=ac` and `B=cb` with
        lists `a`, `b`, `c`, where `c` is nonempty. Then
        :func:`merge_overlapping` returns the pair `(acb, acb)`.

        Suppose a ``key``-function is specified and `A=ac_A` and
        `B=c_Bb`, where the list of keys of the elements of `c_A`
        equals the list of keys of the elements of `c_B`. Then
        :func:`merge_overlapping` returns the pair `(ac_Ab, ac_Bb)`.

        After unsuccessfully merging `A=ac` and `B=cb`,
        a merge of `A=ca` and `B=bc` is tried.

    TESTS::

        sage: from sage.rings.asymptotic.misc import merge_overlapping
        sage: def f(L, s):
        ....:     return list((ell, s) for ell in L)
        sage: key = lambda k: k[0]
        sage: merge_overlapping(f([0..3], 'a'), f([5..7], 'b'), key)
        Traceback (most recent call last):
        ...
        ValueError: Input does not have an overlap.
        sage: merge_overlapping(f([0..2], 'a'), f([4..7], 'b'), key)
        Traceback (most recent call last):
        ...
        ValueError: Input does not have an overlap.
        sage: merge_overlapping(f([4..7], 'a'), f([0..2], 'b'), key)
        Traceback (most recent call last):
        ...
        ValueError: Input does not have an overlap.
        sage: merge_overlapping(f([0..3], 'a'), f([3..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a'), (4, 'b')],
         [(0, 'a'), (1, 'a'), (2, 'a'), (3, 'b'), (4, 'b')])
        sage: merge_overlapping(f([3..4], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'b'), (2, 'b'), (3, 'a'), (4, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b'), (4, 'a')])
        sage: merge_overlapping(f([0..1], 'a'), f([0..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'b'), (3, 'b'), (4, 'b')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b'), (4, 'b')])
        sage: merge_overlapping(f([0..3], 'a'), f([0..1], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'a'), (3, 'a')])
        sage: merge_overlapping(f([0..3], 'a'), f([1..3], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'a'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: merge_overlapping(f([1..3], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: merge_overlapping(f([0..6], 'a'), f([3..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a'), (4, 'a'), (5, 'a'), (6, 'a')],
         [(0, 'a'), (1, 'a'), (2, 'a'), (3, 'b'), (4, 'b'), (5, 'a'), (6, 'a')])
        sage: merge_overlapping(f([0..3], 'a'), f([1..2], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'a'), (1, 'b'), (2, 'b'), (3, 'a')])
        sage: merge_overlapping(f([1..2], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'a'), (2, 'a'), (3, 'b')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: merge_overlapping(f([1..3], 'a'), f([1..3], 'b'), key)
        ([(1, 'a'), (2, 'a'), (3, 'a')],
         [(1, 'b'), (2, 'b'), (3, 'b')])
    """
    if key is None:
        Akeys = A
        Bkeys = B
    else:
        Akeys = tuple(key(a) for a in A)
        Bkeys = tuple(key(b) for b in B)

    def find_overlapping_index(A, B):
        if len(B) > len(A) - 2:
            raise StopIteration
        matches = iter(i for i in xrange(1, len(A) - len(B))
                       if A[i:i+len(B)] == B)
        return next(matches)

    def find_mergedoverlapping_index(A, B):
        """
        Return in index i where to merge two overlapping tuples/lists ``A`` and ``B``.

        Then ``A + B[i:]`` or ``A[:-i] + B`` are the merged tuples/lists.

        Adapted from http://stackoverflow.com/a/30056066/1052778.
        """
        matches = iter(i for i in xrange(min(len(A), len(B)), 0, -1)
                       if A[-i:] == B[:i])
        return next(matches, 0)

    i = find_mergedoverlapping_index(Akeys, Bkeys)
    if i > 0:
        return A + B[i:], A[:-i] + B

    i = find_mergedoverlapping_index(Bkeys, Akeys)
    if i > 0:
        return B[:-i] + A, B + A[i:]

    try:
        i = find_overlapping_index(Akeys, Bkeys)
    except StopIteration:
        pass
    else:
        return A, A[:i] + B + A[i+len(B):]

    try:
        i = find_overlapping_index(Bkeys, Akeys)
    except StopIteration:
        pass
    else:
        return B[:i] + A + B[i+len(A):], B

    raise ValueError('Input does not have an overlap.')
