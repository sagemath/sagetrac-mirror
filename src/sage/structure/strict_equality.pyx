r"""
Strict equality comparison switch

Some elements implement equality in a way that is against python's
recommendations. With a stricter version of equality, which can be enabled
temporarily, the equality operator respects python's recommendations and in
particular such objects can be used as keys in dictionaries.

AUTHORS:

- Julian Rüth (2016-03-22): initial version

Some elements (such as `p`-adics) implement equality in a way that is against
python's recommendations::

    sage: R = Zp(3)
    sage: a = R(0,1); a
    sage: b = R(3,2); b
    sage: c = R(6,2); c
    sage: a == b
    True
    sage: a == c
    True
    sage: b == c
    False

This makes it impossible to use such elements consistently as keys in dictionaries. Therefore, hashing is disabled for them::

    sage: D = {}
    sage: a in D

With a stricter version of equality, the equality operator respects python's recommendations::

    sage: with strict_equality(True):
    ....:     a == b || a == c || b == c
    False

In this way these elements can be used in dictionaries for example for the purpose of caching results::

    sage: with strict_equality(True):
    ....:     a in D
    False

"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from libcpp.stack cimport stack

def strict_equality(use_strict_equality = None):
    r"""
    Enable or disable strict equality checks temporarily (or query the
    currently active strictness.)

    INPUT:

    - ``use_strict_equality`` -- ``True``, ``False`` or ``None`` (default)

    .. WARNING::

        Be careful only to set ``strict_equality`` for the lines of code where it
        is necessary. And make sure to explicitly unset it if you query it (for
        example in a ``_cmp`` or ``__hash__`` method if you call subroutines that
        expect it not to be enabled.) Not doing so might produce caching errors
        which are hard to track down.

    EXAMPLES:

    Without a parameter the function reports whether strict equality checking
    is currently enabled::

        sage: strict_equality()
        False

    With the parameter ``True``, the returned object can be used in a ``with``
    statement and will temporarily enable strict equality cheking. After the
    ``with`` block, the previous state is restored::

        sage: with strict_equality(True):
        ....:     strict_equality()
        True
        sage: strict_equality()
        False

    Equally, ``False`` temporarily disables strict checking::

        sage: with strict_equality(True):
        ....:     with strict_equality(False):
        ....:         strict_equality()
        False
        sage: strict_equality()
        False

    """
    if use_strict_equality is None:
        return get_strict_equality()
    else:
        return StrictEqualityContextManager(use_strict_equality)

cdef stack[bint] strict_equality_history
strict_equality_history.push(False)

cdef void enter_strict_equality(bint use_strict_equality):
    r"""
    Set the current strictness for equality comparisons to
    ``use_strict_equality``.

    EXAMPLES::

        sage: with strict_equality(True): # indirect doctest
        ....:     strict_equality()
        True
    """
    strict_equality_history.push(use_strict_equality)

cdef void leave_strict_equality():
    r"""
    Restore previous strictness of equality comparison.

    EXAMPLES::

        sage: with strict_equality(True): pass # indirect doctest
        sage: strict_equality()
        False
    """
    strict_equality_history.pop()

cdef bint get_strict_equality():
    r"""
    Return whether strict equality checking is currently enabled, i.e., the top
    of ``strict_equality_history``.

    EXAMPLES::

        sage: strict_equality() # indirect doctest
        False
    """
    return strict_equality_history.top()

class StrictEqualityContextManager(object):
    r"""
    A ContextManager which pushes to `strict_equality_history` on enter and
    pops from it on exit.

    EXAMPLES::

        sage: type(strict_equality(True))
        sage: strict_equality()
        False
    """
    def __init__(self, use_strict_equality):
        self.use_strict_equality = use_strict_equality

    def __enter__(self):
        r"""
        TESTS::

        sage: with strict_equality(True):
        ....:     strict_equality()
        True
        """
        enter_strict_equality(self.use_strict_equality)

    def __exit__(self, typ, value, traceback):
        r"""
        TESTS::

        sage: with strict_equality(True):
        ....:     pass
        ....: strict_equality()
        False
        """
        leave_strict_equality()
