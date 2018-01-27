"""
The Unknown truth value

The ``Unknown`` object is used in Sage in several places as return value
in addition to ``True`` and ``False``, in order to signal uncertainty
about or inability to compute the result. ``Unknown`` can be identified
using ``is``, or by catching :class:`UnknownError` from a boolean operation.

.. WARNING::

    Calling ``bool()`` with ``Unknown`` as argument will throw an
    ``UnknownError``. This also means that, unless PEP 335 is accepted,
    in the following cases, ``and``, ``not``, and ``or`` fail or return
    a somewhat wrong value::

        sage: not Unknown         # should return Unknown
        Traceback (most recent call last):
        ...
        UnknownError: Unknown does not evaluate in boolean context
        sage: Unknown and False   # should return False
        Traceback (most recent call last):
        ...
        UnknownError: Unknown does not evaluate in boolean context
        sage: Unknown or False    # should return Unknown
        Traceback (most recent call last):
        ...
        UnknownError: Unknown does not evaluate in boolean context

EXAMPLES::

    sage: [Unknown is val for val in [True, False, Unknown]]
    [False, False, True]

    sage: def func(n): return [True, False, Unknown][n]
    sage: def verbose(n):
    ....:    try:
    ....:        if func(n):
    ....:            return 'yes'
    ....:        else:
    ....:            return 'no'
    ....:    except UnknownError:
    ....:        return 'maybe'
    sage: [verbose(n) for n in range(3)]
    ['yes', 'no', 'maybe']

AUTHORS:

- Florent Hivert (2010): initial version.
"""
from __future__ import print_function

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.richcmp import richcmp_method, rich_to_bool

class UnknownError(TypeError):
    """
    Raised whenever :class:`Unknown` is used in a boolean operation.

    EXAMPLES::

        sage: not Unknown
        Traceback (most recent call last):
        ...
        UnknownError: Unknown does not evaluate in boolean context
    """
    pass

@richcmp_method
class UnknownClass(UniqueRepresentation, SageObject):
    """
    The Unknown truth value

    The ``Unknown`` object is used in Sage in several places as return value
    in addition to ``True`` and ``False``, in order to signal uncertainty
    about or inability to compute the result. ``Unknown`` can be identified
    using ``is``, or by catching :class:`UnknownError` from a boolean
    operation.

    .. WARNING::

        Calling ``bool()`` with ``Unknown`` as argument will throw an
        ``UnknownError``. This also means that applying ``and``, ``not``,
        and ``or`` to ``Unknown`` might fail.

    TESTS::

        sage: TestSuite(Unknown).run()
    """
    def _repr_(self):
        """
        TESTS::

            sage: Unknown
            Unknown
        """
        return "Unknown"

    def __bool__(self):
        """
        When evaluated in a boolean context ``Unknown()`` raises an error.

        EXAMPLES::

            sage: bool(Unknown)
            Traceback (most recent call last):
            ...
            UnknownError: Unknown does not evaluate in boolean context
            sage: not Unknown
            Traceback (most recent call last):
            ...
            UnknownError: Unknown does not evaluate in boolean context
        """
        raise UnknownError('Unknown does not evaluate in boolean context')

    __nonzero__ = __bool__

    def __and__(self, other):
        """
        The ``and`` logical connector.

        .. WARNING::

            This is not used by ``and`` unless PEP 335 is accepted.

        EXAMPLES::

            sage: Unknown & False
            False
            sage: Unknown & Unknown
            Unknown
            sage: Unknown & True
            Unknown

        Compare with::

            sage: Unknown and False    # should return False
            Traceback (most recent call last):
            ...
            UnknownError: Unknown does not evaluate in boolean context
        """
        if other is False:
            return False
        else:
            return self

    def __or__(self, other):
        """
        The ``or`` logical connector.

        .. WARNING::

            This is not used by ``or`` unless PEP 335 is accepted.

        EXAMPLES::

            sage: Unknown | False
            Unknown
            sage: Unknown | Unknown
            Unknown
            sage: Unknown | True
            True

        Compare with::

            sage: Unknown or False    # should return Unknown
            Traceback (most recent call last):
            ...
            UnknownError: Unknown does not evaluate in boolean context
        """
        if other is True:
            return True
        else:
            return self

    def __not__(self):
        """
        The ``not`` logical connector.

        .. WARNING::

            This is not used by ``not`` unless PEP 335 is accepted.

        EXAMPLES::

            sage: Unknown.__not__()
            Unknown

        Compare with::

            sage: not Unknown  # should return Unknown
            Traceback (most recent call last):
            ...
            UnknownError: Unknown does not evaluate in boolean context
        """
        return self

    def __richcmp__(self, other, op):
        """
        Comparison of truth value.

        EXAMPLES::

            sage: l = [False, Unknown, True]
            sage: for a in l: print([a < b for b in l])
            [False, True, True]
            [False, False, True]
            [False, False, False]

            sage: for a in l: print([a <= b for b in l])
            [True, True, True]
            [False, True, True]
            [False, False, True]
        """
        if other is self:
            return rich_to_bool(op, 0)
        if not isinstance(other, bool):
            return NotImplemented
        if other:
            return rich_to_bool(op, -1)
        else:
            return rich_to_bool(op, +1)


Unknown = UnknownClass()
