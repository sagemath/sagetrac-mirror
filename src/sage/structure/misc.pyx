"""
Miscellaneous utilities
"""


cdef class SafeSortable:
    """
    Wrapper for arbitrary objects allowing them to always be ordered in the
    manner of Python 2.x objects.

    When comparing a `SafeSortable` to another `SafeSortable`, it first
    attempts comparing their wrapped objects using the standard ``__lt__``
    comparison.  If that fails, it falls back to comparing them by their type
    name, then by their id, which is consistent with what Python 2.x does with
    a pair of otherwise unorderable types.

    `SafeSortable`s may only be compared with other `SafeSortable`s.

    EXAMPLES::

        sage: int(1) < 'a'  # py3
        Traceback (most recent call last):
        ...
        TypeError: '<' not supported between instances of 'int' and 'str'
        sage: from sage.structure.misc import SafeSortable
        sage: SafeSortable(int(1)) < SafeSortable('a')
        True
        sage: class ExampleClass(object): pass
        sage: e1 = ExampleClass()
        sage: e2 = ExampleClass()
        sage: e1 < e2  # py3
        Traceback (most recent call last):
        ...
        TypeError: '<' not supported between instances of 'ExampleClass' and
        'ExampleClass'
        sage: (SafeSortable(e1) < SafeSortable(e2)) == (id(e1) < id(e2))
        True
        sage: int(1) < SafeSortable('a')  # py3
        Traceback (most recent call last):
        ...
        TypeError: '<' not supported between instances of 'int' and
        'sage.structure.misc.SafeSortable'
    """

    def __cinit__(self, obj):
        self.obj = obj

    def __lt__(left, right):
        cdef object l, r

        if not (type(left) is SafeSortable and type(right) is SafeSortable):
            return NotImplemented

        l = (<SafeSortable>left).obj
        r = (<SafeSortable>right).obj

        try:
            return l < r
        except TypeError:
            return (str(type(l)), id(l)) < (str(type(r)), id(r))

    @classmethod
    def lt(cls, left, right):
        """
        Convenience method for evaluating ``left < right`` with each operand
        wrapped with ``SafeSortable``.

        Just a shortcut for ``SafeSortable(left) < SafeSortable(right)``.

        EXAMPLES::

            sage: from sage.structure.misc import SafeSortable
            sage: SafeSortable.lt(1, 'a')
            True
        """

        return cls(left) < cls(right)

    @classmethod
    def sort(cls, iterable, key=None, reverse=False):
        """
        Convenience method for safe-sorting an iterable by wrapping each
        element in ``SafeSortable``.

        If a ``key`` argument is given, the key is wrapped in ``SafeSortable``,
        so that each object returned by the key function an be compared.

        EXAMPLES::

            sage: from sage.structure.misc import SafeSortable
            sage: L = [int(1), 'a', complex(1+2j), 'b']
            sage: SafeSortable.sort(L)
            [(1+2j), 1, 'a', 'b']
            sage: SafeSortable.sort(L, key=lambda x: x * int(2))
            [(1+2j), 1, 'a', 'b']
        """

        if key is None:
            sort_key = cls
        else:
            sort_key = lambda x: cls(key(x))

        return sorted(iterable, key=sort_key, reverse=reverse)


def is_extension_type(cls):
    """
    INPUT:

    - cls: a class

    Tests whether cls is an extension type (int, list, cython compiled classes, ...)

    EXAMPLES::

        sage: from sage.structure.parent import is_extension_type
        sage: is_extension_type(int)
        True
        sage: is_extension_type(list)
        True
        sage: is_extension_type(ZZ.__class__)
        True
        sage: is_extension_type(QQ.__class__)
        False
    """
    # Robert B claims that this should be robust
    try:
        return cls.__dictoffset__ == 0
    except AttributeError:
        pass
    return False
