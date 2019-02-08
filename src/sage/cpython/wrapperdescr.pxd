from cpython.descr cimport wrapper_descriptor, wrapperbase


cdef wrapperdescr_fastcall(wrapper_descriptor slotwrapper, self, args, kwds)


cdef inline wrapperbase* get_slotdef(wrapper_descriptor slotwrapper) except NULL:
    """
    Given a slot wrapper, return the corresponding ``slotdef``.

    A ``slotdef`` is associated to a specific slot like ``__eq__``
    and does not depend at all on the type. In other words, calling
    ``get_slotdef(t.__eq__)`` will return the same ``slotdef``
    independent of the type ``t`` (provided that the type implements
    rich comparison in C).

    TESTS::

        sage: cython('''
        ....: from sage.cpython.wrapperdescr cimport get_slotdef
        ....: from cpython.long cimport PyLong_FromVoidPtr
        ....: def py_get_slotdef(slotwrapper):
        ....:     return PyLong_FromVoidPtr(get_slotdef(slotwrapper))
        ....: ''')
        sage: py_get_slotdef(object.__init__)  # random
        140016903442416
        sage: py_get_slotdef(bytes.__lt__)  # random
        140016903441800
        sage: py_get_slotdef(bytes.__lt__) == py_get_slotdef(Integer.__lt__)
        True
        sage: py_get_slotdef(bytes.__lt__) == py_get_slotdef(bytes.__gt__)
        False
        sage: class X(object):
        ....:     def __eq__(self, other):
        ....:         return False
        sage: py_get_slotdef(X.__eq__)
        Traceback (most recent call last):
        ...
        TypeError: Cannot convert ... to wrapper_descriptor
    """
    return slotwrapper.d_base
