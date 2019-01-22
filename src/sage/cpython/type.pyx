"""
Type internals
"""

from cpython.object cimport Py_TYPE, PyTypeObject, Py_TPFLAGS_HEAPTYPE
from cpython.ref cimport PyObject, Py_INCREF, Py_XINCREF

cdef PyTypeObject* PyInstance_Type = NULL
try:
    from types import InstanceType
    PyInstance_Type = <PyTypeObject*>InstanceType
except ImportError:
    pass

cdef extern from *:
    """
    #define has_del(obj) ((obj)->tp_del != NULL)
    """
    void PyObject_GC_UnTrack(op)
    PyObject* PyDict_GetItem(PyObject*, key)
    bint has_del(PyTypeObject*)


cpdef bint can_assign_class(obj):
    """
    Can we assign ``obj.__class__``?

    Note that Python 3.5 has experimented with allowing assigning
    ``__class__`` in more cases but this was mostly reverted. In this
    function, we apply the Python 2.7 semantics.

    EXAMPLES::

        sage: class X: pass
        sage: from sage.cpython.type import can_assign_class
        sage: can_assign_class(X())
        True
        sage: class Y(int): pass
        sage: from sage.cpython.type import can_assign_class
        sage: can_assign_class(Y())
        True
        sage: can_assign_class(1)
        False
    """
    cdef PyTypeObject* tp = Py_TYPE(obj)
    if tp is PyInstance_Type:
        return True
    return (tp.tp_flags & Py_TPFLAGS_HEAPTYPE) != 0


cdef void extension_class_dealloc(PyObject* self):
    """
    Deallocator for classes created with :func:`extension_class`.
    """
    origcls = PyDict_GetItem(self.ob_type.tp_dict, '__heap_type')
    if origcls is NULL:
        # This exception is not actually raised, just written out
        raise RuntimeError("No __heap_type attribute found")

    # Change type of self
    Py_XINCREF(origcls)
    self.ob_type = <PyTypeObject*>origcls

    # Delegate to the deallocator of origcls
    (<PyTypeObject*>origcls).tp_dealloc(self)


def extension_class(cls):
    """
    Return a copy of class ``cls`` as extension class.

    If ``cls`` is already an extension class, simply return ``cls``.

    The main use case is to allow Cython optimizations for ``cpdef``
    methods which only work for extension classes.

    .. NOTE::

        Since extension classes do not support ``__del__``, a
        ``NotImplementedError`` is raised when calling
        ``extension_class`` on a class with ``__del__``.

    .. WARNING::

        Extension classes are never deallocated, so this can lead to
        memory leaks if used badly.

    INPUT:

    - ``cls`` -- any class (not a Python 2 old-style class)

    EXAMPLES::

        sage: from sage.cpython.type import extension_class
        sage: class Test(object):
        ....:     def __repr__(self):
        ....:         return "OK!"
        sage: Test()
        OK!
        sage: T = extension_class(Test)
        sage: T
        <... 'Test'>
        sage: T()
        OK!

    We can see the difference between ``Test`` and ``T`` using
    :func:`can_assign_class`::

        sage: from sage.cpython.type import can_assign_class
        sage: can_assign_class(Test())
        True
        sage: can_assign_class(T())
        False

    The metaclass is correctly preserved::

        sage: from sage.misc.fast_methods import Singleton
        sage: S = extension_class(Singleton)
        sage: type(S)
        <... 'sage.misc.classcall_metaclass.ClasscallMetaclass'>

    When used on an extension class, we just get that class back::

        sage: extension_class(list) is list
        True

    Classes with ``__del__`` cannot be supported::

        sage: @extension_class
        ....: class X(object):
        ....:     def __del__(self):
        ....:         pass
        Traceback (most recent call last):
        ...
        NotImplementedError: extension_class() does not support classes with __del__

    TESTS:

    We test various special methods::

        sage: @extension_class
        ....: class Special(object):
        ....:     def __repr__(self):
        ....:         return "<instance of Special>"
        ....:     def __add__(self, other):
        ....:         print("__add__")
        ....:     def __iter__(self):
        ....:         yield 42
        sage: x = Special()
        sage: x
        <instance of Special>
        sage: x + x
        __add__
        sage: list(x)
        [42]

    Check that attribute assignment and deallocation work properly::

        sage: class DeleteMe:
        ....:     def __del__(self):
        ....:         print("deleted!")
        sage: @extension_class
        ....: class Object(object):
        ....:     pass
        sage: x = Object()
        sage: x.attr = DeleteMe()
        sage: del x
        deleted!

    The same example with ``__slots__``::

        sage: @extension_class
        ....: class Object(object):
        ....:     __slots__ = ('attr',)
        sage: x = Object()
        sage: x.attr = DeleteMe()
        sage: del x
        deleted!

    These classes are not tracked by the garbage collector::

        sage: from gc import is_tracked
        sage: is_tracked(Test)
        True
        sage: is_tracked(T)
        False
    """
    if not isinstance(cls, type):
        raise TypeError(f"{cls!r} is not a type")

    if (<PyTypeObject*>cls).tp_flags & Py_TPFLAGS_HEAPTYPE == 0:
        # Already an extension class
        return cls

    # First make a copy of "cls" as ordinary Python class.
    # In the new class, we store the original class as "__heap_type"
    # attribute
    D = dict(cls.__dict__)
    D['__heap_type'] = cls
    # If __slots__ is used, remove the member descriptors for the
    # slots. Copying these doesn't work, they need to be generated
    # from scratch.
    try:
        slots = D['__slots__']
    except KeyError:
        pass
    else:
        for slot in slots:
            del D[slot]
    meta = type(cls)

    # Use tp_new instead of __new__ to bypass metaclass check
    args = (cls.__name__, cls.__bases__, D)
    newcls = (<PyTypeObject*>meta).tp_new(meta, args, <dict>NULL)

    if has_del(<PyTypeObject*>newcls):
        raise NotImplementedError("extension_class() does not support classes with __del__")

    # Extension classes are not refcounted and do not participate in
    # garbage collection. So we remove it from GC and artificially
    # increase the refcount to ensure that our new class will never
    # be deallocated.
    PyObject_GC_UnTrack(newcls)
    Py_INCREF(newcls)

    # There is one extra complication: the default deallocator for
    # instances only works correctly for actual heap types and there is
    # no way around this. So we have a custom deallocator which
    # dynamically changes the type from "newcls" to "cls" and then calls
    # the original deallocator. This probably only works correctly if
    # __del__ is not set, so we checked that above.
    (<PyTypeObject*>newcls).tp_dealloc = extension_class_dealloc

    # Change it to an extension class
    (<PyTypeObject*>newcls).tp_flags &= ~Py_TPFLAGS_HEAPTYPE

    return newcls
