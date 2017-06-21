r"""
Cython-like rich comparisons in Python

With "rich comparisons", we mean the Python 3 comparisons which are
usually implemented in Python using methods like ``__eq__`` and
``__lt__``. Internally in Python, there is only one rich comparison
slot ``tp_richcompare``. The actual operator is passed as an integer
constant (defined in this module as
``op_LT``, ``op_LE``, ``op_EQ``, ``op_NE``, ``op_GT``, ``op_GE``).

Cython exposes rich comparisons in ``cdef`` classes as the
``__richcmp__`` special method. The Sage coercion model also supports
rich comparisons this way: for two instances ``x`` and ``y``
of :class:`~sage.structure.element.Element`, ``x._richcmp_(y, op)``
is called when the user does something like ``x <= y``
(possibly after coercion if ``x`` and ``y`` have different parents).

Various helper functions exist to make it easier to implement rich
comparison: the most important one is the :func:`richcmp` function.
This is analogous to the Python 2 function ``cmp()`` but implements
rich comparison, with the comparison operator (e.g. ``op_GE``) as
third argument. There is also :func:`richcmp_not_equal` which is like
:func:`richcmp` but it is optimized assuming that the compared objects
are not equal.

The functions :func:`rich_to_bool` and :func:`rich_to_bool_sgn` can be
used to convert results of ``cmp()`` (i.e. -1, 0 or 1) to a boolean
``True``/``False`` for rich comparisons.

AUTHORS:

- Jeroen Demeyer
"""

from cpython.object cimport PyObject_RichCompare, Py_TYPE, PyTypeObject

op_LT = Py_LT   # operator <
op_LE = Py_LE   # operator <=
op_EQ = Py_EQ   # operator ==
op_NE = Py_NE   # operator !=
op_GT = Py_GT   # operator >
op_GE = Py_GE   # operator >=


# Slotdefs for richcmp methods (the "bytes" below is arbitrary,
# any type which implements these methods would work)
cdef wrapperbase* richcmp_slotdef[6]
richcmp_slotdef[Py_EQ] = get_slotdef(bytes.__eq__)
richcmp_slotdef[Py_NE] = get_slotdef(bytes.__ne__)
richcmp_slotdef[Py_LT] = get_slotdef(bytes.__lt__)
richcmp_slotdef[Py_GT] = get_slotdef(bytes.__gt__)
richcmp_slotdef[Py_LE] = get_slotdef(bytes.__le__)
richcmp_slotdef[Py_GE] = get_slotdef(bytes.__ge__)


def richcmp(x, y, int op):
    """
    Return the result of the rich comparison of ``x`` and ``y`` with
    operator ``op``.

    INPUT:

    - ``x``, ``y`` -- arbitrary Python objects

    - ``op`` -- comparison operator (one of ``op_LT`, ``op_LE``,
      ``op_EQ``, ``op_NE``, ``op_GT``, ``op_GE``).

    EXAMPLES::

        sage: from sage.structure.richcmp import *
        sage: richcmp(3, 4, op_LT)
        True
        sage: richcmp(x, x^2, op_EQ)
        x == x^2

    The two examples above are completely equivalent to ``3 < 4``
    and ``x == x^2``. For this reason, it only makes sense in practice
    to call ``richcmp`` with a non-constant value for ``op``.

    We can write a custom ``Element`` class which shows a more
    realistic example of how to use this::

        sage: from sage.structure.element import Element
        sage: class MyElement(Element):
        ....:     def __init__(self, parent, value):
        ....:         Element.__init__(self, parent)
        ....:         self.v = value
        ....:     def _richcmp_(self, other, op):
        ....:         return richcmp(self.v, other.v, op)
        sage: P = Parent()
        sage: x = MyElement(P, 3)
        sage: y = MyElement(P, 3)
        sage: x < y
        False
        sage: x == y
        True
        sage: x > y
        False
    """
    return PyObject_RichCompare(x, y, op)


cdef slot_tp_richcompare(self, other, int op):
    """
    Function to put in the ``tp_richcompare`` slot.
    """
    return self.__richcmp__(other, op)


def richcmp_method(cls):
    """
    Class decorator to implement rich comparions using the special
    method ``__richcmp__`` (analogous to Cython) instead of the 6
    methods ``__eq__`` and friends.

    This changes the class in-place and returns the given class.

    EXAMPLES::

        sage: from sage.structure.richcmp import *
        sage: sym = {op_EQ: "==", op_NE: "!=", op_LT: "<", op_GT: ">", op_LE: "<=", op_GE: ">="}
        sage: @richcmp_method
        ....: class A(str):
        ....:     def __richcmp__(self, other, op):
        ....:         print("%s %s %s" % (self, sym[op], other))
        sage: A("left") < A("right")
        left < right
        sage: object() <= A("right")
        right >= <object object at ...>

    We can call this comparison with the usual Python special methods::

        sage: x = A("left"); y = A("right")
        sage: x.__eq__(y)
        left == right
        sage: A.__eq__(x, y)
        left == right

    Everything still works in subclasses::

        sage: class B(A):
        ....:     pass
        sage: x = B("left"); y = B("right")
        sage: x != y
        left != right
        sage: x.__ne__(y)
        left != right
        sage: B.__ne__(x, y)
        left != right

    We can override ``__richcmp__`` with standard Python rich
    comparison methods and conversely::

        sage: class C(A):
        ....:     def __ne__(self, other):
        ....:         return False
        sage: C("left") != C("right")
        False
        sage: C("left") == C("right")  # Calls __eq__ from class A
        left == right

        sage: class Base(object):
        ....:     def __eq__(self, other):
        ....:         return False
        sage: @richcmp_method
        ....: class Derived(Base):
        ....:     def __richcmp__(self, other, op):
        ....:         return True
        sage: Derived() == Derived()
        True

    TESTS::

        sage: richcmp_method(None)
        Traceback (most recent call last):
        ...
        TypeError: None is not a class

        sage: @richcmp_method
        ....: class X(object):
        ....:     def __eq__(self, other):
        ....:         pass
        ....:     def __richcmp__(self, other, op):
        ....:         pass
        Traceback (most recent call last):
        ...
        TypeError: class <class '__main__.X'> defines __eq__ which cannot be combined with @richcmp_method
    """
    if not isinstance(cls, type):
        raise TypeError(f"{cls!r} is not a class")
    cdef PyTypeObject* tp = <PyTypeObject*>cls
    tp.tp_richcompare = slot_tp_richcompare

    # Install slot wrappers in the class dict
    cdef dict D = <dict>(tp.tp_dict)
    cdef wrapperbase* slotdef
    for slotdef in richcmp_slotdef:
        name = <object>slotdef.name_strobj
        if name in D:
            raise TypeError("class %r defines %s which cannot be combined with @richcmp_method" % (cls, name))
        D[name] = PyDescr_NewWrapper(tp, slotdef, <void*>slot_tp_richcompare)

    PyType_Modified(tp)

    return cls
