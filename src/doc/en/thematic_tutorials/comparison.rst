.. -*- coding: utf-8 -*-

.. _comparison:

**************************************
Comparisons in Python, Cython and Sage
**************************************

.. contents::
   :depth: 3

The goal of this tutorial is to explain how the various comparison
methods and functions like ``__lt__``, ``__cmp__``, ``__richcmp__``,
``cmp()``, ... work in Python and Cython and what Sage add on top of
this.

Comparisons are a major difference between Python 2 and Python 3, this
document only treats Python 2.

Comparisons for Python classes
==============================

.. linkall

This section deals with comparisons for new-style Python classes, that
is classes defined using ``class MyClass(object):`` (or indirectly
inheriting from ``object``). Old-style classes like ``class MyClass:``
implement comparison differently and are not discussed here.
For builtin classes or Cython extension types, see :ref:`compare-cython`.

A Python class has 2 largely independent mechanisms to implement
comparison: "three-way comparison" using ``__cmp__`` and
"rich comparison" using the ``__lt__``, ``__le__``, ``__gt__``,
``__ge__``, ``__eq__`` and ``__ne__`` methods.

Three-way comparisons
---------------------

Three-way comparison is implemented using the ``cmp()`` function and
the ``__cmp__`` special method. In early versions of Python, this was
the only comparison available. In Python 3, three-way comparisons are
no longer supported.

The call ``x.__cmp__(y)`` must return -1, 0 or 1 according to whether
`x < y`, `x == y` or `x > y`::

    sage: 3.__cmp__(9)
    -1
    sage: (4/5).__cmp__(4/5)
    0
    sage: (6.28).__cmp__(3.14)
    1

Three-way comparison is usually invoked using ``cmp(x, y)``. This is
what Python does with this call, assuming that ``x`` and ``y`` are
instances of (possibly different) new-style classes.

.. NOTE::

    When we say that ``x.method(y)`` is implemented, we mean that ``x``
    has a method :meth:`method` and that the call ``x.method(y)`` does
    not return the value ``NotImplemented`` (this is not the same as
    raising ``NotImplementedError``, which is treated like other
    exceptions).

1. if ``x is y``, return 0.
2. compare using ``__cmp___``: if ``x`` and ``y`` have the same type
   and that type has a ``__cmp__`` method, then ``x.__cmp__(y)`` is
   called.

   - If ``x.__cmp__(y)`` returns an integer, return the sign
     (-1, 0 or 1) of this integer as the result of ``cmp(x, y)``.
   - If ``x.__cmp__(y)`` is not implemented, call
     ``y.__cmp__(x)`` and handle its result analogously.
   - If neither ``x.__cmp__(y)`` and ``y.__cmp__(x)`` are implemented,
     return ``cmp(id(x), id(y))``.
   - If ``x.__cmp__(y)`` returns something which is not an integer and
     which is not ``NotImplemented``, raise an exception.

3. compare using ``__eq__`` (if ``type(y)`` is a subtype of ``type(x)``,
   reverse the order of the following 2 steps):

   - if ``x.__eq__(y)`` is implemented, convert the result to boolean.
     If true, return 0. If false, go to step 4.
   - if ``y.__eq__(x)`` is implemented, convert the result to boolean.
     If true, return 0. If false, go to step 4.

4. similar as step 3, but calling ``x.__lt__(y)`` and ``y.__gt__(x)``.
5. similar as step 3, but calling ``x.__gt__(y)`` and ``y.__lt__(x)``.
6. if either ``x`` or ``y`` has a ``__cmp__`` method, compare using
   ``__cmp___``, exactly like step 2.
7. if ``coerce(x, y)`` succeeds, let ``(v, w) = coerce(x, y)``.
   If ``v`` and ``w`` have a ``__cmp__`` method, execute step 2 with
   ``v`` and ``w``.
8. if ``x`` and ``y`` are of the same type, return
   ``cmp(id(x), id(y))``.
9. if both ``x`` and ``y`` are "numbers" (i.e. have an ``__int__``
   or ``__float__`` method), return ``cmp(id(type(x)), id(type(y)))``.
9. if neither ``x`` nor ``y`` is a number, return
   ``cmp(type(x).__name__, type(y).__name__)``.
9. if one of ``x`` or ``y`` is a number, the number is considered
   smaller.

If any of the above raise an exception, then that exception will
be propagated to the ``cmp()`` call.

Note that the exact code path is even more complicated than the
above, but assuming that the methods mentioned above are idempotent
(e.g calling ``x.__cmp__(y)`` twice returns the same result twice
without side effects), the answer should be consistent with the above.

.. _compare-cython:

Comparisons in Cython
=====================

While the title of this section refers to Cython, it is really about
comparisons for more general objects than Python classes: it is about
Cython extension types (i.e. a ``cdef class``), built-in Python types
(like ``dict``) and types implemented in C (like the ``numpy`` types).

Rich comparisons
must be implemented using a ``__richcmp__`` method.

The ``__cmp__`` method *requires* that both arguments have the same
Python type::

    sage: 1.__cmp__(1.0)
    Traceback (most recent call last):
    ...
    TypeError: sage.rings.integer.Integer.__cmp__(x,y) requires y to be a 'sage.rings.integer.Integer', not a 'sage.rings.real_mpfr.RealLiteral'

Conversely, when writing a ``__cmp__`` method, you may assume that the
types are the same.

