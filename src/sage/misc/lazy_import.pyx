r"""
Lazy imports

This module allows one to lazily import objects into a namespace,
where the actual import is delayed until the object is actually called
or inspected. This is useful for modules that are expensive to import
or may cause circular references, though there is some overhead in its
use.

EXAMPLES::

    sage: from sage.misc.lazy_import import _lazyimport_
    sage: with _lazyimport_:
    ....:     from sage.rings.all import ZZ
    sage: type(ZZ)
    <type 'sage.misc.lazy_import.LazyImport'>
    sage: ZZ(4.0)
    4

By default, a warning is issued if a lazy import module is resolved
during Sage's startup. In case a lazy import's sole purpose is to
break a circular reference and it is known to be resolved at startup
time, one can use the ``at_startup`` option::

    sage: with _lazyimport_(at_startup=True):
    ....:     from sage.rings.all import ZZ

This option can also be used as an intermediate step toward not
importing by default a module that is used in several places, some of
which can already afford to lazy import the module but not all.

A lazy import that is marked as "at_startup" will print a message if
it is actually resolved after the startup, so that the developer knows
that (s)he can remove the flag::

    sage: ZZ
    Option ``at_startup=True`` for lazy import ZZ not needed anymore
    Integer Ring

.. SEEALSO:: :class:`LazyImportContext`, :class:`LazyImport, :func:`lazy_import`
"""

#*****************************************************************************
#       Copyright (C) 2009 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Keep OLD division semantics for Python 2 compatibility, such that
# lazy imports support old and true division.
from __future__ import absolute_import

cimport cython
from cpython.module cimport PyImport_ImportModuleLevel
from cpython.object cimport PyObject_RichCompare
from cpython.number cimport PyNumber_TrueDivide, PyNumber_Power, PyNumber_Index
from cpython.pystate cimport PyThreadState, PyThreadState_Get
cdef extern from *:
    ctypedef PyThreadState volatile_PyThreadState "volatile PyThreadState"

import os
import sys
from six.moves import builtins, cPickle as pickle
import inspect
from time import sleep

from . import sageinspect
from .lazy_import_cache import get_cache_file


cdef inline obj(x):
    if isinstance(x, LazyImport):
        return (<LazyImport>x)._get_object()
    else:
        return x


# boolean to determine whether Sage is still starting up
cdef bint startup_guard = True


cpdef finish_startup():
    """
    This function must be called exactly once at the end of the Sage
    import process

    TESTS::

        sage: from sage.misc.lazy_import import finish_startup
        sage: finish_startup()
        Traceback (most recent call last):
        ...
        AssertionError: finish_startup() must be called exactly once
    """
    global startup_guard
    assert startup_guard, 'finish_startup() must be called exactly once'
    startup_guard = False

cpdef bint is_during_startup():
    """
    Return whether Sage is currently starting up

    OUTPUT:

    Boolean

    TESTS::

        sage: from sage.misc.lazy_import import is_during_startup
        sage: is_during_startup()
        False
    """
    global startup_guard
    return startup_guard

cpdef test_fake_startup():
    """
    For testing purposes only.

    Switch the startup lazy import guard back on.

    EXAMPLES::

        sage: sage.misc.lazy_import.test_fake_startup()
        sage: from sage.misc.lazy_import import lazy_import
        sage: lazy_import('sage.rings.all', 'ZZ', 'my_ZZ')
        sage: my_ZZ(123)
        Traceback (most recent call last):
          File "/usr/local/src/sage-config/local/lib/python2.7/site-packages/sage/doctest/forker.py", line 498, in _run
            self.compile_and_execute(example, compiler, test.globs)
          File "/usr/local/src/sage-config/local/lib/python2.7/site-packages/sage/doctest/forker.py", line 861, in compile_and_execute
            exec(compiled, globs)
          File "<doctest sage.misc.lazy_import.test_fake_startup[3]>", line 1, in <module>
            my_ZZ(Integer(123))
          File "sage/misc/lazy_import.pyx", line 346, in sage.misc.lazy_import.LazyImport.__call__ (build/cythonized/sage/misc/lazy_import.c:3495)
            return self._get_object()(*args, **kwds)
          File "sage/misc/lazy_import.pyx", line 210, in sage.misc.lazy_import.LazyImport._get_object (build/cythonized/sage/misc/lazy_import.c:2139)
            raise RuntimeError(f"resolving lazy import {self._name} during startup")
        RuntimeError: resolving lazy import ZZ during startup
        sage: sage.misc.lazy_import.finish_startup()
    """
    global startup_guard
    startup_guard = True


cdef find_dict_value_by_id(D, value, guess, default):
    """
    Return the key of the dict ``D`` such that ``D[key] is value``.
    Try ``guess`` first as key.

    Return ``default`` is there is no such key.
    """
    try:
        if D[guess] is value:
            return guess
    except KeyError:
        pass
    for key in D:
        if D[key] is value:
            return key
    return default


@cython.final
cdef class LazyImport(object):
    """
    EXAMPLES::

        sage: from sage.misc.lazy_import import LazyImport
        sage: my_integer = LazyImport('sage.rings.all', 'Integer')
        sage: my_integer(4)
        4
        sage: my_integer('101', base=2)
        5
        sage: my_integer(3/2)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
    """
    cdef readonly _object  # The actual object if imported, None otherwise
    cdef _module
    cdef _name
    cdef _as_name
    cdef _namespace
    cdef bint _at_startup
    cdef _deprecation
    cdef int _level

    def __init__(self, module, name, as_name=NotImplemented,
        at_startup=False, namespace=None, deprecation=None, level=0):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: my_isprime(5)
            True
            sage: my_isprime(55)
            False
        """
        self._object = None
        self._module = module
        self._name = name
        self._as_name = name if as_name is None else as_name
        self._namespace = namespace
        self._at_startup = at_startup
        self._deprecation = deprecation
        self._level = level

    cpdef _get_object(self):
        """
        Return the wrapped object, importing it if necessary.

        OUTPUT:

        - the wrapped object

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_integer_ring = LazyImport('sage.rings.all', 'ZZ')
            sage: my_integer_ring._object is None
            True
            sage: my_integer_ring._get_object()
            Integer Ring
            sage: my_integer_ring._object is None
            False
            sage: my_integer_ring = LazyImport('sage.rings.all', 'ZZ', at_startup=True)
            sage: my_integer_ring
            Option ``at_startup=True`` for lazy import ZZ not needed anymore
            Integer Ring
        """
        if self._object is not None:
            return self._object

        if startup_guard and not self._at_startup:
            raise RuntimeError(f"resolving lazy import {self._name} during startup")
        elif self._at_startup and not startup_guard:
            print('Option ``at_startup=True`` for lazy import {0} not needed anymore'.format(self._name))

        name = self._as_name
        # Special case: if name is NotImplemented, search for it in the namespace
        if name is NotImplemented:
            name = self._name  # Reasonable default
            if self._namespace is not None:
                name = find_dict_value_by_id(self._namespace, self, name, name)

        cdef bint do_replace = 0
        if self._namespace is not None:
            do_replace = (self._namespace.get(name) is self)
            if do_replace:
                # Delete the name from the namespace, otherwise we
                # potentially break the real import below.
                del self._namespace[name]

        # Do a real import of the object
        globals = {} if self._namespace is None else self._namespace
        module = PyImport_ImportModuleLevel(self._module, globals, {},
                (self._name,), self._level)
        self._object = getattr(module, self._name)

        if self._deprecation is not None:
            from sage.misc.superseded import deprecation
            try:
                trac_number, message = self._deprecation
            except TypeError:
                trac_number = self._deprecation
                message = ('\nImporting {name} from here is deprecated. ' +
                    'If you need to use it, please import it directly from' +
                    ' {module_name}').format(name=name, module_name=self._module)
            deprecation(trac_number, message)
        # Replace the lazy import in the namespace by the actual object
        if do_replace:
            self._namespace[name] = self._object
        return self._object

    def _get_deprecation_ticket(self):
        """
        Return the ticket number of the deprecation, or 0 if this lazy
        import is not deprecated.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: H = LazyImport('sage.categories.homsets', 'Homsets')
            sage: H._get_deprecation_ticket()
            0
            sage: H = LazyImport('sage.categories.homsets', 'Homsets', deprecation=10668)
            sage: H._get_deprecation_ticket()
            10668
            sage: H = LazyImport('sage.categories.homsets', 'Homsets', deprecation=(10668, "this is deprecated"))
            sage: H._get_deprecation_ticket()
            10668
        """
        if self._deprecation is None:
            return 0
        try:
            return self._deprecation[0]
        except TypeError:
            return self._deprecation

    def _instancedoc_(self):
        """
        Return the docstring of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: my_isprime.__doc__ is is_prime.__doc__
            True

        TESTS:

        Check that :trac:`19475` is fixed::

            sage: 'A subset of the real line' in RealSet.__doc__
            True

        Check that :trac:`20626` is fixed::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell("sage.combinat.tutorial?")
            Type:            LazyImport
            String form:     <module 'sage.combinat.tutorial' from '/usr/local/src/sage-git/local/lib/python2.7/site-packages/sage/combinat/tutorial.pyc'>
            File:            /usr/local/src/sage-git/src/sage/misc/lazy_import.pyx
            Docstring:
            Introduction to combinatorics in Sage...
        """
        return sageinspect.sage_getdoc_original(self._get_object())

    def _sage_src_(self):
        """
        Returns the source of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: 'def is_prime(' in my_isprime._sage_src_()
            True
        """
        return sageinspect.sage_getsource(self._get_object())

    def _sage_argspec_(self):
        """
        Returns the argspec of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: rm = LazyImport('sage.all', 'random_matrix')
            sage: rm._sage_argspec_()
            ArgSpec(args=['ring', 'nrows', 'ncols', 'algorithm'], varargs='args', keywords='kwds', defaults=(None, 'randomize'))
        """
        return sageinspect.sage_getargspec(self._get_object())

    def __getattr__(self, attr):
        """
        Attribute lookup on self defers to attribute lookup on the
        wrapped object.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_integer = LazyImport('sage.rings.all', 'Integer')
            sage: my_integer.sqrt is Integer.sqrt
            True
        """
        return getattr(self._get_object(), attr)

    # We need to wrap all the slot methods, as they are not forwarded
    # via getattr.

    def __dir__(self):
        """
        Tab completion on self defers to completion on the wrapped
        object.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_ZZ = LazyImport('sage.rings.all', 'ZZ')
            sage: dir(my_ZZ) == dir(ZZ)
            True
        """
        return dir(self._get_object())

    def __call__(self, *args, **kwds):
        """
        Calling self calls the wrapped object.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: my_isprime(12)
            False
            sage: my_isprime(13)
            True
        """
        return self._get_object()(*args, **kwds)

    def __repr__(self):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: type(lazy_ZZ)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: lazy_ZZ
            Integer Ring
            sage: repr(lazy_ZZ)
            'Integer Ring'
        """
        return repr(self._get_object())

    def __str__(self):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: str(lazy_ZZ)
            'Integer Ring'
        """
        return str(self._get_object())

    def __unicode__(self):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: unicode(lazy_ZZ)
            u'Integer Ring'
        """
        return unicode(self._get_object())

    def __nonzero__(self):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: not lazy_ZZ
            True
        """
        return not self._get_object()

    def __hash__(self):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: hash(lazy_ZZ) == hash(1.parent())
            True
        """
        return hash(self._get_object())

    def __cmp__(left, right):
        """
        Removed by :trac:`21247` (for compatibility with Python 3)

        TESTS::

            sage: lazy_import('sage.all', ['ZZ', 'QQ'])
            sage: cmp(ZZ, QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: old-style comparisons are not supported for lazily imported objects (see https://trac.sagemath.org/ticket/21247)
        """
        raise NotImplementedError("old-style comparisons are not supported "
            "for lazily imported objects (see https://trac.sagemath.org/ticket/21247)")

    def __richcmp__(left, right, int op):
        """
        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: lazy_ZZ == RR
            False
            sage: lazy_ZZ == 1.parent()
            True
        """
        return PyObject_RichCompare(obj(left), obj(right), op)

    def __len__(self):
        """
        TESTS::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: len(version_info)
            5
        """
        return len(self._get_object())

    def __get__(self, instance, owner):
        """
        EXAMPLES:

        Here we show how to take a function in a module, and lazy
        import it as a method of a class. For the sake of this
        example, we add manually a function in sage.all::

            sage: def my_method(self): return self
            sage: sage.all.my_method = my_method

        Now we lazy import it as a method of a new class ``Foo``::

            sage: from sage.misc.lazy_import import LazyImport
            sage: class Foo:
            ....:     my_method = LazyImport('sage.all', 'my_method')

        Now we can use it as a usual method::

            sage: Foo().my_method()
            <__main__.Foo instance at ...>
            sage: Foo.my_method
            <unbound method Foo.my_method>
            sage: Foo().my_method
            <bound method Foo.my_method of <__main__.Foo instance at ...>>

        When a :class:`LazyImport` method is a method (or attribute)
        of a class, then extra work must be done to replace this
        :class:`LazyImport` object with the actual object. See the
        documentation of :meth:`_get_object` for an explanation of
        this.

        .. NOTE::

           For a :class:`LazyImport` object that appears in a class
           namespace, we need to do something special. Indeed, the
           class namespace dictionary at the time of the class
           definition is not the one that actually gets used. Thus,
           ``__get__`` needs to manually modify the class dict::

               sage: class Foo(object):
               ....:     lazy_import('sage.all', 'plot')
               sage: class Bar(Foo):
               ....:     pass
               sage: type(Foo.__dict__['plot'])
               <type 'sage.misc.lazy_import.LazyImport'>

           We access the ``plot`` method::

               sage: Bar.plot
               <unbound method Bar.plot>

           Now ``plot`` has been replaced in the dictionary of ``Foo``::

               sage: type(Foo.__dict__['plot'])
               <... 'function'>

        TESTS:

        Check that :trac:`15648` is fixed::

            sage: class A:
            ....:     Associative = LazyImport('sage.categories.magmas', 'Magmas')
            sage: type(A.Associative)
            <type 'sage.misc.classcall_metaclass.ClasscallMetaclass'>
            sage: type(A.__dict__["Associative"])
            <type 'sage.misc.classcall_metaclass.ClasscallMetaclass'>

            sage: class B(object):
            ....:     Associative = LazyImport('sage.categories.magmas', 'Magmas')
            sage: type(B.Associative)
            <type 'sage.misc.classcall_metaclass.ClasscallMetaclass'>
            sage: type(B.__dict__["Associative"])
            <type 'sage.misc.classcall_metaclass.ClasscallMetaclass'>

            sage: class C:
            ....:     Associative = LazyImport('sage.categories.magmas', 'Magmas', 'Associative')
            sage: type(C.Associative)
            <type 'sage.misc.classcall_metaclass.ClasscallMetaclass'>
            sage: type(C.__dict__["Associative"])
            <type 'sage.misc.classcall_metaclass.ClasscallMetaclass'>

            sage: class D:
            ....:     Associative = LazyImport('sage.categories.magmas', 'Magmas', 'WRONG')
            sage: type(D.Associative)
            <type 'sage.misc.classcall_metaclass.ClasscallMetaclass'>
            sage: type(D.__dict__["Associative"])
            <type 'sage.misc.lazy_import.LazyImport'>
        """
        # Don't use the namespace of the class definition
        self._namespace = None
        obj = self._get_object()

        as_name = self._as_name
        name = None
        for cls in inspect.getmro(owner):
            if as_name is NotImplemented:
                name = find_dict_value_by_id(cls.__dict__, self, as_name, None)
            elif cls.__dict__.get(as_name) is self:
                name = as_name
            else:
                continue
            if name is not None:
                setattr(cls, name, obj)
                break

        # Check whether the imported object is itself a descriptor
        try:
            get = obj.__get__
        except AttributeError:
            return obj
        else:
            return get(instance, owner)

    def __getitem__(self, key):
        """
        TESTS::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: version_info[0]
            2
        """
        return self._get_object()[key]

    def __setitem__(self, key, value):
        """
        TESTS::

            sage: sage.all.foo = list(range(10))
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo[1] = 100
            sage: print(foo)
            [0, 100, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        self._get_object()[key] = value

    def __delitem__(self, key):
        """
        TESTS::

            sage: sage.all.foo = list(range(10))
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: del foo[1]
            sage: print(foo)
            [0, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        del self._get_object()[key]

    def __iter__(self):
        """
        TESTS::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: iter(version_info)
            <iterator object at ...>
        """
        return iter(self._get_object())

    def __contains__(self, item):
        """
        TESTS::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: 2 in version_info
            True

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: 2000 not in version_info
            True
        """
        return item in self._get_object()

    def __add__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo + 1
            11
        """
        return obj(left) + obj(right)

    def __sub__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo - 1
            9
        """
        return obj(left) - obj(right)

    def __mul__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo * 2
            20
        """
        return obj(left) * obj(right)

    def __matmul__(left, right):
        """
        TESTS::

            sage: from sympy import Matrix
            sage: sage.all.foo = Matrix([[1,1],[0,1]])
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo.__matmul__(foo)
            Matrix([
            [1, 2],
            [0, 1]])
        """
        return obj(left) @ obj(right)

    def __div__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo / 2
            5
        """
        return obj(left) / obj(right)

    def __floordiv__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo  // 3
            3
        """
        return obj(left) // obj(right)

    def __truediv__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: operator.truediv(foo, 3)
            10/3
        """
        return PyNumber_TrueDivide(obj(left), obj(right))

    def __pow__(left, right, mod):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo ** 2
            100
        """
        return PyNumber_Power(obj(left), obj(right), obj(mod))

    def __mod__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo % 7
            3
        """
        return obj(left) % obj(right)

    def __lshift__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo << 3
            80
        """
        return obj(left) << obj(right)

    def __rshift__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo >> 2
            2
        """
        return obj(left) >> obj(right)

    def __and__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo & 7
            2
        """
        return obj(left) & obj(right)

    def __or__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo | 7
            15
        """
        return obj(left) | obj(right)

    def __xor__(left, right):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: foo ^^ 7
            13
        """
        return obj(left) ^ obj(right)

    def __neg__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: -foo
            -10
        """
        return -self._get_object()

    def __pos__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: +foo
            10
        """
        return +self._get_object()

    def __abs__(self):
        """
        TESTS::

            sage: sage.all.foo = -1000
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: abs(foo)
            1000
        """
        return abs(self._get_object())

    def __invert__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: ~foo
            1/10
        """
        return ~self._get_object()

    def __complex__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: complex(foo)
            (10+0j)
        """
        return complex(self._get_object())

    def __int__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: int(foo)
            10
        """
        return int(self._get_object())

    def __long__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: long(foo)
            10L
        """
        return long(self._get_object())

    def __float__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: float(foo)
            10.0
        """
        return float(self._get_object())

    def __oct__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: oct(foo)
            '12'
        """
        return oct(self._get_object())

    def __hex__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: hex(foo)
            'a'
        """
        return hex(self._get_object())

    def __index__(self):
        """
        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: range(100)[foo]
            10
        """
        return PyNumber_Index(self._get_object())

    def __copy__(self):
        """
        Support copy()

        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: copy(foo)
            10
        """
        return self._get_object()

    def __deepcopy__(self, memo=None):
        """
        Support copy()

        TESTS::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: deepcopy(foo)
            10
        """
        return self._get_object()

    def __instancecheck__(self, x):
        """
        Support ``isinstance()``.

        EXAMPLES::

            sage: lazy_import('sage.rings.rational_field', 'RationalField')
            sage: isinstance(QQ, RationalField)
            True
        """
        return isinstance(x, self._get_object())

    def __subclasscheck__(self, x):
        """
        Support ``issubclass()``.

        EXAMPLES::

            sage: lazy_import('sage.structure.parent', 'Parent')
            sage: issubclass(RationalField, Parent)
            True
        """
        return issubclass(x, self._get_object())


def lazy_import(module, names, as_=None, *,
    at_startup=False, namespace=None, overwrite=None, deprecation=None):
    """
    Create a lazy import object and inject it into the caller's global
    namespace. For the purposes of introspection and calling, this is
    like performing a lazy "from module import name" where the import
    is delayed until the object actually is used or inspected.

    INPUT:

    - ``module`` -- a string representing the module to import

    - ``names`` -- a string or list of strings representing the names to
      import from module

    - ``as_`` -- (optional) a string or list of strings representing the
      names of the objects in the importing module. This is analogous to
      ``from ... import ... as ...``.

    - ``at_startup`` -- a boolean (default: ``False``);
      whether the lazy import is supposed to be resolved at startup time

    - ``namespace`` -- the namespace where importing the names; by default,
      import the names to current namespace

    - ``deprecation`` -- (optional) if not ``None``, a deprecation warning
      will be issued when the object is actually imported;
      ``deprecation`` should be either a trac number (integer) or a
      pair ``(trac_number, message)``

    .. SEEALSO:: :mod:`sage.misc.lazy_import`, :class:`LazyImport`

    EXAMPLES::

        sage: from sage.misc.lazy_import import lazy_import
        sage: lazy_import('sage.rings.all', 'ZZ')
        sage: type(ZZ)
        <type 'sage.misc.lazy_import.LazyImport'>
        sage: ZZ(4.0)
        4
        sage: lazy_import('sage.rings.all', 'RDF', 'my_RDF')
        sage: my_RDF._get_object() is RDF
        True
        sage: my_RDF(1/2)
        0.5

        sage: lazy_import('sage.all', ['QQ', 'RR'], ['my_QQ', 'my_RR'])
        sage: my_QQ._get_object() is QQ
        True
        sage: my_RR._get_object() is RR
        True

    Upon the first use, the object is injected directly into
    the calling namespace::

        sage: lazy_import('sage.all', 'ZZ', 'my_ZZ')
        sage: my_ZZ is ZZ
        False
        sage: my_ZZ(37)
        37
        sage: my_ZZ is ZZ
        True

    We check that :func:`lazy_import` also works for methods::

        sage: class Foo(object):
        ....:     lazy_import('sage.all', 'plot')
        sage: class Bar(Foo):
        ....:     pass
        sage: type(Foo.__dict__['plot'])
        <type 'sage.misc.lazy_import.LazyImport'>
        sage: 'EXAMPLES' in Bar.plot.__doc__
        True
        sage: type(Foo.__dict__['plot'])
        <... 'function'>

    If deprecated then a deprecation warning is issued::

        sage: lazy_import('sage.all', 'Qp', 'my_Qp', deprecation=14275)
        sage: my_Qp(5)
        doctest:...: DeprecationWarning:
        Importing my_Qp from here is deprecated. If you need to use it, please import it directly from sage.all
        See http://trac.sagemath.org/14275 for details.
        5-adic Field with capped relative precision 20

    An example of deprecation with a message::

        sage: lazy_import('sage.all', 'Qp', 'my_Qp_msg', deprecation=(14275, "This is an example."))
        sage: my_Qp_msg(5)
        doctest:...: DeprecationWarning: This is an example.
        See http://trac.sagemath.org/14275 for details.
        5-adic Field with capped relative precision 20
    """
    if overwrite is not None:
        from sage.misc.superseded import deprecation
        deprecation(22755, "lazy_import(overwrite=False) is no longer supported")
    if as_ is None:
        as_ = names
    if isinstance(names, basestring):
        names = [names]
        as_ = [as_]
    else:
        names = list(names)
        as_ = list(as_)
    if namespace is None:
        namespace = inspect.currentframe().f_locals
    if "*" in names:
        ix = names.index("*")
        all = get_star_imports(module)
        names[ix:ix+1] = all
        as_[ix:ix+1] = all
    for name, alias in zip(names, as_):
        namespace[alias] = LazyImport(module, name, alias, at_startup, namespace, deprecation)


star_imports = None

def save_cache_file():
    """
    Used to save the cached * import names.

    TESTS::

        sage: import sage.misc.lazy_import
        sage: sage.misc.lazy_import.save_cache_file()
    """
    from sage.misc.misc import sage_makedirs
    from sage.misc.temporary_file import atomic_write

    global star_imports
    if star_imports is None:
        star_imports = {}
    cache_file = get_cache_file()
    cache_dir = os.path.dirname(cache_file)

    sage_makedirs(cache_dir)
    with atomic_write(cache_file) as f:
        pickle.dump(star_imports, f)


def get_star_imports(module_name):
    """
    Lookup the list of names in a module that would be imported with
    ``import *`` either via a cache or by actually importing.

    EXAMPLES::

        sage: from sage.misc.lazy_import import get_star_imports
        sage: 'get_star_imports' in get_star_imports('sage.misc.lazy_import')
        True
        sage: 'EllipticCurve' in get_star_imports('sage.schemes.all')
        True

    TESTS::

        sage: import os, tempfile
        sage: fd, cache_file = tempfile.mkstemp()
        sage: os.write(fd, 'invalid')
        7
        sage: os.close(fd)
        sage: import sage.misc.lazy_import as lazy
        sage: lazy.get_cache_file = (lambda: cache_file)
        sage: lazy.star_imports = None
        sage: lazy.get_star_imports('sage.schemes.all')
        doctest:...: UserWarning: star_imports cache is corrupted
        [...]
        sage: os.remove(cache_file)
        sage: lazy.get_cache_file = sage.misc.lazy_import_cache.get_cache_file
    """
    global star_imports
    if star_imports is None:
        star_imports = {}
        try:
            with open(get_cache_file()) as cache_file:
                star_imports = pickle.load(cache_file)
        except IOError:        # file does not exist
            pass
        except Exception:  # unpickling failed
            import warnings
            warnings.warn('star_imports cache is corrupted')
    try:
        return star_imports[module_name]
    except KeyError:
        module = PyImport_ImportModuleLevel(module_name, {}, {}, ("*",), 0)
        try:
            all = module.__all__
        except AttributeError:
            all = [key for key in dir(module) if key[0] != "_"]
        star_imports[module_name] = all
        save_cache_file()
        return all


# Add support for _instancedoc_
from sage.docs.instancedoc import instancedoc
instancedoc(LazyImport)


cdef class LazyModule:
    """
    A lazy imported module.

    All attributes of this module become lazy imports.

    INPUT:

    - ``module`` -- the name of the module

    - ``**kwds`` -- additional arguments to pass to :class:`LazyImport`

    EXAMPLES::

        sage: from sage.misc.lazy_import import LazyModule
        sage: mod = LazyModule("sage.all")
        sage: mod
        <lazy imported module 'sage.all'>
        sage: type(mod.ZZ)
        <type 'sage.misc.lazy_import.LazyImport'>

    ::

        sage: deprecated_mod = LazyModule("sage.all", deprecation=14974)
        sage: deprecated_mod.ZZ
        doctest:...: DeprecationWarning:
        Importing ZZ from here is deprecated. If you need to use it, please import it directly from sage.all
        See http://trac.sagemath.org/14974 for details.
        Integer Ring

    An indirect example using :class:`LazyImportContext`::

        sage: from sage.misc.lazy_import import _lazyimport_
        sage: with _lazyimport_ as lazy:
        ....:     mod = lazy.__import__("sage.all", {}, {}, ["ZZ"])
        sage: mod
        <lazy imported module 'sage.all'>
    """
    cdef module
    cdef kwds

    def __init__(self, module, **kwds):
        self.module = module
        self.kwds = kwds

    def __repr__(self):
        return f"<lazy imported module {self.module!r}>"

    def __getattribute__(self, attr):
        if attr == "__all__":
            return get_star_imports(self.module)
        return LazyImport(self.module, name=attr, as_name=NotImplemented, **self.kwds)


# Unique thread which is inside a LazyImportContext(). To avoid
# threading issues, at most one thread can do this. Due to the Python
# GIL, there are no race conditions in accessing this global pointer.
cdef volatile_PyThreadState* lazyimport_thread = NULL

class LazyImportContext(object):
    """
    Context in which all imports become lazy imports.

    .. WARNING::

        You should use this only for module-level imports, not for
        imports inside a function. This is because the namespace where
        this import ends up is assumed to be the globals.

    EXAMPLES::

        sage: from sage.misc.lazy_import import _lazyimport_
        sage: with _lazyimport_:
        ....:     from sage.all import ZZ
        sage: type(ZZ)
        <type 'sage.misc.lazy_import.LazyImport'>
        sage: ZZ("0x2a")
        42
        sage: type(ZZ)
        <type 'sage.rings.integer_ring.IntegerRing_class'>

    You can pass additional keywords to :class:`LazyImport`::

        sage: with _lazyimport_(deprecation=14974):
        ....:     from sage.all import QQ
        sage: QQ
        doctest:...: DeprecationWarning:
        Importing QQ from here is deprecated. If you need to use it, please
        import it directly from sage.all
        See http://trac.sagemath.org/14974 for details.
        Rational Field

    Star imports work::

        sage: with _lazyimport_:
        ....:     from sage.structure.parent import *
        sage: type(Parent)
        <type 'sage.misc.lazy_import.LazyImport'>

    This works also with ``import ... as ...`` in Python and in Cython::

        sage: with _lazyimport_:
        ....:     from sage.all import ZZ as my_ZZ
        sage: type(my_ZZ)
        <type 'sage.misc.lazy_import.LazyImport'>
        sage: my_ZZ("0x2a")
        42
        sage: type(my_ZZ)
        <type 'sage.rings.integer_ring.IntegerRing_class'>

        sage: cython(  # long time
        ....: '''
        ....: from sage.misc.lazy_import import _lazyimport_
        ....: with _lazyimport_:
        ....:     from sage.all import ZZ as my_ZZ
        ....: print(type(my_ZZ))
        ....: fortytwo = my_ZZ("0x2a")
        ....: print(type(my_ZZ))
        ....: ''')
        <type 'sage.misc.lazy_import.LazyImport'>
        <type 'sage.rings.integer_ring.IntegerRing_class'>
    """
    def __init__(self, kwds={}):
        """
        TESTS::

            sage: from sage.misc.lazy_import import LazyImportContext
            sage: ctx = LazyImportContext({"foo":"bar"}); ctx
            <sage.misc.lazy_import.LazyImportContext object at ...>
            sage: ctx.kwds
            {'foo': 'bar'}
        """
        self.kwds = kwds

    def __enter__(self):
        """
        Change the builtin ``__import__`` function to ``self.__import__``.

        EXAMPLES::

            sage: from sage.misc.lazy_import import _lazyimport_
            sage: with _lazyimport_:
            ....:     print(__import__)
            <bound method LazyImportContext.__import__ of <sage.misc.lazy_import.LazyImportContext object at ...>>

        Nesting is illegal::

            sage: with _lazyimport_:
            ....:     with _lazyimport_:
            ....:         pass
            Traceback (most recent call last):
            ...
            RuntimeError: nesting "with _lazyimport_" contexts is not allowed
        """
        global lazyimport_thread
        while lazyimport_thread is not NULL:
            if lazyimport_thread is PyThreadState_Get():
                # Current thread is already doing lazy imports.
                raise RuntimeError('nesting "with _lazyimport_" contexts is not allowed')
            # A different thread is doing lazy imports, yield the GIL
            # for 20ms.
            sleep(0.020)

        lazyimport_thread = PyThreadState_Get()
        self.old_import = builtins.__import__
        builtins.__import__ = self.__import__
        return self

    def __exit__(self, *args):
        """
        Change the builtin ``__import__`` function back to what it was
        before.

        EXAMPLES::

            sage: from sage.misc.lazy_import import _lazyimport_
            sage: with _lazyimport_:
            ....:     pass
            sage: print(__import__)
            <built-in function __import__>
            sage: with _lazyimport_:
            ....:     raise Exception
            Traceback (most recent call last):
            ...
            Exception
            sage: print(__import__)
            <built-in function __import__>

        We get a message if ``__import__`` was changed inside the
        context::

            sage: from six.moves import builtins
            sage: with _lazyimport_:
            ....:     builtins.__import__ = None
            WARNING: __import__ was changed inside a "with _lazyimport_" context to None (possibly by a different thread)
        """
        global lazyimport_thread
        if builtins.__import__ != self.__import__:
            print(f'WARNING: __import__ was changed inside a "with _lazyimport_" context to {builtins.__import__} (possibly by a different thread)')
        builtins.__import__ = self.old_import
        del self.old_import
        lazyimport_thread = NULL

    def __call__(self, **kwds):
        """
        Return a new :class:`LazyImportContext` which will apply
        ``kwds`` when creating a lazy import.

        EXAMPLES::

            sage: from sage.misc.lazy_import import _lazyimport_
            sage: _lazyimport_(foo="bar")(hello="yo")(foo="foo")().kwds
            {'foo': 'foo', 'hello': 'yo'}
        """
        newkwds = dict(self.kwds)
        newkwds.update(kwds)
        return type(self)(newkwds)

    def __import__(self, module, globals=None, locals=None, fromlist=(), level=0):
        """
        Replacement function for ``builtins.__import__``.

        It is not allowed to call this outside of a ``with _lazyimport_``
        context. If one thread does ``with _lazyimport_``, then other
        threads doing imports will fall back to the previous (non-lazy)
        ``__import__`` function.

        TESTS::

            sage: from sage.misc.lazy_import import _lazyimport_
            sage: with _lazyimport_ as lazy:
            ....:     print(lazy.__import__("sage.all", {}, {}, ["ZZ"]))
            <lazy imported module 'sage.all'>
            sage: _lazyimport_.__import__("sage.all", {}, {}, ["ZZ"])
            Traceback (most recent call last):
            ...
            RuntimeError: __import__ called outside "with _lazyimport_" context
            sage: with _lazyimport_:
            ....:     import sage.all
            Traceback (most recent call last):
            ...
            RuntimeError: lazy import only works with 'from ... import ...' imports, not with 'import ...'

        We test 100 threads each doing an import, half of them lazy and
        half of them a real import. Note that this will not cause a
        lot of CPU activity due to the Python GIL. ::

            sage: from sage.misc.lazy_import import LazyImport, _lazyimport_
            sage: from threading import Thread
            sage: from time import sleep
            sage: def do_import():
            ....:     sleep(0.001)  # Yield the GIL
            ....:     from sage.structure.sage_object import SageObject
            ....:     t = type(SageObject)
            ....:     if t is not type:
            ....:         print("SageObject should be type, got {!r}".format(t))
            sage: def do_lazy_import():
            ....:     with _lazyimport_:
            ....:         sleep(0.001)  # Yield the GIL
            ....:         from sage.structure.sage_object import SageObject
            ....:     t = type(SageObject)
            ....:     if t is not LazyImport:
            ....:         print("SageObject should be LazyImport, got {!r}".format(t))
            sage: threads = [Thread(target=do_import if i%2 else do_lazy_import) for i in range(100)]
            sage: for T in threads:
            ....:     T.start()
            sage: for T in threads:
            ....:     T.join()
        """
        # If __import__ was called from a different thread, call
        # old_import() instead.
        if PyThreadState_Get() is not lazyimport_thread:
            try:
                old = self.old_import
            except AttributeError:
                raise RuntimeError('__import__ called outside "with _lazyimport_" context')
            return old(module, globals, locals, fromlist, level)

        if not fromlist:
            raise RuntimeError("lazy import only works with 'from ... import ...' imports, not with 'import ...'")
        return LazyModule(module, namespace=globals, level=level, **self.kwds)


_lazyimport_ = LazyImportContext()
