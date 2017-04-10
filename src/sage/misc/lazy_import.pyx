r"""
Lazy imports

This module allows one to lazily import objects into a namespace,
where the actual import is delayed until the object is actually called
or inspected. This is useful for modules that are expensive to import
or may cause circular references, though there is some overhead in its
use.

EXAMPLES::

    sage: from sage.misc.lazy_import import lazy_import
    sage: lazy_import('sage.rings.all', 'ZZ')
    sage: type(ZZ)
    <class 'sage.misc.lazy_import.LazyImport'>
    sage: ZZ(4.0)
    4

By default, a warning is issued if a lazy import module is resolved
during Sage's startup. In case a lazy import's sole purpose is to
break a circular reference and it is known to be resolved at startup
time, one can use the ``at_startup`` option::

    sage: lazy_import('sage.rings.all', 'ZZ', at_startup=True)

This option can also be used as an intermediate step toward not
importing by default a module that is used in several places, some of
which can already afford to lazy import the module but not all.

A lazy import that is marked as "at_startup" will print a message if
it is actually resolved after the startup, so that the developer knows
that (s)he can remove the flag::

    sage: ZZ
    Option ``at_startup=True`` for lazy import ZZ not needed anymore
    Integer Ring

.. SEEALSO:: :func:`lazy_import`, :class:`LazyImport`

AUTHOR:

 - Robert Bradshaw
"""

#*****************************************************************************
#       Copyright (C) 2009 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

import importlib
import inspect
import os
from six.moves import cPickle as pickle
import operator
import sys
import traceback

from . import sageinspect

from .lazy_import_cache import get_cache_file

from lazy_object_proxy import Proxy


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
        -------------------------------------------------------------------------------
        Resolving lazy import ZZ during startup
        Calling stack:
        ...
        -------------------------------------------------------------------------------
        123
        sage: sage.misc.lazy_import.finish_startup()
    """
    global startup_guard
    startup_guard = True


class LazyImport(Proxy):
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

    TESTS::

        sage: my_integer_ring = LazyImport('sage.rings.all', 'ZZ', at_startup=True)
        sage: my_integer_ring
        Option ``at_startup=True`` for lazy import ZZ not needed anymore
        Integer Ring

        Check that :trac:`19475` is fixed::

            sage: 'A subset of the real line' in RealSet.__doc__
            True

        Test ``__getattr__``::

            sage: my_integer = LazyImport('sage.rings.all', 'Integer')
            sage: my_integer.sqrt is Integer.sqrt
            True

        Test ``__dir__``::

            sage: my_ZZ = LazyImport('sage.rings.all', 'ZZ')
            sage: dir(my_ZZ) == dir(ZZ)
            True

        Test ``__str__``::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: str(lazy_ZZ)
            'Integer Ring'

        Test ``__unicode__``::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: unicode(lazy_ZZ)
            u'Integer Ring'

        Test ``__nonzero__``::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: not lazy_ZZ
            False

        Test ``__hash__``::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: hash(lazy_ZZ) == hash(1.parent())
            True

        Test ``__richcmp__``::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: lazy_ZZ == RR
            False
            sage: lazy_ZZ == 1.parent()
            True

        Test ``__len__``::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: len(version_info)
            5

        Test ``__getitem__``::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: version_info[0]
            2

        Test ``__setitem__``::

            sage: sage.all.foo = list(range(10))
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo[1] = 100
            sage: print(foo)
            [0, 100, 2, 3, 4, 5, 6, 7, 8, 9]

        Test ``__delitem__``::

            sage: sage.all.foo = list(range(10))
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: del foo[1]
            sage: print(foo)
            [0, 2, 3, 4, 5, 6, 7, 8, 9]

        Test ``__iter__``::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: iter(version_info)
            <iterator object at ...>

        Test ``__contains__``::

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: 2 in version_info
            True

            sage: lazy_import('sys', 'version_info')
            sage: type(version_info)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: 2000 not in version_info
            True

        Test ``__add__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo + 1
            11

        Test ``__sub__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo - 1
            9

        Test ``__mul__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo * 2
            20

        Test ``__div__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo / 2
            5

        Test ``__floordiv__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo  // 3
            3

        Test ``__truediv__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: operator.truediv(foo, 3)
            10/3

        Test ``__pow__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo ** 2
            100

        Test ``__mod__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo % 7
            3

        Test ``__lshift__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo << 3
            80

        Test ``__rshift__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo >> 2
            2

        Test ``__and__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo & 7
            2

        Test ``__or__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo | 7
            15

        Test ``__xor__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: foo ^^ 7
            13

        Test ``__neg__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: -foo
            -10

        Test ``__pos__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: +foo
            10

        Test ``__abs__``::

            sage: sage.all.foo = -1000
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: abs(foo)
            1000

        Test ``__invert__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: ~foo
            1/10

        Test ``__complex__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: complex(foo)
            (10+0j)

        Test ``__init__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: int(foo)
            10

        Test ``__long__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: long(foo)
            10L

        Test ``__float__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: float(foo)
            10.0

        Test ``__oct__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: oct(foo)
            '12'

        Test ``__hex__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: hex(foo)
            'a'

        Test ``__index__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: range(100)[foo]
            10

        Test ``__copy__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: isinstance(copy(foo), LazyImport)
            False

        Test ``__deepcopy__``::

            sage: sage.all.foo = 10
            sage: lazy_import('sage.all', 'foo')
            sage: type(foo)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: isinstance(deepcopy(foo), LazyImport)
            False
    """

    __deprecation_ticket__ = None
    __initialized_wrapped__ = False

    def __init__(self, module, name, as_name=None, namespace=None, at_startup=False, deprecation=None):
        try:
            trac_number, message = deprecation
        except TypeError:
            trac_number = deprecation
            message = None

        self.__deprecation_ticket__ = trac_number

        def factory(owner=None):
            # For future calls to self.__factory__ just return the already
            # initalized wrapped object; this is used in particular in __get__
            # below
            if self.__initialized_wrapped__:
                return self.__wrapped__

            if startup_guard and not at_startup:
                print('-' * 79)
                print('Resolving lazy import {0} during startup'.format(name))
                print('Calling stack:')
                traceback.print_stack(None, None, sys.stdout)
                print('-' * 79)
            elif at_startup and not startup_guard:
                print('Option ``at_startup=True`` for lazy import {0} not needed anymore'.format(name))

            try:
                obj = getattr(importlib.import_module(module), name)
            except AttributeError as exc:
                # <module>.<name> might itself be a submodule if <module> is
                # a package
                try:
                    obj = importlib.import_module('%s.%s' % (module, name))
                except ImportError:
                    raise exc

            alias = as_name or name

            if deprecation is not None:
                from sage.misc.superseded import deprecation as _deprecation
                msg = message
                if msg is None:
                    msg = ('\nImporting {name} from here is deprecated. ' +
                        'If you need to use it, please import it directly from' +
                        ' {module_name}').format(name=alias, module_name=module)

                _deprecation(trac_number, msg)

            if owner is None:
                if namespace and namespace[alias] is self:
                    namespace[alias] = obj
            else:
                for cls in inspect.getmro(owner):
                    if cls.__dict__.get(alias, None) is self:
                        setattr(cls, alias, obj)
                        break

            self.__initialized_wrapped__ = True
            return obj

        super(LazyImport, self).__init__(factory)

    @property
    def __class__(self):
        """
        The ``Proxy`` base class overrides the normal behavior of the
        ``__class__`` attribute so it always returns the class of the wrapped
        object.

        This creates a problem, however, especially at import time as there is
        some code that runs at import that performs ``isinstance`` checks on
        objects that are proxied, and the ``isinstance`` implementation
        accesses the ``__class__`` attribute.

        We override ``__class__`` here so that it returns the proxy class
        itself ''until'' the wrapped object has been initialized through any
        other means.  After that it returns the ``__class__`` attribute of the
        wrapped type.

        TESTS::

            sage: real_ZZ = ZZ
            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: lazy_ZZ.__class__
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: isinstance(lazy_ZZ, type(real_ZZ))
            False
            sage: lazy_ZZ
            Integer Ring
            sage: lazy_ZZ.__class__ == type(real_ZZ)
            True
            sage: isinstance(lazy_ZZ, type(real_ZZ))
            True
        """

        if not self.__initialized_wrapped__:
            return type(self)

        return self.__wrapped__.__class__

    def __repr__(self):
        """
        The ``Proxy`` base class's ``__repr__`` does ''not'' resolve the
        wrapped object by default.

        TESTS::

            sage: lazy_import('sage.all', 'ZZ'); lazy_ZZ = ZZ
            sage: type(lazy_ZZ)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: lazy_ZZ
            Integer Ring
            sage: repr(lazy_ZZ)
            'Integer Ring'
        """

        return repr(self.__wrapped__)

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

    def __instancecheck__(self, x):
        """
        Support ``isinstance()``.

        EXAMPLES::

            sage: lazy_import('sage.rings.rational_field', 'RationalField')
            sage: isinstance(QQ, RationalField)
            True
        """

        return isinstance(x, self.__wrapped__)

    def __subclasscheck__(self, x):
        """
        Support ``issubclass()``.

        EXAMPLES::

            sage: lazy_import('sage.structure.parent', 'Parent')
            sage: issubclass(RationalField, Parent)
            True
        """

        return issubclass(x, self.__wrapped__)

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
        """

        wrapped = self.__wrapped__ = self.__factory__(owner)

        if hasattr(wrapped, '__get__'):
            return wrapped.__get__(instance, owner)

        return wrapped

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
        """

        return sageinspect.sage_getdoc_original(self.__wrapped__)

    def _sage_src_(self):
        """
        Returns the source of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: my_isprime = LazyImport('sage.all', 'is_prime')
            sage: 'def is_prime(' in my_isprime._sage_src_()
            True
        """

        return sageinspect.sage_getsource(self.__wrapped__)

    def _sage_argspec_(self):
        """
        Returns the argspec of the wrapped object for introspection.

        EXAMPLES::

            sage: from sage.misc.lazy_import import LazyImport
            sage: rm = LazyImport('sage.all', 'random_matrix')
            sage: rm._sage_argspec_()
            ArgSpec(args=['ring', 'nrows', 'ncols', 'algorithm'], varargs='args', keywords='kwds', defaults=(None, 'randomize'))
        """

        return sageinspect.sage_getargspec(self.__wrapped__)



def lazy_import(module, names, _as=None, namespace=None, bint overwrite=True, at_startup=False, deprecation=None):
    """
    Create a lazy import object and inject it into the caller's global
    namespace. For the purposes of introspection and calling, this is
    like performing a lazy "from module import name" where the import
    is delayed until the object actually is used or inspected.

    INPUT:

    - ``module`` -- a string representing the module to import

    - ``names`` -- a string or list of strings representing the names to
      import from module

    - ``_as`` -- (optional) a string or list of strings representing the
      aliases of the names imported

    - ``namespace`` -- the namespace where importing the names; by default,
      import the names to current namespace

    - ``overwrite`` -- (default: ``True``) if set to ``True`` and a name is
      already in the namespace, overwrite it with the lazy_import-ed name

    - ``at_startup`` -- a boolean (default: ``False``);
      whether the lazy import is supposed to be resolved at startup time

    - ``deprecation`` -- (optional) if not ``None``, a deprecation warning
      will be issued when the object is actually imported;
      ``deprecation`` should be either a trac number (integer) or a
      pair ``(trac_number, message)``

    .. SEEALSO:: :mod:`sage.misc.lazy_import`, :class:`LazyImport`

    EXAMPLES::

        sage: from sage.misc.lazy_import import lazy_import
        sage: lazy_import('sage.rings.all', 'ZZ')
        sage: type(ZZ)
        <class 'sage.misc.lazy_import.LazyImport'>
        sage: ZZ(4.0)
        4
        sage: lazy_import('sage.rings.all', 'RDF', 'my_RDF')
        sage: my_RDF.__wrapped__ is RDF
        True
        sage: my_RDF(1/2)
        0.5

        sage: lazy_import('sage.all', ['QQ', 'RR'], ['my_QQ', 'my_RR'])
        sage: my_QQ.__wrapped__ is QQ
        True
        sage: my_RR.__wrapped__ is RR
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
        <class 'sage.misc.lazy_import.LazyImport'>
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
    if _as is None:
        _as = names
    if isinstance(names, str):
        names = [names]
        _as = [_as]
    if namespace is None:
        namespace = inspect.currentframe().f_locals
    if "*" in names:
        ix = names.index("*")
        names[ix:ix+1] = get_star_imports(module)
        _as[ix:ix+1] = [None] * (len(names) - len(_as) + 1)
    for name, alias in zip(names, _as):
        if not overwrite and (alias or name) in namespace:
            continue
        namespace[alias or name] = LazyImport(module, name, alias, namespace, at_startup, deprecation)


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
    Lookup the list of names in a module that would be imported with "import \\*"
    either via a cache or actually importing.

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
        module = __import__(module_name, {}, {}, ["*"])
        if hasattr(module, "__all__"):
            all = module.__all__
        else:
            all = [key for key in dir(module) if key[0] != "_"]
        star_imports[module_name] = all
        return all


# Add support for _instancedoc_
from sage.docs.instancedoc import instancedoc
instancedoc(LazyImport)
