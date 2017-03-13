"""
Handling Superseded Functionality

The main mechanism in Sage to deal with superseded functionality is to
add a deprecation warning. This will be shown once, the first time
that the deprecated function is called.

Note that all doctests in the following use the trac ticket number
:trac:`13109`, which is where this mandatory argument to
:func:`deprecation` was introduced.

Functions and classes
---------------------
"""



########################################################################
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################
from __future__ import print_function
from six import iteritems

from warnings import warn
import inspect

from sage.misc.lazy_attribute import lazy_attribute


def _check_trac_number(trac_number):
    """
    Check that the argument is likely to be a valid trac issue number.

    INPUT:

    - ``trac_number`` -- anything.

    OUTPUT:

    This function returns nothing. A ``ValueError`` or ``TypeError`` is
    raised if the argument can not be a valid trac number.

    EXAMPLES::

        sage: from sage.misc.superseded import _check_trac_number
        sage: _check_trac_number(1)
        sage: _check_trac_number(0)
        Traceback (most recent call last):
        ...
        ValueError: 0 is not a valid trac issue number
        sage: _check_trac_number(int(10))
        sage: _check_trac_number(long(1000))
        sage: _check_trac_number(10.0)
        Traceback (most recent call last):
        ...
        TypeError: 10.0000000000000 is not a valid trac issue number
        sage: _check_trac_number('10')
        Traceback (most recent call last):
        ...
        TypeError: '10' is not a valid trac issue number
    """
    try:
        trac_number = trac_number.__index__()
    except Exception:
        raise TypeError('%r is not a valid trac issue number'%trac_number)
    if trac_number <= 0:
        raise ValueError('%r is not a valid trac issue number'%trac_number)

def deprecation(trac_number, message, stacklevel=4):
    r"""
    Issue a deprecation warning.

    INPUT:

    - ``trac_number`` -- integer. The trac ticket number where the
      deprecation is introduced.

    - ``message`` -- string. An explanation why things are deprecated
      and by what it should be replaced.

    - ``stack_level`` -- (default: ``4``) an integer. This is passed on to
      :func:`warnings.warn`.

    EXAMPLES::

        sage: def foo():
        ....:  sage.misc.superseded.deprecation(13109, 'the function foo is replaced by bar')
        sage: foo()
        doctest:...: DeprecationWarning: the function foo is replaced by bar
        See http://trac.sagemath.org/13109 for details.

    .. SEEALSO::

        :func:`experimental`,
        :func:`warning`.
    """
    warning(trac_number, message, DeprecationWarning, stacklevel)

def warning(trac_number, message, warning_class=Warning, stacklevel=3):
    r"""
    Issue a warning.

    INPUT:

    - ``trac_number`` -- integer. The trac ticket number where the
      deprecation is introduced.

    - ``message`` -- string. An explanation what is going on.

    - ``warning_class`` -- (default: ``Warning``) a class inherited
      from a Python :class:`~exceptions.Warning`.

    - ``stack_level`` -- (default: ``3``) an integer. This is passed on to
      :func:`warnings.warn`.

    EXAMPLES::

        sage: def foo():
        ....:     sage.misc.superseded.warning(
        ....:         99999,
        ....:         'The syntax will change in future.',
        ....:         FutureWarning)
        sage: foo()
        doctest:...: FutureWarning: The syntax will change in future.
        See http://trac.sagemath.org/99999 for details.

    .. SEEALSO::

        :func:`deprecation`,
        :func:`experimental`,
        :class:`exceptions.Warning`.
    """
    _check_trac_number(trac_number)
    message += '\n'
    message += 'See http://trac.sagemath.org/'+ str(trac_number) + ' for details.'
    # Stack level 3 to get the line number of the code which called
    # the deprecated function which called this function.
    warn(message, warning_class, stacklevel)

def experimental_warning(trac_number, message, stacklevel=4):
    r"""
    Issue a warning that the functionality or class is experimental
    and might change in future.

    INPUT:

    - ``trac_number`` -- an integer. The trac ticket number where the
      experimental functionality was introduced.

    - ``message`` -- a string. An explanation what is going on.

    - ``stack_level`` -- (default: ``4``) an integer. This is passed on to
      :func:`warnings.warn`.

    EXAMPLES::

        sage: def foo():
        ....:    sage.misc.superseded.experimental_warning(
        ....:        66666, 'This function is experimental and '
        ....:               'might change in future.')
        sage: foo()
        doctest:...: FutureWarning: This function is experimental and
        might change in future.
        See http://trac.sagemath.org/66666 for details.

    .. SEEALSO::

        :class:`mark_as_experimental`,
        :func:`warning`,
        :func:`deprecation`.
    """
    warning(trac_number, message, FutureWarning, stacklevel)


class experimental(object):
    def __init__(self, trac_number, stacklevel=4):
        """
        A decorator which warns about the experimental/unstable status of
        the decorated class/method/function.

        INPUT:

        - ``trac_number`` -- an integer. The trac ticket number where this
          code was introduced.

        - ``stack_level`` -- (default: ``4``) an integer. This is passed on to
          :func:`warnings.warn`.

        EXAMPLES::

            sage: @sage.misc.superseded.experimental(trac_number=79997)
            ....: def foo(*args, **kwargs):
            ....:     print("{} {}".format(args, kwargs))
            sage: foo(7, what='Hello')
            doctest:...: FutureWarning: This class/method/function is
            marked as experimental. It, its functionality or its
            interface might change without a formal deprecation.
            See http://trac.sagemath.org/79997 for details.
            (7,) {'what': 'Hello'}

        ::

            sage: class bird(SageObject):
            ....:     @sage.misc.superseded.experimental(trac_number=99999)
            ....:     def __init__(self, *args, **kwargs):
            ....:         print("piep {} {}".format(args, kwargs))
            sage: _ = bird(99)
            doctest:...: FutureWarning: This class/method/function is
            marked as experimental. It, its functionality or its
            interface might change without a formal deprecation.
            See http://trac.sagemath.org/99999 for details.
            piep (99,) {}

        TESTS:

        The following test works together with the doc-test for
        :meth:`__experimental_self_test` to demonstrate that warnings are issued only
        once, even in doc-tests (see :trac:`20601`).
        ::

            sage: from sage.misc.superseded import __experimental_self_test
            sage: _ = __experimental_self_test("A")
            doctest:...: FutureWarning: This class/method/function is
            marked as experimental. It, its functionality or its
            interface might change without a formal deprecation.
            See http://trac.sagemath.org/88888 for details.
            I'm A

        .. SEEALSO::

            :func:`experimental`,
            :func:`warning`,
            :func:`deprecation`.
        """
        self.trac_number = trac_number
        self.stacklevel = stacklevel

    def __call__(self, func):
        """
        Print experimental warning.

        INPUT:

        - ``func`` -- the function to decorate.

        OUTPUT:

        The wrapper to this function.

        TESTS::

            sage: def foo(*args, **kwargs):
            ....:     print("{} {}".format(args, kwargs))
            sage: from sage.misc.superseded import experimental
            sage: ex_foo = experimental(trac_number=99399)(foo)
            sage: ex_foo(3, what='Hello')
            doctest:...: FutureWarning: This class/method/function is
            marked as experimental. It, its functionality or its
            interface might change without a formal deprecation.
            See http://trac.sagemath.org/99399 for details.
            (3,) {'what': 'Hello'}
        """
        from sage.misc.decorators import sage_wraps
        @sage_wraps(func)
        def wrapper(*args, **kwds):
            from sage.misc.superseded import experimental_warning
            if not wrapper._already_issued:
                experimental_warning(self.trac_number,
                            'This class/method/function is marked as '
                            'experimental. It, its functionality or its '
                            'interface might change without a '
                            'formal deprecation.',
                            self.stacklevel)
                wrapper._already_issued = True
            return func(*args, **kwds)
        wrapper._already_issued = False

        return wrapper

from sage.structure.sage_object import SageObject
class __experimental_self_test(SageObject):
    r"""
    This is a class only to demonstrate with a doc-test that the @experimental
    decorator only issues a warning message once (see :trac:`20601`).

    The test below does not issue a warning message because that warning has
    already been issued by a previous doc-test in the @experimental code. Note
    that this behaviour can not be demonstrated within a single documentation
    string: Sphinx will itself supress multiple issued warnings.
    
    TESTS::

        sage: from sage.misc.superseded import __experimental_self_test
        sage: _ = __experimental_self_test("B")
        I'm B
    """
    @experimental(trac_number=88888)
    def __init__(self, x):
        print("I'm " + x)


class DeprecatedFunctionAlias(object):
    """
    A wrapper around methods or functions which automatically print
    the correct deprecation message. See
    :func:`deprecated_function_alias`.

    AUTHORS:

    - Florent Hivert (2009-11-23), with the help of Mike Hansen.
    - Luca De Feo (2011-07-11), printing the full module path when different from old path
    """
    def __init__(self, trac_number, func, module, instance = None, unbound = None):
        r"""
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: g = deprecated_function_alias(13109, number_of_partitions)
            sage: from sage.misc.superseded import deprecated_function_alias
            sage: g.__doc__
            'Deprecated: Use :func:`number_of_partitions` instead.\nSee :trac:`13109` for details.\n\n'
        """
        _check_trac_number(trac_number)
        try:
            self.__dict__.update(func.__dict__)
        except AttributeError:
            pass # Cython classes don't have __dict__
        self.func = func
        self.trac_number  = trac_number
        self.instance = instance # for use with methods
        self.unbound = unbound
        self.__module__ = module
        if isinstance(func, type(deprecation)):
            sphinxrole = "func"
        else:
            sphinxrole = "meth"
        doc = 'Deprecated: '
        doc += 'Use :' + sphinxrole + ':`' + self.func.__name__ + '` instead.\n'
        doc += 'See :trac:`' + str(self.trac_number) + '` for details.\n\n'
        self.__doc__ = doc

    @lazy_attribute
    def __name__(self):
        """
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: g = deprecated_function_alias(13109, number_of_partitions)
            sage: g.__name__
            'g'

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: class cls(object):
            ....:    def new_meth(self): return 42
            ....:    old_meth = deprecated_function_alias(13109, new_meth)
            sage: cls.old_meth.__name__
            'old_meth'
            sage: cls().old_meth.__name__
            'old_meth'

            sage: cython('\n'.join([
            ....:     r"from sage.misc.superseded import deprecated_function_alias",
            ....:     r"cdef class cython_cls(object):",
            ....:     r"    def new_cython_meth(self):",
            ....:     r"        return 1",
            ....:     r"    old_cython_meth = deprecated_function_alias(13109, new_cython_meth)"
            ....: ]))
            ....:
            sage: cython_cls().old_cython_meth.__name__
            'old_cython_meth'
        """
        # first look through variables in stack frames
        for frame in inspect.stack():
            for name, obj in iteritems(frame[0].f_globals):
                if obj is self:
                    return name
        # then search object that contains self as method
        import gc, copy
        gc.collect()
        def is_class(gc_ref):
            if not isinstance(gc_ref, dict):
                return False
            is_python_class = '__module__' in gc_ref or '__package__' in gc_ref
            is_cython_class = '__new__' in gc_ref
            return is_python_class or is_cython_class
        search_for = self if (self.unbound is None) else self.unbound
        for ref in gc.get_referrers(search_for):
            if is_class(ref) and ref is not self.__dict__:
                ref_copy = copy.copy(ref)
                for key, val in iteritems(ref_copy):
                    if val is search_for:
                        return key
        raise AttributeError("The name of this deprecated function can not be determined")

    def __call__(self, *args, **kwds):
        """
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: def bla(): return 42
            sage: blo = deprecated_function_alias(13109, bla)
            sage: blo()
            doctest:...: DeprecationWarning: blo is deprecated. Please use bla instead.
            See http://trac.sagemath.org/13109 for details.
            42
        """
        if self.instance is None and self.__module__ != self.func.__module__:
            other = self.func.__module__ + "." + self.func.__name__
        else:
            other = self.func.__name__

        deprecation(self.trac_number,
                    "%s is deprecated. Please use %s instead."%(self.__name__, other))
        if self.instance is None:
            return self.func(*args, **kwds)
        else:
            return self.func(self.instance, *args, **kwds)

    def __get__(self, inst, cls = None):
        """
        TESTS::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: class cls(object):
            ....:    def new_meth(self): return 42
            ....:    old_meth = deprecated_function_alias(13109, new_meth)
            sage: obj = cls()
            sage: obj.old_meth.instance is obj
            True

        :trac:`19125`::

            sage: from sage.misc.superseded import deprecated_function_alias
            sage: class A:
            ....:    def __init__(self, x):
            ....:        self.x = x
            ....:    def f(self, y):
            ....:        return self.x+y
            ....:    g = deprecated_function_alias(42, f)
            sage: a1 = A(1)
            sage: a2 = A(2)
            sage: a1.g(a2.g(0))
            doctest:...: DeprecationWarning: g is deprecated. Please use f instead.
            See http://trac.sagemath.org/42 for details.
            3
            sage: a1.f(a2.f(0))
            3

        """
        return self if (inst is None) else DeprecatedFunctionAlias(self.trac_number, self.func, self.__module__, instance = inst, unbound = self)


def deprecated_function_alias(trac_number, func):
    """
    Create an aliased version of a function or a method which raise a
    deprecation warning message.

    If f is a function or a method, write
    ``g = deprecated_function_alias(trac_number, f)``
    to make a deprecated aliased version of f.

    INPUT:

    - ``trac_number`` -- integer. The trac ticket number where the
      deprecation is introduced.

    - ``func`` -- the function or method to be aliased

    EXAMPLES::

        sage: from sage.misc.superseded import deprecated_function_alias
        sage: g = deprecated_function_alias(13109, number_of_partitions)
        sage: g(5)
        doctest:...: DeprecationWarning: g is deprecated. Please use sage.combinat.partition.number_of_partitions instead.
        See http://trac.sagemath.org/13109 for details.
        7

    This also works for methods::

        sage: class cls(object):
        ....:    def new_meth(self): return 42
        ....:    old_meth = deprecated_function_alias(13109, new_meth)
        sage: cls().old_meth()
        doctest:...: DeprecationWarning: old_meth is deprecated. Please use new_meth instead.
        See http://trac.sagemath.org/13109 for details.
        42

    :trac:`11585`::

        sage: def a(): pass
        sage: b = deprecated_function_alias(13109, a)
        sage: b()
        doctest:...: DeprecationWarning: b is deprecated. Please use a instead.
        See http://trac.sagemath.org/13109 for details.

    AUTHORS:

     - Florent Hivert (2009-11-23), with the help of Mike Hansen.
     - Luca De Feo (2011-07-11), printing the full module path when different from old path
    """
    _check_trac_number(trac_number)
    module_name = inspect.getmodulename(
        inspect.currentframe(1).f_code.co_filename)
    if module_name is None:
        module_name = '__main__'
    return DeprecatedFunctionAlias(trac_number, func, module_name)


def deprecated_callable_import(trac_number, module_name, globs, locs, fromlist, message=None):
    """
    Imports a list of callables into the namespace from
    which it is called.  These callables however give a deprecation
    warning whenever they are called.  This is primarily used from
    deprecating things from Sage's ``all.py`` files.

    INPUT:

    - ``trac_number`` -- integer. The trac ticket number where the
      deprecation is introduced.

    - ``param module_name`` -- string or ``None``. The name of the
      module from which to import the callables or ``None``.

    - ``globs`` -- dictionary. The ``globals()`` from where this is being called.

    - ``locs`` -- dictionary. The ``locals()`` from where this is being called.

    - ``fromlist`` -- list of strings. The list the names of the
      callables to deprecate

    - ``message`` -- string. Message to display when the deprecated functions are called.

    .. note::

       If ``module_name`` is ``None``, then no importing will be done, and
       it will be assumed that the functions have already been
       imported and are present in ``globs``

    .. warning::

       This should really only be used for functions.

    EXAMPLES::

       sage: from sage.misc.superseded import deprecated_callable_import
       sage: is_prime(3)
       True
       sage: message = "Using %(name)s from here is deprecated."
       sage: deprecated_callable_import(13109, None, globals(), locals(), ['is_prime'], message)
       sage: is_prime(3)
       doctest:...: DeprecationWarning:
       Using is_prime from here is deprecated.
       See http://trac.sagemath.org/13109 for details.
       True
       sage: del is_prime
       sage: deprecated_callable_import(13109, 'sage.arith.all', globals(), locals(), ['is_prime'])
       sage: is_prime(3)
       doctest:...: DeprecationWarning:
       Using is_prime from here is deprecated.  If you need to use it, please import it directly from sage.arith.all.
       See http://trac.sagemath.org/13109 for details.
       True
    """
    _check_trac_number(trac_number)
    if message is None:
        message = '\nUsing %(name)s from here is deprecated. ' + \
            'If you need to use it, please import it directly from %(module_name)s.'
    from functools import partial
    from sage.misc.decorators import sage_wraps
    if module_name is None:
        mod_dict = globs
    else:
        mod_dict = __import__(module_name, globs, locs, fromlist).__dict__
    for name in fromlist:
        func = mod_dict[name]
        def wrapper(func, name, *args, **kwds):
            from sage.misc.superseded import deprecation
            deprecation(trac_number, message%{'name': name, 'module_name': module_name})
            return func(*args, **kwds)
        wrapped_function = sage_wraps(func)(partial(wrapper, func, name))
        wrapped_function.__doc__ = message%{'name': name, 'module_name': module_name}
        globs[name] = wrapped_function
    del name
