"""
Function, method, and module pickling

AUTHORS:

- Yi Qiang (2008-06-17): initial version

- William Stein (2008-06-17): added docstrings

- William Stein (2008-06-19): added call_pickled_function()

- Julian Rueth (2011-10-29): copied pickleMethod(), unpickleMethod(),
  pickleModule(), unpickleModule() from twisted.persisted.style

- Julian Rueth (2011-12-20): added docstrings
"""
#*****************************************************************************
#       Copyright (C) 2008 Yi Qiang <yqiang@gmail.com>
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Julian Rueth <julian.rueth@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import types, copy_reg, cPickle

def pickle_code(co):
    """
    This function provides pickling for code objects.

    INPUT::

    - ``co`` -- a code object to pickle

    OUTPUT::

    Returns data that can be used to unpickle the code object.

    TESTS::

        sage: def f(): return 1
        sage: pickled = pickle_function(f) # indirect doctest
        sage: unpickle_function(pickled)()
        1
        sage: def g(x):
        ...       def h(): return x
        ...       return h
        sage: pickled = pickle_function(g) # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: Cannot pickle code objects from closures
    """
    if co.co_freevars or co.co_cellvars:
        raise ValueError, "Cannot pickle code objects from closures"
    return code_ctor, (co.co_argcount, co.co_nlocals, co.co_stacksize,
                       co.co_flags, co.co_code, co.co_consts, co.co_names,
                       co.co_varnames, co.co_filename, co.co_name,
                       co.co_firstlineno, co.co_lnotab)

def code_ctor(*args):
    """
    This function provides unpickling for code objects pickled with :meth:`pickle_code`.

    TESTS::

        sage: def f(): return 1
        sage: pickled = pickle_function(f) # indirect doctest
        sage: unpickle_function(pickled)()
        1
    """
    return types.CodeType(*args)

copy_reg.pickle(types.CodeType, pickle_code)

def pickle_function(func):
    """
    This function provides pickling for functions. This is not a normal pickle; you
    must use the unpickle_function method to unpickle the pickled function.

    .. note::

        This does not work on all functions, but does work on
        'surprisingly' many functions.  In particular, it does not
        work on functions that includes nested functions.

    INPUT:

    - ``func`` -- a Python function

    OUTPUT:

    Returns the pickled function as a string.

    EXAMPLES::

        sage: def f(n): return n+1
        sage: pickled = pickle_function(f) # indirect doctest
        sage: unpickle_function(pickled)(2)
        3
    """
    return cPickle.dumps(func.func_code)

def unpickle_function(pickled):
    """
    Unpickle a pickled function that was pickled with :meth:`pickle_function`.

    INPUT:

    - ``pickled`` -- the pickled function as a string

    OUTPUT:

    Returns the unpickled function.

    EXAMPLES::

        sage: def f(n,m): return n*m
        sage: pickled = pickle_function(f) # indirect doctest
        sage: unpickle_function(pickled)(2,3)
        6
    """
    recovered = cPickle.loads(pickled)
    return types.FunctionType(recovered, globals())

def pickleMethod(method):
    """
    This function provides pickling for methods. It was originally copied from
    twisted.persisted.styles.  We chose not to import twisted.persisted.styles
    anymore since it takes a long time to import and we hardly use anything it
    provides. (see #11716)

    INPUT:

    - ``method`` -- a method to pickle

    OUTPUT:

    Returns data that can be used to unpickle the method.

    TESTS::

        sage: class A:
        ...       def __init__(self, x):
        ...           self.x = x
        ...       def f(self): return self.x
        ...       @classmethod
        ...       def g(cls): return 0
        ...
        sage: a = A(1)
        sage: from sage.misc.fpickle import pickleMethod
        sage: unpickleMethod, data = pickleMethod(a.f)
        sage: unpickleMethod(*data)()
        1
        sage: unpickleMethod, data = pickleMethod(A.g)
        sage: unpickleMethod(*data)()
        0
    """
    if isinstance(method.__self__,types.ClassType):
        return unpickleMethod, (method.__name__, None, method.__self__)
    else:
        return unpickleMethod, (method.__name__, method.__self__, method.__self__.__class__)

def unpickleMethod(im_name, im_self, im_class):
    """
    This function provides unpickling for methods pickled with help of
    :meth:`pickleMethod`. It was originally copied from twisted.persisted.styles.
    We chose not to import twisted.persisted.styles anymore since it takes a long
    time to import and we hardly use anything it provides. (see #11716)

    INPUT:

    - ``im_name`` -- the method's name

    - ``im_self`` -- the method's class instance object or None if this is a class method

    - `` im_class`` -- the class which defined the method

    OUTPUT:

    Returns the unpickled instance method.

    TESTS::

        sage: class A:
        ...       def __init__(self, x):
        ...           self.x = x
        ...       def f(self): return self.x
        ...       @classmethod
        ...       def g(cls): return 0
        ...
        sage: a = A(1)
        sage: from sage.misc.fpickle import pickleMethod
        sage: unpickleMethod, data = pickleMethod(a.f)
        sage: unpickleMethod(*data)()
        1
        sage: del A.f
        sage: unpickleMethod(*data)
        Traceback (most recent call last):
        ...
        UnpicklingError: no method `f' in class `A'
        sage: unpickleMethod, data = pickleMethod(A.g)
        sage: unpickleMethod(*data)()
        0
        sage: del A.g
        sage: unpickleMethod(*data)()
        Traceback (most recent call last):
        ...
        UnpicklingError: no method `g' in class `A'
    """
    if hasattr(im_class, im_name):
        unbound = getattr(im_class, im_name)
    elif im_self is not None and hasattr(im_self.__class__, im_name):
        # Attempt a common fix before bailing -- if classes have
        # changed around since we pickled this method, we may still be
        # able to get it by looking on the instance's current class.
        unbound = getattr(im_self.__class__, im_name)
    else:
        from pickle import UnpicklingError
        raise UnpicklingError, "no method `%s' in class `%s'"%(im_name, im_class.__name__)

    if im_self is None:
        return unbound
    return types.MethodType(unbound.im_func, im_self)

copy_reg.pickle(types.MethodType, pickleMethod, unpickleMethod)

def pickleModule(module):
    """
    This function provides pickling for modules. It was originally copied from
    twisted.persisted.styles.  We chose not to import twisted.persisted.styles
    anymore since it takes a long time to import and we hardly use anything it
    provides. (see #11716)

    INPUT:

    - ``module`` -- a module to pickle

    OUTPUT:

    Returns data that can be used to unpickle the module.

    TESTS::

        sage: from sage.misc.fpickle import pickleModule
        sage: unpickleModule, data = pickleModule(sage)
        sage: unpickleModule(*data) is sage
        True
    """
    return unpickleModule, (module.__name__,)

def unpickleModule(name):
    """
    This function provides unpickling for modules pickled with the help of
    :meth:`pickleModule`. It was originally copied from twisted.persisted.styles.
    We chose not to import twisted.persisted.styles anymore since it takes a long
    time to import and we hardly use anything it provides. (see #11716)

    INPUT:

    - ``name`` -- the module's name

    OUTPUT:

    Returns the unpickled module.

    TESTS::

        sage: from sage.misc.fpickle import pickleModule
        sage: unpickleModule, data = pickleModule(sage)
        sage: unpickleModule(*data) is sage
        True
    """
    return __import__(name,{},{})

copy_reg.pickle(types.ModuleType, pickleModule, unpickleModule)
