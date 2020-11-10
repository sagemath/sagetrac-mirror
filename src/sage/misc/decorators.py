"""
Decorators

Python decorators for use in Sage.

AUTHORS:

- Tim Dumol (5 Dec 2009) -- initial version.
- Johan S. R. Nielsen (2010) -- collect decorators from various modules.
- Johan S. R. Nielsen (8 apr 2011) -- improve introspection on decorators.
- Simon King (2011-05-26) -- put this file into the reference manual.

"""
#*****************************************************************************
#       Copyright (C) 2009 Tim Dumol
#                     2010,2011 Johan S. R. Nielsen
#                     2011 Simon King <simon.king@uni-jena.de>
#                     2014 Julian Rueth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from functools import (partial, update_wrapper, WRAPPER_ASSIGNMENTS,
                       WRAPPER_UPDATES)
from copy import copy

from decorator import decorate, decorator

from sage.misc.sageinspect import (sage_getsource, sage_getsourcelines,
                                   sage_getargspec)
from inspect import ArgSpec

# Infix operator decorator
class infix_operator(object):
    """
    A decorator for functions which allows for a hack that makes
    the function behave like an infix operator.

    This decorator exists as a convenience for interactive use.

    EXAMPLES:

    An infix dot product operator::

        sage: @infix_operator('multiply')
        ....: def dot(a, b):
        ....:     '''Dot product.'''
        ....:     return a.dot_product(b)
        sage: u = vector([1, 2, 3])
        sage: v = vector([5, 4, 3])
        sage: u *dot* v
        22

    An infix element-wise addition operator::

        sage: @infix_operator('add')
        ....: def eadd(a, b):
        ....:   return a.parent([i + j for i, j in zip(a, b)])
        sage: u = vector([1, 2, 3])
        sage: v = vector([5, 4, 3])
        sage: u +eadd+ v
        (6, 6, 6)
        sage: 2*u +eadd+ v
        (7, 8, 9)

    A hack to simulate a postfix operator::

        sage: @infix_operator('or')
        ....: def thendo(a, b):
        ....:     return b(a)
        sage: x |thendo| cos |thendo| (lambda x: x^2)
        cos(x)^2
    """

    operators = {
        'add': {'left': '__add__', 'right': '__radd__'},
        'multiply': {'left': '__mul__', 'right': '__rmul__'},
        'or': {'left': '__or__', 'right': '__ror__'},
    }

    def __init__(self, precedence):
        """
        INPUT:

        - ``precedence`` -- one of ``'add'``, ``'multiply'``, or ``'or'``
          indicating the new operator's precedence in the order of operations.
        """
        self.precedence = precedence

    def __call__(self, func):
        """Returns a function which acts as an inline operator."""

        left_meth = self.operators[self.precedence]['left']
        right_meth = self.operators[self.precedence]['right']
        wrapper_name = func.__name__
        wrapper_members = {
            'function': staticmethod(func),
            left_meth: _infix_wrapper._left,
            right_meth: _infix_wrapper._right,
            '_sage_src_': lambda: sage_getsource(func)
        }
        for attr in WRAPPER_ASSIGNMENTS:
            try:
                wrapper_members[attr] = getattr(func, attr)
            except AttributeError:
                pass

        wrapper = type(wrapper_name, (_infix_wrapper,), wrapper_members)

        wrapper_inst = wrapper()
        wrapper_inst.__dict__.update(getattr(func, '__dict__', {}))
        return wrapper_inst


class _infix_wrapper(object):
    function = None

    def __init__(self, left=None, right=None):
        """
        Initialize the actual infix object, with possibly a specified left
        and/or right operand.
        """
        self.left = left
        self.right = right

    def __call__(self, *args, **kwds):
        """Call the passed function."""
        return self.function(*args, **kwds)

    def _left(self, right):
        """The function for the operation on the left (e.g., __add__)."""
        if self.left is None:
            if self.right is None:
                new = copy(self)
                new.right = right
                return new
            else:
                raise SyntaxError("Infix operator already has its "
                                  "right argument")
        else:
            return self.function(self.left, right)

    def _right(self, left):
        """The function for the operation on the right (e.g., __radd__)."""
        if self.right is None:
            if self.left is None:
                new = copy(self)
                new.left = left
                return new
            else:
                raise SyntaxError("Infix operator already has its "
                                  "left argument")
        else:
            return self.function(left, self.right)


class suboptions(object):
    def __init__(self, name, **options):
        """
        A decorator for functions which collects all keywords
        starting with ``name+'_'`` and collects them into a dictionary
        which will be passed on to the wrapped function as a
        dictionary called ``name_options``.

        The keyword arguments passed into the constructor are taken
        to be default for the ``name_options`` dictionary.

        EXAMPLES::

            sage: from sage.misc.decorators import suboptions
            sage: s = suboptions('arrow', size=2)
            sage: s.name
            'arrow_'
            sage: s.options
            {'size': 2}
        """
        self.name = name + "_"
        self.options = options

    def __call__(self, func):
        """
        Returns a wrapper around func

        EXAMPLES::

            sage: from sage.misc.decorators import suboptions
            sage: def f(*args, **kwds): print(sorted(kwds.items()))
            sage: f = suboptions('arrow', size=2)(f)
            sage: f(size=2)
            [('arrow_options', {'size': 2}), ('size', 2)]
            sage: f(arrow_size=3)
            [('arrow_options', {'size': 3})]
            sage: f(arrow_options={'size':4})
            [('arrow_options', {'size': 4})]
            sage: f(arrow_options={'size':4}, arrow_size=5)
            [('arrow_options', {'size': 5})]

         Demonstrate that the introspected argument specification of the
         wrapped function is updated (see :trac:`9976`).

            sage: from sage.misc.sageinspect import sage_getargspec
            sage: sage_getargspec(f)
            ArgSpec(args=['arrow_size'], varargs='args', keywords='kwds', defaults=(2,))
        """
        def wrapper(func, *args, **kwds):
            suboptions = copy(self.options)
            suboptions.update(kwds.pop(self.name+"options", {}))

            # Collect all the relevant keywords in kwds
            # and put them in suboptions
            for key, value in list(kwds.items()):
                if key.startswith(self.name):
                    suboptions[key[len(self.name):]] = value
                    del kwds[key]

            kwds[self.name + "options"] = suboptions

            return func(*args, **kwds)
        return decorate(func, wrapper)

        # Add the options specified by @options to the signature of the wrapped
        # function in the Sphinx-generated documentation (Trac 9976), using the
        # special attribute _sage_argspec_ (see e.g. sage.misc.sageinspect)
        def argspec():
            argspec = sage_getargspec(func)

            def listForNone(l):
                return l if not l is None else []
            newArgs = [self.name + opt for opt in self.options.keys()]
            args = (argspec.args if not argspec.args is None else []) + newArgs
            defaults = (argspec.defaults if not argspec.defaults is None else ()) \
                        + tuple(self.options.values())
            # Note: argspec.defaults is not always a tuple for some reason
            return ArgSpec(args, argspec.varargs, argspec.keywords, defaults)
        wrapper._sage_argspec_ = argspec

        return wrapper


class options(object):
    def __init__(self, **options):
        """
        A decorator for functions which allows for default options to be
        set and reset by the end user.  Additionally, if one needs to, one
        can get at the original keyword arguments passed into the
        decorator.

        TESTS::

            sage: from sage.misc.decorators import options
            sage: o = options(rgbcolor=(0,0,1))
            sage: o.options
            {'rgbcolor': (0, 0, 1)}
            sage: o = options(rgbcolor=(0,0,1), __original_opts=True)
            sage: o.original_opts
            True
            sage: loads(dumps(o)).options
            {'rgbcolor': (0, 0, 1)}

        Demonstrate that the introspected argument specification of the wrapped
        function is updated (see :trac:`9976`)::

            sage: from sage.misc.decorators import options
            sage: o = options(rgbcolor=(0,0,1))
            sage: def f(*args, **kwds):
            ....:     print("{} {}".format(args, sorted(kwds.items())))
            sage: f1 = o(f)
            sage: from sage.misc.sageinspect import sage_getargspec
            sage: sage_getargspec(f1)
            ArgSpec(args=['rgbcolor'], varargs='args', keywords='kwds', defaults=((0, 0, 1),))
        """
        self.options = options
        self.original_opts = options.pop('__original_opts', False)

    def __call__(self, func):
        """
        EXAMPLES::

            sage: from sage.misc.decorators import options
            sage: o = options(rgbcolor=(0,0,1))
            sage: def f(*args, **kwds):
            ....:     print("{} {}".format(args, sorted(kwds.items())))
            sage: f1 = o(f)
            sage: f1()
            () [('rgbcolor', (0, 0, 1))]
            sage: f1(rgbcolor=1)
            () [('rgbcolor', 1)]
            sage: o = options(rgbcolor=(0,0,1), __original_opts=True)
            sage: f2 = o(f)
            sage: f2(alpha=1)
            () [('__original_opts', {'alpha': 1}), ('alpha', 1), ('rgbcolor', (0, 0, 1))]

        """
        def wrapper(func, *args, **kwds):
            options = copy(wrapper.options)
            if self.original_opts:
                options['__original_opts'] = kwds
            options.update(kwds)
            return func(*args, **options)
        return decorate(func, wrapper)

        #Add the options specified by @options to the signature of the wrapped
        #function in the Sphinx-generated documentation (Trac 9976), using the
        #special attribute _sage_argspec_ (see e.g. sage.misc.sageinspect)
        def argspec():
            argspec = sage_getargspec(func)
            args = ((argspec.args if not argspec.args is None else []) +
                    list(self.options))
            defaults = (argspec.defaults or ()) + tuple(self.options.values())
            # Note: argspec.defaults is not always a tuple for some reason
            return ArgSpec(args, argspec.varargs, argspec.keywords, defaults)

        wrapper._sage_argspec_ = argspec

        def defaults():
            """
            Return the default options.

            EXAMPLES::

                sage: from sage.misc.decorators import options
                sage: o = options(rgbcolor=(0,0,1))
                sage: def f(*args, **kwds):
                ....:     print("{} {}".format(args, sorted(kwds.items())))
                sage: f = o(f)
                sage: f.options['rgbcolor']=(1,1,1)
                sage: f.defaults()
                {'rgbcolor': (0, 0, 1)}
            """
            return copy(self.options)

        def reset():
            """
            Reset the options to the defaults.

            EXAMPLES::

                sage: from sage.misc.decorators import options
                sage: o = options(rgbcolor=(0,0,1))
                sage: def f(*args, **kwds):
                ....:     print("{} {}".format(args, sorted(kwds.items())))
                sage: f = o(f)
                sage: f.options
                {'rgbcolor': (0, 0, 1)}
                sage: f.options['rgbcolor']=(1,1,1)
                sage: f.options
                {'rgbcolor': (1, 1, 1)}
                sage: f()
                () [('rgbcolor', (1, 1, 1))]
                sage: f.reset()
                sage: f.options
                {'rgbcolor': (0, 0, 1)}
                sage: f()
                () [('rgbcolor', (0, 0, 1))]
            """
            wrapper.options = copy(self.options)

        wrapper.options = copy(self.options)
        wrapper.reset = reset
        wrapper.reset.__doc__ = """
        Reset the options to the defaults.

        Defaults:
        %s
        """ % self.options

        wrapper.defaults = defaults
        wrapper.defaults.__doc__ = """
        Return the default options.

        Defaults:
        %s
        """ % self.options

        return wrapper


def rename_keyword(deprecated=None, deprecation=None, **renames):
    """
    A decorator which renames keyword arguments and optionally
    deprecates the new keyword.

    INPUT:

    - ``deprecation`` -- integer. The trac ticket number where the
        deprecation was introduced.

    - the rest of the arguments is a list of keyword arguments in the
        form ``renamed_option='existing_option'``.  This will have the
        effect of renaming ``renamed_option`` so that the function only
        sees ``existing_option``.  If both ``renamed_option`` and
        ``existing_option`` are passed to the function, ``existing_option``
        will override the ``renamed_option`` value.

    EXAMPLES::

        sage: from sage.misc.decorators import rename_keyword
        sage: r = rename_keyword(color='rgbcolor')
        sage: r.renames
        {'color': 'rgbcolor'}
        sage: loads(dumps(r)).renames
        {'color': 'rgbcolor'}

    To deprecate an old keyword::

        sage: r = rename_keyword(deprecation=13109, color='rgbcolor')

        sage: from sage.misc.decorators import rename_keyword
            sage: r = rename_keyword(color='rgbcolor')
            sage: def f(*args, **kwds):
            ....:     print("{} {}".format(args, kwds))
            sage: f = r(f)
            sage: f()
            () {}
            sage: f(alpha=1)
            () {'alpha': 1}
            sage: f(rgbcolor=1)
            () {'rgbcolor': 1}
            sage: f(color=1)
            () {'rgbcolor': 1}

        We can also deprecate the renamed keyword::

            sage: r = rename_keyword(deprecation=13109, deprecated_option='new_option')
            sage: def f(*args, **kwds):
            ....:     print("{} {}".format(args, kwds))
            sage: f = r(f)
            sage: f()
            () {}
            sage: f(alpha=1)
            () {'alpha': 1}
            sage: f(new_option=1)
            () {'new_option': 1}
            sage: f(deprecated_option=1)
            doctest:...: DeprecationWarning: use the option 'new_option' instead of 'deprecated_option'
            See http://trac.sagemath.org/13109 for details.
            () {'new_option': 1}
    """
    assert deprecated is None, 'Use @rename_keyword(deprecation=<trac_number>, ...)'
    
    def wrapper(func, *args, **kwds):
        for old_name, new_name in renames.items():
            if old_name in kwds and new_name not in kwds:
                if deprecation is not None:
                    from sage.misc.superseded import deprecation as print_deprecation
                    print_deprecation(deprecation, "use the option "
                                "%r instead of %r" % (new_name, old_name))
                kwds[new_name] = kwds[old_name]
                del kwds[old_name]
        return func(*args, **kwds)

    return decorator(wrapper)
