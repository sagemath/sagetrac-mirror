r"""
Symbolic Minimum and Maximum

Sage provides a symbolic maximum and minimum due to the fact that the Python
builtin max and min are not able to deal with variables as users might expect.
These functions wait to evaluate if there are variables.

Here you can see some differences::

   sage: max(x,x^2)
   x
   sage: max_symbolic(x,x^2)
   max(x, x^2)
   sage: f(x) = max_symbolic(x,x^2); f(1/2)
   1/2

This works as expected for more than two entries::

   sage: max(3,5,x)
   5
   sage: min(3,5,x)
   3
   sage: max_symbolic(3,5,x)
   max(x, 5)
   sage: min_symbolic(3,5,x)
   min(x, 3)

"""
###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2010 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################
from __future__ import absolute_import

from sage.symbolic.function import BuiltinFunction
from sage.symbolic.expression import Expression
from sage.symbolic.ring import SR

from six.moves.builtins import max as builtin_max, min as builtin_min

class MinMax_base(BuiltinFunction):
    def eval_helper(self, this_f, builtin_f, initial_val, args):
        """
        EXAMPLES::

            sage: max_symbolic(3,5,x) # indirect doctest
            max(x, 5)
            sage: min_symbolic(3,5,x)
            min(x, 3)
        """
        # __call__ ensures that if args is a singleton, the element is iterable
        arg_is_iter = False
        if len(args) == 1:
            arg_is_iter = True
            args = args[0]

        symb_args = []
        res = initial_val
        num_non_symbolic_args = 0
        for x in args:
            if isinstance(x, Expression):
                symb_args.append(x)
            else:
                num_non_symbolic_args += 1
                res = builtin_f(res, x)

        # if no symbolic arguments, return the result
        if len(symb_args) == 0:
            if res is None:
                # this is a hack to get the function to return None to the user
                # the convention to leave a given symbolic function unevaluated
                # is to return None from the _eval_ function, so we need
                # a trick to indicate that the return value of the function is
                # really None
                # this is caught in the __call__ method, which knows to return
                # None in this case
                raise ValueError("return None")
            return res

        # if all arguments were symbolic return
        if num_non_symbolic_args <= 1 and not arg_is_iter:
            return None

        if res is not None: symb_args.append(res)
        return this_f(*symb_args)

    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: max_symbolic(3,5,x)
            max(x, 5)
            sage: max_symbolic(3,5,x, hold=True)
            max(3, 5, x)
            sage: max_symbolic([3,5,x])
            max(x, 5)

        ::

            sage: min_symbolic(3,5,x)
            min(x, 3)
            sage: min_symbolic(3,5,x, hold=True)
            min(3, 5, x)
            sage: min_symbolic([3,5,x])
            min(x, 3)

        TESTS:

        We get an exception if no arguments are given::

            sage: max_symbolic()
            Traceback (most recent call last):
            ...
            ValueError: number of arguments must be > 0

        Check if we return None, when the builtin function would::

            sage: max_symbolic([None]) is None
            True
            sage: max_symbolic([None, None]) is None
            True
            sage: min_symbolic([None]) is None
            True
            sage: min_symbolic([None, None]) is None
            True

        Check if a single argument which is not iterable works::

            sage: max_symbolic(None)
            Traceback (most recent call last):
            ...
            TypeError: 'NoneType' object is not iterable
            sage: max_symbolic(5)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
            sage: max_symbolic(x)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.symbolic.expression.Expression' object is not iterable
            sage: min_symbolic(5)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
            sage: min_symbolic(x)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.symbolic.expression.Expression' object is not iterable
        """
        if len(args) == 0:
            raise ValueError("number of arguments must be > 0")
        if len(args) == 1:
            try:
                args=(SR._force_pyobject(iter(args[0])),)
            except TypeError as e:
                raise e

        try:
            return BuiltinFunction.__call__(self, *args, **kwds)
        except ValueError as e:
            if e.args[0] == "return None":
                return None

class MaxSymbolic(MinMax_base):
    def __init__(self):
        r"""
        Symbolic `\max` function.

        The Python builtin `\max` function doesn't work as expected when symbolic
        expressions are given as arguments. This function delays evaluation
        until all symbolic arguments are substituted with values.

        EXAMPLES::

            sage: max_symbolic(3, x)
            max(3, x)
            sage: max_symbolic(3, x).subs(x=5)
            5
            sage: max_symbolic(3, 5, x)
            max(x, 5)
            sage: max_symbolic([3,5,x])
            max(x, 5)

        TESTS::

            sage: loads(dumps(max_symbolic(x,5)))
            max(x, 5)
            sage: latex(max_symbolic(x,5))
            \max\left(x, 5\right)
            sage: max_symbolic(x, 5)._sympy_()
            Max(5, x)
        """
        BuiltinFunction.__init__(self, 'max', nargs=0, latex_name="\max",
                                 conversions=dict(sympy='Max'))

    def _eval_(self, *args):
        """
        EXAMPLES::

            sage: t = max_symbolic(x, 5); t
            max(x, 5)
            sage: t.subs(x=3) # indirect doctest
            5
            sage: max_symbolic(5,3)
            5
            sage: u = max_symbolic(*(range(10)+[x])); u
            max(x, 9)
            sage: u.subs(x=-1)
            9
            sage: u.subs(x=10)
            10
            sage: max_symbolic([0,x])
            max(x, 0)

        TESTS::

            sage: max_symbolic()
            Traceback (most recent call last):
            ...
            ValueError: number of arguments must be > 0
        """
        return self.eval_helper(max_symbolic, builtin_max, None, args)

    def _evalf_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: t = max_symbolic(sin(x), cos(x))
            sage: t.subs(x=1).n(200)
            0.84147098480789650665250232163029899962256306079837106567275
            sage: var('y')
            y
            sage: t = max_symbolic(sin(x), cos(x), y)
            sage: u = t.subs(x=1); u
            max(sin(1), cos(1), y)
            sage: u.n()
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate symbolic expression numerically

        ::

            sage: f = max_symbolic(sin(x), cos(x))
            sage: r = integral(f, x, 0, 1)
            sage: r.n()
            0.8739124411567263
        """
        return max_symbolic(args)

max_symbolic = MaxSymbolic()


class MinSymbolic(MinMax_base):
    def __init__(self):
        r"""
        Symbolic `\min` function.

        The Python builtin `\min` function doesn't work as expected when symbolic
        expressions are given as arguments. This function delays evaluation
        until all symbolic arguments are substituted with values.

        EXAMPLES::

            sage: min_symbolic(3, x)
            min(3, x)
            sage: min_symbolic(3, x).subs(x=5)
            3
            sage: min_symbolic(3, 5, x)
            min(x, 3)
            sage: min_symbolic([3,5,x])
            min(x, 3)

        TESTS::

            sage: loads(dumps(min_symbolic(x,5)))
            min(x, 5)
            sage: latex(min_symbolic(x,5))
            \min\left(x, 5\right)
            sage: min_symbolic(x, 5)._sympy_()
            Min(5, x)
        """
        BuiltinFunction.__init__(self, 'min', nargs=0, latex_name="\min",
                                 conversions=dict(sympy='Min'))

    def _eval_(self, *args):
        """
        EXAMPLES::

            sage: t = min_symbolic(x, 5); t
            min(x, 5)
            sage: t.subs(x=3) # indirect doctest
            3
            sage: min_symbolic(5,3)
            3
            sage: u = min_symbolic(*(range(10)+[x])); u
            min(x, 0)
            sage: u.subs(x=-1)
            -1
            sage: u.subs(x=10)
            0
            sage: min_symbolic([3,x])
            min(x, 3)

        TESTS::

            sage: min_symbolic()
            Traceback (most recent call last):
            ...
            ValueError: number of arguments must be > 0
        """
        return self.eval_helper(min_symbolic, builtin_min, float('inf'), args)

    def _evalf_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: t = min_symbolic(sin(x), cos(x))
            sage: t.subs(x=1).n(200)
            0.54030230586813971740093660744297660373231042061792222767010
            sage: var('y')
            y
            sage: t = min_symbolic(sin(x), cos(x), y)
            sage: u = t.subs(x=1); u
            min(sin(1), cos(1), y)
            sage: u.n()
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate symbolic expression numerically
        """
        return min_symbolic(args)

min_symbolic = MinSymbolic()
