#*****************************************************************************
#       Copyright (C) 2014 Sage developers, http://trac.sagemath.org/wiki
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.symbolic.function import GinacFunction, BuiltinFunction
from sage.rings.all import Integer, Rational, RealField, RR, QQ, ZZ, ComplexField
from sage.rings.complex_number import is_ComplexNumber
from sage.symbolic.expression import Expression
from sage.structure.element import get_coercion_model
from sage.functions.other import gamma
from sage.misc.latex import latex
import math


class IntegerFunction(BuiltinFunction):
    
    def __init__(self, name, nargs=1, latex_name=None, conversions=None,
            evalf_params_first=True, alt_name=None):
        """
        TESTS::

        """
        BuiltinFunction.__init__(self, name, nargs, latex_name, conversions,
                evalf_params_first, alt_name = alt_name)


    def _is_numerical(self, x):
        """
        Return True if `x` is a numerical object.

        This is used to determine whether to call the :meth:`_evalf_`
        method instead of the :meth:`_eval_` method.

        This is a non-static method since whether or not an argument is
        considered numerical may depend on the specific function.

        TESTS::

            sage: factorial._is_numerical(5)
            True
            sage: factorial._is_numerical(5.)
            True
            sage: factorial._is_numerical(pi)
            False
            sage: factorial._is_numerical(5r)
            True
            sage: factorial._is_numerical(5.4r)
            True
            sage: factorial._is_numerical(QQ(1/2))
            True
        """
        if isinstance(x, (float, complex, int, Integer, Rational)):
            return True
        return hasattr(x, 'precision')


class IntegerGinacFunction(IntegerFunction, GinacFunction):
    
    def __init__(self, name, nargs=1, latex_name=None, conversions=None,
            evalf_params_first=True):
        """
        TESTS::

            sage: from sage.functions.trig import Function_cot
            sage: c = Function_cot() # indirect doctest
            sage: c(pi/2)
            0
        """
        IntegerFunction.__init__(self, name, nargs, latex_name, conversions,
                evalf_params_first)
        GinacFunction.__init__(self, name, nargs, latex_name=latex_name,
                               conversions=conversions,
                               evalf_params_first=evalf_params_first)


class Function_factorial(IntegerGinacFunction):
    def __init__(self):
        r"""
        Return the factorial of `n`.

        INPUT:

        -  ``n`` - any complex argument (except negative
           integers) or any symbolic expression

        OUTPUT: an integer or symbolic expression

        EXAMPLES::

            sage: x = var('x')
            sage: factorial(0)
            1
            sage: factorial(4)
            24
            sage: factorial(10)
            3628800
            sage: factorial(6) == 6*5*4*3*2
            True
            sage: f = factorial(x + factorial(x)); f
            factorial(x + factorial(x))
            sage: f(x=3)
            362880
            sage: factorial(x)^2
            factorial(x)^2

        To prevent automatic evaluation use the ``hold`` argument::

            sage: factorial(5,hold=True)
            factorial(5)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: factorial(5,hold=True).simplify()
            120

        We can also give input other than nonnegative integers.  For
        other nonnegative numbers, the :func:`gamma` function is used::

            sage: factorial(1/2)
            1/2*sqrt(pi)
            sage: factorial(3/4)
            gamma(7/4)
            sage: factorial(2.3)
            2.68343738195577

        TESTS:

            sage: factorial(6) == 6*5*4*3*2
            True
            sage: factorial(1) == factorial(0)
            True
            sage: factorial(71) == 71* factorial(70)
            True

        We verify that we can convert this function to Maxima and
        bring it back into Sage.::

            sage: z = var('z')
            sage: factorial._maxima_init_()
            'factorial'
            sage: maxima(factorial(z))
            factorial(_SAGE_VAR_z)
            sage: _.sage()
            factorial(z)
            sage: k = var('k')
            sage: factorial(k)
            factorial(k)

            sage: factorial(3.14)
            7.173269190187...

        Test latex typesetting::

            sage: latex(factorial(x))
            x!
            sage: latex(factorial(2*x))
            \left(2 \, x\right)!
            sage: latex(factorial(sin(x)))
            \sin\left(x\right)!
            sage: latex(factorial(sqrt(x+1)))
            \left(\sqrt{x + 1}\right)!
            sage: latex(factorial(sqrt(x)))
            \sqrt{x}!
            sage: latex(factorial(x^(2/3)))
            \left(x^{\frac{2}{3}}\right)!

            sage: latex(factorial)
            {\rm factorial}

        Check that #11539 is fixed::

            sage: (factorial(x) == 0).simplify()
            factorial(x) == 0
            sage: maxima(factorial(x) == 0).sage()
            factorial(x) == 0
            sage: y = var('y')
            sage: (factorial(x) == y).solve(x)
            [factorial(x) == y]

        Test pickling::

            sage: loads(dumps(factorial))
            factorial
        """
        IntegerGinacFunction.__init__(self, "factorial", nargs=1,
                                 latex_name='{\\rm factorial}',
                conversions=dict(maxima='factorial',
                                 mathematica='Factorial',
                                 sympy='factorial'))

    def _evalf_(self, x, **kwds):
        """
        Return the factorial function.
 
        Note that this method overrides the eval method defined in GiNaC
        which calls numeric evaluation on all numeric input. We preserve
        exact results if the input is a rational number.
 
        EXAMPLES::
 
            sage: k = var('k')
            sage: k.factorial()
            factorial(k)
            sage: SR(1/2).factorial()
            1/2*sqrt(pi)
            sage: SR(3/4).factorial()
            gamma(7/4)
            sage: SR(5).factorial()
            120
            sage: SR(3245908723049857203948572398475r).factorial()
            factorial(3245908723049857203948572398475L)
            sage: SR(3245908723049857203948572398475).factorial()
            factorial(3245908723049857203948572398475)
            sage: factorial(5,hold=True).n(200)
            120.000000000000000000000000000...
            sage: parent(factorial(QQ(4)))
            Rational Field
        """
        algorithm = kwds.get('algorithm','gmp')
        hold = kwds.get('hold', False)
#        print('algo:, hold:'), algorithm, kwds.get('hold', False)
        if hold is not True:
            parent = ZZ
            if isinstance(x, Rational) and x.is_integer():
                parent = QQ
                x = ZZ(x)
            if isinstance(x, (Integer, int)):
                if algorithm == 'gmp' or algorithm is None:
#                    print('calling ZZ(x).factorial()')
                    return parent(ZZ(x).factorial())
                elif algorithm == 'pari':
                    from sage.libs.pari.pari_instance import pari
#                    print('calling pari')
                    return parent(pari.factorial(x))
                else:
                    raise ValueError('unknown algorithm')
            else:
                # let gamma deal with other parents
#                print('calling gamma()')
                return gamma(x+1)

factorial = Function_factorial()
