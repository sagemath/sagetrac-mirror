r"""
Arithmetic functions

AUTHORS:

- Jonas Jermann (2015): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013-2015 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import Element
from sage.rings.all import ZZ
from sage.rings.arith import divisors

class ArithmeticFunctionElement(Element):
    r"""
    Element of ArithmeticFunctions
    """

    def __init__(self, parent, f, repr_function=None, repr_latex=None, check=True, convert=True):
        r"""
        Element of a ArithmeticFunctions
        """

        if check or convert:
            from sage.structure.element import parent as getParent
            for _ in range(3):
                res = f(ZZ(ZZ.random_element(1,10)))
                if convert:
                    parent = parent.extend_codomain(getParent(res))
                elif check:
                    if not (res in parent.codomain):
                        raise ValueError("f is not a function from N to {}".format(parent.codomain))

        self.function = f

        if not (repr_function is None):
            self.repr_function = repr_function
        else:
            self.repr_function = lambda var, get_var: "<fun>({})".format(var)

        if not (repr_latex is None):
            self.repr_latex = repr_latex
        elif repr_function is None:
            self.repr_latex = lambda var, get_var: "\\left<fun\\right>\\left({}\\right)".format(var)
        else:
            self.repr_latex = self.repr_function

        super(ArithmeticFunctionElement, self).__init__(parent)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.
        """

        result = self.repr_function(self.parent().get_variable(), self.parent().get_variable)
        self.parent().reset_variables()

        return result

    def _latex_(self):
        r"""
        Return the string representation of ``self``.
        """

        result = self.repr_latex(self.parent().get_variable(), self.parent().get_variable)
        self.parent().reset_variables()

        return result

    def _eq_(self, other):
        return self.function == other.function
    def _neq_(self, other):
        return not (self == other)

    def _add_(self,other):
        r"""
        Return the sum of ``self`` and ``other``.
        """

        def f(n):
            return self.function(n)+other.function(n)
        def repr_f(var, get_var):
            return "({}+{})".format(self.repr_function(var, get_var), other.repr_function(var, get_var))
        def repr_latex(var, get_var):
            return "\\left({}+{}\\right)".format(self.repr_latex(var, get_var), other.repr_latex(var, get_var))

        return self.__class__(self.parent(), f, repr_f, repr_latex)

    def _sub_(self,other):
        r"""
        Return the difference of ``self`` and ``other``.
        """

        def f(n):
            return self.function(n)-other.function(n)
        def repr_f(var, get_var):
            return "({}-{})".format(self.repr_function(var, get_var), other.repr_function(var, get_var))
        def repr_latex(var, get_var):
            return "\\left({}-{}\\right)".format(self.repr_latex(var, get_var), other.repr_latex(var, get_var))

        return self.__class__(self.parent(), f, repr_f, repr_latex)

    def _mul_(self,other):
        r"""
        Return the convolution product of ``self`` and ``other``.
        """

        def f(n):
            return sum([self.function(d)*other.function(n/d) for d in divisors(n)])
        def repr_f(var, get_var):
            var2=get_var()
            return "sum([{}*{} for {} in divisors({})])".format(self.repr_function(var2, get_var), other.repr_function("{}/{}".format(var,var2), get_var), var2, var)
        def repr_latex(var, get_var):
            var2=get_var()
            return "\\sum_{{\\left.{}\\middle|{}\\right.}}\\left({}\\cdot{}\\right)".format(var2, var, self.repr_latex(var2, get_var), other.repr_latex("\\frac{{{}}}{{{}}}".format(var,var2), get_var))

        return self.__class__(self.parent(), f, repr_f, repr_latex)

    def pmul(self,other):
        r"""
        Return the pointwise product of ``self`` and ``other``.
        """

        def f(n):
            return self.function(n)*other.function(n)
        def repr_f(var, get_var):
            return "({}*{})".format(self.repr_function(var, get_var),other.repr_function(var, get_var))
        def repr_latex(var, get_var):
            return "\\left({}\\cdot{}\\right)".format(self.repr_latex(var, get_var),other.repr_latex(var, get_var))

        return self.__class__(self.parent(), f, repr_f, repr_latex)

    def pdiv(self,other):
        r"""
        Return the pointwise division of ``self`` by ``other``.
        """

        def f(n):
            return self.function(n)/other.function(n)
        def repr_f(var, get_var):
            return "({}/{})".format(self.repr_function(var, get_var),other.repr_function(var, get_var))
        def repr_latex(var, get_var):
            return "\\left(\\frac{{{}}}{{{}}}\\right)".format(self.repr_latex(var, get_var),other.repr_latex(var, get_var))

        return self.__class__(self.parent(), f, repr_f, repr_latex)

    def __call__(self, n):
        r"""
        Evaluate the function at the given argument.
        """

        return self.function(n)
