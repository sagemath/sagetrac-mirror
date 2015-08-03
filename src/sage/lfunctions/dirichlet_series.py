"""
Class for manipulation of formal Dirichlet series
=================================================

        EXAMPLES::

            sage: s,t = var('s,t')
            sage: dirichlet_series(1)
            1 + O(20^(-s))
            sage: dirichlet_series(zeta(s))
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + 1/(10^s) + 1/(11^s) + 1/(12^s) + 1/(13^s) + 1/(14^s) + 1/(15^s) + 1/(16^s) + 1/(17^s) + 1/(18^s) + 1/(19^s) + O(20^(-s))
            sage: dirichlet_series([1,0,1,0,1,0,1])
            1 + 1/(3^s) + 1/(5^s) + 1/(7^s) + O(8^(-s))
            sage: dirichlet_series([1,0,1,0,1,0,1], precision=4)
            1 + 1/(3^s) + O(4^(-s))
            sage: dirichlet_series(zeta(s-1))
            1 + 2/(2^s) + 3/(3^s) + 4/(4^s) + 5/(5^s) + 6/(6^s) + 7/(7^s)...
            sage: dirichlet_series(zeta(s)/(1-2^(-s)))
            1 + 1/(3^s) + 1/(5^s) + 1/(7^s) + 1/(9^s) + ...
            sage: dirichlet_series(dirichlet_L(4,2,s))
            1 + -1/(3^s) + 1/(5^s) + -1/(7^s) + 1/(9^s) + -1/(11^s)...

        ::

            sage: D1 = dirichlet_series([1,1,1,1]);
            sage: D1 + D1
            2 + 2/(2^s) + 2/(3^s) + 2/(4^s) + O(5^(-s))
            sage: D1 * D1
            1 + 2/(2^s) + 2/(3^s) + 3/(4^s) + O(5^(-s))
            sage: D2 = dirichlet_series(zeta(s), precision=8); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + O(8^(-s))
            sage: D2^2
            1 + 2/(2^s) + 2/(3^s) + 3/(4^s) + 2/(5^s) + 4/(6^s) + 2/(7^s) + O(8^(-s))
            sage: D2^(-1)
            1 + -1/(2^s) + -1/(3^s) + -1/(5^s) + 1/(6^s) + -1/(7^s) + O(8^(-s))
"""

#*****************************************************************************
#       Copyright (C) 2011 Jonathan Hanke
#       Copyright (C) 2015 Ralf Stephan <ralf@ark.in-berlin.de>
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

from sage.structure.sage_object import SageObject
from sage.structure.element import parent
from sage.interfaces.gp import gp
from sage.libs.pari.pari_instance import pari
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.real_mpfr import RR
from sage.misc.defaults import series_precision
from sage.symbolic.expression import Expression
from sage.symbolic.function import BuiltinFunction
from sage.symbolic.ring import SR
from sage.functions.transcendental import zeta

class DirichletLFunction(BuiltinFunction):
    """
    The Dirichlet L function is implemented as symbol for use
    in :class:`DirichletSeries`.
    """
    def __init__(self):
        """
        Initialize class.
        """
        BuiltinFunction.__init__(self, 'dirichlet_L', nargs=3,
                                 latex_name=r'\operatorname{L}',
                                 conversions={'mathematica': 'DirichletL'})

    def _eval_(self, modulus, index, arg, **kwargs):
        """
        EXAMPLES::

            sage: s = var('s')
            sage: dirichlet_L(3, 2, s)
            dirichlet_L(3, 2, s)
        """
        self._modulus = modulus
        self._index = index
        return None

dirichlet_L = DirichletLFunction()

class DirichletSeries(SageObject):
    """
    Class for manipulation of formal Dirichlet series.
    """
    def __init__(self, arg, precision=series_precision(), base_ring=ZZ):
        """
        Create a Dirichlet series from either a list of coefficients,
        or a generating function (g.f.) expression. The g.f. can
        be given as a product of zeta function, L-function, and other
        expressions.

        INPUT:

        - ``arg`` -- a list with coefficients, or a generating function
            as symbolic expression

        - ``precision`` -- the output precision in case of a generated
            series (default: global ``series_precision()``); equal to the
            list length if the series is created from a list

        - ``base_ring`` -- the coefficient ring

        The series is considered to be exact (has infinite absolute precision)
        if it is created from a generating function.
        If created from a list its output precision is equal to the relative
        precision, i.e., the list length. If not given, the precision defaults
        to the global series precision.

        EXAMPLES::

            sage: s,t = var('s,t')
            sage: dirichlet_series(1)
            1 + O(20^(-s))
            sage: dirichlet_series(zeta(s))
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + 1/(10^s) + 1/(11^s) + 1/(12^s) + 1/(13^s) + 1/(14^s) + 1/(15^s) + 1/(16^s) + 1/(17^s) + 1/(18^s) + 1/(19^s) + O(20^(-s))
            sage: dirichlet_series(zeta(t), precision=3)
            1 + 1/(2^t) + O(3^(-t))
            sage: dirichlet_series([1,0,1,0,1,0,1])
            1 + 1/(3^s) + 1/(5^s) + 1/(7^s) + O(8^(-s))
            sage: dirichlet_series([1,0,1,0,1,0,1], precision=4)
            1 + 1/(3^s) + O(4^(-s))

        Creation from g.f. is limited to products of ``zeta(A*s+B)``,
        ``dirichlet_L(C, D, s)``, and ``(1-E^(-s+F))``, with ``A,B,C,D,E,F`` positive
        integers. This also ensures the multiplicativity of series
        coefficients::

            sage: dirichlet_series(zeta(s-1))
            1 + 2/(2^s) + 3/(3^s) + 4/(4^s) + 5/(5^s) + 6/(6^s) + 7/(7^s)...
            sage: dirichlet_series(zeta(s-2)/zeta(2*s))
            1 + 4/(2^s) + 9/(3^s) + 15/(4^s) + 25/(5^s) + 36/(6^s) + ...
            sage: dirichlet_series(zeta(s)/(1-2^(-s)))
            1 + 1/(3^s) + 1/(5^s) + 1/(7^s) + 1/(9^s) + ...
            sage: L2a = dirichlet_series(dirichlet_L(2,1,s))
            sage: L2b = dirichlet_series(zeta(s)/(1-2^(-s)))
            sage: assert(L2a.list() == L2b.list())
            sage: dirichlet_series(zeta(s)*(1-2^(-s))/(1-2^(-s+1)))
            1 + 1/(3^s) + -1/(4^s) + 1/(5^s) + 1/(7^s) + -2/(8^s) + ...
            sage: dirichlet_series(dirichlet_L(4,2,s))
            1 + -1/(3^s) + 1/(5^s) + -1/(7^s) + 1/(9^s) + -1/(11^s)...

        TESTS::

            sage: dirichlet_series(zeta(5))
            Traceback (most recent call last):
            ...
            ValueError: Generating function must have exactly one variable
            sage: dirichlet_series(zeta(s)*zeta(t))
            Traceback (most recent call last):
            ...
            ValueError: Generating function must have exactly one variable
            sage: dirichlet_series(zeta(s)+1)
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot construct Dirichlet series ...
            sage: dirichlet_series(zeta(s)*(1-2^(-s))/(1-2^(-2*s)))
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot construct Dirichlet series ...
        """
        if not hasattr(base_ring, 'is_ring') or not base_ring.is_ring():
            raise TypeError("The base_ring argument must be a ring!")
        
        self._creation_function_expression = None
        self._var = SR.var('s')
        if not isinstance(arg, list):
            if SR(arg).is_integer():
                self._factors = [(base_ring(arg),1)]
            elif isinstance(arg, Expression):
                vars = arg.variables()
                if len(vars) != 1:
                    raise ValueError('Generating function must have exactly one variable')
                self._var = vars[0]

                # the specific form of the factors returned by Maxima
                # (see ticket #18081) makes masking necessary
                arg = dirichlet_series._mask_exp_factors(arg)
                self._factors = arg.factor_list()
                self._check_construction()
            else:
                raise NotImplementedError
            self._coeffs = []
            self._creation_function_expression = arg
        else:
            self._coeffs = [base_ring(x) for x in arg]
            self._creation_function_expression = None

        self._precision = precision
        self._base_ring = base_ring

    def _check_construction(self):
        for f,exp in self._factors:
            if not (SR(f).is_integer()
                or dirichlet_series._is_zeta_expr(f, self._var)
                or dirichlet_series._is_dirichlet_L_expr(f, self._var)
                or dirichlet_series._is_exp_expr(f, self._var)):
                raise NotImplementedError('Cannot construct Dirichlet series from the factor {}'.format(f))

    def _eval(self, precision):
        """
        Return a list with the coefficients of this series up to
        a given index (if existing).

        EXAMPLES::

            sage: s = var('s')
            sage: dirichlet_series(zeta(s-1))._eval(5)
            [1, 2, 3, 4, 5]
            sage: dirichlet_series(zeta(s-2)/zeta(2*s))._eval(5)
            [1, 4, 9, 15, 25]
            sage: dirichlet_series(zeta(s)/(1-2^(-s)))._eval(9)
            [1, 0, 1, 0, 1, 0, 1, 0, 1]
        """
        if not self.has_infinite_precision():
            return self._coeffs
        product = [self._base_ring(1)] + [0] * (precision-1)
        for f,exp in self._factors:
            if SR(f).is_integer():
                fseries = [self._base_ring(f)] + [0] * (precision-1)
            elif dirichlet_series._is_zeta_expr(f, self._var):
                fseries = dirichlet_series._zeta_expr(f, self._var, precision)
            elif dirichlet_series._is_dirichlet_L_expr(f, self._var):
                fseries = dirichlet_series._dirichlet_L_expr(f, precision, self._base_ring)
            elif dirichlet_series._is_exp_expr(f, self._var):
                fseries = dirichlet_series._exp_expr(f, precision)
            else:
                raise NotImplementedError('Cannot construct exact Dirichlet series from the factor {}'.format(f))
            if exp > 0:
                for i in range(exp):
                    product = pari(product).dirmul(fseries).sage()
            else:
                for i in range(-exp):
                    product = pari(product).dirdiv(fseries).sage()
        return product

    class MaskFunction(BuiltinFunction):
        def __init__(self):
            BuiltinFunction.__init__(self, 'dummy', nargs=2)

    @staticmethod
    def _mask_exp_factors(ex):
        """
        EXAMPLES::

            sage: s = var('s')
            sage: dirichlet_series._mask_exp_factors(1-2^(-s))
            dummy(2, -s)
            sage: dirichlet_series._mask_exp_factors((1-3^(-s))^2*(1-5^(-s+1)))
            dummy(5, -s + 1)*dummy(3, -s)^2
            sage: dirichlet_series._mask_exp_factors(1/(1-2^(-s))^2*(1-3^(-s+2)))
            dummy(3, -s + 2)/dummy(2, -s)^2
        """
        w0 = SR.wild(0); w1 = SR.wild(1)
        ex = ex.subs(1-w0**w1 == dirichlet_series.MaskFunction()(w0, w1))
        ex = ex.subs(w0**w1-1 == -dirichlet_series.MaskFunction()(w0, w1))
        return ex

    @staticmethod
    def _is_zeta_expr(ex, var):
        """
        Return True if ``ex`` is of form ``dirichlet_L(a,b,var)``, with
        ``a>1``,``b>0`` positive integers.

        EXAMPLES::

            sage: dirichlet_series._is_zeta_expr(sin(x),x)
            False
            sage: dirichlet_series._is_zeta_expr(zeta(x),x)
            True
            sage: dirichlet_series._is_zeta_expr(zeta(2*x),x)
            True
            sage: dirichlet_series._is_zeta_expr(zeta(-2*x),x)
            False
            sage: dirichlet_series._is_zeta_expr(zeta(x-1),x)
            True
            sage: dirichlet_series._is_zeta_expr(zeta(2*x-1),x)
            True
            sage: dirichlet_series._is_zeta_expr(zeta(-x+.5),x)
            False

        TESTS::

            sage: dirichlet_series._is_zeta_expr(zeta(x,1),x)
            Traceback (most recent call last):
            ...
            TypeError: Symbolic function zeta takes exactly 1 arguments (2 given)
        """
        if (ex.operator() == zeta and len(ex.operands()) == 1
            and len(ex.variables()) == 1):
            arg = ex.operands()[0]
            if bool(arg == 1) or bool(arg == var):
                return True
            w0 = SR.wild(0); w1 = SR.wild(1)
            d = arg.match(w0*var)
            if d is not None:
                factor = d.values()[0]
                return factor.is_integer() and bool(factor > 1)
            d = arg.match(var+w0)
            if d is not None:
                summand = d.values()[0]
                return summand.is_integer()
            d = arg.match(w0*var+w1)
            if d is not None:
                factor = d.get(w0)
                summand = d.get(w1)
                return factor.is_integer() and bool(factor > 1) and summand.is_integer()
        return False

    @staticmethod
    def _is_dirichlet_L_expr(ex, var):
        """
        Return True if ``ex`` is of form ``dirichlet_L(a,b,var)``, with
        ``a>1``,``b>0`` positive integers.

        EXAMPLES::

            sage: dirichlet_series._is_dirichlet_L_expr(sin(x),x)
            False
            sage: dirichlet_series._is_dirichlet_L_expr(dirichlet_L(2,1,x+1),x)
            False
            sage: dirichlet_series._is_dirichlet_L_expr(dirichlet_L(2,1,x),x)
            True
        """
        # Most checks happened already in class DirichletLFunction
        return ex.operator() == dirichlet_L and bool(ex.operands()[2] == var)

    @staticmethod
    def _is_exp_expr(ex, v):
        from sage.symbolic.operators import add_vararg
        w0 = SR.wild(0); w1 = SR.wild(1)
        dummy = dirichlet_series.MaskFunction()
        d = ex.match(dummy(w0, w1))
        if d is None:
            return False
        var = ex.variables()[0]
        arg2 = d.get(w1)
        return (d is not None
                and (arg2 == -var or arg2.operator() == add_vararg))

    @staticmethod
    def _Lseries_coeff(m, r, prec, R=ZZ):
        """
        Return the coefficients of the Dirichlet series generated
        by the L-function with given modulus ``m`` and representation
        ``r``, with precision ``prec``.

        EXAMPLES::

            sage: dirichlet_series._Lseries_coeff(3, 2, 10)
            [1, -1, 0, 1, -1, 0, 1, -1, 0, 1]
            sage: dirichlet_series._Lseries_coeff(11, 3, 5, CyclotomicField(10))
            [1, zeta10^2, -zeta10, zeta10^3 - zeta10^2 + zeta10 - 1, -zeta10^3]

        REFERENCES:

        .. [MatharTable2010] R. J. Mathar, Table of Dirichlet L-Series and Prime
             Zeta Modulo functions for small moduli, :arxiv:`1008.2547`
        """
        from sage.modular.dirichlet import DirichletGroup
        dg = DirichletGroup(m)
        period = dg.list()[r-1].values()
        l = len(period)
        period = [R(elem) for elem in period]
        period = period[1:] + [R(0)]
        quo,rem = Integer(prec).quo_rem(l)
        return period * quo + period[:rem]

    @staticmethod
    def _zeta_expr(ex, var, prec):
        """
        Return the coefficients of the Dirichlet series generated
        by the zeta function expression with variable ``var`` and
        precision ``prec``.

        EXAMPLES::

            sage: s = var('s')
            sage: dirichlet_series._zeta_expr(zeta(s), s, 5)
            [1, 1, 1, 1, 1]
            sage: dirichlet_series._zeta_expr(zeta(s-1), s, 5)
            [1, 2, 3, 4, 5]
            sage: dirichlet_series._zeta_expr(zeta(2*s-1), s, 5)
            [1, 0, 0, 2, 0]
            sage: dirichlet_series._zeta_expr(zeta(2*s), s, 5)
            [1, 0, 0, 1, 0]
        """
        arg = ex.operands()[0]
        if bool(arg == var):
            return [1] * prec
        w0 = SR.wild(0); w1 = SR.wild(1)
        d = arg.match(w0*var)
        if d is not None:
            factor = d.values()[0]
            paricmd = "direuler(p=1,{0},1/(1-X^{1}))".format(prec, factor)
            return gp(paricmd)
        d = arg.match(var+w0)
        if d is not None:
            summand = d.values()[0]
            paricmd = "direuler(p=1,{0},1/(1-p^{1}*X))".format(prec, -summand)
            return gp(paricmd)
        d = arg.match(w0*var+w1)
        if d is not None:
            factor = d.get(w0)
            summand = d.get(w1)
            paricmd = "direuler(p=1,{0},1/(1-p^{1}*X^{2}))".format(prec, -summand, factor)
            return gp(paricmd)

    @staticmethod
    def _dirichlet_L_expr(ex, prec, R=ZZ):
        """
        Return the coefficients of the Dirichlet series generated
        by the L-function expression with precision ``prec``.

        EXAMPLES::

            sage: s = var('s')
            sage: dirichlet_series._dirichlet_L_expr(dirichlet_L(2,0,s), 9, ZZ)
            [1, 0, 1, 0, 1, 0, 1, 0, 1]
            sage: dirichlet_series._dirichlet_L_expr(dirichlet_L(3,0,s), 9, ZZ)
            [1, -1, 0, 1, -1, 0, 1, -1, 0]
            sage: dirichlet_series._dirichlet_L_expr(dirichlet_L(3,1,s), 9, ZZ)
            [1, 1, 0, 1, 1, 0, 1, 1, 0]
            sage: dirichlet_series._dirichlet_L_expr(dirichlet_L(9,0,s), 9, CyclotomicField(9))
            [1, -zeta9^3, 0, -zeta9^3 - 1, zeta9^3 + 1, 0, zeta9^3, -1, 0]
        """
        m,r,_ = ex.operands()
        return dirichlet_series._Lseries_coeff(m, r, prec, R)

    @staticmethod
    def _exp_expr(ex, prec):
        """
        Return the coefficients of the Dirichlet series generated
        by the dummy expression with precision ``prec``.

        EXAMPLES::

            sage: s = var('s')
            sage: dummy = dirichlet_series.MaskFunction()
            sage: dirichlet_series._exp_expr(dummy(2, -s), 9)
            [1, 1, 0, 1, 0, 0, 0, 1, 0]
            sage: dirichlet_series._exp_expr(dummy(2, -s+1), 9)
            [1, 2, 0, 4, 0, 0, 0, 8, 0]
            sage: dirichlet_series._exp_expr(dummy(3, -s+1), 9)
            [1, 0, 3, 0, 0, 0, 0, 0, 9]
        """
        w0 = SR.wild(0); w1 = SR.wild(1)
        dummy = dirichlet_series.MaskFunction()
        d1 = ex.match(dummy(w0, w1))
        arg1 = d1.get(w0)
        arg2 = d1.get(w1)
        var = ex.variables()[0]
        if arg2 == -var:
            paricmd = "direuler(p=1,{0},if(p=={1},1/(1-X),1))".format(prec, arg1)
            return gp(paricmd)
        d2 = arg2.match(-var+w0)
        summand = d2.values()[0]
        paricmd = "direuler(p=1,{0},if(p=={1},1/(1-p^{2}*X),1))".format(prec, arg1, summand)
        return gp(paricmd)

    def base_ring(self):
        """
        Return the base ring of the Dirichlet series.

        EXAMPLES::

            sage: dirichlet_series(1).base_ring()
            Integer Ring
        """
        return self._base_ring

    def has_infinite_precision(self):
        """
        Return if the Dirichlet series was generated from a generating function.

        EXAMPLES::

            sage: dirichlet_series(1).has_infinite_precision()
            True
            sage: dirichlet_series([1]).has_infinite_precision()
            False
        """
        return len(self._coeffs) == 0

    def list(self, n_prec=series_precision()):
        """
        Return the list of coefficients of the Dirichlet series,
        optionally up to a specific index (default: global
        ``series_precision()``.

        EXAMPLES::

            sage: s = var('s')
            sage: D = dirichlet_series(dirichlet_L(3,2,s))
            sage: D.list()
            [1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1]
            sage: D.list(5)
            [1, -1, 0, 1, -1]
        """
        if self.has_infinite_precision():
            return self._eval(n_prec)
        else:
            l = min(len(self._coeffs), n_prec)
            return self._coeffs[:l]

    def __repr__(self, n_max=Infinity):
        """
        Print first n coefficients of the Dirichlet series,
        using its default variable.
        """
        out_str = ""
        n_prec = min(self._precision, n_max+1) - 1
        if self.has_infinite_precision():
            coeffs = self._eval(n_prec)
        else:
            coeffs = self._coeffs
            n_prec = min(n_prec, len(self._coeffs))

        ## Add the first coefficient
        if n_prec >= 1:
            out_str += str(coeffs[0])

        ## Add all other non-zero coefficients up to the desired precision
        for i in range(1, n_prec):
            if (coeffs[i] != 0):
                out_str += " + " + str(coeffs[i]) + "/(" + str(i+1) + "^" + str(self._var) + ")"
            
        ## Add the error term, if the result has infinite precision.
        out_str += " + O(" + str(n_prec + 1) + "^(-" + str(self._var) + "))"
            
        ## Return the output string
        return out_str

    def __add__(self, other):
        """
        Form the sum of two Dirichlet series.

        Addition makes the precision finite because sums are
        not supported as generating functions at the moment.

        EXAMPLES::

            sage: D1 = dirichlet_series([1,1,1,1]);
            sage: D1 + D1
            2 + 2/(2^s) + 2/(3^s) + 2/(4^s) + O(5^(-s))
            sage: D1 + D1 + D1
            3 + 3/(2^s) + 3/(3^s) + 3/(4^s) + O(5^(-s))
            sage: s = var('s')
            sage: D2 = dirichlet_series(zeta(s), precision=10); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + O(10^(-s))
            sage: D1 + D2
            2 + 2/(2^s) + 2/(3^s) + 2/(4^s) + O(5^(-s))

        TESTS::

            sage: D1 = dirichlet_series([1,1,1])
            sage: D1 + 1
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot add DirichletSeries to 1
        """
        if not isinstance(other, DirichletSeries):
            raise NotImplementedError('cannot add DirichletSeries to {}'.format(other))
        if self.base_ring() != other.base_ring():
            raise NotImplementedError("addition of DirichletSeries having different base rings")

        if self.has_infinite_precision():
            self_prec = Infinity
        else:
            self_prec = self._precision
        if other.has_infinite_precision():
            other_prec = Infinity
        else:
            other_prec = other._precision

        prec = min(self_prec, other_prec)
        if prec == Infinity:
            prec = min(self._precision, other._precision)

        if self.has_infinite_precision():
            self_list = self._eval(prec)
        else:
            self_list = self._coeffs
        if other.has_infinite_precision():
            other_list = other._eval(prec)
        else:
            other_list = other._coeffs

        sum_coeff_list = [self_list[i] + other_list[i] for i
                          in range(min(len(self_list), len(other_list)))]
        D_sum = dirichlet_series(sum_coeff_list)
        return D_sum

    def __sub__(self, other):
        """
        Form the sum of two Dirichlet series.

        EXAMPLES::

            sage: s = var('s')
            sage: D1 = dirichlet_series([1,1,1,1]);
            sage: D1 - D1
            0 + O(5^(-s))
            sage: D2 = dirichlet_series(zeta(s), precision=10); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + 1/(8^s) + 1/(9^s) + O(10^(-s))
            sage: D1 - D2
            0 + O(5^(-s))

        """
        return self + (other * ZZ(-1))

    def __mul__(self, other):
        """
        Define the product of two Dirichlet series, or of a Dirichlet series and a number.

        EXAMPLES::

            sage: s = var('s')
            sage: D1 = dirichlet_series([1,1,1,1]); D1
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + O(5^(-s))
            sage: D1 * D1
            1 + 2/(2^s) + 2/(3^s) + 3/(4^s) + O(5^(-s))
            sage: D2 = dirichlet_series(zeta(s), precision=8); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + O(8^(-s))
            sage: D22 = D2 * D2; D22
            1 + 2/(2^s) + 2/(3^s) + 3/(4^s) + 2/(5^s) + 4/(6^s) + 2/(7^s) + O(8^(-s))

        TESTS::

            sage: D1 = dirichlet_series([1,1,1])
            sage: D1 * 0.5
            Traceback (most recent call last):
            ...
            ValueError: cannot multiply DirichletSeries with 0.500000000000000
        """
        if not isinstance(other, DirichletSeries) and not parent(other) is self.base_ring():
            raise ValueError('cannot multiply DirichletSeries with {}'.format(other))
        if isinstance(other, DirichletSeries):
            if self.base_ring() != other.base_ring():
                raise NotImplementedError("multiplication of DirichletSeries having different base rings")
            if self.has_infinite_precision() and other.has_infinite_precision():
                new_prec = min(self._precision, other._precision)
                return dirichlet_series(self._creation_function_expression
                                        * other._creation_function_expression,
                                        precision=new_prec,
                                        base_ring=self._base_ring)
            else:
                new_coeff_list = pari(self.list()).dirmul(other.list()).sage()
                return dirichlet_series(new_coeff_list)
        else:
            scale_factor = self._base_ring(other)
            if self.has_infinite_precision():
                return dirichlet_series(self._creation_function_expression * scale_factor)
            else:
                return dirichlet_series([scale_factor * self._coeffs[i]
                                     for i in range(len(self._coeffs))])

    def __pow__(self, n):
        """
        Take any integer power of the current Dirichlet series.

        EXAMPLES::

            sage: s = var('s')
            sage: D2 = dirichlet_series(zeta(s), precision=8); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + O(8^(-s))
            sage: D2^2
            1 + 2/(2^s) + 2/(3^s) + 3/(4^s) + 2/(5^s) + 4/(6^s) + 2/(7^s) + O(8^(-s))
            sage: D2^5
            1 + 5/(2^s) + 5/(3^s) + 15/(4^s) + 5/(5^s) + 25/(6^s) + 5/(7^s) + O(8^(-s))
            sage: D2^0
            1 + O(8^(-s))
            sage: D2^(-1)
            1 + -1/(2^s) + -1/(3^s) + -1/(5^s) + 1/(6^s) + -1/(7^s) + O(8^(-s))

        """
        if not SR(n).is_integer():
            raise TypeError("The power must be an integer!")
        if self.has_infinite_precision():
            if n == 0:
                return dirichlet_series(1,
                                    precision=self._precision,
                                    base_ring=self._base_ring)
            else:
                return dirichlet_series(self._creation_function_expression ** n,
                                    precision=self._precision,
                                    base_ring=self._base_ring)
        else:
            if n == 0:
                return dirichlet_series([1] + [0]*(len(self._coeffs)-1))
            elif n > 0:
                tmp_new_series = self
                for i in range(n-1):
                    tmp_new_series = tmp_new_series * self
            elif n < 0:
                inv_series = self.inverse()
                tmp_new_series = inv_series
                for i in range(abs(n)-1):
                    tmp_new_series = tmp_new_series * inv_series
            return tmp_new_series

    def inverse(self):
        """
        Compute the inverse Dirichlet series under Dirichlet multiplication.

        EXAMPLES::

            sage: s = var('s')
            sage: D2 = dirichlet_series(zeta(s), precision=8); D2
            1 + 1/(2^s) + 1/(3^s) + 1/(4^s) + 1/(5^s) + 1/(6^s) + 1/(7^s) + O(8^(-s))
            sage: D2.inverse()
            1 + -1/(2^s) + -1/(3^s) + -1/(5^s) + 1/(6^s) + -1/(7^s) + O(8^(-s))
            sage: D2 * D2.inverse()
            1 + O(8^(-s))

        TESTS::

            sage: D1 = dirichlet_series([2,1,1])
            sage: D1^(-1)
            Traceback (most recent call last):
            ...
            RuntimeError: The leading term is not invertible ...
        """
        from sage.rings.arith import divisors
        R = self.base_ring()

        if self.has_infinite_precision():
            return self**(-1)
        else:
            ## Check if the first coefficient is invertible in R
            a1 = self._coeffs[0]
            a1_inv = a1**(-1)
            if not a1_inv in R:
                raise RuntimeError("The leading term is not invertible in the base ring, so the Dirichlet series is not invertible.")
            new_coeff_list = pari([1]+[0]*len(self._coeffs)).dirdiv(self.list()).sage()
            return dirichlet_series(new_coeff_list)

dirichlet_series = DirichletSeries
