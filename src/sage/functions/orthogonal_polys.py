r"""
Orthogonal Polynomials

This module wraps some of the orthogonal/special functions in the
Maxima package "orthopoly". This package was written by Barton
Willis of the University of Nebraska at Kearney. It is released
under the terms of the General Public License (GPL). Send
Maxima-related bug reports and comments on this module to
willisb@unk.edu. In your report, please include Maxima and specfun
version information.


-  The Chebyshev polynomial of the first kind arises as a solution
   to the differential equation

   .. math::

         (1-x^2)\,y'' - x\,y' + n^2\,y = 0


   and those of the second kind as a solution to

   .. math::

         (1-x^2)\,y'' - 3x\,y' + n(n+2)\,y = 0.


   The Chebyshev polynomials of the first kind are defined by the
   recurrence relation

   .. math::

     T_0(x) = 1 \, T_1(x) = x \, T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x). \,


   The Chebyshev polynomials of the second kind are defined by the
   recurrence relation

   .. math::

     U_0(x) = 1 \, U_1(x) = 2x \, U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x). \,



   For integers `m,n`, they satisfy the orthogonality
   relations

   .. math::

     \int_{-1}^1 T_n(x)T_m(x)\,\frac{dx}{\sqrt{1-x^2}} =\left\{ \begin{matrix} 0 &: n\ne m~~~~~\\ \pi &: n=m=0\\ \pi/2 &: n=m\ne 0 \end{matrix} \right.


   and


   .. math::

     \int_{-1}^1 U_n(x)U_m(x)\sqrt{1-x^2}\,dx =\frac{\pi}{2}\delta_{m,n}.



   They are named after Pafnuty Chebyshev (alternative
   transliterations: Tchebyshef or Tschebyscheff).

-  The Hermite polynomials are defined either by

   .. math::

     H_n(x)=(-1)^n e^{x^2/2}\frac{d^n}{dx^n}e^{-x^2/2}


   (the "probabilists' Hermite polynomials"), or by


   .. math::

     H_n(x)=(-1)^n e^{x^2}\frac{d^n}{dx^n}e^{-x^2}


   (the "physicists' Hermite polynomials"). Sage (via Maxima)
   implements the latter flavor. These satisfy the orthogonality
   relation

   .. math::

     \int_{-\infty}^\infty H_n(x)H_m(x)\,e^{-x^2}\,dx ={n!2^n}{\sqrt{\pi}}\delta_{nm}



   They are named in honor of Charles Hermite.

-  Each *Legendre polynomial* `P_n(x)` is an `n`-th degree polynomial.
   It may be expressed using Rodrigues' formula:

   .. math::

      P_n(x) = (2^n n!)^{-1} {\frac{d^n}{dx^n} } \left[ (x^2 -1)^n \right].

   These are solutions to Legendre's differential equation:

   .. math::

      {\frac{d}{dx}} \left[ (1-x^2) {\frac{d}{dx}} P(x) \right] + n(n+1)P(x) = 0.

   and satisfy the orthogonality relation

   .. math::

      \int_{-1}^{1} P_m(x) P_n(x)\,dx = {\frac{2}{2n + 1}} \delta_{mn}

   The *Legendre function of the second kind* `Q_n(x)` is another
   (linearly independent) solution to the Legendre differential equation.
   It is not an "orthogonal polynomial" however.

   The associated Legendre functions of the first kind
   `P_\ell^m(x)` can be given in terms of the "usual"
   Legendre polynomials by

   .. math::

     \begin{array}{ll} P_\ell^m(x)    &=  (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}P_\ell(x) \\ &=  \frac{(-1)^m}{2^\ell \ell!} (1-x^2)^{m/2}\frac{d^{\ell+m}}{dx^{\ell+m}}(x^2-1)^\ell. \end{array}


   Assuming `0 \le m \le \ell`, they satisfy the orthogonality
   relation:

   .. math::

      \int_{-1}^{1} P_k ^{(m)} P_\ell ^{(m)} dx  = \frac{2 (\ell+m)!}{(2\ell+1)(\ell-m)!}\ \delta _{k,\ell},


   where `\delta _{k,\ell}` is the Kronecker delta.

   The associated Legendre functions of the second kind
   `Q_\ell^m(x)` can be given in terms of the "usual"
   Legendre polynomials by


   .. math::

     Q_\ell^m(x)   =  (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}Q_\ell(x).



   They are named after Adrien-Marie Legendre.

-  Laguerre polynomials may be defined by the Rodrigues formula

   .. math::

      L_n(x)=\frac{e^x}{n!}\frac{d^n}{dx^n}\left(e^{-x} x^n\right).


   They are solutions of Laguerre's equation:


   .. math::

      x\,y'' + (1 - x)\,y' + n\,y = 0\,

   and satisfy the orthogonality relation


   .. math::

      \int_0^\infty L_m(x) L_n(x) e^{-x}\,dx = \delta_{mn}.



   The generalized Laguerre polynomials may be defined by the
   Rodrigues formula:


   .. math::

       L_n^{(\alpha)}(x)   = {\frac{x^{-\alpha} e^x}{n!}}{\frac{d^n}{dx^n}} \left(e^{-x} x^{n+\alpha}\right) .


   (These are also sometimes called the associated Laguerre
   polynomials.) The simple Laguerre polynomials are recovered from
   the generalized polynomials by setting `\alpha =0`.

   They are named after Edmond Laguerre.

-  Jacobi polynomials are a class of orthogonal polynomials. They
   are obtained from hypergeometric series in cases where the series
   is in fact finite:

   .. math::

     P_n^{(\alpha,\beta)}(z) =\frac{(\alpha+1)_n}{n!} \,_2F_1\left(-n,1+\alpha+\beta+n;\alpha+1;\frac{1-z}{2}\right) ,


   where `()_n` is Pochhammer's symbol (for the rising
   factorial), (Abramowitz and Stegun p561.) and thus have the
   explicit expression


   .. math::

     P_n^{(\alpha,\beta)} (z) = \frac{\Gamma (\alpha+n+1)}{n!\Gamma (\alpha+\beta+n+1)} \sum_{m=0}^n {n\choose m} \frac{\Gamma (\alpha + \beta + n + m + 1)}{\Gamma (\alpha + m + 1)} \left(\frac{z-1}{2}\right)^m .



   They are named after Carl Jacobi.

-  Ultraspherical or Gegenbauer polynomials are given in terms of
   the Jacobi polynomials `P_n^{(\alpha,\beta)}(x)` with
   `\alpha=\beta=a-1/2` by


   .. math::

     C_n^{(a)}(x)= \frac{\Gamma(a+1/2)}{\Gamma(2a)}\frac{\Gamma(n+2a)}{\Gamma(n+a+1/2)} P_n^{(a-1/2,a-1/2)}(x).


   They satisfy the orthogonality relation

   .. math::

     \int_{-1}^1(1-x^2)^{a-1/2}C_m^{(a)}(x)C_n^{(a)}(x)\, dx =\delta_{mn}2^{1-2a}\pi \frac{\Gamma(n+2a)}{(n+a)\Gamma^2(a)\Gamma(n+1)} ,


   for `a>-1/2`. They are obtained from hypergeometric series
   in cases where the series is in fact finite:


   .. math::

     C_n^{(a)}(z) =\frac{(2a)^{\underline{n}}}{n!} \,_2F_1\left(-n,2a+n;a+\frac{1}{2};\frac{1-z}{2}\right)


   where `\underline{n}` is the falling factorial. (See
   Abramowitz and Stegun p561)

   They are named for Leopold Gegenbauer (1849-1903).


For completeness, the Pochhammer symbol, introduced by Leo August
Pochhammer, `(x)_n`, is used in the theory of special
functions to represent the "rising factorial" or "upper factorial"

.. math::

         (x)_n=x(x+1)(x+2)\cdots(x+n-1)=\frac{(x+n-1)!}{(x-1)!}.


On the other hand, the "falling factorial" or "lower factorial" is

.. math::

     x^{\underline{n}}=\frac{x!}{(x-n)!} ,


in the notation of Ronald L. Graham, Donald E. Knuth and Oren
Patashnik in their book Concrete Mathematics.

.. note::

   The first call of any of these will usually cost a bit extra
   (it loads "specfun", but I'm not sure if that is the real reason).
   The next call is usually faster but not always.

.. TODO::

    Implement associated Legendre polynomials and Zernike
    polynomials. (Neither is in Maxima.)
    :wikipedia:`Associated_Legendre_polynomials`
    :wikipedia:`Zernike_polynomials`

REFERENCES:

.. [ASHandbook] Abramowitz and Stegun: Handbook of Mathematical Functions,
    http://www.math.sfu.ca/ cbm/aands/

.. :wikipedia:`Chebyshev_polynomials`

.. :wikipedia:`Legendre_polynomials`

.. :wikipedia:`Hermite_polynomials`

.. http://mathworld.wolfram.com/GegenbauerPolynomial.html

.. :wikipedia:`Jacobi_polynomials`

.. :wikipedia:`Laguerre_polynomia`

.. :wikipedia:`Associated_Legendre_polynomials`

.. [EffCheby] Wolfram Koepf: Effcient Computation of Chebyshev Polynomials
    in Computer Algebra
    Computer Algebra Systems: A Practical Guide.
    John Wiley, Chichester (1999): 79-99.

AUTHORS:

- David Joyner (2006-06)
- Stefan Reiterer (2010-)
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdj@usna.edu>
#                     2010 Stefan Reiterer <maldun.finsterschreck@gmail.com>
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

import warnings

from sage.misc.sage_eval import sage_eval
from sage.rings.all import ZZ, QQ, RR, CC
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.real_mpfr import is_RealField
from sage.rings.complex_field import is_ComplexField
from sage.calculus.calculus import maxima


from sage.symbolic.ring import SR, is_SymbolicVariable
from sage.symbolic.function import BuiltinFunction
from sage.symbolic.expression import Expression
from sage.functions.other import factorial, binomial
from sage.structure.coerce import parent

_done = False
def _init():
    """
    Internal function which checks if Maxima has loaded the
    "orthopoly" package.  All functions using this in this
    file should call this function first.

    TEST:

    The global starts ``False``::

        sage: sage.functions.orthogonal_polys._done
        False

    Then after using one of these functions, it changes::

        sage: from sage.functions.orthogonal_polys import laguerre
        sage: laguerre(2,x)
        1/2*x^2 - 2*x + 1
        sage: sage.functions.orthogonal_polys._done
        True


    Note that because here we use a Pynac variable ``x``,
    the representation of the function is different from
    its actual doctest, where a polynomial indeterminate
    ``x`` is used.
    """
    global _done
    if _done:
        return
    maxima.eval('load("orthopoly");')
    # TODO -- make it possible to use the intervals returned
    # instead of just discarding this info!
    maxima.eval('orthopoly_returns_intervals:false;')
    _done = True


class OrthogonalPolynomial(BuiltinFunction):
    """
    Base class for orthogonal polynomials.

    This class is an abstract base class for all orthogonal polynomials since
    they share similar properties. The evaluation as a polynomial
    is either done via maxima, or with pynac.

    Convention: The first argument is always the order of the polynomial,
    the others are other values or parameters where the polynomial is
    evaluated.
    """
    def __init__(self, name, nargs=2, latex_name=None, conversions={}):
        """
        :class:`OrthogonalPolynomial` class needs the same input parameter as
        it's parent class.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import OrthogonalPolynomial
            sage: new = OrthogonalPolynomial('testo_P')
            sage: new
            testo_P
        """
        try:
            self._maxima_name = conversions['maxima']
        except KeyError:
            self._maxima_name = None

        super(OrthogonalPolynomial,self).__init__(name=name, nargs=nargs,
                                 latex_name=latex_name, conversions=conversions)

    def _maxima_init_evaled_(self, *args):
        r"""
        Return a string which represents this function evaluated at
        ``n, x`` in Maxima.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import OrthogonalPolynomial
            sage: P = OrthogonalPolynomial('testo_P')
            sage: P._maxima_init_evaled_(2, 5) is None
            True
        """
        return None

    def eval_formula(self, *args):
        """
        Evaluate this polynomial using an explicit formula.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import OrthogonalPolynomial
            sage: P = OrthogonalPolynomial('testo_P')
            sage: P.eval_formula(1,2.0)
            Traceback (most recent call last):
            ...
            NotImplementedError: no explicit calculation of values implemented
        """
        raise NotImplementedError("no explicit calculation of values implemented")

    def _eval_special_values_(self, *args):
        """
        Evaluate the polynomial explicitly for special values.

        EXAMPLES::

            sage: var('n')
            n
            sage: chebyshev_T(n,-1)
            (-1)^n
        """
        raise ValueError("no special values known")

    def _eval_(self, n, *args):
        """
        The :meth:`_eval_()` method decides which evaluation suits best
        for the given input, and returns a proper value.

        EXAMPLES::

            sage: var('n,x')
            (n, x)
            sage: chebyshev_T(5,x)
            16*x^5 - 20*x^3 + 5*x
        """
        return None

    def __call__(self, n, *args, **kwds):
        """
        This overides the call method from SageObject to avoid problems with coercions,
        since the _eval_ method is able to handle more data types than symbolic functions
        would normally allow.
        Thus we have the distinction between algebraic objects (if n is an integer),
        and else as symbolic function.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-x-1)
            sage: chebyshev_T(5, a)
            16*a^2 + a - 4
        """
        return super(OrthogonalPolynomial,self).__call__(n, *args, **kwds)

class ChebyshevPolynomial(OrthogonalPolynomial):
    """
    Abstract base class for Chebyshev polynomials of the first and second kind.

    EXAMPLES::

        sage: chebyshev_T(3,x)
        4*x^3 - 3*x
    """
    def __call__(self, n, *args, **kwds):
        """
        This overides the call method from SageObject to avoid problems with coercions,
        since the _eval_ method is able to handle more data types than symbolic functions
        would normally allow.
        Thus we have the distinction between algebraic objects (if n is an integer),
        and else as symbolic function.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-x-1)
            sage: chebyshev_T(5, a)
            16*a^2 + a - 4
            sage: chebyshev_T(5,MatrixSpace(ZZ, 2)([1, 2, -4, 7]))
            [-40799  44162]
            [-88324  91687]
            sage: R.<x> = QQ[]
            sage: parent(chebyshev_T(5, x))
            Univariate Polynomial Ring in x over Rational Field
            sage: chebyshev_T(5, 2, hold=True)
            chebyshev_T(5, 2)
            sage: chebyshev_T(1,2,3)
            Traceback (most recent call last):
            ...
            TypeError: Symbolic function chebyshev_T takes exactly 2 arguments (3 given)
        """
        # If n is an integer: consider the polynomial as an algebraic (not symbolic) object
        if n in ZZ and not kwds.get('hold', False):
            try:
                return self._eval_(n, *args)
            except Exception:
                pass

        return super(ChebyshevPolynomial,self).__call__(n, *args, **kwds)

    def _eval_(self, n, x):
        """
        The :meth:`_eval_()` method decides which evaluation suits best
        for the given input, and returns a proper value.

        EXAMPLES::

            sage: var('n,x')
            (n, x)
            sage: chebyshev_T(5,x)
            16*x^5 - 20*x^3 + 5*x
            sage: chebyshev_T(64, x)
            2*(2*(2*(2*(2*(2*x^2 - 1)^2 - 1)^2 - 1)^2 - 1)^2 - 1)^2 - 1
            sage: chebyshev_T(n,-1)
            (-1)^n
            sage: chebyshev_T(-7,x)
            64*x^7 - 112*x^5 + 56*x^3 - 7*x
            sage: chebyshev_T(3/2,x)
            chebyshev_T(3/2, x)
            sage: R.<t> = QQ[]
            sage: chebyshev_T(2,t)
            2*t^2 - 1
            sage: chebyshev_U(2,t)
            4*t^2 - 1
            sage: parent(chebyshev_T(4, RIF(5)))
            Real Interval Field with 53 bits of precision
            sage: RR2 = RealField(5)
            sage: chebyshev_T(100000,RR2(2))
            8.9e57180
            sage: chebyshev_T(5,Qp(3)(2))
            2 + 3^2 + 3^3 + 3^4 + 3^5 + O(3^20)
            sage: chebyshev_T(100001/2, 2)
            doctest:...: RuntimeWarning: mpmath failed, keeping expression unevaluated
            chebyshev_T(100001/2, 2)
            sage: chebyshev_U._eval_(1.5, Mod(8,9)) is None
            True
        """
        # n is an integer => evaluate algebraically (as polynomial)
        if n in ZZ:
            n = ZZ(n)
            # Expanded symbolic expression only for small values of n
            if isinstance(x, Expression) and n.abs() < 32:
                return self.eval_formula(n, x)
            return self.eval_algebraic(n, x)

        if isinstance(x, Expression) or isinstance(n, Expression):
            # Check for known identities
            try:
                return self._eval_special_values_(n, x)
            except ValueError:
                # Don't evaluate => keep symbolic
                return None

        # n is not an integer and neither n nor x is symbolic.
        # We assume n and x are real/complex and evaluate numerically
        try:
            import sage.libs.mpmath.all as mpmath
            return self._evalf_(n, x)
        except mpmath.NoConvergence:
            warnings.warn("mpmath failed, keeping expression unevaluated",
                          RuntimeWarning)
            return None
        except Exception:
            # Numerical evaluation failed => keep symbolic
            return None

    
class Func_chebyshev_T(ChebyshevPolynomial):
    """
    Chebyshev polynomials of the first kind.

    REFERENCE:

    - [ASHandbook]_ 22.5.31 page 778 and 6.1.22 page 256.

    EXAMPLES::

       sage: chebyshev_T(5,x)
       16*x^5 - 20*x^3 + 5*x
       sage: var('k')
       k
       sage: test = chebyshev_T(k,x)
       sage: test
       chebyshev_T(k, x)
    """
    def __init__(self):
        """
        Init method for the chebyshev polynomials of the first kind.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import Func_chebyshev_T
            sage: chebyshev_T2 = Func_chebyshev_T()
            sage: chebyshev_T2(1,x)
            x
        """
        ChebyshevPolynomial.__init__(self, "chebyshev_T", nargs=2,
                                     conversions=dict(maxima='chebyshev_t',
                                                      mathematica='ChebyshevT'))

    def _eval_special_values_(self, n, x):
        """
        Values known for special values of x.
        For details see [ASHandbook]_ 22.4 (p. 777)

        EXAMPLES:

            sage: var('n')
            n
            sage: chebyshev_T(n,1)
            1
            sage: chebyshev_T(n,0)
            1/2*(-1)^(1/2*n)*((-1)^n + 1)
            sage: chebyshev_T(n,-1)
            (-1)^n
            sage: chebyshev_T._eval_special_values_(3/2,x)
            Traceback (most recent call last):
            ...
            ValueError: no special value found
            sage: chebyshev_T._eval_special_values_(n, 0.1)
            Traceback (most recent call last):
            ...
            ValueError: no special value found
        """
        if x == 1:
            return x

        if x == -1:
            return x**n

        if x == 0:
            return (1+(-1)**n)*(-1)**(n/2)/2

        raise ValueError("no special value found")

    def _evalf_(self, n, x, **kwds):
        """
        Evaluates :class:`chebyshev_T` numerically with mpmath.

        EXAMPLES::

            sage: chebyshev_T._evalf_(10,3)
            2.26195370000000e7
            sage: chebyshev_T._evalf_(10,3,parent=RealField(75))
            2.261953700000000000000e7
            sage: chebyshev_T._evalf_(10,I)
            -3363.00000000000
            sage: chebyshev_T._evalf_(5,0.3)
            0.998880000000000
            sage: chebyshev_T(1/2, 0)
            0.707106781186548
            sage: chebyshev_T(1/2, 3/2)
            1.11803398874989
            sage: chebyshev_T._evalf_(1.5, Mod(8,9))
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate chebyshev_T with parent Ring of integers modulo 9

        This simply evaluates using :class:`RealField` or :class:`ComplexField`::

            sage: chebyshev_T(1234.5, RDF(2.1))
            5.48174256255782e735
            sage: chebyshev_T(1234.5, I)
            -1.21629397684152e472 - 1.21629397684152e472*I

        For large values of ``n``, mpmath fails (but the algebraic formula
        still works)::

            sage: chebyshev_T._evalf_(10^6, 0.1)
            Traceback (most recent call last):
            ...
            NoConvergence: Hypergeometric series converges too slowly. Try increasing maxterms.
            sage: chebyshev_T(10^6, 0.1)
            0.636384327171504
        """
        try:
            real_parent = kwds['parent']
        except KeyError:
            real_parent = parent(x)

            if not is_RealField(real_parent) and not is_ComplexField(real_parent):
                # parent is not a real or complex field: figure out a good parent
                if x in RR:
                    x = RR(x)
                    real_parent = RR
                elif x in CC:
                    x = CC(x)
                    real_parent = CC

        if not is_RealField(real_parent) and not is_ComplexField(real_parent):
            raise TypeError("cannot evaluate chebyshev_T with parent {}".format(real_parent))

        from sage.libs.mpmath.all import call as mpcall
        from sage.libs.mpmath.all import chebyt as mpchebyt

        return mpcall(mpchebyt, n, x, parent=real_parent)

    def _maxima_init_evaled_(self, n, x):
        """
        Evaluate the Chebyshev polynomial ``self`` with maxima.

        EXAMPLES::

            sage: var('n, x')
            (n, x)
            sage: chebyshev_T._maxima_init_evaled_(1,x)
            '_SAGE_VAR_x'
            sage: maxima(chebyshev_T(n, chebyshev_T(n, x)))
            chebyshev_t(_SAGE_VAR_n,chebyshev_t(_SAGE_VAR_n,_SAGE_VAR_x))
        """
        return maxima.eval('chebyshev_t({0},{1})'.format(n._maxima_init_(), x._maxima_init_()))

    def eval_formula(self, n, x):
        """
        Evaluate ``chebyshev_T`` using an explicit formula.
        See [ASHandbook]_ 227 (p. 782) for details for the recurions.
        See also [EffCheby]_ for fast evaluation techniques.

        INPUT:

        - ``n`` -- an integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: chebyshev_T.eval_formula(-1,x)
            x
            sage: chebyshev_T.eval_formula(0,x)
            1
            sage: chebyshev_T.eval_formula(1,x)
            x
            sage: chebyshev_T.eval_formula(2,0.1) == chebyshev_T._evalf_(2,0.1)
            True
            sage: chebyshev_T.eval_formula(10,x)
            512*x^10 - 1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1
            sage: chebyshev_T.eval_algebraic(10,x).expand()
            512*x^10 - 1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1
        """
        if n < 0:
            return self.eval_formula(-n, x)
        elif n == 0:
            return parent(x).one()

        res = parent(x).zero()
        for j in xrange(0, n//2+1):
            f = factorial(n-1-j) / factorial(j) / factorial(n-2*j)
            res += (-1)**j * (2*x)**(n-2*j) * f
        res *= n/2
        return res

    def eval_algebraic(self, n, x):
        """
        Evaluate :class:`chebyshev_T` as polynomial, using a recursive
        formula.

        INPUT:

        - ``n`` -- an integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: chebyshev_T.eval_algebraic(5, x)
            2*(2*(2*x^2 - 1)*x - x)*(2*x^2 - 1) - x
            sage: chebyshev_T(-7, x) - chebyshev_T(7,x)
            0
            sage: R.<t> = ZZ[]
            sage: chebyshev_T.eval_algebraic(-1, t)
            t
            sage: chebyshev_T.eval_algebraic(0, t)
            1
            sage: chebyshev_T.eval_algebraic(1, t)
            t
            sage: chebyshev_T(7^100, 1/2)
            1/2
            sage: chebyshev_T(7^100, Mod(2,3))
            2
            sage: n = 97; x = RIF(pi/2/n)
            sage: chebyshev_T(n, cos(x)).contains_zero()
            True
            sage: R.<t> = Zp(2, 8, 'capped-abs')[]
            sage: chebyshev_T(10^6+1, t)
            (2^7 + O(2^8))*t^5 + (O(2^8))*t^4 + (2^6 + O(2^8))*t^3 + (O(2^8))*t^2 + (1 + 2^6 + O(2^8))*t + (O(2^8))
        """
        if n == 0:
            return parent(x).one()
        if n < 0:
            return self._eval_recursive_(-n, x)[0]
        return self._eval_recursive_(n, x)[0]

    def _eval_recursive_(self, n, x, both=False):
        """
        If ``both=True``, compute ``(T(n,x), T(n-1,x))`` using a
        recursive formula.
        If ``both=False``, return instead a tuple ``(T(n,x), False)``.

        EXAMPLES::

            sage: chebyshev_T._eval_recursive_(5, x)
            (2*(2*(2*x^2 - 1)*x - x)*(2*x^2 - 1) - x, False)
            sage: chebyshev_T._eval_recursive_(5, x, True)
            (2*(2*(2*x^2 - 1)*x - x)*(2*x^2 - 1) - x, 2*(2*x^2 - 1)^2 - 1)
        """
        if n == 1:
            return x, parent(x).one()

        assert n >= 2
        a, b = self._eval_recursive_((n+1)//2, x, both or n % 2)
        if n % 2 == 0:
            return 2*a*a - 1, both and 2*a*b - x
        else:
            return 2*a*b - x, both and 2*b*b - 1


    def _eval_numpy_(self, n, x):
        """
        Evaluate ``self`` using numpy.

        EXAMPLES::

            sage: import numpy
            sage: z = numpy.array([1,2])
            sage: z2 = numpy.array([[1,2],[1,2]])
            sage: z3 = numpy.array([1,2,3.])
            sage: chebyshev_T(1,z)
            array([ 1.,  2.])
            sage: chebyshev_T(1,z2)
            array([[ 1.,  2.],
                   [ 1.,  2.]])
            sage: chebyshev_T(1,z3)
            array([ 1.,  2.,  3.])
            sage: chebyshev_T(z,0.1)
            array([ 0.1 , -0.98])
        """
        from scipy.special import eval_chebyt
        return eval_chebyt(n, x)

    def _derivative_(self, n, x, diff_param):
        """
        Return the derivative of :class:`chebyshev_T` in form of the Chebyshev
        polynomial of the second kind :class:`chebyshev_U`.

        EXAMPLES::

            sage: var('k')
            k
            sage: derivative(chebyshev_T(k,x),x)
            k*chebyshev_U(k - 1, x)
            sage: derivative(chebyshev_T(3,x),x)
            12*x^2 - 3
            sage: derivative(chebyshev_T(k,x),k)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative w.r.t. to the index is not supported yet
        """
        if diff_param == 0:
            raise NotImplementedError("derivative w.r.t. to the index is not supported yet")
        elif diff_param == 1:
            return n*chebyshev_U(n-1, x)
        raise ValueError("illegal differentiation parameter {}".format(diff_param))

chebyshev_T = Func_chebyshev_T()

class Func_chebyshev_U(ChebyshevPolynomial):
    """
    Class for the Chebyshev polynomial of the second kind.

    REFERENCE:

    - [ASHandbook]_ 22.8.3 page 783 and 6.1.22 page 256.

    EXAMPLES::

        sage: R.<t> = QQ[]
        sage: chebyshev_U(2,t)
        4*t^2 - 1
        sage: chebyshev_U(3,t)
        8*t^3 - 4*t
    """
    def __init__(self):
        """
        Init method for the chebyshev polynomials of the second kind.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import Func_chebyshev_U
            sage: chebyshev_U2 = Func_chebyshev_U()
            sage: chebyshev_U2(1,x)
            2*x
        """
        ChebyshevPolynomial.__init__(self, "chebyshev_U", nargs=2,
                                     conversions=dict(maxima='chebyshev_u',
                                                      mathematica='ChebyshevU'))

    def eval_formula(self, n, x):
        """
        Evaluate ``chebyshev_U`` using an explicit formula.
        See [ASHandbook]_ 227 (p. 782) for details on the recurions.
        See also [EffCheby]_ for the recursion formulas.

        INPUT:

        - ``n`` -- an integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: chebyshev_U.eval_formula(10, x)
            1024*x^10 - 2304*x^8 + 1792*x^6 - 560*x^4 + 60*x^2 - 1
            sage: chebyshev_U.eval_formula(-2, x)
            -1
            sage: chebyshev_U.eval_formula(-1, x)
            0
            sage: chebyshev_U.eval_formula(0, x)
            1
            sage: chebyshev_U.eval_formula(1, x)
            2*x
            sage: chebyshev_U.eval_formula(2,0.1) == chebyshev_U._evalf_(2,0.1)
            True
        """
        if n < -1:
            return -self.eval_formula(-n-2, x)

        res = parent(x).zero()
        for j in xrange(0, n//2+1):
            f = binomial(n-j, j)
            res += (-1)**j * (2*x)**(n-2*j) * f
        return res

    def eval_algebraic(self, n, x):
        """
        Evaluate :class:`chebyshev_U` as polynomial, using a recursive
        formula.

        INPUT:

        - ``n`` -- an integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: chebyshev_U.eval_algebraic(5,x)
            -2*((2*x + 1)*(2*x - 1)*x - 4*(2*x^2 - 1)*x)*(2*x + 1)*(2*x - 1)
            sage: parent(chebyshev_U(3, Mod(8,9)))
            Ring of integers modulo 9
            sage: parent(chebyshev_U(3, Mod(1,9)))
            Ring of integers modulo 9
            sage: chebyshev_U(-3,x) + chebyshev_U(1,x)
            0
            sage: chebyshev_U(-1,Mod(5,8))
            0
            sage: parent(chebyshev_U(-1,Mod(5,8)))
            Ring of integers modulo 8
            sage: R.<t> = ZZ[]
            sage: chebyshev_U.eval_algebraic(-2, t)
            -1
            sage: chebyshev_U.eval_algebraic(-1, t)
            0
            sage: chebyshev_U.eval_algebraic(0, t)
            1
            sage: chebyshev_U.eval_algebraic(1, t)
            2*t
            sage: n = 97; x = RIF(pi/n)
            sage: chebyshev_U(n-1, cos(x)).contains_zero()
            True
            sage: R.<t> = Zp(2, 6, 'capped-abs')[]
            sage: chebyshev_U(10^6+1, t)
            (2 + O(2^6))*t + (O(2^6))
        """
        if n == -1:
            return parent(x).zero()
        if n < 0:
            return -self._eval_recursive_(-n-2, x)[0]
        return self._eval_recursive_(n, x)[0]

    def _eval_recursive_(self, n, x, both=False):
        """
        If ``both=True``, compute ``(U(n,x), U(n-1,x))`` using a
        recursive formula.
        If ``both=False``, return instead a tuple ``(U(n,x), False)``.

        EXAMPLES::

            sage: chebyshev_U._eval_recursive_(3, x)
            (4*((2*x + 1)*(2*x - 1) - 2*x^2)*x, False)
            sage: chebyshev_U._eval_recursive_(3, x, True)
            (4*((2*x + 1)*(2*x - 1) - 2*x^2)*x, ((2*x + 1)*(2*x - 1) + 2*x)*((2*x + 1)*(2*x - 1) - 2*x))
        """
        if n == 0:
            return parent(x).one(), 2*x

        assert n >= 1
        a, b = self._eval_recursive_((n-1)//2, x, True)
        if n % 2 == 0:
            return (b+a)*(b-a), both and 2*b*(x*b-a)
        else:
            return 2*a*(b-x*a), both and (b+a)*(b-a)

    def _maxima_init_evaled_(self, n, x):
        """
        Uses maxima to evaluate ``self``.

        EXAMPLES::

            sage: var('n, x')
            (n, x)
            sage: maxima(chebyshev_U(5,x))
            32*_SAGE_VAR_x^5-32*_SAGE_VAR_x^3+6*_SAGE_VAR_x
            sage: maxima(chebyshev_U(n,x))
            chebyshev_u(_SAGE_VAR_n,_SAGE_VAR_x)
            sage: maxima(chebyshev_U(2,x))
            4*_SAGE_VAR_x^2-1
        """
        return maxima.eval('chebyshev_u({0},{1})'.format(n._maxima_init_(), x._maxima_init_()))

    def _evalf_(self, n, x, **kwds):
        """
        Evaluate :class:`chebyshev_U` numerically with mpmath.

        EXAMPLES::

            sage: chebyshev_U(5,-4+3.*I)
            98280.0000000000 - 11310.0000000000*I
            sage: chebyshev_U(10,3).n(75)
            4.661117900000000000000e7
            sage: chebyshev_U._evalf_(1.5, Mod(8,9))
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate chebyshev_U with parent Ring of integers modulo 9
        """
        try:
            real_parent = kwds['parent']
        except KeyError:
            real_parent = parent(x)

            if not is_RealField(real_parent) and not is_ComplexField(real_parent):
                # parent is not a real or complex field: figure out a good parent
                if x in RR:
                    x = RR(x)
                    real_parent = RR
                elif x in CC:
                    x = CC(x)
                    real_parent = CC

        if not is_RealField(real_parent) and not is_ComplexField(real_parent):
            raise TypeError("cannot evaluate chebyshev_U with parent {}".format(real_parent))

        from sage.libs.mpmath.all import call as mpcall
        from sage.libs.mpmath.all import chebyu as mpchebyu

        return mpcall(mpchebyu, n, x, parent=real_parent)

    def _eval_special_values_(self, n, x):
        """
        Values known for special values of x.
        See [ASHandbook]_ 22.4 (p.777).

        EXAMPLES::

            sage: var('n')
            n
            sage: chebyshev_U(n,1)
            n + 1
            sage: chebyshev_U(n,0)
            1/2*(-1)^(1/2*n)*((-1)^n + 1)
            sage: chebyshev_U(n,-1)
            (-1)^n*(n + 1)
            sage: chebyshev_U._eval_special_values_(n, 2)
            Traceback (most recent call last):
            ...
            ValueError: no special value found
        """
        if x == 1:
            return x*(n+1)

        if x == -1:
            return x**n*(n+1)

        if x == 0:
            return (1+(-1)**n)*(-1)**(n/2)/2

        raise ValueError("no special value found")

    def _eval_numpy_(self, n, x):
        """
        Evaluate ``self`` using numpy.

        EXAMPLES::

            sage: import numpy
            sage: z = numpy.array([1,2])
            sage: z2 = numpy.array([[1,2],[1,2]])
            sage: z3 = numpy.array([1,2,3.])
            sage: chebyshev_U(1,z)
            array([ 2.,  4.])
            sage: chebyshev_U(1,z2)
            array([[ 2.,  4.],
                   [ 2.,  4.]])
            sage: chebyshev_U(1,z3)
            array([ 2.,  4.,  6.])
            sage: chebyshev_U(z,0.1)
            array([ 0.2 , -0.96])
        """
        from scipy.special import eval_chebyu
        return eval_chebyu(n, x)

    def _derivative_(self, n, x, diff_param):
        """
        Return the derivative of :class:`chebyshev_U` in form of the Chebyshev
        polynomials of the first and second kind.

        EXAMPLES::

            sage: var('k')
            k
            sage: derivative(chebyshev_U(k,x),x)
            ((k + 1)*chebyshev_T(k + 1, x) - x*chebyshev_U(k, x))/(x^2 - 1)
            sage: derivative(chebyshev_U(3,x),x)
            24*x^2 - 4
            sage: derivative(chebyshev_U(k,x),k)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative w.r.t. to the index is not supported yet
        """
        if diff_param == 0:
            raise NotImplementedError("derivative w.r.t. to the index is not supported yet")
        elif diff_param == 1:
            return ((n+1)*chebyshev_T(n+1, x) - x*chebyshev_U(n,x)) / (x*x-1)
        raise ValueError("illegal differentiation parameter {}".format(diff_param))

chebyshev_U = Func_chebyshev_U()


class Func_legendre_P(OrthogonalPolynomial):
    def __init__(self):
        r"""
        Init method for the Legendre polynomials of the first kind.

        EXAMPLES::

            sage: loads(dumps(legendre_P))
            legendre_P
        """
        OrthogonalPolynomial.__init__(self, "legendre_P", nargs=2, latex_name=r"P",
                conversions={'maxima':'legendre_p', 'mathematica':'LegendreP',
                    'maple':'LegendreP'})

    def __call__(self, n, x, *args, **kwds):
        r"""
        Return an evaluation or call super.
        
        EXAMPLES::
            
            sage: legendre_P(1/2, I+1)
            1.05338240025858 + 0.359890322109665*I
            sage: legendre_P(1/2, I+1, hold=True)
            legendre_P(1/2, I + 1)
        """
        algorithm = kwds.get('algorithm', None)
        if algorithm == 'pari':
            return self.eval_pari(n, x)
        
        if (n in ZZ or SR(n).is_real()) and not kwds.get('hold', False):
            ret = self._eval_(n, x, *args, **kwds)
            if ret is not None:
                return ret

        return super(OrthogonalPolynomial,self).__call__(n, x, *args, **kwds)

    def _eval_(self, n, x, *args, **kwds):
        r"""
        Return an evaluation of this Legendre P expression.
        
        EXAMPLES::
            
            sage: legendre_P(4, 2.0, algorithm='pari')
            55.3750000000000
            sage: legendre_P(1, x)
            x
            sage: legendre_P(4, x+1)
            35/8*(x + 1)^4 - 15/4*(x + 1)^2 + 3/8
            sage: legendre_P(1/2, I+1)
            1.05338240025858 + 0.359890322109665*I
        """
        algorithm = kwds.get('algorithm', None)
        if algorithm == 'pari':
            return self.eval_pari(n, x)

        ret = self._eval_special_values_(n, x)
        if ret is not None:
            return ret
        if n in ZZ:
            ret = self.eval_pari(n, x)
            if ret is not None:
                return ret
        if SR(x).is_numeric() and SR(n).is_numeric():
            return self._evalf_(n, x, **kwds)

    def _eval_special_values_(self, n, x):
        """
        Special values known.
        
        EXAMPLES::
            
            sage: legendre_P(0, 0)
            1
            sage: legendre_P(1, x)
            x
        """
        if n == 0 or n == -1:
            return 1
        if n == 1 or n == -2:
            return x
    
    def _evalf_(self, n, x, **kwds):
        """
        EXAMPLES::
            
            sage: legendre_P(4, 2.)
            55.3750000000000
            sage: legendre_P(5.5,1.00001)
            1.00017875754114
            sage: legendre_P(1/2, I+1).n()
            1.05338240025858 + 0.359890322109665*I
            sage: legendre_P(1/2, I+1).n(59)
            1.0533824002585803 + 0.35989032210966537*I
        """
        ret = self._eval_special_values_(n, x)
        if ret is not None:
            return ret
        
        real_parent = kwds.get('parent', None)
        if real_parent is None:
            real_parent = parent(x)

            if not is_RealField(real_parent) and not is_ComplexField(real_parent):
                # parent is not a real or complex field: figure out a good parent
                if x in RR:
                    x = RR(x)
                    real_parent = RR
                elif x in CC:
                    x = CC(x)
                    real_parent = CC

        if not is_RealField(real_parent) and not is_ComplexField(real_parent):
            raise TypeError("cannot evaluate legendre_P with parent {}".format(real_parent))

        if real_parent.prec() <= 53:
            from scipy.special import eval_legendre
            if real_parent is RR:
                return RR(eval_legendre(float(n), float(x)))
            else:
                return real_parent(eval_legendre(float(n), complex(x)))
        else:
            import mpmath
            from sage.libs.mpmath.all import call as mpcall   
            return mpcall(mpmath.legenp, n, 0, x, parent=real_parent, prec=real_parent.prec())

    def eval_pari(self, n, arg):
        """
        Use Pari to evaluate legendre_P for integer, symbolic, and
        polynomial argument.
        
        EXAMPLES::
            
            sage: R.<x> = QQ[]
            sage: legendre_P(4,x)
            35/8*x^4 - 15/4*x^2 + 3/8
            sage: legendre_P(10000,x).coeffs()[1]
            0
            sage: var('t,x')
            (t, x)
            sage: legendre_P(-5,t)
            35/8*t^4 - 15/4*t^2 + 3/8
            sage: legendre_P(4, x+1)
            35/8*(x + 1)^4 - 15/4*(x + 1)^2 + 3/8
            sage: legendre_P(4, sqrt(2)).simplify()
            83/8
            sage: legendre_P(4, I*e)
            35/8*e^4 + 15/4*e^2 + 3/8
        """
        if n<0:
            n = - n - 1
        P = parent(arg)
        if P in (ZZ, QQ, RR, CC, SR):
            from sage.libs.pari.all import pari
            R = PolynomialRing(QQ, 'x')
            pol = R(pari.pollegendre(n))
            pol = sum([b*arg**a for (a,b) in enumerate(pol)])
            return pol
        elif is_PolynomialRing(P):
            from sage.libs.pari.all import pari
            if arg == P.gen():
                return P(pari.pollegendre(n))
            else:
                R = PolynomialRing(QQ, 'x')
                pol = R(pari.pollegendre(n))
                pol = pol.subs({pol.parent().gen():arg})
                pol = pol.change_ring(P.base_ring())
                return pol

    def _derivative_(self, n, x, *args,**kwds):
        """
        Return the derivative of legendre_P.
  
        EXAMPLES:: 
            sage: n = var('n')
            sage: derivative(legendre_P(n,x), x)
            (n*x*legendre_P(n, x) - n*legendre_P(n - 1, x))/(x^2 - 1) 
            sage: derivative(legendre_P(3,x), x) 
            15/2*x^2 - 3/2
            sage: derivative(legendre_P(n,x), n)
            Traceback (most recent call last):
            ...
            NotImplementedError: Derivative w.r.t. to the index is not supported.
        """
        diff_param = kwds['diff_param'] 
        if diff_param == 0:  
            raise NotImplementedError("Derivative w.r.t. to the index is not supported.")
        else:
            return (n*legendre_P(n-1, x) - n*x*legendre_P(n, x))/(1 - x**2)
         
legendre_P = Func_legendre_P()

class Func_legendre_Q(OrthogonalPolynomial):
    def __init__(self):
        r"""
        EXAMPLES::

            sage: loads(dumps(legendre_Q))
            legendre_Q
        """
        OrthogonalPolynomial.__init__(self, "legendre_Q", nargs=2, latex_name=r"Q",
                conversions={'maxima':'legendre_q', 'mathematica':'LegendreQ',
                    'maple':'LegendreQ'})

    def __call__(self, n, x, *args, **kwds):
        r"""
        Return an evaluation or call super.
        
        EXAMPLES::
            
            sage: legendre_Q(1/2, I+1)
            -0.511424110789061 + 1.34356195297194*I
            sage: legendre_Q(1/2, I+1, hold=True)
            legendre_Q(1/2, I + 1)
        """
        algorithm = kwds.get('algorithm', None)
        if algorithm == 'maxima':
            return self._maxima_init_evaled_(n, x)
        
        if (n in ZZ or SR(n).is_real()) and not kwds.get('hold', False):
            ret = self._eval_(n, x, *args, **kwds)
            if ret is not None:
                return ret

        return super(OrthogonalPolynomial,self).__call__(n, x, *args, **kwds)

    def _eval_(self, n, x, *args, **kwds):
        r"""
        Return an evaluation of this Legendre Q expression.
        
        EXAMPLES::
            
            sage: legendre_Q(2,x,algorithm='recursive')
            3/4*x^2*(log(x + 1) - log(-x + 1)) - 3/2*x - 1/4*log(x + 1) + 1/4*log(-x + 1)
            sage: legendre_Q(5,0)
            -8/15
            sage: legendre_Q(2,2*x)
            1/4*(12*x^2 - 1)*(log(2*x + 1) - log(-2*x + 1)) - 3*x
            sage: legendre_Q(1/2, I+1)
            -0.511424110789061 + 1.34356195297194*I

        NOTE::
        
            Maxima (``algorithm='maxima'``) will output the complex
            conjugate of the correct result, see :trac:`16813`.
        """
        algorithm = kwds.get('algorithm', None)
        if algorithm == 'maxima':
            return self._maxima_init_evaled_(n, x)
        if algorithm == 'recursive':
            return self.eval_recursive(n, x)

        ret = self._eval_special_values_(n, x)
        if ret is not None:
            return ret
        if n in ZZ and not SR(x).is_numeric():
            return self.eval_formula(n, x)
        if SR(x).is_numeric() and SR(n).is_numeric():
            return self._evalf_(n, x, **kwds)

    def _maxima_init_evaled_(self, n, x):
        """
        Return a string which represents this function evaluated at
        ``n, x`` in Maxima.

        NOTE::
        
            Maxima will output the complex conjugate of the
            correct result, see :trac:`16813`.

        EXAMPLES::
            
            sage: legendre_Q._maxima_init_evaled_(20,x).coeff(x^10)
            -29113619535/131072*log(-(x + 1)/(x - 1))
        """
        _init()
        return sage_eval(maxima.eval('legendre_q(%s,x)'%ZZ(n)), locals={'x':x}) 
     
    def _eval_special_values_(self, n, x):
        """
        Special values known.
        
        EXAMPLES::
            
            sage: var('n')
            n
            sage: legendre_Q(n,0)
            -1/2*sqrt(pi)*sin(1/2*pi*n)*gamma(1/2*n + 1/2)/gamma(1/2*n + 1)
        """
        if x == 1:
            from sage.symbolic.constants import NaN
            return NaN
        
        if x == -1:
            from sage.symbolic.constants import NaN
            return NaN
        
        if x == 0:
            from sage.functions.other import gamma, sqrt
            from sage.functions.trig import sin
            try:
                return -(sqrt(SR.pi()))/2*sin(SR.pi()/2*n)*\
                       gamma((n+1)/2)/gamma(n/2 + 1)
            except TypeError:
                pass

    def _evalf_(self, n, x, **kwds):
        """
        Float evaluation of Legendre Q(n, x) function.
        
        EXAMPLES::
            
            sage: legendre_Q(4, 2)
            0.00116107583162041 - 86.9828465962674*I
            sage: legendre_Q(1/2, I+1)
            -0.511424110789061 + 1.34356195297194*I
            sage: legendre_Q(1/2, I+1).n(59)
            -0.51142411078906080 + 1.3435619529719394*I
        """
        ret = self._eval_special_values_(n, x)
        if ret is not None:
            return ret
        
        real_parent = kwds.get('parent', None)
        if real_parent is None:
            real_parent = parent(x)

            if not is_RealField(real_parent) and not is_ComplexField(real_parent):
                # parent is not a real or complex field: figure out a good parent
                if x in RR:
                    x = RR(x)
                    real_parent = RR
                elif x in CC:
                    x = CC(x)
                    real_parent = CC

        if not is_RealField(real_parent) and not is_ComplexField(real_parent):
            raise TypeError("cannot evaluate legendre_Q with parent {}".format(real_parent))

        import mpmath
        from sage.libs.mpmath.all import call as mpcall   
        return mpcall(mpmath.legenq, n, 0, x, parent=real_parent, prec=real_parent.prec())

    def eval_recursive(self, n, arg, **kwds):
        """
        Return expanded Legendre Q(n, arg) function expression.
        
        EXAMPLES::
        
            sage: legendre_Q(2,x,algorithm='recursive')
            3/4*x^2*(log(x + 1) - log(-x + 1)) - 3/2*x - 1/4*log(x + 1) + 1/4*log(-x + 1)
            sage: legendre_Q.eval_recursive(20,x).expand().coeff(x^10)
            -29113619535/131072*log(x + 1) + 29113619535/131072*log(-x + 1)
        """
        from sage.functions.log import ln
        if n == 0: 
            return (ln(1+arg)-ln(1-arg))/2
        elif n == 1: 
            return arg/2*(ln(1+arg)-ln(1-arg))-1 
        
        R = PolynomialRing(QQ, 'x,l')
        (x,l) = R.gens()
        help1 = l/2
        help2 = x/2*l-1 
        for j in xrange(1,n): 
            help3 = (2*j+1)*x*help2 - j*help1
            help3 = help3/(j+1)
            help1 = help2
            help2 = help3 

        sum1 = sum([help3.monomial_coefficient(mon)*arg**(mon.exponents()[0][0])
                    for mon in help3.monomials() if not l.divides(mon)])
        sum2 = sum([help3.monomial_coefficient(mon)*arg**(mon.exponents()[0][0])*(ln(1+arg)-ln(1-arg))
                    for mon in help3.monomials() if l.divides(mon)])
        return sum1 + sum2
    
    def eval_formula(self, n, arg, **kwds):
        """
        Return expanded Legendre ``Q(n, arg)`` function expression.
        
        We use the fact that ``Q(n,x)`` is of form
        
        .. math::

            Q_n(x) = \frac{1}{2}P_n(x)(\log(1+x)-\log(1-x))+W_n(x),
        
        with ``W_n(x)`` a holonomic function.

        REFERENCE:
    
        - T. M. Dunster, Legendre and Related Functions, http://dlmf.nist.gov/14.7#E2
        
        EXAMPLES::
        
            sage: legendre_Q.eval_formula(1,x)
            1/2*x*(log(x + 1) - log(-x + 1)) - 1
            sage: legendre_Q.eval_formula(2,x).expand().collect(log(1+x)).collect(log(1-x))
            1/4*(3*x^2 - 1)*log(x + 1) - 1/4*(3*x^2 - 1)*log(-x + 1) - 3/2*x
            sage: legendre_Q.eval_formula(20,x).coeff(x^10)
            -29113619535/131072*log(x + 1) + 29113619535/131072*log(-x + 1)
        """
        from sage.functions.log import ln
        if n == 0: 
            return (ln(1+arg)-ln(1-arg))/2
        elif n == 1: 
            return arg/2*(ln(1+arg)-ln(1-arg))-1 
        
        arg = SR(arg)
        return legendre_P(n, arg)*(ln(1+arg)-ln(1-arg))/2 - self._Wfunc(n, arg)
    
    def _Wfunc(self, n, arg):
        """
        Helper function for ``eval_formula()``.
        
        EXAMPLES::
        
            sage: legendre_Q._Wfunc(2, x)
            3/2*x
            sage: legendre_Q._Wfunc(7, x)
            429/16*x^6 - 275/8*x^4 + 849/80*x^2 - 16/35
        """
        if n == 0:
            return 0
        if n == 1:
            return 1
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        help1 = 0
        help2 = 1 
        for j in xrange(2,n+1): 
            help3 = (2*j-1)*x*help2 - (j-1)*help1
            help3 = help3/j
            help1 = help2
            help2 = help3 

        return sum([b*arg**a for (a,b) in enumerate(help3)])

    def _derivative_(self, n, x, *args,**kwds):
        """
        Return the derivative of legendre_Q.
  
        EXAMPLES:: 
            sage: n = var('n')
            sage: derivative(legendre_Q(n,x), x)
            (n*x*legendre_Q(n, x) - n*legendre_Q(n - 1, x))/(x^2 - 1)
            sage: ex1 = legendre_Q(5,x,hold=True).diff(x).expand().simplify_full()
            sage: ex2 = legendre_Q(5,x).diff(x).expand().simplify_full()
            sage: ex1.subs(x=7).n() == ex2.subs(x=7).n()
            True
            sage: derivative(legendre_Q(n,x), n)
            Traceback (most recent call last):
            ...
            NotImplementedError: Derivative w.r.t. to the index is not supported.
        """
        diff_param = kwds['diff_param'] 
        if diff_param == 0:  
            raise NotImplementedError("Derivative w.r.t. to the index is not supported.")
        else: 
            return (n*x*legendre_Q(n, x) - n*legendre_Q(n-1, x))/(x**2 - 1)

legendre_Q = Func_legendre_Q()

class Func_assoc_legendre_P(OrthogonalPolynomial):
    def __init__(self):
        r"""
        EXAMPLES::

            sage: loads(dumps(gen_legendre_P))
            gen_legendre_P
        """
        OrthogonalPolynomial.__init__(self, "gen_legendre_P", nargs=3, latex_name=r"P",
                conversions={'maxima':'assoc_legendre_p', 'mathematica':'LegendreP',
                    'maple':'LegendreP'})

    def __call__(self, n, m, x, *args, **kwds):
        r"""
        Return an evaluation or call super.
        
        EXAMPLES::
        
            sage: gen_legendre_P(3,2,I+1)
            45.0000000000000 - 15.0000000000000*I
            sage: gen_legendre_P(3,2,I+1,hold=True)
            gen_legendre_P(3, 2, I + 1)
        """
        algorithm = kwds.get('algorithm', None)
        if algorithm == 'maxima':
            return self._maxima_init_evaled_(n, m, x)
        
        if (n in ZZ or SR(n).is_real()) and not kwds.get('hold', False):
            ret = self._eval_(n, m, x, *args, **kwds)
            if ret is not None:
                return ret

        return super(OrthogonalPolynomial,self).__call__(n, m, x, *args, **kwds)

    def _eval_(self, n, m, x, *args, **kwds):
        r"""
        Return an evaluation of this Legendre P(n, m, x) expression.
        
        EXAMPLES::
        
            sage: gen_legendre_P(3,2,2,algorithm='maxima')
            -90
            sage: gen_legendre_P(13/2,2,0)
            2*sqrt(2)*gamma(19/4)/(sqrt(pi)*gamma(13/4))
            sage: gen_legendre_P(3,2,x)
            -15*(x^2 - 1)*x
            sage: gen_legendre_P(3,2,2.)
            -90.0000000000000
        """
        algorithm = kwds.get('algorithm', None)
        if algorithm == 'maxima':
            return self._maxima_init_evaled_(n, m, x)

        ret = self._eval_special_values_(n, m, x)
        if ret is not None:
            return ret
        if (n in ZZ and m in ZZ
            and n >= 0 and m >= 0
            and (x in ZZ or not SR(x).is_numeric())):
            return self.eval_poly(n, m, x)
        if SR(x).is_numeric() and SR(n).is_numeric() and SR(m).is_numeric():
            return self._evalf_(n, m, x, **kwds)

    def _maxima_init_evaled_(self, n, m, x):
        """
        Return a string which represents this function evaluated at
        ``n, m, x`` in Maxima.

        EXAMPLES::
        
            sage: gen_legendre_P._maxima_init_evaled_(20,6,x).expand().coeff(x^10)
            2508866163428625/128
        """
        _init()
        return sage_eval(maxima.eval('assoc_legendre_p(%s,%s,x)'%(ZZ(n),ZZ(m))), locals={'x':x}) 
     
    def _eval_special_values_(self, n, m, x):
        """
        Special values known.
        
        EXAMPLES::
            sage: gen_legendre_P(2,3,4)
            0
            sage: gen_legendre_P(2,0,4)==legendre_P(2,4)
            True
            sage: gen_legendre_P(2,2,4)
            45
            sage: gen_legendre_P(2,2,x)
            3*x^2 - 3
            sage: gen_legendre_P(13/2,2,0)
            2*sqrt(2)*gamma(19/4)/(sqrt(pi)*gamma(13/4))
            sage: (m,n) = var('m,n')
            sage: gen_legendre_P(n,m,0)
            2^m*cos(1/2*pi*(m + n))*gamma(1/2*m + 1/2*n + 1/2)/(sqrt(pi)*gamma(-1/2*m + 1/2*n + 1))
            sage: gen_legendre_P(n,3,0)
            8*cos(1/2*pi*(n + 3))*gamma(1/2*n + 2)/(sqrt(pi)*gamma(1/2*n - 1/2))
            sage: gen_legendre_P(3,m,0)
            2^m*cos(1/2*pi*(m + 3))*gamma(1/2*m + 2)/(sqrt(pi)*gamma(-1/2*m + 5/2))
        """
        if m > n:
            return 0
        if m == 0:
            return legendre_P(n, x)
        if n == m:
            return factorial(2*m)/2**m/factorial(m) * (x**2-1)**(m/2)
        if x == 0:
            from sage.functions.other import gamma, sqrt
            from sage.functions.trig import cos
            if m in QQ and n in QQ:
                return 2**m/sqrt(SR.pi())*cos((n+m)/2*SR.pi())*(gamma(QQ(n+m+1)/2)/gamma(QQ(n-m)/2+1))
            elif isinstance(n, Expression) or isinstance(m, Expression):
                return 2**m/sqrt(SR.pi())*cos((n+m)/2*SR.pi())*(gamma((n+m+1)/2)/gamma((n-m)/2+1))

    def _evalf_(self, n, m, x, **kwds):
        """
        Float evaluation of Legendre P(n, m, x) function.
        
        EXAMPLES::
            
            sage: gen_legendre_P(10,2,3.)
            -7.19496360000000e8
            sage: gen_legendre_P(5/2,2,1+I)
            14.3165258449040 - 12.7850496155152*I
            sage: gen_legendre_P(5/2,2,ComplexField(70)(1+I))
            14.316525844904028532 - 12.785049615515157033*I
        """
        ret = self._eval_special_values_(n, m, x)
        if ret is not None:
            return ret
        
        real_parent = kwds.get('parent', None)
        if real_parent is None:
            real_parent = parent(x)

            if not is_RealField(real_parent) and not is_ComplexField(real_parent):
                # parent is not a real or complex field: figure out a good parent
                if x in RR:
                    x = RR(x)
                    real_parent = RR
                elif x in CC:
                    x = CC(x)
                    real_parent = CC

        if not is_RealField(real_parent) and not is_ComplexField(real_parent):
            raise TypeError("cannot evaluate legendre_Q with parent {}".format(real_parent))

        import mpmath
        from sage.libs.mpmath.all import call as mpcall   
        return mpcall(mpmath.legenp, n, m, x, parent=real_parent, prec=real_parent.prec())

    def eval_poly(self, n, m, arg, **kwds):
        """
        Return the associated Legendre P(n, m, arg) polynomial for integers `n > -1, m > -1`.
        
        EXAMPLES::
        
            sage: gen_legendre_P(7,4,x)
            3465/2*(13*x^3 - 3*x)*(x^2 - 1)^2
            sage: gen_legendre_P(3,1,sqrt(x))
            -3/2*(5*x - 1)*sqrt(-x + 1)

        REFERENCE:
    
        - T. M. Dunster, Legendre and Related Functions, http://dlmf.nist.gov/14.7#E10
        """
        from sage.functions.other import factorial        
        if n < 0 or m < 0:
            return
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        p = (1-x**2)**ZZ(n)
        for i in range(m+n):
            p = p.diff(x)
        ex1 = (1-arg**2)**(QQ(m)/2)/2**n/factorial(ZZ(n))
        ex2 = sum([b*arg**a for (a,b) in enumerate(p)])
        return (-1)**(m+n)*ex1*ex2

    def _derivative_(self, n, m, x, *args,**kwds):
        """
        Return the derivative of ``gen_legendre_P(n,m,x)``.
  
        EXAMPLES::
        
            sage: (m,n) = var('m,n')
            sage: derivative(gen_legendre_P(n,m,x), x)
            -((n + 1)*x*gen_legendre_P(n, m, x) + (m - n - 1)*gen_legendre_P(n + 1, m, x))/(x^2 - 1)
            sage: gen_legendre_P(3,2,x,hold=True).diff(x).expand().simplify_full()
            -45*x^2 + 15
            sage: derivative(gen_legendre_P(n,m,x), n)
            Traceback (most recent call last):
            ...
            NotImplementedError: Derivative w.r.t. to the index is not supported.
        """
        diff_param = kwds['diff_param'] 
        if diff_param == 0:  
            raise NotImplementedError("Derivative w.r.t. to the index is not supported.")
        else: 
            return ((n-m+1)*gen_legendre_P(n+1, m, x) - (n+1)*x*gen_legendre_P(n, m, x))/(x**2 - 1)

gen_legendre_P = Func_assoc_legendre_P()

class Func_assoc_legendre_Q(OrthogonalPolynomial):
    def __init__(self):
        r"""
        EXAMPLES::

            sage: loads(dumps(gen_legendre_Q))
            gen_legendre_Q
        """
        OrthogonalPolynomial.__init__(self, "gen_legendre_Q", nargs=3, latex_name=r"Q",
                conversions={'maxima':'assoc_legendre_q', 'mathematica':'LegendreQ',
                    'maple':'LegendreQ'})

    def __call__(self, n, m, x, *args, **kwds):
        r"""
        Return an evaluation or call super.
        
        EXAMPLES::
        
        """
        algorithm = kwds.get('algorithm', None)
        if algorithm == 'maxima':
            return self._maxima_init_evaled_(n, m, x)
        
        if (n in ZZ or SR(n).is_real()) and not kwds.get('hold', False):
            ret = self._eval_(n, m, x, *args, **kwds)
            if ret is not None:
                return ret

        return super(OrthogonalPolynomial,self).__call__(n, m, x, *args, **kwds)

    def _eval_(self, n, m, x, *args, **kwds):
        r"""
        Return an evaluation of this Legendre Q(n, m, x) expression.
        
        EXAMPLES::
        
            sage: gen_legendre_Q(2,1,3,algorithm='maxima')
            -3*sqrt(-2)*(3*I*pi + 3*log(2)) + 25/4*sqrt(-2)
        """
        algorithm = kwds.get('algorithm', None)
        if algorithm == 'maxima':
            return self._maxima_init_evaled_(n, m, x)

        ret = self._eval_special_values_(n, m, x)
        if ret is not None:
            return ret
        if (n in ZZ and m in ZZ
            and n >= 0 and m >= 0
            and (x in ZZ or not SR(x).is_numeric())):
            return self.eval_recursive(n, m, x)
        if SR(x).is_numeric() and SR(n).is_numeric() and SR(m).is_numeric():
            return self._evalf_(n, m, x, **kwds)

    def _maxima_init_evaled_(self, n, m, x):
        """
        Return a string which represents this function evaluated at
        ``n, m, x`` in Maxima.

        EXAMPLES::
        
            sage: gen_legendre_Q._maxima_init_evaled_(2,1,3)
            -3*sqrt(-2)*(3*I*pi + 3*log(2)) + 25/4*sqrt(-2)
        """
        if m <= n:
            _init()
            return sage_eval(maxima.eval('assoc_legendre_q(%s,%s,x)'%(ZZ(n),ZZ(m))), locals={'x':x})

    def _eval_special_values_(self, n, m, x):
        """
        Special values known.
        
        EXAMPLES::

            sage: n, m = var('n m') 
            sage: gen_legendre_Q(n,m,0)
            -sqrt(pi)*2^(m - 1)*sin(1/2*pi*(m + n))*gamma(1/2*m + 1/2*n + 1/2)/gamma(-1/2*m + 1/2*n + 1)
        """
        if m == 0:
            return legendre_Q(n, x)
        if x.is_zero():
            from sage.functions.other import gamma, sqrt
            from sage.functions.trig import sin
            if m in QQ and n in QQ:
                return -(sqrt(SR.pi()))*sin(SR.pi()/2*(m+n))*gamma(QQ(m+n+1)/2)/gamma(QQ(n-m)/2 + 1)*2**(m-1)
            elif isinstance(n, Expression) or isinstance(m, Expression):
                return -(sqrt(SR.pi()))*sin(SR.pi()/2*(m+n))*gamma((m+n+1)/2)/gamma((n-m)/2 + 1)*2**(m-1)

    def _evalf_(self, n, m, x, **kwds):
        """
        Float evaluation of Legendre Q(n, m, x) function.
        
        EXAMPLES::

            sage: gen_legendre_Q(2,1,3.)
            -39.9859464434253 + 0.0165114736149193*I
            sage: gen_legendre_Q(2,1,ComplexField(70)(3))
            -39.985946443425296223 + 0.016511473614919329585*I
        """
        ret = self._eval_special_values_(n, m, x)
        if ret is not None:
            return ret
        
        real_parent = kwds.get('parent', None)
        if real_parent is None:
            real_parent = parent(x)

            if not is_RealField(real_parent) and not is_ComplexField(real_parent):
                # parent is not a real or complex field: figure out a good parent
                if x in RR:
                    x = RR(x)
                    real_parent = RR
                elif x in CC:
                    x = CC(x)
                    real_parent = CC

        if not is_RealField(real_parent) and not is_ComplexField(real_parent):
            raise TypeError("cannot evaluate legendre_Q with parent {}".format(real_parent))

        import mpmath
        from sage.libs.mpmath.all import call as mpcall   
        return mpcall(mpmath.legenq, n, m, x, parent=real_parent, prec=real_parent.prec())

    def eval_recursive(self, n, m, x, **kwds):
        """
        Return the associated Legendre Q(n, m, arg) function for integers `n > -1, m > -1`.
        
        EXAMPLES::

            sage: gen_legendre_Q(3,4,x)
            48/(x^2 - 1)^2
            sage: gen_legendre_Q(4,5,x)
            -384/((x^2 - 1)^2*sqrt(-x^2 + 1))
            sage: gen_legendre_Q(0,1,x)
            -1/sqrt(-x^2 + 1)
            sage: gen_legendre_Q(0,2,x)
            -1/2*((x + 1)^2 - (x - 1)^2)/(x^2 - 1)
            sage: gen_legendre_Q(2,2,x).subs(x=2).expand()
            9/2*I*pi - 9/2*log(3) + 14/3
        """
        from sage.functions.all import sqrt
        if m == n + 1 or n == 0:
            if m.mod(2).is_zero():
                denom = (1 - x**2)**(m/2)
            else:
                denom = sqrt(1 - x**2)*(1 - x**2)**((m-1)/2)
            if m == n + 1:
                return (-1)**m*(m-1).factorial()*2**n/denom
            else:
                return (-1)**m*(m-1).factorial()*((x+1)**m - (x-1)**m)/(2*denom)
        else:
            return ((n-m+1)*x*gen_legendre_Q(n,m-1,x)-(n+m-1)*gen_legendre_Q(n-1,m-1,x))/sqrt(1-x**2)

    def _derivative_(self, n, m, x, *args,**kwds):
        """
        Return the derivative of ``gen_legendre_Q(n,m,x)``.
  
        EXAMPLES::
        
            sage: (m,n) = var('m,n')
            sage: derivative(gen_legendre_Q(n,m,x), x)
            -((n + 1)*x*gen_legendre_Q(n, m, x) + (m - n - 1)*gen_legendre_Q(n + 1, m, x))/(x^2 - 1)
            sage: ex1=gen_legendre_Q(3,2,x,hold=True).diff(x).expand().simplify_full()
            sage: ex2=gen_legendre_Q(3,2,x).diff(x).expand().simplify_full()
            sage: ex1.subs(x=5).n() == ex2.subs(x=5).n()
            True
            sage: derivative(gen_legendre_Q(n,m,x), n)
            Traceback (most recent call last):
            ...
            NotImplementedError: Derivative w.r.t. to the index is not supported.
        """
        diff_param = kwds['diff_param'] 
        if diff_param == 0:  
            raise NotImplementedError("Derivative w.r.t. to the index is not supported.")
        else: 
            return ((n-m+1)*gen_legendre_Q(n+1, m, x) - (n+1)*x*gen_legendre_Q(n, m, x))/(x**2 - 1)

gen_legendre_Q = Func_assoc_legendre_Q()


def gen_laguerre(n,a,x):
    """
    Returns the generalized Laguerre polynomial for integers `n > -1`.
    Typically, `a = 1/2` or `a = -1/2`.

    REFERENCES:

    - Table on page 789 in [ASHandbook]_.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: gen_laguerre(2,1,x)
        1/2*x^2 - 3*x + 3
        sage: gen_laguerre(2,1/2,x)
        1/2*x^2 - 5/2*x + 15/8
        sage: gen_laguerre(2,-1/2,x)
        1/2*x^2 - 3/2*x + 3/8
        sage: gen_laguerre(2,0,x)
        1/2*x^2 - 2*x + 1
        sage: gen_laguerre(3,0,x)
        -1/6*x^3 + 3/2*x^2 - 3*x + 1
    """
    _init()
    return sage_eval(maxima.eval('gen_laguerre(%s,%s,x)'%(ZZ(n),a)), locals={'x':x})

def hermite(n,x):
    """
    Returns the Hermite polynomial for integers `n > -1`.

    REFERENCE:

    - [ASHandbook]_ 22.5.40 and 22.5.41, page 779.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: hermite(2,x)
        4*x^2 - 2
        sage: hermite(3,x)
        8*x^3 - 12*x
        sage: hermite(3,2)
        40
        sage: S.<y> = PolynomialRing(RR)
        sage: hermite(3,y)
        8.00000000000000*y^3 - 12.0000000000000*y
        sage: R.<x,y> = QQ[]
        sage: hermite(3,y^2)
        8*y^6 - 12*y^2
        sage: w = var('w')
        sage: hermite(3,2*w)
        8*(8*w^2 - 3)*w
    """
    _init()
    return sage_eval(maxima.eval('hermite(%s,x)'%ZZ(n)), locals={'x':x})

def jacobi_P(n,a,b,x):
    r"""
    Returns the Jacobi polynomial `P_n^{(a,b)}(x)` for
    integers `n > -1` and a and b symbolic or `a > -1`
    and `b > -1`. The Jacobi polynomials are actually defined
    for all a and b. However, the Jacobi polynomial weight
    `(1-x)^a(1+x)^b` isn't integrable for `a \leq -1`
    or `b \leq -1`.

    REFERENCE:

    - Table on page 789 in [ASHandbook]_.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: jacobi_P(2,0,0,x)
        3/2*x^2 - 1/2
        sage: jacobi_P(2,1,2,1.2)        # random output of low order bits
        5.009999999999998
    """
    _init()
    return sage_eval(maxima.eval('jacobi_p(%s,%s,%s,x)'%(ZZ(n),a,b)), locals={'x':x})

def laguerre(n,x):
    """
    Return the Laguerre polynomial for integers `n > -1`.

    REFERENCE:

    - [ASHandbook]_ 22.5.16, page 778 and page 789.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: laguerre(2,x)
        1/2*x^2 - 2*x + 1
        sage: laguerre(3,x)
        -1/6*x^3 + 3/2*x^2 - 3*x + 1
        sage: laguerre(2,2)
        -1
    """
    _init()
    return sage_eval(maxima.eval('laguerre(%s,x)'%ZZ(n)), locals={'x':x})

def ultraspherical(n,a,x):
    """
    Returns the ultraspherical (or Gegenbauer) polynomial for integers
    `n > -1`.
 
    Computed using Maxima.
 
    REFERENCE:
 
    - [ASHandbook]_ 22.5.27
 
    EXAMPLES::
 
        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: ultraspherical(2,3/2,x)
        15/2*x^2 - 3/2
        sage: ultraspherical(2,1/2,x)
        3/2*x^2 - 1/2
        sage: ultraspherical(1,1,x)
        2*x
        sage: t = PolynomialRing(RationalField(),"t").gen()
        sage: gegenbauer(3,2,t)
        32*t^3 - 12*t
    """
    _init()
    return sage_eval(maxima.eval('ultraspherical(%s,%s,x)'%(ZZ(n),a)), locals={'x':x})
 
gegenbauer = ultraspherical
