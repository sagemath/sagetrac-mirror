"""
Miscellaneous p-adic functions

p-adic functions from ell_rational_field.py, moved here to reduce
crowding in that file.
"""

######################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
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
######################################################################


import sage.rings.all as rings
import .padic_lseries as plseries
import math
sqrt = math.sqrt


def __check_padic_hypotheses(self, p):
    p = rings.Integer(p)
    if not p.is_prime():
        raise ValueError("p = (%s) must be prime" % p)
    if p == 2:
        raise ValueError("p must be odd")
    if self.conductor() % p == 0 or self.ap(p) % p == 0:
        raise ArithmeticError("p must be a good ordinary prime")
    return p


def padic_lseries(self, p, normalize='L_ratio'):
    r"""
    Return the `p`-adic `L`-series of self at `p`.

    This is an object whose :meth:`approx` method computes
    approximation to the true `p`-adic `L`-series to
    any desired precision.

    INPUT:

    -  ``p`` -- prime

    -  ``normalize`` --  'L_ratio' (default), 'period' or 'none';
       this is describes the way the modular symbols
       are normalized. See modular_symbol for
       more details.

    EXAMPLES::

        sage: J = J0(37)[0]
        sage: L = J.padic_lseries(5); L
        5-adic L-series of Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        sage: type(L)
        <class 'sage.schemes.elliptic_curves.padic_lseries.pAdicLseriesOrdinary'>

    We compute the `3`-adic `L`-series of two curves of
    rank `0` and in each case verify the interpolation property
    for their leading coefficient (i.e., value at 0)::

        sage: e = EllipticCurve('11a')
        sage: ms = e.modular_symbol()
        sage: [ms(1/11), ms(1/3), ms(0), ms(oo)]
        [0, -3/10, 1/5, 0]
        sage: ms(0)
        1/5
        sage: L = e.padic_lseries(3)
        sage: P = L.series(5)
        sage: P(0)
        2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 3^6 + O(3^7)
        sage: alpha = L.alpha(9); alpha
        2 + 3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^8 + O(3^9)
        sage: R.<x> = QQ[]
        sage: f = x^2 - e.ap(3)*x + 3
        sage: f(alpha)
        O(3^9)
        sage: r = e.lseries().L_ratio(); r
        1/5
        sage: (1 - alpha^(-1))^2 * r
        2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 3^6 + 3^7 + O(3^9)
        sage: P(0)
        2 + 3 + 3^2 + 2*3^3 + 2*3^5 + 3^6 + O(3^7)

    Next consider the curve 37b::

        sage: e = EllipticCurve('37b')
        sage: L = e.padic_lseries(3)
        sage: P = L.series(5)
        sage: alpha = L.alpha(9); alpha
        1 + 2*3 + 3^2 + 2*3^5 + 2*3^7 + 3^8 + O(3^9)
        sage: r = e.lseries().L_ratio(); r
        1/3
        sage: (1 - alpha^(-1))^2 * r
        3 + 3^2 + 2*3^4 + 2*3^5 + 2*3^6 + 3^7 + O(3^9)
        sage: P(0)
        3 + 3^2 + 2*3^4 + 2*3^5 + O(3^6)

    We can use eclib to compute the `L`-series::

        sage: e = EllipticCurve('11a')
        sage: L = e.padic_lseries(3,use_eclib=True)
        sage: L.series(5,prec=10)
        1 + 2*3^3 + 3^6 + O(3^7) + (2 + 2*3 + 3^2 + O(3^4))*T + (2 + 3 + 3^2 + 2*3^3 + O(3^4))*T^2 + (2*3 + 3^2 + O(3^3))*T^3 + (3 + 2*3^3 + O(3^4))*T^4 + (1 + 2*3 + 2*3^2 + O(3^4))*T^5 + (2 + 2*3^2 + O(3^3))*T^6 + (1 + 3 + 2*3^2 + 3^3 + O(3^4))*T^7 + (1 + 2*3 + 3^2 + 2*3^3 + O(3^4))*T^8 + (1 + 3 + O(3^2))*T^9 + O(T^10)

    """
    key = (p, normalize)
    try:
        return self._padic_lseries[key]
    except AttributeError:
        self._padic_lseries = {}
    except KeyError:
        pass

#    if self.ap(p) % p != 0:
    Lp = plseries.pAdicLseriesOrdinary(self, p, normalize=normalize)
#    else:
#        Lp = plseries.pAdicLseriesSupersingular(self, p,
#                              normalize = normalize, use_eclib=use_eclib)
    self._padic_lseries[key] = Lp
    return Lp
