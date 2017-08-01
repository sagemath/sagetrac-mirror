r"""
Hyperelliptic curves of genus 3 over a general ring

AUTHORS:

- Soria Ionica, Elisa Lorenzo Garcia, Anna Somoza (2016-01-11): Initial version
"""
from __future__ import absolute_import
#*****************************************************************************
#  Copyright (C) 2016 Sorina Ionica <sorina.ionica@gmail.com>
#                2016 Elisa Lorenzo Garcia <elisa.lorenzo@gmail.com>
#                2016 Anna Somoza <anna.somoza@upc.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import hyperelliptic_generic
from . import invariants


class HyperellipticCurve_g3_generic(hyperelliptic_generic.HyperellipticCurve_generic):
    def shioda_invariants(self):
        r"""
        Return the Shioda invariants `(J2, J3, J4, J5, J6, J7, J8, J9, J10)` of Shioda, [Sh1967]_.

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.invariants`

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^7 - x^4 + 3
            sage: C = HyperellipticCurve(f)
            sage: C.shioda_invariants()
            [1/70,  -9/34300,  1/57624,  -1/672280,  1/33882912,  333534899/6324810240,
            2334744453/1033052339200,  333534907/1859494210560,  778248143/101239129241600]
            sage: C = HyperellipticCurve(x^8 + 3*x^2 + 2*x +1)
            sage: C.shioda_invariants()
            [2,  81/392,  2/3,  27/196,  -2/9,  -47/1568,  -174401/4609920,  17/1344,
            703679/55319040]
            sage: C = HyperellipticCurve(x*(x^6 + 1))
            sage: C.shioda_invariants()
            [-1/4, 0, 1/1536, 0, 1/147456, 0, 1/4587520, 0, 1/440401920]
            sage: C = HyperellipticCurve(x**8+4*x**4+5, x)
            sage: C.shioda_invariants()
            [358/35, 1901607/1097600, 1858655/115248, 5952869/1075648, -447875123/16941456,
            -229108271171/25299240960, -289161500810731/12396628070400, 
            2763886995805/185949421056, 1393753442816093/36446086526976]

        """
        f, h = self.hyperelliptic_polynomials()
        #assert h == 0, 'Argument must be a simplified model of genus 3.'
        if h != 0:
            f = f + h**2/4
        return invariants.shioda_invariants(f)
