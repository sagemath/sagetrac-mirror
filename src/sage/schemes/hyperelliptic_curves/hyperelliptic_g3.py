r"""
Hyperelliptic curves of genus 3 over a general ring

AUTHORS:

- Sorina Ionica, Elisa Lorenzo Garcia, Anna Somoza (2017-01-11): Initial version
"""
#*****************************************************************************
#  Copyright (C) 2017 Sorina Ionica <sorina.ionica@gmail.com>
#                2017 Elisa Lorenzo Garcia <elisa.lorenzo@gmail.com>
#                2017 Anna Somoza <anna.somoza@upc.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import hyperelliptic_generic
from . import invariants


class HyperellipticCurve_g3(hyperelliptic_generic.HyperellipticCurve_generic):
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
            [8/35, -144/8575, 32/7203, -128/84035, 128/1058841, 10673116768/12353145,
            149423644992/1008840175, 85384936192/1815912315, 398463049216/49433168575]
            sage: C = HyperellipticCurve(x^8 + 3*x^2 + 2*x +1)
            sage: C.shioda_invariants()
            [32, 648/49, 512/3, 6912/49, -8192/9, -24064/49, -89293312/36015,
            69632/21, 1441134592/108045]
            sage: C = HyperellipticCurve(x*(x^6 + 1))
            sage: C.shioda_invariants()
            [-4, 0, 1/6, 0, 1/36, 0, 1/70, 0, 1/420]
            sage: C = HyperellipticCurve(x**8+4*x**4+5, x)
            sage: C.shioda_invariants()
            [5728/35, 1901607/17150, 29738480/7203, 95245904/16807, -114656031488/1058841,
            -1832866169368/12353145, -4626584012971696/3026520525, 1415110141852160/363182463,
            713601762721839616/17795940687]

        TESTS:

        Check that the inheritance is correct::

            sage: R.<x> = QQ[]
            sage: f = x^7 - x^4 + 3
            sage: C = HyperellipticCurve(f)
            sage: type(C).mro()
            [<class 'sage.schemes.hyperelliptic_curves.constructor.HyperellipticCurve_g3_RationalField_with_category'>,
             <class 'sage.schemes.hyperelliptic_curves.constructor.HyperellipticCurve_g3_RationalField'>,
             <class 'sage.schemes.hyperelliptic_curves.hyperelliptic_g3.HyperellipticCurve_g3'>,
             <class 'sage.schemes.hyperelliptic_curves.hyperelliptic_rational_field.HyperellipticCurve_rational_field'>,
             <class 'sage.schemes.hyperelliptic_curves.hyperelliptic_generic.HyperellipticCurve_generic'>,
             ...]
        """
        f, h = self.hyperelliptic_polynomials()
        return invariants.shioda_invariants(4*f + h**2)
