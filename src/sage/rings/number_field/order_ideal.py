"""
Ideals of Number Field Orders
"""

from __future__ import absolute_import
from sage.misc.cachefunc import cached_method
from sage.rings.ideal import Ideal_generic
from sage.rings.all import ZZ
# from sage.rings.number_field.number_field_ideal import basis_to_module


class OrderIdeal(Ideal_generic):
    """
    An ideal of a number field order.
    """

    def __init__(self, order, gens, coerce=True):
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        self._order = order
        field = order.number_field()
        Ideal_generic.__init__(self, field, gens, coerce)

    # This should be changed
    def _repr_short(self):
        return '(%s)' % (', '.join(map(str, self.gens())))

    def __repr__(self):
        return "Fractional ideal %s of non-maximal order" % self._repr_short(
        )

    def is_integral(self):
        """
        Return True if this ideal is integral.

        EXAMPLES::

           sage: R.<x> = PolynomialRing(QQ)
           sage: K.<a> = NumberField(x^5-x+1)
           sage: K.ideal(a).is_integral()
           True
           sage: (K.ideal(1) / (3*a+1)).is_integral()
           False
        """
        try:
            return self.__is_integral
        except AttributeError:
            self.__is_integral = all(a in self._order for a in self.gens())
            return self.__is_integral

    # def order(self):
    #     return self._older

    # def number_field(self):
    #     return self.ring()

    # @cached_method
    # def basis(self):
    #     """
    #     Return a basis for self as a ZZ-module
    #     """
    #     K = self.number_field()
    #     order_basis = self.ring().basis()
    #     from itertools import product
    #     ZZ_gens = [x[0] * x[1] for x in product(order_basis, self.gens())]

    #     from sage.matrix.constructor import Matrix
    #     M = Matrix([x.vector() for x in ZZ_gens])
    #     d = M.denominator()
    #     dM = (d * M).change_ring(ZZ)
    #     H = dM.hermite_form(include_zero_rows=False)
    #     Hd = H / d
    #     return [K(x) for x in Hd.rows()]

    # @cached_method
    # def free_module(self):
    #     return basis_to_module(self.basis(), self.number_field())

    # def prime_factors(self):
    #     O = self.order()
    #     OK = O.integral_closure()

    #     aOK = OK.ideal(self.gens())
    #     above = aOK.prime_factors()

    #     return [O.intersection(p) for p in above]
