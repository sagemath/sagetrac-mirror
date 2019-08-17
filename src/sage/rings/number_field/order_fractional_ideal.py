"""
Ideals of Number Field Orders
"""

from __future__ import absolute_import

from sage.misc.cachefunc import cached_method
from sage.rings.all import ZZ
from sage.rings.ideal import Ideal_generic
from sage.structure.richcmp import richcmp


class OrderFractionalIdeal(Ideal_generic):
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
        return "(%s)" % (", ".join(map(str, self.gens())))

    def __repr__(self):
        return "Fractional ideal %s of non-maximal order" % self._repr_short()

    def _richcmp_(self, other, op):
        """
        Compare this ideal with `other` by comparing their underlying free
        module.
        """
        if not isinstance(other, OrderFractionalIdeal):
            return NotImplemented
        if not self.order() == other.order():
            return NotImplemented
        return richcmp(self.free_module(), other.free_module(), op)

    @cached_method
    def is_integral(self):
        """
        Return whether this ideal is integral.

        EXAMPLES::

           sage: R.<x> = PolynomialRing(QQ)
           sage: K.<a> = NumberField(x^5-x+1)
           sage: K.ideal(a).is_integral()
           True
           sage: (K.ideal(1) / (3*a+1)).is_integral()
           False
        """
        return all(a in self._order for a in self.gens())

    def order(self):
        return self._order

    def number_field(self):
        return self.ring()

    @cached_method
    def basis(self):
        """
        Return a basis for self as a `ZZ`-module.
        """
        K = self.number_field()
        order_basis = self.order().basis()
        from itertools import product

        ZZ_gens = [x[0] * x[1] for x in product(order_basis, self.gens())]

        from sage.matrix.constructor import Matrix

        M = Matrix([x.vector() for x in ZZ_gens])
        d = M.denominator()
        dM = (d * M).change_ring(ZZ)
        H = dM.hermite_form(include_zero_rows=False)
        Hd = H / d
        return [K(x) for x in Hd.rows()]

    @cached_method
    def free_module(self):
        r"""
        Return the free `\ZZ`-module contained in the vector space
        associated to the ambient number field, that corresponds
        to this ideal.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(7)
            sage: I = K.factor(11)[0][0]; I
            Fractional ideal (-2*z^4 - 2*z^2 - 2*z + 1)
            sage: A = I.free_module()
            sage: A              # warning -- choice of basis can be somewhat random
            Free module of degree 6 and rank 6 over Integer Ring
            User basis matrix:
            [11  0  0  0  0  0]
            [ 0 11  0  0  0  0]
            [ 0  0 11  0  0  0]
            [10  4  5  1  0  0]
            [ 5  1  1  0  1  0]
            [ 5  6  2  1  1  1]

        However, the actual `\ZZ`-module is not at all random::

            sage: A.basis_matrix().change_ring(ZZ).echelon_form()
            [ 1  0  0  5  1  1]
            [ 0  1  0  1  1  7]
            [ 0  0  1  7  6 10]
            [ 0  0  0 11  0  0]
            [ 0  0  0  0 11  0]
            [ 0  0  0  0  0 11]

        The ideal doesn't have to be integral::

            sage: J = I^(-1)
            sage: B = J.free_module()
            sage: B.echelonized_basis_matrix()
            [ 1/11     0     0  7/11  1/11  1/11]
            [    0  1/11     0  1/11  1/11  5/11]
            [    0     0  1/11  5/11  4/11 10/11]
            [    0     0     0     1     0     0]
            [    0     0     0     0     1     0]
            [    0     0     0     0     0     1]

        This also works for relative extensions::

            sage: K.<a,b> = NumberField([x^2 + 1, x^2 + 2])
            sage: I = K.fractional_ideal(4)
            sage: I.free_module()
            Free module of degree 4 and rank 4 over Integer Ring
            User basis matrix:
            [  4   0   0   0]
            [ -3   7  -1   1]
            [  3   7   1   1]
            [  0 -10   0  -2]
            sage: J = I^(-1); J.free_module()
            Free module of degree 4 and rank 4 over Integer Ring
            User basis matrix:
            [  1/4     0     0     0]
            [-3/16  7/16 -1/16  1/16]
            [ 3/16  7/16  1/16  1/16]
            [    0  -5/8     0  -1/8]

        An example of intersecting ideals by intersecting free modules.::

            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: I = K.factor(2)
            sage: p1 = I[0][0]; p2 = I[1][0]
            sage: N = p1.free_module().intersection(p2.free_module()); N
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [  1 1/2 1/2]
            [  0   1   1]
            [  0   0   2]
            sage: N.index_in(p1.free_module()).abs()
            2

        TESTS:

        Sage can find the free module associated to quite large ideals
        quickly (see :trac:`4627`)::

            sage: y = polygen(ZZ)
            sage: M.<a> = NumberField(y^20 - 2*y^19 + 10*y^17 - 15*y^16 + 40*y^14 - 64*y^13 + 46*y^12 + 8*y^11 - 32*y^10 + 8*y^9 + 46*y^8 - 64*y^7 + 40*y^6 - 15*y^4 + 10*y^3 - 2*y + 1)
            sage: M.ideal(prod(prime_range(6000, 6200))).free_module()
            Free module of degree 20 and rank 20 over Integer Ring
            User basis matrix:
            20 x 20 dense matrix over Rational Field
        """
        return basis_to_module(self.basis(), self.number_field())

    # def prime_factors(self):
    #     O = self.order()
    #     OK = O.integral_closure()

    #     aOK = OK.ideal(self.gens())
    #     above = aOK.prime_factors()

    #     return [O.intersection(p) for p in above]


def basis_to_module(B, K):
    r"""
    Given a basis `B` of elements for a `\ZZ`-submodule of a number
    field `K`, return the corresponding `\ZZ`-submodule.

    EXAMPLES::

        sage: K.<w> = NumberField(x^4 + 1)
        sage: from sage.rings.number_field.number_field_ideal import basis_to_module
        sage: basis_to_module([K.0, K.0^2 + 3], K)
        Free module of degree 4 and rank 2 over Integer Ring
        User basis matrix:
        [0 1 0 0]
        [3 0 1 0]
    """
    V, from_V, to_V = K.absolute_vector_space()
    M = ZZ ** (V.dimension())
    C = [to_V(K(b)) for b in B]
    return M.span_of_basis(C)
