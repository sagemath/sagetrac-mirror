include "sage/libs/linkages/padics/Polynomial_ram.pxi"
include "CA_template.pxi"

cdef class RelativeRamifiedCappedAbsoluteElement(CAElement):
    def _poly_rep(self):
        """
        Return the underlying polynomial representation of this element 
        (which is used for computations).

        For debug purpose.

        EXAMPLES::

            sage: K.<a> = ZqCA(125)
            sage: S.<x> = PolynomialRing(K)
            sage: W.<w> = K.extension(x^3 - 25*x^2 - 5*a*x + 5)
            sage: w._poly_rep()
            x
            sage: P = W(5)._poly_rep(); P
            5

        The coefficients of P are floating point p-adics::

            sage: ring = P.parent().base_ring() 
            sage: ring
            5-adic Unramified Extension Ring in a defined by x^3 + 3*x + 3
            sage: ring._prec_type()
            'floating-point'
        """
        return self.value
