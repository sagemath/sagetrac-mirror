from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_field
from sage.rings.polynomial.polynomial_quotient_ring_element import PolynomialQuotientRingElement
from sage.rings.real_mpfr import RealField
from sage.rings.padics.padic_base_leaves import pAdicFieldCappedRelative
from sage.rings.padics.padic_valuation import pAdicValuation_padic

class QpAlgebraElement(PolynomialQuotientRingElement):
    def __init__(self, parent, polynomial, check=True):
        PolynomialQuotientRingElement.__init__(self, parent, polynomial, check)
        self.abs_prec = parent.abs_prec

    def norm(self):
        """
        Returns the norm of the given element down to Qp.
        """
        K = self.parent().base_ring()
        return K.coerce(self.charpoly('y').subs(y=0))

    def abs(self,prec=None):
        """
        Returns the absolute value of the given element.
        Our normalization is that the absolute value extends
        the p-adic absolute value.

        If ``prec`` is not specified, returns an exact value,
        otherwise returns an element of RealField(prec).

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Qp(3))
            sage: f = x^3 - 3
            sage: from sage.rings.padics.QpAlgebra import QpAlgebra
            sage: A.<a> = QpAlgebra(f)
            sage: a.abs()
            0.693361274350635
            sage: a.abs(prec=10)
            0.69

        """
        r=self.norm().abs()**(1/self.parent().degree())
        if prec:
            return RealField(prec=prec)(r)
        elif self.abs_prec:
            return RealField(prec=self.abs_prec)(r)
        else:
            return r

class QpAlgebra(PolynomialQuotientRing_field):
    """
    This class is for Qp-algebras constructed as a quotient ring
    of Qp[X] by a polynomial ``f`` over Qp. ``name`` is the variable
    name of the canonical generator. ``abs_prec`` is the precision
    to use for the absolute value.

    EXAMPLES::

    First we try an Eisenstein extension::

        sage: R.<x> = PolynomialRing(Qp(3))
        sage: f = x^3 - 3
        sage: from sage.rings.padics.QpAlgebra import QpAlgebra
        sage: A.<a> = QpAlgebra(f)
        sage: norm(a)
        2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + 2*3^12 + 2*3^13 + 2*3^14 + 2*3^15 + 2*3^16 + 2*3^17 + 2*3^18 + 2*3^19 + 2*3^20 + O(3^21)
        sage: a.abs()
        0.693361274350635

    Now we try an extension which is neither totally ramified nor unramified::
        sage: R.<x> = PolynomialRing(Qp(3))
        sage: f = x^3 - 20
        sage: from sage.rings.padics.QpAlgebra import QpAlgebra
        sage: A.<a> = QpAlgebra(f)
        sage: a.abs()
        1.00000000000000
    """

    Element = QpAlgebraElement
    def __init__(self, f, names=None, abs_prec=53, category = None):

        K = f.base_ring()
        if not isinstance(
                K, pAdicFieldCappedRelative):
            raise TypeError("Defining polynomial must be over a p-adic field Qp.")
        if names is None:
            raise ValueError("Must specify name of a generator.")
        if not f.is_irreducible():
            raise NotImplementedError("QpAlgebra not implemented when f is not irreducible.")

        self.abs_prec = abs_prec
        self.defining_polynomial = f
        super(PolynomialQuotientRing_field,self).__init__(f.parent(), f, name=names, category=category)


    # Implement coercion from the base and from fraction fields
    # over a ring that coerces into the base
    def _coerce_map_from_(self,S):
        if self.base().coerce_map_from(S):
            return True
        elif isinstance(S,QpAlgebra) and self:
            return True
        return False

    def _repr_(self):
        return "Quotient of Qp(%s)[X] by f(X)=%s." % (self.base_ring().residue_characteristic(), self.defining_polynomial())

    def _latex_(self):
        return r"\text{}"+latex(self.base_ring())

    def defining_polynomial(self):
        return self.defining_polynomial

    def valuation(self):
        return pAdicValuation_padic(parent=self)

    def domain(self):
        return self.base_ring()

