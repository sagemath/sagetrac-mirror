"""
Unramified Extension Generic

This file implements the shared functionality for unramified extensions.

AUTHORS:

- David Roe
"""
from __future__ import absolute_import

#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from .padic_extension_generic import pAdicExtensionGeneric
from .misc import precprint
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.misc.cachefunc import cached_method

class UnramifiedExtensionGeneric(pAdicExtensionGeneric):
    """
    An unramified extension of Qp or Zp.
    """
    def __init__(self, poly, prec, print_mode, names, element_class):
        """
        Initializes self

        INPUT:

            - poly -- Polynomial defining this extension.
            - prec -- The precision cap
            - print_mode -- a dictionary with print options
            - names -- a 4-tuple, (variable_name, residue_name,
              unramified_subextension_variable_name, uniformizer_name)
            - element_class -- the class for elements of this unramified extension.

        EXAMPLES::

            sage: R.<a> = Zq(27) #indirect doctest
        """
        #base = poly.base_ring()
        #if base.is_field():
        #    self._PQR = pqr.PolynomialQuotientRing_field(poly.parent(), poly, name = names)
        #else:
        #    self._PQR = pqr.PolynomialQuotientRing_domain(poly.parent(), poly, name = names)
        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)
        self._res_field = GF(self.prime_pow.pow_Integer_Integer(poly.degree()), name = names[1], modulus = poly.change_ring(poly.base_ring().residue_field()))

    def _repr_(self, do_latex = False):
        r"""
        Representation.

        EXAMPLES::

            sage: R.<a> = Zq(125); R #indirect doctest
            Unramified Extension in a defined by x^3 + 3*x + 3 with capped relative precision 20 over 5-adic Ring
            sage: latex(R) #indirect doctest
            \mathbf{Z}_{5^{3}}
        """
        if do_latex:
            if self.is_field():
                return "\\mathbf{Q}_{%s^{%s}}" % (self.prime(), self.degree())
            else:
                return "\\mathbf{Z}_{%s^{%s}}" % (self.prime(), self.degree())
        return "Unramified Extension in %s defined by %s %s over %s-adic %s"%(self.variable_name(), self.defining_polynomial(exact=True), precprint(self._prec_type(), self.precision_cap(), self.prime()), self.prime(), "Field" if self.is_field() else "Ring")

    def ramification_index(self, K = None):
        """
        Returns the ramification index of self over the subring K.

        INPUT:

            - K -- a subring (or subfield) of self.  Defaults to the
              base.

        EXAMPLES::

            sage: R.<a> = Zq(125); R.ramification_index()
            1
        """
        if K is None:
            return 1
        elif K is self:
            return 1
        else:
            raise NotImplementedError

    def inertia_degree(self, K = None):
        """
        Returns the inertia degree of self over the subring K.

        INPUT:

            - K -- a subring (or subfield) of self.  Defaults to the
              base.

        EXAMPLES::

            sage: R.<a> = Zq(125); R.inertia_degree()
            3
        """
        if K is None:
            return self.modulus().degree()
        elif K is self:
            return 1
        else:
            raise NotImplementedError

    #def extension(self, *args, **kwds):
    #    raise NotImplementedError

    #def get_extension(self):
    #    raise NotImplementedError

    def residue_class_field(self):
        """
        Returns the residue class field.

        EXAMPLES::

            sage: R.<a> = Zq(125); R.residue_class_field()
            Finite Field in a0 of size 5^3
        """
        #should eventually take advantage of finite field
        #\code{extension} or finite field
        #\code{unramified_extension_of_degree} over the automatic
        #coercion base.
        return self._res_field

    def residue_ring(self, n):
        """
        Return the quotient of the ring of integers by the `n`th power of its maximal ideal.

        EXAMPLES::

            sage: R.<a> = Zq(125)
            sage: R.residue_ring(1)
            Finite Field in a0 of size 5^3

        The following requires implementing more general Artinian rings::

            sage: R.residue_ring(2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if n == 1:
            return self._res_field
        else:
            raise NotImplementedError

    def discriminant(self, K=None):
        """
        Returns the discriminant of self over the subring K.

        INPUT:

            - K -- a subring/subfield (defaults to the base ring).

        EXAMPLES::

            sage: R.<a> = Zq(125)
            sage: R.discriminant()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if K is self:
            return 1
        else:
            raise NotImplementedError

    #def automorphisms(self):
    #    raise NotImplementedError

    #def galois_group(self):
    #    r"""
    #    Returns the Galois group of self's fraction field over Qp.
    #    """
    #    ##
    #    ## If K is a number field, then K.galois_group() can return
    #    ## other variants, i.e. via Pari or KASH. We could consider
    #    ## doing this.
    #    ##
    #    from sage.groups.perm_gps.permgroup import CyclicPermutationGroup
    #    return CyclicPermutationGroup(self.modulus().degree())

    #def is_abelian(self):
    #    return True

    def is_galois(self, K=None):
        """
        Returns True if this extension is Galois.

        Every unramified extension is Galois.

        INPUT:

            - K -- a subring/subfield (defaults to the base ring).

        EXAMPLES::

            sage: R.<a> = Zq(125); R.is_galois()
            True
        """
        if K is None or K is self:
            return True
        raise NotImplementedError

    def gen(self, n=0):
        """
        Returns a generator for this unramified extension.

        This is an element that satisfies the polynomial defining this
        extension.  Such an element will reduce to a generator of the
        corresponding residue field extension.

        EXAMPLES::

            sage: R.<a> = Zq(125); R.gen()
            a + O(5^20)
        """
        if n != 0:
            raise IndexError("only one generator")
        return self([0,1])

    @cached_method
    def _frob_gen(self, arithmetic = True):
        """
        Returns frobenius of the generator for this unramified extension

        EXAMPLES::
            sage: R.<a> = Zq(9)
            sage: R._frob_gen()
            (2*a + 1) + (2*a + 2)*3 + (2*a + 2)*3^2 + (2*a + 2)*3^3 + (2*a + 2)*3^4 + (2*a + 2)*3^5 + (2*a + 2)*3^6 + (2*a + 2)*3^7 + (2*a + 2)*3^8 + (2*a + 2)*3^9 + (2*a + 2)*3^10 + (2*a + 2)*3^11 + (2*a + 2)*3^12 + (2*a + 2)*3^13 + (2*a + 2)*3^14 + (2*a + 2)*3^15 + (2*a + 2)*3^16 + (2*a + 2)*3^17 + (2*a + 2)*3^18 + (2*a + 2)*3^19 + O(3^20)
        """
        p = self.prime()
        exp = p
        a = self.gen()
        if not arithmetic:
            exp = p**(self.degree()-1)
        approx = (self(a.residue()**exp)).lift_to_precision(self.precision_cap()) #first approximation
        f = self.defining_polynomial()
        g = f.derivative()
        while(f(approx) != 0): #hensel lift frobenius(a)
            approx = approx - f(approx)/g(approx)
        return approx

    def uniformizer_pow(self, n):
        """
        Returns the nth power of the uniformizer of self (as an element of self).

        EXAMPLES::

            sage: R.<a> = ZqCR(125)
            sage: R.uniformizer_pow(5)
            5^5 + O(5^25)
        """
        return self(self.prime_pow(n))

    def uniformizer(self):
        """
        Returns a uniformizer for this extension.

        Since this extension is unramified, a uniformizer for the
        ground ring will also be a uniformizer for this extension.

        EXAMPLES::

            sage: R.<a> = ZqCR(125)
            sage: R.uniformizer()
            5 + O(5^21)
        """
        return self(self.ground_ring().uniformizer())

    def _uniformizer_print(self):
        """
        Returns how the uniformizer is supposed to print.

        EXAMPLES::

            sage: R.<a> = Zq(125); R._uniformizer_print()
            '5'
        """
        return self.ground_ring()._uniformizer_print()

    def _unram_print(self):
        """
        Returns how the generator prints.

        EXAMPLES::

            sage: R.<a> = Zq(125); R._unram_print()
            'a'
        """
        return self.variable_name()

    def has_pth_root(self):
        r"""
        Returns whether or not `\ZZ_p` has a primitive `p^{\mbox{th}}` root of unity.

        Since adjoining a `p^{\mbox{th}}` root of unity yields a
        totally ramified extension, self will contain one if and only
        if the ground ring does.

        INPUT:

            - self -- a p-adic ring

        OUTPUT:

            - boolean -- whether self has primitive `p^{\mbox{th}}`
              root of unity.

        EXAMPLES::

            sage: R.<a> = Zq(1024); R.has_pth_root()
            True
            sage: R.<a> = Zq(17^5); R.has_pth_root()
            False
        """
        return self.ground_ring().has_pth_root()

    def has_root_of_unity(self, n):
        """
        Returns whether or not `\ZZ_p` has a primitive `n^{\mbox{th}}`
        root of unity.

        INPUT:

            - self -- a p-adic ring
            - n -- an integer

        OUTPUT:

            - boolean -- whether self has primitive `n^{\mbox{th}}`
              root of unity

        EXAMPLES::

            sage: R.<a> = Zq(37^8)
            sage: R.has_root_of_unity(144)
            True
            sage: R.has_root_of_unity(89)
            True
            sage: R.has_root_of_unity(11)
            False
        """
        if (self.prime() == 2):
            return n.divides(2*(self.residue_class_field().order()-1))
        else:
            return n.divides(self.residue_class_field().order() - 1)
