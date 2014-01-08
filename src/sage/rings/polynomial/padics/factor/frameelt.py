r"""
Representations of polynomials as sums of powers of OM approximations

A FrameEltTerm object represents a coefficient polynomial multiplied by the
approximation ``phi`` from the previous term.  In the first term, ``phi``
is replaced by the uniformizer, presening a valuation-unit representation
of a constant. The coefficient in FrameEltTerms beyond the first is a
FrameElt from the previous Frame.

A FrameElt is a polynominal in the current Frame as a sum of powers of the
current approximation multiplied by polynomial coefficients.. The
representation is a list of FrameEltTerms, each with a power of the
approximation and coefficient polynomial.

AUTHORS:

- Brian Sinclair (2012-02-22): initial version

"""

from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method

class FrameElt(SageObject):
    """
    Polynomials recursively represented by powers of the approximations
    $\phi_t(x)$ to a factor of $\Phi(x)$.

    INPUT:

    - ``frame`` -- The frame (and thus approximation) this refers to.

    - ``a`` -- Polynomial, default None; The polynomial to be represented
      by self.

    - ``this_exp`` -- int, default None; If ``this_exp`` is not None,
      then self is initialized as having a single term in its sum, namely
      ``a`` * ``phi`` ^ ``this_exp``

    EXAMPLES:

    If the FrameElt comes from the first frame, the term must be a constant::

        sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
        sage: from sage.rings.polynomial.padics.factor.frame import Frame
        sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
        sage: f = Frame(x^32+16); f.seed(x)
        sage: FrameElt(f,6)
        [3*2^1]

    Otherwise the we have an error (from FrameEltTerm)::

        sage: FrameElt(f,x+1)
        Traceback (most recent call last):
        ...
        TypeError: not a constant polynomial

    Moving to a higher frame and representing 6x^2 + 1.  Notice that the
    first FrameEltTerm represents 1 and the second (3 * 2^1)*x^2::

        sage: f = f.polygon[0].factors[0].next_frame(); f
        Frame with phi (1 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (1048574 + O(2^20))
        sage: FrameElt(f,6*x^2 + 1)
        [[1*2^0]phi1^0, [3*2^1]phi1^2]

    """

    def __init__(self,frame,a=None,this_exp=None):
        """
        Initializes self.

        See ``FrameElt`` for full documentation.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f = Frame(x^32+16); f.seed(x)
            sage: fe = FrameElt(f, 20)
            sage: fe.terms
            [5*2^2]
            sage: TestSuite(fe).run()
        """
        # deg(a) < deg(frame.phi)*frame.Eplus*frame.Fplus
        self.frame = frame
        if this_exp is not None:
            # initializes a FrameElt of phi_t ^ this_exp * a
            self.terms = [FrameEltTerm(self,a,this_exp)]
        elif a is None:
            self.terms = []
        elif frame.is_first():
            self.terms = [FrameEltTerm(self,a,0)]
        else:
            # Compute the phi-expansion of a
            a = self.frame.Ox(a)
            if self.frame.prev_frame().phi.degree() > a.degree():
                b = [a]
            else:
                b = []
                q, r = a.quo_rem(self.frame.prev_frame().phi)
                b.append(r)
                while q != 0:
                    q, r = q.quo_rem(self.frame.prev_frame().phi)
                    b.append(r)
            self.terms = [FrameEltTerm(self,b[i],i) for i in range(len(b)) if not b[i].is_zero()]

    def __cmp__(self, other):
        """
        Comparison.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f = Frame(x^32+16); f.seed(x)
            sage: fe = FrameElt(f, 20)
            sage: fe2 = FrameElt(f, 20)
            sage: fe == fe2
            True
            sage: fe == 0
            False
            sage: FrameElt(f) == 0
            False
        """
        c = cmp(type(self), type(other))
        if c: return c
        return cmp((self.frame, self.terms), (other.frame, other.terms))

    def is_single_term(self):
        """
        Returns ``True`` if the FrameElt is a sum of only one term, otherwise ``False``

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0,6)

        Since ``fe0`` represents `6`, it has a single term::

            sage: fe0.is_single_term()
            True

        We can create another ``FrameElt`` that has multiple terms::

            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,6*x^2 + 1)
            sage: fe1.is_single_term()
            False
        """
        if len(self.terms) == 1:
            if self.frame.prev is None:
                return True
            else:
                return self.terms[0]._coefficient.is_single_term()
        return False

    def valuation(self,purge=True):
        """
        Returns the valuation of self.

        The valuation of a sum (and thus a FrameElt) is the minimum of the
        valuations of its terms.

        INPUT:

        - ``purge`` -- Boolean, default True; If True, signals the method
          to remove any terms from the sum any terms of valuation higher
          than the minimum.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0,6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,6*x^2 + 1)
            sage: fe0.valuation()
            1
            sage: fe1.valuation(purge=False)
            0
            sage: fe1.terms
            [[1*2^0]phi1^0, [3*2^1]phi1^2]
            sage: fe1.valuation()
            0
            sage: fe1.terms
            [[1*2^0]phi1^0]
        """
        v = min([a.valuation() for a in self.terms]) if len(self.terms) > 0 else self.frame.O.precision_cap()
        if purge:
            self.terms = [a for a in self.terms if a.valuation() == v]
        return v

    def residue(self):
        """
        Returns the residue of this frame element.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import *
            sage: from sage.rings.polynomial.padics.factor.frameelt import *
            sage: k = ZpFM(5,40,'terse'); kx.<x> = k[]
            sage: f = kx([627500,6111375,13000,40,1])
            sage: t = OM_tree(f)[0]
            sage: t.prev.gamma_frameelt.residue()
            4*a0 + 1
            sage: e = t.prev.gamma_frameelt ** 2 + t.prev.gamma_frameelt
            sage: e.residue()
            3*a0
            sage: e.residue().parent()
            Finite Field in a0 of size 5^2
        """
        if not self.is_reduced():
            self = self.reduce()

        if self.frame.is_first():
            if self.terms[0]._exponent == 0:
                # unable to coerce in Zq
                return self.frame.R(self.terms[0]._unit)
            else:
                return self.frame.R(0)

        Eplus = int(self.frame.prev.segment.Eplus)
        if self.frame.prev.Fplus > 1:
            psi = self.frame.prev.segment.psi
            gamma = self.frame.prev.gamma
            residuelist = [gamma**(a._exponent/Eplus)*self.frame.prev.FFbase_elt_to_FF((psi**(a._exponent/Eplus)*a.value()).residue()) for a in self.terms if int(a._exponent) % Eplus == 0]
        else:
            residuelist = [a.value().residue() for a in self.terms if int(a._exponent) == 0]
        if len(residuelist) > 0:
            return sum(residuelist)
        else:
            return self.frame.R(0)

    def reduce(self):
        """
        Uses identities to fix the exponents of self to between
        zero and E+ times F+.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: k = ZpFM(2,40,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^32+16)[0]
            sage: f = t.phi**3
            sage: e = FrameElt(t,f); e
            [[[4294967295*2^8]phi1^1, [4294967295*2^8]phi1^3, [8589934591*2^7]phi1^7]phi2^0,
             [[8589934591*2^7]phi1^1, [17179869183*2^6]phi1^7]phi2^1, [[3*2^5]phi1^2,
             [3*2^5]phi1^4]phi2^2, [[3*2^4]phi1^2]phi2^3, [[274877906941*2^2]phi1^5]phi2^4,
             [[1*2^0]phi1^0]phi2^6]
            sage: e.reduce()
            [[[2147483649*2^9]phi1^1, [4294967295*2^8]phi1^3, [33*2^40]phi1^7]phi2^0,
             [[8589934591*2^7]phi1^1, [8589934593*2^7]phi1^7]phi2^1]
        """
        if self.frame.is_first():
            return self
        Eplus = self.frame.prev.segment.Eplus
        Fplus = self.frame.prev.Fplus
        reduced_elt = FrameElt(self.frame) # zero FrameElt
        for a in self.terms:
            if a._exponent >= Eplus * Fplus:
                b = FrameElt(self.frame)
                b.terms = [FrameEltTerm(b,a.value(),a._exponent - Eplus * Fplus)]
                b *= self.frame.prev.reduce_elt
                reduced_elt += b.reduce()
            elif a._exponent < 0:
                b = FrameElt(self.frame)
                b.terms = [FrameEltTerm(b,a.value(),a._exponent + Eplus * Fplus)]
                b *= (self.frame.prev.reduce_elt)**(-1)
                reduced_elt += b.reduce()
            else:
                summand = FrameElt(self.frame)
                summand.terms = [FrameEltTerm(reduced_elt,a.value().reduce(),a._exponent)]
                reduced_elt += summand
        return reduced_elt

    def is_reduced(self):
        """
        Returns ``True`` if all the exponents of all terms are less
        than E+, otherwise ``False``.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: k = ZpFM(2,40,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^32+16)[0]
            sage: f = t.phi**3
            sage: e = FrameElt(t,f); e.is_reduced()
            False
            sage: e = e.reduce()
            sage: e.is_reduced()
            True
        """
        return all([a.is_reduced() for a in self.terms])

    def find_denominator(self):
        """
        Returns the lowest power (possibly most negative power) of the
        uniformizer in self.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,20*x^2 + 6)
            sage: fe1.find_denominator()
            1
        """
        if self.frame.is_first():
            return self.terms[0]._exponent
        else:
            return min([fet._coefficient.find_denominator() for fet in self.terms])

    def polynomial(self,denominator=False):
        """
        Returns ``self`` as a polynomial, optionally with the power of the
        uniformizer present in its denominator.

        INPUT:

        - ``denominator`` -- Boolean, default False.  If True, returns the
          polynomial of ``self`` multiplied by the uniformizer to the highest
          power it appears in the denominator of ``self`` and that power as
          a tuple.

        OUTPUT:

        - If ``denominator`` is False, ``self`` as a polynomial.

        - If ``denominator`` is True, returns the tuple of ``self`` multiplied
          to not have a denominator and the power of the uniformizer required
          to clear the denominator.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,6*x^2 + 1)
            sage: fe1.polynomial()
            (6 + O(2^20))*x^2 + (0 + O(2^20))*x + (1 + O(2^20))
        """
        if min([a._exponent for a in self.terms]) < 0:
            raise ValueError("Cannot cconstruct polynomial representation of FrameElt with negative exponents")
        if denominator:
            piexp = self.find_denominator()
            if piexp < 0:
                return (self * FrameElt(self.frame,self.frame.Ox(self.frame.O.uniformizer_pow(-piexp)))).polynomial(),-piexp
            else:
                return self.polynomial(),0
        else:
            if self.frame.is_first():
                return self.frame.O.uniformizer()**int(self.terms[0]._exponent)*self.terms[0]._unit
            else:
                return sum([self.frame.prev_frame().phi**int(a._exponent)*a._coefficient.polynomial() for a in self.terms])

    def __neg__(self):
        """
        Negates a FrameElt.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,20*x^2+6); fe1
            [[3*2^1]phi1^0, [5*2^2]phi1^2]
            sage: -fe1
            [[524285*2^1]phi1^0, [262139*2^2]phi1^2]
        """
        if self.frame.is_first():
            return FrameElt(self.frame,-self.polynomial())
        else:
            neg = FrameElt(self.frame)
            neg.terms = [-r for r in self.terms]
            return neg

    def __sub__(self,r):
        """
        Subtraction.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,20*x^2); fe2 = FrameElt(f1,6)
            sage: fe1 - fe2 # indirect doctest
            [[524285*2^1]phi1^0, [5*2^2]phi1^2]
        """
        return self.__add__(r.__neg__())

    def __radd__(self,l):
        """
        Alias for addition.

        TESTS::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,20*x^2); fe2 = FrameElt(f1,6)
            sage: fe1 + fe2 # indirect doctest
            [[3*2^1]phi1^0, [5*2^2]phi1^2]
        """
        return self.__add__(l)

    def __add__(self,r):
        """
        Addition.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,20*x^2); fe2 = FrameElt(f1,6)
            sage: fe1 + fe2 # indirect doctest
            [[3*2^1]phi1^0, [5*2^2]phi1^2]
        """
        # For using sum() command, must be able to be added to int(0)
        if isinstance(r,int) and r == 0:
            return self
        if self.frame.phi != r.frame.phi:
            raise ValueError("Cannot add FrameElts with different values of phi")
        if len(self.terms) == 0:
            return r
        sum = FrameElt(self.frame)
        if self.frame.prev is None:
            v = min(self.terms[0]._exponent,r.terms[0]._exponent)
            pi = self.frame.O.uniformizer()
            usum = self.terms[0]._unit * pi ** (self.terms[0]._exponent - v)
            usum = usum + r.terms[0]._unit * pi ** (r.terms[0]._exponent - v)
            sum.terms = [FrameEltTerm(sum,usum,v)]
        else:
            if self.terms == []:
                for k in range(len(r.terms)):
                    sum.terms.append(FrameEltTerm(sum,r.terms[k].value(),r.terms[k]._exponent))
            elif r.terms == []:
                for k in range(len(self.terms)):
                    sum.terms.append(FrameEltTerm(sum,self.terms[k].value(),self.terms[k]._exponent))
            else:
                i = 0 ; j = 0
                while i < len(self.terms) and j < len(r.terms):
                    if self.terms[i]._exponent < r.terms[j]._exponent:
                        sum.terms.append(FrameEltTerm(sum,self.terms[i].value(),self.terms[i]._exponent))
                        i = i + 1
                    elif self.terms[i]._exponent > r.terms[j]._exponent:
                        sum.terms.append(FrameEltTerm(sum,r.terms[j].value(),r.terms[j]._exponent))
                        j = j + 1
                    elif self.terms[i]._exponent == r.terms[j]._exponent:
                        sum.terms.append(FrameEltTerm(sum,self.terms[i].value()+r.terms[j].value(),self.terms[i]._exponent))
                        i = i + 1; j = j + 1
                if j < len(r.terms):
                    for k in range(j,len(r.terms)):
                        sum.terms.append(FrameEltTerm(sum,r.terms[k].value(),r.terms[k]._exponent))
                elif i < len(self.terms):
                    for k in range(i,len(self.terms)):
                        sum.terms.append(FrameEltTerm(sum,self.terms[k].value(),self.terms[k]._exponent))
        return sum

    def __rmul__(self,l):
        """
        Alias for multiplication.

        TESTS::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,20*x^2); fe2 = FrameElt(f1,6)
            sage: fe1 * fe2 # indirect doctest
            [[15*2^3]phi1^2]
        """
        return self.__mul__(l)

    def __mul__(self,r):
        """
        Multiplication.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,20*x^2); fe2 = FrameElt(f1,6)
            sage: fe1 * fe2 # indirect doctest
            [[15*2^3]phi1^2]
        """
        if isinstance(r,int) and r == 0:
            return self
        if self.frame.depth != r.frame.depth:
            raise ValueError("Cannot multiply FrameElts with different frame depths")
        product = FrameElt(self.frame)
        if self.frame.prev is None:
            v = self.terms[0]._exponent + r.terms[0]._exponent
            uprod = self.terms[0]._unit
            uprod = uprod * r.terms[0]._unit
            product.terms = [FrameEltTerm(product,uprod,v)]
        else:
            for i in range(len(self.terms)):
                summand = FrameElt(self.frame)
                for j in range(len(r.terms)):
                    summand.terms.append(FrameEltTerm(product,self.terms[i].value()*r.terms[j].value(),self.terms[i]._exponent+r.terms[j]._exponent))
                product = product + summand
        return product

    def __pow__(self,n):
        """
        Raise ``self`` to the integer power ``n``.

        Only single term FrameElts with single terms in all recursive FrameElts
        and FrameEltTerms can be raised this way.

        EXAMPLES::

        Building the needed framework and squaring 6 as a FrameElt::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f = Frame(x^32+16); f.seed(x)
            sage: fe = FrameElt(f,6); fe
            [3*2^1]
            sage: fe.polynomial()
            6 + O(2^20)
            sage: fe.__pow__(2)
            [9*2^2]
            sage: fe ** 2
            [9*2^2]
            sage: (fe**2).polynomial()
            36 + O(2^20)

        Moving to a higher frame and squaring 6x^2 as a FrameElt::

            sage: f = f.polygon[0].factors[0].next_frame()
            sage: fe = FrameElt(f,6*x**2); fe
            [[3*2^1]phi1^2]
            sage: fe.polynomial()
            (6 + O(2^20))*x^2 + (0 + O(2^20))*x + (0 + O(2^20))
            sage: fe.__pow__(2)
            [[9*2^2]phi1^4]
            sage: fe**2
            [[9*2^2]phi1^4]
            sage: (fe**2).polynomial()
            (36 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (0 + O(2^20))

        As soon as we are past the first frame, we must take care not to
        try to take powers of non-single-term FrameElts::

            sage: fe = FrameElt(f,6*x^2+1); fe
            [[1*2^0]phi1^0, [3*2^1]phi1^2]
            sage: fe ** 2
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot take a power of a non-single term FrameElt

        """
        if not self.is_single_term():
            raise NotImplementedError("Cannot take a power of a non-single term FrameElt")
        else:
            product = FrameElt(self.frame)
            product.terms = [self.terms[0]**n]
            return product

    def __div__(self,right):
        """
        Division.

        One can only divide by single term FrameElts.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,20*x^2); fe2 = FrameElt(f1,6)
            sage: fe1 / fe2 # indirect doctest
            [[349527*2^1]phi1^2]
        """
        if not right.is_single_term():
            raise NotImplementedError("Cannot divide by a non-single term FrameElt")
        else:
            quotient = FrameElt(self.frame)
            quotient.terms = [a / right.terms[0] for a in self.terms]
            return quotient

    def __getitem__(self,i):
        """
        Get the `i`th FrameEltTerm.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0,6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,6*x^2 + 1)
            sage: fe0[0] # indirect doctest
            3*2^1
            sage: fe1[0]
            [1*2^0]phi1^0
            sage: fe1[1]
            [3*2^1]phi1^2
        """
        return self.terms[i]

    def __repr__(self):
        """
        String representation.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0,6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1,6*x^2 + 1)
            sage: fe1 # indirect doctest
            [[1*2^0]phi1^0, [3*2^1]phi1^2]
        """
        return repr(self.terms)

class FrameEltTerm(SageObject):
    """
    A single term of the sum of powers of OM representations.

    A FrameEltTerm object should be generated and manipulated by a parent
    FrameElt to which it belongs.

    If the parent FrameElt belongs to the first frame, the FrameEltTerm
    holds a constnat, namely a * pi ^ e.  For frames beyond the first,
    a FrameEltTerm contains the exponent of the previous approximation
    and a FrameElt of the previous frame as a coefficient.

    INPUT:

    - ``frelt`` -- FrameElt. The sum to which this term belongs.

    - ``a`` -- The coefficient of this term.

    - ``e`` -- The exponent of this term.

    EXAMPLES::

    If the parent FrameElt comes from the first frame, the term is a constant::

        sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt,FrameEltTerm
        sage: from sage.rings.polynomial.padics.factor.frame import Frame
        sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
        sage: f = Frame(x^32+16); f.seed(x)
        sage: elt = FrameElt(f)
        sage: FrameEltTerm(elt,3,2)
        3*2^2

    If the uniformizer divides a constant, the FrameEltTerm corrects the exponent::

        sage: FrameEltTerm(elt,6,0)
        3*2^1
        sage: FrameEltTerm(elt,4,0)
        1*2^2

    Moving to a higher frame and representing 6x^2 (or 6 * phi ^ 2)::

        sage: f = f.polygon[0].factors[0].next_frame(); f
        Frame with phi (1 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (1048574 + O(2^20))
        sage: elt = FrameElt(f)
        sage: FrameEltTerm(elt,6,2)
        [3*2^1]phi1^2
        sage: FrameElt(f,6*x^2)[0]
        [3*2^1]phi1^2

    """
    def __init__(self,frelt,a,e):
        """
        Initializes self.

        See ``FrameEltTerm`` for full documentation.

        TESTS::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt,FrameEltTerm
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f = Frame(x^32+16); f.seed(x)
            sage: fe = FrameElt(f)
            sage: fet = FrameEltTerm(fe,3,2)
            sage: TestSuite(fet).run()
        """
        self.frameelt = frelt
        self._scalar_flag = (self.frameelt.frame.prev is None)
        self._exponent = int(e)
        self._zero_flag = False

        if a in self.frameelt.frame.Ox or a in self.frameelt.frame.O:
            if isinstance(a,int):
                a = self.frameelt.frame.Ox(a)
            self._zero_flag = a.is_zero() #cannot be replaced by a == 0
            if self._scalar_flag:
                a = self.frameelt.frame.O(a)
                if self._zero_flag:
                    self._unit = a
                else:
                    self._unit = a.unit_part()
                if a.valuation() > 0 and a.valuation():
                    self._exponent += a.valuation()
            else:
                self._coefficient = FrameElt(self.frameelt.frame.prev_frame(),a)
        else:
            self._coefficient = a
            a._zero_flag = False

    def __cmp__(self, other):
        """
        Comparison.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt,FrameEltTerm
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: fet0 = FrameEltTerm(fe0,3,2)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1, 4*x^2 + 2)
            sage: fet1 = FrameEltTerm(fe1,3,2)
            sage: fet0 == fet1
            False
            sage: fet0 == FrameEltTerm(fe0,4,2)
            False
            sage: fet0 == 12
            False
        """
        c = cmp(type(self), type(other))
        if c: return c
        c = cmp((self.frameelt.frame, self._scalar_flag), (other.frameelt.frame, other._scalar_flag))
        if c: return c
        if self._scalar_flag:
            return cmp((self._exponent, self._unit), (other._exponent, other._unit))
        else:
            return cmp((self._exponent, self._coefficient), (other._exponent, other._coefficient))

    @cached_method
    def valuation(self):
        """
        Returns the valuation of this term.

        The valuation is given by the slope of the previous segment
        times the exponent, plus the valuation of the coefficient.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt,FrameEltTerm
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1)
            sage: fet1 = FrameEltTerm(fe1,3,2); fet1
            [3*2^0]phi1^2
            sage: fet1.valuation()
            1/4
            sage: f1.prev.segment.slope
            1/8
        """
        if self.frameelt.frame.prev is None:
            return Rational(self._exponent)
        else:
            return self.frameelt.frame.prev.segment.slope * self._exponent + self._coefficient.valuation()

    def reduce(self):
        """
        Uses identities to fix the exponent of self to be less than
        E+ times F+.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: k = ZpFM(2,40,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^32+16)[0]
            sage: f = t.phi**3
            sage: e = FrameElt(t,f); e.terms
            [[[4294967295*2^8]phi1^1, [4294967295*2^8]phi1^3, [8589934591*2^7]phi1^7]phi2^0,
             [[8589934591*2^7]phi1^1, [17179869183*2^6]phi1^7]phi2^1,
             [[3*2^5]phi1^2, [3*2^5]phi1^4]phi2^2,
             [[3*2^4]phi1^2]phi2^3,
             [[274877906941*2^2]phi1^5]phi2^4,
             [[1*2^0]phi1^0]phi2^6]

        Note that reducing the terms is not the same as reducing the element::

            sage: for a in e.terms: a.reduce()
            sage: e.terms
            [[[4294967295*2^8]phi1^1, [4294967295*2^8]phi1^3, [8589934591*2^7]phi1^7]phi2^0,
             [[8589934591*2^7]phi1^1, [17179869183*2^6]phi1^7]phi2^1,
             [[3*2^7]phi1^7, [3*2^7]phi1^9]phi2^0,
             [[3*2^6]phi1^7]phi2^1,
             [[274877906941*2^6]phi1^15]phi2^0,
             [[1*2^6]phi1^15]phi2^0]

            sage: ee = FrameElt(t, f); ee.reduce()
            [[[2147483649*2^9]phi1^1, [4294967295*2^8]phi1^3, [33*2^40]phi1^7]phi2^0,
             [[8589934591*2^7]phi1^1, [8589934593*2^7]phi1^7]phi2^1]
        """
        if self.frameelt.frame.prev is None:
            return
        Eplus = self.frameelt.frame.prev.segment.Eplus
        Fplus = self.frameelt.frame.prev.Fplus

        if self._exponent >= Eplus * Fplus:
            q,r = Integer(self._exponent).quo_rem(int(Eplus))
            self._exponent = r
            self._coefficient *= (self.frameelt.frame.prev.segment.psi ** (q*Fplus))
        self._coefficient = self._coefficient.reduce()
        return

    def is_reduced(self):
        """
        Returns ``True`` if all the exponents of all terms are less
        than E+, otherwise ``False``.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: k = ZpFM(2,40,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^32+16)[0]
            sage: f = t.phi**3
            sage: e = FrameElt(t,f)
            sage: [a.is_reduced() for a in e.terms]
            [True, True, False, False, False, False]
        """
        if self.frameelt.frame.prev is None:
            return True
        if self._exponent < self.frameelt.frame.prev.segment.Eplus:
            return self._coefficient.is_reduced()
        return False

    def is_scalar(self):
        """
        Returns ``True`` if the term is just a scalar, otherwise ``False``.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt,FrameEltTerm
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: fet0 = FrameEltTerm(fe0,3,2)
            sage: fet0.is_scalar()
            True

        Terms in later frames are not scalars::

            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1)
            sage: fet1 = FrameEltTerm(fe1,3,2)
            sage: fet1.is_scalar()
            False
        """
        return self._scalar_flag

    def is_zero(self):
        """
        Returns ``True`` if the term is equal to 0, otherwise ``False``.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt,FrameEltTerm
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0)
            sage: fet0 = FrameEltTerm(fe0,3,2)
            sage: fet0.is_zero()
            False
            sage: fet1 = FrameEltTerm(fe0,0,4)
            sage: fet1.is_zero()
            True
        """
        return self._zero_flag

    def is_single_term(self):
        """
        Returns ``True`` if the term is does not have more than one term
        at any recursive level, otherwise ``False``.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: k = ZpFM(2,40,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^32+16)[0]
            sage: f = t.phi**3
            sage: e = FrameElt(t,f); e.terms
            [[[4294967295*2^8]phi1^1, [4294967295*2^8]phi1^3, [8589934591*2^7]phi1^7]phi2^0,
             [[8589934591*2^7]phi1^1, [17179869183*2^6]phi1^7]phi2^1,
             [[3*2^5]phi1^2, [3*2^5]phi1^4]phi2^2,
             [[3*2^4]phi1^2]phi2^3,
             [[274877906941*2^2]phi1^5]phi2^4,
             [[1*2^0]phi1^0]phi2^6]
            sage: [a.is_single_term() for a in e.terms]
            [False, False, False, True, True, True]
        """
        if self._scalar_flag:
            return True
        else:
            return self._coefficient.is_single_term()

    def value(self):
        """
        Returns the coeffecient part of the term.  For scalars, this
        is a single number (the unit).  For polynomials, this is the
        FrameElt representing the polynomial coefficient.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: k = ZpFM(2,40,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^32+16)[0]
            sage: f = t.phi**3
            sage: e = FrameElt(t,f)
            sage: et = e.terms[0]; et
            [[4294967295*2^8]phi1^1, [4294967295*2^8]phi1^3, [8589934591*2^7]phi1^7]phi2^0
            sage: et.value()
            [[4294967295*2^8]phi1^1, [4294967295*2^8]phi1^3, [8589934591*2^7]phi1^7]
            sage: et = et.value().terms[0]; et
            [4294967295*2^8]phi1^1
            sage: et.value()
            [4294967295*2^8]
            sage: et = et.value().terms[0]; et
            4294967295*2^8
            sage: et.value()
            4294967295 + O(2^40)
        """
        if self._scalar_flag:
            return self._unit
        else:
            return self._coefficient

    #def __add__(self,right):
    #def __mul__(self,right):
    #    We don't add or multiply on FrameEltTerms directly -- the parent FrameElt does this

    def __pow__(self,n):
        """
        Raise ``self`` to the integer power ``n``.

        Only FrameEltTerms with single terms in all recursive FrameElts
        and FrameEltTerms can be raised this way.

        EXAMPLES::

        Building the needed framework and squaring 12 as a FrameEltTerm::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt,FrameEltTerm
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f = Frame(x^32+16); f.seed(x)
            sage: fet = FrameEltTerm(FrameElt(f),3,2); fet
            3*2^2
            sage: fet.__pow__(2)
            9*2^4
            sage: fet**2
            9*2^4

        Moving to a higher frame and squaring 12x^2 as a FrameEltTerm::

            sage: f = f.polygon[0].factors[0].next_frame()
            sage: fet = FrameEltTerm(FrameElt(f),12,2); fet
            [3*2^2]phi1^2
            sage: fet.__pow__(2)
            [9*2^4]phi1^4
            sage: fet**2
            [9*2^4]phi1^4

        Starting at a depth of 2, we must take care not to try to take
        powers of non-single-term FrameEltTerms::

            sage: f = f.polygon[0].factors[0].next_frame()
            sage: f = f.polygon[0].factors[0].next_frame()
            sage: f.depth
            2
            sage: fet = FrameElt(f,x+1)[0]; fet
            [[1*2^0]phi1^0, [1*2^0]phi1^1]phi2^0
            sage: fet ** 2 # Error
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot take a power of a non-single term FrameEltTerm

        """
        if not self.is_single_term():
            raise NotImplementedError("Cannot take a power of a non-single term FrameEltTerm")
        else:
            n = int(n)
            if self._scalar_flag:
                return FrameEltTerm(self.frameelt, self._unit**n, self._exponent*n)
            else:
                return FrameEltTerm(self.frameelt, self._coefficient**n, self._exponent*n)

    def __div__(self,right):
        """
        Division.

        One can only divide by single term FrameEltTerms.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt
            sage: k = ZpFM(2,40,'terse'); kx.<x> = k[]
            sage: t = OM_tree(x^32+16)[0]
            sage: f = t.phi**3
            sage: e = FrameElt(t,f)
            sage: e.terms[0] / e.terms[-1] # indirect doctest
            [[4294967295*2^8]phi1^1, [4294967295*2^8]phi1^3, [8589934591*2^7]phi1^7]phi2^-6
        """
        if not right.is_single_term():
            raise NotImplementedError("Cannot divide by a non-single term FrameEltTerm")
        else:
            if self._scalar_flag:
                return FrameEltTerm(self.frameelt, self._unit/right._unit, self._exponent-right._exponent)
            else:
                return FrameEltTerm(self.frameelt, self._coefficient/right._coefficient, self._exponent-right._exponent)

    def __neg__(self):
        """
        Negation.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt,FrameEltTerm
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f = Frame(x^32+16); f.seed(x)
            sage: f = f.polygon[0].factors[0].next_frame();
            sage: elt = FrameElt(f)
            sage: fet = FrameEltTerm(elt,6,2); fet
            [3*2^1]phi1^2
            sage: FrameElt(f,6*x^2)[0]
            [3*2^1]phi1^2
            sage: fet.__neg__()
            [524285*2^1]phi1^2
            sage: -fet
            [524285*2^1]phi1^2
            sage: FrameElt(f,-6*x**2)[0]
            [524285*2^1]phi1^2

        """
        if self.frameelt.frame.is_first():
            return FrameEltTerm(self.frameelt,-self._unit,self._exponent)
        else:
            return FrameEltTerm(self.frameelt,-self._coefficient,self._exponent)

    def __repr__(self):
        """
        String representation.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frameelt import FrameElt,FrameEltTerm
            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: f0 = Frame(x^32+16); f0.seed(x)
            sage: fe0 = FrameElt(f0, 6)
            sage: repr(FrameEltTerm(fe0,0,4)) # indirect doctest
            '0'
            sage: fet0 = FrameEltTerm(fe0,3,2)
            sage: repr(fet0)
            '3*2^2'
            sage: f1 = f0.polygon[0].factors[0].next_frame()
            sage: fe1 = FrameElt(f1)
            sage: fet1 = FrameEltTerm(fe1,3,2)
            sage: repr(fet1)
            '[3*2^0]phi1^2'
        """
        if self._zero_flag:
            return "0"
        if self._scalar_flag:
            return repr(self._unit.lift())+'*'+repr(self.frameelt.frame.O.uniformizer().lift())+'^'+repr(self._exponent)
        else:
            return repr(self._coefficient)+'phi'+repr(self.frameelt.frame.depth)+'^'+repr(self._exponent)
