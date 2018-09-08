from sage.structure.element cimport Element
from sage.structure.element cimport MonoidElement
from sage.structure.element cimport CommutativeAlgebraElement

from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ

from sage.structure.element import coerce_binop


cdef class TateAlgebraTerm(MonoidElement):
    def __init__(self, parent, coeff, exponent):
        MonoidElement.__init__(self, parent)
        self._field = parent.base_ring().fraction_field()
        self._coeff = self._field(coeff)
        if isinstance(exponent, (list, tuple)):
            if len(exponent) != parent._ngens:
                raise ValueError("The number of exponents does not match the number of variables")
            self._exponent = tuple([ ZZ(e) for e in exponent ])
        else:
            self._exponent = (ZZ(exponent),) * parent._ngens

    def _repr_(self):
        parent = self._parent
        s = "(%s)" % self._coeff
        for i in range(parent._ngens):
            if self._exponent[i] == 1:
                s += "*%s" % parent._names[i]
            elif self._exponent[i] > 1:
                s += "*%s^%s" % (parent._names[i], self._exponent[i])
        return s

    def coefficient(self):
        return self._coeff

    def exponent(self):
        return self._exponent

    def _mul_(self, other):
        exponent = tuple([ self._exponent[i] + other.exponent()[i] for i in range(self._parent._ngens) ])
        return self._parent(self._coeff * other.coefficient(), exponent)

    def valuation(self):
        return self._coefficient.valuation() - sum(self._exponent[i]*self._log_radii[i] for i in range(self._parent._ngens))

    @coerce_binop
    def gcd(self, other):
        exponent = tuple([ min(self._exponent[i], other.exponent()[i]) for i in range(self._parent._ngens) ])
        val = min(self._coeff.valuation(), other.coefficient().valuation())
        return self._parent(self._field.uniformizer_pow(val), exponent)

    @coerce_binop
    def lcm(self, other):
        exponent = tuple([ max(self._exponent[i], other.exponent()[i]) for i in range(self._parent._ngens) ])
        val = max(self._coeff.valuation(), other.coefficient().valuation())
        return self._parent(self._field.uniformizer_pow(val), exponent)

    @coerce_binop
    def is_divisible_by(self, other):
        if self.valuation() < other.valuation():
            return False
        for i in range(self._parent._ngens):
            if self._exponent[i] < other.exponent()[i]:
                return False
        return True

    def divides(self, other):
        return other.is_divisible_by(self)

    @coerce_binop
    def __floordiv__(self, other):
        parent = self._parent
        if not parent.base_ring().is_field() and self.valuation() > other.valuation():
            raise ValueError("The division is not exact")
        exponent = [ ]
        for i in range(parent._ngens):
            if self._exponent[i] < other._exponent[i]:
                raise ValueError("The division is not exact")
            exponent.append(self._exponent[i] - other.exponent()[i])
        return parent(self._coeff / other.coefficient(), tuple(exponent))


cdef class TateAlgebraElement(CommutativeAlgebraElement):
    r"""
    Define the class for Tate series, elements of Tate algebras.

    Given a complete discrete valuation ring `R`, variables `X_1,\dots,X_k`
    and convergence radii `r_,\dots, r_n` in `\mathbb{R}_{>0}`, a Tate series is an
    element of the Tate algebra `R{X_1,\dots,X_k}`, that is a power series with
    coefficients `a_{i_1,\dots,i_n}` in `R` and such that
    `|a_{i_1,\dots,i_n}|*r_1^{-i_1}*\dots*r_n^{-i_n}` tends to 0 as
    `i_1,\dots,i_n` go towards infinity.

    INPUT:

    - ???


    """
    def __init__(self, parent, x, prec=None, reduce=True):
        r"""
        Initialize a Tate algebra element

        TESTS::

        
        """
        cdef TateAlgebraElement xc
        CommutativeAlgebraElement.__init__(self, parent)
        self._prec = Infinity
        S = parent._polynomial_ring
        if isinstance(x, TateAlgebraElement):
            xc = <TateAlgebraElement>x
            xparent = x.parent()
            if xparent is parent:
                self._poly = xc._poly
                self._prec = xc._prec
            elif xparent.variable_names() == parent.variable_names():
                ratio = parent._base.absolute_e() / xparent.base_ring().absolute_e()
                for i in range(parent.ngens()):
                    if parent.log_radii()[i] > xparent.log_radii()[i] * ratio:
                        raise ValueError("Cannot restrict to a bigger domain")
                self._poly = S(xc._poly)
                if xc._prec is Infinity:
                    self._prec = Infinity
                else:
                    self._prec = (xc._prec * ratio).ceil()
            else:
                raise TypeError("Variable names do not match")
        else:
            self._poly = S(x)
        if prec is not None:
            self._prec = min(self._prec, prec)
        if reduce:
            self._poly = self._poly.map_coefficients(lambda c: c.add_bigoh(self._prec))
        self._coeffs = None

    def _repr_(self):
        r"""
        Return a printable representation of a Tate series

        The terms are ordered with decreasing term order (increasing valuation of the coefficients, then the monomial order of the parent algebra).

        EXAMPLES::
        
            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: x + 2*x^2 + x^3
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: A(x+2*x^2+x^3,prec=5)
            (...00001)*x^3 + (...00001)*x + (...00010)*x^2 + O(2^5)

        """
        
        base = self._parent.base_ring()
        nvars = self._parent.ngens()
        vars = self._parent.variable_names()
        s = ""
        for (e,c) in self._coefficients():
            if c.valuation() >= self._prec:
                continue
            s += " + (%s)" % c
            for i in range(nvars):
                if e[i] == 1:
                    s += "*%s" % vars[i]
                elif e[i] > 1:
                    s += "*%s^%s" % (vars[i], e[i])
        if self._prec is not Infinity:
            # Only works with padics...
            s += " + %s" % base.change(show_prec="bigoh")(0, self._prec)
        if s == "":
            return "0"
        return s[3:]

    def polynomial(self):
        r"""
        Return a polynomial representation of the Tate series.

        For series which are not polynomials, the output corresponds to ???
        """
        # TODO: should we make sure that the result is reduced in some cases?
        return self._poly

    cpdef _add_(self, other):
        r"""
        Add a Tate series to the Tate series

        The precision of the output is adjusted to the minimum of both precisions.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: g = A(x + 2*y^2,prec=5); g
            (...00001)*x + (...00010)*y^2 + O(2^5)
            sage: h=f+g; h # indirect doctest
            (...00001)*x^3 + (...00010)*x^2 + (...00010)*y^2 + (...00010)*x + O(2^5)

        ::

            sage: f.precision_absolute()
            +Infinity
            sage: g.precision_absolute()
            5
            sage: h.precision_absolute()
            5


        """
        # TODO: g = x + 2*y^2 + O(2^5) fails coercing O(2^5) into R, is that normal? If we force conversion with R(O(2^5)) the result is a Tate series with precision infinity... it's understandable but confusing, no? 
        return self._parent(self._poly + other.polynomial(), min(self._prec, other.precision_absolute()))

    cpdef _neg_(self):
        r"""
        Return the opposite of the Tate series

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: -f # indirect doctest
            (...1111111111)*x^3 + (...1111111111)*x + (...11111111110)*x^2


        """

        return self._parent(-self._poly, self._prec, reduce=False)

    cpdef _sub_(self, other):
        r"""
        Return the difference of the Tate series and another

        The precision of the output is adjusted to the minimum of both precisions.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: g = A(x + 2*y^2,prec=5); g
            (...00001)*x + (...00010)*y^2 + O(2^5)
            sage: h=f-g; h # indirect doctest
            (...00001)*x^3 + (...00010)*x^2 + (...00010)*y^2 + (...00010)*x + O(2^5)

        ::

            sage: f.precision_absolute()
            +Infinity
            sage: g.precision_absolute()
            5
            sage: h.precision_absolute()
            5
        """
    
        return self._parent(self._poly - other.polynomial(), min(self._prec, other.precision_absolute()))

    cpdef _mul_(self, other):
        r"""
        Return the product of the Tate series and another.

        The precision is adjusted to match the best precision on the output.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: f = 2*x + 4*x^2 + 2*x^3; f
            (...00000000010)*x^3 + (...00000000010)*x + (...000000000100)*x^2
            g = A(x + 2*x^2, prec=5); g
            (...00001)*x + (...00010)*x^2 + O(2^5)
            sage: h = f*g; h # indirect doctest
            (...001010)*x^4 + (...000010)*x^2 + (...000100)*x^5 + (...001000)*x^3 + O(2^6)

        ::

            sage: f.precision_absolute()
            +Infinity
            sage: g.precision_absolute()
            5
            sage: h.precision_absolute()
            6


        """
        prec = min(self._prec + other.valuation(), other.precision_absolute() + self.valuation())
        return self._parent(self._poly * other.polynomial(), prec)

    cpdef _lmul_(self, Element right):
        r"""
        Multiply the series by an element of the base ring.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: R(2)*f # indirect doctest
            (...00000000010)*x^3 + (...00000000010)*x + (...000000000100)*x^2
            sage: R(6)*f # indirect doctest
            (...00000000110)*x^3 + (...00000000110)*x + (...000000001100)*x^2

        """
        other = self._parent._base(right)
        prec = self._prec + other.valuation()
        return self._parent(self._poly * other, prec)

    def _lshift(self, n):
        r"""
        Multiply the series by the `n`'th power of the uniformizer `\pi`.

        INPUT::

        - ``n`` - an integer (not necessarily non-negative)


        EXAMPLES::

            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: f << 2 # indirect doctest
            (...000000000100)*x^3 + (...000000000100)*x + (...0000000001000)*x^2

        If the base ring is not a field and we're shifting by a negative number of digits, the result is truncated -- that is, the output is the result of the integer division of the Tate series by `\pi^n`.

            sage: f << -1 # indirect doctest
            (...0000000001)*x^2
        """
        parent = self._parent
        base = parent.base_ring()
        if base.is_field():
            poly = self._poly.map_coefficients(lambda c: c << n)
            return parent(poly, self._prec + n, reduce=False)
        else:
            coeffs = { }
            field = base.fraction_field()
            ngens = parent.ngens()
            for exp, c in self._poly.dict().iteritems():
                minval = sum([ exp[i] * parent.log_radii()[i] for i in range(parent.ngens()) ]).ceil()
                coeffs[exp] = field(base(c) >> (minval-n)) << minval
            return parent(coeffs, self._prec + n, reduce=False)

    def __lshift__(self, n):
        r"""
        Multiply the series by the `n`'th power of the uniformizer `\pi`.

        INPUT::

        - ``n`` - an integer (not necessarily non-negative)


        EXAMPLES::

            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: f << 2 # indirect doctest
            (...000000000100)*x^3 + (...000000000100)*x + (...0000000001000)*x^2

        If the base ring is not a field and we're shifting by a negative number of digits, the result is truncated -- that is, the output is the result of the integer division of the Tate series by `\pi^n`.

            sage: f << -1 # indirect doctest
            (...0000000001)*x^2
        """
        
        return self._lshift(n)

    def __rshift__(self, n):
        r"""
        Divide the series by the `n`'th power of the uniformizer `\pi`.

        If the base ring is not a field and ``n`` is positive, the result is truncated -- that is, the output is the result of the integer division of the Tate series by `\pi^n`.


        INPUT:

        - ``n`` - an integer (not necessarily non-negative)


        EXAMPLES::

            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: f >> -2 # indirect doctest
            (...000000000100)*x^3 + (...000000000100)*x + (...0000000001000)*x^2
            sage: f >> 1 # indirect doctest
            (...0000000001)*x^2
        """
        return self._lshift(-n)

    def is_zero(self, prec=None):
        r"""
        Test whether the Tate algebra is 0.

        INPUT:

        - prec - (default: None) an integer which, if specifies, determines the precision of the test.

        EXAMPLES::

            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: f.is_zero()
            False
            sage: g = f << 4; g
            (...00000000010000)*x^3 + (...00000000010000)*x + (...000000000100000)*x^2
            sage: g.is_zero()
            False
            sage: g.is_zero(5)
            False
            sage: g.is_zero(4)
            True
        
        """
        if prec is None:
            # This assumes that self._poly is reduced
            return self._poly == 0
        else:
            return self.valuation() >= prec

    def inverse_of_unit(self):
        r"""
        Returns the inverse of the Tate series if invertible.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Ring with capped relative precision 10
            sage: f = A(1); f
            (...0000000001)
            sage: f.inverse_of_unit()
            (...0000000001) + O(2^10)
            sage: f = 2*x +1; f
            (...0000000001) + (...00000000010)*x
            sage: f.inverse_of_unit()
            (...0000000001) + (...1111111110)*x + (...0000000100)*x^2 + (...1111111000)*x^3 + (...0000010000)*x^4 + (...1111100000)*x^5 + (...0001000000)*x^6 + (...1110000000)*x^7 + (...0100000000)*x^8 + (...1000000000)*x^9 + O(2^10)
            sage: f = 1+x; f
            (...0000000001)*x + (...0000000001)
            sage: f.inverse_of_unit()
            Traceback (most recent call last)
            ...        
            ValueError: This series in not invertible
        
        """
        if not self.is_unit():
            raise ValueError("This series in not invertible")
        exp, c = self._coefficients()[0]
        parent = self._parent
        v, c = c.val_unit()
        cap = parent.precision_cap()
        inv = parent(c.inverse_of_unit(), prec=cap)
        x = self.add_bigoh(cap)
        prec = 1
        while prec < cap:
            prec *= 2
            inv = 2*inv - self*inv*inv
        return inv << v

    def is_unit(self):
        r"""
        Test whether the Tate series is invertible.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Ring with capped relative precision 10
            sage: f = 1; f
            1
            sage: f = A(1); f
            (...0000000001)
            sage: f.is_unit()
            True
            sage: f = 2*x +1; f
            (...0000000001) + (...000000000100)*x
            sage: f.is_unit()
            True
            sage: f = 1+x; f
            (...00000000010)*x + (...0000000001)
            sage: f.is_unit()
            False

        """
        if self.is_zero():
            return False
        exp, c = self._coefficients()[0]
        if max(exp) != 0:
            return False
        base = self.base_ring()
        if not base.is_field() and c.valuation() > 0:
            return False
        return True

    def restriction(self, log_radii):
        r"""
        Restrict the Tate series to a series converging on a smaller domain

        INPUT:

        - ``log_radii`` - the logs of the wanted convergence radii (see the definition of the class for the possible specifications)

        EXAMPLES::

        
        
        """
        # TODO: Seems broken
        parent = self._parent
        from sage.rings.tate_algebra import TateAlgebra
        ring = TateAlgebra(self.base_ring(), parent.variable_names(), log_radii, parent.precision_cap(), parent.term_order())
        return ring(self)

    def dict(self):
        r"""
        Returns a dictionary of the coefficients of the Tate series.

        The keys of the dictionary are the exponents and the values the
        coefficients.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Ring with capped relative precision 10
            sage: f = 1+2*x + 3*x*y^2; f
            (...0000000011)*x*y^2 + (...0000000001) + (...00000000010)*x
            sage: f.dict()
            {(0, 0): ...0000000001, (1, 0): ...00000000010, (1, 2): ...0000000011}

        
        """
        return self._poly.dict()

    def _coefficients(self):
        if self._coeffs is not None:
            return self._coeffs
        if self.is_zero():
            self._coeffs = []
        else:
            parent = self._parent
            self._coeffs = list(self.dict().iteritems())
            log_radii = parent.log_radii()
            ngens = parent.ngens()
            T = parent.term_order()
            sortkey = lambda c: (-c[1].valuation() + sum(c[0][i]*log_radii[i] for i in range(ngens)), T.sortkey(c[0]))
            self._coeffs.sort(key=sortkey, reverse=True)
        return self._coeffs

    def coefficients(self):
        return [ c for (c,exp) in self._coefficients ]

    def add_bigoh(self, n):
        return self._parent(self, prec=n)

    def precision_absolute(self):
        return self._prec

    def valuation(self):
        cap = self._parent.precision_cap()
        if self.is_zero():
            return cap
        else:
            parent = self._parent
            log_radii = parent.log_radii()
            ngens = parent.ngens()
            lt = self._coefficients()[0]
            val = lt[1].valuation() - sum(lt[0][i]*log_radii[i] for i in range(ngens))
            return min(cap, val)

    def precision_relative(self):
        return self._prec - self.valuation()

    def leading_term(self):
        if self.is_zero():
            return self._parent(0)
        lt = self._coefficients()[0]
        poly = self._parent._polynomial_ring.monomial(*lt[0])
        return self._parent(lt[1] * poly)

    def leading_coefficient(self):
        if self.is_zero():
            return self.base_ring(0)
        lt = self._coefficients()[0]
        return lt[1]
        
    def leading_monomial(self):
        if self.is_zero():
            return self._parent(0)
        lt = self._coefficients()[0]
        poly = self._parent._polynomial_ring.monomial(*lt[0])
        return self._parent(poly)

    def weierstrass_degree(self):
        v = self.valuation()
        return self.residue(v+1).degree()

    def degree(self):
        return self.weierstrass_degree()

    def weierstrass_degrees(self):
        v = self.valuation()
        return self.residue(v+1).degrees()

    def degrees(self):
        return self.weierstrass_degrees()

    def residue(self, n=1):
        for r in self._parent.radii():
            if r != 1:
                raise NotImplementedError("Residues are only implemented for radius 1")
        try:
            Rn = self.base_ring().residue_ring(n)
        except (AttributeError, NotImplementedError):
            Rn = self.base_ring().change(field=False, type="fixed-mod", prec=n)
        return self._poly.change_ring(Rn)

    def quo_rem(self, *divisors):
        parent = self._parent
        nvars = parent.ngens()
        ltds = [ d._coefficients()[0] for d in divisors ]
        # TODO: do everything on self._poly
        cap = parent.precision_cap()
        f = self.add_bigoh(cap)
        q = [ parent(0, cap - d.valuation()) for d in divisors ]
        r = parent(0, cap)
        while not f.is_zero():
            lt = f._coefficients()[0]
            found = False; i = 0
            while not found and i < len(divisors):
                ltd = ltds[i]
                is_divisible = ltd[1].valuation() <= lt[1].valuation()
                j = 0
                while is_divisible and j < nvars:
                    is_divisible = ltd[0][j] <= lt[0][j]
                    j += 1
                if is_divisible:
                    found = True
                    break
                i += 1
            if found:
                exponent = ( lt[0][j] - ltd[0][j] for j in range(nvars) )
                factor = (lt[1]//ltd[1]) * parent._polynomial_ring.monomial(*exponent)
                factor = parent(factor)
                f -= factor * divisors[i]
                q[i] += factor
            else:
                lt = f.leading_term()
                f -= lt; r += lt
        return q, r

    def __mod__(self, divisors):
        if not isinstance(divisors, (list, tuple)):
            divisors = [ divisors ]
        _, r = self.quo_rem(*divisors)
        return r

    def __floordiv__(self, divisors):
        if not isinstance(divisors, (list, tuple)):
            q, _ = self.quo_rem(divisors)
            return q[0]
        else:
            q, _ = self.quo_rem(*divisors)
            return q

    def reduce(self, I):
        g = I.groebner_basis()
        return self.quo_rem(g)

    @coerce_binop
    def Spoly(self, other):
        if self.is_zero() or other.is_zero():
            raise ValueError("Cannot compute the S-polynomial of zero")
        parent = self._parent
        S = parent._polynomial_ring
        se, sc = self._coefficients()[0]
        oe, oc = other._coefficients()[0]
        spoly = self._poly
        opoly = other.polynomial()
        for i in range(parent.ngens()):
            if se[i] > oe[i]:
                opoly *= S.gen(i) ** (se[i] - oe[i])
            elif se[i] < oe[i]:
                spoly *= S.gen(i) ** (oe[i] - se[i])
        if sc.valuation() > oc.valuation():
            opoly *= sc // oc
        else:
            spoly *= oc // sc
        return parent(spoly-opoly, min(self._prec, other.precision_absolute()))
