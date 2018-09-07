from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeAlgebra
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationRings
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.tate_algebra_element import TateAlgebraElement

from sage.categories.pushout import pushout

DEFAULT_CAP = 20

class TateAlgebra(CommutativeAlgebra, UniqueRepresentation):
    r"""Create a Tate series ring over a given complete discrete valuation
    ring.

    Given a complete discrete valuation ring `R`, variables `X_1,\dots,X_k`
    and convergence radii `r_,\dots, r_n` in `\mathbb{R}_{>0}`, the Tate
    algebra `R{X_1,\dots,X_k}` is the algebra of power series with
    coefficients `a_{i_1,\dots,i_n}` in `R` and such that
    `|a_{i_1,\dots,i_n}|*r_1^{-i_1}*\dots*r_n^{-i_n}` tends to 0 as
    `i_1,\dots,i_n` go towards infinity.


    INPUT:

    - ``base`` - a complete discrete valuation ring or field

    - ``names`` - names of the indeterminates

    - ``log-radii`` - (default: ``0``) the value(s) `\log(r_i)`. If only
      one number l is given, all `r_i`'s are defined with `\log(r_i)=l`.

    - ``prec`` - the default precision used if an exact object
      must be changed to an approximate object in order to do an
      arithmetic operation. If left as ``None``, it will be set to
      the precision of the base ring, if any. Otherwise,
      it will be set to the cap in precision of the base
      ring, if any. Otherwise, it will be set to the global
      default (20).

    - ``order`` - (default: ``degrevlex``) the monomial ordering 
      used to break ties when comparing terms with the same 
      coefficient valuation

    EXAMPLES::

        sage: R = Zp(2, 10, print_mode='digits'); R
        2-adic Ring with capped relative precision 10

    ::

        sage: A.<x,y> = TateAlgebra(R, order='lex'); A
        Tate Algebra in x, y over 2-adic Ring with capped relative precision 10

    The term ordering is used to determine how series are displayed. Terms
    are compared first according to the valuation of their coefficient, and
    ties are broken using the monomial ordering.

        sage: A.term_order()
        Lexicographic term order
        sage: f = 2 + y^5 + x^2; f
        (...0000000001)*x^2 + (...0000000001)*y^5 + (...00000000010)
        sage: B.<x,y> = TateAlgebra(R); B
        Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
        sage: B.term_order()
        Degree reverse lexicographic term order
        sage: B(f)
        (...0000000001)*y^5 + (...0000000001)*x^2 + (...00000000010)

    By order of priority, the precision cap is taken to be the given value, the
    precision of the base ring, the precision cap of the base ring, or the
    global default (20).  With a given value:

        sage: A.<x,y> = TateAlgebra(R,prec=5); A
        Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
        sage: A.precision_cap()
        5

    With the precision of the base ring:

        sage: A.<x,y> = TateAlgebra(R); A
        Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
        sage: A.precision_cap()
        10

    With the precision cap of the base ring:

        sage: S.<t> = LaurentSeriesRing(QQ,default_prec=15); S
        Laurent Series Ring in t over Rational Field
        sage: B.<x,y> = TateAlgebra(S); B
        Tate Algebra in x, y over Laurent Series Ring in t over Rational Field
        sage: B.precision_cap()
        15

    """

    def __init__(self, base, names, log_radii=QQ(0), prec=None, order='degrevlex'):
        """
        Initialize the Tate algebra

        TESTS::
        
        """
        if base not in CompleteDiscreteValuationRings() and base not in CompleteDiscreteValuationFields():
            raise TypeError("The base ring must be a complete discrete valuation ring or field")
        if isinstance(names, (list, tuple)):
            names = [ str(var) for var in names ]
        else:
            names = [ str(names) ]
        self._ngens = len(names)
        self.element_class = TateAlgebraElement
        CommutativeAlgebra.__init__(self, base, names, category=CommutativeAlgebras(base))
        if not isinstance(log_radii, (list, tuple)):
            self._log_radii = [ QQ(log_radii) ] * self._ngens
        elif len(log_radii) != self._ngens:
            raise ValueError("The number of radii does not match the number of variables")
        else:
            self._log_radii = [ QQ(r) for r in log_radii ]
        field = base.fraction_field()
        self._polynomial_ring = PolynomialRing(field, names, order=order)
        self._names = self._polynomial_ring.variable_names()
        self._order = self._polynomial_ring.term_order()
        self._gens = tuple([ self((field(1) << self._log_radii[i].ceil()) * self._polynomial_ring.gen(i)) for i in range(self._ngens) ])
        if prec is None:
            try:
                self._cap = base.precision_cap()
            except AttributeError:
                try:
                    self._cap = base.default_prec()
                except AttributeError:
                    self._cap = DEFAULT_CAP
        else:
            self._cap = ZZ(prec)

    # def _an_element_(self):
    #     r"""
    #     Return an element of the Tate series algebra

    #     EXAMPLES::
         
    #     """
    #     return self.element_class(0)

    def _coerce_map_from_(self, R):
        r"""
        Test whether the ring R can be coerced into the Tate algebra.

        R can be coerced into the Tate algebra if it can be coerced into its
        base ring, or if R is a Tate algebra whose base ring can be coerced into
        the base ring, which has the same variables with the same monomial
        order, and if the series of R are converging on a larger ball than those
        of the algebra.

        INPUT:

        - ``R`` - the ring to be coerced

        EXAMPLES::
        
            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10

        B can be coerced into A if B can be coerced into the base ring of A

            sage: B = ZZ; B
            Integer Ring
            sage: A.has_coerce_map_from(B) # indirect doctest
            True
            sage: B = GF(2); B
            Finite Field of size 2
            sage: A.has_coerce_map_from(B) # indirect doctest
            False

        B can be coerced into A if B is a Tate algebra over a base ring which
        can be coerced into the base ring of A, with the same variables and term
        orders, and a larger convergence radius:

            sage: S.<a> = Zq(4); S
            2-adic Unramified Extension Ring in a defined by x^2 + x + 1
            sage: B.<x,y> = TateAlgebra(S); B
            Tate Algebra in x, y over 2-adic Unramified Extension Ring in a defined by x^2 + x + 1
            sage: B.has_coerce_map_from(A) # indirect doctest
            True

        If the base ring of B cannot be coerced into the base ring of A,
        coercion cannot happen:
        
            sage: A.has_coerce_map_from(B) # indirect doctest
            False

        If B has different variables than A, coercion cannot happen:

            sage: B.<x,z> = TateAlgebra(R)
            sage: A.has_coerce_map_from(B) # indirect doctest
            False

        If B has a different term order than A, coercion cannot happen:

            sage: B.<x,y> = TateAlgebra(R,order="lex"); B
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: B.has_coerce_map_from(A) # indirect doctest
            False
            sage: B.<y,x> = TateAlgebra(R); B
            Tate Algebra in y, x over 2-adic Ring with capped relative precision 10
            sage: B.has_coerce_map_from(A) # indirect doctest
            False

        """
        base = self._base
        if base.has_coerce_map_from(R):
            return True
        if isinstance(R, TateAlgebra):
            Rbase = R.base_ring()
            if base.has_coerce_map_from(Rbase) and self._names == R.variable_names() and self._order == R.term_order():
                ratio = base.absolute_e() // Rbase.absolute_e()
                for i in range(self._ngens):
                    if self._log_radii[i] != R.log_radii()[i] * ratio:
                        return False
                return True

    def _ideal_class_(self, n):
        r"""
        Return the class of ideals in the Tate algebra

        INPUT:

        - ``n`` - number of generators

        EXAMPLE::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: A._ideal_class_(3)
            <class 'sage.rings.tate_algebra_ideal.TateAlgebraIdeal'>
        
        .. NOTE::

            The argument ``n`` is disregarded in the current implementation.
        """
        from sage.rings.tate_algebra_ideal import TateAlgebraIdeal
        return TateAlgebraIdeal

    #def _pushout_(self, R):
    #    if not isinstance(R, TateAlgebra):  # should we allow PolynomialRing as well?
    #        return None
    #    if self._names != R.variable_names():
    #        return None
    #    if self._order != R.term_order():
    #        return None
    #    Sbase = self._base
    #    Rbase = R.base_ring()
    #    base = pushout(Sbase, Rbase)
    #    Sratio = base.absolute_e() // Sbase.absolute_e()
    #    Rratio = base.absolute_e() // Rbase.absolute_e()
    #    log_radii = tuple([ min(self._log_radii[i] * Sratio, R.log_radii()[i] * Rratio) for i in range(self._ngens) ])
    #    cap = min(self._cap * Sratio, R.precision_cap() * Rratio)
    #    return TateAlgebra(base, self._names, log_radii, cap, self._order)

    def gen(self, n=0):
        r"""
        Return the ``n``'th generator of the algebra

        INPUT:

        - ``n`` - (default: 0) the generator to return

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: A.gen()
            (...0000000001)*x
            sage: A.gen(0)
            (...0000000001)*x
            sage: A.gen(1)
            (...0000000001)*y
            sage: A.gen(2)
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        
        """
        return self._gens[n]

    def gens(self):
        r"""
        Return the list of generators of the algebra

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: A.gens()
            ((...0000000001)*x, (...0000000001)*y)
        
        """
        return self._gens

    def ngens(self):
        """
        Return the number of generators of the algebra

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: A.ngens()
            2

        """
        return self._ngens

    def _repr_(self):
        """
        Return a printable representation of the algebra

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            
        """
        vars = ", ".join(self._names)
        return "Tate Algebra in %s over %s" % (vars, self.base_ring())

    def variable_names(self):
        """
        Return the list of the names of the variables of the algebra

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: A.variable_names()
            ('x', 'y')
        
        """
        return self._names

    def log_radii(self):
        """
        Return the list of the logs of the convergence radii of the series of the
        algebra.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: A.log_radii()
            [0, 0]
            sage: B.<x,y> = TateAlgebra(R,log_radii=1); B
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: B.log_radii()
            [1, 1]
            sage: C.<x,y> = TateAlgebra(R,log_radii=(1,-1)); C
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: C.log_radii()
            [1, -1]

        """
        
        return self._log_radii

    def term_order(self):
        """
        Return the monomial order used in the algebra

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10

        
        """
        return self._order

    def precision_cap(self):
        """
        Return the precision cap of the algebra

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: A.precision_cap()
            10
        
        """
        return self._cap

    def characteristic(self):
        """
        Return the characteristic of the base ring

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10
            sage: S.<t> = LaurentSeriesRing(QQ); S
            Laurent Series Ring in t over Rational Field
            sage: B.<x,y> = TateAlgebra(S); B
            Tate Algebra in x, y over Laurent Series Ring in t over Rational Field
            sage: B.characteristic()
            0
            sage: SS.<t> = LaurentSeriesRing(GF(2)); SS
            Laurent Series Ring in t over Finite Field of size 2
            sage: BB.<x,y> = TateAlgebra(SS); BB
            Tate Algebra in x, y over Laurent Series Ring in t over Finite Field of size 2
            sage: BB.characteristic()
            2

        """
        
        return self.base_ring().characteristic()

    #def ideal(self, gens):
    #    pass
