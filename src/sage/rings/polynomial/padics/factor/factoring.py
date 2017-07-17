r"""
Methods related to factoring local field polynomials using an OM algorithm

AUTHORS:

- Brian Sinclair (2012-02-22): initial version

"""
from sage.rings.padics.factory import ZpFM, Zp
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.padics.factor.frame import Frame
from sage.structure.factorization import Factorization

def pfactor_non_monic(f):
    r"""
    Wrapper around pfactor to handle non-monic polynomial

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.factor.factoring import pfactor_non_monic
        sage: R.<x> = ZpFM(3,7)[]
        sage: pfactor_non_monic(3*x^2 - 15)
        (3 + O(3^7)) * ((1 + O(3^7))*x^2 + (O(3^7))*x + (1 + 3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + O(3^7)))
        sage: pfactor_non_monic(3*x^2 + 15)
        (3 + O(3^7)) * ((1 + O(3^7))*x + (2 + 2*3^2 + 3^4 + 3^6 + O(3^7))) * ((1 + O(3^7))*x + (1 + 2*3 + 2*3^3 + 3^4 + 2*3^5 + 3^6 + O(3^7)))
    """
    def make_monic(f):
        """
        Transform a polynomial over a local ring into a monic one for
        factoring.

        INPUT:

        - ``f`` -- a polynomial over a local ring
        """

        Kx = f.parent()
        K = Kx.base()
        pi = K.uniformizer()

        if f.leading_coefficient() == K(1):
            return f,K(1),0,0

        minval = min([a.valuation() for a in f.coeffs()])
        divi = pi ** minval
        g = Kx([a // divi for a in f.coeffs()])
        multval = -minval

        lead = g.leading_coefficient()
        vlead = lead.valuation()
        seq = g.coeffs()
        exp = -1
        expmult = g.degree()*exp-vlead
        if not g.degree() <= 0:
            while expmult < 0 or min([seq[i].valuation()-exp*i+expmult for i in range(1,len(seq))]) < 0:
                exp += 1
                expmult = g.degree()*exp-vlead
        else:
            exp = 0
            expmult = -vlead

        g = g * pi ** expmult
        g = Kx([g.coeffs()[i] // pi**(exp*i) for i in range(g.degree()+1)])
        multval = multval + expmult

        uni = g.leading_coefficient()
        g = Kx([a // uni for a in g.coeffs()])

        return g,uni,multval,exp

    def undo_monic(f,exp,h):
        """
        Transform back to a monic polynomial
        """
        return h(f.parent().base().uniformizer() ** exp * f.parent().gen())

    def dist_den(facts,multval):
        if len(facts) == 0:
            return []
        Kx = facts[0].parent()
        K = Kx.base_ring()
        Kpi = facts[0].parent().base_ring().uniformizer()
        newfacts = []
        for fact in facts:
            m = min([a.valuation() for a in fact.coeffs()])
            newfacts.append(Kx([a // Kpi**m for a in fact.coeffs()]))
            multval -= m
        if multval > 0:
            # This REALLY should not happen
            raise StandardError, "Power of pi did not distribute over factors of non-monic polynomial."
        return newfacts

    Kx = f.parent()
    g,uni,multval,exp = make_monic(f)
    if multval > 0:
        pi = f.base_ring().uniformizer()
        prec = f.base_ring().precision_cap()
        Ly = PolynomialRing(ZpFM(pi,prec+2*multval),names='y')
        f = Ly(f)
        g,uni,multval,exp = make_monic(f)
    facts = [fact[0] for fact in pfactor(g)]
    facts = [undo_monic(f,exp,fact) for fact in facts]
    facts = dist_den(facts,multval)
    if multval > 0:
        uni = Kx(uni)
        facts = [Kx(fact) for fact in facts]
    if multval < 0:
        return Factorization([(f.base_ring().uniformizer(),-multval)] + [(fact,1) for fact in facts],unit=uni)
    return Factorization([(fact,1) for fact in facts],unit=uni)

def pfactor(Phi):
    r"""
    Return a factorization of Phi

    This is accomplished by constructing a tree of Frames of approximations
    of factors.

    INPUT:

    - ``Phi`` -- squarefree padic polynomial with fixed precision coefficients

    OUTPUT:

    - A Factorization of the polynomial ``Phi``

    EXAMPLES:

    Factoring polynomials over Zp(2)[x]::

        sage: from sage.rings.polynomial.padics.factor.factoring import jorder4,pfactor
        sage: f = ZpFM(2,50,'terse')['x']( (x^32+16)*(x^32+16+2^16*x^2)+2^34 )
        sage: factors = pfactor(f); len(factors) # long time (5.7s)
        2

    See the irreducibility of x^32+16 in Zp(2)[x]::

        sage: pfactor(ZpFM(2)['x'](x^32+16))
        (1 + O(2^20))*x^32 + (O(2^20))*x^31 + (O(2^20))*x^30 + (O(2^20))*x^29 + (O(2^20))*x^28 + (O(2^20))*x^27 + (O(2^20))*x^26 + (O(2^20))*x^25 + (O(2^20))*x^24 + (O(2^20))*x^23 + (O(2^20))*x^22 + (O(2^20))*x^21 + (O(2^20))*x^20 + (O(2^20))*x^19 + (O(2^20))*x^18 + (O(2^20))*x^17 + (O(2^20))*x^16 + (O(2^20))*x^15 + (O(2^20))*x^14 + (O(2^20))*x^13 + (O(2^20))*x^12 + (O(2^20))*x^11 + (O(2^20))*x^10 + (O(2^20))*x^9 + (O(2^20))*x^8 + (O(2^20))*x^7 + (O(2^20))*x^6 + (O(2^20))*x^5 + (O(2^20))*x^4 + (O(2^20))*x^3 + (O(2^20))*x^2 + (O(2^20))*x + (2^4 + O(2^20))

    Test the irreducibility of test polynomial jorder4 for Zp(3)::

        sage: len(pfactor(jorder4(3))) == 1
        True

    Factor jorder4 for Zp(5) and Zp(7) and check that the products return
    the original::

        sage: ff = pfactor(jorder4(5))
        sage: prod(ff) == jorder4(5)
        True
        sage: ff = pfactor(jorder4(7))
        sage: prod(ff) == jorder4(7)
        True

    REFERENCES:

        S. Pauli, Factoring polynomials over local fields II.
        Algorithmic Number Theory, 9th International Symposium, ANTS-IX,
        Nancy, France, July 2010, LNCS 6197, 301-315, Springer Verlag 2010.

    AUTHORS:

    - Brian Sinclair and Sebastian Pauli (2012-02-22): initial version

    """
    if not Phi.is_monic():
        # Call non monic transform wrapper
        return pfactor_non_monic(Phi)
    else:
        # Phi is monic

        # Handle the situation that x is a factor of $\Phi(x)$
        if Phi.constant_coefficient() == 0:
            x_divides = True
            Phi = Phi >> 1
        else:
            x_divides = False

        if Phi == Phi.base_ring().one():
            return Factorization([(Phi.parent().gen(),1)] if x_divides else [])

        # Build an OM Tree for Phi
        tree = OM_tree(Phi)

        # If we only have one leaf, Phi is irreducible, so we do not lift it.
        if len(tree) == 1:
            return Factorization(([(Phi.parent().gen(),1)] if x_divides else []) +
                                 [(Phi,1)])
        # quo_rem is faster than single_factor_lift, so Phi = f*g is specially handled
        if len(tree) == 2:
            if tree[0].phi.degree() < tree[1].phi.degree():
                fact = single_factor_lift(tree[0])
                return Factorization(([(Phi.parent().gen(),1)] if x_divides else []) +
                                     [(fact,1),(Phi.quo_rem(fact)[0],1)])
            else:
                fact = single_factor_lift(tree[1])
                return Factorization(([(Phi.parent().gen(),1)] if x_divides else []) +
                                     [(fact,1),(Phi.quo_rem(fact)[0],1)])
        # Phi has more than two factors, so we lift them all
        return Factorization(([(Phi.parent().gen(),1)] if x_divides else []) +
                             [(single_factor_lift(frame),1) for frame in tree])


def OM_tree(Phi):
    r"""
    Return an tree of OM (Okutsu-Montes/Ore-Mac Lane) representations for Phi.

    INPUT:

    - ``Phi`` -- squarefree, monic padic polynomial with fixed precision
      coefficients

    OUTPUT:

    The leaves of the OM tree of ``Phi`` as a list of Frames.

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.factor.factoring import OM_tree
        sage: Phi = ZpFM(2,20,'terse')['x'](x^32+16)
        sage: OM_tree(Phi)
        [Frame with phi (1 + O(2^20))*x^16 + (0 + O(2^20))*x^15 + (0 + O(2^20))*x^14 + (0 + O(2^20))*x^13 + (0 + O(2^20))*x^12 + (0 + O(2^20))*x^11 + (1048572 + O(2^20))*x^10 + (0 + O(2^20))*x^9 + (1048572 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (1048572 + O(2^20))*x^5 + (4 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (8 + O(2^20))*x^2 + (0 + O(2^20))*x + (4 + O(2^20))]
    """
    from sage.misc.flatten import flatten

    def followsegment(next,Phi):
        """
        Returns next if it corresponds to an irreducible factor of $\Phi$ 
        and follows every branch if not.

        """
        # Handle the unlikely event that our approximation is actually a factor
        if next.is_first() == False and next.phi == next.prev_frame().phi:
            return [next]
        if next.phi_divides_Phi():
            return [next]+[[followsegment(fact.next_frame(fact.rhoexp+1),Phi)
                            for fact in seg.factors] for seg in next.polygon[1:]]
        # With E potentially increased, Check to see if E*F == deg(Phi)
        # (and thus Phi is irreducible)
        if next.E * next.F * next.polygon[0].Eplus == Phi.degree():
            return next
        # With F potentially increaded, Check to see if E*F == deg(Phi)
        # (and thus Phi is irreducible)
        if (next.E * next.polygon[0].Eplus * 
                next.F * next.polygon[0].factors[0].Fplus) == Phi.degree():
            return next
        # Check if we should begin Single Factor Lifting
        if sum([seg.length for seg in next.polygon]) == 1:
            return next
        return [[followsegment(fact.next_frame(fact.rhoexp+1),Phi)
                 for fact in seg.factors] for seg in next.polygon]

    # Construct and initialize the first frame (phi = x)
    next = Frame(Phi)
    next.seed(Phi.parent().gen())

    # With E potentially increased, Check to see if E*F == deg(Phi)
    # (and thus Phi is irreducible)
    if next.E * next.F * next.polygon[0].Eplus == Phi.degree():
        return [next]

    # With F potentially increaded, Check to see if E*F == deg(Phi)
    # (and thus Phi is irreducible)
    if (next.E * next.polygon[0].Eplus * 
            next.F * next.polygon[0].factors[0].Fplus) == Phi.degree():
        return [next]

    # Handle the special case wherein our initial approximation (phi = x) is a factor
    if next.phi_divides_Phi():
        tree = [next] + [[followsegment(fact.next_frame(fact.rhoexp+1),Phi)
                          for fact in seg.factors] for seg in next.polygon[1:]]

    # Construct the next level of the tree by following every factor of the
    # residual polynomials of every Newton polygon segment in our frame
    else:
        tree = [[followsegment(fact.next_frame(fact.rhoexp+1),Phi)
                 for fact in seg.factors] for seg in next.polygon]

    # tree contains the leaves of the tree of frames and each leaf corresponds
    # to an irreducible factor of Phi, so we flatten the list and start lifting
    return flatten(tree)


def jorder4(p):
    r"""
    Produce a particularly complicated example of polynomials for
    factorization over Zp(p).

    This method exists to test p-adic polynomial factorization.

    INPUT:

    - ``p`` -- a prime number

    OUTPUT:

    - A p-adic polynomial over Zp(``p``)

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.factor.factoring import jorder4
        sage: jorder4(3)
        (1 + O(3^20))*x^24 + (24 + O(3^20))*x^23 + (276 + O(3^20))*x^22 + (2048 + O(3^20))*x^21 + (11310 + O(3^20))*x^20 + (51180 + O(3^20))*x^19 + (201652 + O(3^20))*x^18 + (709092 + O(3^20))*x^17 + (2228787 + O(3^20))*x^16 + (6232484 + O(3^20))*x^15 + (15469950 + O(3^20))*x^14 + (34143276 + O(3^20))*x^13 + (67323664 + O(3^20))*x^12 + (119268300 + O(3^20))*x^11 + (190652502 + O(3^20))*x^10 + (275456480 + O(3^20))*x^9 + (359189415 + O(3^20))*x^8 + (420635664 + O(3^20))*x^7 + (438402286 + O(3^20))*x^6 + (400618284 + O(3^20))*x^5 + (313569267 + O(3^20))*x^4 + (203945072 + O(3^20))*x^3 + (105227142 + O(3^20))*x^2 + (38341248 + O(3^20))*x + (10912597 + O(3^20))

    Input must be prime::

        sage: jorder4(4)
        Traceback (most recent call last):
        ...
        ValueError: p must be prime

    """
    K = ZpFM(p,20,print_mode='terse')
    Kx = PolynomialRing(K,names='x')
    x = Kx.gen()
    f1 = (x+1)**3+p;
    f2 = f1**2+p**2*(x+1);
    f3 = f2**2+4*p**3*f1*(x+1)**2;
    f4 = f3**2+20*p**2*f3*f2*(x+1)**2+64*p**9*f1;
    return f4;

def single_factor_lift(frame,prec=2**20):
    r"""
    Lift a Frame to a factor of the polynomial it approximates.

    INPUT:

    - ``frame`` -- a Frame that is the leaf of an OM tree.

    OUTPUT:

    A factor of the polynomial referred to by the input frame

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.factor.factoring import single_factor_lift
        sage: from sage.rings.polynomial.padics.factor.frame import Frame
        sage: Kx.<x> = PolynomialRing(ZpFM(3,20,'terse'))
        sage: f = (x**3+x-1)*(x**2+1)
        sage: fr = Frame(f)
        sage: fr.seed(x)
        sage: fr = fr.polygon[0].factors[0].next_frame()
        sage: fact = single_factor_lift(fr) ; fact
        (1 + O(3^20))*x + (904752403 + O(3^20))
        sage: f % fact
        0

    REFERENCES:

        J. Guardia, E. Nart, S. Pauli. Single-factor lifting and
        factorization of polynomials over local fields.
        J. Symb. Comput. 47(11): 1318-1346 (2012)

    """
    def _reduce(poly,phi,d):
        """ returns poly mod phi and simplifies the denominator of poly """
        poly = phi.parent()(poly) % phi
        if d != 0:
            g = min([d] + [p.valuation() for p in poly])
            if g > 0:
                poly = poly.parent( [p >> g for p in poly] )
                d = d - g
        return poly,d
    if frame.phi_divides_Phi():
        return frame.phi

    LiftRing = ZpFM(frame.O.uniformizer(),2*frame.O.precision_cap())
    Lifty = PolynomialRing(LiftRing,names='y')
    Phi = Lifty(frame.Phi)
    phi = Lifty(frame.phi)
    a0,a1 = frame._phi_expansion_as_elts[0:2]

    Psi = frame.find_psi(-a1.valuation())
    A0 = Psi * a0
    A1 = Psi * a1

    Psi,Psiden = Psi.polynomial(True)
    Psi = Lifty(Psi)

    C1inv = frame.polygon[0].factors[0].lift(1/(A1.residue()))
    C1inv,C1invden = C1inv.polynomial(True)
    C1inv,C1invden = _reduce(C1inv,phi,C1invden)

    A0,A0den = A0.polynomial(True)
    A0,A0den = _reduce(Lifty(A0),phi,A0den)

    C,Cden = _reduce(frame.Ox(A0*C1inv),phi,A0den+C1invden)
    phi = (phi + C)

    h = 2
    oldphi = None
    while h < prec and phi != oldphi:
        oldphi = phi
        C1, C0 = Phi.quo_rem(phi)

        C0,C0den = _reduce((Psi*C0),phi,Psiden)
        C1,C1den = _reduce((Psi*C1),phi,Psiden)

        x1,x1den = _reduce((LiftRing(2)<<(C1den+C1invden))-C1*C1inv,phi,C1den+C1invden)
        C1inv,C1invden = _reduce(C1inv*x1,phi,C1invden+x1den)

        C,Cden = _reduce((C0*C1inv),phi,C0den+C1invden)

        phi = (phi + C)
        h = 2 * h
    return frame.Ox(phi)
