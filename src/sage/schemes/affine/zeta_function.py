r"""
Zeta function of hypersurfaces
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from time import time
from sys import stdout

from sage.arith.all import factorial, binomial, factor
from sage.combinat.integer_lists import IntegerListsLex
from sage.combinat.integer_vector import IntegerVectors
from sage.categories.fields import Fields
from sage.functions.log import log
from sage.functions.other import floor
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.all import Matrix
from sage.matrix.constructor import ones_matrix, identity_matrix, zero_matrix
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Zp
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.parent_gens import localvars


def zeta_function_of_polynomial(f, p=None, N=None, affine=False, verbose=False):
    r"""
    Return the zeta function of `f` up to precision `N`.

    The zeta function is the generating function for the number of points
    on the hypersurface defined by `f` over a finite field.

    The ambient space may be either the algebraic torus
    `(\overline{\GF_q}^\times)^n` (our default), or the
    affine space `(\overline{\GF_q})^n`.

    This implementation assumes that the finite field is of prime order.

    INPUT:

    - ``f'' -- polynomial defined over a finite field or
               ZZ if p is given (currently of prime size q = p)

    - ``p'' -- prime integer (default: ``None'') prime to compute over

    - ``N'' -- positive integer (default: ``None'') precision to
               calculate zeta function with. If ``None'' compute
               sufficient `N` to recover entire zeta function

    - ``affine'' -- boolean (default: ``False``) Flag to compute affine
                    or toric zeta function

    - ``verbose'' -- boolean (default: ``False``) debugging tool

    OUTPUT:

    - ``Z'' -- The zeta function of `f`, `Z(f, T)`

    EXAMPLES::

        sage: R.<x> = PolynomialRing(ZZ,1)
        sage: (1+x^3).zeta_function(7)
        -1/(T^3 - 3*T^2 + 3*T - 1)

        sage: R.<x,y> = PolynomialRing(FiniteField(7))
        sage: (x^3+x+1-y^2).zeta_function(7)
        (7*T^7 - 17*T^6 + 14*T^5 - 12*T^4 + 18*T^3 - 14*T^2 + 5*T - 1)/(7*T - 1)
    """
    R = PolynomialRing(ZZ, "T")
    T = R.gen()

    # ------------------------------------------------------------------
    #                             Input Checks
    # ------------------------------------------------------------------

    # Check that the polynomial is defined over Z or a finite field
    S = f.base_ring()
    if not (S == ZZ or S in Fields().Finite()):
        raise TypeError("the polynomial must be defined over ZZ or a finite field")

    if S.is_finite():
        if not S.is_prime_field():
            raise NotImplementedError("The field must be of prime order")

    # ------------------------------------------------------------------
    #                             Initialization
    # ------------------------------------------------------------------

    if verbose:
        print("Starting...\n")
        # Global variable holding start time of each step
        global start_time

    # Print function at start of each step
    def step_start(step):
        global start_time

        print("Running Step %d..." % step,)
        # Flush output
        stdout.flush()
        # Return time when the step started
        start_time = time()

    # Print function at the end of each step
    def step_end(step):
        global start_time

        print("complete.")

        total_time = time() - start_time
        if total_time < 0.1:
            print("(%.2fms)" % (total_time * 1000))
        elif total_time > 120:
            print("(%.2fm)" % (total_time / 60))
        else:
            print("(%.2fs)" % total_time)

    # Helper function that computes Z(G_m^n, T)
    # the zeta function of the algebraic torus
    def zeta_torus(p, n):
        Z = R.one()
        for i in range(n + 1):
            Z *= (1 - p**i * T)**(ZZ(n).binomial(i) * (-1)**(n + 1 - i))
        return Z

    # Helper function that computes Z(A^n, T)
    # the zeta function of affine space
    def zeta_affine(p, n):
        return R.one() / (1 - p**n * T)

    # ------------------------------------------------------------------
    #                                 Step 0
    # ------------------------------------------------------------------

    # Determine p
    # Use the characteristic if not given
    if p is None:
        p = S.characteristic()
        if p == 0:
            raise ValueError("p must be prime")
    # Otherwise make sure the given value is prime
    p = ZZ(p)
    if not p.is_prime():
        raise ValueError("p must be prime")

    # Determine the number of variables
    n = f.parent().ngens()

    # Form the Newton polytope of f - self.exponents()
    NP = Polyhedron(f.exponents())

    # Handle base cases based on dimension of Newton polytope
    n_tilde = NP.dimension()

    # Empty polytope (f = 0), return Z(G_m^n,T)
    if n_tilde == -1:
        if affine:
            return zeta_affine(p, n)
        else:
            return zeta_torus(p, n)

    # 0-dimensional polytope (f = ax^u, a,u != 0), return 1
    elif n_tilde == 0:
        if affine:
            raise NotImplementedError("Computing affine zeta functions of a"
                                      " single monomial not yet implemented")
        else:
            return PolynomialRing(ZZ, "T")(1)

    # Determine normalized volume to the Newton polytope
    v = factorial(n_tilde) * NP.affine_hull().volume()

    # Determine best guess for N if not given
    if N is None:
        bounds = []
        for i in range(v):
            bounds.append(ZZ(v - 1).binomial(i) * (p**(i * ZZ(n - 1) / 2)).n())
        N = log(max(bounds), p).floor() + 1

    # Update user with relevant information
    if verbose:
        print("+-------------------------")
        print("| f =", f)
        print("| p =", p)
        print("| n =", n)
        print("| v =", v)
        print("| N =", N)
        print("| NP =", NP)
        print("+-------------------------")
        print("")

    # Create polynomial ring with coefficients in ZZ/p^N ZZ,
    # in n+1 variables x = (x0, x1, ..., xn) [w = x0],
    # under lexicographical order
    Rx = PolynomialRing(Zmod(p**N), n + 1, "x", order="lex")
    x = Rx.gens()
    w = x[0]

    # Convert f to be defined over R with w=x0 the dummy variable
    with localvars(f.parent(), ["x%d" % i for i in range(1, n + 1)]):
        self = Rx(f)

    # ------------------------------------------------------------------
    #                                 Step 1
    # ------------------------------------------------------------------

    if verbose:
        step_start(1)

    # Construct the p-adic integers with fixed modulus N
    ZpN = Zp(p, N, "fixed-mod")

    # Take the teichmuller lift of the coefficients of f
    self = self.map_coefficients(ZpN.teichmuller)
    self = self.map_coefficients(Zmod(p**N))

    # Introduce dummy variable w = x0
    self = w * self

    # Calculate the generators of the Jacobian ideal
    Jacobian = [xi * self.derivative(xi) for xi in x]

    # Initialize J with empty matrix
    J = [Matrix(ZZ, sparse=True)]

    # Initialize RwD
    if affine:
        RwD = [[[] for i in range(n + 1)]]
    else:
        RwD = [[[Rx.one()]]]

    # Helper function for multi-index exponentiation:
    # x^u = x1^u1 * x2^u2 * ... * xn^un
    def mexp(x, u):
        return prod(map(pow, x, u))

    # Form Jd: RwD_{d-1}^{n+1} -> RwD_{n} for d = 1, ..., n+1
    for d in range(1, n + 2):

        integral_points = NP.dilation(d).integral_points()

        RwD.append([])
        span_set = []

        if affine:
            for i in range(n + 1):
                # Compute monomials in R[wD]_d^{S-i}
                RwD[d].append([w**d * mexp(x[1:], pt)
                               for pt in integral_points
                               if all((w**d * x[i] * mexp(x[1:], pt)).exponents()[0])])
                span_set.extend([mon * Jacobian[i] for mon in RwD[d - 1][i]])
        else:
            RwD[d].append([w**d * mexp(x[1:], pt) for pt in integral_points])
            span_set = [mon * fi for fi in Jacobian for mon in RwD[d - 1][0]]

        Jd = Matrix(ZZ, len(span_set), len(RwD[d][0]), sparse=True)
        for i in range(len(span_set)):
            for j in range(len(RwD[d][0])):
                Jd[i, j] = span_set[i].monomial_coefficient(RwD[d][0][j])

        J.append(Jd)

    # Initialize V
    V = []

    # Add to V the corresponding monomials of the nonpivots of
    # each Jd matrix for d = 1, ..., n. Then append to each
    # Jd a corresponding pivot row for step 3
    for d in range(1, n + 1):
        nonpivots = J[d].nonpivots()
        V.extend(RwD[d][0][i] for i in nonpivots)
        J[d] = J[d].stack(identity_matrix(len(RwD[d][0]))[nonpivots, :])

    # Compute the left inverse of each Jd for use in step 3
    for d in range(1, n + 2):
        pivot_rows = J[d].change_ring(Zmod(p)).pivot_rows()
        rows, cols = J[d].dimensions()
        if len(pivot_rows) < cols:
            raise ValueError("J%d has no left inverse" % d)
        JI = J[d][pivot_rows, :].change_ring(Zmod(p**N)).inverse()
        ZM = zero_matrix(Zmod(p**N), cols, rows)
        ZM[:, pivot_rows] = JI
        J[d] = ZM

    if verbose:
        step_end(1)

    # ------------------------------------------------------------------
    #                                 Step 2
    # ------------------------------------------------------------------

    if verbose:
        step_start(2)

    # List of splitting function coefficients; ell0 = 1
    ell = [ZZ.one()]

    # Use recurrence relation to calculate splitting
    # function coefficients i = 1, ..., floor((N+1)*p^2/(p-1))
    for i in range(1, (N + 1) * p**2 // (p - 1) + 1):
        if i < p:
            ell.append(ell[i - 1] / ZZ(i))
        else:
            ell.append((ell[i - 1] + ell[i - p]) / ZZ(i))

    # Convert coefficients to ThetaCoefficient object
    # for easier computation
    ell = [_ThetaCoefficient(a, p, N) for a in ell]

    # Exponent matrix U over ZZ
    U = Matrix(ZZ, self.exponents(as_ETuples=False))

    # Determine kernel and number of monomials = s
    s = U.nrows()
    ker = U.change_ring(Zmod(p)).kernel()

    # Compute alpha(v) for v in V and store it in list alpha_V
    alpha_V = []
    for vi in V:
        # Determine d and mu from exponents of v
        dmu = vector(ZZ, vi.exponents()[0])  # 0 since v is a monomial
        d = dmu[0]                          # exponent of x0

        # Compute -A- solution to (d,mu) = e*K  (mod p)
        sol = U.change_ring(Zmod(p)).solve_left(-dmu)

        # Compute alpha(v) from equation [?]
        alpha_v = 0
        for t in ker:
            k = (sol + t).change_ring(ZZ)

            piece1 = mexp(self.coefficients(), k) * mexp(x, (k * U + dmu) / ZZ(p))

            # Limits for possible integer vectors, e, based on p and N
            beta = ZZ(p**2 - p) / ZZ(p**2 - 3 * p + 1)
            limit = beta * (N + 1 - d / ZZ(p)) - sum(k) / ZZ(p)
            outer = [floor((N + 1) * ZZ(p) / ZZ(p - 1) - ki / ZZ(p))
                     for ki in k]

            piece2 = 0
            for e in IntegerListsLex(length=U.nrows(), min_sum=0,
                                     max_sum=limit.floor(), ceiling=outer):
                e = vector(ZZ, e)
                p_exponent = (sum(k) + d) // p + sum(e) - 1

                piece = prod([ell[ki + p * ei] for (ki, ei) in zip(k, e)])
                piece = piece.mul_p(p_exponent)
                piece = piece.as_int()
                piece *= Zmod(p)((-1)**(p_exponent + 1))
                piece *= mexp(x, e * U) * mexp(self.coefficients(), e)

                piece2 += piece

            alpha_v += piece1 * piece2

        alpha_V.append(alpha_v)

    if verbose:
        step_end(2)

    # ------------------------------------------------------------------
    #                                 Step 3
    # ------------------------------------------------------------------

    if verbose:
        step_start(3)

    # Triangulate NP into simplices t,
    # and compute U|_t and U|_t^-1 for each t
    Uts = []

    for t in NP.triangulate():
        t = Polyhedron([NP.Vrepresentation(i) for i in t])
        m = t.vertices_matrix().transpose()
        m = ones_matrix(n_tilde + 1, 1).augment(m)
        IM = identity_matrix(ZZ, n + 1)
        IM[m.pivot_rows(), :] = m
        Uts.append((IM, IM.inverse()))

    # Reduce each G in alpha_V and record transformation matrix A
    A = []
    for G in alpha_V:

        # Vector representing reduced G in cohomology
        Gbar = []

        # Iterate until G is a constant
        while not G.is_constant():

            # Determine exponent of lm(G)
            dmu = vector(ZZ, max(G.exponents()))
            d = dmu[0]

            # Decompose ulm if necessary with Uts
            if d > n + 1:
                for Ut, Uti in Uts:
                    y = dmu * Uti
                    if all(yi >= 0 for yi in y):
                        for z in IntegerVectors(d - (n + 1), n + 1,
                                                outer=map(floor, y)):
                            z = vector(ZZ, z)
                            um0 = (y - z) * Ut
                            if affine:
                                if all(umi > 0 for umi in um0):
                                    break
                            else:
                                break

                        um = dmu - um0
                        break
            else:
                um = vector(ZZ, n + 1)

            # Compute m in RwD_{n+1}
            m = mexp(x, um)

            # Determine what step we're on
            step = min(d, n + 1)

            # Gm = part of G that is "divisible" by m, i.e. the terms in m*RwD
            mRwD = set(m * mon for mon in RwD[step][0])
            Gm = sum(prod(g) for g in G if g[1] in mRwD)

            # xi = Gm/m as vector according to basis
            xi = vector(G.base_ring(), len(RwD[step][0]))
            for i, mon in enumerate(RwD[step][0]):
                Gm = (m * mon).parent()(Gm)
                xi[i] = Gm.monomial_coefficient(m * mon)

            # Compute inverse of xi
            # xii = (eta0, eta1, eta2, ..., etan, [basis part])
            xii = xi * J[step]

            # Extract eta = (eta0, eta1, ..., etan) from xii
            eta = []
            idx = 0
            for i in range(n + 1):
                s = len(RwD[step - 1][i * affine])
                etai = 0
                for j in range(s):
                    etai += xii[j + idx] * RwD[step - 1][i * affine][j]
                eta.append(etai)
                idx += s

            # Record basis pieces
            Gbar += reversed(list(xii[idx:]))

            # Reduce remaining part
            meta = Gm.parent()(m * eta[i])
            x_i = Gm.parent()(x[i])
            Gmbar = -sum([x_i * meta.derivative(x_i) for i in range(n + 1)])

            # Update G by removing Gm but adding reduced piece
            G = G - Gm + Gmbar

        # --------------------------------------- #
        # After G is a reduced to just a constant #
        # --------------------------------------- #

        # Constant MUST be zero
        assert G == 0, "Something is wrong, didn't reduce completely"

        # Insert 0's if we reduced to 0 before going through all the steps
        while len(Gbar) < len(V):
            Gbar.append(0)

        # Reverse Gbar so matches order of basis
        Gbar.reverse()

        # Add to A matrix
        A.append(Gbar)

    # Convert list A into a Sage matrix
    A = Matrix(A)

    if verbose:
        step_end(3)

    # ------------------------------------------------------------------
    #                                 Step 4
    # ------------------------------------------------------------------

    if verbose:
        step_start(4)

    # Function to convert coefficients
    # to lie in [-(p^N-1)/2, (p^N-1)/2]
    def normalize(a):
        if 2 * a > p**N - 1:
            return a - p**N
        else:
            return a

    # Compute the characteristic polynomial of A
    g = A.charpoly("T").reverse()

    # Normalize the coefficients and add (1-T) factor
    g = g.change_ring(ZZ)
    g = g.map_coefficients(normalize)
    g = g(p * T)

    # Multiply by (1 - T) factor if toric
    if not affine:
        g *= 1 - T

    # Compute L(wf) = g^((delta^(n-n_tilde))*(-1)^n)
    L = 1
    for i in range(n - n_tilde + 1):
        L *= g(p**i * T)**(binomial(n - n_tilde, i) * (-1)**i)
    L = L**((-1)**n)

    # Attempt to convert L to lie in ZZ(t) if possible
    try:
        L = R.fraction_field()(L)
    except TypeError:
        # Cannot convert to ZZ(t)
        print("WARNING: Cannot convert L(wf,T) to ZZ(t)")
        print("Very likely precision is too low")

    if affine:
        # Z(f, pT) = Z(A^n, T) * L(wf, T)
        ZpT = zeta_affine(p, n) * L
    else:
        # Z(f, pT) = Z(G_m^n, T) * L(wf, T)
        ZpT = zeta_torus(p, n) * L

    Z = ZpT(T / ZZ(p))

    # Attempt to convert Z to lie in ZZ(t) if possible
    try:
        Z = R.fraction_field()(Z)
    except TypeError:
        # Cannot convert to ZZ(t)
        print("WARNING: Cannot convert Z(f,T) to ZZ(T)")
        print("Very likely precision is too low")

    if verbose:
        step_end(4)

    if verbose:
        print("")
        print("Basis V = %s" % V)

        print("")
        print("A matrix =")
        print(A.change_ring(ZZ).apply_map(normalize))

        print("")
        print("g(T) =", g)
        print("     =", factor(g))

        if affine:
            print("")
            h = R(g(T / ZZ(p)))
            print("g(T/p) =", h)
            print("       =", factor(h))

        else:
            print("")
            h = R((g / (1 - T))(T / ZZ(p)))
            print("(g/(1-T))(T/p) =", h)
            print("               =", factor(h))

        print("")
        print("L(wf,T) =", L)
        print("        =", factor(L))

        print("")
        print("Chosen N = %d" % N)
        print("Minimum Necessary N =", )
        print(floor(log(max(map(abs, h.coefficients())), p) + 1))

    if verbose:
        print("\n...finished.")

    # Output Zeta function Z(f, T)
    return Z


class _ThetaCoefficient:
    """
    Helper class for :func:`zeta_function_of_polynomial`

    Handles splitting function coefficients and possible p's in denominator

    EXAMPLES::

        sage: from sage.schemes.affine.zeta_function import _ThetaCoefficient
        sage: a = _ThetaCoefficient(55/43, 5, 3); a
        (10, 0)
        sage: b = _ThetaCoefficient(23/5, 5, 3); b
        (23, 1)
        sage: a*b
        (105, 1)
        sage: a.as_int()
        10
    """
    def __init__(self, t, p, N):
        self.p = p
        self.N = N

        if isinstance(t, tuple):
            x, v = t
            self.x = x
            self.v = v
        else:
            self.v = max(-t.valuation(p), 0)
            self.x = Zmod(p**N)(t * p**self.v)

    def __str__(self):
        return str((self.x, self.v))

    def __repr__(self):
        return str(self)

    def __mul__(self, other):

        if self.p != other.p:
            raise ValueError("cannot combine coefficients"
                             "modulo different primes")
        if self.N != other.N:
            raise ValueError("cannot combine coefficients"
                             "with different precision")

        x = self.x * other.x
        v = self.v + other.v

        return _ThetaCoefficient((x, v), self.p, self.N)

    def __eq__(self, other):
        return self.x == other.x and self.v == other.v

    def __ne__(self, other):
        return not (self == other)

    def mul_p(self, c):
        if c <= 0:
            v = self.v - c
            x = self.x
        else:
            x = self.x * self.p**max(c - self.v, 0)
            v = min(c - self.v, 0)

        return _ThetaCoefficient((x, v), self.p, self.N)

    def as_int(self):
        if self.v > 0:
            raise ValueError("still has %d p's in denominator" % self.v)
        else:
            return self.x
