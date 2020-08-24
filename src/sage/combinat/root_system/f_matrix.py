"""
F-Matrix Factory for FusionRings
"""
# ****************************************************************************
#  Copyright (C) 2019 Daniel Bump <bump at match.stanford.edu>
#                     Guillermo Aboumrad <gh_willieab>
#                     Travis Scrimshaw <tcscrims at gmail.com>
#                     Galit Anikeeva <physicstravels@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.misc import inject_variable
from sage.matrix.constructor import matrix
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.ideal import Ideal
from sage.combinat.root_system.fusion_ring import FusionRing
import sage.graphs
from sage.graphs.generators.basic import EmptyGraph
from itertools import product

class FMatrix():
    r"""
    Return an F-Matrix factory for a FusionRing.

    INPUT:

    - ``FR`` -- a FusionRing.

    We only undertake to compute the F-matrix if the
    FusionRing is *multiplicity free* meaning that
    the Fusion coefficients `N^{ij}_k` are bounded
    by 1. For Cartan Types `X_r` and level `k`,
    the multiplicity-free cases are given by the
    following table.

+------------------------+----------+
| Cartan Type            | `k`      |
+========================+==========+
| `A_1`                  | any      |
+------------------------+----------+
| `A_r, r\geq 2`         | `\leq 2` |
+------------------------+----------+
| `B_r, r\geq 2`         | `\leq 2` |
+------------------------+----------+
| `C_2`                  | `\leq 2` |
+------------------------+----------+
| `C_r, r\geq 3`         | `\leq 1` |
+------------------------+----------+
| `D_r, r\geq 4`         | `\leq 2` |
+------------------------+----------+
| `G_2,F_4,E_r`          | `\leq 2` |
+------------------------+----------+

    Beyond this limitation, computation of the F-matrix
    can involve very large systems of equations. A
    rule of thumb is that this code can compute the
    F-matrix for systems with `\leq 4` primary fields,
    with the exception of `G_2` at level `2`.

    The FusionRing and its methods capture much
    of the structure of the underlying tensor category.
    But an important aspect that is not encoded in the
    fusion ring is the associator, which is a homomorphism 
    `(A\otimes B)\otimes C\to A\otimes(B\otimes C)`
    requires an additional tool, the F-matrix or 6j-symbol.
    To specify this, we fix a simple object `D`
    and represent the transformation

    .. MATH::

         \text{Hom}(D,(A\otimes B)\otimes C) \to \text{Hom}(D,A\otimes(B\otimes C))

    by a matrix `F^{ABC}_D`. This depends on a pair of
    additional simple objects `X` and `Y`. Indeed, we can
    get a basis for `\text{Hom}(D,(A\otimes B)\otimes C)`
    indexed by simple objects `X` in which the corresponding
    homomorphism factors through `X\otimes C`, and similarly
    `\text{Hom}(D,A\otimes(B\otimes C))` has a basis indexed
    by `Y`, in which the basis vector factors through `A\otimes Y`.

    See [TTWL2009]_ for an introduction to this topic,
    [EGNO2015]_ Section 4.9 for a precise mathematical
    definition and [Bond2007]_ Section 2.5 for a discussion 
    of how to compute the F-matrix.

    The F-matrix is only determined up to a gauge. This
    is a family of embeddings `C\to A\otimes B` for
    simple objects `A,B,C` such that `\text{Hom}(C,A\otimes B)`
    is nonzero. Changing the gauge changes the F-matrix though
    not in a very essential way. By varying the gauge it is
    possible to make the F-matrices unitary, or it is possible
    to make them cyclotomic. We choose the latter.

    Due to the large number of equations we may fail to find a 
    Groebner basis if there are too many variables.

    EXAMPLES::

        sage: I=FusionRing("E8",2,conjugate=True)
        sage: I.fusion_labels(["i0","p","s"],inject_variables=True)
        sage: f = FMatrix(I)
        Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
        sage: f.fmatrix(s,s,s,s)
        [fx10 fx11]
        [fx12 fx13]
        sage: f.f_from(s,s,s,s), f.f_to(s,s,s,s)
        ([i0, p], [i0, p])

    The F-matrix has not been computed at this stage, so
    the F-matrix `F^{sss}_s` is filled with variables
    ``fx10``, ``fx11``, ``fx12``, ``fx13``. These are subject
    pentagon equations reflecting the Pentagon relation
    in the monoidal category. They are also subject to
    Hexagon relations, related to the braiding in the
    category. The last two statments show that the possible
    values of `X` and `Y` when `A=B=C=D=s` are `i_0` and `p`.

    EXAMPLES::

        sage: f.pentagon()[1:3]                                                                                                                
        equations: 41
        [-fx0*fx1 + fx1, -fx1*fx2^2 + fx1]
        sage: f.hexagon()[1:3]
        equations: 14
        [-fx1*fx5 - fx2, (zeta128^32)*fx1*fx2 + (zeta128^32)*fx1]

You may solve these 41+14=55 equations to compute the F-matrix.

    EXAMPLES::

        sage: f.get_solution()
        Setting up hexagons and pentagons...
        equations: 14
        equations: 41
        Finding a Groebner basis...
        Solving...
        Fixing the gauge...
        adding equation... fx1 - 1
        adding equation... fx3 - 1
        Done!
        {(s, s, s, s, i0, i0): -1/2*zeta128^48 + 1/2*zeta128^16,
         (s, s, s, s, i0, p): -1/2*zeta128^48 + 1/2*zeta128^16,
         (s, s, s, s, p, i0): -1/2*zeta128^48 + 1/2*zeta128^16,
         (s, s, s, s, p, p): 1/2*zeta128^48 - 1/2*zeta128^16,
         (s, s, p, i0, p, s): 1,
         (s, s, p, p, i0, s): 1,
         (s, p, s, i0, s, s): 1,
         (s, p, s, p, s, s): -1,
         (s, p, p, s, s, i0): 1,
         (p, s, s, i0, s, p): 1,
         (p, s, s, p, s, i0): 1,
         (p, s, p, s, s, s): -1,
         (p, p, s, s, i0, s): 1,
         (p, p, p, p, i0, i0): 1}
        sage: f.fmatrix(s,s,s,s)
        [-1/2*zeta128^48 + 1/2*zeta128^16 -1/2*zeta128^48 + 1/2*zeta128^16]
        [-1/2*zeta128^48 + 1/2*zeta128^16  1/2*zeta128^48 - 1/2*zeta128^16]
        sage: f.fmat(s,s,s,s,p,p)
        1/2*zeta128^48 - 1/2*zeta128^16

    """
    def __init__(self, fusion_ring, fusion_label="f", var_prefix='fx'):
        self.FR = fusion_ring
        if self.FR._fusion_labels is None:
            self.FR.fusion_labels(fusion_label, inject_variables=True)
            #Set up F-symbols entry by entry
        n_vars = self.findcases()
        self._poly_ring = PolynomialRing(self.FR.field(),n_vars,var_prefix)
        self._poly_ring.inject_variables()
        self._var_to_sextuple, self._fvars = self.findcases(output=True)

        #Initialize set of defining equations
        self.ideal_basis = set()

        #Initialize empty set of solved F-symbols
        self.solved = set()

    def __repr__(self):
        return "F-Matrix factory for %s"%self.FR

    def remaining_vars(self):
        """
        Return a list of unknown F-symbols (reflects current stage of computation)
        """
        return [var for var in self._poly_ring.gens() if var not in self.solved]

    def f_from(self,a,b,c,d):
        r"""
        Return the possible `x` such that there are morphisms
        `d\to x\otimes c\to (a\otimes b)\otimes c`.

        INPUT:

        - ``a,b,c,d`` -- basis elements of the FusionRing.
        """

        return [x for x in self.FR.basis() if self.FR.Nk_ij(a,b,x) != 0 and self.FR.Nk_ij(x,c,d) != 0]

    def f_to(self,a,b,c,d):
        r"""
        Return the possible `y` such that there are morphisms
        `d\to a\otimes y\to a\otimes(b\otimes c)`.

        INPUT:

        - ``a,b,c,d`` -- basis elements of the FusionRing.
        """

        return [y for y in self.FR.basis() if self.FR.Nk_ij(b,c,y) != 0 and self.FR.Nk_ij(a,y,d) != 0]

    def fix_gauge(self):
        """
        Fix the gauge by forcing F-symbols not already fixed to equal 1.
        This method should be used AFTER adding hex and pentagon eqns to ideal_basis
        """
        while len(self.solved) < len(self._poly_ring.gens()):
            #Get a variable that has not been fixed
            #In ascending index order, for consistent results
            for var in self._poly_ring.gens():
                if var not in self.solved:
                    break

            #Fix var = 1, substitute, and solve equations
            self.ideal_basis.add(var-1)
            print("adding equation...", var-1)
            self.ideal_basis = set(Ideal(list(self.ideal_basis)).groebner_basis())
            self.substitute_known_values()
            self.update_equations()

    def substitute_known_values(self, eqns=None):
        if not eqns:
            eqns = self.ideal_basis

        for eq in eqns:
            #Ensure polynomial is degree 1 in a single variable
            if sum(eq.degrees()) == 1 and eq.monomials()[0] not in self.solved:
                var = eq.monomials()[0]
                self._fvars[self._var_to_sextuple[var]] = -eq.constant_coefficient()
                #Add variable to set of known values
                self.solved.add(var)

    def update_equations(self):
        """
        Update ideal_basis equations by plugging in known values
        """
        special_values = { known : self._fvars[self._var_to_sextuple[known]] for known in self.solved }
        self.ideal_basis = set(eq.specialization(special_values) for eq in self.ideal_basis)
        self.ideal_basis.discard(0)

    def get_solution(self, equations=None, factor=False):
        """
        Solve the the hexagon and pentagon relations to evaluate the F-matrix.

        INPUT:

        - ``equations`` -- (optional) a set of equations to be
          solved. Defaults to the hexagon and pentagon equations.
        - ``factor`` -- (default: ``False``). Set true to use
          the sreduce method to simplify the hexagon and pentagon
          equations before solving them.

        """
        if equations is None:
            print("Setting up hexagons and pentagons...")
            equations = self.hexagon(factor=factor)+self.pentagon(factor=factor)
        print("Finding a Groebner basis...")
        self.ideal_basis = set(Ideal(equations).groebner_basis())
        print("Solving...")
        self.substitute_known_values()
        print("Fixing the gauge...")
        self.fix_gauge()
        print("Done!")
        return self._fvars

    def add_equations(self, eqns):
        #TODO: consider replacing set union by ideal intersection. (study computational cost)
        self.ideal_basis = self.ideal_basis.union(set(eqns))

    def clear_equations(self):
        """
        Clear the set of equations to be solved.
        """
        self.ideal_basis = set()

    def clear_vars(self):
        """
        Clear the set of variables.
        """
        self._fvars = { self._var_to_sextuple[key] : key for key in self._var_to_sextuple }

    def fmat(self, a, b, c, d, x, y, data=True):
        """
        Return the F-Matrix coefficient `(F^{a,b,c}_d)_{x,y}`
        """
        #Determine if fusion tree is admissible
        admissible = self.FR.Nk_ij(a,b,x) * self.FR.Nk_ij(x,c,d) * self.FR.Nk_ij(b,c,y) * self.FR.Nk_ij(a,y,d)

        if admissible == 0:
            return 0

        #Some known zero F-symbols
        if a == self.FR.one():
            if x == b and y == d:
                return 1
            else:
                return 0
        if b == self.FR.one():
            if x == a and y == c:
                return 1
            else:
                return 0
        if c == self.FR.one():
            if x == d and y == b:
                return 1
            else:
                return 0
        if data:
            return self._fvars.get((a,b,c,d,x,y),0)
        else:
            return (a,b,c,d,x,y)
        
    def fmatrix(self,a,b,c,d):
        """
        Return the F-Matrix `F^{a,b,c}_d`.
        """

        X = self.f_from(a,b,c,d)
        Y = self.f_to(a,b,c,d)
        return matrix([[self.fmat(a,b,c,d,x,y) for y in Y] for x in X])
    
    def findcases(self,output=False):
        """
        Find the unknown F-matrix entries. If run with output=True,
        this returns two dictionaries; otherwise it just returns the
        number of unknown values.
        """
        i = 0
        if output:
            idx_map = dict()
            ret = dict()
        for (a,b,c,d) in list(product(self.FR.basis(), repeat=4)):
            for x in self.f_from(a, b, c, d):
                for y in self.f_to(a, b, c, d):
                    fm = self.fmat(a, b, c, d, x, y, data=False)
                    if fm is not None and fm not in [0,1]:
                        if output:
                            v = self._poly_ring.gens()[i]
                            ret[(a,b,c,d,x,y)] = v
                            idx_map[v] = (a, b, c, d, x, y)
                        i += 1
        if output:
            return idx_map, ret
        else:
            return i

    def sreduce(self, expr, nonzeros=None):
        """
        From an equation, discard the leading coefficient and
        factors of degree one.
        """
        if nonzeros is None:
            nonzeros = self._poly_ring.gens()
        ret = 1
        for (a,e) in expr.factor()._Factorization__x:
            if a not in nonzeros:
                ret *= a**e
        return ret

    def feq(self, a, b, c, d, e, f, g, k, l):
        """
        Return True if the Pentagon axiom (Bonderson (2.77)) is satisfied.
        """
        lhs = self.fmat(f,c,d,e,g,l)*self.fmat(a,b,l,e,f,k)
        rhs = sum(self.fmat(a,b,c,g,f,h)*self.fmat(a,h,d,e,g,k)*self.fmat(b,c,d,k,h,l) for h in self.FR.basis())
        return lhs - rhs

    def req(self, a, b, c, d, e, g):
        """
        Return True if the Hexagon axiom (Bonderson (2.78)) is satisfied.
        """
        lhs = self.FR.r_matrix(a,c,e)*self.fmat(a,c,b,d,e,g)*self.FR.r_matrix(b,c,g)
        rhs = sum(self.fmat(c,a,b,d,e,f)*self.FR.r_matrix(f,c,d)*self.fmat(a,b,c,d,f,g) for f in self.FR.basis())
        return lhs-rhs

    def hexagon(self, verbose=False, output=True, factor=False):
        """
        Return generators of the ideal of solutions to the Hexagon equations.

        INPUT:

        - ``verbose`` -- (optional) set True for verbose. Default False
        - ``output`` -- (optional) set True to output a set of equations. Default True
        - ``factor`` -- (optional) set False for sreduce simplified equations.

        EXAMPLES::

            sage: f = FMatrix(FusionRing("A1",2),fusion_label="a")
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
            sage: f.hexagon()[-3:]
            equations: 14
            [(zeta32^8)*fx11*fx12 + (zeta32^8)*fx12, -fx8*fx12 - fx11, -fx13^2 + fx13]
            sage: f.hexagon(factor=True)[-3:]
            equations: 14
            [fx11 + 1, fx8*fx12 + fx11, fx13 - 1]

        """
        if output:
            ret = []
        for (a,b,c,d,e,g) in list(product(self.FR.basis(), repeat=6)):
            rd = self.req(a,b,c,d,e,g)
            if rd != 0:
                if factor:
                    rd = self.sreduce(rd)
                if output:
                    ret.append(rd)
                if verbose:
                    print ("%s,%s,%s,%s,%s,%s : %s"%(a,b,c,d,e,g,rd.factor()))
        print ("equations: %s"%len(ret))
        if output:
            return ret

    def pentagon(self, verbose=False, output=True, factor=None):
        """
        Return generators of the ideal of Pentagon equations.

        - ``verbose`` -- (optional) set True for verbose. Default False
        - ``output`` -- (optional) set True to output a set of equations. Default True
        - ``factor`` -- (optional) set False for sreduce simplified equations.

        In contrast with the hexagon equations, where setting ``factor`` True
        is a big improvement, for the pentagon equations this option produces
        little or no simplification.

        EXAMPLES::

            sage: p = FMatrix(FusionRing("A2",1),fusion_label="c")
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7
            sage: p.pentagon()[-3:]                                               
            equations: 16
            [-fx5*fx6 + fx1, -fx4*fx6*fx7 + fx2, -fx5*fx7^2 + fx3*fx6]
            sage: p.pentagon(factor=True)[-3:]
            equations: 16
            [fx5*fx6 - fx1, fx4*fx6*fx7 - fx2, fx5*fx7^2 - fx3*fx6]

        """
        if output:
            ret = []
        for (a,b,c,d,e,f,g,k,l) in list(product(self.FR.basis(), repeat=9)):
            pd = self.feq(a,b,c,d,e,f,g,k,l)
            if pd != 0:
                if factor:
                    pd = self.sreduce(pd)
                if output:
                    ret.append(pd)
                if verbose:
                    print ("%s,%s,%s,%s,%s,%s,%s,%s,%s : %s"%(a,b,c,d,e,f,g,k,l,pd))
        print ("equations: %s"%len(ret))
        if output:
            return ret

    def equation_graph(self, equations):
        """
        Return a graph whose vertices are variables
        and whose edges correspond to equations
        relating two variables.
        """
        G = sage.graphs.generators.basic.EmptyGraph()
        for e in equations:
            s = [v for v in e.variables()]
            for x in s:
                for y in s:
                    if y!=x:
                        G.add_edge(x,y)
        return(G)
