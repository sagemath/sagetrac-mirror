from sage.structure.element import ModuleElement
from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from sage.rings.integer_ring import ZZ
from manin_map import ManinMap
import operator
from sage.misc.cachefunc import cached_method
from sage.rings.padics.factory import Qp
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.arith import next_prime
from sage.rings.infinity import infinity
from sage.misc.misc import verbose

from sage.categories.action import Action

from fund_domain import M2ZSpace, M2Z, Id
minusproj = M2Z([1,0,0,-1])

class PSModSymAction(Action):
    def __init__(self, actor, MSspace):
        Action.__init__(self, actor, MSspace, False, operator.mul)

    def _call_(self, sym, g):
        return sym.__class__(sym._map * g, sym.parent(), construct=True)

class PSModularSymbolElement(ModuleElement):
    def __init__(self, map_data, parent, construct=False):
        ModuleElement.__init__(self, parent)
        if construct:
            self._map = map_data
        else:
            self._map = ManinMap(parent._coefficients, parent._source, map_data)

    def _repr_(self):
        r"""
        Returns the print representation.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi._repr_()
            'Modular symbol with values in Sym^0 Q^2'
        """
        return "Modular symbol with values in %s"%(self.parent().coefficient_module())

    def dict(self):
        """
        Returns dictionary on self, where keys are generators and values are the corresponding values of self on generators

        EXAMPLES::
        
            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.dict()
            {[1 0]
            [0 1]: -1/5, [ 0 -1]
            [ 1  3]: 3/2, [-1 -1]
            [ 3  2]: -1/2}
            
    
        """
        D = {}
        for g in self.parent().source().gens():
            D[g] = self._map[g]
        return D

    def weight(self):
        """
        Returns the weight of this Pollack-Stevens modular symbol.

        This is `k-2`, where `k` is the usual notion of weight for modular
        forms!

        EXAMPLES::
            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)  
            sage: phi.weight()
            0

        """
        return self.parent().weight()

    def values(self):
        """
        Returns the values of self on our chosen generators (generators are listed in self.dict().keys())

        EXAMPLES::
        
             sage: E = EllipticCurve('11a')
             sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
             sage: phi = ps_modsym_from_elliptic_curve(E)
             sage: phi.values()
             [-1/5, 3/2, -1/2]
             sage: phi.dict().keys()
             [
             [1 0]  [ 0 -1]  [-1 -1]
             [0 1], [ 1  3], [ 3  2]
             ]
             sage: phi.values() == phi.dict().values()
             True
        
        """
        return [self._map[g] for g in self.parent().source().gens()]

    def _normalize(self):
        for val in self._map:
            val.normalize()
        return self

    def __cmp__(self, other):
        """
        Checks if self == other

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi == phi
            True
            sage: phi == 2*phi
            False
            sage: psi = ps_modsym_from_elliptic_curve(EllipticCurve('37a'))
            sage: psi == phi
            False
        
        """
        gens = self.parent().source().gens()
        for g in gens:
            c = cmp(self._map[g], other._map[g])
            if c: return c
        return 0

    def _add_(self, right):
        """
        Returns self + right

        EXAMPLES::
        
            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi + phi
            Modular symbol with values in Sym^0 Q^2
            sage: (phi + phi).values()
            [-2/5, 3, -1]
        """
        return self.__class__(self._map + right._map, self.parent(), construct=True)

    def _lmul_(self, right):
        """
        Returns self * right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: 2*phi
            Modular symbol with values in Sym^0 Q^2
            sage: (2*phi).values()
            [-2/5, 3, -1]
        """
        return self.__class__(self._map * right, self.parent(), construct=True)

    def _rmul_(self, right):
        """
        Returns self * right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi*2
            Modular symbol with values in Sym^0 Q^2
            sage: (phi*2).values()
            [-2/5, 3, -1]
        """
        return self.__class__(self._map * right, self.parent(), construct=True)

    def _sub_(self, right):
        """
        Returns self - right

        EXAMPLES:;

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi - phi
            Modular symbol with values in Sym^0 Q^2
            sage: (phi - phi).values()
            [0, 0, 0]
        """
        return self.__class__(self._map - right._map, self.parent(), construct=True)

    def plus_part(self):
        r"""
        Returns the plus part of self -- i.e. self + self | [1,0,0,-1].

        Note that we haven't divided by 2.  Is this a problem?

        OUTPUT:

        - self + self | [1,0,0,-1]

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: (phi.plus_part()+phi.minus_part()) == 2 * phi
            True
        """
        return self * minusproj + self

    def minus_part(self):
        r"""
        Returns the minus part of self -- i.e. self - self | [1,0,0,-1]

        Note that we haven't divided by 2.  Is this a problem?

        OUTPUT:

        - self - self | [1,0,0,-1]

        EXAMPLES::

            sage: E = EllipticCurve('11a')  
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: (phi.plus_part()+phi.minus_part()) == phi * 2
            True
        """
        return self - self * minusproj

    def hecke(self, ell, algorithm="prep"):
        r"""
        Returns self | `T_{\ell}` by making use of the precomputations in
        self.prep_hecke()

        INPUT:

        - ``ell`` -- a prime

        - ``algorithm`` -- a string, either 'prep' (default) or
          'naive'

        OUTPUT:

        - The image of this element under the hecke operator
          `T_{\ell}`

        ALGORITHMS:

        - If ``algorithm == 'prep'``, precomputes a list of matrices
          that only depend on the level, then uses them to speed up
          the action.

        - If ``algorithm == 'naive'``, just acts by the matrices
          defining the Hecke operator.  That is, it computes
          sum_a self | [1,a,0,ell] + self | [ell,0,0,1],
          the last term occurring only if the level is prime to ell.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E); phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi.hecke(2) == phi * E.ap(2)
            True
            sage: phi.hecke(3) == phi * E.ap(3)
            True
            sage: phi.hecke(5) == phi * E.ap(5)
            True
            sage: phi.hecke(101) == phi * E.ap(101)
            True

            sage: all([phi.hecke(p, algorithm='naive') == phi * E.ap(p) for p in [2,3,5,101]])
            True
        """
        return self.__class__(self._map.hecke(ell, algorithm), self.parent(), construct=True)

    def valuation(self, p):
        r"""
        Returns the valuation of self at `p`.

        Here the valuation if the exponent of the largest power of `p`
        which divides all of the coefficients of all values of self.

        INPUT:

        - ``p`` - prime

        OUTPUT:

        - The valuation of self at `p`

        EXAMPLES::

           sage: E = EllipticCurve('11a')
           sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
           sage: phi = ps_modsym_from_elliptic_curve(E)
           sage: phi.values()
           [-1/5, 3/2, -1/2]
           sage: phi.valuation(2)
           -1
           sage: phi.valuation(3)
           0
           sage: phi.valuation(5)
           -1
           sage: phi.valuation(7)
           0
        """
        return min([val.valuation(p) for val in self._map])

    def diagonal_valuation(self, p):
        return min([val.diagonal_valuation(p) for val in self._map])

    @cached_method
    def is_Tq_eigensymbol(self,q,p=None,M=None):
        r"""
        Determines if self is an eigenvector for `T_q` modulo `p^M`

        INPUT:

        - ``q`` -- prime of the Hecke operator
        - ``p`` -- prime we are working modulo
        - ``M`` -- degree of accuracy of approximation

        OUTPUT:

        - True/False

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi_ord = phi.p_stabilize(p = 3, ap = E.ap(3), M = 10, ordinary = True)
            sage: phi_ord.is_Tq_eigensymbol(2,3,10)
            True
            sage: phi_ord.is_Tq_eigensymbol(2,3,100)
            True
            sage: phi_ord.is_Tq_eigensymbol(2,3,1000)
            True
            sage: phi_ord.is_Tq_eigensymbol(3,3,10)
            True
            sage: phi_ord.is_Tq_eigensymbol(3,3,100)
            False
        """
        try:
            aq = self.Tq_eigenvalue(q, p, M)
            return True
        except ValueError:
            return False


    # what happens if a cached method raises an error?  Is it recomputed each time?
    @cached_method
    def Tq_eigenvalue(self, q, p=None, M=None, check=True):
        r"""
        Eigenvalue of `T_q` modulo `p^M`

        INPUT:

        - ``q`` -- prime of the Hecke operator
        - ``p`` -- prime we are working modulo
        - ``M`` -- degree of accuracy of approximation

        OUTPUT:

        - Constant `c` such that `self|T_q - c * self` has valuation greater than
          or equal to `M` (if it exists), otherwise raises ValueError

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi_ord = phi.p_stabilize(p = 3, ap = E.ap(3), M = 10, ordinary = True)
            sage: phi_ord.Tq_eigenvalue(2,3,10)
            -2
            sage: phi_ord.Tq_eigenvalue(2,3,100)
            -2
            sage: phi_ord.Tq_eigenvalue(2,3,1000)
            -2


            sage: phi_ord.Tq_eigenvalue(3,3,10)
            -95227/47611
            sage: phi_ord.Tq_eigenvalue(3,3,100)
            Traceback (most recent call last):
            ...
            ValueError: not a scalar multiple
        """
        f = self.hecke(q)
        gens = self.parent().source().gens()
        g = gens[0]
        if p is None:
            p = self.parent().prime()
        aq = self._map[g].find_scalar(f._map[g], p, M, check)
        if check:
            if p is None or M is None:
                for g in gens[1:]:
                    if f._map[g] != aq * self._map[g]:
                        raise ValueError("not a scalar multiple")
            elif (f - aq * self).valuation(p) < M:
                raise ValueError("not a scalar multiple")
        return aq

class PSModularSymbolElement_symk(PSModularSymbolElement):
    def p_stabilize(self, p=None, M=None, alpha=None, ap=None, new_base_ring=None, ordinary=True, check=True):
        if check:
            pp = self.parent().prime()
            ppp = (alpha is not None) and alpha.parent().hasattr('prime') and alpha.parent().prime()
            p = ZZ(p) or pp or ppp
            if not p:
                raise ValueError("you must specify a prime")
            if (pp and p != pp) or (ppp and p != ppp):
                raise ValueError("inconsistent prime")
            if M is None:
                M = self.parent().precision_cap()
        k = self.parent().weight()
        if alpha is None:
            if ap is None:
                ap = self.Tq_eigenvalue(p, check=check)
            if check and ap % p == 0:
                raise ValueError("p is not ordinary")
            if new_base_ring is None:
                disc = ap**2 - 4*p**(k+1)
                if M is None:
                    if disc.is_square():
                        new_base_ring = disc.parent()
                    else:
                        poly = PolynomialRing(disc.parent(), 'x')([-disc, 0, 1])
                        new_base_ring = disc.parent().extension(poly, 'a')
                elif p == 2:
                    new_base_ring = Qp(2, M+1)
                else:
                    new_base_ring = Qp(p, M)
            sdisc = new_base_ring(disc).sqrt()
            v0 = (new_base_ring(ap) + sdisc) / 2
            v1 = (new_base_ring(ap) - sdisc) / 2
            if v0.valuation(p) > 0:
                v0, v1 = v1, v0
            if ordinary:
                alpha = v0
            else:
                alpha = v1
        else:
            if new_base_ring is None:
                new_base_ring = alpha.parent()
            if check:
                if ap is None:
                    ap = alpha + p**(k+1)/alpha
                elif alpha**2 - ap * alpha + p**(k+1) != 0:
                    raise ValueError("alpha must be a root of x^2 - a_p*x + p^(k+1)")
                if self.hecke(p) != ap * self:
                    raise ValueError("alpha must be a root of x^2 - a_p*x + p^(k+1)")
        V = self.parent()._p_stabilize_parent_space(p, new_base_ring)
        return self.__class__(self._map.p_stabilize(p, alpha, V), V, construct=True)

    def completions(self, p, M):
        r"""
        If `K` is the base_ring of self, this function takes all maps
        `K-->Q_p` and applies them to self return a vector of
        <modular symbol,map: `K-->Q_p`> as map varies over all such maps.

        NOTE: This only returns all completions when `p` splits completely in `K`

        INPUT:

        - ``p`` -- prime
        - ``M`` -- precision

        OUTPUT:

        - A vector of <modular symbol,map: `K-->Q_p`> as map varies over all such maps

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
            sage: D = ModularSymbols(67,2,1).cuspidal_submodule().new_subspace().decomposition()[1]
            sage: f = ps_modsym_from_simple_modsym_space(D)
            sage: f.completions(41,10)
            [(Modular symbol with values in Space of 41-adic distributions with k=0 action and precision cap 1, Ring morphism:
              From: Number Field in alpha with defining polynomial x^2 + 3*x + 1
              To:   41-adic Field with capped relative precision 10
              Defn: alpha |--> 5 + 22*41 + 19*41^2 + 10*41^3 + 28*41^4 + 22*41^5 + 9*41^6 + 25*41^7 + 40*41^8 + 8*41^9 + O(41^10)), (Modular symbol with values in Space of 41-adic distributions with k=0 action and precision cap 1, Ring morphism:
              From: Number Field in alpha with defining polynomial x^2 + 3*x + 1
              To:   41-adic Field with capped relative precision 10
              Defn: alpha |--> 33 + 18*41 + 21*41^2 + 30*41^3 + 12*41^4 + 18*41^5 + 31*41^6 + 15*41^7 + 32*41^9 + O(41^10))]
        
        """
        K = self.base_ring()
        f = K.defining_polynomial()
        R = Qp(p,M+10)['x']
        x = R.gen()
        v = R(f).roots()
        if len(v) == 0:
            raise ValueError, "No coercion possible -- no prime over p has degree 1"
        else:
            roots = [r[0] for r in v]
            ans = []
            V = self.parent().change_ring(Qp(p, M))
            Dist = V.coefficient_module()
            for r in roots:
                psi = K.hom([r],Qp(p,M))
                embedded_sym = self.__class__(self._map.apply(psi, codomain=Dist, to_moments=True), V, construct=True)
                ans.append((embedded_sym,psi))
            return ans

    def lift(self, p=None, M=None, new_base_ring=None, algorithm = None, eigensymbol = None):
        r"""
        """
        if p is None:
            p = self.parent().prime()
            if p is None:
                raise ValueError("must specify a prime")
        elif self.parent().prime() is not None and p != self.parent().prime():
            raise ValueError("inconsistent prime")
        if M is None:
            M = self.parent().precision_cap() + 1
        if new_base_ring is None:
            if isinstance(self.parent().base_ring(), pAdicGeneric):
                new_base_ring = self.parent().base_ring()
            else:
                # should eventually be a completion
                new_base_ring = Qp(p, M)
        if algorithm is None:
            raise NotImplementedError
        if algorithm == 'stevens':
            if eigensymbol:
                return self._lift_to_OMS_eigen(p, M, new_base_ring)
            else:
                return self._lift_to_OMS(p, M, new_base_ring)
        else:
            return self._lift_greenberg(p, M, new_base_ring)

    def _lift_to_OMS(self, p, M, new_base_ring):
        """
        Returns a (`p`-adic) overconvergent modular symbol with `M` moments which lifts self up to an Eisenstein error

        Here the Eisenstein error is a symbol whose system of Hecke eigenvalues equals `ell+1` for `T_ell` when `ell`
        does not divide `Np` and 1 for `U_q` when `q` divides `Np`.

        INPUT:
            - ``p`` -- prime
            - ``M`` -- integer equal to the number of moments
            - ``new_base_ring`` -- new base ring

        OUTPUT: 

        An overconvergent modular symbol whose specialization equals self up to some Eisenstein error.

        EXAMPLES:

        """
        
        D = {}
        manin = self.parent().source()
        MSS = self.parent()._lift_parent_space(p, M, new_base_ring)
        half = ZZ(1) / ZZ(2)
        for g in manin.gens()[1:]:
            twotor = g in manin.reps_with_two_torsion
            threetor = g in manin.reps_with_three_torsion
            if twotor:
                # See [PS] section 4.1
                gam = manin.two_torsion[g]
                mu = self._map[g].lift(p, M, new_base_ring)
                D[g] = (mu * gam - mu) * half
            elif threetor:
                # See [PS] section 4.1
                gam = manin.three_torsion[g]
                mu = self._map[g].lift(p, M, new_base_ring)
                D[g] = (2 * mu - mu * gam - mu * (gam**2)) * half
            else:
                # no two or three torsion
                D[g] = self._map[g].lift(p, M, new_base_ring)

        t = self.parent().coefficient_module().lift(p, M, new_base_ring).zero_element()
        for h in manin[2:]:
            R = manin.relations(h)
            if len(R) == 1:
                c, A, g = R[0]
                if c == 1:
                    t += self._map[h].lift(p, M, new_base_ring)
                elif A is not Id:
                    # rules out extra three torsion terms
                    t += c * self._map[g].lift(p, M, new_base_ring) * A
        D[manin.gen(0)] = t.solve_diff_eqn()
        return MSS(D)

    def _lift_to_OMS_eigen(self, p, M, new_base_ring):
        r"""
        Returns Hecke-eigensymbol OMS lifting self -- self must be a  
        `p`-ordinary eigensymbol

        INPUT:
            - ``p`` -- prime
            - ``M`` -- integer equal to the number of moments  
            - ``new_base_ring`` -- new base ring

        OUTPUT:

        EXAMPLES:
        
        
        """
        
        ap = self.Tq_eigenvalue(p, p, M)
        apinv = ~ap
        k = self.parent().weight()
        Phi = self._lift_to_OMS(p, M, new_base_ring)
        s = - Phi.diagonal_valuation(p)
        if s > 0:
            verbose("scaling by %s^%s"%(p, s))
            Phi = p**s * Phi
            need_unscaling = True
        else:
            need_unscaling = False
        Phi._normalize()
        verbose("Applying Hecke")
        Phi = apinv * Phi.hecke(p)
        verbose("Killing eisenstein part")
        if ap % p**M != 1:
            Phi = 1 / (1 - ap) * (Phi - Phi.hecke(p))
            e = (1 - ap).valuation(p)
            if e > 0:
                verbose("change precision to %s"%(M - e))
                Phi = Phi.reduce_precision(M - e)
        else:
            q = ZZ(2)
            aq = self.Tq_eigenvalue(q, p, M)
            while q != p and (aq - q**(k+1) - 1) % p**M == 0:
                q = next_prime(q)
                aq = self.Tq_eigenvalue(q, p, M)
            Phi = (q**(k+1) + 1 - aq) * ((q**(k+1) + 1) * Phi - Phi.hecke(q))
            e = (q**(k+1) + 1 - aq).valuation(p)
            if e > 0:
                verbose("change precision to %s"%(M - e))
                Phi = Phi.reduce_precision(M - e)
        verbose("Iterating U_p")
        Psi = apinv * Phi.hecke(p)
        err = (Psi - Phi).diagonal_valuation(p)
        Phi = Psi
        while err < M:
            if need_unscaling and Phi.valuation(p) >= s:
                Phi *= (1 / p**s)
                # Can't we get this to better precision....
                Phi = Phi.reduce_precision(Phi.precision_absolute() - s)._normalize()
                verbose("unscaling by %s^%s"%(p, s))
                need_unscaling = False
            Psi = Phi.hecke(p) * 1/ap
            err = (Psi - Phi).diagonal_valuation(p)
            verbose("error is zero modulo p^%s"%(err))
            Phi = Psi
        return Phi._normalize()

class PSModularSymbolElement_dist(PSModularSymbolElement):


    def reduce_precision(self, M):
        r"""
        Only holds on to `M` moments of each value of self
        """
        return self.__class__(self._map.reduce_precision(M), self.parent(), construct=True)

    def precision_absolute(self):
        r"""
        Returns the number of moments of each value of self
        """
        return min([a.precision_absolute() for a in self._map])

    def specialize(self, new_base_ring=None):
        r"""
        Returns the underlying classical symbol of weight `k` -- i.e.,
        applies the canonical map `D_k --> Sym^k` to all values of
        self.

        EXAMPLES::

            sage: D = Distributions(0, 5, 10);  M = PSModularSymbolSpace(Gamma0(2), D); M
            Space of overconvergent modular symbols for Congruence Subgroup Gamma0(2) with sign 0 and values in Space of 5-adic distributions with k=0 action and precision cap 10
            sage: f = M(1)
            sage: f.specialize()
            Modular symbol with values in Space of 5-adic distributions with k=0 action and precision cap 1
            sage: f.specialize().values()
            [1 + O(5^10), 1 + O(5^10)]
            sage: f.values()
            [1, 1]
            sage: f.specialize().parent()
            Space of overconvergent modular symbols for Congruence Subgroup Gamma0(2) with sign 0 and values in Space of 5-adic distributions with k=0 action and precision cap 1
            sage: f.specialize().parent().coefficient_module()
            Space of 5-adic distributions with k=0 action and precision cap 1
            sage: f.specialize().parent().coefficient_module().is_symk()
            True

            sage: f.specialize(QQ)
            Modular symbol with values in Sym^0 Q^2
            sage: f.specialize(QQ).values()
            [1, 1]
            sage: f.specialize(QQ).parent().coefficient_module()
            Sym^0 Q^2            
        """
        if new_base_ring is None:
            new_base_ring = self.base_ring()
        return self.__class__(self._map.specialize(new_base_ring),
                              self.parent()._specialize_parent_space(new_base_ring), construct=True)

    def _consistency_check(self):
        """
        Check that the map really does satisfy the Manin relations loop (for debugging).
        """
        rels = self.parent()._grab_relations()
        # TODO: no clue how to do this until this object fully works again...
        raise NotImplementedError

