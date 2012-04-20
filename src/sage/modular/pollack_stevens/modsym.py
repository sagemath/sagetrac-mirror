from sage.structure.element import ModuleElement
from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from sage.rings.integer_ring import ZZ
from manin_map import ManinMap
import operator

from sage.categories.action import Action

from fund_domain import M2ZSpace, M2Z
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
            self._map = ManinMap(parent._coefficients, parent._manin_relations, map_data)

    def _repr_(self):
        return "A modular symbol with values in %s"%(self.parent().coefficient_module())

    def dict(self):
        D = {}
        for g in self.parent().source().gens():
            D[g] = self._map[g]
        return D

    def values(self):
        return [self._map[g] for g in self.parent().source().gens()]

    def __cmp__(self, other):
        gens = self.parent().source().gens()
        for g in gens:
            c = cmp(self._map[g], other._map[g])
            if c: return c
        return 0

    def _add_(self, right):
        return self.__class__(self._map + right._map, self.parent(), construct=True)

    def _lmul_(self, right):
        return self.__class__(self._map * right, self.parent(), construct=True)

    def _rmul_(self, right):
        return self.__class__(self._map * right, self.parent(), construct=True)

    def _sub_(self, right):
        return self.__class__(self._map - right._map, self.parent(), construct=True)

    def plus_part(self):
        r"""
        Returns the plus part of self -- i.e. self + self | [1,0,0,-1].

        Note that we haven't divided by 2.  Is this a problem?

        OUTPUT:

        self + self | [1,0,0,-1]

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: (phi.plus_part()+phi.minus_part()) == phi.scale(2)
        True
        """
        return self * minusproj + self

    def minus_part(self):
        r"""
        Returns the minus part of self -- i.e. self - self | [1,0,0,-1]

        Note that we haven't divided by 2.  Is this a problem?

        OUTPUT:

        self - self | [1,0,0,-1]

        EXAMPLES:

        ::
        
        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: (phi.plus_part()+phi.minus_part()) == phi.scale(2)
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
            sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
            sage: phi = form_modsym_from_elliptic_curve(E); phi
            [-1/5, 3/2, -1/2]
            sage: ell=2
            sage: phi.hecke(ell) == phi.scale(E.ap(ell))
            True
            sage: ell=3
            sage: phi.hecke(ell) == phi.scale(E.ap(ell))
            True
            sage: ell=5
            sage: phi.hecke(ell) == phi.scale(E.ap(ell))
            True
            sage: ell=101
            sage: phi.hecke(ell) == phi.scale(E.ap(ell))
            True

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
            sage: phi = form_modsym_from_elliptic_curve(E); phi
            [-1/5, 3/2, -1/2]
            sage: ell = 2
            sage: phi.hecke_from_defn(ell) == phi.scale(E.ap(ell))
            True
            sage: ell = 3
            sage: phi.hecke_from_defn(ell) == phi.scale(E.ap(ell))
            True
            sage: ell=5
            sage: phi.hecke_from_defn(ell) == phi.scale(E.ap(ell))
            True
            sage: ell=101
            sage: phi.hecke_from_defn(ell) == phi.scale(E.ap(ell))
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

        The valuation of self at `p`

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
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
        return min([val for val in self._map])

    def change_ring(self,R):
        r"""
        Changes the base ring of self to `R`
        """
        pass

    def is_Tq_eigensymbol(self,q,p,M):
        r"""
        Determines if self is an eigenvector for `T_q` modulo `p^M` 

        INPUT:
            - ``q`` -- prime of the Hecke operator
            - ``p`` -- prime we are working modulo
            - ``M`` -- degree of accuracy of approximation

        OUTPUT:
        True/False

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi_ord = phi.p_stabilize_ordinary(3,E.ap(3),10)
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
        pass
    
    def Tq_eigenvalue(self,q,p,M):
        r"""
        Eigenvalue of `T_q` modulo `p^M`

        INPUT:
        - ``q`` -- prime of the Hecke operator
        - ``p`` -- prime we are working modulo
        - ``M`` -- degree of accuracy of approximation     
        
        OUTPUT:

        Constant `c` such that `self|T_q - c * self` has valuation greater than
        or equal to `M` (if it exists), otherwise raises ValueError

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi_ord = phi.p_stabilize_ordinary(3,E.ap(3),10)
        sage: phi_ord.Tq_eigenvalue(2,3,10)
        -2
        sage: phi_ord.Tq_eigenvalue(2,3,100)
        -2
        sage: phi_ord.Tq_eigenvalue(2,3,1000)
        -2
        sage: phi_ord.Tq_eigenvalue(3,3,10)
        ...
        ValueError: No eigenvalue exists modulo 3^10
        """
        pass
    
    def lift(self, algorithm = None, eigensymbol = None):
        r"""
        """
        pass
    
class PSModularSymbolElement_symk(PSModularSymbolElement):
    
    def p_stabilize(alpha = None, ap = None, ordinary = True):
        pass
    
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

        EXAMPLES:
        """
        K = self.base_ring()
        f = K.defining_polynomial()
        R = pAdicField(p,M+10)['x']
        x = R.gen()
        v = R(f).roots()
        if len(v) == 0:
            raise ValueError, "No coercion possible -- no prime over p has degree 1"
        else:
            roots = [r[0] for r in v]
            ans = []
            for r in root:
                psi = K.hom([root],pAdicField(p,M))
                ans.append([self.map(psi),psi])
            return ans
        
        
class PSModularSymbolElement_dist(PSModularSymbolElement):

    
    def reduce_precision(self, M):
        r"""
        Only holds on to `M` moments of each value of self
        """
        sd = self._dict
        for val in sd.itervalues():
            val.reduce_precision(M)
        return self
    
    def precision_absolute(self):
        r"""
        Returns the number of moments of each value of self
        """
        return self.precision_cap()
    
    def specialize(self):
        r"""
        Returns the underlying classical symbol of weight `k` -- i.e.,
        applies the canonical map `D_k --> Sym^k` to all values of
        self
        """
        sd = self._dict
        for val in sd.itervalues():
            val.specialize()
        return self
