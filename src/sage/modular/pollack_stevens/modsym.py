from sage.structure.element import ModuleElement
from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from manin_map import ManinMap
import operator

from sage.categories.action import Action

M2Z = MatrixSpace_ZZ_2x2()
minusproj = M2Z([1,0,0,-1])

class PSModSymAction(Action):
    def __init__(self, MSspace):
        Action.__init__(self, M2Z, MSspace, False, operator.mul)

    def _call_(self, sym, g):
        return sym.__class__(sym._map * g, sym.parent())

class PSModularSymbolElement(ModuleElement):
    def __init__(self, map_data, parent, construct=False):
        ModuleElement.__init__(self, parent)
        if construct:
            self._map = map_data
        else:
            self._map = ManinMap(parent._coefficients, parent._manin_relations, map_data)

    def _add_(self, right):
        return self.__class__(self._map + right._map, self.parent())

    def _lmul_(self, right):
        return self.__class__(self._map * right, self.parent())

    def _sub_(self, right):
        return self.__class__(self._map - right._map, self.parent())

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
    
## the methods below need the new dict (currently using old code with indices)

    def hecke(self, ell):
        r"""
        Returns self | T_ell by making use of the precomputations in
        self.prep_hecke()

        INPUT:
            - ``ell`` - a prime

        OUTPUT:
        
        self | T_ell 

        EXAMPLES:

        ::

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
        """
            
        if self._data == 0:
            self.compute_full_data_from_gen_data()
            self.normalize_full_data()
        psi = self.zero()   ## psi will denote self | T_ell

        ## v is a long list of lists of lists with the property that the value
        ## of self | T_ell on the m-th generator is given by
        ## sum_j sum_r self(j-th coset rep) | v[m][j][r]
        ## where j runs thru all coset reps and r runs thru all entries
        ## of v[m][j]
        
        v = self.prep_hecke(ell)

        ## This loop computes (self | T_ell)(m-th generator) 

        for m in range(self.parent().ngens()):
            for j in range(self.parent().ncoset_reps()):
                for r in range(len(v[m][j])):
                     psi._data[m] = psi.data(m) + self.full_data(j).act_right(v[m][j][r])

        return psi.normalize()

    def hecke_from_defn(self,ell):
        r"""
        Computes self | T_ell directly from the definition of acting by double
        coset reps.

        That is, it computes sum_a self | [1,a,0,ell] + self | [ell,0,0,1]
        where the last term occurs only if the level is prime to ell.

        INPUT:
             - ``ell`` - prime

        OUTPUT:

        self | T_ell

        EXAMPLES:

        ::

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
        # If needed, precomputes the value of self on all coset reps
        if self.full_data() == 0:
            self.compute_full_data_from_gen_data()
            
        psi = self.zero() + [self.act_right(M2Z([1,a,0,ell])) for a in range(ell)]
        if self.parent().level() % ell != 0:
            psi = psi + self.act_right(M2Z([ell,0,0,1]))
        return psi.normalize()
   

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
        sd = self._map._dict
        return min([val.valuation(p) for  ky,val in sd.iteritems()])

    def change_ring():
        r"""
        """

    def is_Tq_eigensymbol():
        r"""
        """

    def Tq_eigenvalue():
        r"""
        """

    def lift(self, algorithm = ..., eigensymbol = ...):
        r"""
        """
        
class PSModularSymbolElement_symk(PSModularSymbolElement):
    def __init__():

    def p_stabilize(alpha = None, ap = None, ordinary = True):

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
    def __init__():

    def reduce_precision(self, M):
        r"""
        Only holds on to `M` moments of each value of self
        """

    def precision_absolute( ):

    def specialize( ):
