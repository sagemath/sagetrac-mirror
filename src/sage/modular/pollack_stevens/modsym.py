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
    def __init__(self, map_data, parent):
        ModuleElement.__init__(self, parent)
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


    def act_right(self, gamma):
        r"""
        Returns self | gamma

        This action is defined by (self | gamma)(D) = self(gamma D)|gamma

        For this to work gamma must normalize Gamma_0(N) and be able to act on
        the values of self.  However, it can also be used to define Hecke
        operators.  Even if each individual self | gamma is not really defined
        on Gamma_0(N), the sum over acting by the appropriate double coset reps
        will be defined over Gamma_0(N).

        INPUT:
            - ``gamma`` - 2 x 2 matrix which acts on the values of self

        OUTPUT:

        self | gamma

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.act_right(Matrix(ZZ,2,2,[1,0,0,-1]))
        [-1/5, -1/2, 3/2]
        """
        v = []  ## This will be the data defining self | gamma

        ##this loop runs over all generators

        #for j in range(0,self.ngens()):
        for j in range(self.parent().ngens()):
            rj = self.manin().generator_indices(j)
            v = v + [self.eval(gamma * self.manin().coset_reps(rj)).act_right(gamma)]

        C = type(self)
        ans = C(v, self.manin())
        ans.normalize()

        return ans

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
        if self.full_data() == 0:
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

        for m in range(self.ngens()):
            for j in range(self.ncoset_reps()):
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
        ## If needed, precomputes the value of self on all coset reps
        if self.full_data() == 0:
            self.compute_full_data_from_gen_data()

        psi = self.zero()
        for a in range(ell):
            psi = psi + self.act_right(M2Z([1,a,0,ell]))
        if self.level()%ell != 0:
            psi = psi + self.act_right(M2Z([ell,0,0,1]))

        return psi.normalize()

    def normalize(self):
        r"""
        Normalizes all of the values of the symbol self.
        """
        #for j in range(self.ngens()):
        for j in range(self.parent().ngens()):
            self._data[j].normalize()

        if self.full_data() != 0:
            for j in range(self.parent().ncoset_reps()):
            #for j in range(self.ncoset_reps()):
                self._full_data[j].normalize()

        return self
