r"""
    Represents maps from a set of right coset representatives to a coefficient module.
    
    This is a class that represents maps from a set of right coset
    representatives to a coefficient module.  This is a basic building
    block for implementing modular symbols, and provides basic arithmetic
    and right action of matrices.    
    
    Lots and lots of examples.
    """

#*****************************************************************************
#       Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.arith import convergents
from sage.misc.misc import verbose
from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2, Matrix_integer_2x2
from fund_domain import M2Z, t00, t10, t01, t11, Id, basic_hecke_matrix
from distributions import Distributions

def unimod_matrices_to_infty(r, s):
    r"""
    Return a list of matrices whose associated unimodular paths connect `0` to ``r/s``.

    INPUT:

    - ``r``, ``s`` -- rational numbers

    OUTPUT:

    - a list of matrices in `SL_2(\ZZ)`

    EXAMPLES::

        sage: v = sage.modular.pollack_stevens.manin_map.unimod_matrices_to_infty(19,23); v
        [
        [1 0]  [ 0  1]  [1 4]  [-4  5]  [ 5 19]
        [0 1], [-1  1], [1 5], [-5  6], [ 6 23]
        ]
        sage: [a.det() for a in v]
        [1, 1, 1, 1, 1]

        sage: sage.modular.pollack_stevens.manin_map.unimod_matrices_to_infty(11,25)
        [
        [1 0]  [ 0  1]  [1 3]  [-3  4]  [ 4 11]
        [0 1], [-1  2], [2 7], [-7  9], [ 9 25]
        ]

    ALGORITHM:

    This is Manin's continued fraction trick, which gives an expression
    `{0,r/s} = {0,\infty} + ... + {a,b} + ... + {*,r/s}`, where each `{a,b}` is
    the image of `{0,\infty}` under a matrix in `SL_2(\ZZ)`.

    """
    if s == 0:
        return []
    # the function contfrac_q in
    # https://github.com/williamstein/psage/blob/master/psage/modform/rational/modular_symbol_map.pyx
    # is very, very relevant to massively optimizing this.
    L = convergents(r / s)
    # Computes the continued fraction convergents of r/s
    v = [M2Z([1, L[0].numerator(), 0, L[0].denominator()])]
    # Initializes the list of matrices
    for j in range(0, len(L)-1):
        a = L[j].numerator()
        c = L[j].denominator()
        b = L[j + 1].numerator()
        d = L[j + 1].denominator()
        v.append(M2Z([(-1)**(j + 1) * a, b, (-1)**(j + 1) * c, d]))
        # The matrix connecting two consecutive convergents is added on
    return v


def unimod_matrices_from_infty(r, s):
    r"""
    Return a list of matrices whose associated unimodular paths connect `\infty` to ``r/s``.

    INPUT:

    - ``r``, ``s`` -- rational numbers

    OUTPUT:

    - a list of `SL_2(\ZZ)` matrices

    EXAMPLES::

        sage: v = sage.modular.pollack_stevens.manin_map.unimod_matrices_from_infty(19,23); v
        [
        [ 0  1]  [-1  0]  [-4  1]  [-5 -4]  [-19   5]
        [-1  0], [-1 -1], [-5  1], [-6 -5], [-23   6]
        ]
        sage: [a.det() for a in v]
        [1, 1, 1, 1, 1]
        
        sage: sage.modular.pollack_stevens.manin_map.unimod_matrices_from_infty(11,25)
        [
        [ 0  1]  [-1  0]  [-3  1]  [-4 -3]  [-11   4]
        [-1  0], [-2 -1], [-7  2], [-9 -7], [-25   9]
        ]
        
    ALGORITHM:
        
    This is Manin's continued fraction trick, which gives an expression
    `{\infty,r/s} = {\infty,0} + ... + {a,b} + ... + {*,r/s}`, where each
    `{a,b}` is the image of `{0,\infty}` under a matrix in `SL_2(\ZZ)`.
        
    """
    if s != 0:
        L = convergents(r / s)
        # Computes the continued fraction convergents of r/s
        v = [M2Z([-L[0].numerator(), 1, -L[0].denominator(), 0])]
        # Initializes the list of matrices
        # the function contfrac_q in https://github.com/williamstein/psage/blob/master/psage/modform/rational/modular_symbol_map.pyx
        # is very, very relevant to massively optimizing this.
        for j in range(0, len(L) - 1):
            a = L[j].numerator()
            c = L[j].denominator()
            b = L[j + 1].numerator()
            d = L[j + 1].denominator()
            v.append(M2Z([-b, (-1)**(j + 1) * a, -d, (-1)**(j + 1) * c]))
            # The matrix connecting two consecutive convergents is added on
        return v
    else:
        return []

class ManinMap(object):
    r"""
    Map from a set of right coset representatives of `\Gamma_0(N)` in
    `SL_2(\ZZ)` to a coefficient module that satisfies the Manin
    relations.
    """
    def __init__(self, codomain, manin_relations, defining_data, check=True):
        """
        INPUT:

        - ``codomain`` -- coefficient module
        - ``manin_relations`` -- a ManinRelations object
        - ``defining_data`` -- a dictionary whose keys are a superset of
          manin_relations.gens() and a subset of manin_relations.reps(),
          and whose values are in the codomain.
        - ``check`` -- do numerous (slow) checks and transformations to
          ensure that the input data is perfect.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap, Distributions
            sage: D = Distributions(0, 5, 10); D
            Space of 5-adic distributions with k=0 action and precision cap 10
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data); f
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 5-adic distributions with k=0 action and precision cap 10
            sage: f(M2Z([1,0,0,1]))
            (1, 2)
            
        """
        self._codomain = codomain
        self._manin = manin_relations
        if check:
            self._dict = {}
            if isinstance(defining_data, (list, tuple)):
                if len(defining_data) != manin_relations.ngens():
                    raise ValueError("length of defining data must be the same as number of manin generators")
                for i in range(len(defining_data)):
                    self._dict[manin_relations.gen(i)] = defining_data[i]
            elif isinstance(defining_data, dict):
                for ky, val in defining_data.iteritems():
                    if not isinstance(ky, Matrix_integer_2x2):
                        # should eventually check that ky is actually a coset rep,
                        # handle elements of P^1, make sure that we cover all cosets....
                        ky = M2Z(ky)
                    self._dict[ky] = val
            else:
                # constant function
                try:
                    c = codomain(defining_data)
                except TypeError:
                    raise TypeError("unrecognized type for defining_data")
                g = manin_relations.gens()
                self._dict = dict(zip(g, [c]*len(g)))
        else:
            self._dict = defining_data
            
    def _compute_image_from_gens(self, B):
        r"""
        Compute image of ``B`` under ``self``.

        INPUT:
            
        - ``B`` --  generator of Manin relations.
            
        OUTPUT:
            
        - an element in the codomain of self (e.g. a distribution), the image of ``B`` under ``self``.
        
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap, Distributions
            sage: D = Distributions(0, 5, 10)
            sage: MR = ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, MR, data)
            sage: f._compute_image_from_gens(MR.reps()[1])
            (24, 0)
        
        """
        L = self._manin.relations(B)
        # could raise KeyError if B is not a generator
        if len(L) == 0:
            t = self._codomain.zero_element()
        else:
            c, A, g = L[0]
            A=M2Z(A)
            A.set_immutable()
            g1 = (self._dict[self._manin.reps(g)] * A)
            t = g1 * c
            for c, A, g in L[1:]:
                A=M2Z(A)
                A.set_immutable()
                g1 = (self._dict[self._manin.reps(g)] * A)
                t += g1 * c
        return t

    def __getitem__(self, B):
        r"""
        
        Compute image of ``B`` under ``self``.
        
        INPUT:
            
        - ``B`` -- coset representative of Manin relations.
            
        OUTPUT:
            
        - an element in the codomain of self (e.g. a distribution), the image of ``B`` under ``self``.
            
        EXAMPLES::
            
            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: S = Symk(0,QQ)
            sage: MR = ManinRelations(37); MR.gens()
            [
            [1 0]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -3]  [-3 -1]  [-1 -4]  [-4 -3]
            [0 1], [ 1  4], [ 4  3], [ 3  5], [ 5  7], [ 7  2], [ 2  7], [ 7  5],
            
            [-2 -3]
            [ 3  4]
            ]
            sage: data  = {M2Z([-2,-3,5,7]): 0, M2Z([1,0,0,1]): 0, M2Z([-1,-2,3,5]): 0, M2Z([-1,-4,2,7]): 1, M2Z([0,-1,1,4]): 1, M2Z([-3,-1,7,2]): -1, M2Z([-2,-3,3,4]): 0, M2Z([-4,-3,7,5]): 0, M2Z([-1,-1,4,3]): 0}
            sage: f = ManinMap(D,MR,data)
            sage: f.__getitem__(MR.gens()[1])
            1
            sage: f.__getitem__(MR.gens()[3])
            0
            sage: f.__getitem__(MR.gens()[5])
            -1
            
        """
        try:
            return self._dict[B]
        except KeyError:
            # To prevent memory overflow
            return self._compute_image_from_gens(B)
            # self._dict[B] = self._compute_image_from_gens(B)
            # return self._dict[B]

    def clear_cache(self):
        r"""
        Clear the cache of ``self``.
            
        EXAMPLES::
            
            sage: D = Distributions(2,23,40)
            sage: MR = ManinRelations(13)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D,MR,data)
            sage: f._dict == {}
            False
            sage: f.clear_cache
            <bound method ManinMap.clear_cache of Map from the set of right cosets of Gamma0(13) in SL_2(Z) to Space of 23-adic distributions with k=2 action and precision cap 40>
            sage: f._dict == {}
            True
            
        """
        self._dict = {}
        self.compute_full_data()

    def compute_full_data(self):
        r"""
        Computes the values of self on all coset reps from its values on our generating set.
        """
        for B in self._manin.reps():
            if not self._dict.has_key(B):
                self._dict[B] = self._compute_image_from_gens(B)

    def __add__(self, right):
        r"""
        Return sum self + right, where self and right are
        assumed to have identical codomains and Manin relations.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap, Distributions
            sage: D = Distributions(0, 5, 10); D
            Space of 5-adic distributions with k=0 action and precision cap 10
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data); f
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 5-adic distributions with k=0 action and precision cap 10
            sage: f(M2Z([1,0,0,1]))
            (1, 2)
            sage: f+f
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 5-adic distributions with k=0 action and precision cap 10
            sage: (f+f)(M2Z([1,0,0,1]))
            (2, 4)
        """
        D = {}
        sd = self._dict
        rd = right._dict
        for ky, val in sd.iteritems():
            if ky in rd:
                D[ky] = val + rd[ky]
        return self.__class__(self._codomain, self._manin, D, check=False)

    def __sub__(self, right):
        """
        Return difference self - right, where self and right are
        assumed to have identical codomains and Manin relations.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap
            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(0, 5, 10); D
            Space of 5-adic distributions with k=0 action and precision cap 10
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11) 
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data); f
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 5-adic distributions with k=0 action and precision cap 10
            sage: f(M2Z([1,0,0,1]))
            (1, 2)
            sage: f-f
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 5-adic distributions with k=0 action and precision cap 10
            sage: (f-f)(M2Z([1,0,0,1]))
            (0, 0)
        
        """
        D = {}
        sd = self._dict
        rd = right._dict
        for ky, val in sd.iteritems():
            if ky in rd:
                D[ky] = val - rd[ky]
        return self.__class__(self._codomain, self._manin, D, check=False)

    def __mul__(self, right):
        """
        Return scalar multiplication self * right, where right is in the
        base ring of the codomain.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap, Distributions
            sage: D = Distributions(0, 5, 10)
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data)
            sage: f(M2Z([1,0,0,1]))
            (1, 2)
            sage: f*2
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 5-adic distributions with k=0 action and precision cap 10
            sage: (f*2)(M2Z([1,0,0,1]))
            (2, 4)
        
        """
        if isinstance(right, Matrix_integer_2x2):
            return self._right_action(right)
        D = {}
        sd = self._dict
        for ky, val in sd.iteritems():
            D[ky] = val * right
        return self.__class__(self._codomain, self._manin, D, check=False)

    def __repr__(self):
        """
        Return print representation of self.

        EXAMPLES::
 
            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap, Distributions
            sage: D = Distributions(0, 5, 10); D
            Space of 5-adic distributions with k=0 action and precision cap 10
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data)
            sage: f.__repr__()
            'Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 5-adic distributions with k=0 action and precision cap 10'
            
        """
        return "Map from the set of right cosets of Gamma0(%s) in SL_2(Z) to %s"%(
            self._manin.level(), self._codomain)
    
    def _eval_sl2(self, A):
        r"""
        Return the value of self on the unimodular divisor corresponding to `A`.

        Note that `A` must be in `SL_2(Z)` for this to work.
        
        INPUT:
            
        - ``A`` - an element of `SL_2(Z)`

        OUTPUT:

        - The value of self on the divisor corresponding to `A` -- i.e. on the divisor `{A(0)} - {A(\infty)}`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap, Distributions
            sage: D = Distributions(0, 5, 10)
            sage: MR = ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, MR, data)   
            sage: A = MR.reps()[1]
            sage: f._eval_sl2(A)
            (15, 0)
            
        """

        B = self._manin.equivalent_rep(A)
        gaminv = M2Z(B * A._invert_unit())
        gaminv.set_immutable()
        return self[B] * gaminv

    def __call__(self, A):
        """
        Evaluate self at A.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.manin_map import M2Z, ManinMap, Distributions
            sage: D = Distributions(0, 5, 10); D
            Space of 5-adic distributions with k=0 action and precision cap 10
            sage: manin = sage.modular.pollack_stevens.fund_domain.ManinRelations(11)
            sage: data  = {M2Z([1,0,0,1]):D([1,2]), M2Z([0,-1,1,3]):D([3,5]), M2Z([-1,-1,3,2]):D([1,1])}
            sage: f = ManinMap(D, manin, data); f
            Map from the set of right cosets of Gamma0(11) in SL_2(Z) to Space of 5-adic distributions with k=0 action and precision cap 10
            sage: f(M2Z([1,0,0,1]))
            (1, 2)
        
        """
        a = A[t00]
        b = A[t01]
        c = A[t10]
        d = A[t11]
        # v1: a list of unimodular matrices whose divisors add up to {b/d} - {infty}
        v1 = unimod_matrices_to_infty(b,d)
        # v2: a list of unimodular matrices whose divisors add up to {a/c} - {infty}
        v2 = unimod_matrices_to_infty(a,c)
        # ans: the value of self on A
        ans = self._codomain.zero_element()
        # This loop computes self({b/d}-{infty}) by adding up the values of self on elements of v1
        for B in v1:
            ans = ans + self._eval_sl2(B)

        # This loops subtracts away the value self({a/c}-{infty}) from ans by subtracting away the values of self on elements of v2
        # and so in the end ans becomes self({b/d}-{a/c}) = self({A(0)} - {A(infty)}
        for B in v2:
            ans = ans - self._eval_sl2(B)
        return ans

    def apply(self, f, codomain=None, to_moments=False):
        r"""
        Return Manin map given by `x \mapsto f(self(x))`, where `f` is
        anything that can be called with elements of the coefficient
        module.

        This might be used to normalize, reduce modulo a prime, change
        base ring, etc.
        """
        D = {}
        sd = self._dict
        if codomain is None:
            codomain = self._codomain
        for ky, val in sd.iteritems():
            if to_moments:
                D[ky] = codomain([f(val.moment(a)) for a in range(val.precision_absolute())])
            else:
                D[ky] = f(val)
        return self.__class__(codomain, self._manin, D, check=False)

    def __iter__(self):
        r"""
        Return iterator over the values of this map on the reduced
        representatives.

        This might be used to compute the valuation.
        """
        for A in self._manin.gens():
            yield self._dict[A]

    def _right_action(self, gamma):
        r"""
        Return self | gamma, where gamma is a 2x2 integer matrix.

        The action is defined by `(self | gamma)(D) = self(gamma D)|gamma`

        For the action by a single element gamma to be well defined,
        gamma must normalize `\Gamma_0(N)`.  However, this right action
        can also be used to define Hecke operators, in which case each
        individual self | gamma is not a modular symbol on `\Gamma_0(N)`,
        but the sum over acting by the appropriate double coset
        representatives is.

        INPUT:

        - ``gamma`` - 2 x 2 matrix which acts on the values of self

        OUTPUT:

        - ManinMap
            
        EXAMPLES
        """
        D = {}
        sd = self._dict
        # we should eventually replace the for loop with a call to apply_many
        keys = [ky for ky in sd.iterkeys()]
        for ky in keys:
            D[ky] = self(gamma*ky) * gamma
        return self.__class__(self._codomain, self._manin, D, check=False)

    def normalize(self):
        r"""
        Normalize every value of self -- e.g., reduces each value's
        `j`-th moment modulo `p^(N-j)`
        """
        sd = self._dict
        for val in sd.itervalues():
            val.normalize()
        return self

    def reduce_precision(self, M):
        r"""
        
        """
        D = {}
        sd = self._dict
        for ky, val in sd.iteritems():
            D[ky] = val.reduce_precision(M)
        return self.__class__(self._codomain, self._manin, D, check=False)

    def specialize(self, new_base_ring):
        D = {}
        sd = self._dict
        for ky, val in sd.iteritems():
            D[ky] = val.specialize(new_base_ring)
        return self.__class__(self._codomain.specialize(new_base_ring), self._manin, D, check=False)

    def hecke(self, ell, algorithm = 'prep'):
        r"""
        Return the image of this Manin map under the Hecke operator `T_{\ell}`.

        INPUT:

        - ``ell`` -- a prime

        - ``algorithm`` -- a string, either 'prep' (default) or
          'naive'

        OUTPUT:

        - The image of this ManinMap under the Hecke operator
          `T_{\ell}`

        EXAMPLES:

        ::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi.values()
            [-1/5, 3/2, -1/2]
            sage: phi.is_Tq_eigensymbol(7,7,10)
            True
            sage: phi.hecke(7).values()
            [2/5, -3, 1]
            sage: phi.Tq_eigenvalue(7,7,10)
            -2
        """
        self.compute_full_data()
        self.normalize()
        M = self._manin
        if algorithm == 'prep':
            ## psi will denote self | T_ell
            psi = {}
            for g in M.gens():
                ## v is a dictionary so that the value of self | T_ell
                ## on g is given by
                ## sum_h sum_A self(h) * A
                ## where h runs over all coset reps and A runs over
                ## the entries of v[h] (a list)
                # verbose("prepping for T_%s: %s"%(ell, g), level = 2)
                v = M.prep_hecke_on_gen(ell, g)
                psi[g] = self._codomain.zero_element()
                for h in M:
                    for A in v[h]:
                        psi[g] += self[h] * A
                psi[g].normalize()
            return self.__class__(self._codomain, self._manin, psi, check=False)
        elif algorithm == 'naive':
            psi = self._right_action(M2Z([1,0,0,ell]))
            for a in range(1, ell):
                psi += self._right_action(M2Z([1,a,0,ell]))
            if self._manin.level() % ell != 0:
                psi += self._right_action(M2Z([ell,0,0,1]))
            return psi.normalize()

    def p_stabilize(self, p, alpha, V):
        manin = V.source()
        pmat = M2Z([p,0,0,1])
        D = {}
        scalar = 1/alpha
        one = scalar.parent()(1)
        for g in manin.gens():
            # we use scale here so that we don't need to define a
            # construction functor in order to scale by something
            # outside the base ring.
            D[g] = self._eval_sl2(g).scale(one) - (self(pmat * g) * pmat).scale(1/alpha)
        return self.__class__(self._codomain.change_ring(scalar.parent()), manin, D, check=False)
