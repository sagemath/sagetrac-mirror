"""
Manin Maps

This is a class that represents maps from a set of right coset
representatives to a coefficient module.  This is a basic building
block for implementing modular symbols, and provides basic arithmetic
and right action of matrices.


"""
from sage.rings.arith import convergents
from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2, Matrix_integer_2x2

from distributions import Distributions
from fund_domain import M2Z, t00, t10, t01, t11, Id, unimod_matrices_to_infty

def unimod_matrices_to_infty(r, s):
    r"""
    Returns a list of matrices whose associated unimodular paths
    connect `0` to `r/s`.  This is Manin's continued fraction trick, which
    gives an expression `{0,r/s} = {0,oo} + ... + {a,b} + ... + {*,r/s}`,
    where each `{a,b}` is the image of `{0,oo}` under a matrix in `SL_2(ZZ)`.

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
    Returns a list of matrices whose associated unimodular paths
    connect `0` to `r/s`.  This is Manin's continued fraction trick, which
    gives an expression `{oo,r/s} = {oo,0} + ... + {a,b} + ... + {*,r/s}`,
    where each `{a,b}` is the image of `{0,oo}` under a matrix in `SL_2(ZZ)`.

    INPUT:

    - ``r``, ``s`` -- rational numbers

    OUTPUT:

    - a list of `SL_2(Z)` matrices

    EXAMPLES:

        sage: v = sage.modular.pollack_stevens.manin_map.unimod_matrices_from_infty(19,23); v
        [
        [ 0  1]  [-1  0]  [-4  1]  [-5 -4]  [-19   5]
        [-1  0], [-1 -1], [-5  1], [-6 -5], [-23   6]
        ]
        sage: [a.det() for a in v]
        [1, 1, 1, 1, 1]    
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

def basic_hecke_matrix(a, ell):
    r"""
    Returns the matrix [1, a, 0, ell] (if `a < ell`) and [ell, 0, 0, 1] if `a \geq ell`

    INPUT:

    - ``a`` -- an integer or Infinity
    - ``ell`` -- a prime

    OUTPUT:

    - a 2 x 2 matrix of determinant `ell`

    EXAMPLES:

        sage: sage.modular.pollack_stevens.manin_map.basic_hecke_matrix(0, 7)
        [1 0]
        [0 7]
        sage: sage.modular.pollack_stevens.manin_map.basic_hecke_matrix(5, 7)
        [1 5]
        [0 7]
        sage: sage.modular.pollack_stevens.manin_map.basic_hecke_matrix(7, 7)
        [7 0]
        [0 1]
        sage: sage.modular.pollack_stevens.manin_map.basic_hecke_matrix(19, 7)
        [7 0]
        [0 1]    
    """
    # TODO: probably a bottleneck.
    if a < ell:
        return M2Z([1, a, 0, ell])
    else:
        return M2Z([ell, 0, 0, 1])

class ManinMap(object):
    """
    Map from a set of right coset representatives of `\Gamma_0(N)` in
    `SL_2(Z)` to a coefficient module that satisfies the Manin
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
            sage: f(M2Z([0,-1,1,3]))
            (3, 5)
            sage: f(M2Z([-1,-1,3,2]))
            (1, 1)
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
                raise TypeError("unrecognized type for defining_data")
        else:
            self._dict = defining_data

    def _compute_image_from_gens(self, B):
        L = self._manin.relations(B)
        # could raise KeyError if B is not a coset rep
        if len(L) == 0:
            t = self._codomain.zero_element()
        else:
            c, A, g = L[0]
            t = (self._dict[g] * A) * c
            for c, A, g in L[1:]:
                t += (self._dict[g] * A) * c
        return t

    def __getitem__(self, B):
        try:
            return self._dict[B]
        except KeyError:
            self._dict[B] = self._compute_image_from_gens(B)
            return self._dict[B]

    def compute_full_data(self):
        r"""
        Computes the values of self on all coset reps from its values on our generating set.
        """
        for B in self._manin.reps():
            if not self._dict.has_key(B):
                self._dict[B] = self._compute_image_from_gens(B)

    def __add__(self, right):
        """
        Return sum self + right, where self and right are
        assumed to have identical codomains and Manin relations.
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
        Return scalar multiplication self*right, where right is in the
        base ring of the codomain.
        """
        if isinstance(right, Matrix_integer_2x2):
            return self._right_action(right)
        D = {}
        sd = self._dict
        for ky, val in sd.iteritems():
            D[ky] = val * right
        return self.__class__(self._codomain, self._manin, D, check=False)

    def __repr__(self):
        return "Map from the set of right cosets of Gamma0(%s) in SL_2(Z) to %s"%(
            self._manin.level(), self._codomain)

    def _eval_sl2(self, A):
        r"""
        Returns the value of self on the unimodular divisor corresponding to `A`.

        Note that `A` must be in `SL_2(Z)` for this to work.

        INPUT:
            - ``A`` - an element of `SL_2(Z)`

        OUTPUT:

        The value of self on the divisor corresponding to `A` -- i.e. on the divisor `{A(0)} - {A(\infty)}`.

        EXAMPLES:

        ::

        sage: E = EllipticCurve('11a')
        sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
        sage: phi = form_modsym_from_elliptic_curve(E); phi
        [-1/5, 3/2, -1/2]
        sage: phi.manin().generator_indices()
        [0, 2, 3]
        sage: A = phi.manin().coset_reps(3); A
        [-1 -1]
        [ 3  2]
        sage: phi.eval_sl2(A) == phi.data(2)
        True
        sage: sig = Matrix(2,2,[0,1,-1,0])
        sage: A = Matrix(ZZ,2,2,[8,29,3,11]); A
        [ 8 29]
        [ 3 11]
        sage: A.determinant()
        1
        sage: (phi.eval_sl2(A)+phi.eval_sl2(A*sig)) == phi.zero_elt()
        True
        sage: tau = Matrix(2,2,[0,-1,1,-1])
        sage: (phi.eval_sl2(A)+phi.eval_sl2(A*tau)+phi.eval_sl2(A*tau*tau)) == phi.zero_elt()
        True
        """
        
        B = self._manin.equivalent_rep(A)
        gaminv = B * A._invert_unit()
        return self[B] * gaminv
    
    def __call__(self, A):
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

    def apply(self, f):
        """
        Returns Manin map given by `x |--> f(self(x))`, where `f` is
        anything that can be called with elements of the coefficient
        module.

        This might be used to normalize, reduce modulo a prime, change
        base ring, etc.
        """
        D = {}
        sd = self._dict
        for ky, val in sd.iteritems():
            D[ky] = f(val)
        return self.__class__(self._codomain, self._manin, D, check=False)

    def __iter__(self):
        """
        Returns iterator over the values of this map on the reduced
        representatives.

        This might be used to compute the valuation.
        """
        for A in self._manin.gens():
            yield self._dict[A]

    def _right_action(self, gamma):
        """
        Returns self | gamma, where gamma is a 2x2 integer matrix.

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
        """
        D = {}
        sd = self._dict
        # we should eventually replace the for loop with a call to apply_many
        for ky in sd:
            D[ky] = self(gamma*ky) * gamma
        return self.__class__(self._codomain, self._manin, D, check=False)

    def normalize(self):
        r"""
        Normalizes every value of self -- e.g., reduces each value's
        `j`-th moment modulo `p^(N-j)`
        """
        sd = self._dict
        for val in sd.itervalues():
            val.normalize()
        return self

    def hecke(self, ell, algorithm = 'prep'):
        """
        Returns the image of this Manin map under the Hecke operator `T_{\ell}`.

        INPUT:

        - ``ell`` -- a prime

        - ``algorithm`` -- a string, either 'prep' (default) or
          'naive'

        OUTPUT:

        - The image of this ManinMap under the Hecke operator
          `T_{\ell}`

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: from sage.modular.overconvergent.pollack.modsym_symk import form_modsym_from_elliptic_curve
            sage: phi = form_modsym_from_elliptic_curve(E); phi
            [-1/5, 3/2, -1/2]
            sage: phi.prep_hecke(2)
            [[[[1 0]
            [0 2], [1 1]
            [0 2], [2 0]
            [0 1]], [], [], [], [], [], [], [], [], [], [], [[ 1 -1]
            [ 0  2]]], [[[1 2]
            [0 2], [1 1]
            [0 2], [2 1]
            [0 1]], [[ 1 -2]
            [ 0  2], [ 1 -1]
            [ 0  2], [ 2 -1]
            [ 0  1]], [], [[-4 -2]
            [11  5], [-8 -3]
            [22  8]], [], [], [], [], [], [[1 0]
            [0 2]], [[-5 -2]
            [11  4], [-1  1]
            [ 0 -2]], []], [[[1 2]
            [0 2], [1 1]
            [0 2], [2 1]
            [0 1]], [[1 0]
            [0 2], [ 1 -1]
            [ 0  2], [2 0]
            [0 1]], [], [], [[-6 -4]
            [11  7]], [[-7 -4]
            [11  6], [-1  1]
            [ 0 -2], [-2  0]
            [ 0 -1]], [], [], [[-1  0]
            [ 0 -2]], [[1 0]
            [0 2]], [], [[-5 -2]
            [11  4]]]]

        WARNING: changed from this (which disagreed with .sage file)::

            sage: phi.prep_hecke(2)
            [[[[1 0]
            [0 2], [1 1]
            [0 2], [2 0]
            [0 1]], [], [], [], [], [], [[ 1 -1]
            [ 0  2]], [], [], [], [], []], [[[1 2]
            [0 2], [1 1]
            [0 2], [2 1]
            [0 1]], [[ 1 -2]
            [ 0  2], [ 1 -1]
            [ 0  2], [ 2 -1]
            [ 0  1]], [], [[-4 -2]
            [11  5], [-8 -3]
            [22  8]], [], [], [], [[-5 -2]
            [11  4], [-1  1]
            [ 0 -2]], [[1 0]
            [0 2]], [], [], []], [[[1 2]
            [0 2], [1 1]
            [0 2], [2 1]
            [0 1]], [[1 0]
            [0 2], [ 1 -1]
            [ 0  2], [2 0]
            [0 1]], [], [], [[-6 -4]
            [11  7]], [[-7 -4]
            [11  6], [-1  1]
            [ 0 -2], [-2  0]
            [ 0 -1]], [[-5 -2]
            [11  4]], [], [[1 0]
            [0 2]], [[-1  0]
            [ 0 -2]], [], []]]
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


