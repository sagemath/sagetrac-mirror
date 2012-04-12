from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.arith import binomial, bernoulli
from sage.modules.free_module_element import vector, zero_vector
from sage.matrix.matrix cimport Matrix
from sage.matrix.all import matrix, MatrixSpace
from sage.misc.prandom import random
from sage.functions.other import floor
from sage.structure.element cimport RingElement, Element
import operator
#from sage.modular.overconvergent.pollack.S0p import S0
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.padics.padic_capped_absolute_element cimport pAdicCappedAbsoluteElement
from sage.rings.integer cimport Integer

M2Z = MatrixSpace(ZZ,2,2)

include "../../../ext/stdsage.pxi"

def get_dist_classes(p, prec_cap, base):
    if isinstance(base, pAdicGen
    if 7*p**(prec_cap) < ZZ(2)**(4*sizeof(long)-1):
        return Dist_long, WeightKAction_long
    else:
        return Dist_vector, WeightKAction_vector
    

cdef class Dist(ModuleElement):
    cpdef normalize(self):
        raise NotImplementedError

    def scale(self,left):
        return self * left

    def valuation(self):
        """
	returns the highest power of p which divides all moments of the distribution
	"""
        p = self.parent()._p
        return min([self.moment(a).valuation(p) for a in range(self.num_moments())])

    def specialize(self):
        """
	Specializes to weight k -- i.e. projects to Sym^k
	"""
        from sage.modular.overconvergent.pollack.symk import symk
        k=self.parent()._k
        assert k >= 0,"negative weight"
        # R.<X,Y>=PolynomialRing(QQ,2)
        R = PolynomialRing(QQ,('X','Y'))
        X,Y = R.gens()
        P=0
        for j in range(0,k+1):
            P=P+binomial(k,j)*((-1)**j)*self.moment(j)*(X**j)*(Y**(k-j))
        return symk(k,P)

    def act_right(self,gam):
        """
	Return self|gam
	"""
        return self.parent()._act.act(self, gam)

cdef class Dist_vector(Dist):
    def __init__(self,space,moments,check=True):
        """
	A distribution is stored as a vector whose j-th entry is the j-th moment of the distribution.
        The j-th entry is stored modulo p^(N-j) where N is the total number of moments.
        (This is the accuracy that is maintained after acting by Gamma_0(p).)

        INPUTS:

        - ``space`` -- a :class:`distributions.Distributions_Zp` instance
	- ``moments`` -- the list of moments given as a vector
        - ``check`` -- boolean, whether to coerce the vector into the appropriate module
	"""
        Dist.__init__(self,space)
        if check:
            base = space.base_ring()
            try:
                M = len(moments)
            except TypeError:
                M = 1
                moments = [moments]
            moments = space.approx_module(M)(moments)
        self.moments = moments

    cdef Dist_vector _new_c(self):
        cdef Dist_vector ans = PY_NEW(Dist_vector)
        ans._parent = self._parent
        return ans

    def _repr_(self):
        """
	Displays the moments of the distribution
	"""
        self.normalize()
        return repr(self.moments)

    def moment(self,n):
        """
	Returns the n-th moment
	"""
        return self.moments[n]

    cpdef ModuleElement _add_(self, ModuleElement _right):
        cdef Dist_vector ans = self._new_c()
        cdef Dist_vector right = _right
        if len(self.moments) == len(right.moments):
            ans.moments = self.moments + right.moments
        elif len(self.moments) < len(right.moments):
            ans.moments = self.moments + right.moments[:len(self.moments)]
        else:
            ans.moments = self.moments[:len(right.moments)] + right.moments
        return ans

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        cdef Dist_vector ans = self._new_c()
        cdef Dist_vector right = _right
        if len(self.moments) == len(right.moments):
            ans.moments = self.moments - right.moments
        elif len(self.moments) < len(right.moments):
            ans.moments = self.moments - right.moments[:len(self.moments)]
        else:
            ans.moments = self.moments[:len(right.moments)] - right.moments
        return ans

    cpdef ModuleElement _lmul_(self, RingElement right):
        cdef Dist_vector ans = self._new_c()
        ans.moments = self.moments * right
        return ans

    def num_moments(self):
        """
	Returns the number of moments of the distribution
	"""
        return len(self.moments)

    cdef int _cmp_c_impl(left, Element right) except -2:
        return cmp(left.moments, right.moments)

    def zero(self):
        cdef Dist_vector ans = self._new_c()
        ans.moments = 0 * self.moments # could make this faster
        return ans

    cpdef normalize(self):
        """
	reduced modulo Fil^N -- that is the i-th moments is reduced modulo p^(N-i)
	"""
        p=self.parent()._p
        assert self.valuation() >= 0, "moments not integral in normalization"
        V = self.moments.parent()
        n = self.num_moments()
        self.moments = V([self.moment(i)%(p**(n-i)) for i in range(n)])
        return self

    def change_precision(self, M):
        """
	Only holds on to M moments
	"""
        assert M<=self.num_moments(),"not enough moments"

        cdef Dist_vector ans = self._new_c()
        ans.moments = self.moments[:M]
        return ans

    def solve_diff_eqn(self):
        # assert self.moments[0][0]==0, "not total measure zero"
        # print "result accurate modulo p^",self.moment(0).valuation(self.p)
        #v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
        M = self.num_moments()
        K = self.parent().base_ring().fraction_field()
        V = self.moments.parent()
        v = [K(0) for i in range(M)]
        for m in range(1,M):
            scalar = K(self.moment(m) / m) # division can take us out of the base ring.
            for j in range(m-1,M):
                v[j] += binomial(j,m-1) * bernoulli(j-m+1) * scalar
        cdef Dist_vector ans = self._new_c()
        ans.moments = V(v)
        return ans

    def lift(self):
        """
	Increases the number of moments by 1
	"""
        cdef Dist_vector ans = self._new_c()
        n = len(self.moments)
        ans.moments = self.parent().approx_module(n+1)(list(self.moments) + [0])
        return ans

cdef class Dist_long(Dist):
    def __init__(self, space, moments, check=True):
        Dist.__init__(self, space)
        p = space._p
        cdef int i
        if check:
            if len(moments) > 100 or 7*p**len(moments) > ZZ(2)**(4*sizeof(long) - 1): # 6 is so that we don't overflow on gathers
                raise ValueError("moments too long")
        for i in range(len(moments)):
            self.moments[i] = int(moments[i])
        self.prec = len(moments)
        #gather = 2**(4*sizeof(long)-1) // p**len(moments)
        #if gather >= len(moments):
        #    gather = 0
        #self._gather = gather

    cdef Dist_long _new_c(self):
        cdef Dist_long ans = PY_NEW(Dist_long)
        ans._parent = self._parent
        return ans

    def _repr_(self):
        self.normalize()
        return "(" + ", ".join([repr(self.moments[i]) for i in range(self.prec)]) + ")"

    def moment(self, _n):
        cdef int n = _n
        if n < 0:
            n += self.prec
        if n < 0 or n >= self.prec:
            raise IndexError("list index out of range")
        return self.moments[n]

    cpdef ModuleElement _add_(self, ModuleElement _right):
        cdef Dist_long ans = self._new_c()
        cdef Dist_long right = _right
        ans.prec = self.prec if self.prec < right.prec else right.prec
        cdef int i
        for i in range(ans.prec):
            ans.moments[i] = self.moments[i] + right.moments[i]
        return ans

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        cdef Dist_long ans = self._new_c()
        cdef Dist_long right = _right
        ans.prec = self.prec if self.prec < right.prec else right.prec
        cdef int i
        for i in range(ans.prec):
            ans.moments[i] = self.moments[i] - right.moments[i]
        return ans

    cpdef ModuleElement _lmul_(self, RingElement _right):
        cdef Dist_long ans = self._new_c()
        ans.prec = self.prec
        cdef long scalar
        if PY_TYPE_CHECK(_right, Integer):
            

cdef class WeightKAction_vector(Action):
    def __init__(self, Dk, character):
        self._k = Dk._k
        self._character = character
        self._p = Dk._p
        if character is None:
            self._Np = Dk._p
        else:
            self._Np = Dk._p * character.conductor()
        self._actmat = {}
        self._maxprecs = {}
        Action.__init__(self, M2Z, Dk, False, operator.mul)

    def clear_cache(self):
        self._actmat = {}
        self._maxprecs = {}

    cpdef acting_matrix(self, g, M):
        g.set_immutable() # a bit sketchy
        if not self._maxprecs.has_key(g):
            self._maxprecs[g] = M
            A = self._compute_acting_matrix(g, M)
            self._actmat[g] = {M:A}
            return A
        else:
            mats = self._actmat[g]
            if mats.has_key(M):
                return mats[M]
            maxprec = self._maxprecs[g]
            if M < maxprec:
                A = mats[maxprec][:M,:M] # submatrix; might want to reduce precisions
                mats[M] = A
                return A
            if M < 2*maxprec:
                maxprec = 2*maxprec
            else:
                maxprec = M
            self._maxprecs[g] = maxprec
            mats[maxprec] = self._compute_acting_matrix(g, maxprec) # could lift from current maxprec
            if M == maxprec:
                return mats[maxprec]
            A = mats[maxprec][:M,:M] # submatrix; might want to reduce precisions
            mats[M] = A
            return A

    cpdef _compute_acting_matrix(self, g, M):
        a = g[0,0]
        b = g[0,1]
        c = g[1,0]
        d = g[1,1]
        if a*d == b*c:
            raise ValueError("zero determinant")
        if self._p.divides(a):
            raise ValueError("p divides a")
        if not self._Np.divides(c):
            raise ValueError("Np does not divide c")
        k = self._k
        K = self.S.base_ring().fraction_field()
        padic = isinstance(K, pAdicGeneric)
        if padic:
            R = PowerSeriesRing(K, 'y', default_prec = M)
            OK = K.integer_ring()
        else:
            R = PowerSeriesRing(QQ, 'y', default_prec = M)
            OK = ZZ
        y = R.gen()
        scale = (b+d*y)/(a+c*y)
        t = (a+c*y)**k # will already have precision M
        cdef Matrix B = matrix(OK,M,M)
        q = self.S._p**M
        for c in range(M):
            for r in range(M):
                if padic:
                    B.set_unsafe(r, c, t[r].add_bigoh(M))
                else:
                    B.set_unsafe(r, c, t[r] % q)
            t *= scale
        return B

    cpdef _call_(self, _v, g):
        # if g is a matrix it needs to be immutable
        # hashing on arithmetic_subgroup_elements is by str
        cdef Dist_vector v = <Dist_vector?>_v
        cdef Dist_vector ans = v._new_c()
        ans.moments = v.moments * self.acting_matrix(g, len(v.moments))
        return ans




# @cached_function
def form_acting_matrix_on_dist(p,M,k,a,b,c,d):
    """
    Forms a large M x M matrix say G such that if v is the vector of
    moments of a distribution mu, then v*G is the vector of moments of
    mu|[a,b;c,d]
    """

    print("Checking...")
    print(a,b,c,d)
    print(p)

    assert (a%p != 0) and (c%p == 0), "acting by bad matrix"

    R=PowerSeriesRing(QQ,'y',default_prec=M)
    y=R.gen()

    scale=(b+d*y)/(a+c*y)
    t=((a+c*y)**k).truncate(M)

    A = []
    for i in range(0,M):
        temp1=t.list();
        d=len(temp1)
        for j in range(d,M):
            temp1 = temp1 + [0]
        #while len(temp1)>M:
        #    temp1.pop()
        A = A + [temp1]
        t=(t*scale).truncate(M)
    q=p**M
    B=Matrix(QQ,A).transpose()
    for r in range(0,M):
        for c in range(0,M):
            #B[r,c]=B[r,c]%(p**(M-c))
            B[r,c]=B[r,c]%(q)

    return B

#@cached_function
#def eta(Dk, i, M):
#    """
#    Helper function in solving the difference equation -- see Lemma 4.4 of [PS]
#    """
#    v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
#    return dist(p,k,vector(v))
#
#def random_dist(p,k,M):
#    """
#    Returns a random distribution with prime p, weight k, and M moments
#    """
#    moments=vector([ZZ(floor(random()*(p**M))) for i in range(1,M+1)])
#
#    return dist(p,k,moments)


