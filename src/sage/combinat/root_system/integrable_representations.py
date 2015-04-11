"""
Integrable Representations of Affine Lie Algebras
"""

#*****************************************************************************
#  Copyright (C) 2014 Daniel Bump <bump at match.stanford.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.category_object import CategoryObject
#from sage.structure.parent import Parent
from sage.categories.modules import Modules
from sage.combinat.root_system.cartan_type import CartanType
from sage.rings.all import ZZ, QQ
from sage.misc.all import cached_method
from root_space import RootSpace
from weight_space import WeightSpace
from sage.functions.other import floor

# TODO: Make this a proper parent and implement actions
class IntegrableRepresentation(CategoryObject, UniqueRepresentation):
    r"""
    An irreducible highest weight representation of an affine Lie algebra.

    INPUT:

    - ``Lam`` -- a dominant weight in an extended weight lattice
      of affine type

    OPTIONAL:

    - ``depth`` -- a parameter indicating how far to push computations

    REFERENCES:

    .. [Kac] Kac, *Infinite-dimensional Lie algebras*, Third Edition.
       Cambridge, 1990.

    .. [KMPS] Kass, Moody, Patera and Slansky, Affine Lie algebras,
       weight multiplicities, and branching rules. Vols. 1, 2. University of
       California Press, Berkeley, CA, 1990.

    .. [KacPeterson] Kac and Peterson. Infinite-dimensional Lie algebras, theta
       functions and modular forms. Adv. in Math. 53 (1984), no. 2, 125-264.

    If `\Lambda` is a dominant integral weight for an affine root system,
    there exists a unique integrable representation of highest weight
    `\Lambda`. If `\mu` is another weight such that `\Lambda - \mu` is in
    the root lattice, then multiplicity of `\mu` in this representation will
    be denoted `m(\mu)`.

    Let `\delta` be the nullroot. Then for fixed `\mu` the function
    `m(\mu - k\delta)` is a monotone increasing function of `\mu`. It is
    useful to take `\mu` to be such that this function is nonzero if and
    only if `k \geq 0`. Therefore we make the following definition.
    If `\mu` is such that `m(\mu) \neq 0` but `m(\mu + \delta) = 0`
    then `\mu` is called *maximal*.

    Since `\delta` is fixed under the action of the affine Weyl group,
    and since the weight multiplicities are Weyl group invariant, the
    *string function* `k \mapsto m(\mu-k\delta)` is unchanged if `\mu`
    is replaced by an equivalent weight. Therefore in tabulating the
    string functions, we may assume that `\mu` is dominant. There are
    only a finite number of dominant maximal weights.

    Since every nonzero weight multiplicity appears in the string
    `\mu - k\delta` for one of the finite number of dominant maximal
    weights `\mu`, it is important to be able to compute these. We may
    do this as follows.

    EXAMPLE::

         sage: Lambda = RootSystem(['A',3,1]).weight_lattice(extended=true).fundamental_weights()
         sage: IntegrableRepresentation(Lambda[1]+Lambda[2]+Lambda[3]).strings()
         3*Lambda[2] - delta: 3 21 107 450 1638 5367 16194 45687 121876 310056 757056 1783324
         2*Lambda[0] + Lambda[2]: 4 31 161 665 2380 7658 22721 63120 166085 417295 1007601 2349655
         Lambda[1] + Lambda[2] + Lambda[3]: 1 10 60 274 1056 3601 11199 32354 88009 227555 563390 1343178
         Lambda[0] + 2*Lambda[3]: 2 18 99 430 1593 5274 16005 45324 121200 308829 754884 1779570
         Lambda[0] + 2*Lambda[1]: 2 18 99 430 1593 5274 16005 45324 121200 308829 754884 1779570
         sage: Lambda = RootSystem(['D',4,1]).weight_lattice(extended=true).fundamental_weights()
         sage: IntegrableRepresentation(Lambda[0]+Lambda[1]).strings()                        # long time
         Lambda[3] + Lambda[4] - delta: 3 25 136 590 2205 7391 22780 65613 178660 463842 1155717 2777795
         Lambda[0] + Lambda[1]: 1 10 62 293 1165 4097 13120 38997 109036 289575 735870 1799620

    In this example, we construct the extended weight lattice of Cartan
    type `A_3^{(1)}`, then define ``Lambda`` to be the fundamental
    weights. We find there are 5 maximal dominant weights in irreducible
    representation of highest weight ``Lambda[1]+Lambda[2]+Lambda[3]``,
    and we determine their string functions.

    It was shown by Kac and Peterson that each string function is the
    set of Fourier coefficients of a modular form.

    Every weight `\mu` such that the weight multiplicity `m(\mu)` is
    nonzero has the form

      .. MATH::

          \Lambda - n_0 \alpha_0 - n_1 \alpha_1 - \cdots,

    where the `n_i` are nonnegative integers. This is represented internally
    as a tuple ``(n0, n1, n2, ...)``. If you want an individual multiplicity
    you use the method ``m`` and supply it with this tuple::

        sage: Lambda = RootSystem(['C',2,1]).weight_lattice(extended=true).fundamental_weights()
        sage: v = IntegrableRepresentation(2*Lambda[0]); v
        Integrable representation of ['C', 2, 1] with highest weight 2*Lambda[0]
        sage: v.m((3,5,3))
        18

    The :class:`IntegrableRepresentation` class has methods :meth:`to_weight`
    and :meth:`from_weight` to convert between this internal representation
    and the weight lattice::

        sage: delta = v._Q.null_root()
        sage: v.to_weight((4,3,2))
        -3*Lambda[0] + 6*Lambda[1] - Lambda[2] - 4*delta
        sage: v.from_weight(-3*Lambda[0] + 6*Lambda[1] - Lambda[2] - 4*delta)
        (4, 3, 2)

    To get more values, use the depth parameter::

        sage: L0 = RootSystem(["A",1,1]).weight_lattice(extended=true).fundamental_weight(0); L0
        Lambda[0]
        sage: IntegrableRepresentation(4*L0).strings(depth=20)
        4*Lambda[1] - 2*delta: 1 2 6 11 23 41 75 126 215 347 561 878 1368 2082 3153 4690 6936 10121 14677 21055
        2*Lambda[0] + 2*Lambda[1] - delta: 1 2 5 10 20 36 66 112 190 310 501 788 1230 1880 2850 4256 6303 9222 13396 19262
        4*Lambda[0]: 1 1 3 6 13 23 44 75 131 215 354 561 889 1368 2097 3153 4712 6936 10151 14677
    """
    def __init__(self, Lam):
        """
        Initialize ``self``.
        """
        CategoryObject.__init__(self, base=ZZ, category=Modules(ZZ))

        self._Lam = Lam
        self._P = Lam.parent()
        self._Q = self._P.root_system.root_lattice()

        self._cartan_matrix = self._P.root_system.cartan_matrix()
        self._cartan_type = self._P.root_system.cartan_type()
        self._classical_rank = self._cartan_type.classical().rank()
        self._index_set = self._P.index_set()
        self._index_set_classical = self._cartan_type.classical().index_set()

        cmi = self._cartan_type.classical().cartan_matrix().inverse()
        self._cminv = {}
        for i in self._index_set_classical:
            for j in self._index_set_classical:
                self._cminv[(i,j)] = cmi[i-1][j-1]

        self._smat = {}
        for i in self._index_set:
            for j in self._index_set:
                self._smat[(i,j)] = -self._cartan_matrix[i,j]
            self._smat[(i,i)] += 1

        self._shift = {}
        alphacheck = self._Q.simple_coroots()
        for i in self._index_set:
            self._shift[i] = self._Lam.scalar(alphacheck[i])

        self._ddict = {}
        self._mdict = {tuple(0 for i in self._index_set): 1}

        Lam_rho = self._Lam + self._P.rho()
        self._den0 = Lam_rho.symmetric_form(Lam_rho)

        # Coerce a classical root into Q
        from_cl_root = lambda h: self._Q._from_dict(h._monomial_coefficients)
        self._classical_roots = [from_cl_root(al)
                                 for al in self._Q.classical().roots()]
        self._classical_positive_roots = [from_cl_root(al)
                                          for al in self._Q.classical().positive_roots()]

        self._eps = {i: self._cartan_type.a()[i] / self._cartan_type.dual().a()[i]
                     for i in self._index_set}

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return "Integrable representation of %s with highest weight %s"%(self._cartan_type, self._Lam)

    def inner_qq(self, qelt1, qelt2):
        """
        Symmetric form between two elements of the root lattice.

        .. WARNING:

            If ``qelt1`` or ``qelt1`` accidentally gets coerced into
            the extended weight lattice, this will return an answer,
            and it will be wrong. To make this code robust, parents
            should be checked. This is not done since in the application
            the parents are known, so checking would unnecessarily slow
            us down.
        """
        mc1 = qelt1.monomial_coefficients()
        mc2 = qelt2.monomial_coefficients()
        zero = ZZ.zero()
        return sum(mc1.get(i, zero) * mc2.get(j, zero)
                   * self._cartan_matrix[i,j] / self._eps[i]
                   for i in self._index_set for j in self._index_set)

    def inner_pq(self, pelt, qelt):
        """
        Symmetric form between an element of the weight and root lattices

        .. WARNING:

            If ``qelt`` accidentally gets coerced into the extended weight
            lattice, this will return an answer, and it will be wrong. To make
            this code robust, parents should be checked. This is not done
            since in the application the parents are known, so checking would
            unnecessarily slow us down.
        """
        mcp = pelt.monomial_coefficients()
        mcq = qelt.monomial_coefficients()
        zero = ZZ.zero()
        return sum(mcp.get(i, zero) * mcq.get(i, zero) / self._eps[i]
                   for i in self._index_set)

    def to_weight(self, n):
        """
        Return the weight associated to the tuple ``n``.

        If ``n`` is the tuple `(n_1, n_2, \ldots)`, then the associated
        weight is `\Lambda - \sum_i n_i alpha_i`, where `\Lambda` is the
        weight of the representation.
        """
        alpha = self._P.simple_roots()
        I = self._index_set
        return self._Lam - self._P.sum(ZZ(val) * alpha[I[i]]
                                       for i,val in enumerate(n))

    # This recieves a significant number of calls.
    # TODO: Try to reduce the number of such calls.
    def _from_weight_helper(self, mu):
        """
        It is assumeed that ``mu`` is in the root lattice.
        returns ``(n[0], n[1], ...)`` such that ``mu = \sum n[i]*alpha[i]``
        """
        mu = self._P(mu)
        n0 = mu.monomial_coefficients().get('delta', 0)
        mu0 = mu - n0*self._P.simple_root(self._cartan_type.special_node())
        ret = [ZZ(n0)]
        mc_mu0 = mu0.monomial_coefficients()
        zero = ZZ.zero()
        for i in self._index_set_classical:
            ret.append( ZZ.sum(self._cminv[(i,j)] * mc_mu0.get(j, zero)
                               for j in self._index_set_classical) )
        return tuple(ret)
    
    def from_weight(self, mu):
        """
        Return ``(n[0], n[1], ...)`` such that ``mu = Lam - \sum n[i]*alpha[i]``
        """
        return self._from_weight_helper(self._Lam - mu)

    def s(self, n, i):
        """
        Implements the `i`-th simple reflection in the internal
        representation of weights by tuples.
        """
        ret = list(n) # This makes a copy
        ret[i] += self._Lam.monomial_coefficients().get(i, ZZ.zero())
        row = self._cartan_matrix[i]
        ret[i] -= sum(n[j] * row[j] for j in self._index_set)
        return tuple(ret)

    def to_dominant(self, n):
        if self._ddict.has_key(n):
            return self._ddict[n]
        mc = self.to_weight(n).monomial_coefficients()
        for i in self._index_set:
            if mc.get(i, 0) < 0:
                m = self.s(n, i)
                v = self.to_dominant(m)
                self._ddict[n] = v
                return v
        return n

    def freudenthal_roots_imaginary(self, nu):
        r"""
        It is assumed that ``nu`` is in `Q`. Returns the set of imaginary roots
        `\alpha \in \Delta^+` such that `nu - \alpha \in Q^+`.

        Lambda = self._P.fundamental_weights()
        kp = min(ZZ(nu.symmetric_form(Lambda[i])) for i in self._index_set)
        delta = self._Q.null_root()
        return [u*delta for u in [1..kp]]
        """
        l = self._from_weight_helper(nu)
        a = self._cartan_type.a()
        kp = min(floor(l[i] / a[i]) for i in self._index_set)
        delta = self._Q.null_root()
        return [u*delta for u in range(1, kp+1)]

    def freudenthal_roots_real(self, nu):
        r"""
        It is assumed that ``nu`` is in `Q`. Returns the set of real positive
        roots `\alpha \in \Delta^+` such that `nu - \alpha \in Q^+`.
        """
        ret = []
        for al in self._classical_positive_roots:
            if all(x >= 0 for x in self._from_weight_helper(nu-al)):
                ret.append(al)
        for al in self._classical_roots:
            for ir in self.freudenthal_roots_imaginary(nu-al):
                ret.append(al+ir)
        return ret

    def freudenthal_roots(self, nu):
        r"""
        It is assumed that ``nu`` is in `Q`. Returns the set of real roots
        `\alpha \in \Delta^+` such that `nu - \alpha \in Q^+`.

        This code is not called in the main algorithm.
        """
        ret = []
        for alpha in self.freudenthal_roots_imaginary(nu):
            ret.append(alpha)
        for alpha in self.freudenthal_roots_real(nu):
            ret.append(alpha)
        return ret

    # TODO: (mostly to self) try to optimize starting here
    def _freudenthal_accum(self, nu, al):
        ret = 0
        k = 1
        while all(val >= 0 for val in self._from_weight_helper(self._Lam - nu - k*al)):
            mk = self.m( self.from_weight(nu + k*al) )
            ip = self.inner_pq(nu + k*al, al)
            ret += 2*mk*ip
            k += 1
        return ret

    def m_freudenthal(self, n):
        """
        Compute the weight multiplicity using the Freudenthal multiplicity formula.
        """
        if min(n) < 0:
            return 0
        mu = self.to_weight(n)
        I = self._index_set
        al = self._Q._from_dict({I[i]: val for i,val in enumerate(n) if val},
                                remove_zeros=False)
        den = 2*self.inner_pq(self._Lam+self._P.rho(), al) - self.inner_qq(al,al)
        num = 0
        for al in self.freudenthal_roots_real(self._Lam - mu):
            num += self._freudenthal_accum(mu, al)
        for al in self.freudenthal_roots_imaginary(self._Lam - mu):
            num += self._classical_rank*self._freudenthal_accum(mu, al)
        if den == 0 and num == 0:
            print "m_freudenthal","m: n=%s, num=den=0"%n.__repr__()
        try:
            return ZZ(num / den)
        except:
            return None

    def m(self, n):
        """
        Return the multiplicity.
        """
        # TODO: Make this non-recursive by implementing our own stack
        if self._mdict.has_key(n):
            return self._mdict[n]
        elif self._ddict.has_key(n):
            self._mdict[n] = self.m(self._ddict[n])
        m = self.to_dominant(n)
        if self._mdict.has_key(m):
            return self._mdict[m]
        ret = self.m_freudenthal(m)
        if ret is not None:
            self._mdict[n] = ret
        else:
            print "m: error - failed to compute m%s"%n.__repr__()
        return ret

    # FIXME: Make this generate itself without having needing to be called by string()
    #@lazy_attribute
    def dominant_maximal(self):
        """
        Return the finite set of dominant maximal weights.
        """
        ret = set()
        delta = self._Q.null_root()
        for x in self._ddict.values():
            if self.m(x) > 0:
                if min(x) == 0:
                    ret.add(x)
                else:
                    y = self.from_weight(self.to_weight(x) + delta)
                    if self.m(y) == 0:
                        ret.add(x)
        return [self.to_weight(x) for x in ret]

    def string(self, max_weight, depth=12):
        """
        Return the list of multiplicities `m(\Lambda - k \delta)` where
        `\Lambda` is ``max_weight`` and `k` runs from `0` to ``depth``.

        INPUT:

        - ``max_weight`` -- a dominant maximal weight
        - ``depth`` -- (default: 12) the maximum value of `k`
            """
        ret = []
        delta = self._Q.null_root()
        for k in range(depth):
            ret.append(self.m( self.from_weight(max_weight - k*delta) ))
        return ret

    def strings(self, depth=12):
        """
        Return the set of dominant maximal weights of ``self``, together
        with the string coefficients for each.
        """
        # FIXME: This call to string should not be necessary as it is
        #   highly redundant to generate the data for dominant_maximal
        self.string(self._Lam, depth)
        for max_weight in self.dominant_maximal():
            s = self.string(max_weight, depth)
            print "%s:"%max_weight,
            for j in s:
                print j,
            print

