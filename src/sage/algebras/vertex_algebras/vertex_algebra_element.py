r"""
Vertex Algebra Element Class

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element_wrapper import ElementWrapper
from sage.functions.other import binomial, factorial
from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.structure.element import parent
from sage.combinat.partition_tuple import PartitionTuples_level
from sage.combinat.partition import Partition


class VertexAlgebraElement(ElementWrapper):
    """
    Vertex Algebra base element class.
    """
    pass

class UniversalEnvelopingVertexAlgebraElement(VertexAlgebraElement):
    """
    Universal Enveloping Vertex Algebra base element class.
    """
    def _repr_(self):
        p = self.parent()
        if self == p.zero():
            return "0";
        coeff = list(self.value.monomial_coefficients().items())
        ret = ""
        for i in range(len(coeff)):
            if i > 0  and coeff[i][1] > 0:
                ret += "+"
            if coeff[i][1] < 0:
                ret += "-"

            if abs(coeff[i][1]) != 1:
                ret += "{}*".format(abs(coeff[i][1]))
            for idx,j in enumerate(coeff[i][0]):
                for k in j.to_list():
                    ret += "{}_".format(p._lca._repr_generator(p._lca.gen(idx)))
                    if p.is_graded():
                        ret+= "{}".format(1 -k-p._lca.gen(idx).degree())
                    else:
                        ret+= "({})".format(-k)
            ret+="|0>"
        return ret

    def _add_(self,right):
        return type(self)(self.parent(), self.value + right.value)

    def __nonzero__(self):
        return bool(self.value)

    def _sub_(self, right):
        return type(self)(self.parent(), self.value - right.value)

    def _neg_(self):
        return type(self)(self.parent(), -self.value)

    def _bracket_(self, other):
        # sum over monomials for other. Use non-commutative Wick formula to
        # Reduce to other being a generator or a derivative.
        # Use Skew-Symmetry to reduce to self being a generator or a
        # derivative.
        p = self.parent()
        ret = {}
        svmc = other.value.monomial_coefficients()
        for k in svmc.keys():
            c = svmc[k]
            i = next((j for j,x in enumerate(k) if x), None)
            if i is None:
                continue
            l = [Partition([]),]*p.ngens()
            l[i] = Partition([k[i].get_part(0)])
            if k == l:
                svmc2 = self.value.monomial_coefficients()
                for k2 in svmc2.keys():
                    c2 = svmc2[k2]
                    i2 = next((j for j,x in enumerate(k2) if x), None)
                    if i2 is None:
                        continue
                    l2 = [Partition([]),]*p.ngens()
                    l2[i2] = Partition([k2[i2].get_part(0)])
                    if k2 == l2:
                        #Here we only use T from LieConformalAlgebras
                        br = p._lca.gen(i2).T(l2[i2].size()-1) \
                            ._bracket_(p._lca.gen(i).T(l[i].size()-1))
                        for j in br.keys():
                            rec = ret.get(j,p.zero())
                            ret[j] = rec + factorial(l[i].size()-1)**(-1)*\
                            factorial(l2[i2].size()-1)**(-1)*c*c2*p(br[j])
                    else:
                        #skew-symmetry
                        #Here is the only place we use T from VertexAlgebra
                        #We need T^j (other_(n) self)
                        #Note however that these terms are lower in the PBW
                        #filtration.
                        br = p(k)._bracket_(p(k2))
                        for cl in range(max(br.keys())+1):
                            rec = ret.get(cl, p.zero())
                            ret[cl] = rec+sum(c*c2*(-1)**(j+1)\
                                      *p.base_ring()(factorial(j-cl))\
                                      .inverse_of_unit()\
                                      *br[j].T(j-cl)
                                      for j in br.keys() if j >= cl)
            else:
                #non-commutative Wick Formula
                sf = k.to_list()
                sf[i] = Partition(Partition(sf[i]).to_list()[1:])
                br = self._bracket_(p(l))
                for j in br.keys():
                    rec = ret.get(j, p.zero())
                    ret[j] = rec + c*br[j]._mul_(p(sf))
                    #integral term
                    br2 = br[j]._bracket_(p(sf))
                    for m in br2.keys():
                        rec = ret.get(j+m+1 ,p.zero())
                        ret[j+m+1] = rec + c*binomial(j+m+1,j)*br2[m]
                br = self._bracket_(p(sf))
                for j in br.keys():
                    rec = ret.get(j, p.zero())
                    ret[j] = rec + c*p(l)._mul_(br[j])
        return {k:ret[k] for k in ret.keys() if ret[k]}

    def _mul_(self,right):
        r"""
        returns the normally ordered product :`self``right`:
        This is where the magic happens
        """
        #Use quasi-associativity to reduce to self being a generator. This
        #uses T and bracket from Vertex Algebra, but for terms lower in the
        #PBW filtration
        p = self.parent()
        ret = p.zero()
        svmc = self.value.monomial_coefficients()
        for k in svmc.keys():
            c = svmc[k]
            i = next((j for j,x in enumerate(k) if x), None)
            if i is None:
                ret += c*right
            else:
                l = [Partition([]),]*p.ngens()
                l[i] = Partition([k[i].get_part(0)])
                if k == l:
                    ni = k[i].get_part(0)
                    rvmc = right.value.monomial_coefficients()
                    for k2 in rvmc.keys():
                        c2 = rvmc[k2]
                        i2 = next((j for j,x in enumerate(k2) if x), None)
                        if i2 is None:
                            ret += c*c2*p(l)
                        elif  i2 < i or (i2==i and k2[i2].get_part(0)>ni):
                            ni2 = k2[i2].get_part(0)
                            l2 = [Partition([]),]*p.ngens()
                            l2[i2] = Partition([ni2])
                            sf = k2.to_list()
                            sf[i2] = Partition(Partition(sf[i2]).to_list()[1:])
                            rec = p(l)._mul_(p(sf))
                            ret += c*c2*p(l2)._mul_(rec)
                            #Here's the only place where we need the
                            #structure constants. We are using that
                            #T^(n)a_{(-1)} = a_{(-1-n)} and the commutator
                            #formula. Given the implementation it's faster
                            #to compute the bracket in the other direction
                            #We only use T from LieConformalAlgebra and we
                            #only use the bracket in LieConformalAlgebra
                            br = p._lca.gen(i2).bracket(p._lca.gen(i))
                            if br:
                                ret -= c*c2*sum( binomial(-ni2,j)*\
                                factorial(ni+ni2+j-1)**(-1)*p(br[j].\
                                T(ni+ni2+j-1))._mul_(p(sf)) for j in br.keys() )
                        else:
                            sf = k2.to_list()
                            l2 = Partition(k2[i]).to_list()
                            l2.insert(0,ni)
                            sf[i] = Partition(l2)
                            ret += c*c2*p(sf)
                else:
                    #quasi-associativity
                    #This uses bracket in VerteAlgebra and T in Vertex
                    #algebra but it computes T in terms lower in the PBW
                    #filtration.

                    l = p(l)
                    sf = k.to_list()
                    sf[i] = Partition(Partition(sf[i]).to_list()[1:])
                    sf = p(sf)
                    rec = sf._mul_(right)
                    ret+= c*l._mul_(rec)
                    sfbr = sf.bracket(right)
                    lbr = l.bracket(right)
                    if sfbr:
                        ret+=c*sum(prod(Integer(1)/i for i in range(1 ,j+2 ))*
                            l.T(j+1)._mul_(sfbr[j]) for j in sfbr.keys())
                    if lbr:
                        ret+=c*sum(prod(Integer(1)/i for i in range(1 ,j+2 ))*
                            sf.T(j+1)._mul_(lbr[j]) for j in lbr.keys())
        return ret

    def T(self,n=1):
        r"""The `n`-th derivative of this element. If ``n`` is not provided
        it defaults to `1`.

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.inject_variables()
            Defining L
            sage: L.T()
            L_-3|0>
            sage: v = V.vacuum(); v.T()
            0
            sage: L.T(3)
            6*L_-5|0>
            sage: (L*L).T()
            2*L_-3L_-2|0>+L_-5|0>

        """
        if n==0:
            return self
        if n > 1:
            return self.T().T(n-1)

        p = self.parent()
        if self.is_zero() or self == p.vacuum():
            return p.zero()

        selfmon = self.monomials()
        if len(selfmon) > 1:
            return sum( m.T() for m in selfmon)

        #Now we just have a monomial to compute.
        k,c = list(self.value.monomial_coefficients().items())[0]
        kl = k.to_list()
        i = next((i for i,x in enumerate(k) if x))
        g = kl[i].pop(0)
        PT = PartitionTuples_level(k.level())
        return factorial(g-1)**(-1)*c*p(p._lca.gen(i).T(g))._mul_(
                    p(PT(kl)))+ factorial(g-1)**(-1)*c*\
                    p(p._lca.gen(i).T(g-1))._mul_(p(PT(kl)).T())


    def __getitem__(self, i):
        return self.value.__getitem__(i)

    def _acted_upon_(self, scalar, self_on_left=False):
        scalar_parent = parent(scalar)
        if scalar_parent != self.parent().base_ring():
            if self.parent().base_ring() \
                    .has_coerce_map_from(scalar_parent):
                scalar = self.parent().base_ring()( scalar )
            else:
                return None
        if self_on_left:
            return type(self)(self.parent(), self.value * scalar)
        return type(self)(self.parent(), scalar * self.value)

    def monomial_coefficients(self):
        r"""
        Return the monomial coefficients of ``self`` as a dictionary.
        """
        c = self.value.monomial_coefficients()
        p = self.parent()
        return { p(k) : c[k] for k in c.keys() }

    def is_monomial(self):
        """Return whether this element is a monomial"""
        return (len(self.monomial_coefficients()) == 1 or self.is_zero())

    def monomials(self):
        """
        The tuple of monomials in this element

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0; w = (L*L)*L;
            sage: w.monomials()
            (2*L_-3L_-3|0>, 4*L_-4L_-2|0>, 1/2*L_-6|0>, L_-2L_-2L_-2|0>)

        """
        return tuple(v[1]*v[0] for v in
                    self.monomial_coefficients().items())

    def index(self):
        r"""If this element is a monomial it returns the index parametrizing
        this element in the basis of this vertex algebra. If it is not a
        monomial it raises an error.

        EXAMPLES::

            sage: V = AffineVertexAlgebra(QQ, 'A1', 1); e = V.gen(0); f = V.gen(2);
            sage: e.index()
            ([1], [], [])
            sage: f.index()
            ([], [], [1])
            sage: e.T()*e*f
            E(alpha[1])_-3E(alphacheck[1])_-1|0>+E(alpha[1])_-4|0>+E(alpha[1])_-2E(alpha[1])_-1E(-alpha[1])_-1|0>
            sage: _.index()
            Traceback (most recent call last):
            ...
            ValueError: index can only be computed for monomials
            sage: e.T()*(e*f)
            E(alpha[1])_-2E(alpha[1])_-1E(-alpha[1])_-1|0>
            sage: _.index()
            ([2, 1], [], [1])

        """
        if self.is_zero():
            return None
        if not self.is_monomial():
            raise ValueError ("index can only be computed for"
                     " monomials")
        return list(self.value.monomial_coefficients().keys())[0]

    def weight(self):
        r"""The conformal weight o this element. It raises an error if the
        element is not homogeneous

        EXAMPLES::

            sage: V = AffineVertexAlgebra(QQ, 'A1', 1); e = V.gen(0); f = V.gen(2);
            sage: (e*f).weight()
            2
            sage: V.vacuum().weight()
            0
            sage: V.zero().weight()
            +Infinity
            sage: f*e
            E(alpha[1])_-1E(-alpha[1])_-1|0>-E(alphacheck[1])_-2|0>
            sage: (f*e).weight()
            2
            sage: (e*f + e).weight()
            Traceback (most recent call last):
            ...
            ValueError: E(alpha[1])_-1E(-alpha[1])_-1|0>+E(alpha[1])_-1|0> is not homogeneous!

        """
        p = self.parent()
        if not p.is_graded():
            raise ValueError("Weight is only defined in H-graded"
                " vertex algebras")
        if self.is_zero():
            return Infinity
        ls = []
        for idx in self.value.monomial_coefficients().keys():
            ret = sum(len(idx[i])*(p._lca.gen(i).degree()-1) + idx[i].size()
                for i in range(p.ngens()))
            ls.append(ret)
        if ls[1:] == ls[:-1]:
            return ls[0]
        raise ValueError("{} is not homogeneous!".format(self))

    def _li_filtration_monomial_degree(self):
        r"""The Li filtration degree of this non-zero monomial"""
        p = list(self.value.monomial_coefficients().keys())[0]
        return p.size() - sum(j.length() for j in p.components())

    def li_filtration_leading_terms(self):
        r"""Returns the leading terms of this element with respect to the Li
        filtration

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2);
            sage: V.find_singular(6)
            [L_-2L_-2L_-2|0>-33/8*L_-4L_-2|0>+93/64*L_-3L_-3|0>-27/16*L_-6|0>]
            sage: v = _[0]
            sage: v
            L_-2L_-2L_-2|0>-33/8*L_-4L_-2|0>+93/64*L_-3L_-3|0>-27/16*L_-6|0>
            sage: v.li_filtration_leading_terms()
            L_-2L_-2L_-2|0>

        """
        if self.is_zero():
            return self
        lt = [(m._li_filtration_monomial_degree(),m) for m in
              self.monomials()]
        lideg = min(k[0] for k in lt)
        return sum(k[1] for k in lt if k[0]==lideg)

    def li_filtration_degree(self):
        r"""Let `F_pV` denote the Li filtration of this vertex algebra. This
        method returns the maximum `p` such that this element belongs in
        `F_pV`.

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.inject_variables()
            Defining L
            sage: V.find_singular(6)
            [L_-2L_-2L_-2|0>-33/8*L_-4L_-2|0>+93/64*L_-3L_-3|0>-27/16*L_-6|0>]
            sage: v = _[0]
            sage: v
            L_-2L_-2L_-2|0>-33/8*L_-4L_-2|0>+93/64*L_-3L_-3|0>-27/16*L_-6|0>
            sage: v.li_filtration_leading_terms()
            L_-2L_-2L_-2|0>
            sage: v.li_filtration_degree()
            0
            sage: w = (L*L).T(); w
            2*L_-3L_-2|0>+L_-5|0>
            sage: w.li_filtration_degree()
            1
            sage: (w*v).li_filtration_degree()
            1

        """
        if self.is_zero():
            return Infinity
        return min(m._li_filtration_monomial_degree() for m in
                   self.monomials())

        def is_singular(self):
            """
            Return whether this vector is a singular vector

            If `a \in V` is a vector in a finitely generated H-Graded
            vertex algebra, then `a` is singular if for each homogeneous
            vector `v in V` we have `v_n a = 0` whenever `n > 0`.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
                sage: v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                sage: v.is_singular()
                True
                sage: V = AffineVertexAlgebra(QQ, 'A1', 2); E = V.0
                sage: (E*E*E).is_singular()
                True

            """
            p = self.parent()
            try:
                weight = self.weight()
            except ValueError:
                raise ValueError("Couldn't compute weight of {}, "\
                                 "it's not homogeneous?".format(self))
            #it's cheaper to compute the whole bracket first
            for g in p.gens():
                gw = g.weight()
                br = g._bracket_(self)
                if not all(br.get(n+gw-1, p.zero()).is_zero() for
                           n in range(1,weight+2)):
                    return False
            return True

        def _action_from_partition_tuple(self,p,negative=True):
            """
            helper function to apply elements of a vertex algebra
            constructed from partitions

            INPUT:

            - ``p`` -- `PartitionTuple`. The level of ``p`` needs to
              equal the number of generators of the vertex algebra
            - ``negative`` -- boolean (default: `True`);

            OUTPUT: the result of repeatedly applying
            modes determined by ``p`` of the generators of
            this vertex algebra to the vector ``self``. By default
            negative modes are applied. Thus if `p = [[1]]` and `L` is
            the unique generator of `V`, the mode `L_{-1}` will be
            applied. If ``negative`` is `False`, non-negative modes are
            applied instead. Thus in the example above `L_0` will be
            applied.

            EXAMPLES:


                sage: V = AffineVertexAlgebra(QQ,'A2',1); f=V.2; f
                E(alpha[1] + alpha[2])_-1|0>
                sage: V.ngens()
                8
                sage: f._action_from_partition_tuple(PartitionTuple([[2,1],[3],[],[],[1],[2],[],[]]))
                E(alpha[2])_-2E(alpha[2])_-1E(alpha[1])_-3E(alpha[1] + alpha[2])_-1E(alphacheck[2])_-1E(-alpha[2])_-2|0>+E(alpha[2])_-2E(alpha[2])_-1E(alpha[1])_-3E(alpha[1] + alpha[2])_-2E(-alpha[2])_-2|0>+E(alpha[2])_-2E(alpha[2])_-1E(alpha[1])_-3E(alpha[1])_-3E(alphacheck[2])_-1|0>-E(alpha[2])_-2E(alpha[2])_-1E(alpha[1])_-4E(alpha[1])_-3|0>

            """
            ngens = self.parent().ngens()
            if len(p) != ngens:
                raise ValueError("p has to be a partition tuple of "
                                 "level {0}, got {1}".format(ngens,p))
            ret = self
            p = p.to_list()
            p.reverse()
            for j in range(ngens):
                p[j].reverse()
                g = self.parent()(self.parent().gen(ngens-j-1))
                for n in p[j]:
                    if negative:
                        ret = g.nmodeproduct(ret,-n)
                    else:
                        ret = g.nmodeproduct(ret, n-1 )
            return ret

    def filtered_degree(self):
        """
        The smallest space `F^p` in the Li filtration of this vertex
        algebra containing this element

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
            sage: L.li_filtration_degree()
            0
            sage: (L.T(2)*L.T()).li_filtration_degree()
            3

        """
        return max(m.weight() for m in self.monomial())

    def is_homogeneous(self):
        """
        Whether this element is homogeneous with respect to conformal
        weight

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
            sage: L.is_homogeneous()
            True
            sage: (L + L.T()).is_homogeneous()
            False

        """
        try:
            self.weight()
        except ValueError:
            return False
        return True

    def degrevlex_lm(self):
        r"""
        Returns the leading monomials with respect to the degrevlex
        monomial order. Let this vertex algebra `V` be generated by fields
        `a^1,\dots,a^n` of conformal weights `\Delta_1,\dots,\Delta_n`. A
        typical monomial element of `V` is written as a product of elements
        `a^i_{-j - \Delta_i}` for `i=1,\dots,n` and `j \geq 0`. We define a
        lexycografical order on
        these generators in by degree first, and then by the order in the
        list of generators, namely
        `a^i_{-j-\Delta_i} < a^{i'}_{-j' - \Delta_{i'}}` if
        `j+\Delta_i < j' + \Delta_{i'}` or if
        `j+\Delta_i = j' + \Delta_{i'}` and `i < i'`.

        This method then returns the leading monomial with respect to the
        weighted degree reverse lexycographical order.

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2);  V.inject_variables()
            Defining L
            sage: v = L*L*L + L.T()*L.T() + L.T(4)
            sage: v
            3*L_-3L_-3|0>+4*L_-4L_-2|0>+49/2*L_-6|0>+L_-2L_-2L_-2|0>
            sage: v.degrevlex_lm()
            49/2*L_-6|0>
            sage: V = AffineVertexAlgebra(QQ, 'A1', 1); e = V.gen(0); h = V.gen(1); f = V.gen(2);
            sage: v = e.T() + f*f; v
            E(alpha[1])_-2|0>+E(-alpha[1])_-1E(-alpha[1])_-1|0>
            sage: v.degrevlex_lm()
            E(-alpha[1])_-1E(-alpha[1])_-1|0>
            sage: v = f.T() + e*e; v.degrevlex_lm() == e*e
            False
            sage: v.degrevlex_lm()
            E(-alpha[1])_-2|0>

        """
        def dgrlcmp(a):
            pa = list(a.value.monomial_coefficients().keys())[0]
            bigest = max(p.get_part(0) for p in pa)
            pa = pa.to_exp(bigest)
            ret = [a.weight(),]
            for j in range(bigest):
                for l in range(len(pa)):
                    ret.append(pa[l][bigest-j-1])
            return ret
        if self.is_zero():
            return self
        return sorted(self.monomials(),key=dgrlcmp)[0]

    def PBW_filtration_degree(self):
        r"""If `G_pV` is the increasing PBW filtration in this universal
        enveloping vertex algebra, this method returns the minimal `p` such
        that this element belongs to `G_pV`

        EXAMPLES::

            sage: V = AffineVertexAlgebra(QQ, 'A1', 1); e = V.gen(0); h = V.gen(1); f = V.gen(2);
            sage: e.PBW_filtration_degree()
            1
            sage: (e*f.T()).PBW_filtration_degree()
            2
            sage: (e*f - f*e).PBW_filtration_degree()
            1
            sage: V.vacuum().PBW_filtration_degree()
            0

        """
        p = self.parent()
        if self.is_zero():
            return -1
        ls = []
        for p in self.value.monomial_coefficients().keys():
            ret = sum(j.length() for j in p.components())
            ls.append(ret)
        return max(ls)

    def _pbw_one_less(self):
        if self.is_zero():
            return (None, None, None)
        svmc = self.value.monomial_coefficients()
        if len(svmc) != 1:
            raise ValueError("{} is not a monomial.".format(self))
        k = list(svmc.keys())[0]
        c = svmc[k]
        p = self.parent()
        i = next((j for j,x in enumerate(k) if x), None)
        if i is None:
            return (p.vacuum(), p.vacuum(), c)
        l = [Partition([]), ]*p.ngens()
        l[i] = Partition([k[i].get_part(0)])
        sf = k.to_list()
        sf[i] = Partition(Partition(sf[i]).to_list()[1:])
        return (p(l),p(sf),c)


