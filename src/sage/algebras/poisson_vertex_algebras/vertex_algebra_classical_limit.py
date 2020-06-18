r"""
Vertex Algebra Classical Limit

AUTHORS:

- Reimundo Heluani (06-15-2020): Initial implementation.

"""


#******************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.poisson_vertex_algebras import PoissonVertexAlgebras
from sage.categories.vertex_algebras import VertexAlgebras
from sage.modules.with_basis.indexed_element import IndexedFreeModuleElement
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method
from sage.combinat.family import Family
from sage.misc.misc import repr_lincomb
from sage.combinat.partition import Partition
from sage.rings.all import QQ, ZZ
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.morphism import RingHomomorphism_im_gens
from sage.categories.homset import Hom

class SingularSupportCoverMorphism(RingHomomorphism_im_gens):
    def __init__(self, domain, codomain):
        H = Hom(domain,codomain)
        RingHomomorphism_im_gens.__init__(self, H, codomain.gens())

    def kernel(self, deg):
        source = self.domain().get_weight(deg)
        target = self.codomain().get_weight(deg)
        B = source.basis()
        linear = source.module_morphism(on_basis=lambda v: target.retract(
                                        self(B[v])),codomain=target)
        return self.domain().submodule([v.lift() for v in\
                                        linear.kernel_basis()])

class ClassicalLimitElement(IndexedFreeModuleElement):

    def _repr_(self):
        """
        A

        EXAMPLES::

            sage: FreeFermionsVertexAlgebra(AA,ngens=2).inject_variables()
            Defining psi_0, psi_1
            sage: psi_0*psi_1.T(3)
            6*psi_0_-1/2psi_1_-7/2|0>
            sage: VirasoroVertexAlgebra(QQ,1/2).inject_variables()
            Defining L
            sage: L.T(5)*L.T(2) +L*L
            240*L_-7L_-4|0> + L_-2L_-2|0>

        For non H-graded vertex algebras we use normal modes::

            sage: WeylVertexAlgebra(QQ,ngens=2).inject_variables()
            Defining alpha0, alpha1
            sage: alpha1.T()*(alpha1.T(3)*alpha0)
            6*alpha0_(-1)alpha1_(-4)alpha1_(-2)|0>
            sage: alpha1*alpha1*alpha0
            alpha0_(-1)alpha1_(-1)alpha1_(-1)|0> - 2*alpha1_(-2)|0>
        """
        if self.is_zero():
            return "0";
        p = self.parent()
        terms = [("".join(["".join([ "{}{}^{}*".format(p.variable_names()[j],\
                 i+p.gen(j).weight(),r) if r > 1\
            else "{}{}*".format(p.variable_names()[j],i+p.gen(j).weight(),r)\
            for i,r in reversed(list(enumerate(mu))) if r])\
            for j,mu in enumerate(k.to_exp())]), v) for k,v in \
            sorted(self.monomial_coefficients().items())]
        terms = [(m[:-1],v) for m,v in terms]
        terms = [(k,v) if k or v != 1 else ('1',1) for k,v in terms]
        return repr_lincomb(terms,strip_one = True)

    def _latex_(self):
        """
        A visual representation of this element.

        EXAMPLES::

            sage: V = WeylVertexAlgebra(QQ,ngens=2); V.inject_variables()
            Defining alpha0, alpha1
            sage: latex(alpha1*alpha1*alpha0)
            alpha0_{(-1)}alpha1_{(-1)}alpha1_{(-1)}|0\rangle> - 2alpha1_{(-2)}|0\rangle>
        """
        if self.is_zero():
            return "0";
        p = self.parent()
        terms = [("".join(["".join(["{}_{{{}}}^{{{}}}".format(
                p._latex_names()[j],-1+i+p.gen(j).weight(),r) if r >1\
                else "{}_{{{}}}".format(p._latex_names()[j],
               -1+i+p.gen(j).weight(),r)  for i,r in reversed(list(enumerate(
               mu))) if r]) for j,mu in enumerate(k)]), v)\
              for k,v in sorted(self.monomial_coefficients().items())]
        terms = [(k,v) if k or v != 1 else ('1',1) for k,v in terms]
        return repr_lincomb(terms,strip_one = True, is_latex=True)

    def index(self):
        return self.lift().index()

    def weight(self):
        if self.is_zero():
            from sage.rings.infinity import Infinity
            return +Infinity
        weights = [k.energy() for k,v in self._monomial_coefficients.items()]
        if weights[1:] == weights[:-1]:
            return weights[0]
        raise ValueError("{} is not homogeneous".format(self))

    def is_homogeneous(self):
        weights = [k.energy() for k,v in self._monomial_coefficients.items()]
        return weights[1:] == weights[:-1]

    def homogeneous_terms(self):
        if self.is_zero():
            return tuple([self])
        S = {}
        p = self.parent()
        for m in self.terms():
            w = m.weight()
            S[w] = S.get(w,p.zero()) + m
        return tuple(S.values())

    def degree(self):
        return max(m.weight() for m in self.monomials())

    def _li_filtration_monomial_degree(self):
        return self.lift()._li_filtration_monomial_degree()

    def lift(self):
        return self.parent()._ambient._from_dict(self._monomial_coefficients)

    def T(self, n=1):
        #This implementation now only works for universal enveloping
        if self.is_zero() or n == 0:
            return self
        if n > 1:
            return self.T().T(n-1)
        if self.is_monomial():
            p = self.parent()
            ds = self._li_filtration_monomial_degree()
            a,b,c = self._pbw_one_less()
            if b == p.one():
                if a == p.one():
                    return p.zero()
                pt = a.index().to_list()
                i = next((j for j,x in enumerate(pt) if x))
                c = c*(pt[i][0])
                pt[i][0] += 1
                pt = p._ambient.indices()(pt)
                if p._ambient in VertexAlgebras(p.base_ring()).Quotients():
                    ret = c*p._ambient.cover_algebra()(pt)
                    ret = p._ambient.retract(ret)
                    ret = sum([m for m in ret.terms() if\
                            m._li_filtration_monomial_degree() == ds+1],
                            p._ambient.zero())
                else:
                    ret = c*p._ambient(pt)
                return p._from_dict(ret._monomial_coefficients)
            return c*(a.T()*b + a*b.T())
        return sum(m.T() for m in self.terms())

    def _mul_(self,other):
        p = self.parent()
        if self.is_zero() or other.is_zero():
            return p.zero()
        if self.is_monomial() and other.is_monomial():
            a,b,c = self._pbw_one_less()
            if b == p.one():
                if a == p.one():
                    return c*other
                idx = self.index().to_list()
                odx,c2 = next(iter(other._monomial_coefficients.items()))
                i = next((j for j,x in enumerate(idx) if x))
                n = idx[i][0]
                sgn = p._parity_list[i]
                sgn = sum(p._parity_list[i]*p._parity_list[j]*len(odx[j])\
                          for j in range(i))
                sgn += p._parity_list[i]*len([j for j in odx[i] if j > n])
                sgn = (-1)**(sgn)
                odx = other.index().to_exp(n)
                odx[i][n-1] += 1
                if p._ambient in VertexAlgebras(p.base_ring()).Quotients():
                    pt = p._ambient.cover_algebra().indices()(tuple([
                                        Partition(exp=m) for m in odx]))
                    ret = sgn*c*c2*p._ambient.cover_algebra()(pt)
                    ret = p._ambient.retract(ret)
                else:
                    pt = p._ambient.indices()(tuple([Partition(exp=m)\
                                                     for m in odx]))
                    ret = sgn*c*c2*p._ambient(pt)

                ds = self._li_filtration_monomial_degree()
                do = other._li_filtration_monomial_degree()
                ret = sum([m for m in ret.terms() if\
                            m._li_filtration_monomial_degree() == ds+do],
                            p._ambient.zero())
                return p._from_dict(ret._monomial_coefficients)
            return c*a._mul_(b._mul_(other))
        return sum(i._mul_(j) for i in self.terms() for j in other.terms())

    def _bracket_(self,other):
        p = self.parent()
        if self.is_zero() or other.is_zero():
            return {}
        if self.is_monomial() and other.is_monomial():
            ds = self._li_filtration_monomial_degree()
            do = other._li_filtration_monomial_degree()
            br = self.lift()._bracket_(other.lift())
            ret = {k:sum([p._from_dict(t._monomial_coefficients) for t in
                    v.terms() if t._li_filtration_monomial_degree() ==\
                    ds + do - k],p.zero()) for k,v in br.items()}
            return {k:v for k,v in ret.items() if v}

        diclist = [i._bracket_(j) for i in self.terms() for j in other.terms()]
        ret = {}
        for d in diclist:
            for k in d:
                ret[k] = ret.get(k,p.zero()) + d[k]
        return {k:v for k,v in ret.items() if v}

    def is_monomial(self):
        return len(self._monomial_coefficients) == 1 or self.is_zero()

    def _pbw_one_less(self):
        """
        Expresses this element as a product of two elements lower in
        the PBW filtration.

        OUTPUT: a triple ``(a,b,c)`` such that ``self=c*a*b`` with
        ``c`` a scalar and ``a`` a monomial of PBW degree ``1``.
        """
        if self.is_zero():
            return (None, None, None)
        svmc = self.monomial_coefficients()
        if len(svmc) != 1:
            raise ValueError("{} is not a monomial.".format(self))
        k = next(iter(svmc))
        c = svmc[k]
        p = self.parent()
        i = next((j for j,x in enumerate(k) if x), None)
        if i is None:
            return (p.one(), p.one(), c)
        l = [[],]*len(k)
        l[i] = [k[i].get_part(0)]
        sf = k.to_list()
        sf[i] = sf[i][1:]
        return (p(l),p(sf),c)

    def _to_polynomial(self,ord=None, termorder='wdegrevlex'):
        p = self.parent()
        if ord == None or ord < self.degree():
            ord = self.degree()
        try:
            PR = p.jet_algebra(ord,termorder)
        except NotImplementedError:
            raise NotImplementedError("_to_polynomial is not implemented for "\
                                      "elements of {}".format(p))
        ret = PR.zero()
        weights = tuple(ZZ(g.weight()) for g in p.gens())
        for m,c in self.monomial_coefficients().items():
            lexp = [part.to_exp(ord - w +1) for part,w in zip(m,weights)]
            tuplexp = [(i,j,item) for i,part in enumerate(lexp) for j,item in\
                       enumerate(part)]
            if termorder == "wdegrevlex":
                k = tuple(t[2] for t in sorted(tuplexp, 
                          key = lambda x : (weights[x[0]]+x[1],x[0])))
            elif termorder == "wdeglex":
                k = tuple(t[2] for t in sorted(tuplexp, 
                          key = lambda x: (ord-weights[x[0]]-x[1], x[0])))
            elif termorder == "revlexwdeg":
                k = tuple(t[2] for t in sorted(tuplexp, 
                                               key = lambda x: (x[0],x[1])))
            elif termorder == "lexwdeg":
                k = tuple(t[2] for t in sorted(tuplexp,
                          key = lambda x : (x[0], ord - x[1])))
            else:
                raise NotImplementedError("termorder not implemented")

            ret += PR({k:c})
        return ret

    def _im_gens_(self, codomain, im_gens, base_map=None):
        if self.is_monomial():
            k,c = next(iter(self._monomial_coefficients.items()))
            from sage.functions.other import factorial
            from sage.misc.misc_c import prod
            return c*prod(g.T(n-1)/factorial(n-1) for i,g in enumerate(im_gens)\
                          for n in k[i])
        return sum(m._im_gens_(codomain, im_gens, base_map) for m in self.terms())

class VertexAlgebraClassicalLimit(CombinatorialFreeModule):
    def __init__(self, R, V, category=None):
        assert V in VertexAlgebras(R).Graded().FinitelyGenerated().WithBasis()
        from sage.algebras.vertex_algebras.universal_enveloping_vertex_algebra\
            import UniversalEnvelopingVertexAlgebra as UEA
        if not isinstance(V, UEA):
            if not V in VertexAlgebras(R).Quotients():
                raise NotImplementedError("classical limit is not implemented"\
                                          "for {}".format(V))
        defcat = PoissonVertexAlgebras(R).Graded().FinitelyGenerated().\
                 WithBasis()
        category = defcat.or_subcategory(category)
        try:
            names = V.variable_names()
        except ValueError:
            names = ['L%d'%i for i in range(V.ngens())]

        try:
            self._latex_names = V._latex_names
        except AttributeError:
            self._latex_names = names

        CombinatorialFreeModule.__init__(self, R, basis_keys=V._indices,
                                         element_class=ClassicalLimitElement,
                                         category=category, names=names)
        self._ambient = V
        self._parity_list = tuple(g.is_even_odd() for g in V.gens())

    def _element_constructor_(self,x):

        try:
            return CombinatorialFreeModule._element_constructor_(self,x)
        except (ValueError,TypeError):
            pass

        if x in self._ambient:
            return self._from_dict(x.monomial_coefficients())

        raise ValueError("Do not know how to convert {} into an element of {}"\
                         .format(x,self))

    def __contains__(self,x):
        if self.has_coerce_map_from(x.parent()):
            return True
        return super(VertexAlgebraClassicalLimit, self).__contains__(x)

    def basis(self):
        return Family(self._indices, self.monomial)

    def _repr_(self):
        return "The classical limit of {}".format(self._ambient)

    def dimension_at_weight(self,n):
        return self._ambient.dimension_at_weight(n)

    @cached_method
    def one(self):
        return self.monomial(self._ambient.vacuum().index())

    @cached_method
    def gens(self):
        return tuple(self._from_dict(g._monomial_coefficients) for g in \
                     self._ambient.gens())

    def ideal(self,gens,check=True):
        from .poisson_vertex_algebra_ideal import PoissonVertexAlgebraIdeal
        return PoissonVertexAlgebraIdeal(self,gens,check=check)

    def quotient(self, I):
        from .poisson_vertex_algebra_quotient import \
                                                   PoissonVertexAlgebraQuotient
        return PoissonVertexAlgebraQuotient(I,
                                           category=self.category().Quotients())

    def hilbert_series(self,ord=None):
        try:
            R = self.jet_algebra(ord)
        except NotImplementedError:
            from sage.arith.functions import lcm
            from sage.functions.other import floor
            weights = [g.weight() for g in self.gens()]
            if any([w not in QQ or w < 0 for w in weights]):
                raise NotImplementedError("hilbert_series is not "\
                                  "implemented for {}".format(self))
            if ord not in QQ or ord < 0:
                raise ValueError("ord must be a positive rational "\
                                 "number")
            l = lcm([g.weight().denominator() for g in self.gens()])
            if l==1:
                from sage.rings.power_series_ring import\
                                                    PowerSeriesRing
                q = PowerSeriesRing(ZZ,'q', default_prec=ord).gen()
                return sum(self.dimension_at_weight(n)*q**n for\
                           n in range(floor(ord))).O(floor(ord))
            else:
                from sage.rings.puiseux_series_ring import\
                                                  PuiseuxSeriesRing
                q = PuiseuxSeriesRing(ZZ,'q').gen()
                ord = floor(ord*l)
                f = sum(self.dimension_at_weight(n/l)*q**(n/l) for\
                        n in range(ord))
                return f.add_bigoh(ord/l)
        
        weights = tuple(g.degree() for g in R.gens())
        from sage.rings.power_series_ring import PowerSeriesRing
        q = PowerSeriesRing(ZZ, 'q', default_prec=ord).gen()
        return R.ideal(0).hilbert_series(grading=weights)(q)

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        from sage.algebras.vertex_algebras.universal_enveloping_vertex_algebra\
            import UniversalEnvelopingVertexAlgebra as UEA
        if not isinstance(self._ambient, UEA):
           raise NotImplementedError("morphisms are not implemented for "\
                                     "{}".format(self))
        return True

    @lazy_attribute
    def singular_support_cover(self):
        if self._ambient in VertexAlgebras(self.base_ring()).Quotients():
            P = self._ambient.singular_support()
            return SingularSupportCoverMorphism(P,self)

        from sage.categories.homset import End
        return End(self).identity()

    def jet_algebra(self,ord, termorder='wdegrevlex'):
        if self.is_super():
            raise NotImplementedError("jet_algebra is not implemented for "\
                                      "super Poisson vertex algebras.")

        from sage.algebras.vertex_algebras.universal_enveloping_vertex_algebra\
            import UniversalEnvelopingVertexAlgebra as UEA
        if not isinstance(self._ambient, UEA):
           raise NotImplementedError("jet_algebra is not implemented for "\
                                     "{}".format(self))

        if any (g.weight() not in ZZ or g.weight() < 0 for g in self.gens()):
            raise NotImplementedError("generators must have integral weight")

        names = self.variable_names()
        pgens = {g:(i,names[i]) for i,g in enumerate(self.gens())}
        
        vardict = {(v,j) : "{}_{}".format(pgens[v][1],j)
                   for v in self.gens() for j in range(v.weight(),ord+1)}

        if termorder == "wdegrevlex":
            varlist = tuple((vardict[v],v[1]) for v in sorted(vardict,
                key=lambda x:(x[1],pgens[x[0]][0])))
        elif termorder == "wdeglex":
            varlist = tuple((vardict[v],v[1]) for v in sorted(vardict,
                key=lambda x:(ord-x[1],pgens[x[0]][0])))
        elif termorder == "revlexwdeg":
            varlist = tuple((vardict[v],v[1]) for v in sorted(vardict,
                key=lambda x:(pgens[x[0]][0],x[1])))
            termorder = "wdegrevlex"
        elif termorder == "lexwdeg":
            varlist = tuple((vardict[v],v[1]) for v in sorted(vardict,
                key=lambda x:(pgens[x[0]][0],ord-x[1])))
            termorder = "wdeglex"
        else:
            raise ValueError("termorder needs to be one of 'wdegrevlex', "\
                             "'wdeglex', 'revlexwdeg' or 'lexwdeg'")

        from sage.rings.polynomial.polynomial_ring_constructor import\
                                                                PolynomialRing
        from sage.rings.polynomial.term_order import TermOrder
        varnames = tuple(a[0] for a in varlist)
        vardegs = tuple(a[1] for a in varlist)
        return PolynomialRing(self.base_ring(), len(varnames), varnames, 
                              order=TermOrder(termorder,vardegs))
   
