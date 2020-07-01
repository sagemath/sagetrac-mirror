r"""
Vertex Algebra Element

AUTHORS:

- Reimundo Heluani (2019-08-09): Initial implementation.
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

from sage.functions.other import binomial, factorial
from sage.rings.infinity import Infinity
from sage.modules.with_basis.indexed_element import IndexedFreeModuleElement
from sage.misc.misc import repr_lincomb

class UniversalEnvelopingVertexAlgebraElement(IndexedFreeModuleElement):
    """
    Universal Enveloping Vertex Algebra base element class.
    """
    def _repr_(self):
        """
        String representation of this element.

        EXAMPLES::

            sage: vertex_algebras.FreeFermions(AA,ngens=2).inject_variables()
            Defining psi_0, psi_1
            sage: psi_0*psi_1.T(3)
            6*psi_0_-1/2psi_1_-7/2|0>
            sage: vertex_algebras.Virasoro(QQ,1/2).inject_variables()
            Defining L
            sage: L.T(5)*L.T(2) +L*L
            L_-2L_-2|0> + 240*L_-7L_-4|0>

        For non H-graded vertex algebras we use normal modes::

            sage: vertex_algebras.Weyl(QQ,ngens=2).inject_variables()
            Defining alpha0, alpha1
            sage: alpha1.T()*(alpha1.T(3)*alpha0)
            6*alpha0_(-1)alpha1_(-4)alpha1_(-2)|0>
            sage: alpha1*alpha1*alpha0
            alpha0_(-1)alpha1_(-1)alpha1_(-1)|0> - 2*alpha1_(-2)|0>
        """
        if self.is_zero():
            return "0";
        p = self.parent()
        if p.is_graded():
            terms = [("".join(["".join(["{}_{}".format(p._lca.gen(j),\
            1-i-p._lca.gen(j).degree()) for i in mu])\
            for j,mu in enumerate(k)]), v) for k,v in \
            sorted(self.monomial_coefficients().items())]
        else:
            terms = [("".join(["".join(["{}_({})".format(p._lca.gen(j),-i)\
            for i in mu]) for j,mu in enumerate(k)]), v) for k,v in \
            self.monomial_coefficients().items()]
        return repr_lincomb(terms,strip_one = True,
                            repr_monomial = lambda s:s+"|0>")

    def _latex_(self):
        """
        A visual representation of this element.

        EXAMPLES::

            sage: V = vertex_algebras.Weyl(QQ,ngens=2); V.inject_variables()
            Defining alpha0, alpha1
            sage: latex(alpha1*alpha1*alpha0)
            alpha0_{(-1)}alpha1_{(-1)}alpha1_{(-1)}|0\rangle> - 2alpha1_{(-2)}|0\rangle>
        """
        if self.is_zero():
            return "0";
        p = self.parent()
        if p.is_graded():
            terms = [("".join(["".join(["{}_{{{}}}".format(p._lca.gen(j),\
            1-i-p._lca.gen(j).degree()) for i in mu]) for j,mu in enumerate(k)]\
            ), v) for k,v in sorted(self.monomial_coefficients().items())]
        else:
            terms = [("".join(["".join(["{}_{{({})}}".format(p._lca.gen(j),-i)\
            for i in mu]) for j,mu in enumerate(k)]), v) for k,v in \
            self.monomial_coefficients().items()]
        return repr_lincomb(terms,strip_one = True, is_latex=True,
                            repr_monomial = lambda s:s+r"|0\rangle>")

    def is_even_odd(self):
        """
        Return ``0`` if this element is `even` or ``1`` if it is `odd`.

        EXAMPLES::

            sage: V = vertex_algebras.NeveuSchwarz(QQbar,1/2); V
            The Neveu-Schwarz super vertex algebra of central charge 1/2 over Algebraic Field
            sage: V.inject_variables()
            Defining L, G
            sage: G.is_even_odd()
            1
            sage: G.is_odd()
            True
            sage: L.is_odd()
            False

        TESTS::

            sage: V = vertex_algebras.NeveuSchwarz(QQbar,1/2); v = V.an_element(); v
            |0> + 2*G_-3/2|0> + 3*L_-2|0> + L_-2G_-3/2|0>
            sage: v.is_odd()
            Traceback (most recent call last):
            ...
            ValueError: |0> + 2*G_-3/2|0> + 3*L_-2|0> + L_-2G_-3/2|0> is not homogeneous
        """
        p = self.parent()
        if self.is_zero() or self.monomials() == [p.vacuum()]:
            return 0
        if self.is_monomial():
            pt = self.index()
            return sum([len(mu) for i,mu in enumerate(pt)\
                        if p._lca.gen(i).is_even_odd()]) % 2
        parity_list = [m.is_even_odd() for m in self.monomials()]
        if parity_list[1:] != parity_list[:-1]:
            raise ValueError("{} is not homogeneous".format(self))
        return parity_list[0]

    def is_even(self):
        """
        Whether this element is `even`.

        EXAMPLES::

            sage: V = vertex_algebras.NeveuSchwarz(QQbar,1/2); V
            The Neveu-Schwarz super vertex algebra of central charge 1/2 over Algebraic Field
            sage: V.inject_variables()
            Defining L, G
            sage: G.is_even_odd()
            1
            sage: G.is_odd()
            True
            sage: V = vertex_algebras.Virasoro(AA,1); v = V.an_element(); v.is_odd()
            False
        """
        return not self.is_even_odd()

    def is_odd(self):
        """
        Whether this element is `odd`.

        EXAMPLES::

            sage: V = vertex_algebras.NeveuSchwarz(QQbar,1/2); V
            The Neveu-Schwarz super vertex algebra of central charge 1/2 over Algebraic Field
            sage: V.inject_variables()
            Defining L, G
            sage: G.is_even_odd()
            1
            sage: G.is_odd()
            True
            sage: V = vertex_algebras.Virasoro(AA,1); v = V.an_element(); v.is_odd()
            False
        """
        return not not self.is_even_odd()

    def _bracket_(self, other):
        """
        The OPE of these two elements.

        INPUT:

        - ``other`` -- an element of this vertex algebra

        OUTPUT:

        A dictionary with non-negative integer entries as keys.
        The value of the key ``n`` is the `n-th` product of ``self``
        and ``other``.

        EXAMPLES:

        We consider the bosonization of the charged Fermions::

            sage: V = vertex_algebras.FreeFermions(QQ, gram_matrix=Matrix([[0,1],[1,0]])); V
            The Free Fermions super vertex algebra with generators (psi_0_-1/2|0>, psi_1_-1/2|0>) over Rational Field
            sage: V.inject_variables()
            Defining psi_0, psi_1
            sage: alpha = psi_0*psi_1; alpha.bracket(alpha)
            {1: |0>}

        We check the standard conformal vector of the Free Fermions::

            sage: V = vertex_algebras.FreeFermions(QQ); V
            The Free Fermions super vertex algebra with generators (psi_-1/2|0>,) over Rational Field
            sage: V.inject_variables()
            Defining psi
            sage: L = 1/2*psi.T()*psi
            sage: Family(L.bracket(L))
            Finite family {0: psi_-5/2psi_-1/2|0>, 1: psi_-3/2psi_-1/2|0>, 3: 1/4*|0>}
            sage: L.bracket(L) == {0:L.T(), 1: 2*L, 3:1/4*V.vacuum()}
            True

        A similar computation for the free Boson::

            sage: V = vertex_algebras.FreeBosons(QQ); V
            The Free Bosons vertex algebra with generators (alpha_-1|0>,) over Rational Field
            sage: V.inject_variables()
            Defining alpha
            sage: L = 1/2*alpha*alpha
            sage: Family(L.bracket(L))
            Finite family {0: alpha_-2alpha_-1|0>, 1: alpha_-1alpha_-1|0>, 3: 1/2*|0>}
            sage: L.bracket(L) == {0:L.T(), 1: 2*L, 3: 1/2*V.vacuum()}
            True

        The topological twist of the `N=2` super vertex algebra::

            sage: V = vertex_algebras.N2(QQ,3/2); V
            The N=2 super vertex algebra of central charge 3/2 over Rational Field
            sage: V.inject_variables()
            Defining L, J, G1, G2
            sage: Lm = L - 1/2*J.T(); Family(Lm.bracket(Lm))
            Finite family {0: L_-3|0> - J_-3|0>, 1: 2*L_-2|0> - J_-2|0>}
            sage: Lm.bracket(Lm) == {0:Lm.T(), 1: 2*Lm}
            True

        TESTS::

            sage: V = vertex_algebras.Affine(QQ,'A1',1)
            sage: V.0.bracket(V.vacuum())
            {}
            sage: V.vacuum().bracket(V.an_element())
            {}
            sage: V.vacuum().bracket(V(0))
            {}
        """
        p = self.parent()
        pz = p.zero()

        if self.is_zero() or other.is_zero() or \
            self.monomials() == [p.vacuum()] or \
            other.monomials()== [p.vacuum()]:
            return {}

        if self.is_monomial() and other.is_monomial():
            a2,b2,c2 = other._pbw_one_less()
            if b2 == p.vacuum():
                a1,b1,c1 = self._pbw_one_less()
                if b1 == p.vacuum():
                    i1 = next((j for j,x in enumerate(a1.index()) if x))
                    n1 = a1.index()[i1][0]
                    i2 = next((j for j,x in enumerate(a2.index()) if x))
                    n2 = a2.index()[i2][0]
                    ret = {k: c1*c2*v.lift()/factorial(n1-1)/factorial(n2-1)\
                            for k,v in p._lca.gen(i1).T(n1-1).bracket(
                                                p._lca.gen(i2).T(n2-1)).items()}
                    return {k:v for k,v in ret.items() if v}
                #Use Skew-Symmetry
                if self.is_even_odd()*other.is_even_odd():
                    parsgn = -1
                else:
                    parsgn = 1
                br = other._bracket_(self)
                ret = {n: -sum([parsgn*(-1)**(k)*v.T(k-n)/factorial(k-n)\
                    for k,v in br.items() if k >= n]) for n in range(max(br,
                    default=0)+1)}
                return {k:v for k,v in ret.items() if v}
            #Use Non-Commutative Wick-Formula
            br1 = self._bracket_(a2)
            ret1 = {k: c2*v._mul_(b2) for k,v in br1.items()}
            br2 = self._bracket_(b2)
            parsgn = (-1)**(self.is_even_odd()*a2.is_even_odd())
            ret2 = {k: parsgn*c2*a2._mul_(v) for k,v in br2.items()}
            intdic = {k : v._bracket_(b2) for k,v in br1.items()}
            intdic = {k:v for k,v in intdic.items() if v}
            max1 = max(br1, default=0)
            max2 = max([max(v,default=0) for k,v in intdic.items()],default=0)
            ret3 = {n : sum(c2*binomial(n,k)*intdic[k][n-k-1] for k in intdic\
                    if n-k-1 in intdic[k]) for n in range(1,max1+max2+2)}
            ret = {}
            for d in [ret1,ret2,ret3]:
                for k in d:
                    ret[k] = ret.get(k,pz) + d[k]
            return {k:v for k,v in ret.items() if v}

        diclist = [i._bracket_(j) for i in self.terms() for j in other.terms()]
        ret = {}
        for d in diclist:
            for k in d:
                ret[k] = ret.get(k,pz) + d[k]
        return {k:v for k,v in ret.items() if v}

    def _mul_(self,right):
        r"""
        The normally ordered product of these two elements.

        INPUT:

        - ``right`` -- an element of this vertex algebra

        .. WARNING::

            The normally ordered product on a vertex algebra is not
            associative. If parenthesis are not added by default they
            are nested to the *left*, that is ``a*b*c == (a*b)*c``.

        EXAMPLES:

        The square of a Free Fermion vanishes::

            sage: V = vertex_algebras.FreeFermions(RR); V.inject_variables()
            Defining psi
            sage: psi*psi
            0

        This is not the case in the N=1 super vertex algebra::

            sage: V = vertex_algebras.NeveuSchwarz(QQ, 1); V.inject_variables()
            Defining L, G
            sage: G*G
            L_-3|0>

        The product is not associative on the Free Bosons::

            sage: V = vertex_algebras.FreeBosons(QQ); V.inject_variables()
            Defining alpha
            sage: (alpha*alpha)*alpha - alpha*(alpha*alpha)
            2*alpha_-3|0>

        The Sugawara construction for the affine Kac-Moody vertex
        algebra of `\mathfrak{sl}_2` at level ``1``::

            sage: V = vertex_algebras.Affine(QQ, 'A1', 1, names = ('e','h', 'f')); V
            The universal affine vertex algebra of CartanType ['A', 1] at level 1 over Rational Field
            sage: V.inject_variables()
            Defining e, h, f
            sage: L = 1/6*(e*f + h*h/2 + f*e); L.bracket(L)
            {0: 1/3*e_-2f_-1|0> + 1/6*h_-2h_-1|0> + 1/3*e_-1f_-2|0> - 1/3*h_-3|0>,
             1: 1/6*h_-1h_-1|0> + 2/3*e_-1f_-1|0> - 1/3*h_-2|0>,
             3: 1/2*|0>}
            sage: L.bracket(L) == {0: L.T(), 1:2*L, 3:V.vacuum()/2}
            True
        """
        p = self.parent()
        pz = p.zero()
        if self.is_zero() or right.is_zero():
            return pz
        if self.is_monomial() and right.is_monomial():
            a,b,c = self._pbw_one_less()
            if a == p.vacuum():
                return c*right
            if b == p.vacuum():
                a2,b2,c2 = right._pbw_one_less()
                if a2 == p.vacuum():
                    return c*c2*a
                i1 = next((j for j,x in enumerate(a.index()) if x))
                n1 = a.index()[i1][0]
                i2 = next((j for j,x in enumerate(a2.index()) if x))
                n2 = a2.index()[i2][0]
                if i1 < i2 or i1==i2 and n1 > n2 or a == a2 and a.is_even():
                    pt = right.index().to_list()
                    pt[i1].insert(0,n1)
                    return c*c2*p(pt)
                #Need to commute a with a2. If a = a2 is odd is special
                br = a._bracket_(a2)
                if a == a2:
                    ff = -sum([(-1)**(k+1)*v.T(k+1)/factorial(k+1) for k,v in\
                                br.items()],pz)/2
                    return c*c2*ff._mul_(b2)
                ff = sum([(-1)**(k+1)*v.T(k+1)/factorial(k+1) for k,v in\
                                br.items()],pz)
                ret = (-1)**(a.is_even_odd()*a2.is_even_odd())*\
                        a2._mul_(a._mul_(b2)) - p(ff)._mul_(b2)
                return c*c2*ret
            #quasi-associativity
            br1 = a._bracket_(right)
            br2 = b._bracket_(right)
            return c*a._mul_(b._mul_(right))+\
                    c*sum([a.T(k+1)._mul_(v)/factorial(k+1) for k,v in\
                    br2.items()],pz) + (-1)**(a.is_even_odd()*b.is_even_odd())*\
                    c*sum([b.T(k+1)._mul_(v)/factorial(k+1) for k,v in\
                    br1.items()],pz)
        return sum(i._mul_(j) for i in self.terms() for j in right.terms())

    def T(self,n=1):
        r"""
        The ``n``-th derivative of this element.

        INPUT:

        - ``n`` -- a non-negative integer (default: ``1``); the number
          of derivatives to apply.

        EXAMPLES::

            sage: V = vertex_algebras.Affine(QQ, 'A1', 1, names = ('e','h', 'f'));
            sage: V.0.T(3)
            6*e_-4|0>

            sage: V = vertex_algebras.Virasoro(QQ,1/2); V.inject_variables()
            Defining L
            sage: (L*L).T()
            2*L_-3L_-2|0> + L_-5|0>

            sage: V = vertex_algebras.NeveuSchwarz(QQ, 1); V.inject_variables()
            Defining L, G
            sage: (G.T()*G).T()
            2*G_-7/2G_-3/2|0> + L_-5|0>
            sage: G.T()
            G_-5/2|0>
        """
        if n==0:
            return self
        if n > 1:
            return self.T().T(n-1)

        p = self.parent()
        if self.is_monomial():
            if self.is_zero() or self.monomials()[0] == p.vacuum():
                return p.zero()

            #we cannot simply use Leibniz T(ab) = T(a)b + aT(b)
            #because the product uses T()
            (a,b,c) = self._pbw_one_less()
            if b == p.vacuum():
                i = next((j for j,x in enumerate(a.index()) if x))
                n = a.index()[i][0]
                return (c*p._lca.gen(i).T(n)/factorial(n-1)).lift()
            return c*a.T()._mul_(b) + c*a._mul_(b.T())

        return sum(m.T() for m in self.terms())

    def weight(self):
        """
        The conformal weight o this element.

        EXAMPLES::

            sage: V = vertex_algebras.Affine(QQ, 'A1', 1, names = ('e','h', 'f'));
            sage: V.inject_variables()
            Defining e, h, f
            sage: (f.T(3)*e).weight()
            5
            sage: V.vacuum().weight()
            0

        The zero element has infinite weight::

            sage: V(0).weight()
            +Infinity

        weight is only defined for homogeneous elements::

            sage: (e*f + e).weight()
            Traceback (most recent call last):
            ...
            ValueError: e_-1|0> + e_-1f_-1|0> is not homogeneous

        weight is only defined for graded vertex algebras::

            sage: V = vertex_algebras.Weyl(QQ); V.inject_variables()
            Defining alpha0, alpha1
            sage: alpha0.weight()
            Traceback (most recent call last):
            ...
            ValueError: Weight is only defined in H-graded vertex algebras
        """
        p = self.parent()
        if not p.is_graded():
            raise ValueError("Weight is only defined in H-graded vertex "\
                             "algebras")
        if self.is_zero():
            return Infinity

        if self.is_monomial():
            return self.index().energy()

        weightlist = [m.weight() for m in self.monomials()]
        if weightlist[1:] != weightlist[:-1]:
            raise ValueError("{} is not homogeneous".format(self))
        return weightlist[0]

    def _li_filtration_monomial_degree(self):
        """
        The Li filtration degree of this non-zero monomial.

        OUTPUT:

        The maximal `p` such that this element belongs
        to `F_pV`, where `F_pV` is the Li filtration of this vertex
        algebra.

        EXAMPLES::

            sage: V = vertex_algebras.Weyl(QQ); V.inject_variables()
            Defining alpha0, alpha1
            sage: v = alpha0*(alpha1*alpha1); v
            alpha0_(-1)alpha1_(-1)alpha1_(-1)|0>
            sage: v._li_filtration_monomial_degree()
            0
            sage: v = alpha0.T(2)*(alpha1.T()*alpha1); v
            2*alpha0_(-3)alpha1_(-2)alpha1_(-1)|0>
            sage: v._li_filtration_monomial_degree()
            3

        TESTS::

            sage: V.zero()._li_filtration_monomial_degree()
            +Infinity
            sage: V.vacuum()._li_filtration_monomial_degree()
            0
        """
        if self.is_zero():
            return Infinity
        if self.is_monomial():
            idx = self.index()
            return idx.size() - sum(len(j) for j in idx)
        raise ValueError("_li_filtration_monomial_degree is only defined"\
                         " for monomials, got {}".format(self))


    def pbw_filtration_degree(self):
        r"""
        The minimal `p` such that this element is contained in
        the `p`-th filtered piece of the vertex algebra with respect
        to the PBW filtration.

        EXAMPLES::

            sage: V = vertex_algebras.Affine(QQ, 'A1', 1, names=('e','h','f'))
            sage: V.inject_variables()
            Defining e, h, f
            sage: e.pbw_filtration_degree()
            1
            sage: (e*f.T()).pbw_filtration_degree()
            2
            sage: V([[4,2,1],[5,5],[]]).pbw_filtration_degree()
            5
            sage: (e*f - f*e).pbw_filtration_degree()
            1
            sage: V.vacuum().pbw_filtration_degree()
            0
            sage: V.zero().pbw_filtration_degree()
            -Infinity
        """
        if self.is_zero():
            return -Infinity
        return max(sum(j.length() for j in m.index()) for m in self.monomials())

    def _pbw_one_less(self):
        """
        Express this element as a product of two elements lower in
        the PBW filtration.

        OUTPUT:

        A triple ``(a,b,c)`` such that ``self=c*a*b`` with
        ``c`` a scalar and ``a`` a monomial of PBW degree ``1``.

        EXAMPLES::

            sage: V = vertex_algebras.Affine(QQ, 'A1', 1, names = ('e','h', 'f')); V.inject_variables(); V.register_lift()
            Defining e, h, f
            sage: v = e.T(2)*(e*(h.T()*f)); v
            2*e_-3e_-1h_-2f_-1|0>
            sage: v._pbw_one_less()
            (e_-3|0>, e_-1h_-2f_-1|0>, 2)
            sage: v = 3*V.vacuum(); v._pbw_one_less()
            (|0>, |0>, 3)
            sage: V = vertex_algebras.Virasoro(QQ);
            sage: v = V([[3,2,1,1]]); v
            L_-4L_-3L_-2L_-2|0>
            sage: v._pbw_one_less()
            (L_-4|0>, L_-3L_-2L_-2|0>, 1)
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
            return (p.vacuum(), p.vacuum(), c)
        l = [[],]*len(k)
        l[i] = [k[i].get_part(0)]
        sf = k.to_list()
        sf[i] = sf[i][1:]
        return (p(l),p(sf),c)


