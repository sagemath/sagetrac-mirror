"""
Lie Conformal Algebra Element Class.

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
from sage.functions.other import factorial
from sage.misc.misc_c import prod
from sage.combinat.partition import Partition
from sage.misc.misc import repr_lincomb
from sage.misc.latex import latex

class LieConformalAlgebraElementWrapper(ElementWrapper):

    def _add_(self, right):
        """
        Add ``self`` and ``rhs``.
        """
        return type(self)(self.parent(), self.value + right.value)

    def _sub_(self, right):
        """
        Subtract ``self`` and ``rhs``.
        """
        return type(self)(self.parent(), self.value - right.value)

    def _mul_(left, right):
        return left.lift()*right.lift()

    def __neg__(self):
        """
        Return the negation of ``self``.
        """
        return type(self)(self.parent(), -self.value)

    def __getitem__(self, i):
        """
        Redirect the ``__getitem__()`` to the wrapped element.
        """
        return self.value.__getitem__(i)

    def _acted_upon_(self, scalar, self_on_left=False):
        """
        Return the product ``scalar`` times this element of the Lie
        conformal algebra.
        """
        return type(self)(self.parent(), scalar * self.value)

    def monomial_coefficients(self):
        """
        Return the monomial coefficients of this element as a
        dictionary.

        The keys are elements of the Lie conformal algebra.

        EXAMPLES::

            sage: V = VirasoroLieConformalAlgebra(QQ); V.inject_variables()
            Defining L, C
            sage: v = L + 2*L.T() + 3/2*C; Family(v.monomial_coefficients())
            Finite family {L: 1, TL: 2, C: 3/2}

            sage: R = WeylLieConformalAlgebra(AA)
            sage: R.inject_variables()
            Defining alpha0, alpha1, K
            sage: v = alpha0.T(3)+2*K-4*alpha1.T(2); Family(v.monomial_coefficients())
            Finite family {T^(3)alpha0: 6, K: 2, T^(2)alpha1: -8}

        TESTS::

            sage: R = FreeBosonsLieConformalAlgebra(ZZ,ngens=3)
            sage: R.zero().monomial_coefficients()
            {}
        """
        p = self.parent()
        return { p.monomial(k):v for k,v in
                self.value.monomial_coefficients().items() }

    def monomials(self):
        """
        The tuple monomials in this element.

        EXAMPLES::

            sage: R = NeveuSchwarzLieConformalAlgebra(QQ); R.inject_variables()
            Defining L, G, C
            sage: (L + G.T(2)).monomials()
            (L, 2*T^(2)G)

        TESTS::

            sage: R = FreeFermionsLieConformalAlgebra(ZZ,ngens=3); R.zero().monomials()
            ()
            sage: R.inject_variables()
            Defining psi_0, psi_1, psi_2, K
            sage: (psi_0.T(2) + R.zero()).monomials()
            (2*T^(2)psi_0,)
        """
        coefs = self.monomial_coefficients()
        return tuple(coefs[k]*k for k in coefs.keys())

class LCAWithGeneratorsElement(LieConformalAlgebraElementWrapper):
        """
        The element class of a Lie conformal algebra with a
        preferred set of generators.
        """
        def T(self,n=1):
            r"""
            The n-th derivative of this element.

            INPUT:

            - ``n`` -- integer (default:`1`); How many times
              to apply `T` to this element.

            We use the notation `T^{(j)} = \frac{T^j}{j!}`.

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ)
                sage: L = Vir.0; C = Vir.1
                sage: L.T()
                TL
                sage: L.T(3)
                6*T^(3)L
                sage: C.T()
                0

                sage: R = NeveuSchwarzLieConformalAlgebra(QQbar); R.inject_variables()
                Defining L, G, C
                sage: (L + 2*G.T() + 4*C).T(2)
                2*T^(2)L + 12*T^(3)G

            TESTS::

                sage: R = FreeBosonsLieConformalAlgebra(QQ); R.zero().T()
                0
                sage: (R.zero() + R.0).T()
                Talpha
            """
            if n ==  0 or self.is_zero():
                return self
            #it's faster to sum than to use recursion
            if self.is_monomial():
                p = self.parent()
                a,m = self.index()
                coef = self.value.monomial_coefficients()[(a,m)]
                if (a,m+n) in p._indices:
                    return coef*prod(j for j in range(m+1,m+n+1))\
                            *p.monomial((a,m+n))
                else:
                    return p.zero()

            return sum(mon.T(n) for mon in self.monomials())

        def lift(self):
            r"""
            Returns the image of this element under the canonical lift
            to the universal enveloping vertex algebra.

            .. WARNING::

                The universal enveloping algebra needs to be constructed
                first for this morphism to be defined.

                This morphism is registered as a coercion between this
                Lie conformal algebra and its universal enveloping
                vertex algebra upon creation. Since we consider central
                quotients of the universal enveloping vertex algebras
                by fixed central parameters, each time a different
                universal enveloping vertex algebra is constructed, this
                lift morphism is changed. See the examples below and
                also :meth:`register_lift(\
                )<sage.algebras.vertex_algebras.vertex_algebra.\
                UniversalEnvelopingVertexAlgebra.register_lift>`.

            EXAMPLES:

            We lift to the universal enveloping vertex algebra of the
            Virasoro Lie conformal algebra with central charge `0`::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: V = Vir.universal_enveloping_algebra()
                sage: L.lift()
                L_-2|0>
                sage: L.lift().__class__
                <class 'sage.algebras.vertex_algebras.universal_enveloping_vertex_algebra.UniversalEnvelopingVertexAlgebra_with_category.element_class'>
                sage: L.lift().parent()
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field

            Notice that the target of the ``lift`` morphism changes when
            we construct another universal enveloping vertex algebra::

                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field

                sage: V = VirasoroVertexAlgebra(QQ,1/2);
                sage: V.register_lift()
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra of central charge 1/2 over Rational Field

            Notice that recreation may not re-establish the right
            coercion depending on the method of construction::

                sage: Vir = VirasoroLieConformalAlgebra(QQ)
                sage: cp = Family({Vir.1:1/3}); V = Vir.universal_enveloping_algebra(cp)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field
                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field
                sage: V.register_lift()
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra of central charge 1/2 over Rational Field
            """
            p = self.parent()
            if not hasattr(p, 'lift'):
                raise NotImplementedError(
                    "In order to lift an element first need to "\
                    "construct the universal enveloping vertex "\
                    "algebra")
            V = p.lift.codomain()
            ret = V.zero()
            for c in self.value.monomial_coefficients().items():
                if p.monomial(c[0]) in p.central_elements():
                    ret += c[1]*V.central_parameters()[p.monomial(c[0])]*\
                                                                    V.vacuum()
                else:
                    l = [Partition([])]*V.ngens()
                    l[p._index_to_pos[c[0][0]]] = Partition([c[0][1]+1])
                    ret += c[1]*V(l)
            return ret

        def is_monomial(self):
            """
            Whether this element is a monomial.

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: (L + L.T()).is_monomial()
                False
                sage: L.T().is_monomial()
                True

            TESTS::

                sage: R = FreeFermionsLieConformalAlgebra(AA,ngens=3); R.inject_variables()
                Defining psi_0, psi_1, psi_2, K
                sage: R.zero().is_monomial()
                True
                sage: (psi_0.T(2) + R.zero()).is_monomial()
                True
            """
            return (len(self.value.monomial_coefficients()) == 1 \
                    or self.is_zero())

        def index(self):
            """
            The index parametrizing this monomial element.

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: L.index()
                ('L', 0)
                sage: L.T(4).index()
                ('L', 4)
                sage: (L + L.T()).index()
                Traceback (most recent call last):
                ...
                ValueError: index can only be computed for monomials
            """
            if self.is_zero():
                return tuple()
            if not self.is_monomial():
                raise ValueError ("index can only be computed for monomials")
            return list(self.value.monomial_coefficients().keys())[0]


class LCAStructureCoefficientsElement(LCAWithGeneratorsElement):
    """
    An element of a Lie conformal algebra given by structure
    coefficients.
    """

    def is_even_odd(self):
        """
        Return ``0`` if this element is even or ``1`` if it is odd.

        EXAMPLES::

            sage: R = lie_conformal_algebras.NeveuSchwarz(QQ); R.inject_variables()
            Defining L, G, C
            sage: L.is_even_odd()
            0
            sage: G.is_odd()
            True

        TESTS::

            sage: R = N2LieConformalAlgebra(QQ); R.inject_variables()
            Defining L, J, G1, G2, C
            sage: (G1 + R.zero()).is_odd()
            True
            sage: R.zero().is_even_odd()
            0
            sage: (J + R.zero()).is_odd()
            False
        """
        if self.is_zero():
            return 0
        p = self.parent()
        coefs = self.monomial_coefficients()
        paritylist = [p._parity[p.monomial((k.index()[0],0))] \
                      for k in coefs.keys()]
        if paritylist[1:] == paritylist[:-1]:
            return paritylist[0]
        raise ValueError("{} is not homogeneous".format(self))

    def _bracket_(self, right):
        """
        The lambda bracket of these two elements.

        The result is a dictionary with non-negative integer keys.
        The value corresponding to the entry `j` is ``self_{(j)}right``.

        EXAMPLES::

            sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
            sage: L.bracket(L)
            {0: TL, 1: 2*L, 3: 1/2*C}
            sage: L.T().bracket(L)
            {1: -TL, 2: -4*L, 4: -2*C}

            sage: R = lie_conformal_algebras.Affine(QQbar, 'A1', names=('e','h','f')); R
            The affine Lie conformal algebra of type ['A', 1] over Algebraic Field
            sage: R.inject_variables()
            Defining e, h, f, K
            sage: e.bracket(f)
            {0: h, 1: K}
            sage: h.bracket(h.T())
            {2: 4*K}

        When the Lie conformal algebra is actually a vertex algebra,
        elements are denoted in a different fashion::

            sage: V = VirasoroVertexAlgebra(QQ,1/2)
            sage: V.inject_variables()
            Defining L
            sage: Family(L.bracket(L*L))
            Finite family {0: 2*L_-3L_-2|0> + L_-5|0>,  1: 4*L_-2L_-2|0>,  2: 3*L_-3|0>,  3: 17/2*L_-2|0>,  5: 3/2*|0>}

        TESTS::

            sage: R = N2LieConformalAlgebra(QQbar); R.inject_variables()
            Defining L, J, G1, G2, C
            sage: G1.bracket(G2)
            {0: L + 1/2*TJ, 1: J, 2: 1/3*C}
            sage: G2.bracket(G1)
            {0: L - 1/2*TJ, 1: -J, 2: 1/3*C}
            sage: G1.bracket(C)
            {}
            sage: G1.bracket(R.zero())
            {}
            sage: G1.bracket(0)
            {}
            sage: R.zero().bracket(0)
            {}
        """
        p = self.parent()
        if self.is_monomial() and right.is_monomial():
            if self.is_zero() or right.is_zero():
                return {}
            s_coeff = p._s_coeff
            a,k = self.index()
            coefa = self.value.monomial_coefficients()[(a,k)]
            b,m = right.index()
            coefb = right.value.monomial_coefficients()[(b,m)]
            try:
                mbr = dict(s_coeff[(a,b)])
            except KeyError:
                return {}
            pole = max(mbr.keys())
            ret =  {l: coefa*coefb*(-1)**k/factorial(k)*sum(factorial(l)\
                    /factorial(m+k+j-l)/factorial(l-k-j)/factorial(j)*\
                    mbr[j].T(m+k+j-l) for j in mbr.keys() if j >= l-m-k and\
                    j <= l-k) for l in range(m+k+pole+1)}
            return {k:v for k,v in ret.items() if v}

        diclist = [ i._bracket_(j) for i in self.monomials() for
                j in right.monomials() ]
        ret = {}
        pz = p.zero()
        for d in diclist:
            for k in d.keys():
                ret[k] = ret.get(k,pz) + d[k]
        return ret

    def __getitem__(self, i):
        """
        Return the coefficient of the basis element indexed by ``i``.

        ``i`` must be a basis element of this Lie conformal algebra.

        EXAMPLES::

            sage: V = VirasoroLieConformalAlgebra(QQ); V.inject_variables()
            Defining L, C
            sage: v = L + L.T(2) + 3/2*C; v.monomials()
            (L, 2*T^(2)L, 3/2*C)
            sage: v[L]
            1
            sage: v[L.T(2)/2]
            2
            sage: v[L.T(2)]
            Traceback (most recent call last):
            ...
            KeyError: 2*T^(2)L
       """
        return self.monomial_coefficients()[i]

    def _repr_(self):
        r"""
        A visual representation of this element.

        For a free generator `L`, the element `\frac{T^{j}}{j!}L` is
        denoted by ``T^(j)L``.

        EXAMPLES::

            sage: V = lie_conformal_algebras.Virasoro(QQ); V.inject_variables()
            Defining L, C
            sage: v = L.T(5).nproduct(L,6); v
            -1440*L
            sage: L.T(2) + L + C
            2*T^(2)L + L + C
            sage: L.T(4)
            24*T^(4)L
        """
        if self.is_zero():
            return "0";
        p = self.parent()
        if p._names:
            terms = [("T^({0}){1}".format(k[1],
                        p._names[p._index_to_pos[k[0]]]),v) if k[1] > 1 \
                    else("T{}".format(p._names[p._index_to_pos[k[0]]]),v) \
                    if k[1] == 1 \
                    else ("{}".format(p._names[p._index_to_pos[k[0]]]),v)\
                        for k,v in self.value.monomial_coefficients().items()]
        else:
            terms = [("T^({0}){1}".format(k[1], k[0]),v) if k[1] > 1 \
                      else("T{}".format(k[0]),v) if k[1] == 1 \
                        else ("{}".format(k[0]),v)\
                        for k,v in self.value.monomial_coefficients().items()]

        return repr_lincomb(terms, strip_one = True)

    def _latex_(self):
        r"""
        A visual representation of this element.

        For a free generator `L`, the element `\frac{T^{j}}{j!}L` is
        denoted by ``T^(j)L``.

        EXAMPLES::

            sage: V = lie_conformal_algebras.Virasoro(QQ); V.inject_variables()
            Defining L, C
            sage: latex(L.T(2))
            2T^{(2)}L

            sage: R = lie_conformal_algebras.Affine(QQbar, 'A1', names=('e','h','f')); R.inject_variables()
            Defining e, h, f, K
            sage: latex(e.bracket(f))
            \left\{0 : h, 1 : K\right\}
            sage: latex(e.T(3))
            6T^{(3)}e

            sage: R = lie_conformal_algebras.Affine(QQbar, 'A1')
            sage: latex(R.0.bracket(R.2))
            \left\{0 : \alpha^\vee_{1}, 1 : \text{\texttt{K}}\right\}

            sage: R = WeylLieConformalAlgebra(QQ,ngens=4); R.inject_variables()
            Defining alpha0, alpha1, alpha2, alpha3, K
            sage: latex(alpha1.T(4))
            24T^{(4)}\alpha_{1}

            sage: R = AffineLieConformalAlgebra(QQ, 'A1'); latex(R.0.T(3))
            6T^{(3)}\alpha_{1}

            sage: R = BosonicGhostsLieConformalAlgebra(QQ,ngens=4); R.inject_variables()
            Defining beta0, beta1, gamma0, gamma1, K
            sage: latex(beta1.T(3))
            6T^{(3)}\beta_{1}
        """
        if self.is_zero():
            return "0";
        p = self.parent()
        try:
            names = p.latex_variable_names()
        except ValueError:
            names = None
        if names:
            terms = [("T^{{({0})}}{1}".format(k[1],
                        names[p._index_to_pos[k[0]]]),v) if k[1] > 1 \
                else("T{}".format(names[p._index_to_pos[k[0]]]),v)\
                if k[1] == 1\
                else ("{}".format(names[p._index_to_pos[k[0]]]),v)\
                        for k,v in self.value.monomial_coefficients().items()]
        else:
            terms = [("T^{{({0})}}{1}".format(k[1], latex(k[0])),v) if k[1] > 1 \
                      else("T{}".format(latex(k[0])),v) if k[1] == 1 \
                        else ("{}".format(latex(k[0])),v)\
                        for k,v in self.value.monomial_coefficients().items()]

        return repr_lincomb(terms, is_latex=True, strip_one = True)


