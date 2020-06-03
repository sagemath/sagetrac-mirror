r"""
Lie conformal algebra element class
AUTHORS

- Reimundo Heluani (08-09-2019): Initial implementation
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
from sage.structure.element import parent
from sage.rings.integer import Integer
from sage.combinat.partition import Partition

class LieConformalAlgebraElementWrapper(ElementWrapper):

    def _repr_(self):
        """
        Return a string representation of ``self``.
       """
        return repr(self.value)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.
       """
        from sage.misc.latex import latex
        return latex(self.value)

    def _ascii_art_(self):
        """
        Return an ascii art representation of ``self``.
       """
        from sage.typeset.ascii_art import ascii_art
        return ascii_art(self.value)

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.
       """
        from sage.typeset.unicode_art import unicode_art
        return unicode_art(self.value)

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

    # Need to bypass the coercion model
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
        Return the product ``scalar`` times this element of the Lie conformal
        algebra
        """
        return type(self)(self.parent(), scalar * self.value)

    def monomial_coefficients(self):
        """
        Return the monomial coefficients of this element as a dictionary.

        The keys are element of the Lie conformal algebra

        EXAMPLES::

            sage: V = VirasoroLieConformalAlgebra(QQ); V.inject_variables()
            Defining L, C
            sage: v = L + 2*L.T() + 3/2*C; v.monomial_coefficients()
            {TL: 2, L: 1, C: 3/2}

        """
        p = self.parent()
        return { p.monomial(k):v for k,v in 
                self.value.monomial_coefficients().items() }

    def monomials(self):
        """
        The monomials in this element

        EXAMPLES::

            sage: R = NeveuSchwarzLieConformalAlgebra(QQ); R.inject_variables()
            Defining L, G, C
            sage: (L + G.T(2)).monomials()
            (L, 2*T^(2)G)
        """
        coefs = self.monomial_coefficients()
        return tuple(coefs[k]*k for k in coefs.keys())

class LCAWithGeneratorsElement(LieConformalAlgebraElementWrapper):
        """
        The element class of a Lie conformal algebra with a 
        preferred set of generators
        """
        def T(self,n=1):
            r"""
            The n-th derivative of this element

            INPUT:

            - ``n`` -- integer (default:`1`); How many times
            to apply `T` to this element. 

            We use the notation `T^{(j)} = \frac{T^j}{j!}` 

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ)
                sage: L = Vir.0; C = Vir.1
                sage: L.T()
                TL
                sage: L.T(3)
                6*T^(3)L
                sage: C.T()
                0
            """
            if n == 0:
                return self
            coef = self.value.monomial_coefficients()
            p = self.parent()
            ret = p.zero()
            for k in coef.keys():
                if (k[0],k[1]+1) in p._indices:
                    ret += prod(range(k[1]+1,k[1]+n+1))*coef[k]\
                                *p.monomial((k[0],k[1]+n))
            return ret

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
                <class 'sage.algebras.vertex_algebras.vertex_algebra.UniversalEnvelopingVertexAlgebra_with_category.element_class'>
                sage: L.lift().parent()
                The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.

            Notice that the target of the ``lift`` morphism changes when
            we construct another universal enveloping vertex algebra::

                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.
                sage: V = VirasoroVertexAlgebra(QQ,1/2);
                sage: V.register_lift()
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra at central charge 1/2

            Notice that recreation may not re-establish the right 
            coercion depending on the method of construction::

                sage: Vir = VirasoroLieConformalAlgebra(QQ)
                sage: cp = Family({Vir.1:1/3}); V = Vir.universal_enveloping_algebra(cp)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.
                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.
                sage: V.register_lift()
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra at central charge 1/2                

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
                    ret += c[1]*V.central_parameters()[
                                p.monomial(c[0])]*V.vacuum()
                else:
                    l = [Partition([])]*V.ngens()
                    l[p._index_to_pos[c[0][0]]] = Partition(
                                                    [c[0][1]+1])
                    ret += c[1]*V(l)
            return ret

        def is_monomial(self):
            r"""
            Whether this element is a monomial

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: (L + L.T()).is_monomial()
                False
                sage: L.T().is_monomial()
                True

            """
            return (len(self.value.monomial_coefficients()) == 1  or self.is_zero())

        def index(self):
            r"""
            The index parametrizing this monomial element

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
    An element of a Lie conformal algebra given by structure coefficients.
    """

    def _bracket_(self, right):
        """
        Return the lambda bracket of these two elements.

        The result is a dictionary with non-negative integer keys. The value
        corresponding to the entry `j` is ``self_{(j)}right``

        EXAMPLES::

            sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
            sage: L.bracket(L)
            {0: TL, 1: 2*L, 3: 1/2*C}
            sage: L.T().bracket(L)
            {1: -TL, 2: -4*L, 4: -2*C}

        When the Lie conformal algebra is actually a vertex algebra, elements
        are denoted in a different fashion::

            sage: V = VirasoroVertexAlgebra(QQ,1/2)
            sage: V.inject_variables()
            Defining L
            sage: L.bracket(L*L)
            {0: L_-5|0>+2*L_-3L_-2|0>,
             1: 4*L_-2L_-2|0>,
             2: 3*L_-3|0>,
             3: 17/2*L_-2|0>,
             5: 3/2*|0>}

        """
        p = self.parent()
        if self.is_monomial() and right.is_monomial():
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

        ``i`` must be a basis element of this Lie conformal algebra

        EXAMPLES::

            sage: V = VirasoroLieConformalAlgebra(QQ); V.inject_variables()
            Defining L, C
            sage: v = L + 2*L.T() + 3/2*C; v.monomial_coefficients()
            {TL: 2, L: 1, C: 3/2}
            sage: v[L]
            1
            sage: v[C]
            3/2

       """
        return self.monomial_coefficients()[i]

    def _repr_(self):
        r"""
        A visual representation of this element

        For a free generator `L`, the element `\frac{T^{j}}{j!}L` is denoted by
        ``T^(j)L``.

        EXAMPLES::

            sage: V = VirasoroLieConformalAlgebra(QQ); V.inject_variables()
            Defining L, C
            sage: v = L.T(5).nproduct(L,6); v
            -1440*L
            sage: L.T(2) + L + C
            C+2*T^(2)L+L
            sage: L.T(4)
            24*T^(4)L

        """

        if self == self.parent().zero():
            return "0";
        coeff = list(self.value.monomial_coefficients().items())
        ret = ""
        p = self.parent()
        try:
            names = p.variable_names()
        except ValueError:
            names = None

        for i in range(len(coeff)):    
            #TODO: deal with this without using an ordering of the ring. 
            if i > 0  and coeff[i][1] > 0 :
                ret += "+"
            if coeff[i][1] < 0 :
                ret += "-"
    
            if abs(coeff[i][1]) != 1 :
                ret += "{}*".format(abs(coeff[i][1]))
            if coeff[i][0][1] == 1 :
                ret += "T"
            elif coeff[i][0][1] > 1 :
                ret += "T^({})".format(coeff[i][0][1])
            if names:
                idx = p._index_to_pos[coeff[i][0][0]]
                ret += "{}".format(names[idx])
            else:
                ret += "{}".format(coeff[i][0][0])
        return ret
    
    def _latex_(self):
        r"""
        A visual representation of this element

        For a free generator `L`, the element `\frac{T^{j}}{j!}L` is denoted by
        ``T^(j)L``.

        EXAMPLES::

            sage: V = VirasoroLieConformalAlgebra(QQ); V.inject_variables()
            Defining L, C
            sage: v = L.T(5).nproduct(L,6); latex(v)
            -1440*L
            sage: latex(L.T(3))
            6*T^{(3)}L
            sage: latex(L.T(2)  + C)
            C+2*T^{(2)}L

        """
        if self == self.parent().zero():
            return "0";
        coeff = self.value.monomial_coefficients().items()
        ret = ""
        p = self.parent()
        try:
            names = p.variable_names()
        except ValueError:
            names = None

        for i in range(len (coeff)):    
            #TODO: deal with this without using an ordering of the ring. 
            if i > 0  and coeff[i][1] > 0 :
                ret += "+"
            if coeff[i][1] < 0 :
                ret += "-"
    
            if abs(coeff[i][1]) != 1 :
                ret += "{}*".format(abs(coeff[i][1]))
            if coeff[i][0][1] == 1 :
                ret += "T"
            elif coeff[i][0][1] > 1 :
                ret += "T^{("+"{}".format(coeff[i][0][1])+")}"
            if names:
                idx = p._index_to_pos[coeff[i][0][0]]
                ret += "{}".format(names[idx])
            else:
                ret += "{}".format(coeff[i][0][0])
        return ret

 

