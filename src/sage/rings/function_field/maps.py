r"""
Function Field Morphisms

AUTHORS:

- William Stein (2010): initial version

- Julian Rueth (2011-09-14, 2014-06-23): refactored class hierarchy; added
  derivation classes

EXAMPLES::

    sage: K.<x> = FunctionField(QQ); R.<y> = K[]
    sage: K.hom(1/x)
    Function Field endomorphism of Rational function field in x over Rational Field
      Defn: x |--> 1/x
    sage: L.<y> = K.extension(y^2 - x)
    sage: K.hom(y)
    Function Field morphism:
      From: Rational function field in x over Rational Field
      To:   Function field in y defined by y^2 - x
      Defn: x |--> y
    sage: L.hom([y,x])
    Function Field endomorphism of Function field in y defined by y^2 - x
      Defn: y |--> y
            x |--> x
    sage: L.hom([x,y])
    Traceback (most recent call last):
    ...
    ValueError: invalid morphism
"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2011-2014 Julian Rueth <julian.rueth@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.morphism import Morphism
from sage.categories.map import Map
from sage.rings.morphism import RingHomomorphism

class FunctionFieldDerivation(Map):
    r"""
    A base class for derivations on function fields.

    A derivation on `R` is map `R\to R` with
    `D(\alpha+\beta)=D(\alpha)+D(\beta)` and `D(\alpha\beta)=\beta
    D(\alpha)+\alpha D(\beta)` for all `\alpha,\beta\in R`.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: d = K.derivation()
        sage: isinstance(d, sage.rings.function_field.maps.FunctionFieldDerivation)
        True

    """
    def __init__(self, K):
        r"""
        Initialize a derivation from ``K`` to ``K``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation() # indirect doctest

        """
        from .function_field import is_FunctionField
        if not is_FunctionField(K):
            raise ValueError("K must be a function field")
        self.__field = K
        from sage.categories.homset import Hom
        from sage.categories.sets_cat import Sets
        Map.__init__(self, Hom(K,K,Sets()))

    def _repr_type(self):
        r"""
        Return the type of this map (a derivation), for the purposes of printing out self.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d._repr_type()
            'Derivation'

        """
        return "Derivation"

    def is_injective(self):
        r"""
        Return whether this derivation is injective.

        OUTPUT:

        Returns ``False`` since derivations are never injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d.is_injective()
            False

        """
        return False

class FunctionFieldDerivation_rational(FunctionFieldDerivation):
    """
    A derivation on a rational function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: d = K.derivation()
        sage: isinstance(d, sage.rings.function_field.maps.FunctionFieldDerivation_rational)
        True
    """
    def __init__(self, K, u):
        """
        Initialize a derivation of ``K`` which sends the generator of ``K`` to
        ``u``.

        INPUT:

        - ``K`` -- a rational function field

        - ``u`` -- an element of ``K``, the image of the generator of ``K`` under
          the derivation

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation() # indirect doctest
            sage: type(d)
            <class 'sage.rings.function_field.maps.FunctionFieldDerivation_rational'>
        """
        FunctionFieldDerivation.__init__(self, K)

        self._u = u

    def _call_(self, x):
        """
        Compute the derivation of ``x``.

        INPUT:

        - ``x`` -- an element of the rational function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d(x) # indirect doctest
            1
            sage: d(x^3)
            3*x^2
            sage: d(1/x)
            -1/x^2
        """
        f = x.numerator()
        g = x.denominator()

        numerator = f.derivative() * g - f * g.derivative()
        if numerator.is_zero():
            return self.codomain().zero()
        else:
            return self._u * self.codomain()(numerator / g**2)

class FunctionFieldDerivation_separable(FunctionFieldDerivation):
    """
    The unique extension of the derivation ``d`` to ``L``.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: d = L.derivation()
    """
    def __init__(self, L, d):
        """
        Initialization.

        INPUT:

        - ``L`` -- a function field which is a separable extension of the domain of
          ``d``

        - ``d`` -- a derivation on the base function field of ``L``

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation() # indirect doctest
            sage: type(d)
            <class 'sage.rings.function_field.maps.FunctionFieldDerivation_separable'>
        """
        FunctionFieldDerivation.__init__(self, L)

        f = self.domain().polynomial()
        if not f.gcd(f.derivative()).is_one():
            raise ValueError("L must be a separable extension of its base field.")

        x = self.domain().gen()

        self._d = d
        self._gen_image = - f.map_coefficients(lambda c:d(c))(x) / f.derivative()(x)

    def _call_(self, x):
        r"""
        Evaluate the derivation on ``x``.

        INPUT:

        - ``x`` -- an element of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: d(x) # indirect doctest
            1
            sage: d(y)
            (-1/2/-x)*y
            sage: d(y^2)
            1
        """
        if x.is_zero():
            return self.codomain().zero()

        x = x._x
        y = self.domain().gen()

        return x.map_coefficients(self._d) + x.derivative()(y) * self._gen_image

    def _repr_defn(self):
        """
        Return the string representation of the map.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L.derivation() # indirect doctest
            Derivation map:
              From: Function field in y defined by y^2 - x
              To:   Function field in y defined by y^2 - x
              Defn: y |--> (-1/2/-x)*y

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.derivation()
            Derivation map:
              From: Function field in z defined by z^2 - y
              To:   Function field in z defined by z^2 - y
              Defn: y |--> (-1/2/-x)*y
                    z |--> 1/4/x*z
        """
        base = self._d._repr_defn()
        ret = '{} |--> {}'.format(self.domain().gen(), self._gen_image)
        if base:
            return base + '\n' + ret
        else:
            return ret

class FunctionFieldIsomorphism(Morphism):
    r"""
    A base class for isomorphisms between function fields and
    vector spaces.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space()
        sage: isinstance(f, sage.rings.function_field.maps.FunctionFieldIsomorphism)
        True
    """
    def _repr_type(self):
        """
        Return the type of this map (an isomorphism), for the purposes of
        printing this map.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"

    def is_injective(self):
        """
        Return True, since this isomorphism is injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        Return True, since this isomorphism is surjective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.is_surjective()
            True
        """
        return True

class MapVectorSpaceToFunctionField(FunctionFieldIsomorphism):
    r"""
    An isomorphism from a vector space to a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space(); f
        Isomorphism morphism:
          From: Vector space of dimension 2 over Rational function field in x over Rational Field
          To:   Function field in y defined by y^2 - x*y + 4*x^3
    """
    def __init__(self, V, K):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space(); type(f)
            <class 'sage.rings.function_field.maps.MapVectorSpaceToFunctionField'>
        """
        self._V = V
        self._K = K
        self._R = K.polynomial_ring()
        from sage.categories.homset import Hom
        FunctionFieldIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        """
        Map ``v`` to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f(x*V.0 + (1/x^3)*V.1) # indirect doctest
            1/x^3*y + x

        TESTS:

        Test that this map is a bijection for some random inputs::

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y - x)
            sage: for F in [K,L,M]:
            ....:     for base in F._intermediate_fields(K):
            ....:         V, f, t = F.vector_space(base)
            ....:         for i in range(100):
            ....:             a = F.random_element()
            ....:             assert(f(t(a)) == a)

        """
        fields = self._K._intermediate_fields(self._V.base_field())
        fields.pop()
        degrees = [k.degree() for k in fields]
        gens = [k.gen() for k in fields]

        # construct the basis composed of powers of the generators of all the
        # intermediate fields, i.e., 1, x, y, x*y, ...
        from sage.misc.misc_c import prod
        from itertools import product
        exponents = product(*[range(d) for d in degrees])
        basis = [prod(g**e for g,e in zip(gens,es)) for es in exponents]

        # multiply the entries of v with the values in basis
        coefficients = self._V(v).list()
        ret = sum([c*b for (c,b) in zip(coefficients,basis)])
        return self._K(ret)

    def domain(self):
        """
        Return the vector space which is the domain of this isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.domain()
            Vector space of dimension 2 over Rational function field in x over Rational Field
        """
        return self._V

    def codomain(self):
        """
        Return the function field which is the codomain of this isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.codomain()
            Function field in y defined by y^2 - x*y + 4*x^3
        """
        return self._K

class MapFunctionFieldToVectorSpace(FunctionFieldIsomorphism):
    """
    An isomorphism from a function field to a vector space.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space(); t
        Isomorphism morphism:
          From: Function field in y defined by y^2 - x*y + 4*x^3
          To:   Vector space of dimension 2 over Rational function field in x over Rational Field
    """
    def __init__(self, K, V):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space(); type(t)
            <class 'sage.rings.function_field.maps.MapFunctionFieldToVectorSpace'>
        """
        self._V = V
        self._K = K
        self._zero = K.base_ring()(0)
        self._n = K.degree()
        from sage.categories.homset import Hom
        FunctionFieldIsomorphism.__init__(self, Hom(K, V))

    def domain(self):
        """
        Return the function field which is the domain of this isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t.domain()
            Function field in y defined by y^2 - x*y + 4*x^3
        """
        return self._K

    def codomain(self):
        """
        Return the vector space which is the domain of this isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t.codomain()
            Vector space of dimension 2 over Rational function field in x over Rational Field
        """
        return self._V

    def _repr_type(self):
        """
        Return the type of this map (an isomorphism), for the purposes of
        printing this map.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"

    def _call_(self, x):
        """
        Map ``x`` to the vector space.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t(x + (1/x^3)*y) # indirect doctest
            (x, 1/x^3)

        TESTS:

        Test that this map is a bijection for some random inputs::

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y - x)
            sage: for F in [K,L,M]:
            ....:     for base in F._intermediate_fields(K):
            ....:         V, f, t = F.vector_space(base)
            ....:         for i in range(100):
            ....:             a = V.random_element()
            ....:             assert(t(f(a)) == a)

        """
        ret = [x]
        fields = self._K._intermediate_fields(self._V.base_field())
        fields.pop()
        from itertools import chain
        for k in fields:
            ret = chain.from_iterable([y.list() for y in ret])
        ret = list(ret)
        assert all([t.parent() is self._V.base_field() for t in ret])
        return self._V(ret)

class FunctionFieldMorphism(RingHomomorphism):
    r"""
    Base class for morphisms between function fields.
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: f = K.hom(1/x); f
            Function Field endomorphism of Rational function field in x over Rational Field
              Defn: x |--> 1/x
            sage: isinstance(f, sage.rings.function_field.maps.FunctionFieldMorphism)
            True
        """
        RingHomomorphism.__init__(self, parent)

        self._im_gen = im_gen
        self._base_morphism = base_morphism

    def _repr_type(self):
        r"""
        Return the type of this map (a morphism of function fields), for the
        purposes of printing this map.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2)
            sage: f._repr_type()
            'Function Field'

        """
        return "Function Field"

    def _repr_defn(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2)
            sage: f._repr_defn()
            'y |--> 2*y'
        """
        a = '%s |--> %s'%(self.domain().gen(), self._im_gen)
        if self._base_morphism is not None:
            a += '\n' + self._base_morphism._repr_defn()
        return a

class FunctionFieldMorphism_polymod(FunctionFieldMorphism):
    """
    Morphism from a finite extension of a function field to a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: f = L.hom(-y); f
        Function Field endomorphism of Function field in y defined by y^2 - x
          Defn: y |--> -y
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2); f
            Function Field endomorphism of Function field in y defined by y^3 + 6*x^3 + x
              Defn: y |--> 2*y
            sage: type(f)
            <class 'sage.rings.function_field.maps.FunctionFieldMorphism_polymod'>
            sage: factor(L.polynomial())
            y^3 + 6*x^3 + x
            sage: f(y).charpoly('y')
            y^3 + 6*x^3 + x
        """
        FunctionFieldMorphism.__init__(self, parent, im_gen, base_morphism)
        # Verify that the morphism is valid:
        R = self.codomain()['X']
        v = parent.domain().polynomial().list()
        if base_morphism is not None:
            v = [base_morphism(a) for a in v]
        f = R(v)
        if f(im_gen):
            raise ValueError("invalid morphism")

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x); f = L.hom(y*2)
            sage: f(y/x + x^2/(x+1))            # indirect doctest
            2/x*y + x^2/(x + 1)
            sage: f(y)
            2*y
        """
        v = x.list()
        if self._base_morphism is not None:
            v = [self._base_morphism(a) for a in v]
        f = v[0].parent()['X'](v)
        return f(self._im_gen)

class FunctionFieldMorphism_rational(FunctionFieldMorphism):
    """
    Morphism from a rational function field to a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: f = K.hom(1/x); f
        Function Field endomorphism of Rational function field in x over Rational Field
          Defn: x |--> 1/x
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: f = K.hom(1/x); f
            Function Field endomorphism of Rational function field in x over Finite Field of size 7
              Defn: x |--> 1/x
            sage: type(f)
            <class 'sage.rings.function_field.maps.FunctionFieldMorphism_rational'>
        """
        FunctionFieldMorphism.__init__(self, parent, im_gen, base_morphism)

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: f = K.hom(1/x); f
            Function Field endomorphism of Rational function field in x over Finite Field of size 7
              Defn: x |--> 1/x
            sage: f(x+1)                          # indirect doctest
            (x + 1)/x
            sage: 1/x + 1
            (x + 1)/x

        You can specify a morphism on the base ring::

            sage: Qi = GaussianIntegers().fraction_field()
            sage: i = Qi.gen()
            sage: K.<x> = FunctionField(Qi)
            sage: phi1 = Qi.hom([CC.gen()])
            sage: phi2 = Qi.hom([-CC.gen()])
            sage: f = K.hom(CC.pi(),phi1)
            sage: f(1+i+x)
            4.14159265358979 + 1.00000000000000*I
            sage: g = K.hom(CC.pi(),phi2)
            sage: g(1+i+x)
            4.14159265358979 - 1.00000000000000*I
        """
        a = x.element()
        if self._base_morphism is None:
            return a.subs({a.parent().gen():self._im_gen})
        else:
            f = self._base_morphism
            num = a.numerator()
            den = a.denominator()
            R = self._im_gen.parent()['X']
            num = R([f(c) for c in num.list()])
            den = R([f(c) for c in den.list()])
            return num.subs(self._im_gen) / den.subs(self._im_gen)

class FunctionFieldConversionToConstantBaseField(Map):
    r"""
    Conversion map from the function field to its constant base field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: QQ.convert_map_from(K)
        Conversion map:
          From: Rational function field in x over Rational Field
          To:   Rational Field

    """
    def __init__(self, parent):
        r"""
        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: f = QQ.convert_map_from(K)
            sage: from sage.rings.function_field.maps import FunctionFieldConversionToConstantBaseField
            sage: isinstance(f, FunctionFieldConversionToConstantBaseField)
            True

        """
        Map.__init__(self, parent)

    def _repr_type(self):
        r"""
        Return the type of this map (a conversion), for the purposes of printing out self.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: QQ.convert_map_from(K) # indirect doctest
            Conversion map:
              From: Rational function field in x over Rational Field
              To:   Rational Field

        """
        return "Conversion"

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: QQ(K(1)) # indirect doctest
            1

        """
        return x.parent()._to_constant_base_field(x)
