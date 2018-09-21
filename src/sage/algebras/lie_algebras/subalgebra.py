r"""
Subalgebras and ideals of Lie algebras

AUTHORS:

- Eero Hakavuori (2018-08-29): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                 https://www.gnu.org/licenses/
# ****************************************************************************

from sage.algebras.lie_algebras.lie_algebra_element import LieSubalgebraElementWrapper
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.modules.free_module_element import vector
from sage.sets.family import Family
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class LieSubalgebra_finite_dimensional_with_basis(Parent, UniqueRepresentation):
    r"""
    A Lie subalgebra of a finite dimensional Lie algebra with basis.

    INPUT:

    - ``ambient`` -- the Lie algebra containing the subalgebra
    - ``gens`` -- a list of generators of the subalgebra
    - ``ideal`` -- (default: ``False``) a boolean; if ``True``, then ``gens``
      is interpreted as the generating set of an ideal instead of a subalgebra
    - ``category`` -- (optional) a subcategory of subobjects of finite
      dimensional Lie algebras with basis

    EXAMPLES:

    Subalgebras and ideals are defined by giving a list of generators::

        sage: L = lie_algebras.Heisenberg(QQ, 1)
        sage: X, Y, Z = L.basis()
        sage: S =  L.subalgebra([X, Z]); S
        Subalgebra generated by (p1, z) of Heisenberg algebra of rank 1 over Rational Field
        sage: I =  L.ideal([X, Z]); I
        Ideal (p1, z) of Heisenberg algebra of rank 1 over Rational Field

    An ideal is in general larger than the subalgebra with the same generators::

        sage: S = L.subalgebra(Y)
        sage: S.basis()
        Family (q1,)
        sage: I = L.ideal(Y)
        sage: I.basis()
        Family (q1, z)

    The zero dimensional subalgebra can be created by giving 0 as a generator
    or with an empty list of generators::

        sage: L.<X,Y,Z> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}})
        sage: S1 = L.subalgebra(0)
        sage: S2 = L.subalgebra([])
        sage: S1 is S2
        True
        sage: S1.basis()
        Family ()

    Elements of the ambient Lie algebra can be reduced modulo an
    ideal or subalgebra::

        sage: L.<X,Y,Z> = LieAlgebra(SR, {('X','Y'): {'Z': 1}})
        sage: I = L.ideal(Y)
        sage: I.reduce(X + 2*Y + 3*Z)
        X
        sage: S = L.subalgebra(Y)
        sage: S.reduce(X + 2*Y + 3*Z)
        X + 3*Z

    The reduction gives elements in a fixed complementary subspace.
    When the base ring is a field, the complementary subspace is spanned by
    those basis elements which are not leading supports of the basis::

        sage: I =  L.ideal(X + Y)
        sage: I.basis()
        Family (X + Y, Z)
        sage: el = var('x')*X + var('y')*Y + var('z')*Z; el
        x*X + y*Y + z*Z
        sage: I.reduce(el)
        (x-y)*X

    A subalgebra of a subalgebra is a subalgebra of the original::

        sage: sc = {('X','Y'): {'Z': 1}, ('X','Z'): {'W': 1}}
        sage: L.<X,Y,Z,W> = LieAlgebra(QQ, sc)
        sage: S1 = L.subalgebra([Y, Z, W]); S1
        Subalgebra generated by (Y, Z, W) of Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        sage: S2 = S1.subalgebra(S1.basis()[1:]); S2
        Subalgebra generated by (Z, W) of Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        sage: S3 = S2.subalgebra(S2.basis()[1:]); S3
        Subalgebra generated by (W) of Lie algebra on 4 generators (X, Y, Z, W) over Rational Field

    An ideal of an ideal is not necessarily an ideal of the original::

        sage: I = L.ideal(Y); I
        Ideal (Y) of Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        sage: J = I.ideal(Z); J
        Ideal (Z) of Ideal (Y) of Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        sage: J.basis()
        Family (Z,)
        sage: J.is_ideal(L)
        False
        sage: K = L.ideal(J.basis().list())
        sage: K.basis()
        Family (Z, W)

    TESTS:

    Test suites::

        sage: S =  L.subalgebra(X + Y)
        sage: TestSuite(S).run()
        sage: I =  L.ideal(X + Y)
        sage: TestSuite(I).run()

    Verify that subalgebras and ideals of nilpotent Lie algebras are nilpotent::

        sage: L = LieAlgebra(QQ, 3, step=4)
        sage: x,y,z = L.homogeneous_component_basis(1)
        sage: S = L.subalgebra([x, y])
        sage: S in LieAlgebras(QQ).Nilpotent()
        True
        sage: S.step()
        4
        sage: I = L.ideal(z)
        sage: I in LieAlgebras(QQ).Nilpotent()
        True
        sage: I.step()
        3

    Test computation for a nested ideal::

        sage: sc = {('X','Y'): {'Z': 1}, ('X','Z'): {'W': 1}}
        sage: L.<X,Y,Z,W> = LieAlgebra(QQ, sc)
        sage: I = L.ideal(Y)
        sage: J = I.ideal(Z)
        sage: J.reduce(I(Z) + I(W))
        W
    """

    @staticmethod
    def __classcall_private__(cls, ambient, gens, ideal=False, category=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES:

        Various ways to input one generator::

            sage: L.<X,Y> = LieAlgebra(QQ, {('X','Y'): {'X': 1}})
            sage: S1 = L.subalgebra(X)
            sage: S2 = L.subalgebra((X,))
            sage: S3 = L.subalgebra([X])
            sage: S1 is S2 and S2 is S3
            True

        Zero generators are ignored::

            sage: S1 = L.subalgebra(X)
            sage: S2 = L.subalgebra((X, 0))
            sage: S3 = L.subalgebra([X, 0, 0])
            sage: S1 is S2 and S2 is S3
            True
            sage: T1 = L.subalgebra(0)
            sage: T2 = L.subalgebra([])
            sage: T3 = L.subalgebra([0, 0])
            sage: T1 is T2 and T2 is T3
            True
        """
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        gens = tuple(ambient(gen) for gen in gens if not gen.is_zero())

        if not ideal and isinstance(ambient,
                            LieSubalgebra_finite_dimensional_with_basis):
            # a nested subalgebra is a subalgebra
            gens = tuple(ambient.lift(gen) for gen in gens)
            ambient = ambient.ambient()

        cat = LieAlgebras(ambient.base_ring()).FiniteDimensional().WithBasis()
        category = cat.Subobjects().or_subcategory(category)
        if ambient in LieAlgebras(ambient.base_ring()).Nilpotent():
            category = category.Nilpotent()

        sup = super(LieSubalgebra_finite_dimensional_with_basis, cls)
        return sup.__classcall__(cls, ambient, gens, ideal, category)

    def __init__(self, ambient, gens, ideal, category=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: L.<X,Y,Z> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}})
            sage: S = L.subalgebra(X)
            sage: TestSuite(S).run()
            sage: I = L.ideal(X)
            sage: TestSuite(I).run()
        """
        self._ambient = ambient
        self._gens = gens
        self._is_ideal = ideal

        sup = super(LieSubalgebra_finite_dimensional_with_basis, self)
        sup.__init__(ambient.base_ring(), category=category)

        # register a coercion to the ambient Lie algebra
        H = Hom(self, ambient)
        f = SetMorphism(H, self.lift)
        ambient.register_coercion(f)

    def __contains__(self, x):
        r"""
        Return ``True`` if ``x`` is an element of ``self``.

        EXAMPLES:

        Elements of the ambient Lie algebra are contained in the subalgebra
        if they are iterated brackets of the generators::

            sage: sc = {('x','y'): {'z': 1}, ('x','z'): {'w': 1}}
            sage: L.<x,y,z,w,u> = LieAlgebra(QQ, sc)
            sage: S = L.subalgebra([x, y])
            sage: z in S
            True
            sage: w in S
            True
            sage: u in S
            False

        TESTS::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z': 1}})
            sage: I = L.subalgebra(x)
            sage: I(x) in I
            True

        """
        if x in self.ambient():
            x = self.ambient()(x)
            return x.to_vector() in self.module()
        sup = super(LieSubalgebra_finite_dimensional_with_basis, self)
        return sup.__contains__(x)

    def __getitem__(self, x):
        r"""
        If `x` is a pair `(a, b)`, return the Lie bracket `[a, b]`.
        Otherwise try to return the `x`-th element of ``self``.

        This replicates the convenience syntax for Lie brackets of Lie algebras.

        EXAMPLES::

            sage: L.<x,y, z> = LieAlgebra(QQ, {('x','y'): {'z': 1}})
            sage: S = L.subalgebra([x, y])
            sage: a = S(x); b = S(y)
            sage: S[a, b]
            z
            sage: S[a, a + S[a,b]]
            0
        """
        if isinstance(x, tuple) and len(x) == 2:
            return self(x[0])._bracket_(self(x[1]))
        sup = super(LieSubalgebra_finite_dimensional_with_basis, self)
        return sup.__getitem__(x)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L.<X,Y> = LieAlgebra(QQ, abelian=True)
            sage: L.subalgebra([X, Y])
            Subalgebra generated by (X, Y) of Abelian Lie algebra on 2 generators (X, Y) over Rational Field
            sage: L.ideal([X, Y])
            Ideal (X, Y) of Abelian Lie algebra on 2 generators (X, Y) over Rational Field
        """
        gens = self.gens()
        if len(gens) == 1:
            gens = gens[0]

        if self._is_ideal:
            basestr = "Ideal"
        else:
            basestr = "Subalgebra generated by"

        return "%s %s of %s" % (basestr, self._repr_short(), self.ambient())

    def _repr_short(self):
        """
        Represent the list of generators.

        EXAMPLES::

            sage: L.<X,Y> = LieAlgebra(QQ, abelian=True)
            sage: L.ideal([X, Y])._repr_short()
            '(X, Y)'
            sage: L.ideal(X)._repr_short()
            '(X)'
            sage: L.subalgebra(X)._repr_short()
            '(X)'
        """
        return '(%s)' % (', '.join(str(X) for X in self.gens()))

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: L.<X,Y> = LieAlgebra(QQ, abelian=True)
            sage: S = L.subalgebra([X, Y])
            sage: S._an_element_()
            X
        """
        return self.element_class(self, self.lie_algebra_generators()[0])

    def _element_constructor_(self, x):
        """
        Convert ``x`` into ``self``.

        EXAMPLES:

        Elements of subalgebras are created directly from elements
        of the ambient Lie algebra::

            sage: L.<x,y,z,w> = LieAlgebra(ZZ, {('x','y'): {'z': 1}})
            sage: S = L.subalgebra([x, y])
            sage: S(y)
            y
            sage: S(y).parent()
            Subalgebra generated by (x, y) of Lie algebra on 4 generators (x, y, z, w) over Integer Ring

        A vector contained in the module corresponding to the subalgebra is
        interpreted as a coordinate vector::

            sage: S.module()
            Free module of degree 4 and rank 3 over Integer Ring
            User basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            sage: S(vector(ZZ, [2, 3, 5, 0]))
            2*x + 3*y + 5*z

        A list of 2 elements is interpreted as a Lie bracket::

            sage: S([S(x), S(y)])
            z
            sage: S([S(x), S(y)]) == S(L[x, y])
            True
        """
        try:
            P = x.parent()
            if P is self:
                return x
            if P == self.ambient():
                return self.retract(x)
        except AttributeError:
            pass

        if x in self.module():
            return self.from_vector(x)

        if isinstance(x, list) and len(x) == 2:
            return self(x[0])._bracket_(self(x[1]))

        sup = super(LieSubalgebra_finite_dimensional_with_basis, self)
        return sup._element_constructor_(x)

    # for submodule computations, the order of the basis is reversed so that
    # the pivot elements in the echelon form are the leading terms
    def _to_m(self, X):
        r"""
        Return the reversed vector of an element of the ambient Lie algebra.

        INPUT:

        - ``X`` -- an element of the ambient Lie algebra

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z': 1}})
            sage: I = L.ideal([x, z])
            sage: el = x + 2*y + 3*z
            sage: el.to_vector()
            (1, 2, 3)
            sage: I._to_m(el)
            (3, 2, 1)
        """
        return vector(self.ambient().base_ring(), reversed(X.to_vector()))

    def _from_m(self, v):
        r"""
        Return the element of the ambient Lie algebra from a reversed vector.

        INPUT:

        - ``v`` -- a vector

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z': 1}})
            sage: I = L.ideal([x, z])
            sage: I._from_m([3, 2, 1])
            x + 2*y + 3*z
        """
        R = self.ambient().base_ring()
        return self.ambient().from_vector(vector(R, reversed(v)))

    @lazy_attribute
    def _indices(self):
        r"""
        Return the set of indices for the basis of ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, abelian=True)
            sage: S = L.subalgebra([x, y])
            sage: S._indices
            {0, 1}
            sage: [S.basis()[k] for k in S._indices]
            [x, y]
        """
        return FiniteEnumeratedSet(self.basis().keys())

    @cached_method
    def zero(self):
        r"""
        Return the element `0`.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: S = L.subalgebra(x)
            sage: S.zero()
            0
            sage: S.zero() == S(L.zero())
            True
        """
        return self.element_class(self, self.ambient().zero())

    def ambient(self):
        r"""
        Return the ambient Lie algebra of ``self``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: S = L.subalgebra(x)
            sage: S.ambient() is L
            True
        """
        return self._ambient

    def lift(self, X):
        r"""
        Coerce an element ``X`` of ``self`` into the ambient Lie algebra.

        INPUT:

        - ``X`` -- an element of ``self``

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: S = L.subalgebra(x)
            sage: sx = S(x); sx
            x
            sage: sx.parent()
            Subalgebra generated by (x) of Abelian Lie algebra on 2 generators (x, y) over Rational Field
            sage: a = S.lift(sx); a
            x
            sage: a.parent()
            Abelian Lie algebra on 2 generators (x, y) over Rational Field
        """
        return X.value

    def retract(self, X):
        r"""
        Retract ``X`` to ``self``.

        INPUT:

        - ``X`` -- an element of the ambient Lie algebra

        EXAMPLES:

        Retraction to a subalgebra of a free nilpotent Lie algebra::

            sage: L = LieAlgebra(QQ, 3, step=2)
            sage: L.inject_variables()
            Defining X_1, X_2, X_3, X_12, X_13, X_23
            sage: S = L.subalgebra([X_1, X_2])
            sage: el = S.retract(2*X_1 + 3*X_2 + 5*X_12); el
            2*X_1 + 3*X_2 + 5*X_12
            sage: el.parent()
            Subalgebra generated by (X_1, X_2) of Free Nilpotent Lie algebra on
            6 generators (X_1, X_2, X_3, X_12, X_13, X_23) over Rational Field

        Retraction raises an error if the element is not contained in the
        subalgebra::

            sage: S.retract(X_3)
            Traceback (most recent call last):
            ...
            ValueError: the element X_3 is not in Subalgebra generated
            by (X_1, X_2) of Free Nilpotent Lie algebra on 6 generators
            (X_1, X_2, X_3, X_12, X_13, X_23) over Rational Field
        """
        if X not in self:
            raise ValueError("the element %s is not in %s" % (X, self))

        return self.element_class(self, X)

    def gens(self):
        r"""
        Return the generating set of ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z': 1}})
            sage: S = L.subalgebra(x)
            sage: S.gens()
            (x,)
        """
        return self._gens

    def lie_algebra_generators(self):
        r"""
        Return the generating set of ``self`` as a Lie algebra.

        EXAMPLES:

        The Lie algebra generators of a subalgebra are the original generators::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z': 1}})
            sage: S = L.subalgebra(x)
            sage: S.lie_algebra_generators()
            (x,)

        The Lie algebra generators of an ideal is usually a larger set::

            sage: I = L.ideal(x)
            sage: I.lie_algebra_generators()
            Family (x, z)
        """
        if self._is_ideal:
            return self.basis()
        return self._gens

    @cached_method
    def basis(self):
        r"""
        Return a basis of ``self``.

        EXAMPLES:

        A basis of a subalgebra::

            sage: sc = {('x','y'): {'z': 1}, ('x','z'): {'w': 1}}
            sage: L.<x,y,z,w> = LieAlgebra(QQ, sc)
            sage: L.subalgebra([x + y, z + w]).basis()
            Family (x + y, z, w)

        A basis of an ideal::

            sage: sc = {('x','y'): {'z': 1}, ('x','z'): {'w': 1}}
            sage: L.<x,y,z,w> = LieAlgebra(QQ, sc)
            sage: L.ideal([x + y + z + w]).basis()
            Family (x + y, z, w)
        """

        L = self.ambient()
        B = [self._to_m(X) for X in L.basis()]

        # use ambient module in case L is an ideal or subalgebra
        m = L.module().ambient_module()

        sm = m.submodule([self._to_m(X) for X in self.gens()])
        d = 0

        while sm.dimension() > d:
            d = sm.dimension()
            SB = sm.basis()
            if not self._is_ideal:
                B = SB

            brackets = [self._to_m(L.bracket(self._from_m(v), self._from_m(w)))
                        for v in B for w in SB]
            sm = m.submodule(sm.basis() + brackets)

        return Family(reversed([self.element_class(self, self._from_m(v))
                                for v in sm.echelonized_basis()]))

    def from_vector(self, v):
        r"""
        Return the element of ``self`` corresponding to the vector ``v``

        INPUT:

        - ``v`` -- a vector in ``self.module()`` or ``self.ambient().module()``

        EXAMPLES:

        An element from a vector of the intrinsic module::

            sage: L.<X,Y,Z> = LieAlgebra(ZZ, abelian=True)
            sage: L.dimension()
            3
            sage: S = L.subalgebra([X, Y])
            sage: S.dimension()
            2
            sage: el = S.from_vector([1, 2]); el
            X + 2*Y
            sage: el.parent() == S
            True

        An element from a vector of the ambient module

            sage: el = S.from_vector([1, 2, 0]); el
            X + 2*Y
            sage: el.parent() == S
            True
        """
        if len(v) == self.ambient().dimension():
            return self.retract(self.ambient().from_vector(v))

        sup = super(LieSubalgebra_finite_dimensional_with_basis, self)
        return sup.from_vector(v)

    def basis_matrix(self):
        r"""
        Return the basis matrix of ``self`` as a submodule
        of the ambient Lie algebra.

        EXAMPLES::

            sage: L.<X,Y,Z> = LieAlgebra(ZZ, {('X','Y'): {'Z': 3}})
            sage: S1 = L.subalgebra([4*X + Y, Y])
            sage: S1.basis_matrix()
            [ 4  0  0]
            [ 0  1  0]
            [ 0  0 12]
            sage: K.<X,Y,Z> = LieAlgebra(QQ, {('X','Y'): {'Z': 3}})
            sage: S2 = K.subalgebra([4*X + Y, Y])
            sage: S2.basis_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return self.module().basis_matrix()

    @cached_method
    def module(self, sparse=False):
        r"""
        Return the submodule of the ambient Lie algebra
        corresponding to ``self``.

        EXAMPLES::

            sage: L.<X,Y,Z> = LieAlgebra(ZZ, {('X','Y'): {'Z': 3}})
            sage: S = L.subalgebra([X, Y])
            sage: S.module()
            Free module of degree 3 and rank 3 over Integer Ring
            User basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 3]
        """
        try:
            m = self.ambient().module(sparse=sparse)
        except TypeError:
            m = self.ambient().module()
        ambientbasis = [self.lift(X).to_vector() for X in self.basis()]
        return m.submodule_with_basis(ambientbasis)

    @cached_method
    def is_ideal(self, A):
        """
        Return if ``self`` is an ideal of ``A``.

        EXAMPLES:

        Some subalgebras are ideals::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z': 1}})
            sage: S1 = L.subalgebra([x])
            sage: S1.is_ideal(L)
            False
            sage: S2 = L.subalgebra([x, y])
            sage: S2.is_ideal(L)
            True
            sage: S3 = L.subalgebra([y, z])
            sage: S3.is_ideal(L)
            True

        All ideals are ideals::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x': 1}})
            sage: I = L.ideal(x)
            sage: I.is_ideal(L)
            True
            sage: I.is_ideal(I)
            True

        TESTS::

            sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x': 1}})
            sage: I = L.ideal(x)
            sage: L.is_ideal(I)
            False
        """
        if A == self.ambient() and self._is_ideal:
            return True

        sup = super(LieSubalgebra_finite_dimensional_with_basis, self)
        return sup.is_ideal(A)

    def reduce(self, X):
        r"""
        Reduce an element of the ambient Lie algebra modulo the
        ideal ``self``.

        INPUT:

        - ``X`` -- an element of the ambient Lie algebra

        OUTPUT:

        An element `Y` of the ambient Lie algebra that is contained in a fixed
        complementary submodule `V` to ``self`` such that `X = Y` mod ``self``.

        When the base ring of ``self`` is a field, the complementary submodule
        `V` is spanned by the elements of the basis that are not the leading
        supports of the basis of ``self``.

        EXAMPLES:

        An example reduction in a 6 dimensional Lie algebra::

            sage: sc = {('a','b'): {'d': 1}, ('a','c'): {'e': 1},
            ....:       ('b','c'): {'f': 1}}
            sage: L.<a,b,c,d,e,f> = LieAlgebra(QQ, sc)
            sage: I =  L.ideal(c)
            sage: I.reduce(a + b + c + d + e + f)
            a + b + d

        The reduction of an element is zero if and only if the
        element belongs to the subalgebra::

            sage: I.reduce(c + e)
            0
            sage: c + e in I
            True

        Over non-fields, the complementary submodule may not be spanned by
        a subset of the basis of the ambient Lie algebra::

            sage: L.<X,Y,Z> = LieAlgebra(ZZ, {('X','Y'): {'Z': 3}})
            sage: I = L.ideal(Y)
            sage: I.basis()
            Family (Y, 3*Z)
            sage: I.reduce(3*Z)
            0
            sage: I.reduce(Y + 14*Z)
            2*Z
        """
        R = self.base_ring()
        for Y in self.basis():
            Y = self.lift(Y)
            k, c = Y.leading_item()

            if R.is_field():
                X = X - X[k] / c * Y
            else:
                try:
                    q, r = X[k].quo_rem(c)
                    X = X - q * Y
                except AttributeError:
                    pass

        return X

    class Element(LieSubalgebraElementWrapper):
        pass

