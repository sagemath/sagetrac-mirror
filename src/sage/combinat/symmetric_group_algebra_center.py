r"""
The center of the group algebra of the symmetric group

For any nonnegative integer `n`, the center of the group algebra of
the symmetric group `S_n` has two standard bases indexed by the
integer partitions of `n`:

- The conjugacy class basis, where the basis element `C_\lambda`
  is the sum of all permutations in `S_n` with cycle type `\lambda`.

- The orthogonal idempotent basis, where the basis element
  `F_\lambda` is the orthogonal projection operator onto the
  irreducible representation of `S_n` indexed by `\lambda`.

This module implements these two bases, together with coercions and
conversions between them and other associated rings. In particular,
there is a surjective ring map from the ring of symmetric functions
to the center of the group algebra of the symmetric group, given
by evaluation at the Jucys-Murphy elements.

AUTHORS:

- Mathieu Guay-Paquet, Travis Scrimshaw (2013-10-08): initial version


EXAMPLES:

We begin by creating our bases of the center of the symmetric group
algebra and doing some basic computations::

    sage: ZS = SymmetricGroupAlgebraCenter(QQ, 4)
    sage: C = ZS.conjugacy_class_basis()
    sage: F = ZS.orthogonal_idempotent_basis()
    sage: C[2, 1, 1] * (3*C[4] + C[2, 2])
    C[2, 1, 1] + 12*C[2, 2] + 9*C[3, 1] + 2*C[4]
    sage: (F[3, 1] + F[4]) * (F[2, 2] - 4*F[3, 1])
    -4*F[3, 1]

We can convert between the two bases::

    sage: F(C[3, 1])
    8*F[1, 1, 1, 1] - 4*F[2, 2] + 8*F[4]
    sage: C(F[3, 1])
    3/8*C[1, 1, 1, 1] + 1/8*C[2, 1, 1] - 1/8*C[2, 2] - 1/8*C[4]

The center is of course embedded in the symmetric group algebra
itself::

    sage: ZS = SymmetricGroupAlgebraCenter(QQ, 3)
    sage: SGA = SymmetricGroupAlgebra(QQ, 3)
    sage: C = ZS.conjugacy_class_basis()
    sage: F = ZS.orthogonal_idempotent_basis()
    sage: SGA(C[2, 1])
    [1, 3, 2] + [2, 1, 3] + [3, 2, 1]
    sage: SGA(F[3] + F[1, 1, 1])
    1/3*[1, 2, 3] + 1/3*[2, 3, 1] + 1/3*[3, 1, 2]


COERCIONS:

There is always a coercion from the conjugacy class basis to the
orthonormal idempotent basis, but depending on the base ring, there
may be either a coercion or a conversion in the other direction
(from the orthogonal idempotent basis to the conjugacy class basis).

This is because the change of basis matrix in the first direction
has integer entries, but its inverse has rational entries with
non-trivial denominators in general. So, if the base ring doesn't
contain the rational numbers, the conversion may fail for certain
elements. In that case, we register a conversion from the F basis
to the C basis instead of a coercion.

An example of base ring where coercion works both ways::

    sage: ZS_coerce = SymmetricGroupAlgebraCenter(QQ, 3)
    sage: C = ZS_coerce.conjugacy_class_basis()
    sage: F = ZS_coerce.orthogonal_idempotent_basis()
    sage: x = C[2, 1]
    sage: y = 3*F[3] + 3*F[1, 1, 1]
    sage: z = F[2, 1] + 5*F[3]
    sage: F(x)
    -3*F[1, 1, 1] + 3*F[3]
    sage: C(y)
    C[1, 1, 1] + C[3]
    sage: C(z)
    3/2*C[1, 1, 1] + 5/6*C[2, 1] + 1/2*C[3]

For mixed arithmetic, the result lives in the left parent::

    sage: x + y
    C[1, 1, 1] + C[2, 1] + C[3]
    sage: y + x
    6*F[3]

And now for a base ring where conversion may fail::

    sage: ZS_convert = SymmetricGroupAlgebraCenter(ZZ, 3)
    sage: C = ZS_convert.conjugacy_class_basis()
    sage: F = ZS_convert.orthogonal_idempotent_basis()
    sage: x = C[2, 1]
    sage: y = 3*F[3] + 3*F[1, 1, 1]
    sage: z = F[2, 1] + 5*F[3]
    sage: F(x)
    -3*F[1, 1, 1] + 3*F[3]
    sage: C(y) # known bug
    C[1, 1, 1] + C[3]
    sage: C(z) # known bug
    Traceback (most recent call last):
    ...
    SomeError: conversion failed.

For mixed arithmetic, the result always lives in the F basis::

    sage: x + y
    6*F[3]
    sage: y + x
    6*F[3]

There is also a non-trivial coercion from symmetric functions to
the symmetric group algebra center given by evaluation at the
Jucys-Murphy elements. The `k`th Jucys-Murphy element is the sum
of all transpositions that swap the element `k` with a smaller
element::

    sage: QS5 = SymmetricGroupAlgebra(QQ, 5)
    sage: QS5.jucys_murphy(4)
    [1, 2, 4, 3, 5] + [1, 4, 3, 2, 5] + [4, 2, 3, 1, 5]

The Jucys-Murphy elements are elements of the symmetric group
algebra which commute with each other, but they are not actually
in the center. However, symmetric polynomials in the Jucys-Murphy
elements are in the center, and in fact the evaluation map is
surjective.

For example::

    sage: sym = SymmetricFunctions(QQ)
    sage: e = sym.e()
    sage: h = sym.h()
    sage: m = sym.m()
    sage: p = sym.p()
    sage: s = sym.s()
    sage: ZS = SymmetricGroupAlgebraCenter(QQ, 4)
    sage: C = ZS.C()
    sage: F = ZS.F()
    sage: C(e[2]) # known bug
    'TODO'
    sage: F(e[2]) # known bug
    'TODO'
    sage: C(h[1]) # known bug
    'TODO'
    sage: F(h[1]) # known bug
    'TODO'
    sage: C(m[4, 2]) # known bug
    'TODO'
    sage: F(m[4, 2]) # known bug
    'TODO'
    sage: C(p[3, 1]) # known bug
    'TODO'
    sage: F(p[3, 1]) # known bug
    'TODO'
    sage: C(s[5]) # known bug
    'TODO'
    sage: F(s[5]) # known bug
    'TODO'

REFERENCES:

 - A.-A. A. Jucys, Symmetric polynomials and the center of the
   symmetric group ring, Reports on Mathematical Physics, volume 5,
   number 1, 1974, pp. 107--112.

 - G. E. Murphy, A new construction of Young's seminormal
   representation of the symmetric groups, Journal of Algebra,
   volume 69, number 2, 1981, pp. 287--297.

 - P. Diaconis and C. Greene, Applications of Murphy's elements,
   Stanford Technical Report 335, 1989.

 - S. Corteel and A. Goupil and G. Schaeffer, Content evaluation and
   class symmetric functions, Advances in Mathematics, volume 188,
   number 2, 2004, pp. 315--336.

"""

#*****************************************************************************
#       Copyright (C) 2013 Mathieu Guay-Paquet <mathieu.guaypaquet@gmail.com>
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.structure.all import UniqueRepresentation, Parent
from sage.categories.all import Rings, AlgebrasWithBasis, Hom
from sage.categories.realizations import Realizations, Category_realization_of_parent
from sage.categories.morphism import SetMorphism
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition import Partition, Partitions
from sage.combinat.permutation import Permutations
from sage.combinat.sf.all import SymmetricFunctions
from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
from sage.rings.all import ZZ, QQ


########################
#   Helper functions   #
########################

dimension_cache = {}
def dimension(partition):
    r"""
    Return the dimension of the irreducible representation indexed by ``partition``.

    The dimensions of all irreducible representations for the appropriate size
    symmetric group are computed at once, then cached.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra_center import dimension
        sage: ps = Partitions(5)
        sage: ps.list()
        [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
        sage: [dimension(p) for p in ps]
        [1, 4, 5, 6, 5, 4, 1]
    """
    partition = Partition(partition)
    try:
        return dimension_cache[partition]
    except KeyError:
        sym = SymmetricFunctions(QQ)
        s = sym.schur()
        p = sym.powersum()
        dimensions = s(p(Partition([1] * partition.size())))
        dimension_cache.update(dimensions)
        return dimension_cache[partition]


###############
#   Parents   #
###############

class SymmetricGroupAlgebraCenter(UniqueRepresentation, Parent):
    r"""
    The center of the group algebra of a symmetric group.

    There are two important bases of the center, namely the conjugacy
    class basis and the orthogonal idempotent basis.

    .. RUBRIC: Conjugacy class basis

    The center of the group algebra of the symmetric group `S_n` has
    a basis indexed by all partitions of `n`, where the basis element
    indexed by a partition corresponds to the formal sum of all
    permutations in `S_n` with this cycle type.

    This basis is called the conjugacy class basis, and can be
    accessed via ``self.conjugacy_class_basis()`` or ``self.C()``.

    .. RUBRIC: Orthogonal idempotent basis

    The center of the group algebra of the symmetric group `S_n` has
    another basis indexed by all partitions of `n`, where the basis
    element indexed by a partition corresponds to a scaled version of
    the irreducible character of `S_n` indexed by this partition; the
    scaling is chosen so that the basis elements are idempotent.
    Equivalently, the basis elements are the orthogonal projection
    operators onto the irreducible representations of `S_n`.

    This basis is called the orthogonal idempotent basis, and can be
    accessed via ``self.orthogonal_idempotent_basis()`` or
    ``self.F()``.

    .. NOTE::

        Currently, this is only implemented for base rings with
        characteristic zero. Also, if the base ring doesn't contain
        the integers, then the algebra spanned by the orthogonal
        idempotent basis may be larger than the algebra spanned by
        the conjugacy class basis.

    INPUT:

    - ``base_ring`` -- base ring
    - ``n`` -- the order of the symmetric group

    OUTPUT:

    The center of the group algebra of the symmetric group `S_n`
    over ``base_ring``, in no particular basis.

    """
    def __init__(self, base_ring, n):
        r"""
        See ``SymmetricGroupAlgebraCenter`` for full documentation.

        TESTS::

            sage: ZS = SymmetricGroupAlgebraCenter(QQ, 5)
            sage: TestSuite(ZS).run()
            sage: ZS = SymmetricGroupAlgebraCenter(ZZ, 5) # known bug
            sage: TestSuite(ZS).run()

        """
        if base_ring not in Rings():
            raise TypeError('base_ring must be a ring.')
        if not base_ring.has_coerce_map_from(ZZ):
            # This tests what we actually use about characteristic zero
            raise NotImplementedError('base_ring must have characteristic zero.')

        self._n = n
        self._category = AlgebrasWithBasis(base_ring)
        Parent.__init__(
            self,
            base=base_ring,
            category=self._category.WithRealizations(),
            )

        # Register coercions between bases of self
        C = self.conjugacy_class_basis()
        F = self.orthogonal_idempotent_basis()
        C.module_morphism(
            on_basis=C._coercion_to_F,
            codomain=F,
            ).register_as_coercion()
        if base_ring.has_coerce_map_from(QQ):
            # In this case we know coercion always works
            F.module_morphism(
                on_basis=F._coercion_to_C,
                codomain=C,
                ).register_as_coercion()
        else:
            # Otherwise coercion may fail, so we register a conversion instead
            C.register_conversion(
                SetMorphism(
                    parent=Hom(F, C),
                    function=F._conversion_to_C,
                    ),
                )

#        # Register conversions from symmetric functions to self via the
#        # Jucys-Murphy evaluation map, also known as content evaluation
#        sym = SymmetricFunctions(base_ring)
#        sym.e().module_morphism(F._conversion_from_e, codomain=F).register_as_conversion()
#        sym.h().module_morphism(F._conversion_from_h, codomain=F).register_as_conversion()
#        sym.m().module_morphism(F._conversion_from_m, codomain=F).register_as_conversion()
#        sym.p().module_morphism(F._conversion_from_p, codomain=F).register_as_conversion()
#        sym.s().module_morphism(F._conversion_from_s, codomain=F).register_as_conversion()

        # Register coercion to the group algebra of the symmetric group
        SGA = SymmetricGroupAlgebra(base_ring, n)
        C.module_morphism(
            on_basis=C._coercion_to_sga,
            codomain=SGA,
            ).register_as_coercion()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: SymmetricGroupAlgebraCenter(QQ, 5)
            Center of the symmetric group algebra of order 5 over Rational Field

        """
        return "Center of the symmetric group algebra of order {} over {}".format(self._n, self.base_ring())

    def a_realization(self):
        r"""
        Return a realization of ``self``.

        EXAMPLES::

            sage: SymmetricGroupAlgebraCenter(QQ, 5).a_realization()
            Center of the symmetric group algebra of order 5 over Rational Field in the orthogonal idempotent basis

        """
        return self.F()

    # See sage.categories.sets_cat.Sets.WithRealizations.ParentMethods.inject_shorthands
    _shorthands = ('C', 'F')

    class conjugacy_class_basis(CombinatorialFreeModule, BindableClass):
        r"""
        The conjugacy class basis for the symmetric group algebra center.

        See ``SymmetricGroupAlgebraCenter`` for full documentation.

        TESTS::

            sage: SymmetricGroupAlgebraCenter(QQ, 5).conjugacy_class_basis()
            Center of the symmetric group algebra of order 5 over Rational Field in the conjugacy class basis

        """
        def __init__(self, ZS):
            r"""
            See ``SymmetricGroupAlgebraCenter`` for full documentation.

            TESTS::

                sage: ZS = SymmetricGroupAlgebraCenter(QQ, 5)
                sage: C = ZS.conjugacy_class_basis()
                sage: TestSuite(C).run()

            """
            self._basis_name = "conjugacy class"
            CombinatorialFreeModule.__init__(
                self,
                R=ZS.base_ring(),
                basis_keys=Partitions(ZS._n),
                prefix='C',
                bracket=False,
                category=SymmetricGroupAlgebraCenterBases(ZS),
                )

        @cached_method
        def _coercion_to_F(self, partition):
            r"""
            Return `C_{\lambda}` in the orthogonal idempotent basis.

            TESTS::

                sage: C = SymmetricGroupAlgebraCenter(QQ, 4).C()
                sage: C._coercion_to_F(Partition([3, 1]))
                8*F[1, 1, 1, 1] - 4*F[2, 2] + 8*F[4]

            """
            F = self.realization_of().F()
            P = Partitions(self.realization_of()._n)
            sym = SymmetricFunctions(self.base_ring())
            s = sym.schur()
            p = sym.powersum()

            ccs = partition.conjugacy_class_size()
            return F._from_dict(dict( (P(la), c * ccs / dimension(la)) for la,c in s(p[partition]) ))

        @cached_method
        def _coercion_to_sga(self, partition):
            r"""
            Lift `C_{\lambda}` to the symmetric group algebra.

            TESTS::

                sage: C = SymmetricGroupAlgebraCenter(QQ, 3).C()
                sage: C._coercion_to_sga(Partition([3]))
                [2, 3, 1] + [3, 1, 2]
                sage: C._coercion_to_sga(Partition([2,1]))
                [1, 3, 2] + [2, 1, 3] + [3, 2, 1]
                sage: C._coercion_to_sga(Partition([1,1,1]))
                [1, 2, 3]

            """
            n = self.realization_of()._n
            SGA = SymmetricGroupAlgebra(self.base_ring(), n)
            return SGA._from_dict(dict((p, 1) for p in Permutations(n) if p.cycle_type() == partition))

        def product(self, left, right):
            r"""
            Compute the product of two elements of the ring.

            This is done by coercing the elements to the orthogonal idempotent
            basis, computing the product there, and coercing the result back.

            INPUT:

            - ``left`` -- an element of the ring
            - ``right`` -- an element of the ring

            OUTPUT:

            The product of the two elements.

            EXAMPLES::

                sage: C = SymmetricGroupAlgebraCenter(QQ, 5).C()
                sage: C.product(C[3,1,1], C[3,1,1])
                20*C[1, 1, 1, 1, 1] + 8*C[2, 2, 1] + 7*C[3, 1, 1] + 5*C[5]
                sage: C[3,1,1] * C[3,1,1]
                20*C[1, 1, 1, 1, 1] + 8*C[2, 2, 1] + 7*C[3, 1, 1] + 5*C[5]
                sage: C[3,1,1] * C[2,2,1]
                4*C[2, 2, 1] + 6*C[3, 1, 1] + 5*C[5]
                sage: C[5] * C[2,1,1,1]
                6*C[3, 2] + 4*C[4, 1]

            """
            F = self.realization_of().orthogonal_idempotent_basis()
            return self(F(left) * F(right))

        def one_basis(self):
            r"""
            Return the index of the identity element.

            EXAMPLES::

                sage: C = SymmetricGroupAlgebraCenter(QQ, 5).C()
                sage: C.one_basis()
                [1, 1, 1, 1, 1]
                sage: C.one_basis().parent()
                Partitions of the integer 5

            """
            n = self.realization_of()._n
            P = Partitions(n)
            return P([1]*n)

    C = conjugacy_class_basis

    class orthogonal_idempotent_basis(CombinatorialFreeModule, BindableClass):
        r"""
        The orthogonal idempotent basis for the symmetric group algebra center.

        See ``SymmetricGroupAlgebraCenter`` for full documentation.

        TESTS::

            sage: SymmetricGroupAlgebraCenter(QQ, 5).orthogonal_idempotent_basis()
            Center of the symmetric group algebra of order 5 over Rational Field in the orthogonal idempotent basis

        """
        def __init__(self, ZS):
            r"""
            See ``SymmetricGroupAlgebraCenter`` for full documentation.

            TESTS::

                sage: ZS = SymmetricGroupAlgebraCenter(QQ, 5)
                sage: F = ZS.orthogonal_idempotent_basis()
                sage: TestSuite(F).run()

            """
            self._basis_name = "orthogonal idempotent"
            CombinatorialFreeModule.__init__(
                self,
                R=ZS.base_ring(),
                basis_keys=Partitions(ZS._n),
                prefix='F',
                bracket=False,
                category=SymmetricGroupAlgebraCenterBases(ZS),
                )

        def _conversion_to_C(self, element):
            r"""
            Return the ``element`` of ``self`` in the conjugacy class basis.

            This method is used when the base ring for the symmetric
            group algebra center doesn't contain the rational numbers.
            Some elements which have only integer coefficients in the
            orthogonal idempotent basis have coefficients with
            non-trivial denominators in the conjugacy class basis,
            and the conversion fails in those cases.

            Since the conversion may fail, it cannot be used as a
            coercion, only as a conversion.

            .. SEE ALSO::

                :meth:`~SymmetricGroupAlgebraCenter.orthogonal_idempotent_basis._coercion_to_C`

            TESTS::

                sage: ZS = SymmetricGroupAlgebraCenter(ZZ, 4)
                sage: C = ZS.conjugacy_class_basis()
                sage: F = ZS.orthogonal_idempotent_basis()
                sage: x = 2*F[1, 1, 1, 1] + 2*F[4]
                sage: y = F[3, 1]
                sage: F._conversion_to_C(x) # known bug
                'TODO'
                sage: C(x) # known bug
                'TODO'
                sage: F._conversion_to_C(y) # known bug
                'TODO'
                sage: C(y) # known bug
                'TODO'

            """
            raise NotImplementedError()

        @cached_method
        def _coercion_to_C(self, partition):
            r"""
            Return `F_{\lambda}` in the conjugacy class basis.

            This method is used when the base ring for the symmetric
            group algebra center contains the rational numbers.
            When this is the case, the conversion is guaranteed to
            work, so it can be used as a coercion.

            .. SEE ALSO::

                :meth:`~SymmetricGroupAlgebraCenter.orthogonal_idempotent_basis._conversion_to_C`

            TESTS::

                sage: ZS = SymmetricGroupAlgebraCenter(QQ, 4)
                sage: C = ZS.conjugacy_class_basis()
                sage: F = ZS.orthogonal_idempotent_basis()
                sage: x = 2*F[1, 1, 1, 1] + 2*F[4]
                sage: y = F[3, 1]
                sage: C(x)
                1/6*C[1, 1, 1, 1] + 1/6*C[2, 2] + 1/6*C[3, 1]
                sage: C(y)
                3/8*C[1, 1, 1, 1] + 1/8*C[2, 1, 1] - 1/8*C[2, 2] - 1/8*C[4]
                sage: F._coercion_to_C(Partition([3, 1]))
                3/8*C[1, 1, 1, 1] + 1/8*C[2, 1, 1] - 1/8*C[2, 2] - 1/8*C[4]

            """
            C = self.realization_of().C()
            P = Partitions(self.realization_of()._n)
            sym = SymmetricFunctions(self.base_ring())
            s = sym.schur()
            p = sym.powersum()
            dim = dimension(partition)
            return C._from_dict(dict((P(la), c * dim / la.conjugacy_class_size()) for la,c in p(s[partition])))

        def product(self, left, right):
            r"""
            Compute the product of two elements of the ring.

            The basis consists of orthogonal idempotents, so the product of
            distinct basis elements is zero, and the product of a basis element
            with itself is itself.

            INPUT:

            - ``left`` -- an element of the ring
            - ``right`` -- an element of the ring

            OUTPUT:

            The product of the two elements.

            EXAMPLES::

                sage: F = SymmetricGroupAlgebraCenter(QQ, 5).F()
                sage: F.product(F[3, 1, 1], F[3, 1, 1])
                F[3, 1, 1]
                sage: F[3, 1, 1] * F[3, 1, 1]
                F[3, 1, 1]
                sage: F[3, 1, 1] * F[2, 2, 1]
                0
                sage: (3*F[3, 1, 1] - F[4, 1]) * (F[5] + 2*F[4, 1])
                -2*F[4, 1]

            """
            return self._from_dict(dict((i, c * right[i]) for (i, c) in left if c * right[i] != 0))

        @cached_method
        def one(self):
            r"""
            Return the identity element.

            EXAMPLES::

                sage: F = SymmetricGroupAlgebraCenter(QQ, 3).F()
                sage: F.one()
                F[1, 1, 1] + F[2, 1] + F[3]

            """
            n = self.realization_of()._n
            return self._from_dict(dict((la, 1) for la in Partitions(self.realization_of()._n)))

    F = orthogonal_idempotent_basis


##################
#   Categories   #
##################

class SymmetricGroupAlgebraCenterBases(Category_realization_of_parent):
    r"""
    The category of bases of ``SymmetricGroupAlgebraCenter``.

    INPUT:

    - ``base`` -- the center of the a symmetric group algebra

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra_center import SymmetricGroupAlgebraCenterBases
        sage: ZS = SymmetricGroupAlgebraCenter(QQ, 5)
        sage: bases = SymmetricGroupAlgebraCenterBases(ZS)

    """
    def __init__(self, base):
        r"""
        See ``SymmetricGroupAlgebraCenterBases`` for full documentation.

        TESTS::

            sage: from sage.combinat.symmetric_group_algebra_center import SymmetricGroupAlgebraCenterBases
            sage: ZS = SymmetricGroupAlgebraCenter(QQ, 5)
            sage: bases = SymmetricGroupAlgebraCenterBases(ZS)
            sage: ZS.C() in bases
            True
            sage: ZS.F() in bases
            True
            sage: TestSuite(bases).run()

        """
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        TESTS::

            sage: from sage.combinat.symmetric_group_algebra_center import SymmetricGroupAlgebraCenterBases
            sage: ZS = SymmetricGroupAlgebraCenter(QQ, 5)
            sage: SymmetricGroupAlgebraCenterBases(ZS)
            Category of bases of Center of the symmetric group algebra of order 5 over Rational Field

        """
        return "Category of bases of {}".format(self.base())

    def super_categories(self):
        r"""
        Return the list of super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.symmetric_group_algebra_center import SymmetricGroupAlgebraCenterBases
            sage: ZS = SymmetricGroupAlgebraCenter(QQ, 5)
            sage: bases = SymmetricGroupAlgebraCenterBases(ZS)
            sage: bases.super_categories()
            [Category of algebras with basis over Rational Field,
             Category of realizations of Center of the symmetric group algebra of order 5 over Rational Field]

        """
        return [self.base()._category, Realizations(self.base())]

    class ParentMethods:
        def _repr_(self):
            r"""
            Return a string representation of ``self``.

            TESTS::

                sage: ZS = SymmetricGroupAlgebraCenter(QQ, 5)
                sage: ZS.C()
                Center of the symmetric group algebra of order 5 over Rational Field in the conjugacy class basis
                sage: ZS.F()
                Center of the symmetric group algebra of order 5 over Rational Field in the orthogonal idempotent basis

            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def __getitem__(self, partition):
            r"""
            Return the basis element indexed by ``partition``.

            INPUT:

            - ``partition`` -- a partition

            EXAMPLES::

                sage: C = SymmetricGroupAlgebraCenter(QQ, 5).C()
                sage: C[5]
                C[5]
                sage: C[3, 1, 1]
                C[3, 1, 1]
                sage: C[[3, 1, 1]]
                C[3, 1, 1]
                sage: C[Partition([3, 1, 1])]
                C[3, 1, 1]

            """
            P = Partitions(self.realization_of()._n)
            try:
                key = P(partition)
            except ValueError:
                key = P([partition])
            return self.monomial(key)

        def is_field(self, proof=True):
            r"""
            Return whether ``self`` is a field.

            EXAMPLES::

                sage: ZS = SymmetricGroupAlgebraCenter(QQ, 5)
                sage: ZS.C().is_field()
                False
                sage: ZS.F().is_field()
                False

            """
            return False

        def is_commutative(self):
            r"""
            Return whether ``self`` is commutative.

            EXAMPLES::

                sage: ZS = SymmetricGroupAlgebraCenter(QQ, 5)
                sage: ZS.C().is_commutative()
                True
                sage: ZS.F().is_commutative()
                True

            """
            return True


