r"""
Unital algebras
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.rings import Rings
from sage.categories.magmas import Magmas
from sage.categories.magmatic_algebras import MagmaticAlgebras

class UnitalAlgebras(CategoryWithAxiom_over_base_ring):
    """
    The category of non-associative algebras over a given base ring.

    A non-associative algebra over a ring `R` is a module over `R`
    which s also a unital magma.

    .. WARNING::

        Until :trac:`15043` is implemented, :class:`Algebras` is the
        category of associative unital algebras; thus, unlike the name
        suggests, :class:`UnitalAlgebras` is not a subcategory of
        :class:`Algebras` but of
        :class:`~.magmatic_algebras.MagmaticAlgebras`.

    EXAMPLES::

        sage: from sage.categories.unital_algebras import UnitalAlgebras
        sage: C = UnitalAlgebras(ZZ); C
        Category of unital algebras over Integer Ring

    TESTS::

        sage: from sage.categories.magmatic_algebras import MagmaticAlgebras
        sage: C is MagmaticAlgebras(ZZ).Unital()
        True
        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (MagmaticAlgebras, "Unital")

    class ParentMethods:
        def from_base_ring(self, r):
            """
            Return the canonical embedding of ``r`` into ``self``.

            INPUT:

            - ``r`` -- an element of ``self.base_ring()``

            EXAMPLES::

                sage: A = AlgebrasWithBasis(QQ).example(); A
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: A.from_base_ring(1)
                B[word: ]
            """
            return self.one()._lmul_(r)

        def __init_extra__(self):
            """
            Declare the canonical coercion from ``self.base_ring()``
            to ``self``, if there has been none before.

            EXAMPLES::

                sage: A = AlgebrasWithBasis(QQ).example(); A
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: coercion_model = sage.structure.element.get_coercion_model()
                sage: coercion_model.discover_coercion(QQ, A)
                (Generic morphism:
                  From: Rational Field
                  To:   An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field, None)
                sage: A(1)          # indirect doctest
                B[word: ]

            """
            # If self has an attribute _no_generic_basering_coercion
            # set to True, then this declaration is skipped.
            # This trick, introduced in #11900, is used in
            # sage.matrix.matrix_space.py and
            # sage.rings.polynomial.polynomial_ring.
            # It will hopefully be refactored into something more
            # conceptual later on.
            if getattr(self, '_no_generic_basering_coercion', False):
                return

            base_ring = self.base_ring()
            if base_ring is self:
                # There are rings that are their own base rings. No need to register that.
                return
            if self._is_coercion_cached(base_ring):
                # We will not use any generic stuff, since a (presumably) better conversion
                # has already been registered.
                return

            # This could be a morphism of Algebras(self.base_ring()); however, e.g., QQ is not in Algebras(QQ)
            H = Hom(base_ring, self, Rings()) # TODO: non associative ring!

            # We need to register a coercion from the base ring to self.
            #
            # There is a generic method from_base_ring(), that just does
            # multiplication with the multiplicative unit. However, the
            # unit is constructed repeatedly, which is slow.
            # So, if the unit is available *now*, then we can create a
            # faster coercion map.
            #
            # This only applies for the generic from_base_ring() method.
            # If there is a specialised from_base_ring(), then it should
            # be used unconditionally.
            generic_from_base_ring = self.category().parent_class.from_base_ring
            if type(self).from_base_ring != generic_from_base_ring:
                # Custom from_base_ring()
                use_from_base_ring = True
            if isinstance(generic_from_base_ring, lazy_attribute):
                # If the category implements from_base_ring() as lazy
                # attribute, then we always use it.
                # This is for backwards compatibility, see Trac #25181
                use_from_base_ring = True
            else:
                try:
                    one = self.one()
                    use_from_base_ring = False
                except (NotImplementedError, AttributeError, TypeError):
                    # The unit is not available, yet. But there are cases
                    # in which it will be available later. So, we use
                    # the generic from_base_ring() after all.
                    use_from_base_ring = True

            mor = None
            if use_from_base_ring:
                mor = SetMorphism(function=self.from_base_ring, parent=H)
            else:
                # We have the multiplicative unit, so implement the
                # coercion from the base ring as multiplying with that.
                #
                # But first we check that it actually works. If not,
                # then the generic implementation of from_base_ring()
                # would fail as well so we don't use it.
                try:
                    if one._lmul_(base_ring.an_element()) is not None:
                        # There are cases in which lmul returns None,
                        # which means that it's not implemented.
                        # One example: Hecke algebras.
                        mor = SetMorphism(function=one._lmul_, parent=H)
                except (NotImplementedError, AttributeError, TypeError):
                    pass
            if mor is not None:
                try:
                    self.register_coercion(mor)
                except AssertionError:
                    pass

    class WithBasis(CategoryWithAxiom_over_base_ring):

        class ParentMethods:

            @abstract_method(optional = True)
            def one_basis(self):
                """
                When the one of an algebra with basis is an element of
                this basis, this optional method can return the index of
                this element. This is used to provide a default
                implementation of :meth:`.one`, and an optimized default
                implementation of :meth:`.from_base_ring`.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.one_basis()
                    word:
                    sage: A.one()
                    B[word: ]
                    sage: A.from_base_ring(4)
                    4*B[word: ]
                """

            @cached_method
            def one_from_one_basis(self):
                """
                Return the one of the algebra, as per
                :meth:`Monoids.ParentMethods.one()
                <sage.categories.monoids.Monoids.ParentMethods.one>`

                By default, this is implemented from
                :meth:`.one_basis`, if available.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.one_basis()
                    word:
                    sage: A.one_from_one_basis()
                    B[word: ]
                    sage: A.one()
                    B[word: ]

                TESTS:

                Try to check that :trac:`5843` Heisenbug is fixed::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: B = AlgebrasWithBasis(QQ).example(('a', 'c'))
                    sage: A == B
                    False
                    sage: Aone = A.one_from_one_basis
                    sage: Bone = B.one_from_one_basis
                    sage: Aone is Bone
                    False

               Even if called in the wrong order, they should returns their
               respective one::

                    sage: Bone().parent() is B
                    True
                    sage: Aone().parent() is A
                    True
                """
                return self.monomial(self.one_basis()) #.

            @lazy_attribute
            def one(self):
                r"""
                Return the multiplicative unit element.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.one_basis()
                    word:
                    sage: A.one()
                    B[word: ]
                """
                if self.one_basis is NotImplemented:
                    return NotImplemented
                return self.one_from_one_basis

            @lazy_attribute
            def from_base_ring(self):
                """
                TESTS::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.from_base_ring(3)
                    3*B[word: ]
                """
                if self.one_basis is NotImplemented:
                    return NotImplemented
                return self.from_base_ring_from_one_basis

            def from_base_ring_from_one_basis(self, r):
                """
                Implement the canonical embedding from the ground ring.

                INPUT:

                - ``r`` -- an element of the coefficient ring

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.from_base_ring_from_one_basis(3)
                    3*B[word: ]
                    sage: A.from_base_ring(3)
                    3*B[word: ]
                    sage: A(3)
                    3*B[word: ]
                """
                return self.term(self.one_basis(), r)
