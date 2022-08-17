"""
Cartesian products

AUTHORS:

- Nicolas Thiery (2010-03): initial version
"""
# ****************************************************************************
#       Copyright (C) 2008-2014 Nicolas Thiery <nthiery at users.sf.net>,
#                     2008      Mike Hansen <mhansen@gmail.com>,
#                     2008      Florent Hivert <Florent.Hivert@univ-rouen.fr>
#                     2014      Nathann Cohen
#                     2014      Peter Bruin
#                     2015      Travis Scrimshaw
#                     2015      Jori Mäntysalo
#                     2015      Daniel Krenn
#                     2015      David Roe
#                     2015-2019 Vincent Delecroix
#                     2016      Johan S. R. Nielsen
#                     2022      Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import numbers

from sage.misc.call import attrcall
from sage.misc.cachefunc import cached_method

from sage.categories.sets_cat import Sets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

from sage.structure.parent import Parent
from sage.structure.unique_representation import WithPicklingByInitArgs, UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapperCheckWrappedClass

from sage.sets.set_from_iterator import EnumeratedSetFromIterator

from sage.categories.rings import Rings
_Rings = Rings()


class CartesianProduct_base(WithPicklingByInitArgs, Parent):
    """
    A class implementing a raw data structure for Cartesian products
    of sets (and elements thereof). See :obj:`cartesian_product` for
    how to construct full fledged Cartesian products.

    EXAMPLES::

        sage: G = cartesian_product([GF(5), Permutations(10)])
        sage: G.cartesian_factors()
        (Finite Field of size 5, Standard permutations of 10)
        sage: G.cardinality()
        18144000
        sage: G.random_element()    # random
        (1, [4, 7, 6, 5, 10, 1, 3, 2, 8, 9])
        sage: G.category()
        Join of Category of finite monoids
            and Category of Cartesian products of monoids
            and Category of Cartesian products of finite enumerated sets

    .. automethod:: CartesianProduct._cartesian_product_of_elements
    """

    def __init__(self, sets, category, flatten=False):
        r"""
        INPUT:

         - ``sets`` -- a tuple of parents
         - ``category`` -- a subcategory of ``Sets().CartesianProducts()``
         - ``flatten`` -- a boolean (default: ``False``)

        ``flatten`` is current ignored, and reserved for future use.

        No other keyword arguments (``kwargs``) are accepted.

        TESTS::

            sage: from sage.sets.cartesian_product import CartesianProduct
            sage: C = CartesianProduct((QQ, ZZ, ZZ), category = Sets().CartesianProducts())
            sage: C
            The Cartesian product of (Rational Field, Integer Ring, Integer Ring)
            sage: C.an_element()
            (1/2, 1, 1)
            sage: TestSuite(C).run()
            sage: cartesian_product([ZZ, ZZ], blub=None)
            Traceback (most recent call last):
            ...
            TypeError: ...__classcall_private__() got an unexpected keyword argument 'blub'
        """
        self._sets = tuple(sets)
        Parent.__init__(self, category=category)

    def _element_constructor_(self,x):
        r"""
        Construct an element of a Cartesian product from a list or iterable

        INPUT:

        - ``x`` -- a list (or iterable)

        Each component of `x` is converted to the corresponding
        Cartesian factor.

        EXAMPLES::

            sage: C = cartesian_product([GF(5), GF(3)])
            sage: x = C((1,3)); x
            (1, 0)
            sage: x.parent()
            The Cartesian product of (Finite Field of size 5, Finite Field of size 3)
            sage: x[0].parent()
            Finite Field of size 5
            sage: x[1].parent()
            Finite Field of size 3

        An iterable is also accepted as input::

            sage: C(i for i in range(2))
            (0, 1)

        TESTS::

            sage: C((1,3,4))
            Traceback (most recent call last):
            ...
            ValueError: (1, 3, 4) should be of length 2

            sage: R = ZZ.cartesian_product(ZZ)
            sage: R(0)
            (0, 0)
            sage: R(-5)
            (-5, -5)
        """
        # NOTE: should we more generally allow diagonal embedding
        # if we have a conversion?
        if self in _Rings and isinstance(x, numbers.Integral):
            return x * self.one()

        from builtins import zip
        x = tuple(x)

        if len(x) != len(self._sets):
            raise ValueError(
                "{} should be of length {}".format(x, len(self._sets)))
        x = tuple(c(xx) for c, xx in zip(self._sets, x))
        return self.element_class(self, x)

    def _repr_(self):
        """
        EXAMPLES::

            sage: cartesian_product([QQ, ZZ, ZZ]) # indirect doctest
            The Cartesian product of (Rational Field, Integer Ring, Integer Ring)
        """
        return "The Cartesian product of %s"%(self._sets,)

    def __contains__(self, x):
        """
        Check if ``x`` is contained in ``self``.

        EXAMPLES::

            sage: C = cartesian_product([list(range(5)), list(range(5))])
            sage: (1, 1) in C
            True
            sage: (1, 6) in C
            False
        """
        if isinstance(x, self.Element):
            if x.parent() == self:
                return True
        elif not isinstance(x, tuple):
            return False
        return ( len(x) == len(self._sets)
                 and all(elt in self._sets[i] for i,elt in enumerate(x)) )

    def cartesian_factors(self):
        """
        Return the Cartesian factors of ``self``.

        .. SEEALSO::

            :meth:`Sets.CartesianProducts.ParentMethods.cartesian_factors()
            <sage.categories.sets_cat.Sets.CartesianProducts.ParentMethods.cartesian_factors>`.

        EXAMPLES::

            sage: cartesian_product([QQ, ZZ, ZZ]).cartesian_factors()
            (Rational Field, Integer Ring, Integer Ring)
        """
        return self._sets

    def _sets_keys(self):
        """
        Return the indices of the Cartesian factors of ``self``
        as per
        :meth:`Sets.CartesianProducts.ParentMethods._sets_keys()
        <sage.categories.sets_cat.Sets.CartesianProducts.ParentMethods._sets_keys>`.

        EXAMPLES::

            sage: cartesian_product([QQ, ZZ, ZZ])._sets_keys()
            {0, 1, 2}
            sage: cartesian_product([ZZ]*100)._sets_keys()
            {0, ..., 99}
        """
        from sage.sets.integer_range import IntegerRange
        return IntegerRange(len(self._sets))

    @cached_method
    def cartesian_projection(self, i):
        """
        Return the natural projection onto the `i`-th Cartesian
        factor of ``self`` as per
        :meth:`Sets.CartesianProducts.ParentMethods.cartesian_projection()
        <sage.categories.sets_cat.Sets.CartesianProducts.ParentMethods.cartesian_projection>`.

        INPUT:

        - ``i`` -- the index of a Cartesian factor of ``self``

        EXAMPLES::

            sage: C = Sets().CartesianProducts().example(); C
            The Cartesian product of (Set of prime numbers (basic implementation), An example of an infinite enumerated set: the non negative integers, An example of a finite enumerated set: {1,2,3})
            sage: x = C.an_element(); x
            (47, 42, 1)
            sage: pi = C.cartesian_projection(1)
            sage: pi(x)
            42

            sage: C.cartesian_projection('hey')
            Traceback (most recent call last):
            ...
            ValueError: i (=hey) must be in {0, 1, 2}
        """
        if i not in self._sets_keys():
            raise ValueError("i (={}) must be in {}".format(i, self._sets_keys()))
        return attrcall("cartesian_projection", i)

    def _cartesian_product_of_elements(self, elements):
        """
        Return the Cartesian product of the given ``elements``.

        This implements :meth:`Sets.CartesianProducts.ParentMethods._cartesian_product_of_elements`.
        INPUT:

        - ``elements`` -- an iterable (e.g. tuple, list) with one element of
          each Cartesian factor of ``self``

        .. WARNING::

            This is meant as a fast low-level method. In particular,
            no coercion is attempted. When coercion or sanity checks
            are desirable, please use instead ``self(elements)`` or
            ``self._element_constructor_(elements)``.

        EXAMPLES::

            sage: S1 = Sets().example()
            sage: S2 = InfiniteEnumeratedSets().example()
            sage: C = cartesian_product([S2, S1, S2])
            sage: C._cartesian_product_of_elements([S2.an_element(), S1.an_element(), S2.an_element()])
            (42, 47, 42)
        """
        elements = tuple(elements)
        assert len(elements) == len(self._sets)
        return self.element_class(self, elements)

    def construction(self):
        r"""
        Return the construction functor and its arguments for this
        Cartesian product.

        OUTPUT:

        A pair whose first entry is a Cartesian product functor and
        its second entry is a list of the Cartesian factors.

        EXAMPLES::

            sage: cartesian_product([ZZ, QQ]).construction()
            (The cartesian_product functorial construction,
             (Integer Ring, Rational Field))
        """
        from sage.categories.cartesian_product import CartesianProductFunctor
        return CartesianProductFunctor(self.category()), self.cartesian_factors()

    def _coerce_map_from_(self, S):
        r"""
        Return ``True`` if ``S`` coerces into this Cartesian product.

        TESTS::

            sage: Z = cartesian_product([ZZ])
            sage: Q = cartesian_product([QQ])
            sage: Z.has_coerce_map_from(Q)  # indirect doctest
            False
            sage: Q.has_coerce_map_from(Z)  # indirect doctest
            True
        """
        if isinstance(S, CartesianProduct_base):
            S_factors = S.cartesian_factors()
            R_factors = self.cartesian_factors()
            if len(S_factors) == len(R_factors):
                if all(r.has_coerce_map_from(s) for r, s in zip(R_factors, S_factors)):
                    return True
        return super()._coerce_map_from_(S)

    an_element = Sets.CartesianProducts.ParentMethods.an_element


class CartesianProduct_with_element_wrapper(CartesianProduct_base):

    class Element(ElementWrapperCheckWrappedClass):

        wrapped_class = tuple

        def cartesian_projection(self, i):
            r"""
            Return the projection of ``self`` on the `i`-th factor of
            the Cartesian product, as per
            :meth:`Sets.CartesianProducts.ElementMethods.cartesian_projection()
            <sage.categories.sets_cat.Sets.CartesianProducts.ElementMethods.cartesian_projection>`.

            INPUT:

            - ``i`` -- the index of a factor of the Cartesian product

            EXAMPLES::

                sage: C = Sets().CartesianProducts().example(); C
                The Cartesian product of (Set of prime numbers (basic implementation), An example of an infinite enumerated set: the non negative integers, An example of a finite enumerated set: {1,2,3})
                sage: x = C.an_element(); x
                (47, 42, 1)
                sage: x.cartesian_projection(1)
                42
            """
            return self.value[i]

        __getitem__ = cartesian_projection

        def __iter__(self):
            r"""
            Iterate over the components of an element.

            EXAMPLES::

                sage: C = Sets().CartesianProducts().example(); C
                The Cartesian product of
                (Set of prime numbers (basic implementation),
                 An example of an infinite enumerated set: the non negative integers,
                 An example of a finite enumerated set: {1,2,3})
                sage: c = C.an_element(); c
                (47, 42, 1)
                sage: for i in c:
                ....:     print(i)
                47
                42
                1
            """
            return iter(self.value)

        def __len__(self):
            r"""
            Return the number of factors in the cartesian product from which ``self`` comes.

            EXAMPLES::

                sage: C = cartesian_product([ZZ, QQ, CC])
                sage: e = C.random_element()
                sage: len(e)
                3
            """
            return len(self.value)

        def cartesian_factors(self):
            r"""
            Return the tuple of elements that compose this element.

            EXAMPLES::

                sage: A = cartesian_product([ZZ, RR])
                sage: A((1, 1.23)).cartesian_factors()
                (1, 1.23000000000000)
                sage: type(_)
                <... 'tuple'>
            """
            return self.value


class CartesianProduct_eq_by_factors(CartesianProduct_with_element_wrapper):
    def __eq__(self, other):
        if isinstance(other, CartesianProduct_base):
            # No flattening, hence we are equal if and only if our factors are equal
            return self.cartesian_factors() == other.cartesian_factors()
        return super().__eq__(other)

    def __hash__(self):
        # No flattening, hence we are equal if and only if our factors are equal
        return hash(self.cartesian_factors())


class CartesianProduct_unique(UniqueRepresentation, CartesianProduct_with_element_wrapper):
    r"""
    A Cartesian product with :class:`~sage.structure.unique_representation.UniqueRepresentation` behavior.
    """
    pass


class CartesianProduct(CartesianProduct_base):

    @staticmethod
    def __classcall_private__(cls, sets, category, flatten=False):
        r"""
        Suppress :class:`~sage.structure.unique_representation.UniqueRepresentation` behavior for certain infinite factors.

        Two :class:`~sage.sets.set_from_iterator.EnumeratedSetFromIterator` objects that are not known to be finite
        cannot be reliably tested for equality. Therefore, we do not put such objects
        in the :class:`~sage.structure.unique_representation.UniqueRepresentation` cache.

        EXAMPLES::

            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: F = EnumeratedSetFromIterator(lambda: iter([1, 2]))
            sage: F.category()
            Category of facade enumerated sets
            sage: cartesian_product([F, F]) == cartesian_product([F, F])
            True
            sage: cartesian_product([F, F]) is cartesian_product([F, F])
            False
            sage: G = EnumeratedSetFromIterator(lambda: iter([1, 2]),
            ....:                               category=FiniteEnumeratedSets())
            sage: G.category()
            Category of facade finite enumerated sets
            sage: cartesian_product([G, G]) == cartesian_product([G, G])
            True
            sage: cartesian_product([G, G]) is cartesian_product([G, G])
            True
        """
        assert cls == CartesianProduct
        # Trac #19195: EnumeratedSetFromIterator instances are not safe to be passed
        # to UniqueRepresentation because EnumeratedSetFromIterator.__eq__ resorts
        # to a semi-decision procedure for equality.
        if not any(isinstance(set, EnumeratedSetFromIterator) and set not in FiniteEnumeratedSets()
                   for set in sets):
            # UniqueRepresentation is safe to use.
            return CartesianProduct_unique(sets, category)
        else:
            return CartesianProduct_eq_by_factors(sets, category)
