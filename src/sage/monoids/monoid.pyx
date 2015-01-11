r"""
Monoids
"""

from sage.structure.parent cimport Parent
from sage.misc.cachefunc import cached_method


def is_Monoid(x):
    r"""
    Returns True if ``x`` is of type ``Monoid_class``.

    EXAMPLES::

        sage: from sage.monoids.monoid import is_Monoid
        sage: is_Monoid(0)
        False
        sage: is_Monoid(ZZ)   # The technical math meaning of monoid has
        ...                   # no bearing whatsoever on the result: it's
        ...                   # a typecheck which is not satisfied by ZZ
        ...                   # since it does not inherit from Monoid_class.
        False
        sage: is_Monoid(sage.monoids.monoid.Monoid_class(('a','b')))
        True
        sage: F.<a,b,c,d,e> = FreeMonoid(5)
        sage: is_Monoid(F)
        True
    """
    return isinstance(x, Monoid_class)


class Monoid_class(Parent):
    def __init__(self,names):
        r"""
        EXAMPLES::

            sage: from sage.monoids.monoid import Monoid_class
            sage: Monoid_class(('a','b'))
            <class 'sage.monoids.monoid.Monoid_class_with_category'>

        TESTS::

            sage: F.<a,b,c,d,e> = FreeMonoid(5)
            sage: TestSuite(F).run()
        """
        from sage.categories.monoids import Monoids
        Parent.__init__(self, base=self,names=names,category=Monoids())

    @cached_method
    def gens(self):
        r"""
        Returns the generators for ``self``.

        EXAMPLES::

            sage: F.<a,b,c,d,e> = FreeMonoid(5)
            sage: F.gens()
            (a, b, c, d, e)
        """
        return tuple(self.gen(i) for i in range(self.ngens()))


cdef class Monoid(Parent):
    """
    Base class for all monoids

    TESTS::

        sage: from sage.monoids.monoid import Monoid
        sage: G = Monoid()
        sage: TestSuite(G).run(skip = ["_test_an_element",\
                                       "_test_associativity",\
                                       "_test_elements",\
                                       "_test_elements_eq_reflexive",\
                                       "_test_elements_eq_symmetric",\
                                       "_test_elements_eq_transitive",\
                                       "_test_elements_neq",\
                                       "_test_inverse",\
                                       "_test_one",\
                                       "_test_pickling",\
                                       "_test_prod",\
                                       "_test_some_elements"])
    """
    def __init__(self, base=None, gens=None, category=None):
        """
        The Python constructor

        TESTS::

            sage: from sage.monoids.monoid import Monoid
            sage: G = Monoid()
            sage: G.category()
            Category of monoids
            sage: G = Monoid(category=Monoids()) # todo: do the same test with some subcategory of Monoids when there will exist one
            sage: G.category()
            Category of monoids
            sage: G = Monoid(category = CommutativeAdditiveMonoids())
            Traceback (most recent call last):
            ...
            ValueError: (Category of commutative additive monoids,) is not a subcategory of Category of monoids
            sage: G._repr_option('element_is_atomic')
            False
        """
        from sage.categories.monoids import Monoids
        if category is None:
            category = Monoids()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(Monoids()) for cat in category):
                raise ValueError("%s is not a subcategory of %s"%(category, Monoids()))
        Parent.__init__(self, base=base, gens=gens, category=category)


    def __contains__(self, x):
        r"""
        Test whether `x` defines a monoid element.

        INPUT:

        - ``x`` -- anything.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.monoids.monoid import Monoid
            sage: G = Monoid()
            sage: 4 in G               #indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        try:
            self(x)
        except TypeError:
            return False
        return True


    def is_abelian(self):
        """
        Test whether this monoid is abelian.

        EXAMPLES::

            sage: from sage.monoids.monoid import Monoid
            sage: G = Monoid()
            sage: G.is_abelian()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def is_commutative(self):
        r"""
        Test whether this monoid is commutative.

        This is an alias for is_abelian, largely to make monoids work
        well with the Factorization class.

        (Note for developers: Derived classes should override is_abelian, not
        is_commutative.)

        EXAMPLE::

            sage: SL(2, 7).is_commutative()
            False
        """
        return self.is_abelian()


    def order(self):
        """
        Returns the number of elements of this monoid, which is either a
        positive integer or infinity.

        EXAMPLES::

            sage: from sage.monoids.monoid import Monoid
            sage: G = Monoid()
            sage: G.order()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def is_finite(self):
        """
        Returns True if this monoid is finite.

        EXAMPLES::

            sage: from sage.monoids.monoid import Monoid
            sage: G = Monoid()
            sage: G.is_finite()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        from sage.rings.infinity import infinity
        return self.order() != infinity


    def is_multiplicative(self):
        """
        Returns True if the monoid operation is given by \* (rather than
        +).

        Override for additive monoids.

        EXAMPLES::

            sage: from sage.monoids.monoid import Monoid
            sage: G = Monoid()
            sage: G.is_multiplicative()
            True
        """
        return True


    def _an_element_(self):
        """
        Return an element

        OUTPUT:

        An element of the monoid.

        EXAMPLES:

            sage: from sage.monoids.monoid import AbelianMonoid
            sage: G = AbelianMonoid([2,3,4,5])
            sage: G.an_element()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def random_element(self, bound=None):
        """
        Return a random element of this monoid.

        EXAMPLES::

            sage: from sage.monoids.monoid import Monoid
            sage: G = Monoid()
            sage: G.random_element()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


cdef class AbelianMonoid(Monoid):
    """
    Generic abelian monoid.
    """
    def is_abelian(self):
        """
        Return True.

        EXAMPLES::

            sage: from sage.monoids.monoid import AbelianMonoid
            sage: G = AbelianMonoid()
            sage: G.is_abelian()
            True
        """
        return True


cdef class FiniteMonoid(Monoid):
    """
    Generic finite monoid.
    """
    def __init__(self, base=None, gens=None, category=None):
        """
        The Python constructor

        TESTS::

            sage: from sage.monoids.monoid import FiniteMonoid
            sage: G = FiniteMonoid()
            sage: G.category()
            Category of finite monoids
        """
        from sage.categories.finite_monoids import FiniteMonoids
        if category is None:
            category = FiniteMonoids()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(FiniteMonoids()) for cat in category):
                raise ValueError("%s is not a subcategory of %s"%(category, FiniteMonoids()))
        Parent.__init__(self, base=base, gens=gens, category=category)


    def is_finite(self):
        """
        Return True.

        EXAMPLES::

            sage: from sage.monoids.monoid import FiniteMonoid
            sage: G = FiniteMonoid()
            sage: G.is_finite()
            True
        """
        return True
