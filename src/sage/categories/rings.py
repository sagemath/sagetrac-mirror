r"""
Rings
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.rngs import Rngs
from sage.structure.element import Element
from functools import reduce


class Rings(CategoryWithAxiom):
    """
    The category of rings

    Associative rings with unit, not necessarily commutative

    EXAMPLES::

        sage: Rings()
        Category of rings
        sage: sorted(Rings().super_categories(), key=str)
        [Category of rngs, Category of semirings]

        sage: sorted(Rings().axioms())
        ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse',
         'AdditiveUnital', 'Associative', 'Distributive', 'Unital']

        sage: Rings() is (CommutativeAdditiveGroups() & Monoids()).Distributive()
        True
        sage: Rings() is Rngs().Unital()
        True
        sage: Rings() is Semirings().AdditiveInverse()
        True

    TESTS::

        sage: TestSuite(Rings()).run()

    .. TODO::

        (see: http://trac.sagemath.org/sage_trac/wiki/CategoriesRoadMap)

        - Make Rings() into a subcategory or alias of Algebras(ZZ);

        - A parent P in the category ``Rings()`` should automatically be
          in the category ``Algebras(P)``.
    """

    _base_category_class_and_axiom = (Rngs, "Unital")

    class MorphismMethods:
        @cached_method
        def is_injective(self):
            """
            Return whether or not this morphism is injective.

            EXAMPLES::

                sage: R.<x,y> = QQ[]
                sage: R.hom([x, y^2], R).is_injective()
                True
                sage: R.hom([x, x^2], R).is_injective()
                False
                sage: S.<u,v> = R.quotient(x^3*y)
                sage: R.hom([v, u], S).is_injective()
                False
                sage: S.hom([-u, v], S).is_injective()
                True
                sage: S.cover().is_injective()
                False

            If the domain is a field, the homomorphism is injective::

                sage: K.<x> = FunctionField(QQ)
                sage: L.<y> = FunctionField(QQ)
                sage: f = K.hom([y]); f
                Function Field morphism:
                  From: Rational function field in x over Rational Field
                  To:   Rational function field in y over Rational Field
                  Defn: x |--> y
                sage: f.is_injective()
                True

            Unless the codomain is the zero ring::

                sage: codomain = Integers(1)
                sage: f = QQ.hom([Zmod(1)(0)], check=False)
                sage: f.is_injective()
                False

            Homomorphism from rings of characteristic zero to rings of positive
            characteristic can not be injective::

                sage: R.<x> = ZZ[]
                sage: f = R.hom([GF(3)(1)]); f
                Ring morphism:
                  From: Univariate Polynomial Ring in x over Integer Ring
                  To:   Finite Field of size 3
                  Defn: x |--> 1
                sage: f.is_injective()
                False

            A morphism whose domain is an order in a number field is injective if
            the codomain has characteristic zero::

                sage: K.<x> = FunctionField(QQ)
                sage: f = ZZ.hom(K); f
                Composite map:
                  From: Integer Ring
                  To:   Rational function field in x over Rational Field
                  Defn:   Conversion via FractionFieldElement_1poly_field map:
                          From: Integer Ring
                          To:   Fraction Field of Univariate Polynomial Ring in x over Rational Field
                        then
                          Isomorphism:
                          From: Fraction Field of Univariate Polynomial Ring in x over Rational Field
                          To:   Rational function field in x over Rational Field
                sage: f.is_injective()
                True

            A coercion to the fraction field is injective::

                sage: R = ZpFM(3)
                sage: R.fraction_field().coerce_map_from(R).is_injective()
                True

            """
            if self.domain().is_zero():
                return True
            if self.codomain().is_zero():
                # the only map to the zero ring that is injective is the map from itself
                return False

            from sage.categories.fields import Fields
            if self.domain() in Fields():
                # A ring homomorphism from a field to a ring is injective
                # (unless the codomain is the zero ring.) Note that ring
                # homomorphism must send the 1 element to the 1 element
                return True

            try:
                ker = self.kernel()
            except (NotImplementedError, AttributeError):
                pass
            else:
                return ker.is_zero()

            if self.domain().characteristic() == 0:
                if self.codomain().characteristic() != 0:
                    return False
                else:
                    from sage.categories.integral_domains import IntegralDomains
                    if self.domain() in IntegralDomains():
                        # if all elements of the domain are algebraic over ZZ,
                        # then the homomorphism must be injective (in
                        # particular if the domain is ZZ)
                        from sage.categories.number_fields import NumberFields
                        if self.domain().fraction_field() in NumberFields():
                            return True

            if self._is_coercion:
                try:
                    K = self.domain().fraction_field()
                except (TypeError, AttributeError, ValueError):
                    pass
                else:
                    if K is self.codomain():
                        return True

            try:
                if self.domain().cardinality() > self.codomain().cardinality():
                    return False
            except AttributeError:
                pass

            raise NotImplementedError

        def _is_nonzero(self):
            r"""
            Return whether this is not the zero morphism.

            .. NOTE::

                We can not override ``is_zero()`` from the category framework
                and we can not implement ``__nonzero__`` because it is a
                special method. That this is why this has a cumbersome name.

            EXAMPLES::

                sage: ZZ.hom(ZZ)._is_nonzero()
                True
                sage: ZZ.hom(Zmod(1))._is_nonzero()
                False

            """
            return bool(self.codomain().one())

        def extend_to_fraction_field(self):
            r"""
            Return the extension of this morphism to fraction fields of
            the domain and the codomain.

            EXAMPLES::

                sage: S.<x> = QQ[]
                sage: f = S.hom([x+1]); f
                Ring endomorphism of Univariate Polynomial Ring in x over Rational Field
                    Defn: x |--> x + 1

                sage: g = f.extend_to_fraction_field(); g
                Ring endomorphism of Fraction Field of Univariate Polynomial Ring in x over Rational Field
                    Defn: x |--> x + 1
                sage: g(x)
                x + 1
                sage: g(1/x)
                1/(x + 1)

            If this morphism is not injective, it does not extend to the fraction
            field and an error is raised::

                sage: f = GF(5).coerce_map_from(ZZ)
                sage: f.extend_to_fraction_field()
                Traceback (most recent call last):
                ...
                ValueError: the morphism is not injective

            TESTS::

                sage: A.<x> = RR[]
                sage: phi = A.hom([x+1])
                sage: phi.extend_to_fraction_field()
                Ring endomorphism of Fraction Field of Univariate Polynomial Ring in x over Real Field with 53 bits of precision
                  Defn: x |--> x + 1.00000000000000
            """
            from sage.rings.morphism import RingHomomorphism_from_fraction_field
            if self.domain().is_field() and self.codomain().is_field():
                return self
            try:
                if not self.is_injective():
                    raise ValueError("the morphism is not injective")
            except (NotImplementedError, TypeError):   # we trust the user
                pass
            domain = self.domain().fraction_field()
            codomain = self.codomain().fraction_field()
            parent = domain.Hom(codomain)   # category = category=self.category_for() ???
            return RingHomomorphism_from_fraction_field(parent, self)


    class SubcategoryMethods:

        def NoZeroDivisors(self):
            r"""
            Return the full subcategory of the objects of ``self`` having
            no nonzero zero divisors.

            A *zero divisor* in a ring `R` is an element `x \in R` such
            that there exists a nonzero element `y \in R` such that
            `x \cdot y = 0` or `y \cdot x = 0`
            (see :wikipedia:`Zero_divisor`).

            EXAMPLES::

                sage: Rings().NoZeroDivisors()
                Category of domains

            .. NOTE::

                This could be generalized to
                :class:`MagmasAndAdditiveMagmas.Distributive.AdditiveUnital`.

            TESTS::

                sage: TestSuite(Rings().NoZeroDivisors()).run()
                sage: Algebras(QQ).NoZeroDivisors.__module__
                'sage.categories.rings'
            """
            return self._with_axiom('NoZeroDivisors')

        def Division(self):
            """
            Return the full subcategory of the division objects of ``self``.

            A ring satisfies the *division axiom* if all non-zero
            elements have multiplicative inverses.

            .. NOTE::

                This could be generalized to
                :class:`MagmasAndAdditiveMagmas.Distributive.AdditiveUnital`.

            EXAMPLES::

                sage: Rings().Division()
                Category of division rings
                sage: Rings().Commutative().Division()
                Category of fields

            TESTS::

                sage: TestSuite(Rings().Division()).run()
                sage: Algebras(QQ).Division.__module__
                'sage.categories.rings'
            """
            return self._with_axiom('Division')


    NoZeroDivisors = LazyImport('sage.categories.domains', 'Domains', at_startup=True)
    Division       = LazyImport('sage.categories.division_rings', 'DivisionRings', at_startup=True)
    Commutative    = LazyImport('sage.categories.commutative_rings', 'CommutativeRings', at_startup=True)

    class ParentMethods:
        def is_ring(self):
            """
            Return True, since this in an object of the category of rings.

            EXAMPLES::

                sage: Parent(QQ,category=Rings()).is_ring()
                True
            """
            return True

        def is_zero(self):
            """
            Return ``True`` if this is the zero ring.

            EXAMPLES::

                sage: Integers(1).is_zero()
                True
                sage: Integers(2).is_zero()
                False
                sage: QQ.is_zero()
                False
                sage: R.<x> = ZZ[]
                sage: R.quo(1).is_zero()
                True
                sage: R.<x> = GF(101)[]
                sage: R.quo(77).is_zero()
                True
                sage: R.quo(x^2+1).is_zero()
                False
            """
            return self.one() == self.zero()

        def bracket(self, x, y):
            """
            Returns the Lie bracket `[x, y] = x y - y x` of `x` and `y`.

            INPUT:

             - ``x``, ``y`` -- elements of ``self``

            EXAMPLES::

                sage: F = AlgebrasWithBasis(QQ).example()
                sage: F
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: a,b,c = F.algebra_generators()
                sage: F.bracket(a,b)
                B[word: ab] - B[word: ba]

            This measures the default of commutation between `x` and `y`.
            `F` endowed with the bracket operation is a Lie algebra;
            in particular, it satisfies Jacobi's identity::

                sage: F.bracket( F.bracket(a,b), c) + F.bracket(F.bracket(b,c),a) + F.bracket(F.bracket(c,a),b)
                0
            """
            return x*y - y*x

        def _Hom_(self, Y, category):
            r"""
            Returns the homset from ``self`` to ``Y`` in the category ``category``

            INPUT:

            - ``Y`` -- a ring
            - ``category`` -- a subcategory of :class:`Rings`() or None

            The sole purpose of this method is to construct the homset
            as a :class:`~sage.rings.homset.RingHomset`. If
            ``category`` is specified and is not a subcategory of
            :class:`Rings`, a ``TypeError`` is raised instead

            This method is not meant to be called directly. Please use
            :func:`sage.categories.homset.Hom` instead.

            EXAMPLES::

                sage: H = QQ._Hom_(QQ, category = Rings()); H
                Set of Homomorphisms from Rational Field to Rational Field
                sage: H.__class__
                <class 'sage.rings.homset.RingHomset_generic_with_category'>

            TESTS::

                sage: Hom(QQ, QQ, category = Rings()).__class__
                <class 'sage.rings.homset.RingHomset_generic_with_category'>

                sage: Hom(CyclotomicField(3), QQ, category = Rings()).__class__
                <class 'sage.rings.number_field.homset.CyclotomicFieldHomset_with_category'>

                sage: TestSuite(Hom(QQ, QQ, category = Rings())).run() # indirect doctest

            """
            if category is not None and not category.is_subcategory(Rings()):
                raise TypeError("%s is not a subcategory of Rings()"%category)
            if Y not in Rings():
                raise TypeError("%s is not a ring"%Y)
            from sage.rings.homset import RingHomset
            return RingHomset(self, Y, category = category)

        # this is already in sage.rings.ring.Ring,
        # but not all rings descend from that class,
        # e.g., matrix spaces.
        def _mul_(self, x, switch_sides=False):
            """
            Multiplication of rings with, e.g., lists.

            NOTE:

            This method is used to create ideals. It is
            the same as the multiplication method for
            :class:`~sage.rings.ring.Ring`. However, not
            all parents that belong to the category of
            rings also inherits from the base class of
            rings. Therefore, we implemented a ``__mul__``
            method for parents, that calls a ``_mul_``
            method implemented here. See :trac:`7797`.

            INPUT:

            - `x`, an object to multiply with.
            - `switch_sides` (optional bool): If ``False``,
              the product is ``self*x``; if ``True``, the
              product is ``x*self``.

            EXAMPLES:

            As we mentioned above, this method is called
            when a ring is involved that does not inherit
            from the base class of rings. This is the case,
            e.g., for matrix algebras::

                sage: MS = MatrixSpace(QQ,2,2)
                sage: isinstance(MS,Ring)
                False
                sage: MS in Rings()
                True
                sage: MS*2     # indirect doctest
                Left Ideal
                (
                  [2 0]
                  [0 2]
                )
                 of Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            In the next example, the ring and the other factor switch sides
            in the product::

                sage: [MS.2]*MS
                Right Ideal
                (
                  [0 0]
                  [1 0]
                )
                 of Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            AUTHOR:

            - Simon King (2011-03-22)

            """
            try:
                if self.is_commutative():
                    return self.ideal(x)
            except (AttributeError,NotImplementedError):
                pass
            try:
                side = x.side()
            except AttributeError:
                return self.ideal(x, side='right' if switch_sides else 'left')
            # presumably x is an ideal...
            try:
                x = x.gens()
            except (AttributeError, NotImplementedError):
                pass # ... not an ideal
            if switch_sides:
                if side in ['right','twosided']:
                    return self.ideal(x,side=side)
                elif side=='left':
                    return self.ideal(x,side='twosided')
            else:
                if side in ['left','twosided']:
                    return self.ideal(x,side=side)
                elif side=='right':
                    return self.ideal(x,side='twosided')
            # duck typing failed
            raise TypeError("Don't know how to transform %s into an ideal of %s"%(x,self))

        def __pow__(self, n):
            """
            Return the free module of rank `n` over this ring.  If n is a tuple of
            two elements, creates a matrix space.

            EXAMPLES::

                sage: QQ^5
                Vector space of dimension 5 over Rational Field
                sage: Integers(20)^1000
                Ambient free module of rank 1000 over Ring of integers modulo 20

                sage: QQ^(2,3)
                Full MatrixSpace of 2 by 3 dense matrices over Rational Field
            """
            if isinstance(n, tuple):
                m, n = n
                from sage.matrix.matrix_space import MatrixSpace
                return MatrixSpace(self, m, n)
            else:
                from sage.modules.free_module import FreeModule
                return FreeModule(self, n)

        @cached_method
        def ideal_monoid(self):
            """
            The monoid of the ideals of this ring.

            NOTE:

            The code is copied from the base class of rings.
            This is since there are rings that do not inherit
            from that class, such as matrix algebras.  See
            :trac:`7797`.

            EXAMPLES::

                sage: MS = MatrixSpace(QQ,2,2)
                sage: isinstance(MS,Ring)
                False
                sage: MS in Rings()
                True
                sage: MS.ideal_monoid()
                Monoid of ideals of Full MatrixSpace of 2 by 2 dense matrices
                over Rational Field

            Note that the monoid is cached::

                sage: MS.ideal_monoid() is MS.ideal_monoid()
                True

            """
            try:
                from sage.rings.ideal_monoid import IdealMonoid
                return IdealMonoid(self)
            except TypeError:
                from sage.rings.noncommutative_ideals import IdealMonoid_nc
                return IdealMonoid_nc(self)

        def characteristic(self):
            """
            Return the characteristic of this ring.

            EXAMPLES::

                sage: QQ.characteristic()
                0
                sage: GF(19).characteristic()
                19
                sage: Integers(8).characteristic()
                8
                sage: Zp(5).characteristic()
                0
            """
            from sage.rings.infinity import infinity
            from sage.rings.integer_ring import ZZ
            order_1 = self.one().additive_order()
            return ZZ.zero() if order_1 is infinity else order_1

        def _test_characteristic(self, **options):
            """
            Run generic tests on the method :meth:`characteristic`.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: ZZ._test_characteristic()
            """
            tester = self._tester(**options)
            try:
                characteristic = self.characteristic()
            except AttributeError:
                return # raised when self.one() does not have a additive_order()
            except NotImplementedError:
                return

            # test that #12988 is fixed
            from sage.rings.integer import Integer
            tester.assertIsInstance(characteristic, Integer)

        def ideal(self, *args, **kwds):
            """
            Create an ideal of this ring.

            NOTE:

            The code is copied from the base class
            :class:`~sage.rings.ring.Ring`. This is
            because there are rings that do not inherit
            from that class, such as matrix algebras.
            See :trac:`7797`.

            INPUT:

            - An element or a list/tuple/sequence of elements.
            - ``coerce`` (optional bool, default ``True``):
              First coerce the elements into this ring.
            - ``side``, optional string, one of ``"twosided"``
              (default), ``"left"``, ``"right"``: determines
              whether the resulting ideal is twosided, a left
              ideal or a right ideal.

            EXAMPLES::

                sage: MS = MatrixSpace(QQ,2,2)
                sage: isinstance(MS,Ring)
                False
                sage: MS in Rings()
                True
                sage: MS.ideal(2)
                Twosided Ideal
                (
                  [2 0]
                  [0 2]
                )
                 of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
                sage: MS.ideal([MS.0,MS.1],side='right')
                Right Ideal
                (
                  [1 0]
                  [0 0],
                <BLANKLINE>
                  [0 1]
                  [0 0]
                )
                 of Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            """
            if 'coerce' in kwds:
                coerce = kwds['coerce']
                del kwds['coerce']
            else:
                coerce = True

            from sage.rings.ideal import Ideal_generic
            from types import GeneratorType
            if len(args) == 0:
                gens = [self(0)]
            else:
                gens = args
                while isinstance(gens, (list, tuple, GeneratorType)) and len(gens) == 1:
                    first = gens[0]
                    if isinstance(first, Ideal_generic):
                        R = first.ring()
                        m = self.convert_map_from(R)
                        if m is not None:
                            gens = [m(g) for g in first.gens()]
                            coerce = False
                        else:
                            m = R.convert_map_from(self)
                            if m is not None:
                                raise NotImplementedError
                            else:
                                raise TypeError
                        break
                    elif isinstance(first, (list, tuple, GeneratorType)):
                        gens = first
                    else:
                        try:
                            if self.has_coerce_map_from(first):
                                gens = first.gens() # we have a ring as argument
                            elif isinstance(first, Element):
                                gens = [first]
                            else:
                                raise ArithmeticError("There is no coercion from %s to %s"%(first,self))
                        except TypeError: # first may be a ring element
                            pass
                        break
            if coerce:
                gens = [self(g) for g in gens]
            from sage.categories.principal_ideal_domains import PrincipalIdealDomains
            if self in PrincipalIdealDomains():
                # Use GCD algorithm to obtain a principal ideal
                g = gens[0]
                if len(gens) == 1:
                    try:
                        g = g.gcd(g) # note: we set g = gcd(g, g) to "canonicalize" the generator: make polynomials monic, etc.
                    except (AttributeError, NotImplementedError):
                        pass
                else:
                    for h in gens[1:]:
                        g = g.gcd(h)
                gens = [g]
            if 'ideal_class' in kwds:
                C = kwds['ideal_class']
                del kwds['ideal_class']
            else:
                C = self._ideal_class_(len(gens))
            if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
                gens = gens[0]
            return C(self, gens, **kwds)

        def _ideal_class_(self,n=0):
            """
            Return the class that is used to implement ideals of this ring.

            NOTE:

            We copy the code from :class:`~sage.rings.ring.Ring`. This is
            necessary because not all rings inherit from that class, such
            as matrix algebras.

            INPUT:

            - ``n`` (optional integer, default 0): The number of generators
              of the ideal to be created.

            OUTPUT:

            The class that is used to implement ideals of this ring with
            ``n`` generators.

            NOTE:

            Often principal ideals (``n==1``) are implemented via a different
            class.

            EXAMPLES::

                sage: MS = MatrixSpace(QQ,2,2)
                sage: MS._ideal_class_()
                <class 'sage.rings.noncommutative_ideals.Ideal_nc'>

            We don't know of a commutative ring in Sage that does not inherit
            from the base class of rings. So, we need to cheat in the next
            example::

                sage: super(Ring,QQ)._ideal_class_.__module__
                'sage.categories.rings'
                sage: super(Ring,QQ)._ideal_class_()
                <class 'sage.rings.ideal.Ideal_generic'>
                sage: super(Ring,QQ)._ideal_class_(1)
                <class 'sage.rings.ideal.Ideal_principal'>
                sage: super(Ring,QQ)._ideal_class_(2)
                <class 'sage.rings.ideal.Ideal_generic'>

            """
            from sage.rings.noncommutative_ideals import Ideal_nc
            try:
                if not self.is_commutative():
                    return Ideal_nc
            except (NotImplementedError,AttributeError):
                return Ideal_nc
            from sage.rings.ideal import Ideal_generic, Ideal_principal
            if n == 1:
                return Ideal_principal
            else:
                return Ideal_generic

        ##
        # Quotient rings
        # Again, this is defined in sage.rings.ring.pyx
        def quotient(self, I, names=None, **kwds):
            """
            Quotient of a ring by a two-sided ideal.

            INPUT:

            - ``I``: A twosided ideal of this ring.
            - ``names``: a list of strings to be used as names
              for the variables in the quotient ring.
            - further named arguments that may be passed to the
              quotient ring constructor.

            EXAMPLES:

            Usually, a ring inherits a method :meth:`sage.rings.ring.Ring.quotient`.
            So, we need a bit of effort to make the following example work with the
            category framework::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: from sage.rings.noncommutative_ideals import Ideal_nc
                sage: from itertools import product
                sage: class PowerIdeal(Ideal_nc):
                ....:  def __init__(self, R, n):
                ....:      self._power = n
                ....:      Ideal_nc.__init__(self, R, [R.prod(m) for m in product(R.gens(), repeat=n)])
                ....:  def reduce(self, x):
                ....:      R = self.ring()
                ....:      return add([c*R(m) for m,c in x if len(m) < self._power], R(0))
                sage: I = PowerIdeal(F,3)
                sage: Q = Rings().parent_class.quotient(F, I); Q
                Quotient of Free Algebra on 3 generators (x, y, z) over Rational Field by the ideal (x^3, x^2*y, x^2*z, x*y*x, x*y^2, x*y*z, x*z*x, x*z*y, x*z^2, y*x^2, y*x*y, y*x*z, y^2*x, y^3, y^2*z, y*z*x, y*z*y, y*z^2, z*x^2, z*x*y, z*x*z, z*y*x, z*y^2, z*y*z, z^2*x, z^2*y, z^3)
                sage: Q.0
                xbar
                sage: Q.1
                ybar
                sage: Q.2
                zbar
                sage: Q.0*Q.1
                xbar*ybar
                sage: Q.0*Q.1*Q.0
                0
            """
            from sage.rings.quotient_ring import QuotientRing
            return QuotientRing(self, I, names=names, **kwds)

        def quo(self, I, names=None, **kwds):
            """
            Quotient of a ring by a two-sided ideal.

            NOTE:

            This is a synonym for :meth:`quotient`.

            EXAMPLES::

                sage: MS = MatrixSpace(QQ,2)
                sage: I = MS*MS.gens()*MS

            ``MS`` is not an instance of :class:`~sage.rings.ring.Ring`.

            However it is an instance of the parent class of the
            category of rings. The quotient method is inherited from
            there::

                sage: isinstance(MS,sage.rings.ring.Ring)
                False
                sage: isinstance(MS,Rings().parent_class)
                True
                sage: MS.quo(I,names = ['a','b','c','d'])
                Quotient of Full MatrixSpace of 2 by 2 dense matrices over Rational Field by the ideal
                (
                  [1 0]
                  [0 0],
                <BLANKLINE>
                  [0 1]
                  [0 0],
                <BLANKLINE>
                  [0 0]
                  [1 0],
                <BLANKLINE>
                  [0 0]
                  [0 1]
                )

            """
            return self.quotient(I,names=names,**kwds)

        def quotient_ring(self, I, names=None, **kwds):
            """
            Quotient of a ring by a two-sided ideal.

            NOTE:

            This is a synonyme for :meth:`quotient`.

            EXAMPLES::

                sage: MS = MatrixSpace(QQ,2)
                sage: I = MS*MS.gens()*MS

            ``MS`` is not an instance of :class:`~sage.rings.ring.Ring`,
            but it is an instance of the parent class of the category of
            rings. The quotient method is inherited from there::

                sage: isinstance(MS,sage.rings.ring.Ring)
                False
                sage: isinstance(MS,Rings().parent_class)
                True
                sage: MS.quotient_ring(I,names = ['a','b','c','d'])
                Quotient of Full MatrixSpace of 2 by 2 dense matrices over Rational Field by the ideal
                (
                  [1 0]
                  [0 0],
                <BLANKLINE>
                  [0 1]
                  [0 0],
                <BLANKLINE>
                  [0 0]
                  [1 0],
                <BLANKLINE>
                  [0 0]
                  [0 1]
                )

            """
            return self.quotient(I,names=names, **kwds)

        def __truediv__(self, I):
            """
            Since assigning generator names would not work properly,
            the construction of a quotient ring using division syntax
            is not supported.

            EXAMPLES::

                sage: MS = MatrixSpace(QQ,2)
                sage: I = MS*MS.gens()*MS
                sage: MS/I
                Traceback (most recent call last):
                ...
                TypeError: Use self.quo(I) or self.quotient(I) to construct the quotient ring.
            """
            raise TypeError("Use self.quo(I) or self.quotient(I) to construct the quotient ring.")

        def __getitem__(self, arg):
            """
            Extend this ring by one or several elements to create a polynomial
            ring, a power series ring, or an algebraic extension.

            This is a convenience method intended primarily for interactive
            use.

            .. SEEALSO::

                :func:`~sage.rings.polynomial.polynomial_ring_constructor.PolynomialRing`,
                :func:`~sage.rings.power_series_ring.PowerSeriesRing`,
                :meth:`~sage.rings.ring.Ring.extension`,
                :meth:`sage.rings.integer_ring.IntegerRing_class.__getitem__`,
                :meth:`sage.rings.matrix_space.MatrixSpace.__getitem__`,
                :meth:`sage.structure.parent.Parent.__getitem__`

            EXAMPLES:

            We create several polynomial rings::

                sage: ZZ['x']
                Univariate Polynomial Ring in x over Integer Ring
                sage: QQ['x']
                Univariate Polynomial Ring in x over Rational Field
                sage: GF(17)['abc']
                Univariate Polynomial Ring in abc over Finite Field of size 17
                sage: GF(17)['a,b,c']
                Multivariate Polynomial Ring in a, b, c over Finite Field of size 17
                sage: GF(17)['a']['b']
                Univariate Polynomial Ring in b over Univariate Polynomial Ring in a over Finite Field of size 17

            We can create Ore polynomial rings::

                sage: k.<t> = GF(5^3)
                sage: Frob = k.frobenius_endomorphism()
                sage: k['x', Frob]
                Ore Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5

                sage: R.<t> = QQ[]
                sage: der = R.derivation()
                sage: R['d', der]
                Ore Polynomial Ring in d over Univariate Polynomial Ring in t over Rational Field twisted by d/dt

            We can also create power series rings by using double brackets::

                sage: QQ[['t']]
                Power Series Ring in t over Rational Field
                sage: ZZ[['W']]
                Power Series Ring in W over Integer Ring

                sage: ZZ[['x,y,z']]
                Multivariate Power Series Ring in x, y, z over Integer Ring
                sage: ZZ[['x','T']]
                Multivariate Power Series Ring in x, T over Integer Ring

            Use :func:`~sage.rings.fraction_field.Frac` or
            :meth:`~sage.rings.ring.CommutativeRing.fraction_field` to obtain
            the fields of rational functions and Laurent series::

                sage: Frac(QQ['t'])
                Fraction Field of Univariate Polynomial Ring in t over Rational Field
                sage: Frac(QQ[['t']])
                Laurent Series Ring in t over Rational Field
                sage: QQ[['t']].fraction_field()
                Laurent Series Ring in t over Rational Field

            Note that the same syntax can be used to create number fields::

                sage: QQ[I]
                Number Field in I with defining polynomial x^2 + 1 with I = 1*I
                sage: QQ[I].coerce_embedding()
                Generic morphism:
                  From: Number Field in I with defining polynomial x^2 + 1 with I = 1*I
                  To:   Complex Lazy Field
                  Defn: I -> 1*I

            ::

                sage: QQ[sqrt(2)]
                Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?
                sage: QQ[sqrt(2)].coerce_embedding()
                Generic morphism:
                  From: Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?
                  To:   Real Lazy Field
                  Defn: sqrt2 -> 1.414213562373095?

            ::

                sage: QQ[sqrt(2),sqrt(3)]
                Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field

            and orders in number fields::

                sage: ZZ[I]
                Order in Number Field in I0 with defining polynomial x^2 + 1 with I0 = 1*I
                sage: ZZ[sqrt(5)]
                Order in Number Field in sqrt5 with defining polynomial x^2 - 5 with sqrt5 = 2.236067977499790?
                sage: ZZ[sqrt(2)+sqrt(3)]
                Order in Number Field in a with defining polynomial x^4 - 10*x^2 + 1 with a = 3.146264369941973?

            Embeddings are found for simple extensions (when that makes sense)::

                sage: QQi.<i> = QuadraticField(-1, 'i')
                sage: QQ[i].coerce_embedding()
                Generic morphism:
                  From: Number Field in i with defining polynomial x^2 + 1 with i = 1*I
                  To:   Complex Lazy Field
                  Defn: i -> 1*I

            TESTS:

            A few corner cases::

                sage: QQ[()]
                Multivariate Polynomial Ring in no variables over Rational Field

                sage: QQ[[]]
                Traceback (most recent call last):
                ...
                TypeError: power series rings must have at least one variable

            These kind of expressions do not work::

                sage: QQ['a,b','c']
                Traceback (most recent call last):
                ...
                ValueError: variable name 'a,b' is not alphanumeric
                sage: QQ[['a,b','c']]
                Traceback (most recent call last):
                ...
                ValueError: variable name 'a,b' is not alphanumeric

                sage: QQ[[['x']]]
                Traceback (most recent call last):
                ...
                TypeError: expected R[...] or R[[...]], not R[[[...]]]

            Extension towers are built as follows and use distinct generator names::

                sage: K = QQ[2^(1/3), 2^(1/2), 3^(1/3)]
                sage: K
                Number Field in a with defining polynomial x^3 - 2 over its base field
                sage: K.base_field()
                Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field
                sage: K.base_field().base_field()
                Number Field in b with defining polynomial x^3 - 3

            Embeddings::

                sage: a = 10^100; expr = (2*a + sqrt(2))/(2*a^2-1)
                sage: QQ[expr].coerce_embedding() is None
                False
                sage: QQ[sqrt(5)].gen() > 0
                True
                sage: expr = sqrt(2) + I*(cos(pi/4, hold=True) - sqrt(2)/2)
                sage: QQ[expr].coerce_embedding()
                Generic morphism:
                  From: Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?
                  To:   Real Lazy Field
                  Defn: a -> 1.414213562373095?
            """
            def normalize_arg(arg):
                if isinstance(arg, (tuple, list)):
                    # Allowing arbitrary iterables would create confusion, but we
                    # may want to support a few more.
                    return tuple(arg)
                elif isinstance(arg, str):
                    return tuple(arg.split(','))
                else:
                    return (arg,)

            # 1. If arg is a list, try to return a power series ring.

            if isinstance(arg, list):
                if not arg:
                    raise TypeError("power series rings must have at least one variable")
                elif len(arg) == 1:
                    # R[["a,b"]], R[[(a,b)]]...
                    if isinstance(arg[0], list):
                        raise TypeError("expected R[...] or R[[...]], not R[[[...]]]")
                    elts = normalize_arg(arg[0])
                else:
                    elts = normalize_arg(arg)
                from sage.rings.power_series_ring import PowerSeriesRing
                return PowerSeriesRing(self, elts)

            if isinstance(arg, tuple):
                from sage.categories.morphism import Morphism
                from sage.rings.derivation import RingDerivation
                if len(arg) == 2 and isinstance(arg[1], (Morphism, RingDerivation)):
                    from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
                    return OrePolynomialRing(self, arg[1], names=arg[0])

            # 2. Otherwise, if all specified elements are algebraic, try to
            #    return an algebraic extension

            elts = normalize_arg(arg)

            try:
                minpolys = [a.minpoly() for a in elts]
            except (AttributeError, NotImplementedError, ValueError, TypeError):
                minpolys = None

            if minpolys:
                # how to pass in names?
                names = tuple(_gen_names(elts))
                if len(elts) == 1:
                    from sage.rings.cif import CIF
                    elt = elts[0]
                    try:
                        iv = CIF(elt)
                    except (TypeError, ValueError):
                        emb = None
                    else:
                        # First try creating an ANRoot manually, because
                        # extension(..., embedding=CLF(expr)) (or
                        # ...QQbar(expr)) would normalize the expression in
                        # QQbar, which currently is VERY slow in many cases.
                        # This may fail when minpoly has close roots or elt is
                        # a complicated symbolic expression.
                        # TODO: Rewrite using #19362 and/or #17886 and/or
                        # #15600 once those issues are solved.
                        from sage.rings.qqbar import AlgebraicNumber, ANRoot
                        try:
                            elt = AlgebraicNumber(ANRoot(minpolys[0], iv))
                        except ValueError:
                            pass
                        # Force a real embedding when possible, to get the
                        # right ordered ring structure.
                        from sage.rings.real_lazy import CLF, RLF
                        if (iv.imag().is_zero() or iv.imag().contains_zero()
                                                   and elt.imag().is_zero()):
                            emb = RLF(elt)
                        else:
                            emb = CLF(elt)
                        return self.extension(minpolys[0], names[0], embedding=emb)
                try:
                    # Doing the extension all at once is best, if possible...
                    return self.extension(minpolys, names)
                except (TypeError, ValueError):
                    # ...but we can also construct it iteratively
                    return reduce(lambda R, ext: R.extension(*ext), zip(minpolys, names), self)

            # 2. Otherwise, try to return a polynomial ring

            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            return PolynomialRing(self, elts)

        def free_module(self, base=None, basis=None, map=True):
            """
            Return a free module `V` over the specified subring together with maps to and from `V`.

            The default implementation only supports the case that the base ring is the ring itself.

            INPUT:

            - ``base`` -- a subring `R` so that this ring is isomorphic
              to a finite-rank free `R`-module `V`

            - ``basis`` -- (optional) a basis for this ring over the base

            - ``map`` -- boolean (default ``True``), whether to return
              `R`-linear maps to and from `V`

            OUTPUT:

            - A finite-rank free `R`-module `V`

            - An `R`-module isomorphism from `V` to this ring
              (only included if ``map`` is ``True``)

            - An `R`-module isomorphism from this ring to `V`
              (only included if ``map`` is ``True``)

            EXAMPLES::

                sage: R.<x> = QQ[[]]
                sage: V, from_V, to_V = R.free_module(R)
                sage: v = to_V(1+x); v
                (1 + x)
                sage: from_V(v)
                1 + x
                sage: W, from_W, to_W = R.free_module(R, basis=(1-x))
                sage: W is V
                True
                sage: w = to_W(1+x); w
                (1 - x^2)
                sage: from_W(w)
                1 + x + O(x^20)
            """
            if base is None:
                base = self.base_ring()
            if base is self:
                V = self**1
                if not map:
                    return V
                if basis is not None:
                    if isinstance(basis, (list, tuple)):
                        if len(basis) != 1:
                            raise ValueError("Basis must have length 1")
                        basis = basis[0]
                    basis = self(basis)
                    if not basis.is_unit():
                        raise ValueError("Basis element must be a unit")
                from sage.modules.free_module_morphism import BaseIsomorphism1D_from_FM, BaseIsomorphism1D_to_FM
                Hfrom = V.Hom(self)
                Hto = self.Hom(V)
                from_V = Hfrom.__make_element_class__(BaseIsomorphism1D_from_FM)(Hfrom, basis=basis)
                to_V = Hto.__make_element_class__(BaseIsomorphism1D_to_FM)(Hto, basis=basis)
                return V, from_V, to_V
            else:
                if not self.has_coerce_map_from(base):
                    raise ValueError("base must be a subring of this ring")
                raise NotImplementedError

    class ElementMethods:
        def is_unit(self):
            r"""
            Return whether this element is a unit in the ring.

            .. NOTE::

                This is a generic implementation for (non-commutative) rings
                which only works for the one element, its additive inverse, and
                the zero element.  Most rings should provide a more specialized
                implementation.

            EXAMPLES::

                sage: MS = MatrixSpace(ZZ, 2)
                sage: MS.one().is_unit()
                True
                sage: MS.zero().is_unit()
                False
                sage: MS([1,2,3,4]).is_unit()
                False
            """
            if self.is_one() or (-self).is_one():
                return True
            if self.is_zero(): # now 0 != 1
                return False
            raise NotImplementedError

        def inverse_of_unit(self):
            r"""
            Return the inverse of this element if it is a unit.

            OUTPUT:

            An element in the same ring as this element.

            EXAMPLES::

                sage: R.<x> = ZZ[]
                sage: S = R.quo(x^2 + x + 1)
                sage: S(1).inverse_of_unit()
                1

            This method fails when the element is not a unit::

                sage: 2.inverse_of_unit()
                Traceback (most recent call last):
                ...
                ArithmeticError: inverse does not exist

            The inverse returned is in the same ring as this element::

                sage: a = -1
                sage: a.parent()
                Integer Ring
                sage: a.inverse_of_unit().parent()
                Integer Ring

            Note that this is often not the case when computing inverses in other ways::

                sage: (~a).parent()
                Rational Field
                sage: (1/a).parent()
                Rational Field

            """
            try:
                if not self.is_unit():
                    raise ArithmeticError("element is not a unit")
            except NotImplementedError:
                # if an element does not implement is_unit, we just try to
                # invert it anyway; if the result is in the ring again, it was
                # a unit
                pass

            inverse = ~self
            if inverse not in self.parent():
                raise ArithmeticError("element is not a unit")

            # return the inverse (with the correct parent)
            return self.parent()(inverse)

        def _divide_if_possible(self, y):
            """
            Divide ``self`` by ``y`` if possible and raise a
            ``ValueError`` otherwise.

            EXAMPLES::

                sage: 4._divide_if_possible(2)
                2
                sage: _.parent()
                Integer Ring

            ::

                sage: 4._divide_if_possible(3)
                Traceback (most recent call last):
                ...
                ValueError: 4 is not divisible by 3
            """
            q, r = self.quo_rem(y)
            if r != 0:
                raise ValueError("%s is not divisible by %s"%(self, y))
            return q

def _gen_names(elts):
    r"""
    Used to find a name for a generator when rings are created using the
    ``__getitem__`` syntax, e.g. ``ZZ['x']``, ``ZZ[sqrt(2)]``.

    EXAMPLES::

        sage: from sage.categories.rings import _gen_names
        sage: list(_gen_names([sqrt(5)]))
        ['sqrt5']
        sage: list(_gen_names([sqrt(-17), 2^(1/3)]))
        ['a', 'b']
        sage: list(_gen_names((1..27)))[-1]
        'aa'
    """
    import re
    from sage.structure.category_object import certify_names
    from sage.combinat.words.words import Words
    it = iter(Words("abcdefghijklmnopqrstuvwxyz", infinite=False))
    next(it) # skip empty word
    for x in elts:
        name = str(x)
        m = re.match(r'^sqrt\((\d+)\)$', name)
        if m:
            name = "sqrt%s" % m.groups()[0]
        try:
            certify_names([name])
        except ValueError:
            name = next(it).string_rep()
        yield name
