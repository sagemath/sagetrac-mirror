r"""
Algebras With Basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from functools import partial
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.constant_function import ConstantFunction
from sage.misc.misc import compose
from sage.categories.all import ModulesWithBasis, Algebras
from sage.categories.tensor import TensorProductsCategory, tensor
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from unital_algebras import UnitalAlgebras

class AlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of algebras with a distinguished basis

    EXAMPLES::

        sage: C = AlgebrasWithBasis(QQ); C
        Category of algebras with basis over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of algebras over Rational Field,
         Category of unital algebras with basis over Rational Field]

    We construct a typical parent in this category, and do some computations with it::

        sage: A = C.example(); A
        An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field

        sage: A.category()
        Category of algebras with basis over Rational Field

        sage: A.one_basis()
        word:
        sage: A.one()
        B[word: ]

        sage: A.base_ring()
        Rational Field
        sage: A.basis().keys()
        Words over {'a', 'b', 'c'}

        sage: (a,b,c) = A.algebra_generators()
        sage: a^3, b^2
        (B[word: aaa], B[word: bb])
        sage: a*c*b
        B[word: acb]

        sage: A.product
        <bound method FreeAlgebra_with_category._product_from_product_on_basis_multiply of An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field>
        sage: A.product(a*b,b)
        B[word: abb]

        sage: TestSuite(A).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_category() . . . pass
        running ._test_characteristic() . . . pass
        running ._test_distributivity() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_nonzero_equal() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_zero() . . . pass
        sage: A.__class__
        <class 'sage.categories.examples.algebras_with_basis.FreeAlgebra_with_category'>
        sage: A.element_class
        <class 'sage.combinat.free_module.FreeAlgebra_with_category.element_class'>

    Please see the source code of `A` (with ``A??``) for how to
    implement other algebras with basis.

    TESTS::

        sage: TestSuite(AlgebrasWithBasis(QQ)).run()
    """

    def example(self, alphabet = ('a','b','c')):
        """
        Returns an example of algebra with basis::

            sage: AlgebrasWithBasis(QQ).example()
            An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field

        An other set of generators can be specified as optional argument::

            sage: AlgebrasWithBasis(QQ).example((1,2,3))
            An example of an algebra with basis: the free algebra on the generators (1, 2, 3) over Rational Field
        """
        from sage.categories.examples.algebras_with_basis import Example
        return Example(self.base_ring(), alphabet)

    FiniteDimensional = LazyImport('sage.categories.finite_dimensional_algebras_with_basis', 'FiniteDimensionalAlgebrasWithBasis')
    Graded = LazyImport('sage.categories.graded_algebras_with_basis', 'GradedAlgebrasWithBasis')

    class ParentMethods:

        # For backward compatibility
        one = UnitalAlgebras.WithBasis.ParentMethods.one

        # Backward compatibility temporary cruft to help migrating form CombinatorialAlgebra
        def _product_from_combinatorial_algebra_multiply(self,left,right):
            """
            Returns left\*right where left and right are elements of self.
            product() uses either _multiply or _multiply basis to carry out
            the actual multiplication.

            EXAMPLES::

                sage: s = SymmetricFunctions(QQ).schur()
                sage: a = s([2])
                sage: s._product_from_combinatorial_algebra_multiply(a,a)
                s[2, 2] + s[3, 1] + s[4]
                sage: s.product(a,a)
                s[2, 2] + s[3, 1] + s[4]
            """
            A = left.parent()
            BR = A.base_ring()
            z_elt = {}

            #Do the case where the user specifies how to multiply basis elements
            if hasattr(self, '_multiply_basis'):
                for (left_m, left_c) in left._monomial_coefficients.iteritems():
                    for (right_m, right_c) in right._monomial_coefficients.iteritems():
                        res = self._multiply_basis(left_m, right_m)
                        #Handle the case where the user returns a dictionary
                        #where the keys are the monomials and the values are
                        #the coefficients.  If res is not a dictionary, then
                        #it is assumed to be an element of self
                        if not isinstance(res, dict):
                            if isinstance(res, self._element_class):
                                res = res._monomial_coefficients
                            else:
                                res = {res: BR(1)}
                        for m in res:
                            if m in z_elt:
                                z_elt[ m ] = z_elt[m] + left_c * right_c * res[m]
                            else:
                                z_elt[ m ] = left_c * right_c * res[m]

            #We assume that the user handles the multiplication correctly on
            #his or her own, and returns a dict with monomials as keys and
            #coefficients as values
            else:
                m = self._multiply(left, right)
                if isinstance(m, self._element_class):
                    return m
                if not isinstance(m, dict):
                    z_elt = m.monomial_coefficients()
                else:
                    z_elt = m

            #Remove all entries that are equal to 0
            BR = self.base_ring()
            zero = BR(0)
            del_list = []
            for m, c in z_elt.iteritems():
                if c == zero:
                    del_list.append(m)
            for m in del_list:
                del z_elt[m]

            return self._from_dict(z_elt)

        #def _test_product(self, **options):
        #    tester = self._tester(**options)
        #    tester.assert_(self.product is not None)
        #    could check that self.product is in Hom( self x self, self)

        def _product_on_basis_tuple(self, tup):
            return self.product_on_basis(*tup)

        def _product_on_tensor_basis(self):
            r"""
            Return the function whose sole input is an index for the tensor square of ``self``,
            and whose output is the product of monomials indexed by the corresponding pair of basis elements.

            EXAMPLES::

                sage: W=WeylGroup("A2",prefix="s")
                sage: A = W.algebra(ZZ)
                sage: AA = tensor([A,A])
                sage: k = AA.basis().keys().an_element(); k
                (s1*s2*s1, s1*s2*s1)
                sage: A._product_on_tensor_basis()(k)
                B[1]
                sage: AAAA = tensor([AA,AA])
                sage: kk = AAAA.basis().keys().an_element(); kk
                (s1*s2*s1, s1*s2*s1, s1*s2*s1, s1*s2*s1)
                sage: AA._product_on_tensor_basis()(kk)
                B[1] # B[1]

            """
            # The algebra tensor product has functions which handle basis indices
            AA = tensor([self,self])
            return compose(self._product_on_basis_tuple, AA.index_to_indices())

        def _product_morphism(self):
            r"""
            The multiplication map as a module morphism from the twofold tensor of ``self``, to ``self``.
            
            EXAMPLES::

                sage: W = WeylGroup(CartanType(['A',2]),prefix="s")
                sage: A = W.algebra(ZZ); A.rename("A")
                sage: mA = A._product_morphism(); mA
                Generic morphism:
                From: A # A
                To:   A
                sage: a = A.monomial(W.from_reduced_word([1,2])); a
                B[s1*s2]
                sage: b = A.an_element(); b
                B[s1*s2*s1] + 3*B[s1*s2] + 3*B[s2*s1]
                sage: ab = tensor([a,b]); ab
                B[s1*s2] # B[s1*s2*s1] + 3*B[s1*s2] # B[s1*s2] + 3*B[s1*s2] # B[s2*s1]
                sage: mA(mA.domain()(ab))
                3*B[s2*s1] + B[s2] + 3*B[1]
                sage: a*b
                3*B[s2*s1] + B[s2] + 3*B[1]
                sage: AA = tensor([A,A])
                sage: mAA = AA._product_morphism()
                sage: mAA(mAA.domain()(tensor([a,b,b,b])))
                27*B[s2*s1] # B[s1*s2] + 27*B[s2*s1] # B[s2*s1] + 18*B[s2*s1] # B[s1] + 18*B[s2*s1] # B[s2] + 57*B[s2*s1] # B[1] + 9*B[s2] # B[s1*s2] + 9*B[s2] # B[s2*s1] + 6*B[s2] # B[s1] + 6*B[s2] # B[s2] + 19*B[s2] # B[1] + 27*B[1] # B[s1*s2] + 27*B[1] # B[s2*s1] + 18*B[1] # B[s1] + 18*B[1] # B[s2] + 57*B[1] # B[1]

            """
            module_category = ModulesWithBasis(self.base_ring())
            AA = tensor([self, self], category=module_category)
            return AA.module_morphism(on_basis = self._product_on_tensor_basis(), codomain=self, category=module_category)

        def _unit_morphism(self):
            r"""
            The unit map as a module morphism from the base ring to ``self``.

            Note that the domain is the tensor unit module,
            the base ring viewed as a free module.

            EXAMPLES::

                sage: W = WeylGroup(CartanType(['A',2]),prefix="s")
                sage: A = W.algebra(ZZ); A.rename("A")
                sage: unit = A._unit_morphism(); unit
                Generic morphism:
                From: The unit object in Category of tensor products of hopf algebras with basis over Integer Ring
                To:   A
                sage: unit(3*unit.domain().one())
                3*B[1]

            """
            U = self.tensor_unit()
            one_func = ConstantFunction(self.one())
            return U.module_morphism(on_basis = one_func, codomain=self, category=ModulesWithBasis(self.base_ring()))

    class ElementMethods:

        def __invert__(self):
            """
            Returns the inverse of self if self is a multiple of one,
            and one is in the basis of this algebra. Otherwise throws
            an error.

            Caveat: this generic implementation is not complete; there
            may be invertible elements in the algebra that can't be
            inversed this way. It is correct though for graded
            connected algebras with basis.

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ).example()
                sage: x = C(2); x
                2*B[word: ]
                sage: ~x
                1/2*B[word: ]
                sage: a = C.algebra_generators().first(); a
                B[word: a]
                sage: ~a
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= B[word: a])
            """
            # FIXME: make this generic
            mcs = self._monomial_coefficients
            one = self.parent().one_basis()
            if len(mcs) == 1 and one in mcs:
                return self.parent()( ~mcs[ one ] )
            else:
                raise ValueError, "cannot invert self (= %s)"%self


    class CartesianProducts(CartesianProductsCategory):
        """
        The category of algebras with basis, constructed as cartesian products of algebras with basis

        Note: this construction give the direct products of algebras with basis.
        See comment in :class:`Algebras.CartesianProducts
        <sage.categories.algebras.Algebras.CartesianProducts>`
        """

        def extra_super_categories(self):
            """
            A cartesian product of algebras with basis is endowed with
            a natural algebra with basis structure.

            EXAMPLES::

                sage: AlgebrasWithBasis(QQ).CartesianProducts().extra_super_categories()
                [Category of algebras with basis over Rational Field]
                sage: AlgebrasWithBasis(QQ).CartesianProducts().super_categories()
                [Category of algebras with basis over Rational Field, Category of Cartesian products of algebras over Rational Field, Category of Cartesian products of modules with basis over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            @cached_method   # todo: reinstate once #5843 is fixed
            def one_from_cartesian_product_of_one_basis(self):
                """
                Returns the one of this cartesian product of algebras, as per ``Monoids.ParentMethods.one``

                It is constructed as the cartesian product of the ones of the
                summands, using their :meth:`.one_basis` methods.

                This implementation does not require multiplication by
                scalars nor calling cartesian_product. This might help keeping
                things as lazy as possible upon initialization.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: A.one_basis()
                    word:

                    sage: B = cartesian_product((A, A, A))
                    sage: B.one_from_cartesian_product_of_one_basis()
                    B[(0, word: )] + B[(1, word: )] + B[(2, word: )]
                    sage: B.one()
                    B[(0, word: )] + B[(1, word: )] + B[(2, word: )]

                    sage: cartesian_product([SymmetricGroupAlgebra(QQ, 3), SymmetricGroupAlgebra(QQ, 4)]).one()
                    B[(0, [1, 2, 3])] + B[(1, [1, 2, 3, 4])]
                """
                return self.sum_of_monomials( zip( self._sets_keys(), (set.one_basis() for set in self._sets)) )

            @lazy_attribute
            def one(self):
                """
                TESTS::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: B = cartesian_product((A, A, A))
                    sage: B.one()
                    B[(0, word: )] + B[(1, word: )] + B[(2, word: )]
                """
                if all(hasattr(module, "one_basis") for module in self._sets):
                    return self.one_from_cartesian_product_of_one_basis
                else:
                    return NotImplemented

            #def product_on_basis(self, t1, t2):
            # would be easy to implement, but without a special
            # version of module morphism, this would not take
            # advantage of the bloc structure


    class TensorProducts(TensorProductsCategory):
        """
        The category of algebras with basis constructed by tensor product of algebras with basis
        """

        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: AlgebrasWithBasis(QQ).TensorProducts().extra_super_categories()
                [Category of algebras with basis over Rational Field]
                sage: AlgebrasWithBasis(QQ).TensorProducts().super_categories()
                [Category of algebras with basis over Rational Field, Category of tensor products of algebras over Rational Field, Category of tensor products of modules with basis over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            """
            Implements operations on tensor products of algebras with basis
            """

            @abstract_method
            def factors(self):
                r"""
                The tensor factor algebras of `self`.

                EXAMPLES::

                    sage: W = WeylGroup("A2",prefix="s")
                    sage: A = W.algebra(ZZ); A.rename("A")
                    sage: A2 = tensor([A,A])
                    sage: A22 = tensor([A2,A2])
                    sage: A22.factors()
                    (A # A, A # A)
                    sage: A4 = tensor([A,A,A,A])
                    sage: A4.factors()
                    (A, A, A, A)

                """
                pass

            @abstract_method
            def index_to_indices(self):
                r"""
                The function which maps a basis index of ``self`` to a tuple of basis indices for ``self.factors()``.
                """
                pass

            @abstract_method
            def indices_to_index(self):
                r"""
                The function which maps a tuple of basis indices for ``self.factors()`` to an index of ``self.basis()``.
                """
                pass

            @cached_method
            def one_basis(self):
                """
                Returns the index of the one of this tensor product of
                algebras, as per ``AlgebrasWithBasis.ParentMethods.one_basis``

                It is the tuple whose operands are the indices of the
                ones of the operands, as returned by their
                :meth:`.one_basis` methods.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: A.one_basis()
                    word:
                    sage: B = tensor((A, A, A))
                    sage: B.one_basis()
                    (word: , word: , word: )
                    sage: B.one()
                    B[word: ] # B[word: ] # B[word: ]


                """
                # FIXME: this method should be conditionally defined,
                # so that .one_basis returns NotImplemented if not
                # all modules provide one_basis
                if all(hasattr(module, "one_basis") for module in self.factors()):
                    return self.indices_to_index()(*tuple(module.one_basis() for module in self.factors()))
                else:
                    raise NotImplementedError

            def product_on_basis(self, t1, t2):
                """
                The product of the algebra on the basis, as per
                ``AlgebrasWithBasis.ParentMethods.product_on_basis``.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: (a,b,c) = A.algebra_generators()

                    sage: x = tensor( (a, b, c) ); x
                    B[word: a] # B[word: b] # B[word: c]
                    sage: y = tensor( (c, b, a) ); y
                    B[word: c] # B[word: b] # B[word: a]
                    sage: x*y
                    B[word: ac] # B[word: bb] # B[word: ca]

                    sage: x = tensor( ((a+2*b), c) )    ; x
                    B[word: a] # B[word: c] + 2*B[word: b] # B[word: c]
                    sage: y = tensor( (c,       a) ) + 1; y
                    B[word: ] # B[word: ] + B[word: c] # B[word: a]
                    sage: x*y
                    B[word: a] # B[word: c] + B[word: ac] # B[word: ca] + 2*B[word: b] # B[word: c] + 2*B[word: bc] # B[word: ca]

                TODO: optimize this implementation!

                """
                return tensor( (module.monomial(x1)*module.monomial(x2) for (module, x1, x2) in zip(self.factors(), self.index_to_indices()(t1), self.index_to_indices()(t2)))) #.

            def _monomial_almost_one(self, i, x):
                ids = [algebra.one_basis() for algebra in self.factors()]
                ids[i] = x
                return self.monomial(tuple(ids))

            @cached_method
            def factor_embedding(self, i):
                r"""
                The algebra embedding from the `i`-th tensor factor of ``self`` to ``self``.

                EXAMPLES::

                    sage: W = WeylGroup("A2",prefix="s")
                    sage: A = W.algebra(ZZ)
                    sage: A3 = tensor([A,A,A])
                    sage: A3.factor_embedding(1)(A.monomial(W.from_reduced_word([2,1])))
                    B[1] # B[s2*s1] # B[1]

                """
                if i not in range(len(self.factors())):
                    raise ValueError, "index is out of range"
                return self.factors()[i].module_morphism(on_basis=partial(self._monomial_almost_one,i),codomain=self)

        class ElementMethods:
            """
            Implements operations on elements of tensor products of algebras with basis
            """
            pass
