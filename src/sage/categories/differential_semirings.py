r"""
Differential Semirings
"""
#*****************************************************************************
#  Copyright (C) 2014 Miguel Marco, John Palmieri, Travis Scrimshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_singleton
from sage.categories.category_singleton import Category_singleton
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory
from sage.categories.semirings import Semirings
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_class_attribute

# Should this be CategoryWithAxiom_singleton?
class DifferentialSemiringsCategory(RegressiveCovariantConstructionCategory, Category_singleton):
    def __init__(self, base_category):
        """
        EXAMPLES::

            sage: C = Semirings().Differential()
            sage: C
            Category of graded algebras over Rational Field
            sage: C.base_category()
            Category of algebras over Rational Field
            sage: sorted(C.super_categories(), key=str)
            [Category of algebras over Rational Field,
             Category of graded modules over Rational Field]

            sage: AlgebrasWithBasis(QQ).Differential().base_ring()
            Rational Field
            sage: GradedHopfAlgebrasWithBasis(QQ).base_ring()
            Rational Field
        """
        super(DifferentialSemiringsCategory, self).__init__(base_category)

    _functor_category = "Differential"

    @lazy_class_attribute
    def _base_category_class(cls):
        """
        Recover the class of the base category.

        OUTPUT:

        A *tuple* whose first entry is the base category class.

        .. WARNING::

            This is only used for graded categories that are not
            implemented as nested classes, and won't work otherwise.

        .. SEEALSO:: :meth:`__classcall__`

        EXAMPLES::

            sage: DifferentialAlgebras._base_category_class
            (<class 'sage.categories.modules.Modules'>,)
            sage: DifferentialAlgebrasWithBasis._base_category_class
            (<class 'sage.categories.algebras_with_basis.AlgebrasWithBasis'>,)

        The reason for wrapping the base category class in a tuple is
        that, often, the base category class implements a
        :meth:`__classget__` method which would get in the way upon
        attribute access::

                sage: F = DifferentialAlgebras
                sage: F._foo = F._base_category_class[0]
                sage: F._foo
                Traceback (most recent call last):
                ...
                AssertionError: base category class for <...Algebras'> mismatch;
                expected <...Algebras'>, got <...DifferentialAlgebras'>
        """
        module_name = cls.__module__.replace("differential_","")
        import sys
        name   = cls.__name__.replace("Differential","")
        __import__(module_name)
        module = sys.modules[module_name]
        return (module.__dict__[name],)

    @staticmethod
    def __classcall__(cls, category, *args):
        """
        Magic support for putting Differential categories in their own file.

        EXAMPLES::

            sage: from sage.categoires.differential_algebras import DifferentialAlgebras
            sage: DifferentialAlgebras(ZZ)   # indirect doctest
            Category of graded modules over Integer Ring
            sage: Algebras(ZZ).Differential()
            Category of graded modules over Integer Ring
            sage: DifferentialAlgebras(ZZ) is Algebras(ZZ).Differential()
            True

        .. TODO::

            Generalize this support for all other functorial
            constructions if at some point we have a category ``Blah`` for
            which we want to implement the construction ``Blah.Foo`` in a
            separate file like we do for e.g. :class:`DifferentialAlgebras`,
            :class:`GradedAlgebras`, ...

        .. SEEALSO:: :meth:`_base_category_class`
        """
        base_category_class = cls._base_category_class[0]
        if isinstance(category, base_category_class):
            return super(DifferentialSemiringsCategory, cls).__classcall__(cls, category, *args)
        return base_category_class(category, *args).Differential()

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Algebra(QQ).Differnetial()  # indirect doctest
            Category of differential algebras over Rational Field
        """
        return "differential {}".format(self.base_category()._repr_object_names())

class DifferentialSemirings(DifferentialSemiringsCategory):
    r"""
    The category of differential semirings.

    A differential semiring `S` is a semiring with a specified derivation
    `\partial : S \to S`, that is a linear map which satisfies the
    Leibniz rule:

    .. MATH::

        \partial(xy) = \partial(x) y + x \partial(y).

    .. SEEALSO: :wikipedia:`Differential_algebra`

    EXAMPLES::

        sage: C = Semirings().Differential()
        Category of differential semirings
        sage: C.super_categories()
        [Category of associative additive commutative additive associative additive unital distributive magmas and additive magmas,
         Category of monoids]

        sage: C.AdditiveInverse()
        Category of differential rings

    TESTS::

        sage: TestSuite(Semirings().Differential()).run()
    """
    _base_category_class_and_axiom = (Semirings, "")

    class ParentMethods:

        @abstract_method
        def differential(self, x):
            """
            Return the differential of ``x`` in ``self``.

            It is recommended that this should be a ``@lazy_attribute``
            which is the differential as a morphism.
            """

        @abstract_method
        def degree_of_differential(self):
            """
            The degree of the differential in ``self``.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
                sage: G = AdditiveAbelianGroup((0,0))
                sage: x = G.gen(0); y = G.gen(1)
                sage: A = CommutativeDGA(degrees=(x, y), differential=({(1,1): 1}, None))
                sage: A.degree_of_differential()
                (0, 1)
            """

        def _bdry_matrix(self, n):
            """
            Matrix representing the differential with domain the
            degree ``n`` part of ``self``.

            In particular, the `j`-th column gives the boundary of the
            `jth` element in the basis of ``self``, expanded in terms of
            the basis in degree `n+d`, where `d` is the degree of the
            differential.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
                sage: E = CommutativeDGA(degrees=(1,1,1,1,1), differential=(None, None, None, {(1,1,0,0,0): 1}, {(0,1,1,0,0): 1}), generator_names=['a', 'b', 'c', 'x', 'y'])

                sage: E.homogeneous_component(1).dimension()
                5
                sage: E.homogeneous_component(2).dimension()
                10
                sage: E._bdry_matrix(1).nrows()
                5
                sage: E._bdry_matrix(1).ncols()
                10

                sage: E._bdry_matrix(2)
                [ 0  0  0 -1  0  0  0  1  0  0]
                [ 0  0  0  0  0  0  0  0  0  0]
                [ 0  0  0  0  0  0  0  0  0 -1]
                [ 0  0  0  0  0  0  0  0  0  0]
                [ 0  0  0  0  0  0  0  0  0  0]
                [ 0  0  0  0  0  0  0  0  0  0]
                [ 0  0  0  0  0  0  0  0  0 -1]
                [ 0  0  0  0  0  0  0  0  0  0]
                [ 0  0  0  0  0  0  0  0  0  0]
                [ 0  0  0  0  0  0  0  0  0  0]
            """
            from sage.matrix.constructor import matrix
            basis = self.basis(n)
            # this gives the standard basis (e_i) for B in degree n:
            std_basis = [self[n].basis()[m].to_vector() for m in basis]
            deg = self.degree_of_differential()
            if len(basis) > 0:
                basis_codomain = self.basis(n+deg)
                std_basis_codomain = [self[n+deg].basis()[m].to_vector()
                                      for m in basis_codomain]
                bdries = [self.differential_on_basis(v) for v in basis]

                return matrix(self.base_ring(), [self[n+deg].sum_of_terms(list(b.monomial_coefficients().items())).to_vector() for b in bdries])
            else:
                return matrix(self.base_ring(), 0, len(self.basis(n+deg)))

        def _raw_cycles(self, n):
            """
            Vector space of cycles in degree ``n``.

            This is presented as a subspace of an abstract vector space
            isomorphic to ``self`` in degree `n`. It is precisely the
            kernel of the boundary matrix in this degree.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
                sage: E = CommutativeDGA(degrees=(1,1,1,1,1), differential=(None, None, None, {(1,1,0,0,0): 1}, {(0,1,1,0,0): 1}), generator_names=['a', 'b', 'c', 'x', 'y'])
                sage: E._raw_cycles(2)
                Vector space of degree 10 and dimension 8 over Rational Field
                Basis matrix:
                [ 0  1  0  0  0  0  0  0  0  0]
                [ 0  0  1  0  0  0 -1  0  0  0]
                [ 0  0  0  1  0  0  0  0  0  0]
                [ 0  0  0  0  1  0  0  0  0  0]
                [ 0  0  0  0  0  1  0  0  0  0]
                [ 0  0  0  0  0  0  0  1  0  0]
                [ 0  0  0  0  0  0  0  0  1  0]
                [ 0  0  0  0  0  0  0  0  0  1]
            """
            return self._bdry_matrix(n).left_kernel()

        def cycles(self, n):
            """
            Vector space of cycles in degree ``n``.

            This is presented as a "subspace" of ``self`` in degree `n`,
            with a basis of polynomials in terms of the generators of this
            algebra. Compare to ``self._raw_cycles(n)``.

            Note that it is not really a subspace, since subspaces are not
            implemented for :class:`CombinatorialFreeModule`. But the
            elements of this vector space are represented as elements of
            the algebra.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
                sage: E = CommutativeDGA(degrees=(1,1,1,1,1), differential=(None, None, None, {(1,1,0,0,0): 1}, {(0,1,1,0,0): 1}), generator_names=['a', 'b', 'c', 'x', 'y'])
                sage: E._raw_cycles(2)
                Vector space of degree 10 and dimension 8 over Rational Field
                Basis matrix:
                [ 0  1  0  0  0  0  0  0  0  0]
                [ 0  0  1  0  0  0 -1  0  0  0]
                [ 0  0  0  1  0  0  0  0  0  0]
                [ 0  0  0  0  1  0  0  0  0  0]
                [ 0  0  0  0  0  1  0  0  0  0]
                [ 0  0  0  0  0  0  0  1  0  0]
                [ 0  0  0  0  0  0  0  0  1  0]
                [ 0  0  0  0  0  0  0  0  0  1]
                sage: E.cycles(2)
                Free module generated by {c_1 * y_1, c_1 * x_1 - a_1 * y_1, b_1 * y_1, b_1 * x_1, b_1 * c_1, a_1 * x_1, a_1 * c_1, a_1 * b_1} over Rational Field

            """
            raw_kernel_basis = self._raw_cycles(n).basis()
            return CombinatorialFreeModule(self.base_ring(),
                                               [self.sum_of_terms(list(self[n].from_vector(v).monomial_coefficients().items())) for v in raw_kernel_basis])

        def _raw_boundaries(self, n):
            """
            Vector space of boundaries in degree ``n``.

            This is presented as a subspace of an abstract vector space
            isomorphic to ``self`` in degree `n`. It is precisely the
            image of the boundary matrix in this degree.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
                sage: E = CommutativeDGA(degrees=(1,1,1,1,1), differential=(None, None, None, {(1,1,0,0,0): 1}, {(0,1,1,0,0): 1}), generator_names=['a', 'b', 'c', 'x', 'y'])
                sage: E._raw_boundaries(2)
                Vector space of degree 10 and dimension 2 over Rational Field
                Basis matrix:
                [0 0 0 0 0 1 0 0 0 0]
                [0 0 0 0 0 0 0 0 0 1]
            """
            return self._bdry_matrix(n-self.degree_of_differential()).row_space()

        def boundaries(self, n):
            """
            Vector space of boundaries in degree ``n``.

            This is presented as a "subspace" of ``self`` in degree `n`,
            with a basis of polynomials in terms of the generators of this
            algebra. Compare to ``self._raw_boundaries(n)``.

            Note also that it is not really a subspace, since subspaces
            are not implemented for :class:`CombinatorialFreeModule`. But
            the elements of this vector space are represented as elements
            of the algebra.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
                sage: E = CommutativeDGA(degrees=(1,1,1,1,1), differential=(None, None, None, {(1,1,0,0,0): 1}, {(0,1,1,0,0): 1}), generator_names=['a', 'b', 'c', 'x', 'y'])
                sage: E._raw_boundaries(2)
                Vector space of degree 10 and dimension 2 over Rational Field
                Basis matrix:
                [0 0 0 0 0 1 0 0 0 0]
                [0 0 0 0 0 0 0 0 0 1]
                sage: E.boundaries(2)
                Free module generated by {b_1 * c_1, a_1 * b_1} over Rational Field

            """
            raw_bdry_basis = self._raw_boundaries(n).basis()
            return CombinatorialFreeModule(self.base_ring(),
                                               [self.sum_of_terms(list(self[n].from_vector(v).monomial_coefficients().items())) for v in raw_bdry_basis])

        def _raw_homology(self, n):
            """
            The quotient cycles mod boundaries in degree ``n``.

            This is computed as
            ``self._raw_cycles(n) / self._raw_boundaries(n)``.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
                sage: E = CommutativeDGA(degrees=(1,1,1,1,1), differential=(None, None, None, {(1,1,0,0,0): 1}, {(0,1,1,0,0): 1}), generator_names=['a', 'b', 'c', 'x', 'y'])
                sage: H = E._raw_homology(2); H
                Vector space quotient V/W of dimension 6 over Rational Field where
                V: Vector space of degree 10 and dimension 8 over Rational Field
                Basis matrix:
                [ 0  1  0  0  0  0  0  0  0  0]
                [ 0  0  1  0  0  0 -1  0  0  0]
                [ 0  0  0  1  0  0  0  0  0  0]
                [ 0  0  0  0  1  0  0  0  0  0]
                [ 0  0  0  0  0  1  0  0  0  0]
                [ 0  0  0  0  0  0  0  1  0  0]
                [ 0  0  0  0  0  0  0  0  1  0]
                [ 0  0  0  0  0  0  0  0  0  1]
                W: Vector space of degree 10 and dimension 2 over Rational Field
                Basis matrix:
                [0 0 0 0 0 1 0 0 0 0]
                [0 0 0 0 0 0 0 0 0 1]
                sage: H.dimension()
                6
            """
            return self._raw_cycles(n)/self._raw_boundaries(n)

        def homology(self, n):
            """
            Vector space of homology in degree ``n``.

            This is a subquotient of ``self`` in degree `n`. It is
            presented as a "subspace" of self in degree `n`, with a basis
            of polynomials in terms of the generators of this algebra,
            obtained by lifting a vector space basis from the quotient
            cycles/boundaries to the vector space of cycles. Compare to
            ``self._raw_homology(n)``.

            Of course mathematically this is not a subspace; it is not one
            from the computer point of view, either, since subspaces are
            not implemented for :class:`CombinatorialFreeModule`. But the
            elements of this vector space are represented as elements of
            the algebra.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import CommutativeDGA
                sage: D = CommutativeDGA(degrees=(1,2), differential=(None, {(1,1): 1}))
                sage: D.homology(1)
                Free module generated by {a_1} over Rational Field
                sage: sum(D.homology(i).dimension() for i in range(2,10))
                0
                sage: E.<a,b,c,x,y> = CommutativeDGA(degrees=(1,1,1,1,1), differential=(None, None, None, {(1,1,0,0,0): 1}, {(0,1,1,0,0): 1}))
                sage: E.homology(2)
                Free module generated by {c_1 * y_1, c_1 * x_1 - a_1 * y_1, b_1 * y_1, b_1 * x_1, a_1 * x_1, a_1 * c_1} over Rational Field
            """
            H = self._raw_homology(n)
            H_basis = [H.lift(H.basis()[i]) for i in range(H.dimension())]
            return CombinatorialFreeModule(self.base_ring(),
                                           [self.sum_of_terms(self[n].from_vector(v).monomial_coefficients().items())
                                            for v in H_basis])

        def differential_matrix(self, n, D):
            """
            Return the matrix that gives the differential on
            the ``n``-th degree.

            INPUT:

            - ``n`` -- integer

            - ``D`` -- the differential

            EXAMPLES::

                sage: A.<x,y,z,t> = CDGAlgebra(GF(5), degrees=(2, 3, 2, 4))
                sage: D = A.differential({t: x*y, x: y, z: y})
                sage: A.differential_matrix(4, D)
                [0 1]
                [2 0]
                [1 1]
                [0 2]
                sage: A.homogeneous_part(4)
                [t, z^2, x*z, x^2]
                sage: A.homogeneous_part(5)
                [y*z, x*y]
                sage: D(t)
                x*y
                sage: D(z^2)
                2*y*z
                sage: D(x*z)
                x*y + y*z
                sage: D(x^2)
                2*x*y
            """
            dom = self.homogeneous_part(n)
            cod = self.homogeneous_part(n+1)
            cokeys = [a.lift().dict().keys()[0] for a in cod]
            m = matrix(self.base_ring(), len(dom), len(cod))
            for i in range(len(dom)):
                im = dom[i].differential(D)
                dic = im.lift().dict()
                for j in dic.keys():
                    k = cokeys.index(j)
                    m[i,k] = dic[j]
            return m

        def cohomology(self, n, differential):
            """
            Return the n'th cohomology group of the algebra. This is
            in fact a vector space over the base ring.

            INPUT:

            - ``n`` -- the degree in which the cohomology is computed

            - ``differential`` -- the differential to be used

            EXAMPLES::

                sage: A.<x,y,z,t> = CDGAlgebra(QQ, degrees = (2, 3, 2, 4))
                sage: D = A.differential({t: x*y, x: y, z: y})
                sage: A.cohomology(4, D)
                Vector space quotient V/W of dimension 2 over Rational Field where
                V: Vector space of degree 4 and dimension 2 over Rational Field
                Basis matrix:
                [   1    0    0 -1/2]
                [   0    1   -2    1]
                W: Vector space of degree 4 and dimension 0 over Rational Field
                Basis matrix:
                []
            """
            N = self.differential_matrix(n, differential)
            if n == 1:
                M = matrix(self.base_ring(), 0, N.nrows())
            else:
                M = self.differential_matrix(n-1, differential)
            V0 = self.base_ring()**M.nrows()
            V1 = self.base_ring()**M.ncols()
            V2 = self.base_ring()**N.ncols()
            h1 = V0.Hom(V1)(M)
            h2 = V1.Hom(V2)(N)
            return h2.kernel().quotient(h1.image())

        def _test_differential(self, **options):
            """
            Verify that the differential of ``self`` satisfies the
            Leibniz rule.

            INPUT:

            - ``options`` -- any keyword arguments accepted
              by :meth:`_tester`
            """
            tester = self._tester(**options)
            D = self.differential
            S = self.some_elements()
            for x in S:
                for y in S:
                    tester.assertEquals(D(x*y), D(x)*y + x*D(y))

# This is the graded version
#        def _test_differential(self, **options):
#            """
#            Verify that the differential of ``self`` satisfies the
#            graded Leibniz rule.
#
#            INPUT:
#
#            - ``options`` -- any keyword arguments accepted
#              by :meth:`_tester`
#            """
#            tester = self._tester(**options)
#            D = self.differential
#            S = self.some_elements()
#            for x in S:
#                for y in S:
#                    tester.assertEquals(D(x*y), D(x)*y + (-1)**x.degree() * x*D(y))

    class ElementMethods:
        def differential(self):
            """
            Return the differential of ``self``.
            """
            return self.parent().differential(self)

        def is_boundary(self):
            """
            True iff ``self`` is a boundary.

            Requires ``self`` to be homogeneous.

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import *
                sage: A.<a,b,c> = CommutativeDGA(degrees=(1,2,2))
                sage: B.<x,y,z> = CommutativeDGA(A, differential=(0, a*c, 0))
                sage: x.is_boundary()
                False
                sage: (x*z).is_boundary()
                True
                sage: (x*z+x*y).is_boundary()
                False
                sage: (x*z+y**2).is_boundary()
                Traceback (most recent call last):
                ...
                ValueError: This element is not homogeneous.
            """
            try:
                assert self.is_homogeneous()
            except AssertionError:
                raise ValueError('this element is not homogeneous')
            # To avoid taking the degree of 0, we special-case it.
            if self == 0:
                return True
            deg = self.degree()
            A = self.parent()
            return A.homogeneous_component(deg)._from_dict(self.monomial_coefficients()).to_vector() in A._raw_boundaries(deg)

        def is_homologous_to(self, other):
            """
            Return True if ``self`` is homologous to ``other``.

            INPUT:

            - ``other``, another element of this DGA

            EXAMPLES::

                sage: from sage.algebras.differential_graded_algebras.commutative_dga import *
                sage: A.<a,b,c,d> = CommutativeDGA(degrees=(1,1,1,1))
                sage: B.<w,x,y,z> = CommutativeDGA(A, differential=(b*c-c*d, 0, 0, 0))
                sage: (x*y).is_homologous_to(y*z)
                True
                sage: (x*y).is_homologous_to(x*z)
                False
                sage: (x*y).is_homologous_to(x*y)
                True

            Two elements in different degrees are homologous if and
            only if they are both boundaries::

                sage: w.is_homologous_to(y*z)
                False
                sage: (x*y-y*z).is_homologous_to(x*y*z)
                True
                sage: (x*y*z).is_homologous_to(B(0)) # make sure 0 works
                True
            """
            try:
                assert isinstance(other, DGA_generic.Element) and self.parent() is other.parent()
            except AssertionError:
                raise ValueError('the element {} does not lie in this DGA'.format(other))
            if (self - other).is_homogeneous():
                return (self - other).is_boundary()

            return (self.is_boundary() and other.is_boundary())

