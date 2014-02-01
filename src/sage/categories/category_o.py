r"""
Category `\mathcal{O}`

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""
#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.category_types import Category_over_base
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.sets_cat import Sets
from sage.categories.homset import Hom
from sage.categories.morphism import Morphism
from sage.structure.sage_object import have_same_parent
from sage.graphs.dot2tex_utils import have_dot2tex

# Use SetsWithAction as common basepoint for this and crystals (finite_semigroups-nt.patch)
class CategoryO(Category_over_base):
    r"""
    The category `\mathcal{O}`.

    The category `\mathcal{O}` is the category of weight modules `V`
    over a Kac-Moody algebra `\mathfrak{g}` with finite dimensional weight
    spaces for which there exists `\lambda_1, \ldots, \lambda_s \in
    \mathfrak{h}^*` such that

    .. MATH::

        \mathrm{wt}(V) \subset D(\lambda_1) \cup \cdots \cup D(\lambda_s)

    where `D(\lambda) = \{\mu \in \mathfrak{h}^* \mid \lambda - \mu \in Q^+\}`
    and `Q^+` is the positive root lattice of `\mathfrak{g}`.

    EXAMPLES::

        sage: C = CategoryO(QQ); C
        Category O over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of algebras over Rational Field,
         Category of non associative algebras with basis over Rational Field,
         Category of non unital algebras with basis over Rational Field]

    We construct a typical parent in this category, and do some computations
    with it::

        sage: A = C.example(); A
        An example of a Lie algebra: the abelian Lie algebra on the generators ('a', 'b', 'c') over Rational Field

        sage: A.category()
        Category of algebras with basis over Rational Field

        sage: TestSuite(A).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
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
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_zero() . . . pass
        sage: A.__class__
        <class 'sage.categories.examples.lie_algebras.AbelianLieAlgebra_with_category'>
        sage: A.element_class
        <class 'sage.algeras.lie_algebras.AbelianLieAlgebra_with_category.element_class'>

    Please see the source code of `A` (with ``A??``) for how to
    implement other objects in category `\mathcal{O}`.

    TESTS::

        sage: TestSuite(CategoryO(QQ)).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: CategoryO(QQ).super_categories()
            [Category of modules with basis over Rational Field]
        """
        return [ModulesWithBasis(self.base_ring())]

    def example(self):
        """
        Return an example of algebra with basis as per
        :meth:`Category.example <sage.categories.category.Category.example>`.

        ::

            sage: CategoryO(QQ).example()
        """
        raise NotImplementedError

    def base_ring(self):
        """
        Return the base ring of ``self`` which is the base ring of the base
        of ``self``.
        """
        return self.base().base_ring()

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):
        def extra_super_categories(self):
            """
            Implements the fact that a finite dimensional representations over
            a finite ring is finite.

            EXAMPLES::

                sage: Modules(IntegerModRing(4)).FiniteDimensional().extra_super_categories()
                [Category of finite sets]
                sage: Modules(ZZ).FiniteDimensional().extra_super_categories()
                []
                sage: Modules(GF(5)).FiniteDimensional().is_subcategory(Sets().Finite())
                True
                sage: Modules(ZZ).FiniteDimensional().is_subcategory(Sets().Finite())
                False
            """
            if self.base_ring() in Sets().Finite():
                return [Sets().Finite()]
            return []

    class ParentMethods:
        def __iter__(self, index_set=None, max_depth=float('inf')):
            """
            Iterate over ``self``.

            INPUT:

            - ``index_set`` -- (Default: ``None``) The index set; if ``None``
              then use the index set of the crystal

            - ``max_depth`` -- (Default: infinity) The maximum depth to build

            EXAMPLES::
            """
            if index_set is None:
                index_set = self.index_set()
            current_level = (self.highest_weight_vector(),)
            known = set(current_level)
            depth = 0
            while len(current_level) > 0 and depth <= max_depth:
                next_level = set()
                for x in current_level:
                    yield x
                    for i in index_set:
                        y = x.f(i)
                        if y == 0:
                            continue
                        # TODO: Fix this HACK of knowing the element is an ElementWrapper subclass
                        y = y.normalize()
                        if y in known:
                            continue
                        next_level.add(y)
                        known.add(y)
                current_level = next_level
                depth += 1

        def _latex_(self, **options):
            r"""
            Return the crystal graph as a latex string. This can be exported
            to a file with ``self.latex_file('filename')``.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2],shape=[1])
                sage: T._latex_()   #optional - dot2tex
                '...tikzpicture...'
                sage: view(T, pdflatex = True, tightpage = True) #optional - dot2tex graphviz

            One can for example also color the edges using the following options::

                sage: T = CrystalOfTableaux(['A',2],shape=[1])
                sage: T._latex_(color_by_label = {0:"black", 1:"red", 2:"blue"})   #optional - dot2tex graphviz
                '...tikzpicture...'
            """
            if not have_dot2tex():
                print "dot2tex not available.  Install after running \'sage -sh\'"
                return
            G = self.digraph()
            G.set_latex_options(**options)
            return G._latex_()

        def index_set(self):
            """
            Return the index set of ``self``.
            """
            return self.base().cartan_type().index_set()

        def digraph(self, subset=None, index_set=None):
            """
            Return the DiGraph associated to ``self``.

            INPUT:

            - ``subset`` -- (Optional) A subset of vertices for
              which the digraph should be constructed

            - ``index_set`` -- (Optional) The index set to draw arrows

            EXAMPLES::
            """
            from sage.graphs.all import DiGraph
            from sage.categories.highest_weight_crystals import HighestWeightCrystals
            d = {}
            if self in HighestWeightCrystals:
                f = lambda (u,v,label): ({})
            else:
                f = lambda (u,v,label): ({"backward":label ==0})

            # Parse optional arguments
            if subset is None:
                subset = self
            if index_set is None:
                index_set = self.index_set()

            for x in subset:
                d[x] = {}
                for i in index_set:
                    child = x.f(i)
                    if child == 0:
                        continue
                    # TODO: Fix this HACK of knowing the element is an ElementWrapper subclass
                    child = child.normalize()
                    if child not in subset:
                        continue
                    d[x][child] = i
            G = DiGraph(d)
            if have_dot2tex():
                G.set_latex_options(format="dot2tex",
                                    edge_labels = True,
                                    color_by_label = self.base().cartan_type()._index_set_coloring,
                                    edge_options = f)
            return G

        def plot(self, **options):
            """
            Returns the plot of self as a directed graph.

            EXAMPLES::
            """
            return self.digraph().plot(edge_labels=True, vertex_size=0, **options)

        @abstract_method(optional=True)
        def e_on_basis(self, i, b):
            """
            Apply the action of `e_i` to the basis indexed by ``b``.
            """

        @abstract_method(optional=True)
        def f_on_basis(self, i, b):
            """
            Apply the action of `f_i` to the basis indexed by ``b``.
            """

        @abstract_method(optional=True)
        def h_on_basis(self, i, b):
            """
            Apply the action of `h_i` to the basis indexed by ``b``.
            """

        @abstract_method(optional=True)
        def d_on_basis(self, i, b):
            """
            Apply the action of `d_i` to the basis indexed by ``b``.
            Not needed for finite types since there are no `d_i` generators.
            """

    class ElementMethods:
        def e(self, i):
            """
            Apply the action of `e_i` to ``self``.
            """
            P = self.parent()
            return P.sum(c * P.e_on_basis(i, b) for b,c in self)

        def f(self, i):
            """
            Apply the action of `f_i` to ``self``.
            """
            P = self.parent()
            return P.sum(c * P.f_on_basis(i, b) for b,c in self)

        def h(self, i):
            """
            Apply the action of `h_i` to ``self``.
            """
            P = self.parent()
            return P.sum(c * P.h_on_basis(i, b) for b,c in self)

        def d(self, i):
            """
            Apply the action of `d_i` to ``self``.
            Not needed for finite types since there are no `d_i` generators.
            """
            P = self.parent()
            return P.sum(c * P.d_on_basis(i, b) for b,c in self)

#        @abstract_method
#        def cartan(self, h):
#            """
#            Apply the action of a `h \in \mathfrak{h}` to ``self``.
#            """

        def is_highest_weight(self, index_set=None):
            """
            Check if ``self`` is highest weight with respect to the given
            index set.
            """
            if index_set is None:
                index_set = self.parent()._g.cartan_type().index_set()
            return all(self.e(i) == 0 for i in index_set)

# TODO: Should have a category for highest weight representations
class CategoryOInt(Category_over_base):
    r"""
    The category `\mathcal{O}_{\mathrm{int}}`.

    The category `\mathcal{O}_{\mathrm{int}}` is the category of integrable
    objects in `\mathcal{O}`. A weight module `V` is called integrable if
    all `e_i` and `f_i` are locally nilpotent on `V`. We say an element `x`
    is locally nilpotent if for any `v \in V`, there exists an `N > 0` such
    that `x^N \cdot v = 0`.

    EXAMPLES::

        sage: C = CategoryOInt(QQ); C
        Category O int over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of algebras over Rational Field,
         Category of non associative algebras with basis over Rational Field,
         Category of non unital algebras with basis over Rational Field]

    We construct a typical parent in this category, and do some computations with it::

        sage: A = C.example(); A
        An example of a Lie algebra: the abelian Lie algebra on the generators ('a', 'b', 'c') over Rational Field

        sage: A.category()
        Category of algebras with basis over Rational Field

        sage: TestSuite(A).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
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
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_zero() . . . pass
        sage: A.__class__
        <class 'sage.categories.examples.lie_algebras.AbelianLieAlgebra_with_category'>
        sage: A.element_class
        <class 'sage.algeras.lie_algebras.AbelianLieAlgebra_with_category.element_class'>

    Please see the source code of `A` (with ``A??``) for how to
    implement other objects in category `\mathcal{O}`.

    TESTS::

        sage: TestSuite(CategoryO(QQ)).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: CategoryOIntegral(QQ).super_categories()
            [Category of modules with basis over Rational Field]
        """
        return [CategoryO(self.base())]

    class ParentMethods:
        pass

    class ElementMethods:
        # This should be in the highest weight category
        def to_highest_weight(self, index_set=None):
            r"""
            Return a path to a highest weight element from ``self``.

            Let `v \in V`, this returns the highest weight element `u` and
            a list `[i_1, \ldots, i_k]` such that
            `v = f_{i_1} \cdots f_{i_k} u`, where `i_1, \ldots, i_k` are
            elements in ``index_set``. By default the index set is assumed
            to be the full index set of ``self``.
            """
            from sage.categories.highest_weight_crystals import HighestWeightCrystals
            if index_set is None:
                index_set = self.index_set()
            for i in index_set:
                next = self.e(i)
                if next is not None:
                    hw = next.to_highest_weight(index_set = index_set)
                    return [hw[0], [i] + hw[1]]
            return [self, []]

