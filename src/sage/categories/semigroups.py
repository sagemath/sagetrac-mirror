r"""
Semigroups
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Florent Hivert <florent.hivert at univ-rouen.fr>
#                2008-2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.misc.misc_c import prod
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.subquotients import SubquotientsCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.quotients import QuotientsCategory
from sage.categories.magmas import Magmas
from sage.structure.element import generic_power

class Semigroups(CategoryWithAxiom):
    """
    The category of (multiplicative) semigroups, i.e. sets with an
    associative operation ``*``.

    EXAMPLES::

        sage: C = Semigroups(); C
        Category of semigroups
        sage: C.super_categories()
        [Category of magmas]
        sage: C.all_super_categories()
        [Category of semigroups, Category of magmas,
         Category of sets, Category of sets with partial maps, Category of objects]
        sage: C.axioms()
        frozenset(['Associative'])
        sage: C.example()
        An example of a semigroup: the left zero semigroup

    TESTS::

        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = [Magmas, "Associative"]

    def example(self, choice="leftzero", **kwds):
        r"""
        Returns an example of a semigroup, as per
        :meth:`Category.example()
        <sage.categories.category.Category.example>`.

        INPUT::

         - ``choice`` -- str [default: 'leftzero']. Can be either 'leftzero'
           for the left zero semigroup, or 'free' for the free semigroup.
         - ``**kwds`` -- keyword arguments passed onto the constructor for the
           chosen semigroup.

        EXAMPLES::

            sage: Semigroups().example(choice='leftzero')
            An example of a semigroup: the left zero semigroup
            sage: Semigroups().example(choice='free')
            An example of a semigroup: the free semigroup generated by ('a', 'b', 'c', 'd')
            sage: Semigroups().example(choice='free', alphabet=('a','b'))
            An example of a semigroup: the free semigroup generated by ('a', 'b')

        """
        import sage.categories.examples.semigroups as examples
        if choice == "leftzero":
            return examples.LeftZeroSemigroup(**kwds)
        else:
            return examples.FreeSemigroup(**kwds)

    class ParentMethods:

        def _test_associativity(self, **options):
            r"""
            Test associativity for (not necessarily all) elements of this
            semigroup.

            INPUT::

             - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

            EXAMPLES:

            By default, this method tests only the elements returned by
            ``self.some_elements()``::

                sage: L = Semigroups().example(choice='leftzero')
                sage: L._test_associativity()

            However, the elements tested can be customized with the
            ``elements`` keyword argument::

                sage: L._test_associativity(elements = (L(1), L(2), L(3)))

            See the documentation for :class:`TestSuite` for more information.

            """
            tester = self._tester(**options)
            S = tester.some_elements()
            from sage.combinat.cartesian_product import CartesianProduct
            for x,y,z in tester.some_elements(CartesianProduct(S,S,S)):
                tester.assert_((x * y) * z == x * (y * z))

        def prod(self, args):
            r"""
            Returns the product of the list of elements ``args`` inside ``self``.

            EXAMPLES::

                sage: S = Semigroups().example("free")
                sage: S.prod([S('a'), S('b'), S('c')])
                'abc'
                sage: S.prod([])
                Traceback (most recent call last):
                ...
                AssertionError: Cannot compute an empty product in a semigroup
            """
            assert len(args) > 0, "Cannot compute an empty product in a semigroup"
            return prod(args[1:], args[0])

        def cayley_graph(self, side="right", simple=False, elements = None, generators = None, connecting_set = None):
            r"""

            Returns the Cayley graph for this finite semigroup.

            INPUT::

             - ``side`` -- "left", "right", or "twosided":
               the side on which the generators act (default:"right")
             - ``simple`` -- boolean (default:False):
               if True, returns a simple graph (no loops, no labels,
               no multiple edges)
             - ``generators`` -- a list, tuple, or family of elements of ``self``
               (default: ``self.semigroup_generators()``)
             - ``connecting_set`` -- alias for ``generators``; deprecated
             - ``elements`` -- a list (or iterable) of elements of ``self``

            OUTPUT::

             - :class:`DiGraph`

            EXAMPLES:

            We start with the (right) Cayley graphs of some classical groups::

                sage: D4 = DihedralGroup(4); D4
                Dihedral group of order 8 as a permutation group
                sage: G = D4.cayley_graph()
                sage: show(G, color_by_label=True, edge_labels=True)
                sage: A5 = AlternatingGroup(5); A5
                Alternating group of order 5!/2 as a permutation group
                sage: G = A5.cayley_graph()
                sage: G.show3d(color_by_label=True, edge_size=0.01, edge_size2=0.02, vertex_size=0.03)
                sage: G.show3d(vertex_size=0.03, edge_size=0.01, edge_size2=0.02, vertex_colors={(1,1,1):G.vertices()}, bgcolor=(0,0,0), color_by_label=True, xres=700, yres=700, iterations=200) # long time (less than a minute)
                sage: G.num_edges()
                120

                sage: w = WeylGroup(['A',3])
                sage: d = w.cayley_graph(); d
                Digraph on 24 vertices
                sage: d.show3d(color_by_label=True, edge_size=0.01, vertex_size=0.03)

            Alternative generators may be specified::

                sage: G = A5.cayley_graph(generators=[A5.gens()[0]])
                sage: G.num_edges()
                60
                sage: g=PermutationGroup([(i+1,j+1) for i in range(5) for j in range(5) if j!=i])
                sage: g.cayley_graph(generators=[(1,2),(2,3)])
                Digraph on 120 vertices

            If ``elements`` is specified, then only the subgraph
            induced and those elements is returned. Here we use it to
            display the Cayley graph of the free monoid truncated on
            the elements of length at most 3::

                sage: M = Monoids().example(); M
                An example of a monoid: the free monoid generated by ('a', 'b', 'c', 'd')
                sage: elements = [ M.prod(w) for w in sum((list(Words(M.semigroup_generators(),k)) for k in range(4)),[]) ]
                sage: G = M.cayley_graph(elements = elements)
                sage: G.num_verts(), G.num_edges()
                (85, 84)
                sage: G.show3d(color_by_label=True, edge_size=0.001, vertex_size=0.01)

            We now illustrate the ``side`` and ``simple`` options on a semigroup::

                sage: S = FiniteSemigroups().example(alphabet=('a','b'))
                sage: g = S.cayley_graph(simple=True)
                sage: g.vertices()
                ['a', 'ab', 'b', 'ba']
                sage: g.edges()
                [('a', 'ab', None), ('b', 'ba', None)]

            ::

                sage: g = S.cayley_graph(side="left", simple=True)
                sage: g.vertices()
                ['a', 'ab', 'b', 'ba']
                sage: g.edges()
                [('a', 'ba', None), ('ab', 'ba', None), ('b', 'ab', None),
                ('ba', 'ab', None)]

            ::

                sage: g = S.cayley_graph(side="twosided", simple=True)
                sage: g.vertices()
                ['a', 'ab', 'b', 'ba']
                sage: g.edges()
                [('a', 'ab', None), ('a', 'ba', None), ('ab', 'ba', None),
                ('b', 'ab', None), ('b', 'ba', None), ('ba', 'ab', None)]

            ::

                sage: g = S.cayley_graph(side="twosided")
                sage: g.vertices()
                ['a', 'ab', 'b', 'ba']
                sage: g.edges()
                [('a', 'a', (0, 'left')), ('a', 'a', (0, 'right')), ('a', 'ab', (1, 'right')), ('a', 'ba', (1, 'left')), ('ab', 'ab', (0, 'left')), ('ab', 'ab', (0, 'right')), ('ab', 'ab', (1, 'right')), ('ab', 'ba', (1, 'left')), ('b', 'ab', (0, 'left')), ('b', 'b', (1, 'left')), ('b', 'b', (1, 'right')), ('b', 'ba', (0, 'right')), ('ba', 'ab', (0, 'left')), ('ba', 'ba', (0, 'right')), ('ba', 'ba', (1, 'left')), ('ba', 'ba', (1, 'right'))]

            ::

                sage: s1 = SymmetricGroup(1); s = s1.cayley_graph(); s.vertices()
                [()]

            TESTS::

                sage: SymmetricGroup(2).cayley_graph(side="both")
                Traceback (most recent call last):
                ...
                ValueError: option 'side' must be 'left', 'right' or 'twosided'

            TODO:

             - Add more options for constructing subgraphs of the
               Cayley graph, handling the standard use cases when
               exploring large/infinite semigroups (a predicate,
               generators of an ideal, a maximal length in term of the
               generators)

             - Specify good default layout/plot/latex options in the graph

             - Generalize to combinatorial modules with module generators / operators

            AUTHORS:

             - Bobby Moretti (2007-08-10)
             - Robert Miller (2008-05-01): editing
             - Nicolas M. Thiery (2008-12): extension to semigroups,
               ``side``, ``simple``, and ``elements`` options, ...
            """
            from sage.graphs.digraph import DiGraph
            from groups import Groups
            if not side in ["left", "right", "twosided"]:
                raise ValueError, "option 'side' must be 'left', 'right' or 'twosided'"
            if elements is None:
                assert self.is_finite(), "elements should be specified for infinite semigroups"
                elements = list(self)
            elements_set = set(elements)
            if simple or self in Groups():
                result = DiGraph()
            else:
                result = DiGraph(multiedges = True, loops = True)
            result.add_vertices(elements)

            if connecting_set is not None:
                generators = connecting_set
            if generators is None:
                generators = self.semigroup_generators()
            if isinstance(generators, (list, tuple)):
                generators = dict((self(g), self(g)) for g in generators)
            left  = (side == "left"  or side == "twosided")
            right = (side == "right" or side == "twosided")
            def add_edge(source, target, label, side_label):
                """
                Skips edges whose targets are not in elements
                Return an appropriate edge given the options
                """
                if target not in elements_set: return
                if simple:
                    result.add_edge([source, target])
                elif side == "twosided":
                    result.add_edge([source, target, (label, side_label)])
                else:
                    result.add_edge([source, target, label])
            for x in elements:
                for i in generators.keys():
                    if left:
                        add_edge(x, generators[i]*x, i, "left" )
                    if right:
                        add_edge(x, x*generators[i], i, "right")
            return result

    class ElementMethods:

        def _pow_(self, n):
            """
            Returns self to the $n^{th}$ power.

            INPUT::

             - ``n`` -- a positive integer

            EXAMPLES::

                sage: S = Semigroups().example("leftzero")
                sage: x = S("x")
                sage: x^1, x^2, x^3, x^4, x^5
                ('x', 'x', 'x', 'x', 'x')
                sage: x^0
                Traceback (most recent call last):
                ...
                AssertionError

            TESTS::

                sage: x._pow_(17)
                'x'

            """
            assert n > 0
            return generic_power(self, n)

        __pow__ = _pow_


    Finite = LazyImport('sage.categories.finite_semigroups', 'FiniteSemigroups', 'Finite', at_startup=True)
    Unital = LazyImport('sage.categories.monoids', 'Monoids', 'Unital', at_startup=True)

    #######################################
    class Subquotients(SubquotientsCategory):
        r"""
        The category of sub/quotient semi-groups.

        EXAMPLES::

            sage: Semigroups().Subquotients().all_super_categories()
            [Category of subquotients of semigroups,
             Category of semigroups,
             Category of subquotients of magmas,
             Category of magmas,
             Category of subquotients of sets,
             Category of sets,
             Category of sets with partial maps,
             Category of objects]

            [Category of subquotients of semigroups,
             Category of semigroups,
             Category of subquotients of magmas,
             Category of magmas,
             Category of subquotients of sets,
             Category of sets,
             Category of sets with partial maps,
             Category of objects]
        """

        def example(self):
            """
            Returns an example of sub quotient of a semigroup, as per
            :meth:`Category.example()
            <sage.categories.category.Category.example>`.

            EXAMPLES::

                sage: Semigroups().Subquotients().example()
                An example of a (sub)quotient semigroup: a quotient of the left zero semigroup
            """
            from sage.categories.examples.semigroups import QuotientOfLeftZeroSemigroup
            return QuotientOfLeftZeroSemigroup(category = self.Subquotients())

    class Quotients(QuotientsCategory):

        def example(self):
            r"""
            Returns an example of quotient of a semigroup, as per
            :meth:`Category.example()
            <sage.categories.category.Category.example>`.

            EXAMPLES::

                sage: Semigroups().Quotients().example()
                An example of a (sub)quotient semigroup: a quotient of the left zero semigroup
            """
            from sage.categories.examples.semigroups import QuotientOfLeftZeroSemigroup
            return QuotientOfLeftZeroSemigroup()

        class ParentMethods:

            def semigroup_generators(self):
                r"""
                Returns semigroup generators for ``self`` by
                retracting the semigroup generators of the ambient
                semigroup.

                EXAMPLES::

                    sage: S = FiniteSemigroups().Quotients().example().semigroup_generators() # todo: not implemented
                """
                return self.ambient().semigroup_generators().map(self.retract)

    class CartesianProducts(CartesianProductsCategory):

        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Semigroups().CartesianProducts().extra_super_categories()
                [Category of semigroups]
                sage: Semigroups().CartesianProducts().super_categories()
                [Category of semigroups, Category of Cartesian products of magmas]
            """
            return [Semigroups()]

    class Algebras(AlgebrasCategory):
        """
        TESTS::

            sage: TestSuite(Semigroups().Algebras(QQ)).run()
            sage: TestSuite(Semigroups().Finite().Algebras(QQ)).run()
        """

        def extra_super_categories(self):
            """
            Implements the fact that the algebra of a semigroup is indeed a (not necessarily unital) algebra

            EXAMPLES::

                sage: Semigroups().Algebras(QQ).extra_super_categories()
                [Category of semigroups]
                sage: Semigroups().Algebras(QQ).super_categories()
                [Category of associative algebras over Rational Field,
                 Category of magma algebras over Rational Field]
            """
            return [Semigroups()]

        class ParentMethods:

            @cached_method
            def algebra_generators(self):
                r"""
                The generators of this algebra, as per
                :meth:`MagmaticAlgebras.ParentMethods.algebra_generators()
                <.magmatic_algebras.MagmaticAlgebras.ParentMethods.algebra_generators>`.

                They correspond to the generators of the semigroup.

                EXAMPLES::

                    sage: A = FiniteSemigroups().example().algebra(ZZ)
                    sage: A.algebra_generators()
                    Finite family {0: B['a'], 1: B['b'], 2: B['c'], 3: B['d']}
                """
                return self.basis().keys().semigroup_generators().map(self.monomial)

            def product_on_basis(self, g1, g2):
                r"""
                Product, on basis elements, as per
                :meth:`MagmaticAlgebras.WithBasis.ParentMethods.product_on_basis()
                <.magmatic_algebras.MagmaticAlgebras.WithBasis.ParentMethods.product_on_basis>`.

                The product of two basis elements is induced by the
                product of the corresponding elements of the group.

                EXAMPLES::

                    sage: S = FiniteSemigroups().example(); S
                    An example of a finite semigroup: the left regular band generated by ('a', 'b', 'c', 'd')
                    sage: A = S.algebra(QQ)
                    sage: a,b,c,d = A.algebra_generators()
                    sage: a * b + b * d * c * d
                    B['ab'] + B['bdc']
                """
                return self.monomial(g1 * g2)

