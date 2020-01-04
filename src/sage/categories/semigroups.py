r"""
Semigroups
"""
from __future__ import absolute_import
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Florent Hivert <florent.hivert at univ-rouen.fr>
#                2008-2015 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.misc.misc_c import prod
from sage.categories.category_with_axiom import CategoryWithAxiom, all_axioms
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.subquotients import SubquotientsCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.quotients import QuotientsCategory
from sage.categories.magmas import Magmas
from sage.arith.power import generic_power


all_axioms += ("HTrivial", "Aperiodic", "LTrivial", "RTrivial", "JTrivial")

class Semigroups(CategoryWithAxiom):
    """
    The category of (multiplicative) semigroups.

    A *semigroup* is an associative :class:`magma <Magmas>`, that is a
    set endowed with a multiplicative binary operation `*` which is
    associative (see :wikipedia:`Semigroup`).

    The operation `*` is not required to have a neutral element. A
    semigroup for which such an element exists is a :class:`monoid
    <sage.categories.monoids.Monoids>`.

    EXAMPLES::

        sage: C = Semigroups(); C
        Category of semigroups
        sage: C.super_categories()
        [Category of magmas]
        sage: C.all_super_categories()
        [Category of semigroups, Category of magmas,
         Category of sets, Category of sets with partial maps, Category of objects]
        sage: C.axioms()
        frozenset({'Associative'})
        sage: C.example()
        An example of a semigroup: the left zero semigroup

    TESTS::

        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (Magmas, "Associative")

    def example(self, choice="leftzero", **kwds):
        r"""
        Returns an example of a semigroup, as per
        :meth:`Category.example()
        <sage.categories.category.Category.example>`.

        INPUT:

        - ``choice`` -- str (default: 'leftzero'). Can be either 'leftzero'
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

            INPUT:

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`

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
            from sage.misc.misc import some_tuples
            for x, y, z in some_tuples(S, 3, tester._max_runs):
                tester.assertEqual((x * y) * z, x * (y * z))

        @abstract_method(optional=True)
        def semigroup_generators(self):
            """
            Return distinguished semigroup generators for ``self``.

            OUTPUT: a family

            This method is optional.

            EXAMPLES::

                sage: S = Semigroups().example("free"); S
                An example of a semigroup: the free semigroup generated by ('a', 'b', 'c', 'd')
                sage: S.semigroup_generators()
                Family ('a', 'b', 'c', 'd')
            """

        def magma_generators(self):
            """
            An alias for :meth:`semigroup_generators`.

            EXAMPLES::

                sage: S = Semigroups().example("free"); S
                An example of a semigroup: the free semigroup generated by ('a', 'b', 'c', 'd')
                sage: S.magma_generators()
                Family ('a', 'b', 'c', 'd')
                sage: S.semigroup_generators()
                Family ('a', 'b', 'c', 'd')
            """
            return self.semigroup_generators()

        def prod(self, args):
            r"""
            Return the product of the list of elements ``args``
            inside ``self``.

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
            Return the Cayley graph for this finite semigroup.

            INPUT:

            - ``side`` -- "left", "right", or "twosided":
              the side on which the generators act (default:"right")
            - ``simple`` -- boolean (default:False):
              if True, returns a simple graph (no loops, no labels,
              no multiple edges)
            - ``generators`` -- a list, tuple, or family of elements
              of ``self`` (default: ``self.semigroup_generators()``)
            - ``connecting_set`` -- alias for ``generators``; deprecated
            - ``elements`` -- a list (or iterable) of elements of ``self``

            OUTPUT:

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

            We now illustrate the ``side`` and ``simple`` options on
            a semigroup::

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

            .. TODO::

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
            from .monoids import Monoids
            from .groups import Groups
            if not side in ["left", "right", "twosided"]:
                raise ValueError("option 'side' must be 'left', 'right' or 'twosided'")
            if elements is None:
                assert self.is_finite(), "elements should be specified for infinite semigroups"
                elements = self
            else:
                elements = set(elements)
            if simple or self in Groups():
                result = DiGraph()
            else:
                result = DiGraph(multiedges = True, loops = True)
            result.add_vertices(elements)

            if connecting_set is not None:
                generators = connecting_set
            if generators is None:
                if self in Monoids and hasattr(self, "monoid_generators"):
                    generators = self.monoid_generators()
                else:
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
                if (elements is not self and
                    target not in elements):
                    return
                if simple:
                    if source != target:
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

        def subsemigroup(self, generators, one=None, category=None):
            r"""
            Return the multiplicative subsemigroup generated by ``generators``.

            INPUT:

            - ``generators`` -- a finite family of elements of
              ``self``, or a list, iterable, ... that can be converted
              into one (see :class:`Family`).

            - ``one`` -- a unit for the subsemigroup, or ``None``.

            - ``category`` -- a category

            This implementation lazily constructs all the elements of
            the semigroup, and the right Cayley graph relations
            between them, and uses the latter as an automaton.

            See :class:`~sage.sets.monoids.AutomaticSemigroup` for details.

            EXAMPLES::

                sage: R = IntegerModRing(15)
                sage: M = R.subsemigroup([R(3),R(5)]); M
                A subsemigroup of (Ring of integers modulo 15) with 2 generators
                sage: M.list()
                [3, 5, 9, 0, 10, 12, 6]

            By default, `M` is just in the category of subsemigroups::

                sage: M in Semigroups().Subobjects()
                True

            In the following example, we specify that `M` is a
            submonoid of the finite monoid `R` (it shares the same
            unit), and a group by itself::

                sage: M = R.subsemigroup([R(-1)],
                ....:     category=Monoids().Finite().Subobjects() & Groups()); M
                A submonoid of (Ring of integers modulo 15) with 1 generators
                sage: M.list()
                [1, 14]
                sage: M.one()
                1

            In the following example `M` is a group; however its unit
            does not coincide with that of `R`, so `M` is only a
            subsemigroup, and we need to specify its unit explicitly::

                sage: M = R.subsemigroup([R(5)],
                ....:     category=Semigroups().Finite().Subobjects() & Groups()); M
                Traceback (most recent call last):
                ...
                ValueError: For a monoid which is just a subsemigroup, the unit should be specified

                sage: M = R.subsemigroup([R(5)], one=R(10),
                ....:     category=Semigroups().Finite().Subobjects() & Groups()); M
                A subsemigroup of (Ring of integers modulo 15) with 1 generators
                sage: M in Groups()
                True
                sage: M.list()
                [10, 5]
                sage: M.one()
                10

            TESTS::

                sage: TestSuite(M).run()
                Failure in _test_inverse:
                Traceback (most recent call last):
                ...
                The following tests failed: _test_inverse

            .. TODO::

                - Fix the failure in TESTS by providing a default
                  implementation of ``__invert__`` for finite groups
                  (or even finite monoids).
                - Provide a default implementation of ``one`` for a
                  finite monoid, so that we would not need to specify
                  it explicitly?
            """
            from sage.monoids.automatic_semigroup import AutomaticSemigroup
            return AutomaticSemigroup(generators, ambient=self, one=one,
                                      category=category)

        def trivial_representation(self, base_ring=None, side="twosided"):
            r"""
            Return the trivial representation of ``self`` over ``base_ring``.

            INPUT:

            - ``base_ring`` -- (optional) the base ring; the default is `\ZZ`
            - ``side`` -- ignored

            EXAMPLES::

                sage: G = groups.permutation.Dihedral(4)
                sage: G.trivial_representation()
                Trivial representation of Dihedral group of order 8
                 as a permutation group over Integer Ring
            """
            if base_ring is None:
                from sage.rings.all import ZZ
                base_ring = ZZ
            from sage.modules.with_basis.representation import TrivialRepresentation
            return TrivialRepresentation(self, base_ring)

        def regular_representation(self, base_ring=None, side="left"):
            """
            Return the regular representation of ``self`` over ``base_ring``.

            - ``side`` -- (default: ``"left"``) whether this is the
              ``"left"`` or ``"right"`` regular representation

            EXAMPLES::

                sage: G = groups.permutation.Dihedral(4)
                sage: G.regular_representation()
                Left Regular Representation of Dihedral group of order 8
                 as a permutation group over Integer Ring
            """
            if base_ring is None:
                from sage.rings.all import ZZ
                base_ring = ZZ
            from sage.modules.with_basis.representation import RegularRepresentation
            return RegularRepresentation(self, base_ring, side)

    class ElementMethods:

        def _pow_int(self, n):
            """
            Return ``self`` to the `n^{th}` power.

            INPUT:

            - ``n`` -- a positive integer

            EXAMPLES::

                sage: S = Semigroups().example("leftzero")
                sage: x = S("x")
                sage: x^1, x^2, x^3, x^4, x^5
                ('x', 'x', 'x', 'x', 'x')
                sage: x^0
                Traceback (most recent call last):
                ...
                ArithmeticError: only positive powers are supported in a semigroup

            TESTS::

                sage: x._pow_int(17)
                'x'
            """
            if n <= 0:
                raise ArithmeticError("only positive powers are supported in a semigroup")
            return generic_power(self, n)


    class SubcategoryMethods:

        @cached_method
        def LTrivial(self):
            r"""
            Return the full subcategory of the `L`-trivial objects of ``self``.

            Let `S` be (multiplicative) :class:`semigroup <Semigroups>`.
            The `L`-*preorder* `\leq_L` on `S` is defined by:

            .. MATH::

                x\leq_L y \qquad \Longleftrightarrow \qquad x \in Sy

            The `L`-*classes* are the equivalence classes for the
            associated equivalence relation. The semigroup `S` is
            `L`-*trivial* if all its `L`-classes are trivial (that is
            of cardinality `1`), or equivalently if the `L`-preorder is
            in fact a partial order.

            EXAMPLES::

                sage: C = Semigroups().LTrivial(); C
                Category of l trivial semigroups

            A `L`-trivial semigroup is `H`-trivial::

                sage: sorted(C.axioms())
                ['Associative', 'HTrivial', 'LTrivial']

            .. SEEALSO::

                - :wikipedia:`Green's_relations`
                - :class:`Semigroups.SubcategoryMethods.RTrivial`
                - :class:`Semigroups.SubcategoryMethods.JTrivial`
                - :class:`Semigroups.SubcategoryMethods.HTrivial`

            TESTS::

                sage: TestSuite(C).run()
                sage: Rings().LTrivial.__module__
                'sage.categories.semigroups'
                sage: C                 # todo: not implemented
                Category of L-trivial semigroups
            """
            return self._with_axiom('LTrivial')

        @cached_method
        def RTrivial(self):
            r"""
            Return the full subcategory of the `R`-trivial objects of ``self``.

            Let `S` be (multiplicative) :class:`semigroup <Semigroups>`.
            The `R`-*preorder* `\leq_R` on `S` is defined by:

            .. MATH::

                x\leq_R y \qquad \Longleftrightarrow \qquad x \in yS

            The `R`-*classes* are the equivalence classes for the
            associated equivalence relation. The semigroup `S` is
            `R`-*trivial* if all its `R`-classes are trivial (that is
            of cardinality `1`), or equivalently if the `R`-preorder is
            in fact a partial order.

            EXAMPLES::

                sage: C = Semigroups().RTrivial(); C
                Category of r trivial semigroups

            An `R`-trivial semigroup is `H`-trivial::

                sage: sorted(C.axioms())
                ['Associative', 'HTrivial', 'RTrivial']

            .. SEEALSO::

                - :wikipedia:`Green's_relations`
                - :class:`Semigroups.SubcategoryMethods.LTrivial`
                - :class:`Semigroups.SubcategoryMethods.JTrivial`
                - :class:`Semigroups.SubcategoryMethods.HTrivial`

            TESTS::

                sage: TestSuite(C).run()
                sage: Rings().RTrivial.__module__
                'sage.categories.semigroups'
                sage: C                 # todo: not implemented
                Category of R-trivial semigroups
            """
            return self._with_axiom('RTrivial')

        @cached_method
        def JTrivial(self):
            r"""
            Return the full subcategory of the `J`-trivial objects of ``self``.

            Let `S` be (multiplicative) :class:`semigroup <Semigroups>`.
            The `J`-*preorder* `\leq_J` on `S` is defined by:

            .. MATH::

                x\leq_J y \qquad \Longleftrightarrow \qquad x \in SyS

            The `J`-*classes* are the equivalence classes for the
            associated equivalence relation. The semigroup `S` is
            `J`-*trivial* if all its `J`-classes are trivial (that is
            of cardinality `1`), or equivalently if the `J`-preorder is
            in fact a partial order.

            EXAMPLES::

                sage: C = Semigroups().JTrivial(); C
                Category of j trivial semigroups

            A semigroup is `J`-trivial if and only if it is
            `L`-trivial and `R`-trivial::

                sage: sorted(C.axioms())
                ['Associative', 'HTrivial', 'JTrivial', 'LTrivial', 'RTrivial']
                sage: Semigroups().LTrivial().RTrivial()
                Category of j trivial semigroups

            For a commutative semigroup, all three axioms are
            equivalent::

                sage: Semigroups().Commutative().LTrivial()
                Category of commutative j trivial semigroups
                sage: Semigroups().Commutative().RTrivial()
                Category of commutative j trivial semigroups

            .. SEEALSO::

                - :wikipedia:`Green's_relations`
                - :class:`Semigroups.SubcategoryMethods.LTrivial`
                - :class:`Semigroups.SubcategoryMethods.RTrivial`
                - :class:`Semigroups.SubcategoryMethods.HTrivial`

            TESTS::

                sage: TestSuite(C).run()
                sage: Rings().JTrivial.__module__
                'sage.categories.semigroups'
                sage: C                 # todo: not implemented
                Category of J-trivial semigroups
            """
            return self._with_axiom('JTrivial')

        @cached_method
        def HTrivial(self):
            r"""
            Return the full subcategory of the `H`-trivial objects of ``self``.

            Let `S` be (multiplicative) :class:`semigroup <Semigroups>`.
            Two elements of `S` are in the same `H`-class if they are
            in the same `L`-class and in the same `R`-class.

            The semigroup `S` is `H`-*trivial* if all its `H`-classes
            are trivial (that is of cardinality `1`).

            EXAMPLES::

                sage: C = Semigroups().HTrivial(); C
                Category of h trivial semigroups
                sage: Semigroups().HTrivial().Finite().example()
                NotImplemented

            .. SEEALSO::

                - :wikipedia:`Green's_relations`
                - :class:`Semigroups.SubcategoryMethods.RTrivial`
                - :class:`Semigroups.SubcategoryMethods.LTrivial`
                - :class:`Semigroups.SubcategoryMethods.JTrivial`
                - :class:`Semigroups.SubcategoryMethods.Aperiodic`

            TESTS::

                sage: TestSuite(C).run()
                sage: Rings().HTrivial.__module__
                'sage.categories.semigroups'
                sage: C                 # todo: not implemented
                Category of H-trivial semigroups
            """
            return self._with_axiom('HTrivial')

        @cached_method
        def Aperiodic(self):
            r"""
            Return the full subcategory of the aperiodic objects of ``self``.

            A (multiplicative) :class:`semigroup <Semigroups>` `S` is
            *aperiodic* if for any element `s\in S`, the sequence
            `s,s^2,s^3,...` eventually stabilizes.

            In terms of variety, this can be described by the equation
            `s^\omega s = s`.

            EXAMPLES::

                sage: Semigroups().Aperiodic()
                Category of aperiodic semigroups

            An aperiodic semigroup is `H`-trivial::

                sage: Semigroups().Aperiodic().axioms()
                frozenset({'Aperiodic', 'Associative', 'HTrivial'})

            In the finite case, the two notions coincide::

                sage: Semigroups().Aperiodic().Finite() is Semigroups().HTrivial().Finite()
                True

            TESTS::

                sage: C = Monoids().Aperiodic().Finite()
                sage: TestSuite(C).run()

            .. SEEALSO::

                - :wikipedia:`Aperiodic_semigroup`
                - :class:`Semigroups.SubcategoryMethods.RTrivial`
                - :class:`Semigroups.SubcategoryMethods.LTrivial`
                - :class:`Semigroups.SubcategoryMethods.JTrivial`
                - :class:`Semigroups.SubcategoryMethods.Aperiodic`

            TESTS::

                sage: TestSuite(C).run()
                sage: Rings().Aperiodic.__module__
                'sage.categories.semigroups'
            """
            return self._with_axiom('Aperiodic')

    Finite = LazyImport('sage.categories.finite_semigroups', 'FiniteSemigroups', at_startup=True)
    FinitelyGeneratedAsMagma = LazyImport('sage.categories.finitely_generated_semigroups', 'FinitelyGeneratedSemigroups')
    Unital = LazyImport('sage.categories.monoids', 'Monoids', at_startup=True)
    LTrivial = LazyImport('sage.categories.l_trivial_semigroups', 'LTrivialSemigroups')
    RTrivial = LazyImport('sage.categories.r_trivial_semigroups', 'RTrivialSemigroups')
    JTrivial = LazyImport('sage.categories.j_trivial_semigroups', 'JTrivialSemigroups')
    HTrivial = LazyImport('sage.categories.h_trivial_semigroups', 'HTrivialSemigroups')
    Aperiodic = LazyImport('sage.categories.aperiodic_semigroups', 'AperiodicSemigroups')

    #######################################
    class Subquotients(SubquotientsCategory):
        r"""
        The category of subquotient semi-groups.

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
            Returns an example of subquotient of a semigroup, as per
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
            Return an example of quotient of a semigroup, as per
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
                Return semigroup generators for ``self`` by
                retracting the semigroup generators of the ambient
                semigroup.

                EXAMPLES::

                    sage: S = FiniteSemigroups().Quotients().example().semigroup_generators() # todo: not implemented
                """
                return self.ambient().semigroup_generators().map(self.retract)

    class CartesianProducts(CartesianProductsCategory):

        def extra_super_categories(self):
            """
            Implement the fact that a Cartesian product of semigroups is a
            semigroup.

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
            Implement the fact that the algebra of a semigroup is indeed
            a (not necessarily unital) algebra.

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

                    sage: M = FiniteSemigroups().example(); M
                    An example of a finite semigroup:
                    the left regular band generated by ('a', 'b', 'c', 'd')
                    sage: M.semigroup_generators()
                    Family ('a', 'b', 'c', 'd')
                    sage: M.algebra(ZZ).algebra_generators()
                    Finite family {0: B['a'], 1: B['b'], 2: B['c'], 3: B['d']}
                """
                return self.basis().keys().semigroup_generators().map(self.monomial)

            # Once there will be some guarantee on the consistency between
            # gens / monoid/group/*_generators, these methods could possibly
            # be removed in favor of aliases gens -> xxx_generators in
            # the Algebras.FinitelyGenerated hierachy
            def gens(self):
                r"""
                Return the generators of ``self``.

                EXAMPLES::

                    sage: a, b = SL2Z.algebra(ZZ).gens(); a, b
                    ([ 0 -1]
                     [ 1  0],
                     [1 1]
                     [0 1])
                    sage: 2*a + b
                    2*[ 0 -1]
                      [ 1  0]
                    +
                    [1 1]
                    [0 1]
                """
                return tuple(self.monomial(g) for g in self.basis().keys().gens())

            def ngens(self):
                r"""
                Return the number of generators of ``self``.

                EXAMPLES::

                    sage: SL2Z.algebra(ZZ).ngens()
                    2
                    sage: DihedralGroup(4).algebra(RR).ngens()
                    2
                """
                return self.basis().keys().ngens()

            def gen(self, i=0):
                r"""
                Return the ``i``-th generator of ``self``.

                EXAMPLES::

                    sage: A = GL(3, GF(7)).algebra(ZZ)
                    sage: A.gen(0)
                    [3 0 0]
                    [0 1 0]
                    [0 0 1]
                """
                return self.monomial(self.basis().keys().gen(i))

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

            def trivial_representation(self, side="twosided"):
                """
                Return the trivial representation of ``self``.

                INPUT:

                - ``side`` -- ignored

                EXAMPLES::

                    sage: G = groups.permutation.Dihedral(4)
                    sage: A = G.algebra(QQ)
                    sage: V = A.trivial_representation()
                    sage: V == G.trivial_representation(QQ)
                    True
                """
                S = self.basis().keys()
                return S.trivial_representation(self.base_ring())

            def regular_representation(self, side="left"):
                """
                Return the regular representation of ``self``.

                INPUT:

                - ``side`` -- (default: ``"left"``) whether this is the
                  ``"left"`` or ``"right"`` regular representation

                EXAMPLES::

                    sage: G = groups.permutation.Dihedral(4)
                    sage: A = G.algebra(QQ)
                    sage: V = A.regular_representation()
                    sage: V == G.regular_representation(QQ)
                    True
                """
                S = self.basis().keys()
                return S.regular_representation(self.base_ring(), side)

