r"""
Finite Coxeter Groups
"""
#*****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#  Copyright (C) 2009    Nicolas Borie <nicolas dot borie at math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_with_axiom import CategoryWithAxiom, axiom
from sage.categories.coxeter_groups import CoxeterGroups
from sage.categories.complex_reflection_groups import ComplexReflectionGroups

class FiniteCoxeterGroups(CategoryWithAxiom):
    r"""
    The category of finite Coxeter groups.

    EXAMPLES::

        sage: FiniteCoxeterGroups()
        Category of finite coxeter groups
        sage: FiniteCoxeterGroups().super_categories()
        [Category of coxeter groups,
         Category of finite groups,
         Category of finite finitely generated semigroups]

        sage: G = FiniteCoxeterGroups().example()
        sage: G.cayley_graph(side = "right").plot()
        Graphics object consisting of 40 graphics primitives

    Here are some further examples::

        sage: FiniteWeylGroups().example()
        The symmetric group on {0, ..., 3}

        sage: WeylGroup(["B", 3])
        Weyl Group of type ['B', 3] (as a matrix group acting on the ambient space)

    Those other examples will eventually be also in this category::

        sage: SymmetricGroup(4)
        Symmetric group of order 4! as a permutation group
        sage: DihedralGroup(5)
        Dihedral group of order 10 as a permutation group
    """
    _base_category_class_and_axiom = (CoxeterGroups, "Finite")

    @cached_method
    def extra_super_categories(self):
        r"""
        EXAMPLES::

            sage: ComplexReflectionGroups().super_categories()
            [Category of groups]
        """
        return [ComplexReflectionGroups().Finite().WellGenerated()]

    class SubcategoryMethods:
        Irreducible = axiom("Irreducible")

    class ParentMethods:
        """
        Ambiguity resolution: the implementation of ``some_elements``
        is preferable to that of :class:`FiniteGroups`. The same holds
        for ``__iter__``, although a breath first search would be more
        natural; at least this maintains backward compatibility after
        :trac:`13589`.

        TESTS::

            sage: W = FiniteCoxeterGroups().example(3)

            sage: W.some_elements.__module__
            'sage.categories.coxeter_groups'
            sage: W.__iter__.__module__
            'sage.categories.coxeter_groups'

            sage: W.some_elements()
            [(1,), (2,), (), (1, 2)]
            sage: list(W)
            [(), (1,), (1, 2), (1, 2, 1), (2,), (2, 1)]
        """
        some_elements = CoxeterGroups.ParentMethods.__dict__["some_elements"]
        __iter__      = CoxeterGroups.ParentMethods.__dict__["__iter__"]

        @lazy_attribute
        def w0(self):
            r"""
            Return the longest element of ``self``.

            This attribute is deprecated.

            EXAMPLES::

                sage: D8 = FiniteCoxeterGroups().example(8)
                sage: D8.w0
                (1, 2, 1, 2, 1, 2, 1, 2)
                sage: D3 = FiniteCoxeterGroups().example(3)
                sage: D3.w0
                (1, 2, 1)
            """
            return self.long_element()

        def long_element(self, index_set = None):
            r"""

            INPUT:

            - ``index_set`` - a subset (as a list or iterable) of the
              nodes of the Dynkin diagram; (default: all of them)

            Returns the longest element of ``self``, or of the
            parabolic subgroup corresponding to the given ``index_set``.

            Should this method be called maximal_element? longest_element?

            EXAMPLES::

                sage: D10 = FiniteCoxeterGroups().example(10)
                sage: D10.long_element()
                (1, 2, 1, 2, 1, 2, 1, 2, 1, 2)
                sage: D10.long_element([1])
                (1,)
                sage: D10.long_element([2])
                (2,)
                sage: D10.long_element([])
                ()

                sage: D7 = FiniteCoxeterGroups().example(7)
                sage: D7.long_element()
                (1, 2, 1, 2, 1, 2, 1)

            """
            if index_set is None:
                index_set = self.index_set()
            w = self.one()
            while True:
                i = w.first_descent(index_set = index_set, positive = True)
                if i is None:
                    return w
                else:
                    w = w.apply_simple_reflection(i)

        @cached_method
        def bruhat_poset(self, facade = False):
            """
            Returns the Bruhat poset of ``self``.

            EXAMPLES::

                sage: W = WeylGroup(["A", 2])
                sage: P = W.bruhat_poset()
                sage: P
                Finite poset containing 6 elements
                sage: P.show()

            Here are some typical operations on this poset::

                sage: W = WeylGroup(["A", 3])
                sage: P = W.bruhat_poset()
                sage: u = W.from_reduced_word([3,1])
                sage: v = W.from_reduced_word([3,2,1,2,3])
                sage: P(u) <= P(v)
                True
                sage: len(P.interval(P(u), P(v)))
                10
                sage: P.is_join_semilattice()
                False

            By default, the elements of `P` are aware that they belong
            to `P`::

                sage: P.an_element().parent()
                Finite poset containing 24 elements

            If instead one wants the elements to be plain elements of
            the Coxeter group, one can use the ``facade`` option::

                sage: P = W.bruhat_poset(facade = True)
                sage: P.an_element().parent()
                Weyl Group of type ['A', 3] (as a matrix group acting on the ambient space)

            .. see also:: :func:`Poset` for more on posets and facade posets.

            TESTS::

                sage: [len(WeylGroup(["A", n]).bruhat_poset().cover_relations()) for n in [1,2,3]]
                [1, 8, 58]

            .. todo::

                - Use the symmetric group in the examples (for nicer
                  output), and print the edges for a stronger test.
                - The constructed poset should be lazy, in order to
                  handle large / infinite Coxeter groups.
            """
            from sage.combinat.posets.posets import Poset
            covers = tuple([u, v] for v in self for u in v.bruhat_lower_covers() )
            return Poset((self, covers), cover_relations = True, facade=facade)

        @cached_method
        def weak_poset(self, side = "right", facade = False):
            """
            INPUT:

            - ``side`` -- "left", "right", or "twosided" (default: "right")
            - ``facade`` -- a boolean (default: ``False``)

            Returns the left (resp. right) poset for weak order.  In
            this poset, `u` is smaller than `v` if some reduced word
            of `u` is a right (resp. left) factor of some reduced word
            of `v`.

            EXAMPLES::

                sage: W = WeylGroup(["A", 2])
                sage: P = W.weak_poset()
                sage: P
                Finite lattice containing 6 elements
                sage: P.show()

            This poset is in fact a lattice::

                sage: W = WeylGroup(["B", 3])
                sage: P = W.weak_poset(side = "left")
                sage: P.is_lattice()
                True

            so this method has an alias :meth:`weak_lattice`::

                sage: W.weak_lattice(side = "left") is W.weak_poset(side = "left")
                True

            As a bonus feature, one can create the left-right weak
            poset::

                sage: W = WeylGroup(["A",2])
                sage: P = W.weak_poset(side = "twosided")
                sage: P.show()
                sage: len(P.hasse_diagram().edges())
                8

            This is the transitive closure of the union of left and
            right order. In this poset, `u` is smaller than `v` if
            some reduced word of `u` is a factor of some reduced word
            of `v`. Note that this is not a lattice::

                sage: P.is_lattice()
                False

            By default, the elements of `P` are aware of that they
            belong to `P`::

                sage: P.an_element().parent()
                Finite poset containing 6 elements

            If instead one wants the elements to be plain elements of
            the Coxeter group, one can use the ``facade`` option::

                sage: P = W.weak_poset(facade = True)
                sage: P.an_element().parent()
                Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)

            .. see also:: :func:`Poset` for more on posets and facade posets.

            TESTS::

                sage: [len(WeylGroup(["A", n]).weak_poset(side = "right").cover_relations()) for n in [1,2,3]]
                [1, 6, 36]
                sage: [len(WeylGroup(["A", n]).weak_poset(side = "left" ).cover_relations()) for n in [1,2,3]]
                [1, 6, 36]

            .. todo::

                - Use the symmetric group in the examples (for nicer
                  output), and print the edges for a stronger test.
                - The constructed poset should be lazy, in order to
                  handle large / infinite Coxeter groups.
            """
            from sage.combinat.posets.posets import Poset
            from sage.combinat.posets.lattices import LatticePoset
            if side == "twosided":
                covers = tuple([u, v] for u in self for v in u.upper_covers(side="left")+u.upper_covers(side="right") )
                return Poset((self, covers), cover_relations = True, facade = facade)
            else:
                covers = tuple([u, v] for u in self for v in u.upper_covers(side=side) )
                return LatticePoset((self, covers), cover_relations = True, facade = facade)

        weak_lattice = weak_poset

        @cached_method
        def cambrian_lattice(self, c):
            """
            Return the c-Cambrian lattice on delta sequences (see arXiv:1503.00710 and arXiv:math/0611106).
            Delta sequences are certain 2-colored minimal factorizations of ``c`` into reflections.

            INPUT:

            - ``c`` -- a standard Coxeter element in ``self`` (as a tuple, or as an element of ``self``)

            EXAMPLES::

                sage: CoxeterGroup(["A", 2]).cambrian_lattice((1,2))
                Finite lattice containing 5 elements

                sage: CoxeterGroup(["B", 2]).cambrian_lattice((1,2))
                Finite lattice containing 6 elements

                sage: CoxeterGroup(["G", 2]).cambrian_lattice((1,2))
                Finite lattice containing 8 elements

            """
            from sage.combinat.posets.lattices import LatticePoset
            return self.m_cambrian_lattice(c,1)

        def inversion_sequence(self, word):
            """
            Return the inversion sequence corresponding to ``word``.  If
            ``word``=`[w_0,w_1,...w_k]`, then the output is `[w_0,w_0w_1w_0,\ldots,w_0w_1\cdots w_k \cdots w_1 w_0]`.

            INPUT:

            - ``word`` -- a word in the simple generators of ``self``

            EXAMPLES::

                sage: CoxeterGroup(["A", 2]).inversion_sequence([1,2,1])
                [
                [-1  1]  [ 0 -1]  [ 1  0]
                [ 0  1], [-1  0], [ 1 -1]
                ]

                sage: [t.reduced_word() for t in CoxeterGroup(["A",3]).inversion_sequence([2,1,3,2,1,3])]
                [[2], [1, 2, 1], [2, 3, 2], [1, 2, 3, 2, 1], [3], [1]]

            """
            return [self.from_reduced_word(word[:i+1]+list(reversed(word[:i]))) for i in range(len(word))]

        def reflections_from_w0(self):
            """
            Return the reflections of ``self`` using the inversion set of ``w_0``.

            EXAMPLES::

                sage: WeylGroup(['A',2]).reflections_from_w0()
                [
                [0 1 0]  [0 0 1]  [1 0 0]
                [1 0 0]  [0 1 0]  [0 0 1]
                [0 0 1], [1 0 0], [0 1 0]
                ]

                sage: WeylGroup(['A',3]).reflections_from_w0()
                [
                [-1  1  0]  [ 0 -1  1]  [ 1  0  0]  [ 0  0 -1]  [ 1  0  0]  [ 1  0  0]
                [ 0  1  0]  [-1  0  1]  [ 1 -1  1]  [-1  1 -1]  [ 1  0 -1]  [ 0  1  0]
                [ 0  0  1], [ 0  0  1], [ 0  0  1], [-1  0  0], [ 1 -1  0], [ 0  1 -1]
                ]

            """
            return self.long_element().inversions_as_reflections()

        @cached_method
        def m_cambrian_lattice(self, c, m = 1):
            """
            Return the m-Cambrian lattice on ``m``-delta sequences (see arXiv:1503.00710 and arXiv:math/0611106).
            ``m``-delta sequences are certain ``m``-colored minimal factorizations of ``c`` into reflections.

            INPUT:

            - ``c`` -- a Coxeter element of ``self`` (as a tuple, or as an element of ``self``)
            - ``m`` -- a positive integer (default: 1)

            EXAMPLES::

                sage: CoxeterGroup(["A",2]).m_cambrian_lattice((1,2))
                Finite lattice containing 5 elements

                sage: CoxeterGroup(["A",2]).m_cambrian_lattice((1,2),2)
                Finite lattice containing 12 elements

            """
            from sage.combinat.posets.posets import Poset
            from sage.combinat.posets.lattices import LatticePoset
            if hasattr(c,"reduced_word"):
               c = c.reduced_word()
            elif not isinstance(c,list):
               c = list(c)
            inv_woc = self.inversion_sequence(self.long_element().coxeter_sorting_word(c))
            S = self.simple_reflections()
            T = self.reflections_from_w0()
            Twords = {t:t.reduced_word() for t in T}#PhiP=T
            id = sorted([[s,0] for s in S])
            elements = []
            covers = []
            new = [id]
            while new != []:
                for new_element in new:
                    new.remove(new_element)
                    elements.append(new_element)
                    for t in new_element:
                        if t[1]<m:
                            cov_element = [s for s in new_element if s!=t]
                            cov_element.append([t[0],t[1]+1])
                            for t_conj in [[i,t[1]] for i in inv_woc[inv_woc.index(t[0]):]]+[[i,t[1]+1] for i in inv_woc[:inv_woc.index(t[0])]]:
                                if t_conj in cov_element:
                                    cov_element.remove(t_conj)
                                    tmp = t[0]*t_conj[0]*t[0]
                                    invs = self.inversion_sequence(Twords[t[0]]+Twords[t_conj[0]])
                                    plus_or_minus = invs.count(tmp)
                                    if plus_or_minus % 2 == 1:
                                        cov_element.append([tmp,t_conj[1]])
                                    else:
                                        cov_element.append([tmp,t_conj[1]-1])
                            cov_element = sorted(cov_element)
                            if cov_element not in elements and cov_element not in new:
                                new.append(cov_element)
                            covers.append([tuple(map(tuple,new_element)),tuple(map(tuple,cov_element))])
            return LatticePoset([[tuple(map(tuple,e)) for e in elements],covers],cover_relations=True)

    class ElementMethods:

        @cached_in_parent_method
        def bruhat_upper_covers(self):
            r"""
            Returns all the elements that cover ``self`` in Bruhat order.

            EXAMPLES::

                sage: W = WeylGroup(["A",4])
                sage: w = W.from_reduced_word([3,2])
                sage: print([v.reduced_word() for v in w.bruhat_upper_covers()])
                [[4, 3, 2], [3, 4, 2], [2, 3, 2], [3, 1, 2], [3, 2, 1]]

                sage: W = WeylGroup(["B",6])
                sage: w = W.from_reduced_word([1,2,1,4,5])
                sage: C = w.bruhat_upper_covers()
                sage: len(C)
                9
                sage: print([v.reduced_word() for v in C])
                [[6, 4, 5, 1, 2, 1], [4, 5, 6, 1, 2, 1], [3, 4, 5, 1, 2, 1], [2, 3, 4, 5, 1, 2],
                [1, 2, 3, 4, 5, 1], [4, 5, 4, 1, 2, 1], [4, 5, 3, 1, 2, 1], [4, 5, 2, 3, 1, 2],
                [4, 5, 1, 2, 3, 1]]
                sage: ww = W.from_reduced_word([5,6,5])
                sage: CC = ww.bruhat_upper_covers()
                sage: print([v.reduced_word() for v in CC])
                [[6, 5, 6, 5], [4, 5, 6, 5], [5, 6, 4, 5], [5, 6, 5, 4], [5, 6, 5, 3], [5, 6, 5, 2],
                [5, 6, 5, 1]]

            Recursive algorithm: write `w` for ``self``. If `i` is a
            non-descent of `w`, then the covers of `w` are exactly
            `\{ws_i, u_1s_i, u_2s_i,..., u_js_i\}`, where the `u_k`
            are those covers of `ws_i` that have a descent at `i`.
            """

            i = self.first_descent(positive=True)
            if i is not None:
                wsi = self.apply_simple_reflection(i)
                return [u.apply_simple_reflection(i) for u in wsi.bruhat_upper_covers() if u.has_descent(i)] + [wsi]
            else:
                return []

        def coxeter_knuth_neighbor(self, w):
            r"""
            Return the Coxeter-Knuth (oriented) neighbors of the reduced word `w` of ``self``.

            INPUT:

            - ``w`` -- reduced word of ``self``

            The Coxeter-Knuth relations are given by `a a+1 a \sim a+1 a a+1`, `abc \sim acb`
            if `b<a<c` and `abc \sim bac` if `a<c<b`. This method returns all neighbors of
            ``w`` under the Coxeter-Knuth relations oriented from left to right.

            EXAMPLES::

                sage: W = WeylGroup(['A',4], prefix='s')
                sage: word = [1,2,1,3,2]
                sage: w = W.from_reduced_word(word)
                sage: w.coxeter_knuth_neighbor(word)
                {(1, 2, 3, 1, 2), (2, 1, 2, 3, 2)}

                sage: word = [1,2,1,3,2,4,3]
                sage: w = W.from_reduced_word(word)
                sage: w.coxeter_knuth_neighbor(word)
                {(1, 2, 1, 3, 4, 2, 3), (1, 2, 3, 1, 2, 4, 3), (2, 1, 2, 3, 2, 4, 3)}

            TESTS::

                sage: W = WeylGroup(['B',4], prefix='s')
                sage: word = [1,2]
                sage: w = W.from_reduced_word(word)
                sage: w.coxeter_knuth_neighbor(word)
                Traceback (most recent call last):
                ...
                NotImplementedError: This has only been implemented in finite type A so far!
            """
            C = self.parent().cartan_type()
            if not C[0] == 'A':
                raise NotImplementedError("This has only been implemented in finite type A so far!")
            d = []
            for i in range(2,len(w)):
                v = [j for j in w]
                if w[i-2] == w[i]:
                    if w[i] == w[i-1] - 1:
                        v[i-2] = w[i-1]
                        v[i] = w[i-1]
                        v[i-1] = w[i]
                        d += [tuple(v)]
                elif w[i-1]<w[i-2] and w[i-2]<w[i]:
                    v[i] = w[i-1]
                    v[i-1] = w[i]
                    d += [tuple(v)]
                elif w[i-2]<w[i] and w[i]<w[i-1]:
                    v[i-2] = w[i-1]
                    v[i-1] = w[i-2]
                    d += [tuple(v)]
            return set(d)

        def coxeter_knuth_graph(self):
            r"""
            Return the Coxeter-Knuth graph of type `A`.

            The Coxeter-Knuth graph of type `A` is generated by the Coxeter-Knuth relations which are
            given by `a a+1 a \sim a+1 a a+1`, `abc \sim acb` if `b<a<c` and `abc \sim bac` if `a<c<b`.

            EXAMPLES::

                sage: W = WeylGroup(['A',4], prefix='s')
                sage: w = W.from_reduced_word([1,2,1,3,2])
                sage: D = w.coxeter_knuth_graph()
                sage: D.vertices()
                [(1, 2, 1, 3, 2),
                (1, 2, 3, 1, 2),
                (2, 1, 2, 3, 2),
                (2, 1, 3, 2, 3),
                (2, 3, 1, 2, 3)]
                sage: D.edges()
                [((1, 2, 1, 3, 2), (1, 2, 3, 1, 2), None),
                ((1, 2, 1, 3, 2), (2, 1, 2, 3, 2), None),
                ((2, 1, 2, 3, 2), (2, 1, 3, 2, 3), None),
                ((2, 1, 3, 2, 3), (2, 3, 1, 2, 3), None)]

                sage: w = W.from_reduced_word([1,3])
                sage: D = w.coxeter_knuth_graph()
                sage: D.vertices()
                [(1, 3), (3, 1)]
                sage: D.edges()
                []

            TESTS::

                sage: W = WeylGroup(['B',4], prefix='s')
                sage: w = W.from_reduced_word([1,2])
                sage: w.coxeter_knuth_graph()
                Traceback (most recent call last):
                ...
                NotImplementedError: This has only been implemented in finite type A so far!
            """
            from sage.graphs.all import Graph
            R = [tuple(v) for v in self.reduced_words()]
            G = Graph()
            G.add_vertices(R)
            G.add_edges([v,vp] for v in R for vp in self.coxeter_knuth_neighbor(v))
            return G

    class Irreducible(CategoryWithAxiom):
        r"""
        The category of finite irreducible Coxeter groups.
        """
        class ParentMethods:
            pass

        class ElementMethods:
            pass
