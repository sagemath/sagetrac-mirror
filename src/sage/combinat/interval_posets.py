r"""
Tamari Interval-posets

This module implements Tamari interval-posets: combinatorial objects which
represent intervals of the Tamari order. They have been introduced in [PCh2013]_
and allow for many combinatorial operations on Tamari intervals. In particular,
they are linked to :class:`DyckWords` and :class:`BinaryTrees`. An introduction
into Tamari interval-posets is given in Chapter 7 of [Pons2013]_.

The Tamari lattice can be defined as a lattice structure on either of several
classes of Catalan objects, especially binary trees and Dyck paths
[TamBrack1962]_ [HuTAss1972]_ [StaEC2]_ . An interval can be seen as a pair of
comparable elements. The number of intervals has been given in
[ChapTamari08]_ .

REFERENCES:

.. [PCh2013] Gregory Chatel, Viviane Pons,
   *Counting smaller trees in the Tamari order*,
   FPSAC, 2013, :arxiv:`1212.0751v1`.
.. [Pons2013] Viviane Pons,
   *Combinatoire algebrique liee aux ordres sur les permutations*,
   thesis, :arxiv:`1310.1805v1`.
.. [TamBrack1962] Dov Tamari,
   *The algebra of bracketings and their enumeration*,
   Nieuw Arch. Wisk., 1962.
.. [HuTAss1972] Samuel Huang and Dov Tamari,
   *Problems of associativity: A simple proof for the lattice property
   of systems ordered by a semi-associative law*,
   J. Combinatorial Theory Ser. A, 1972,
   http://www.sciencedirect.com/science/article/pii/0097316572900039 .
.. [StaEC2] Richard P. Stanley,
   *Enumerative combinatorics, vol. 2*,
   Cambridge University Press, 1999.
.. [ChapTamari08] Frederic Chapoton,
   *Sur le nombre d'intervalles dans les treillis de Tamari*,
   Sem. Lothar. Combin. 2008.
   :arxiv:`math/0602368v1`.

**AUTHORS:**

- Viviane Pons 2014: initial implementation
- Frederic Chapoton 2014: review
"""
#*****************************************************************************
#       Copyright (C) 2013 Viviane Pons <viviane.pons@univie.ac.at>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.finite_posets import FinitePosets
from sage.categories.finite_posets import Posets
from sage.combinat.binary_tree import BinaryTrees
from sage.combinat.binary_tree import LabelledBinaryTrees
from sage.combinat.dyck_word import DyckWords
from sage.combinat.permutation import Permutation
from sage.combinat.posets.posets import Poset
from sage.combinat.posets.posets import FinitePoset
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.integer import Integer
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.structure.element import Element
from sage.structure.global_options import GlobalOptions
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

TamariIntervalPosetOptions = GlobalOptions(name="Tamari Interval-posets",
    doc=r"""
    Set and display the global options for Tamari interval-posets. If no parameters
    are set, then the function returns a copy of the options dictionary.

    The ``options`` to Tamari interval-posets can be accessed as the method
    :obj:`TamariIntervalPosets.global_options` of :class:`TamariIntervalPosets` and
    related parent classes.
    """,
    end_doc=r"""
    EXAMPLES::

        sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
        sage: ip.latex_options()["color_decreasing"]
        'red'
        sage: TamariIntervalPosets.global_options(latex_color_decreasing='green')
        sage: ip.latex_options()["color_decreasing"]
        'green'
        sage: TamariIntervalPosets.global_options.reset()
        sage: ip.latex_options()["color_decreasing"]
        'red'

    """,
    latex_tikz_scale=dict(default=1,
                          description='The default value for the tikz scale when latexed',
                          checker=lambda x: True),  # More trouble than it's worth to check
    latex_line_width_scalar=dict(default=0.5,
                                 description='The default value for the line width as a'
                                             'multiple of the tikz scale when latexed',
                                 checker=lambda x: True),  # More trouble than it's worth to check
    latex_color_decreasing=dict(default="red",
                                description='The default color of decreasing relations when latexed',
                                checker=lambda x: True),  # More trouble than it's worth to check
    latex_color_increasing=dict(default="blue",
                                description='The default color of increasing relations when latexed',
                                checker=lambda x: True),  # More trouble than it's worth to check
    latex_hspace=dict(default=1,
                      description='The default difference between horizontal coordinates of vertices when latexed',
                      checker=lambda x: True),  # More trouble than it's worth to check
    latex_vspace=dict(default=1,
                      description='The default difference between vertical coordinates of vertices when latexed',
                      checker=lambda x: True),  # More trouble than it's worth to check
)


class TamariIntervalPoset(Element):
    r"""
    The class of Tamari interval-posets.

    An interval-poset is a labelled poset of size `n`, with labels
    `1, 2, \dots, n`, satisfying the following conditions:

    - if `a<c` (as integers) and `a` precedes `c` in the poset, then,
      for all `b` such that `a<b<c`, `b` precedes `c`,

    - if `a<c` (as integers) and `c` precedes `a` in the poset, then,
      for all `b` such that `a<b<c`, `b` precedes `a`.

    We use the word "precedes" here to distinguish the poset order and
    the natural order on numbers. "Precedes" means "is smaller than
    with respect to the poset structure"; this does not imply a
    covering relation. 

    Interval-posets of size `n` are in bijection with intervals of
    the Tamari lattice of binary trees of size `n`. Specifically, if
    `P` is an interval-poset of size `n`, then the set of linear
    extensions of `P` (as permutations in `S_n`) is an interval in the
    right weak order (see
    :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`),
    and is in fact the preimage of an interval in the Tamari lattice (of
    binary trees of size `n`) under the operation which sends a
    permutation to its right-to-left binary search tree
    (:meth:`~sage.combinat.permutation.Permutation.binary_search_tree`
    with the ``left_to_right`` variable set to ``False``)
    without its labelling.

    INPUT:

    - ``size`` -- an integer, the size of the interval-posets (number of
      vertices)

    - ``relations`` -- a list (or tuple) of pairs ``(a,b)`` (themselves
      lists or tuples), each representing a relation of the form
      '`a` precedes `b`' in the poset.

    - ``check`` -- (default: ``True``) whether to check the interval-poset
      condition or not.

    .. WARNING::

        The ``relations`` input can be a list or tuple, but not an
        iterator (nor should its entries be iterators).

    NOTATION:

    Here and in the following, the signs `<` and `>` always refer to
    the natural ordering on integers, whereas the word "precedes" refers
    to the order of the interval-poset. "Minimal" and "maximal" refer
    to the natural ordering on integers.

    The *increasing relations* of an interval-poset `P` mean the pairs
    `(a, b)` of elements of `P` such that `a < b` as integers and `a`
    precedes `b` in `P`. The *initial forest* of `P` is the poset
    obtained by imposing (only) the increasing relations on the ground
    set of `P`. It is a sub-interval poset of `P`, and is a forest with
    its roots on top. This forest is usually given the structure of a
    planar forest by ordering brother nodes by their labels; it then has
    the property that if its nodes are traversed in post-order
    (see :meth:~sage.combinat.abstract_tree.AbstractTree.post_order_traversal`,
    and traverse the trees of the forest from left to right as well),
    then the labels encountered are `1, 2, \ldots, n` in this order.

    The *decreasing relations* of an interval-poset `P` mean the pairs
    `(a, b)` of elements of `P` such that `b < a` as integers and `a`
    precedes `b` in `P`. The *final forest* of `P` is the poset
    obtained by imposing (only) the decreasing relations on the ground
    set of `P`. It is a sub-interval poset of `P`, and is a forest with
    its roots on top. This forest is usually given the structure of a
    planar forest by ordering brother nodes by their labels; it then has
    the property that if its nodes are traversed in pre-order
    (see :meth:~sage.combinat.abstract_tree.AbstractTree.pre_order_traversal`,
    and traverse the trees of the forest from left to right as well),
    then the labels encountered are `1, 2, \ldots, n` in this order.

    EXAMPLES::

        sage: TamariIntervalPoset(0,[])
        The tamari interval of size 0 induced by relations []
        sage: TamariIntervalPoset(3,[])
        The tamari interval of size 3 induced by relations []
        sage: TamariIntervalPoset(3,[(1,2)])
        The tamari interval of size 3 induced by relations [(1, 2)]
        sage: TamariIntervalPoset(3,[(1,2),(2,3)])
        The tamari interval of size 3 induced by relations [(1, 2), (2, 3)]
        sage: TamariIntervalPoset(3,[(1,2),(2,3),(1,3)])
        The tamari interval of size 3 induced by relations [(1, 2), (2, 3)]
        sage: TamariIntervalPoset(3,[(1,2),(3,2)])
        The tamari interval of size 3 induced by relations [(1, 2), (3, 2)]
        sage: TamariIntervalPoset(3,[[1,2],[2,3]])
        The tamari interval of size 3 induced by relations [(1, 2), (2, 3)]
        sage: TamariIntervalPoset(3,[[1,2],[2,3],[1,2],[1,3]])
        The tamari interval of size 3 induced by relations [(1, 2), (2, 3)]

        sage: TamariIntervalPoset(3,[(3,4)])
        Traceback (most recent call last):
        ...
        ValueError: The relations do not correspond to the size of the poset.

        sage: TamariIntervalPoset(2,[(2,1),(1,2)])
        Traceback (most recent call last):
        ...
        ValueError: Hasse diagram contains cycles.

        sage: TamariIntervalPoset(3,[(1,3)])
        Traceback (most recent call last):
        ...
        ValueError: This does not satisfy the Tamari interval-poset condition.

    It is also possible to transform a poset directly into an interval-poset::

        sage: TIP = TamariIntervalPosets()
        sage: p = Poset( ([1,2,3], [(1,2)]))
        sage: TIP(p)
        The tamari interval of size 3 induced by relations [(1, 2)]
        sage: TIP(Poset({1: []}))
        The tamari interval of size 1 induced by relations []
        sage: TIP(Poset({}))
        The tamari interval of size 0 induced by relations []
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        Ensure that interval-posets created by the enumerated sets and directly
        are the same and that they are instances of :class:`TamariIntervalPoset`.

        TESTS::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.parent()
            Interval-posets
            sage: type(ip)
            <class 'sage.combinat.interval_posets.TamariIntervalPosets_all_with_category.element_class'>

            sage: ip2 = TamariIntervalPosets()(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip2.parent() is ip.parent()
            True
            sage: type(ip) is type(ip2)
            True

            sage: ip3 = TamariIntervalPosets(4)([(2,4),(3,4),(2,1),(3,1)])
            sage: ip3.parent() is ip.parent()
            True
            sage: type(ip3) is type(ip)
            True
        """
        return TamariIntervalPosets_all().element_class(*args, **opts)

    def __init__(self, size, relations, check=True):
        r"""
        TESTS::

            sage: TamariIntervalPoset(3,[(1,2),(3,2)]).parent()
            Interval-posets
        """
        parent = TamariIntervalPosets()
        self._size = size
        self._poset = Poset((range(1, size + 1), relations))
        if self._poset.cardinality() != size:
            # This can happen as the Poset constructor automatically adds
            # in elements from the relations.
            raise ValueError("The relations do not correspond to the size of the poset.")

        if check and not TamariIntervalPosets.check_poset(self._poset):
            raise ValueError("This does not satisfy the Tamari interval-poset condition.")

        Element.__init__(self, parent)

        self._cover_relations = tuple(self._poset.cover_relations())

        self._latex_options = dict()

    def set_latex_options(self, D):
        r"""
        Set the latex options for use in the ``_latex_`` function.  The
        default values are set in the ``__init__`` function.

        - ``tikz_scale`` -- (default: 1) scale for use with the tikz package.

        - ``line_width`` -- (default: 1*``tikz_scale``) value representing the
          line width.

        - ``color_decreasing`` -- (default: red) the color for decreasing relations.

        - ``color_increasing`` -- (default: blue) the color for increasing relations.

        - ``hspace`` -- (default: 1) the difference between horizontal coordinates of adjacent vertices

        - ``vspace`` -- (default: 1) the difference between vertical coordinates of adjacent vertices

        INPUT:

        - ``D`` -- a dictionary with a list of latex parameters to change

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.latex_options()["color_decreasing"]
            'red'
            sage: ip.set_latex_options({"color_decreasing":'green'})
            sage: ip.latex_options()["color_decreasing"]
            'green'
            sage: ip.set_latex_options({"color_increasing":'black'})
            sage: ip.latex_options()["color_increasing"]
            'black'

        To change the options for all interval-posets, use the parent's
        global latex options::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.latex_options()["color_decreasing"]
            'red'
            sage: ip2.latex_options()["color_decreasing"]
            'red'
            sage: TamariIntervalPosets.global_options(latex_color_decreasing='green')
            sage: ip.latex_options()["color_decreasing"]
            'green'
            sage: ip2.latex_options()["color_decreasing"]
            'green'
            sage: TamariIntervalPosets.global_options.reset()
        """
        for opt in D:
            self._latex_options[opt] = D[opt]

    def latex_options(self):
        r"""
        Return the latex options for use in the ``_latex_`` function as a
        dictionary. The default values are set using the global options.

        - ``tikz_scale`` -- (default: 1) scale for use with the tikz package.

        - ``line_width`` -- (default: 1*``tikz_scale``) value representing the
          line width.

        - ``color_decreasing`` -- (default: red) the color for decreasing relations.

        - ``color_increasing`` -- (default: blue) the color for increasing relations.

        - ``hspace`` -- (default: 1) the difference between horizontal coordinates of adjacent vertices

        - ``vspace`` -- (default: 1) the difference between vertical coordinates of adjacent vertices

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.latex_options()['color_decreasing']
            'red'
            sage: ip.latex_options()['hspace']
            1
        """
        d = self._latex_options.copy()
        if "tikz_scale" not in d:
            d["tikz_scale"] = self.parent().global_options["latex_tikz_scale"]
        if "line_width" not in d:
            d["line_width"] = self.parent().global_options["latex_line_width_scalar"] * d["tikz_scale"]
        if "color_decreasing" not in d:
            d["color_decreasing"] = self.parent().global_options["latex_color_decreasing"]
        if "color_increasing" not in d:
            d["color_increasing"] = self.parent().global_options["latex_color_increasing"]
        if "hspace" not in d:
            d["hspace"] = self.parent().global_options["latex_hspace"]
        if "vspace" not in d:
            d["vspace"] = self.parent().global_options["latex_vspace"]
        return d

    def _latex_(self):
        r"""
        A latex representation of ``self`` using the tikzpicture package.

        If `x` precedes `y`, then `y` will always be placed on top of `x`
        and/or to the right of `x`.
        Decreasing relations are drawn vertically and increasing relations
        horizontally.
        The algorithm tries to avoid superposition but on big
        interval-posets, it might happen.

        You can use ``self.set_latex_options()`` to change default latex
        options. Or you can use the parent's global options.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: print ip._latex_()
            \begin{tikzpicture}[scale=1]
            \node(T2) at (0,-1) {2};
            \node(T3) at (1,-2) {3};
            \node(T1) at (1,0) {1};
            \node(T4) at (2,-1) {4};
            \draw[line width = 0.5, color=blue] (T2) -- (T4);
            \draw[line width = 0.5, color=red] (T2) -- (T1);
            \draw[line width = 0.5, color=blue] (T3) -- (T4);
            \draw[line width = 0.5, color=red] (T3) -- (T1);
            \end{tikzpicture}
        """
        latex.add_package_to_preamble_if_available("tikz")
        latex_options = self.latex_options()
        start = "\\begin{tikzpicture}[scale=" + str(latex_options['tikz_scale']) + "]\n"
        end = "\\end{tikzpicture}"

        def draw_node(j, x, y):
            r"""
            Internal method to draw vertices
            """
            return "\\node(T" + str(j) + ") at (" + str(x) + "," + str(y) + ") {" + str(j) + "};\n"

        def draw_increasing(i, j):
            r"""
            Internal method to draw increasing relations
            """
            return "\\draw[line width = " + str(latex_options["line_width"]) + ", color=" + latex_options["color_increasing"] + "] (T" + str(i) + ") -- (T" + str(j) + ");\n"

        def draw_decreasing(i, j):
            r"""
            Internal method to draw decreasing relations
            """
            return "\\draw[line width = " + str(latex_options["line_width"]) + ", color=" + latex_options["color_decreasing"] + "] (T" + str(i) + ") -- (T" + str(j) + ");\n"
        if self.size() == 0:
            nodes = "\\node(T0) at (0,0){$\emptyset$};"
            relations = ""
        else:
            nodes = ""  # latex for node decraltions
            relations = ""  # latex for drawing relations
            to_draw = []
            to_draw.append((1, 0))  # a pilo of nodes to be drawn with their y position

            current_parent = [self.increasing_parent(1)]  # a pilo for the current increasing parents
            parenty = [0]  # a pilo for the current parent y positions
            if current_parent[-1] is not None:
                relations += draw_increasing(1, current_parent[-1])
            vspace = latex_options["vspace"]
            hspace = latex_options["hspace"]
            x = 0
            y = 0

            # the idea is that we draw the nodes from left to right and save their y position
            for i in xrange(2, self.size() + 1):
                # at each, we draw all possible nodes and add the current node to the to_draw pilo
                decreasing_parent = self.decreasing_parent(i)
                increasing_parent = self.increasing_parent(i)
                while len(to_draw) > 0 and (decreasing_parent is None or decreasing_parent < to_draw[-1][0]):
                    # we draw all the nodes which can be placed at x
                    # we know these nodes won't have any more decreasing children (so their horizontal position is fixed)
                    n = to_draw.pop()
                    nodes += draw_node(n[0], x, n[1])
                if i != current_parent[-1]:
                    #i is not the current increasing parent
                    if (not self.le(i, i - 1) and decreasing_parent is not None):
                        # there is no decreasing relation between i and i-1
                        #they share a decreasing parent and are placed alongside horizontally
                        x += hspace
                        if current_parent[-1] is not None:
                            y -= vspace
                    else:
                        #otherwise, they are placed alongside vertically
                        y -= vspace
                    if increasing_parent != current_parent[-1]:
                        current_parent.append(increasing_parent)
                        parenty.append(y)
                    nodey = y
                else:
                    # i is the current increasing parent so it takes the current vertical position
                    current_parent.pop()
                    x += hspace
                    nodey = parenty.pop()
                    if len(current_parent) == 0 or increasing_parent != current_parent[-1]:
                        current_parent.append(increasing_parent)
                        parenty.append(nodey)
                to_draw.append((i, nodey))
                if increasing_parent is not None:
                    relations += draw_increasing(i, increasing_parent)
                if decreasing_parent is not None:
                    relations += draw_decreasing(i, decreasing_parent)
            for n in to_draw:
                # we draw all remaining nodes
                nodes += draw_node(n[0], x, n[1])
        return start + nodes + relations + end

    def poset(self):
        r"""
        Return ``self`` as a labelled poset.

        An interval-poset is indeed constructed from a labelled poset which 
        is stored internally. This method allows to access the poset and 
        all the associated methods.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(3,2),(2,4),(3,4)])
            sage: pos = ip.poset(); pos
            Finite poset containing 4 elements
            sage: pos.maximal_chains()
            [[3, 2, 4], [1, 2, 4]]
            sage: pos.maximal_elements()
            [4]
            sage: pos.is_lattice()
            False
        """
        return self._poset

    @cached_method
    def increasing_cover_relations(self):
        r"""
        Return the cover relations of the initial forest of ``self``
        (the poset formed by keeping only the relations of the form
        `a` precedes `b` with `a<b`).

        The initial forest of ``self`` is a forest with its roots
        being on top. It is also called the increasing poset of ``self``.

        .. WARNING::

            This method computes the cover relations of the initial
            forest. This is not identical with the cover relations of
            ``self`` which happen to be increasing!

        .. SEEALSO::

            :meth:`initial_forest`.

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(1,2),(3,2),(2,4),(3,4)]).increasing_cover_relations()
            [(1, 2), (2, 4), (3, 4)]
            sage: TamariIntervalPoset(3,[(1,2),(1,3),(2,3)]).increasing_cover_relations()
            [(1, 2), (2, 3)]
        """
        relations = []
        size = self.size()
        for i in xrange(1, size):
            for j in xrange(i + 1, size + 1):
                if self.le(i, j):
                    relations.append((i, j))
                    break
        return relations

    def increasing_roots(self):
        r"""
        Return the root vertices of the initial forest of ``self``,
        i.e., the vertices `a` of ``self`` such that there is no
        `b>a` with `a` precedes `b`.

        OUTPUT:

        The list of all roots of the initial forest of ``self``, in
        decreasing order.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.increasing_roots()
            [6, 5, 2]
            sage: ip.initial_forest().increasing_roots()
            [6, 5, 2]
        """
        size = self.size()
        if size == 0:
            return []
        roots = [size]
        root = size
        for i in xrange(size - 1, 0, -1):
            if not self.le(i, root):
                roots.append(i)
                root = i
        return roots

    def increasing_children(self, v):
        r"""
        Return the children of ``v`` in the initial forest of ``self``.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self`` (between 1 and ``size``)

        OUTPUT:

        The list of all children of ``v`` in the initial forest of
        ``self``, in decreasing order.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.increasing_children(2)
            [1]
            sage: ip.increasing_children(5)
            [4, 3]
            sage: ip.increasing_children(1)
            []
        """
        children = []
        root = None
        for i in xrange(v - 1, 0, -1):
            if not self.le(i, v):
                break
            if root is None or not self.le(i, root):
                children.append(i)
                root = i
        return children

    def increasing_parent(self, v):
        r"""
        Return the vertex parent of ``v`` in the initial forest of ``self``.

        This is the minimal `b>v` such that `v` precedes `b`. If there
        is no such vertex (that is, `v` is an increasing root), then
        ``None`` is returned.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self`` (between 1 and ``size``)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.increasing_parent(1)
            2
            sage: ip.increasing_parent(3)
            5
            sage: ip.increasing_parent(4)
            5
            sage: ip.increasing_parent(5) is None
            True
        """
        parent = None
        for i in xrange(self.size(), v, -1):
            if self.le(v, i):
                parent = i
        return parent

    @cached_method
    def decreasing_cover_relations(self):
        r"""
        Return the cover relations of the final forest of ``self``
        (the poset formed by keeping only the relations of the form
        `a` precedes `b` with `a>b`).

        The final forest of ``self`` is a forest with its roots
        being on top. It is also called the decreasing poset of ``self``.

        .. WARNING::

            This method computes the cover relations of the final
            forest. This is not identical with the cover relations of
            ``self`` which happen to be decreasing!

        .. SEEALSO::

            :meth:`final_forest`.

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(2,1),(3,2),(3,4),(4,2)]).decreasing_cover_relations()
            [(4, 2), (3, 2), (2, 1)]
            sage: TamariIntervalPoset(4,[(2,1),(4,3),(2,3)]).decreasing_cover_relations()
            [(4, 3), (2, 1)]
            sage: TamariIntervalPoset(3,[(2,1),(3,1),(3,2)]).decreasing_cover_relations()
            [(3, 2), (2, 1)]
        """
        relations = []
        for i in xrange(self.size(), 1, -1):
            for j in xrange(i - 1, 0, -1):
                if self.le(i, j):
                    relations.append((i, j))
                    break
        return relations

    def decreasing_roots(self):
        r"""
        Return the root vertices of the final forest of ``self``,
        i.e., the vertices `b` such that there is no `a<b` with `b`
        preceding `a`.

        OUTPUT:

        The list of all roots of the final forest of ``self``, in
        increasing order.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.decreasing_roots()
            [1, 2]
            sage: ip.final_forest().decreasing_roots()
            [1, 2]
        """
        if self.size() == 0:
            return []
        roots = [1]
        root = 1
        for i in xrange(2, self.size() + 1):
            if not self.le(i, root):
                roots.append(i)
                root = i
        return roots

    def decreasing_children(self, v):
        r"""
        Return the children of ``v`` in the final forest of ``self``.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self`` (between 1 and ``size``)

        OUTPUT:

        The list of all children of ``v`` in the final forest of ``self``,
        in increasing order.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.decreasing_children(2)
            [3, 5]
            sage: ip.decreasing_children(3)
            [4]
            sage: ip.decreasing_children(1)
            []
        """
        children = []
        root = None
        for i in xrange(v + 1, self.size() + 1):
            if not self.le(i, v):
                break
            if root is None or not self.le(i, root):
                children.append(i)
                root = i
        return children

    def decreasing_parent(self, v):
        r"""
        Return the vertex parent of ``v`` in the final forest of ``self``.
        This is the maximal `a < v` such that ``v`` precedes ``a``. If
        there is no such vertex (that is, `v` is a decreasing root), then
        ``None`` is returned.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self`` (between
          1 and ``size``)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.decreasing_parent(4)
            3
            sage: ip.decreasing_parent(3)
            2
            sage: ip.decreasing_parent(5)
            2
            sage: ip.decreasing_parent(2) is None
            True
        """
        parent = None
        for i in xrange(1, v):
            if self.le(v, i):
                parent = i
        return parent

    def le(self, e1, e2):
        r"""
        Return whether ``e1`` precedes or equals ``e2`` in ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.le(1,2)
            True
            sage: ip.le(1,3)
            True
            sage: ip.le(2,3)
            True
            sage: ip.le(3,4)
            False
            sage: ip.le(1,1)
            True
        """
        return self._poset.le(e1, e2)

    def lt(self, e1, e2):
        r"""
        Return whether ``e1`` strictly precedes ``e2`` in ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.lt(1,2)
            True
            sage: ip.lt(1,3)
            True
            sage: ip.lt(2,3)
            True
            sage: ip.lt(3,4)
            False
            sage: ip.lt(1,1)
            False
        """
        return self._poset.lt(e1, e2)

    def ge(self, e1, e2):
        r"""
        Return whether ``e2`` precedes or equals ``e1`` in ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.ge(2,1)
            True
            sage: ip.ge(3,1)
            True
            sage: ip.ge(3,2)
            True
            sage: ip.ge(4,3)
            False
            sage: ip.ge(1,1)
            True
        """
        return self._poset.ge(e1, e2)

    def gt(self, e1, e2):
        r"""
        Return whether ``e2`` strictly precedes ``e1`` in ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.gt(2,1)
            True
            sage: ip.gt(3,1)
            True
            sage: ip.gt(3,2)
            True
            sage: ip.gt(4,3)
            False
            sage: ip.gt(1,1)
            False
        """
        return self._poset.gt(e1, e2)

    def size(self):
        r"""
        Return the size (number of vertices) of the interval-poset.

        EXAMPLES::

            sage: TamariIntervalPoset(3,[(2,1),(3,1)]).size()
            3
        """
        return self._size

    def complement(self):
        r"""
        Return the complement of the interval-poset ``self``.

        If `P` is a Tamari interval-poset of size `n`, then the
        *complement* of `P` is defined as the interval-poset `Q` whose
        base set is `[n] = \{1, 2, \ldots, n\}` (just as for `P`), but
        whose order relation has `a` precede `b` if and only if
        `n + 1 - a` precedes `n + 1 - b` in `P`.

        In terms of the Tamari lattice, the *complement* is the symmetric
        of ``self``. It is formed from the left-right symmeterized of
        the binary trees of the interval (switching left and right
        subtrees, see
        :meth:`~sage.combinat.binary_tree.BinaryTree.left_right_symmetry`).
        In particular, initial intervals are sent to final intervals and
        vice-versa.

        EXAMPLES::

            sage: TamariIntervalPoset(3, [(2, 1), (3, 1)]).complement()
            The tamari interval of size 3 induced by relations [(1, 3), (2, 3)]
            sage: TamariIntervalPoset(0, []).complement()
            The tamari interval of size 0 induced by relations []
            sage: ip = TamariIntervalPoset(4, [(1, 2), (2, 4), (3, 4)])
            sage: ip.complement() == TamariIntervalPoset(4, [(2, 1), (3, 1), (4, 3)])
            True
            sage: ip.lower_binary_tree() == ip.complement().upper_binary_tree().left_right_symmetry()
            True
            sage: ip.upper_binary_tree() == ip.complement().lower_binary_tree().left_right_symmetry()
            True
            sage: ip.is_initial_interval()
            True
            sage: ip.complement().is_final_interval()
            True
        """
        N = self._size + 1
        new_covers = [[N - i[0], N - i[1]] for i in self._poset.cover_relations_iterator()]
        return TamariIntervalPoset(N - 1, new_covers)

    def _repr_(self):
        r"""
        TESTS::

            sage: TamariIntervalPoset(3,[(2,1),(3,1)])
            The tamari interval of size 3 induced by relations [(3, 1), (2, 1)]
            sage: TamariIntervalPoset(3,[(3,1),(2,1)])
            The tamari interval of size 3 induced by relations [(3, 1), (2, 1)]
            sage: TamariIntervalPoset(3,[(2,3),(2,1)])
            The tamari interval of size 3 induced by relations [(2, 3), (2, 1)]
        """
        return "The tamari interval of size %s induced by relations %s" % (self.size(), str(self.increasing_cover_relations() + self.decreasing_cover_relations()))

    def __eq__(self, other):
        r"""
        TESTS::

            sage: TamariIntervalPoset(0,[]) == TamariIntervalPoset(0,[])
            True
            sage: TamariIntervalPoset(1,[]) == TamariIntervalPoset(0,[])
            False
            sage: TamariIntervalPoset(3,[(1,2),(3,2)]) == TamariIntervalPoset(3,[(3,2),(1,2)])
            True
            sage: TamariIntervalPoset(3,[(1,2),(3,2)]) == TamariIntervalPoset(3,[(1,2)])
            False
        """
        if (not isinstance(other, TamariIntervalPoset)):
            return False
        return self.size() == other.size() and self._cover_relations == other._cover_relations

    def __le__(self, el2):
        r"""
        TESTS::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip1 <= ip2
            True
            sage: ip1 <= ip1
            True
            sage: ip2 <= ip1
            False
        """
        return self.parent().le(self, el2)

    def __lt__(self, el2):
        r"""
        TESTS::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip1 < ip2
            True
            sage: ip1 < ip1
            False
            sage: ip2 < ip1
            False
        """
        return self.parent().lt(self, el2)

    def __ge__(self, el2):
        r"""
        TESTS::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip1 >= ip2
            False
            sage: ip1 >= ip1
            True
            sage: ip2 >= ip1
            True
        """
        return self.parent().ge(self, el2)

    def __gt__(self, el2):
        r"""
        TESTS::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip1 > ip2
            False
            sage: ip1 > ip1
            False
            sage: ip2 > ip1
            True
        """
        return self.parent().gt(self, el2)

    def contains_interval(self, other):
        r"""
        Return whether the interval represented by ``other`` is contained
        in ``self`` as an interval of the Tamari lattice.

        In terms of interval-posets, it means that all relations of ``self``
        are relations of ``other``.

        INPUT:

        - ``other`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(2,3)])
            sage: ip2.contains_interval(ip1)
            True
            sage: ip3 = TamariIntervalPoset(4,[(2,1)])
            sage: ip2.contains_interval(ip3)
            False
            sage: ip4 = TamariIntervalPoset(3,[(2,3)])
            sage: ip2.contains_interval(ip4)
            False
        """
        if other.size() != self.size():
            return False
        for (i, j) in self._cover_relations:
            if not other.le(i, j):
                return False
        return True

    def lower_contains_interval(self, other):
        r"""
        Return whether the interval represented by ``other`` is contained
        in ``self`` as an interval of the Tamari lattice and if they share
        the same lower bound.

        As interval-posets, it means that ``other`` contains the relations
        of ``self`` plus some extra increasing relations.

        INPUT:

        - ``other`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)]);
            sage: ip2 = TamariIntervalPoset(4,[(4,3)])
            sage: ip2.lower_contains_interval(ip1)
            True
            sage: ip2.contains_interval(ip1) and ip2.lower_binary_tree() == ip1.lower_binary_tree()
            True
            sage: ip3 = TamariIntervalPoset(4,[(4,3),(2,1)])
            sage: ip2.contains_interval(ip3)
            True
            sage: ip2.lower_binary_tree() == ip3.lower_binary_tree()
            False
            sage: ip2.lower_contains_interval(ip3)
            False
        """
        if not self.contains_interval(other):
            return False
        for (i, j) in other.decreasing_cover_relations():
            if not self.le(i, j):
                return False
        return True

    def upper_contains_interval(self, other):
        r"""
        Return whether the interval represented by ``other`` is contained
        in ``self`` as an interval of the Tamari lattice and if they share
        the same upper bound.

        As interval-posets, it means that ``other`` contains the relations
        of ``self`` plus some extra decreasing relations.

        INPUT:

        - ``other`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip2.upper_contains_interval(ip1)
            True
            sage: ip2.contains_interval(ip1) and ip2.upper_binary_tree() == ip1.upper_binary_tree()
            True
            sage: ip3 = TamariIntervalPoset(4,[(1,2),(2,3),(3,4)])
            sage: ip2.upper_contains_interval(ip3)
            False
            sage: ip2.contains_interval(ip3)
            True
            sage: ip2.upper_binary_tree() == ip3.upper_binary_tree()
            False
        """
        if not self.contains_interval(other):
            return False
        for (i, j) in other.increasing_cover_relations():
            if not self.le(i, j):
                return False
        return True

    def is_linear_extension(self, perm):
        r"""
        Return whether the permutation ``perm`` is a linear extension
        of ``self``.

        INPUT:

        - ``perm`` -- a permutation of the size of ``self``

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip.is_linear_extension([1,4,2,3])
            True
            sage: ip.is_linear_extension(Permutation([1,4,2,3]))
            True
            sage: ip.is_linear_extension(Permutation([1,4,3,2]))
            False
        """
        return self._poset.is_linear_extension(perm)

    def contains_binary_tree(self, binary_tree):
        r"""
        Return whether the interval represented by ``self`` contains
        the binary tree ``binary_tree``.

        INPUT:

        - ``binary_tree`` -- a binary tree

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.contains_binary_tree(BinaryTree([[None,[None,[]]],None]))
            True
            sage: ip.contains_binary_tree(BinaryTree([None,[[[],None],None]]))
            True
            sage: ip.contains_binary_tree(BinaryTree([[],[[],None]]))
            False
            sage: ip.contains_binary_tree(ip.lower_binary_tree())
            True
            sage: ip.contains_binary_tree(ip.upper_binary_tree())
            True
            sage: all([ip.contains_binary_tree(bt) for bt in ip.binary_trees()])
            True

        """
        return self.is_linear_extension(binary_tree.to_132_avoiding_permutation())

    def contains_dyck_word(self, dyck_word):
        r"""
        Return whether the interval represented by ``self`` contains
        the Dyck word ``dyck_word``.

        INPUT:

        - ``dyck_word`` -- a Dyck word

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.contains_dyck_word(DyckWord([1,1,1,0,0,0,1,0]))
            True
            sage: ip.contains_dyck_word(DyckWord([1,1,0,1,0,1,0,0]))
            True
            sage: ip.contains_dyck_word(DyckWord([1,0,1,1,0,1,0,0]))
            False
            sage: ip.contains_dyck_word(ip.lower_dyck_word())
            True
            sage: ip.contains_dyck_word(ip.upper_dyck_word())
            True
            sage: all([ip.contains_dyck_word(bt) for bt in ip.dyck_words()])
            True
        """
        return self.contains_binary_tree(dyck_word.to_binary_tree_tamari())

    def intersection(self, other):
        r"""
        Return the interval-poset formed by combining the relations from
        both ``self`` and ``other``. It corresponds to the intersection
        of the two corresponding intervals of the Tamari lattice.

        INPUT:

        - ``other`` -- an interval-poset of the same size as ``self``

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip2 = TamariIntervalPoset(4,[(4,3)])
            sage: ip1.intersection(ip2)
            The tamari interval of size 4 induced by relations [(1, 2), (2, 3), (4, 3)]
            sage: ip3 = TamariIntervalPoset(4,[(2,1)])
            sage: ip1.intersection(ip3)
            Traceback (most recent call last):
            ...
            ValueError: This intersection is empty, it does not correspond to an interval-poset.
            sage: ip4 = TamariIntervalPoset(3,[(2,3)])
            sage: ip2.intersection(ip4)
            Traceback (most recent call last):
            ...
            ValueError: Intersections are only possible on interval-posets of the same size.
        """
        if other.size() != self.size():
            raise ValueError("Intersections are only possible on interval-posets of the same size.")
        try:
            return TamariIntervalPoset(self.size(), self._cover_relations + other._cover_relations)
        except ValueError:
            raise ValueError("This intersection is empty, it does not correspond to an interval-poset.")

    def initial_forest(self):
        r"""
        Return the initial forest of ``self``, i.e., the interval-poset
        formed from only the increasing relations of ``self``.

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(1,2),(3,2),(2,4),(3,4)]).initial_forest()
            The tamari interval of size 4 induced by relations [(1, 2), (2, 4), (3, 4)]
            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.initial_forest() == ip
            True
        """
        return TamariIntervalPoset(self.size(), self.increasing_cover_relations())

    def final_forest(self):
        r"""
        Return the final forest of ``self``, i.e., the interval-poset
        formed with only the decreasing relations of ``self``.

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(2,1),(3,2),(3,4),(4,2)]).final_forest()
            The tamari interval of size 4 induced by relations [(4, 2), (3, 2), (2, 1)]
            sage: ip = TamariIntervalPoset(3,[(2,1),(3,1)])
            sage: ip.final_forest() == ip
            True
        """
        return TamariIntervalPoset(self.size(), self.decreasing_cover_relations())

    def is_initial_interval(self):
        r"""
        Return if ``self`` corresponds to an initial interval of the Tamari
        lattice, i.e. if its lower end is the smallest element of the lattice.
        It consists of checking that ``self`` does not contain any decreasing
        relations.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4, [(1, 2), (2, 4), (3, 4)])
            sage: ip.is_initial_interval()
            True
            sage: ip.lower_dyck_word()
            [1, 0, 1, 0, 1, 0, 1, 0]
            sage: ip = TamariIntervalPoset(4, [(1, 2), (2, 4), (3, 4), (3, 2)])
            sage: ip.is_initial_interval()
            False
            sage: ip.lower_dyck_word()
            [1, 0, 1, 1, 0, 0, 1, 0]
            sage: all([ DyckWord([1,0,1,0,1,0]).tamari_interval(dw).is_initial_interval() for dw in DyckWords(3)])
            True
        """
        return self.decreasing_cover_relations() == []

    def is_final_interval(self):
        r"""
        Return if ``self`` corresponds to a final interval of the Tamari
        lattice, i.e. if its upper end is the largest element of the lattice.
        It consists of checking that ``self`` does not contain any increasing
        relations.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4, [(4, 3), (3, 1), (2, 1)])
            sage: ip.is_final_interval()
            True
            sage: ip.upper_dyck_word()
            [1, 1, 1, 1, 0, 0, 0, 0]
            sage: ip = TamariIntervalPoset(4, [(4, 3), (3, 1), (2, 1), (2, 3)])
            sage: ip.is_final_interval()
            False
            sage: ip.upper_dyck_word()
            [1, 1, 0, 1, 1, 0, 0, 0]
            sage: all([ dw.tamari_interval(DyckWord([1, 1, 1, 0, 0, 0])).is_final_interval() for dw in DyckWords(3)])
            True
        """
        return self.increasing_cover_relations() == []

    def lower_binary_tree(self):
        r"""
        Return the lowest binary tree in the interval of the Tamari
        lattice represented by ``self``.

        This is a binary tree. It is the shape of the unique binary
        search tree whose left-branch ordered forest (i.e., the result
        of applying
        :meth:`~sage.combinat.binary_tree.BinaryTree.to_ordered_tree_left_branch`
        and cutting off the root) is the final forest of ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.lower_binary_tree()
            [[., .], [[., [., .]], [., .]]]
            sage: TamariIntervalPosets.final_forest(ip.lower_binary_tree()) == ip.final_forest()
            True
            sage: ip == TamariIntervalPosets.from_binary_trees(ip.lower_binary_tree(),ip.upper_binary_tree())
            True
        """
        return self.min_linear_extension().binary_search_tree_shape(left_to_right=False)

    def lower_dyck_word(self):
        r"""
        Return the lowest Dyck word in the interval of the Tamari lattice
        represented by ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.lower_dyck_word()
            [1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0]
            sage: TamariIntervalPosets.final_forest(ip.lower_dyck_word()) == ip.final_forest()
            True
            sage: ip == TamariIntervalPosets.from_dyck_words(ip.lower_dyck_word(),ip.upper_dyck_word())
            True
        """
        return self.lower_binary_tree().to_dyck_word_tamari()

    def upper_binary_tree(self):
        r"""
        Return the highest binary tree in the interval of the Tamari
        lattice represented by ``self``.

        This is a binary tree. It is the shape of the unique binary
        search tree whose right-branch ordered forest (i.e., the result
        of applying
        :meth:`~sage.combinat.binary_tree.BinaryTree.to_ordered_tree_right_branch`
        and cutting off the root) is the initial forest of ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.upper_binary_tree()
            [[., .], [., [[., .], [., .]]]]
            sage: TamariIntervalPosets.initial_forest(ip.upper_binary_tree()) == ip.initial_forest()
            True
            sage: ip == TamariIntervalPosets.from_binary_trees(ip.lower_binary_tree(),ip.upper_binary_tree())
            True
        """
        return self.max_linear_extension().binary_search_tree_shape(left_to_right=False)

    def upper_dyck_word(self):
        r"""
        Return the highest Dyck word in the interval of the Tamari lattice
        represented by ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.upper_dyck_word()
            [1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0]
            sage: TamariIntervalPosets.initial_forest(ip.upper_dyck_word()) == ip.initial_forest()
            True
            sage: ip == TamariIntervalPosets.from_dyck_words(ip.lower_dyck_word(),ip.upper_dyck_word())
            True
        """
        return self.upper_binary_tree().to_dyck_word_tamari()

    def sub_poset(self, start, end):
        r"""
        Return the renormalized sub-poset of ``self`` consisting solely
        of integers from ``start`` (inclusive) to ``end`` (not inclusive).

        "Renormalized" means that these integers are relabelled
        `1,2,\ldots,k` in the obvious way (i.e., by subtracting
        ``start - 1``).

        INPUT:

        - ``start`` -- an integer, the starting vertex (inclusive)
        - ``end`` -- an integer, the ending vertex (not inclusive)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.sub_poset(1,3)
            The tamari interval of size 2 induced by relations [(1, 2)]
            sage: ip.sub_poset(1,4)
            The tamari interval of size 3 induced by relations [(1, 2), (3, 2)]
            sage: ip.sub_poset(1,5)
            The tamari interval of size 4 induced by relations [(1, 2), (4, 3), (3, 2)]
            sage: ip.sub_poset(1,7) == ip
            True
            sage: ip.sub_poset(1,1)
            The tamari interval of size 0 induced by relations []
        """
        if start < 1 or start > end or end > self.size() + 1:
            raise ValueError("Invalid starting or ending value, accepted: 1 <= start <= end <= size+1")
        if start == end:
            return TamariIntervalPoset(0, [])
        relations = [(i - start + 1, j - start + 1) for (i, j) in self.increasing_cover_relations() if i >= start and j < end]
        relations.extend([(j - start + 1, i - start + 1) for (j, i) in self.decreasing_cover_relations() if i >= start and j < end])
        return TamariIntervalPoset(end - start, relations)

    def min_linear_extension(self):
        r"""
        Return the minimal permutation for the right weak order which is
        a linear extension of ``self``.

        This is also the minimal permutation in the sylvester
        class of ``self.lower_binary_tree()`` and is a 312-avoiding
        permutation.

        The right weak order is also known as the right permutohedron
        order. See
        :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`
        for its definition.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip.min_linear_extension()
            [1, 2, 4, 3]
            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)])
            sage: ip.min_linear_extension()
            [1, 4, 3, 6, 5, 2]
            sage: ip = TamariIntervalPoset(0,[])
            sage: ip.min_linear_extension()
            []
            sage: ip = TamariIntervalPoset(5, [(1, 4), (2, 4), (3, 4), (5, 4)]); ip
            The tamari interval of size 5 induced by relations [(1, 4), (2, 4), (3, 4), (5, 4)]
            sage: ip.min_linear_extension()
            [1, 2, 3, 5, 4]

        """
        # The min linear extension is build by postfix-reading the
        # final forest of ``self``.
        def add(perm, i):
            r"""
            Internal recursive method to compute the min linear extension.
            """
            for j in self.decreasing_children(i):
                add(perm, j)
            perm.append(i)
        perm = []
        for i in self.decreasing_roots():
            add(perm, i)
        return Permutation(perm)

    def max_linear_extension(self):
        r"""
        Return the maximal permutation for the right weak order which is
        a linear extension of ``self``.

        This is also the maximal permutation in the sylvester
        class of ``self.upper_binary_tree()`` and is a 132-avoiding
        permutation.

        The right weak order is also known as the right permutohedron
        order. See
        :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`
        for its definition.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip.max_linear_extension()
            [4, 1, 2, 3]
            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.max_linear_extension()
            [6, 4, 5, 3, 1, 2]
            sage: ip = TamariIntervalPoset(0,[]); ip
            The tamari interval of size 0 induced by relations []
            sage: ip.max_linear_extension()
            []
            sage: ip = TamariIntervalPoset(5, [(1, 4), (2, 4), (3, 4), (5, 4)]); ip
            The tamari interval of size 5 induced by relations [(1, 4), (2, 4), (3, 4), (5, 4)]
            sage: ip.max_linear_extension()
            [5, 3, 2, 1, 4]

        """
        # The max linear extension is build by right-to-left
        # postfix-reading the initial forest of ``self``. The
        # right-to-leftness here is ensured by the fact that
        # :meth:`increasing_children` and :meth:`increasing_roots`
        # output their results in decreasing order.
        def add(perm, i):
            r"""
            Internal recursive method to compute the max linear extension.
            """
            for j in self.increasing_children(i):
                add(perm, j)
            perm.append(i)
        perm = []
        for i in self.increasing_roots():
            add(perm, i)
        return Permutation(perm)

    def linear_extensions(self):
        r"""
        Return an iterator on the permutations which are linear
        extensions of ``self``.

        They form an interval of the right weak order (also called the
        right permutohedron order -- see
        :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`
        for a definition).

        EXAMPLES::

            sage: ip = TamariIntervalPoset(3,[(1,2),(3,2)])
            sage: list(ip.linear_extensions())
            [[3, 1, 2], [1, 3, 2]]
            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: list(ip.linear_extensions())
            [[4, 1, 2, 3], [1, 4, 2, 3], [1, 2, 4, 3]]
        """
        for ext in self._poset.linear_extensions():
            yield Permutation(ext)

    def lower_contained_intervals(self):
        r"""
        If ``self`` represents the interval `[t_1, t_2]` of the Tamari
        lattice, return an iterator on all intervals `[t_1,t]` with
        `t \leq t_2` for the Tamari lattice.

        In terms of interval-posets, it corresponds to adding all possible
        relations of the form `n` precedes `m` with `n<m`.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.lower_contained_intervals())
            [The tamari interval of size 4 induced by relations [(2, 4), (3, 4), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(1, 4), (2, 4), (3, 4), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(2, 3), (3, 4), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(1, 4), (2, 3), (3, 4), (3, 1), (2, 1)]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: len(list(ip.lower_contained_intervals()))
            14
        """
        size = self._size
        yield self
        r"""
        we try to add links recursively in this order :
        1 -> 2
        2 -> 3
        1 -> 3
        3 -> 4
        2 -> 4
        1 -> 4
        ...
        ("Link" means "relation of the poset".)
        """
        def add_relations(poset, n, m):
            r"""
            Internal recursive method to generate all possible intervals.
            At every step during the iteration, we have n < m and every
            i satisfying n < i < m satisfies that i precedes m in the
            poset ``poset`` (except when m > size).
            """
            if n <= 0:
                #if n<=0, then we go to the next m
                n = m
                m += 1
            if m > size:
                #if m>size, it's finished
                return

            if poset.le(n, m):
                #there is already a link n->m, so we go to the next n
                for pos in add_relations(poset, n - 1, m):
                    yield pos
            elif poset.le(m, n):
                #there is an inverse link m->n, we know we won't be able
                #to create a link i->m with i<=n, so we go to the next m
                for pos in add_relations(poset, m, m + 1):
                    yield pos
            else:
                #there is no link n->m
                #first option : we don't create the link and go to the next m
                #(since the lack of a link n->m forbids any links i->m
                #with i<n)
                for pos in add_relations(poset, m, m + 1):
                    yield pos
                #second option : we create the link
                #(this is allowed because links i->m already exist for all
                #n<i<m, or else we wouldn't be here)
                poset = TamariIntervalPoset(poset.size(), poset._cover_relations + ((n, m),))
                yield poset
                #and then, we go to the next n
                for pos in add_relations(poset, n - 1, m):
                    yield pos

        for inter in add_relations(self, 1, 2):
            yield inter

    def interval_cardinality(self):
        r"""
        Return the cardinality of the interval, i.e., the number of elements
        (binary trees or Dyck words) in the interval represented by ``self``.

        Not to be confused with :meth:`size` which is the number of
        vertices.

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)]).interval_cardinality()
            4
            sage: TamariIntervalPoset(4,[]).interval_cardinality()
            14
            sage: TamariIntervalPoset(4,[(1,2),(2,3),(3,4)]).interval_cardinality()
            1
        """
        return len(list(self.lower_contained_intervals()))

    def binary_trees(self):
        r"""
        Return an iterator on all the binary trees in the interval
        represented by ``self``.

        EXAMPLES::

            sage: list(TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)]).binary_trees())
            [[., [[., [., .]], .]],
             [[., [., [., .]]], .],
             [., [[[., .], .], .]],
             [[., [[., .], .]], .]]
            sage: set(TamariIntervalPoset(4,[]).binary_trees()) == set(BinaryTrees(4))
            True
        """
        for ip in self.lower_contained_intervals():
            yield ip.upper_binary_tree()

    def dyck_words(self):
        r"""
        Return an iterator on all the Dyck words in the interval
        represented by ``self``.

        EXAMPLES::

            sage: list(TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)]).dyck_words())
            [[1, 1, 1, 0, 0, 1, 0, 0],
             [1, 1, 1, 0, 0, 0, 1, 0],
             [1, 1, 0, 1, 0, 1, 0, 0],
             [1, 1, 0, 1, 0, 0, 1, 0]]
            sage: set(TamariIntervalPoset(4,[]).dyck_words()) == set(DyckWords(4))
            True
        """
        for ip in self.lower_contained_intervals():
            yield ip.upper_dyck_word()

    def maximal_chain_intervals(self):
        r"""
        Return an iterator on the upper contained intervals of one
        maximal chain of ``self``.

        If ``self`` represents the interval `[T_1,T_2]` of the Tamari
        lattice, this returns intervals `[T',T_2]` with `T'` following
        one maximal chain between `T_1` and `T_2`. This chain actually
        has the property that the number of decreasing relations
        increases by `1` at each step.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.maximal_chain_intervals())
            [The tamari interval of size 4 induced by relations [(2, 4), (3, 4), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(2, 4), (3, 4), (4, 1), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(2, 4), (3, 4), (4, 1), (3, 2), (2, 1)]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: list(ip.maximal_chain_intervals())
            [The tamari interval of size 4 induced by relations [],
             The tamari interval of size 4 induced by relations [(2, 1)],
             The tamari interval of size 4 induced by relations [(3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(4, 1), (3, 1), (2, 1)],
             The tamari interval of size 4 induced by relations [(4, 1), (3, 2), (2, 1)],
             The tamari interval of size 4 induced by relations [(4, 2), (3, 2), (2, 1)],
             The tamari interval of size 4 induced by relations [(4, 3), (3, 2), (2, 1)]]
        """
        yield self
        size = self.size()
        rel = list(self._cover_relations)
        ti = self
        # we add relations in this order
        # 2 -> 1
        # 3 -> 1
        # 4 -> 1
        # 3 -> 2
        # 4 -> 2
        # 4 -> 3
        for i in xrange(1, size):
            for j in xrange(i + 1, size + 1):
                # Loop invariant: the poset ti has a relation k->i
                # for every k satisfying i < k < j.
                if ti.le(j, i):
                    #the relation j->i is already there, we go to the next j
                    continue
                if ti.le(i, j):
                    #there is a relation i->j which forbids any (>=j)->i
                    # we go to the next i
                    break
                # there is no j->i or i->j, so we add j->i
                rel.append((j, i))
                ti = TamariIntervalPoset(size, rel)
                yield ti

    def maximal_chain_binary_trees(self):
        r"""
        Return an iterator on the binary trees forming a maximal chain of
        ``self`` (regarding ``self`` as an interval of the Tamari
        lattice).

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.maximal_chain_binary_trees())
            [[[., [[., .], .]], .], [., [[[., .], .], .]], [., [[., [., .]], .]]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: list(ip.maximal_chain_binary_trees())
            [[[[[., .], .], .], .],
             [[[., [., .]], .], .],
             [[., [[., .], .]], .],
             [., [[[., .], .], .]],
             [., [[., [., .]], .]],
             [., [., [[., .], .]]],
             [., [., [., [., .]]]]]
        """
        for it in self.maximal_chain_intervals():
            yield it.lower_binary_tree()

    def maximal_chain_dyck_words(self):
        r"""
        Return an iterator on the Dyck words forming a maximal chain of
        ``self`` (regarding ``self`` as an interval of the Tamari
        lattice).

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.maximal_chain_dyck_words())
            [[1, 1, 0, 1, 0, 0, 1, 0], [1, 1, 0, 1, 0, 1, 0, 0], [1, 1, 1, 0, 0, 1, 0, 0]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: list(ip.maximal_chain_dyck_words())
            [[1, 0, 1, 0, 1, 0, 1, 0],
             [1, 1, 0, 0, 1, 0, 1, 0],
             [1, 1, 0, 1, 0, 0, 1, 0],
             [1, 1, 0, 1, 0, 1, 0, 0],
             [1, 1, 1, 0, 0, 1, 0, 0],
             [1, 1, 1, 0, 1, 0, 0, 0],
             [1, 1, 1, 1, 0, 0, 0, 0]]
        """
        for it in self.maximal_chain_intervals():
            yield it.lower_dyck_word()

    def length_of_maximal_chain(self):
        r"""
        Return the length of a maximal chain of ``self`` (regarding
        ``self`` as an interval of the Tamari lattice).

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.length_of_maximal_chain()
            3
            sage: ip = TamariIntervalPoset(4,[])
            sage: ip.length_of_maximal_chain()
            7
            sage: ip = TamariIntervalPoset(3,[])
            sage: ip.length_of_maximal_chain()
            4
        """
        return len(list(self.maximal_chain_intervals()))

    def initial_rise(self):
        r"""
        Return the initial rise of ``self``. The initial rise of an 
        interval-poset is defined as being the initial rise of its upper
        Dyck path [BFP]_.
        
        It can be computed directly from interval-poset [CCP]_: it is the first 
        vertex `k` such that `k` precedes `k+1`. If there is not such vertex,
        then it is the size of the interval-poset.
        
        REFERENCES:

        .. [BFP] The Number of intervals in the m-Tamari lattices, M. Bousquet-Melou,
         E. Fusy, L.-F. Preville-Ratelle
        .. [CCP] Two bijections on Tamari intervals, G. Chatel, V. Pons

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.initial_rise()
            3
            sage: ip.upper_dyck_word().number_of_initial_rises()
            3
            sage: all( [ip.initial_rise() == ip.upper_dyck_word().number_of_initial_rises() for ip in TamariIntervalPosets(3)])
            True
        """
        for i in xrange(1,self.size()):
            if self.le(i,i+1):
                return i
        return self.size()
        
    def number_of_contacts(self):
        r"""
        Return the number of contacts of ``self``. The number of contacts
        of a Tamari interval is defined to be the number of contacts of its
        lower Dyck path [BFP]_. We don't include the initial contact.

        In terms of interval-posets, this is the number of decreasing roots,
        i.e., the number of connected compentends of the final forest [CCP]_.

        REFERENCES:

        .. [BFP] The Number of intervals in the m-Tamari lattices, M. Bousquet-Melou,
         E. Fusy, L.-F. Preville-Ratelle
        .. [CCP] Two bijections on Tamari intervals, G. Chatel, V. Pons

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.number_of_contacts()
            2
            sage: ip.lower_dyck_word().number_of_touch_points()
            2
            sage: all( [ip.number_of_contacts() == ip.lower_dyck_word().number_of_touch_points() for ip in TamariIntervalPosets(3)])
            True

        """
        return len(self.decreasing_roots())

    def lower_contacts_composition(self, ip2, p=None):
        r"""
        Return the lower contacts composition of ``self`` and ``ip2`` as 
        described in [CCP]_. If ``p`` is not ``None``, then the partial 
        composition is performed.

        If ``self.size()`` is `n` and ``ip2.size()`` is `m`, the complete lower 
        contacts composition of ``self`` and ``ip2`` consists of a set of 
        interval-posets of size `n+m+1`. If ``self.number_of_contacts()`` is `c_1` 
        and ``ip2.number_of_contacts()`` is `c_2`, then the number of interval-posets 
        of the composition is `c_2+1` and the number of contacts of the 
        resulting interval-posets are `c_1 + c_2 + 1, c_1 + c_2, c_1 + c_2 - 1, \dots, c_1 + 1`.

        If ``p`` is not ``None``, then only a partial composition is performed 
        and a unique interval-poset is returned with its number of contacts
        equal to `c_1 + 1 + p`.

        For a detailed description of the composition, please refer to [CCP]_.

        REFERENCES:

        .. [CCP] Two bijections on Tamari intervals, G. Chatel, V. Pons

        INPUT:

        - ``ip2`` -- an interval-poset of size `m`
        - ``p`` -- (default: ``None``) a Integer between 0 and ``ip2.number_of_contacts()``

        OUTPUT:

        A list of ``ip2.number_of_contacts()`` interval-posets of size `n+m+1`
        if ``p`` is ``None`` or a unique interval-poset of size `n+m+1` if 
        ``p`` is a number.

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(2,1),(3,1),(2,4),(3,4)])
            sage: ip2 = TamariIntervalPoset(5,[(1,2),(3,2),(5,4)])
            sage: ip1.lower_contacts_composition(ip2)
            [The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (8, 7), (3, 1), (2, 1)],
             The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (8, 7), (6, 5), (3, 1), (2, 1)],
             The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (8, 7), (7, 5), (6, 5), (3, 1), (2, 1)],
             The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (9, 5), (8, 7), (7, 5), (6, 5), (3, 1), (2, 1)]]
            sage: ip1.number_of_contacts()
            2
            sage: ip2.number_of_contacts()
            3    
            sage: [ip.number_of_contacts() for ip in ip1.lower_contacts_composition(ip2)]
            [6, 5, 4, 3]

        partial composition::

            sage: ip1.lower_contacts_composition(ip2,2)
            The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (8, 7), (7, 5), (6, 5), (3, 1), (2, 1)]
            sage: ip1.lower_contacts_composition(ip2,0)
            The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (8, 7), (3, 1), (2, 1)]
            sage: ip1.lower_contacts_composition(ip2,5)
            Traceback (most recent call last):
            ...
            ValueError: Invalid composition parameter
            sage: [ip1.lower_contacts_composition(ip2,p) for p in xrange(ip2.number_of_contacts()+1)] == ip1.lower_contacts_composition(ip2)
            True


        composition with an empty interval-poset::

            sage: ip0 = TamariIntervalPoset(0,[])
            sage: ip1.lower_contacts_composition(ip0)
            [The tamari interval of size 5 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (3, 1), (2, 1)]]
            sage: ip1.lower_contacts_composition(ip0,0)
            The tamari interval of size 5 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (3, 1), (2, 1)]
            sage: ip0.lower_contacts_composition(ip2)
            [The tamari interval of size 6 induced by relations [(2, 3), (6, 5), (4, 3)],
             The tamari interval of size 6 induced by relations [(2, 3), (6, 5), (4, 3), (2, 1)],
             The tamari interval of size 6 induced by relations [(2, 3), (6, 5), (4, 3), (3, 1), (2, 1)],
             The tamari interval of size 6 induced by relations [(2, 3), (6, 5), (5, 1), (4, 3), (3, 1), (2, 1)]]
            sage: ip0.lower_contacts_composition(ip2,2)
            The tamari interval of size 6 induced by relations [(2, 3), (6, 5), (4, 3), (3, 1), (2, 1)]
            sage: ip0.lower_contacts_composition(ip0)
            [The tamari interval of size 1 induced by relations []]
        """
        return TamariIntervalPosets.lower_contacts_composition(self,ip2,p)

    def lower_contacts_decomposition(self):
        r"""
        Return the the lower contacts decompostion of ``self`` into two
        interval-posets of summed sizes equal to ``self.size()-1``. The 
        decomposition algorithm is described in [CCP]_. 

        It is the inverse operation of the composition performed in 
        :method:`TamariIntervalPoset.lower_contacts_composition`. It consists of selecting
        the minimal increasing root of ``self`` and splitting ``self`` into
        two sub posets relatively to the root.

        REFERENCES:

        .. [CCP] Two bijections on Tamari intervals, G. Chatel, V. Pons

        OUTPUT:

        A tuple of 3 elements:

        - an interval-poset ``ip1`` which is the left part of the composition
        - an interval-poset ``ip2`` which is the right part of the composition
        - a parameter ``p``
        
        It is the unique triplet such that ``self == self.lower_contacts_composition(ip1,ip2,p)``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.lower_contacts_decomposition()
            (The tamari interval of size 0 induced by relations [],
             The tamari interval of size 3 induced by relations [(1, 3), (2, 3)],
             2)
            sage: ip == TamariIntervalPosets.lower_contacts_composition(*(ip.lower_contacts_decomposition()))
            True
            sage: ip1 = TamariIntervalPoset(4,[(2,1),(3,1),(2,4),(3,4)])
            sage: ip2 = TamariIntervalPoset(5,[(1,2),(3,2),(5,4)])
            sage: ip1.lower_contacts_composition(ip2,2).lower_contacts_decomposition() == (ip1,ip2,2)
            True

            sage: u = TamariIntervalPoset(1,[])
            sage: u.lower_contacts_decomposition()
            (The tamari interval of size 0 induced by relations [],
             The tamari interval of size 0 induced by relations [],
             0)
            sage: ip0 = TamariIntervalPoset(0,[])
            sage: ip0.lower_contacts_decomposition()
            Traceback (most recent call last):
            ...
            ValueError: The empty interval-poset cannot be decomposed

        TESTS::

            sage: all([ip == TamariIntervalPosets.lower_contacts_composition(*(ip.lower_contacts_decomposition())) for ip in TamariIntervalPosets(4)])
            True

        """
        if self.size()==0:
            raise ValueError, "The empty interval-poset cannot be decomposed"%()
        root = self.increasing_roots()[-1]
        relations_ip1 = [(i,j) for (i,j) in self._cover_relations if i <root and j<root]
        relations_ip2 = [(i-root,j-root) for (i,j) in self._cover_relations if i>root and j>root]
        return (TamariIntervalPoset(root-1,relations_ip1),TamariIntervalPoset(self.size()-root,relations_ip2),len(self.decreasing_children(root)))

    def initial_rise_composition(self,ip2,p=None):
        r"""
        Return the initial rise composition of ``self`` and ``ip2`` as 
        described in [CCP]_. If ``p`` is not ``None``, then the partial 
        composition is performed.

        If ``self.size()`` is `n` and ``ip2.size()`` is `m`, the complete initial 
        rise composition of ``self`` and ``ip2`` consists of a set of 
        interval-posets of size `n+m+1`. If ``self.initial_rise()`` is `ir_1` 
        and ``ip2.initial_rise()`` is `ir_2`, then the number of interval-posets 
        of the composition is `ir_2+ir_1` and the initial rises of the 
        resulting interval-posets are `ir_1 + ir_2 + 1, ir_1 + ir_2, ir_1 + ir_2 - 1, \dots, ir_1 + 1`.
        
        If ``p`` is not ``None``, then only a partial composition is performed 
        and a unique interval-poset is returned with its initial rise
        equal to `ir_1 + 1 + p`.

        For a detailed description of the composition, please refer to [CCP]_.

        REFERENCES:

        .. [CCP] Two bijections on Tamari intervals, G. Chatel, V. Pons

        INPUT:

        - ``ip2`` -- an interval-poset of size `m`
        - ``p`` -- (default: ``None``) a Integer between 0 and ``ip2.initial_rise()``

        OUTPUT:

        A list of ``ip2.initial_rise()`` interval-posets of size `n+m+1`
        if ``p`` is ``None`` or a unique interval-poset of size `n+m+1` if 
        ``p`` is a number.

        EXAMPLES::
        
            sage: ip1 = TamariIntervalPoset(3,[(1,2),(3,2)])
            sage: ip2 = TamariIntervalPoset(5,[(2,3),(3,5),(4,5),(2,1),(3,1),(4,3)])
            sage: ip1.initial_rise_composition(ip2)
            [The tamari interval of size 9 induced by relations [(1, 8), (2, 8), (3, 5), (4, 5), (5, 7), (6, 7), (7, 8), (9, 8), (7, 1), (6, 5), (5, 1), (4, 2), (3, 2), (2, 1)],
             The tamari interval of size 9 induced by relations [(1, 8), (2, 8), (3, 4), (4, 5), (5, 7), (6, 7), (7, 8), (9, 8), (7, 1), (6, 5), (5, 1), (4, 2), (3, 2), (2, 1)],
             The tamari interval of size 9 induced by relations [(1, 8), (2, 3), (3, 8), (4, 5), (5, 7), (6, 7), (7, 8), (9, 8), (7, 1), (6, 5), (5, 3), (4, 3), (3, 1), (2, 1)]]
            sage: ip1.initial_rise()
            1
            sage: ip2.initial_rise()
            2
            sage: [ip.initial_rise() for ip in ip1.initial_rise_composition(ip2)]
            [4, 3, 2]

        partial composition::

            sage: ip1.initial_rise_composition(ip2,2)
            The tamari interval of size 9 induced by relations [(1, 8), (2, 3), (3, 8), (4, 5), (5, 7), (6, 7), (7, 8), (9, 8), (7, 1), (6, 5), (5, 3), (4, 3), (3, 1), (2, 1)]
            sage: ip1.initial_rise_composition(ip2,0)
            The tamari interval of size 9 induced by relations [(1, 8), (2, 8), (3, 5), (4, 5), (5, 7), (6, 7), (7, 8), (9, 8), (7, 1), (6, 5), (5, 1), (4, 2), (3, 2), (2, 1)]
            sage: ip1.initial_rise_composition(ip2,4)
            Traceback (most recent call last):
            ...
            ValueError: Invalid composition parameter
            sage: [ip1.initial_rise_composition(ip2,p) for p in xrange(ip2.initial_rise()+1)] == ip1.initial_rise_composition(ip2)
            True

        composition with an empty interval-poset::

            sage: ip0 = TamariIntervalPoset(0,[])
            sage: ip1.initial_rise_composition(ip0)
            [The tamari interval of size 4 induced by relations [(1, 3), (2, 3), (4, 3), (2, 1)]]
            sage: ip1.initial_rise_composition(ip0,0)
            The tamari interval of size 4 induced by relations [(1, 3), (2, 3), (4, 3), (2, 1)]
            sage: ip0.initial_rise_composition(ip2)
            [The tamari interval of size 6 induced by relations [(2, 4), (3, 4), (4, 6), (5, 6), (5, 4), (3, 1), (2, 1)],
             The tamari interval of size 6 induced by relations [(2, 3), (3, 4), (4, 6), (5, 6), (5, 4), (3, 1), (2, 1)],
             The tamari interval of size 6 induced by relations [(1, 2), (3, 4), (4, 6), (5, 6), (5, 4), (4, 2), (3, 2)]]
            sage: ip0.initial_rise_composition(ip2,2)
            The tamari interval of size 6 induced by relations [(1, 2), (3, 4), (4, 6), (5, 6), (5, 4), (4, 2), (3, 2)]
            sage: ip0.initial_rise_composition(ip0)
            [The tamari interval of size 1 induced by relations []]
        """
        return TamariIntervalPosets.initial_rise_composition(self,ip2,p)
        
    def initial_rise_decomposition(self):
        r"""
        Return the the initial rise decompostion of ``self`` into two
        interval-posets of summed sizes equal to ``self.size()-1``. The 
        decomposition algorithm is described in [CCP]_. 

        It is the inverse operation of the composition performed in 
        :method:`TamariIntervalPoset.initial_rise_composition`.

        REFERENCES:

        .. [CCP] Two bijections on Tamari intervals, G. Chatel, V. Pons

        OUTPUT:

        A tuple of 3 elements:

        - an interval-poset ``ip1`` which is the left part of the composition
        - an interval-poset ``ip2`` which is the right part of the composition
        - a parameter ``p``
        
        It is the unique triplet such that ``self == self.initial_rise_composition(ip1,ip2,p)``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(3,[(1,2),(3,2)])
            sage: ip.initial_rise_decomposition()
            (The tamari interval of size 0 induced by relations [],
             The tamari interval of size 2 induced by relations [(2, 1)],
             2)
            sage: ip == TamariIntervalPosets.initial_rise_composition(*(ip.initial_rise_decomposition()))
            True
            sage: ip1 = TamariIntervalPoset(3,[(1,2),(3,2)])
            sage: ip2 = TamariIntervalPoset(5,[(2,3),(3,5),(4,5),(2,1),(3,1),(4,3)])
            sage: ip1.initial_rise_composition(ip2,1).initial_rise_decomposition() == (ip1,ip2,1)
            True

            sage: u = TamariIntervalPoset(1,[])
            sage: u.initial_rise_decomposition()
            (The tamari interval of size 0 induced by relations [],
             The tamari interval of size 0 induced by relations [],
             0)
            sage: ip0 = TamariIntervalPoset(0,[])
            sage: ip0.initial_rise_decomposition()
            Traceback (most recent call last):
            ...
            ValueError: The empty interval-poset cannot be decomposed

        TESTS::

            sage: all([ip == TamariIntervalPosets.initial_rise_composition(*(ip.initial_rise_decomposition())) for ip in TamariIntervalPosets(4)])
            True
        """

        if self.size()==0:
            raise ValueError, "The empty interval-poset cannot be decomposed"%()

        ### separating ip1 and ip2 ##

        relations_ip1 = []
        relations_ip2 = []

        ir = self.initial_rise()
        a = self.decreasing_parent(ir)
        while a!=None:
            b = self.increasing_parent(a)
            if (b is not None and self.le(b-1,a)) or self.le(self.size(),a):
                break
            a = self.decreasing_parent(a)

        if a!=None:
            if b is None:
                b = self.size()+1
            n2 = b-a-2
        else:
            a=0
            b = self.size()+1
            n2 = self.size()-1
            
        ## increasing relations ##

        for (i,j) in self.increasing_cover_relations():
            if j<=a:
                relations_ip1.append((i,j))
            if i<=a and j>=b:
                relations_ip1.append((i,j-n2-1))
            if i>=b:
                relations_ip1.append((i-n2-1,j-n2-1))
            if i>a and j<b:
                if i<ir:
                    relations_ip2.append((i-a,j-a-1))
                if i>ir:
                    relations_ip2.append((i-a-1,j-a-1))
        
        ## decreasing relations ##

        # for ip1
        for (j,i) in self.decreasing_cover_relations():
            if i>=b:
                relations_ip1.append((j-n2-1,i-n2-1))
            if i<=a and j>=b:
                relations_ip1.append((j-n2-1,i))
            if j<=a:
                relations_ip1.append((j,i))
        
        # for ip2
        children = []
        for i in xrange(b-1,ir,-1):
            while len(children)>0:
                j = children.pop()
                if self.le(j,i):
                    relations_ip2.append((j-1-a,i-1-a))
                else:
                    children.append(j)
                    break 
            children.append(i)
        new_children = [i-1 for i in children]
        children.append(ir)
        for i in xrange(ir-1,a,-1):
            while len(new_children)>0:
                j = children.pop()
                nj = new_children.pop()
                if self.le(j,i):
                    relations_ip2.append((nj-a,i-a))
                else:
                    children.append(j)
                    new_children.append(nj)
                    break
            children.append(i)
            new_children.append(i)
        
        ip1 = TamariIntervalPoset(self.size()-n2-1,relations_ip1)
        ip2 = TamariIntervalPoset(n2,relations_ip2)
        return ip1, ip2, a + ip2.initial_rise() + 1 - ir
        
    def initial_rise_involution(self):
        r"""
        Return the image of ``self`` by the initial rise involution. This 
        involution is described in [CCP]_ : it consists of a recursive 
        decomposition-recomposition of ``self`` through the methods :method:`self.lower_contacts_decomposition`
        and :method:`TamariIntervalPosets.initial_rise_composition`. 

        The main purpose of the involution is to exchange the initial rise
        and number of contacts statistics.

        REFERENCES:

        .. [CCP] Two bijections on Tamari intervals, G. Chatel, V. Pons

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: image = ip.initial_rise_involution(); image
            The tamari interval of size 4 induced by relations [(2, 3), (4, 3)]
            sage: ip.number_of_contacts(), ip.initial_rise()
            (2, 3)
            sage: image.number_of_contacts(), image.initial_rise()
            (3, 2)
            sage: image.initial_rise_involution() == ip
            True

        TESTS::

            sage: all([ip == ip.initial_rise_involution().initial_rise_involution() for ip in TamariIntervalPosets(4)])
            True
            sage: all([ip.number_of_contacts()==ip.initial_rise_involution().initial_rise() for ip in TamariIntervalPosets(4)])
            True
        """
        if self.size()==0 or self.size() ==1:
            return self
        dec = self.lower_contacts_decomposition()
        return TamariIntervalPosets.initial_rise_composition(dec[0].initial_rise_involution(), dec[1].initial_rise_involution(),dec[2])

# Abstract class to serve as a Factory ; no instances are created.
class TamariIntervalPosets(UniqueRepresentation, Parent):
    r"""
    Factory for interval-posets.

    INPUT:

    - ``size`` -- (optional) an integer

    OUTPUT:

    - the set of all interval-posets (of the given ``size`` if specified)

    EXAMPLES::

        sage: TamariIntervalPosets()
        Interval-posets

        sage: TamariIntervalPosets(2)
        Interval-posets of size 2

    .. NOTE:: this is a factory class whose constructor returns instances of
              subclasses.
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        r"""
        TESTS::

            sage: from sage.combinat.interval_posets import TamariIntervalPosets_all, TamariIntervalPosets_size
            sage: isinstance(TamariIntervalPosets(2), TamariIntervalPosets)
            True
            sage: isinstance(TamariIntervalPosets(), TamariIntervalPosets)
            True
            sage: TamariIntervalPosets(2) is TamariIntervalPosets_size(2)
            True
            sage: TamariIntervalPosets() is TamariIntervalPosets_all()
            True
        """
        if n is None:
            return TamariIntervalPosets_all()
        else:
            if not (isinstance(n, (Integer, int)) and n >= 0):
                raise ValueError("n must be a non negative integer")
            return TamariIntervalPosets_size(Integer(n))

    @staticmethod
    def check_poset(poset):
        r"""
        Check if the given poset ``poset`` is a interval-poset, that is,
        if it satisfies the following properties:

        - Its labels are exactly `1,\dots,n` where `n` is its size.
        - If `a<c` (as numbers) and `a` precedes `c`, then `b` precedes
          `c` for all `b` such that `a<b<c`.
        - If `a<c` (as numbers) and `c` precedes `a`, then `b` precedes
          `a` for all `b` such that `a<b<c`.

        INPUT:

        - ``poset`` -- a finite labeled poset

        EXAMPLES::

            sage: p = Poset(([1,2,3],[(1,2),(3,2)]))
            sage: TamariIntervalPosets.check_poset(p)
            True
            sage: p = Poset(([2,3],[(3,2)]))
            sage: TamariIntervalPosets.check_poset(p)
            False
            sage: p = Poset(([1,2,3],[(3,1)]))
            sage: TamariIntervalPosets.check_poset(p)
            False
            sage: p = Poset(([1,2,3],[(1,3)]))
            sage: TamariIntervalPosets.check_poset(p)
            False
        """
        if not set(poset._elements) == set(range(1, poset.cardinality() + 1)):
            return False

        for i in xrange(1, poset.cardinality() + 1):
            stop = False
            for j in xrange(i - 1, 0, -1):
                if not poset.le(j, i):
                    stop = True  # j does not precede i so no j'<j should
                elif stop:
                    return False
            stop = False
            for j in xrange(i + 1, poset.cardinality() + 1):
                if not poset.le(j, i):
                    stop = True  # j does not precede i so no j'>j should
                elif stop:
                    return False
        return True

    @staticmethod
    def final_forest(element):
        r"""
        Return the final forest of a binary tree, an interval-poset or a
        Dyck word.

        A final forest is an interval-poset corresponding to a final
        interval of the Tamari lattice, i.e., containing only decreasing
        relations.

        It can be constructed from a binary tree by its binary
        search tree labeling with the rule: `b` precedes
        `a` in the final forest iff `b` is in the right subtree of `a`
        in the binary search tree.

        INPUT:

        - ``element`` -- a binary tree, a Dyck word or an interval-poset

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: TamariIntervalPosets.final_forest(ip)
            The tamari interval of size 4 induced by relations [(1, 2), (2, 3)]

        from binary trees::

            sage: bt = BinaryTree(); bt
            .
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 0 induced by relations []
            sage: bt = BinaryTree([]); bt
            [., .]
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 1 induced by relations []
            sage: bt = BinaryTree([[],None]); bt
            [[., .], .]
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 2 induced by relations []
            sage: bt = BinaryTree([None,[]]); bt
            [., [., .]]
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 2 induced by relations [(2, 1)]
            sage: bt = BinaryTree([[],[]]); bt
            [[., .], [., .]]
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 3 induced by relations [(3, 2)]
            sage: bt = BinaryTree([[None,[[],None]],[]]); bt
            [[., [[., .], .]], [., .]]
            sage: TamariIntervalPosets.final_forest(bt)
            The tamari interval of size 5 induced by relations [(5, 4), (3, 1), (2, 1)]

        from Dyck words::

            sage: dw = DyckWord([1,0])
            sage: TamariIntervalPosets.final_forest(dw)
            The tamari interval of size 1 induced by relations []
            sage: dw = DyckWord([1,1,0,1,0,0,1,1,0,0])
            sage: TamariIntervalPosets.final_forest(dw)
            The tamari interval of size 5 induced by relations [(5, 4), (3, 1), (2, 1)]
        """
        if isinstance(element, TamariIntervalPoset):
            return element.initial_forest()
        elif element in DyckWords():
            binary_tree = element.to_binary_tree_tamari()
        elif element in BinaryTrees() or element in LabelledBinaryTrees():
            binary_tree = element
        else:
            raise ValueError("Do not know how to construct the initial forest of %s") % str(element)

        def get_relations(bt, start=1):
            r"""
            Recursive method to get the binary tree final forest relations
            with only one recursive reading of the tree.

            The vertices are being labelled with integers starting with
            ``start``.

            OUTPUT:

            - the indexes of the nodes on the left border of the tree
              (these become the roots of the forest)
            - the relations of the final forest (as a list of tuples)
            - the next available index for a node (size of tree +
              ``start``)
            """
            if not bt:
                return [], [], start  # leaf
            roots, relations, index = get_relations(bt[0], start=start)
            rroots, rrelations, rindex = get_relations(bt[1], start=index + 1)
            roots.append(index)
            relations.extend(rrelations)
            relations.extend([(j, index) for j in rroots])
            return roots, relations, rindex

        roots, relations, index = get_relations(binary_tree)
        return TamariIntervalPoset(index - 1, relations)

    @staticmethod
    def initial_forest(element):
        r"""
        Return the inital forest of a binary tree, an interval-poset or
        a Dyck word.

        An initial forest is an interval-poset corresponding to an initial
        interval of the Tamari lattice, i.e., containing only increasing
        relations.

        It can be constructed from a binary tree by its binary
        search tree labeling with the rule: `a` precedes `b` in the
        initial forest iff `a` is in the left subtree of `b` in the
        binary search tree.

        INPUT:

        - ``element`` -- a binary tree, a Dyck word or an interval-poset

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: TamariIntervalPosets.initial_forest(ip)
            The tamari interval of size 4 induced by relations [(1, 2), (2, 3)]

        with binary trees::

            sage: bt = BinaryTree(); bt
            .
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 0 induced by relations []
            sage: bt = BinaryTree([]); bt
            [., .]
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 1 induced by relations []
            sage: bt = BinaryTree([[],None]); bt
            [[., .], .]
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 2 induced by relations [(1, 2)]
            sage: bt = BinaryTree([None,[]]); bt
            [., [., .]]
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 2 induced by relations []
            sage: bt = BinaryTree([[],[]]); bt
            [[., .], [., .]]
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 3 induced by relations [(1, 2)]
            sage: bt = BinaryTree([[None,[[],None]],[]]); bt
            [[., [[., .], .]], [., .]]
            sage: TamariIntervalPosets.initial_forest(bt)
            The tamari interval of size 5 induced by relations [(1, 4), (2, 3), (3, 4)]

        from Dyck words::

            sage: dw = DyckWord([1,0])
            sage: TamariIntervalPosets.initial_forest(dw)
            The tamari interval of size 1 induced by relations []
            sage: dw = DyckWord([1,1,0,1,0,0,1,1,0,0])
            sage: TamariIntervalPosets.initial_forest(dw)
            The tamari interval of size 5 induced by relations [(1, 4), (2, 3), (3, 4)]
        """
        if isinstance(element, TamariIntervalPoset):
            return element.initial_forest()
        elif element in DyckWords():
            binary_tree = element.to_binary_tree_tamari()
        elif element in BinaryTrees() or element in LabelledBinaryTrees():
            binary_tree = element
        else:
            raise ValueError("Do not know how to construct the initial forest of %s") % str(element)

        def get_relations(bt, start=1):
            r"""
            Recursive method to get the binary tree initial forest
            relations with only one recursive reading of the tree.

            The vertices are being labelled with integers starting with
            ``start``.

            OUTPUT:

            - the indexes of the nodes on the right border of the tree
              (these become the roots of the forest)
            - the relations of the initial forest (as a list of tuples)
            - the next available index for a node (size of tree +
              ``start``)
            """
            if not bt:
                return [], [], start  # leaf
            lroots, lrelations, index = get_relations(bt[0], start=start)
            roots, relations, rindex = get_relations(bt[1], start=index + 1)
            roots.append(index)
            relations.extend(lrelations)
            relations.extend([(j, index) for j in lroots])
            return roots, relations, rindex

        roots, relations, index = get_relations(binary_tree)
        return TamariIntervalPoset(index - 1, relations)

    @staticmethod
    def from_binary_trees(tree1, tree2):
        r"""
        Return the interval-poset corresponding to the interval
        [``tree1``,``tree2``] of the Tamari lattice. Raise an exception if
        ``tree1`` is not `\leq` ``tree2`` in the Tamari lattice.

        INPUT:

        - ``tree1`` -- a binary tree
        - ``tree2`` -- a binary tree greater or equal than ``tree1`` for
          the Tamari lattice.

        EXAMPLES::

            sage: tree1 = BinaryTree([[],None])
            sage: tree2 = BinaryTree([None,[]])
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree2)
            The tamari interval of size 2 induced by relations []
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree1)
            The tamari interval of size 2 induced by relations [(1, 2)]
            sage: TamariIntervalPosets.from_binary_trees(tree2,tree2)
            The tamari interval of size 2 induced by relations [(2, 1)]

            sage: tree1 = BinaryTree([[],[[None,[]],[]]])
            sage: tree2 = BinaryTree([None,[None,[None,[[],[]]]]])
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree2)
            The tamari interval of size 6 induced by relations [(4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]

            sage: tree3 = BinaryTree([None,[None,[[],[None,[]]]]])
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree3)
            Traceback (most recent call last):
            ...
            ValueError: The two binary trees are not comparable on the Tamari lattice.
            sage: TamariIntervalPosets.from_binary_trees(tree1,BinaryTree())
            Traceback (most recent call last):
            ...
            ValueError: The two binary trees are not comparable on the Tamari lattice.
        """
        initial_forest = TamariIntervalPosets.initial_forest(tree2)
        final_forest = TamariIntervalPosets.final_forest(tree1)
        try:
            return initial_forest.intersection(final_forest)
        except:
            raise ValueError("The two binary trees are not comparable on the Tamari lattice.")

    @staticmethod
    def from_dyck_words(dw1, dw2):
        r"""
        Return the interval-poset corresponding to the interval
        [``dw1``,``dw2``] of the Tamari lattice. Raise an exception if the
        two Dyck words ``dw1`` and ``dw2`` do not satisfy
        ``dw1`` `\leq` ``dw2`` in the Tamari lattice.

        INPUT:

        - ``dw1`` -- a Dyck word
        - ``dw2`` -- a Dyck word greater or equal than ``dw1`` for the Tamari lattice.

        EXAMPLES::

            sage: dw1 = DyckWord([1,0,1,0])
            sage: dw2 = DyckWord([1,1,0,0])
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw2)
            The tamari interval of size 2 induced by relations []
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw1)
            The tamari interval of size 2 induced by relations [(1, 2)]
            sage: TamariIntervalPosets.from_dyck_words(dw2,dw2)
            The tamari interval of size 2 induced by relations [(2, 1)]

            sage: dw1 = DyckWord([1,0,1,1,1,0,0,1,1,0,0,0])
            sage: dw2 = DyckWord([1,1,1,1,0,1,1,0,0,0,0,0])
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw2)
            The tamari interval of size 6 induced by relations [(4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]

            sage: dw3 = DyckWord([1,1,1,0,1,1,1,0,0,0,0,0])
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw3)
            Traceback (most recent call last):
            ...
            ValueError: The two Dyck words are not comparable on the Tamari lattice.
            sage: TamariIntervalPosets.from_dyck_words(dw1,DyckWord([1,0]))
            Traceback (most recent call last):
            ...
            ValueError: The two Dyck words are not comparable on the Tamari lattice.
        """
        tree1 = dw1.to_binary_tree_tamari()
        tree2 = dw2.to_binary_tree_tamari()
        try:
            return TamariIntervalPosets.from_binary_trees(tree1, tree2)
        except:
            raise ValueError("The two Dyck words are not comparable on the Tamari lattice.")

    @staticmethod
    def lower_contacts_composition(ip1, ip2, p = None):
        r"""
        Return the lower contacts composition of ``ip1`` and ``ip2`` as 
        described in [CCP]_. If ``p`` is not ``None``, then the partial 
        composition is performed.

        If ``ip1.size()`` is `n` and ``ip2.size()`` is `m`, the complete lower 
        contacts composition of ``ip1`` and ``ip2`` consists of a set of 
        interval-posets of size `n+m+1`. If ``ip1.number_of_contacts()`` is `c_1` 
        and ``ip2.number_of_contacts()`` is `c_2`, then the number of interval-posets 
        of the composition is `c_2+1` and the number of contacts of the 
        resulting interval-posets are `c_1 + c_2 + 1, c_1 + c_2, c_1 + c_2 - 1, \dots, c_1 + 1`.
        
        If ``p`` is not ``None``, then only a partial composition is performed 
        and a unique interval-poset is returned with its number of contacts
        equal to `c_1 + 1 + p`.

        For a detailed description of the composition, please refer to [CCP]_.

        REFERENCES:

        .. [CCP] Two bijections on Tamari intervals, G. Chatel, V. Pons

        INPUT:

        - ``ip1`` -- an interval-poset of size `n`
        - ``ip2`` -- an interval-poset of size `m`
        - ``p`` -- (default: ``None``) a Integer between 0 and ``ip2.number_of_contacts()``

        OUTPUT:

        A list of ``ip2.number_of_contacts()`` interval-posets of size `n+m+1`
        if ``p`` is ``None`` or a unique interval-poset of size `n+m+1` if 
        ``p`` is a number.

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(2,1),(3,1),(2,4),(3,4)])
            sage: ip2 = TamariIntervalPoset(5,[(1,2),(3,2),(5,4)])
            sage: TamariIntervalPosets.lower_contacts_composition(ip1,ip2)
            [The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (8, 7), (3, 1), (2, 1)],
             The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (8, 7), (6, 5), (3, 1), (2, 1)],
             The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (8, 7), (7, 5), (6, 5), (3, 1), (2, 1)],
             The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (9, 5), (8, 7), (7, 5), (6, 5), (3, 1), (2, 1)]]
            sage: ip1.number_of_contacts()
            2
            sage: ip2.number_of_contacts()
            3
            sage: [ip.number_of_contacts() for ip in TamariIntervalPosets.lower_contacts_composition(ip1,ip2)]
            [6, 5, 4, 3]
            
        partial composition::
    
            sage: TamariIntervalPosets.lower_contacts_composition(ip1,ip2,2)
            The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (8, 7), (7, 5), (6, 5), (3, 1), (2, 1)]
            sage: TamariIntervalPosets.lower_contacts_composition(ip1,ip2,0)
            The tamari interval of size 10 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (6, 7), (10, 9), (8, 7), (3, 1), (2, 1)]
            sage: TamariIntervalPosets.lower_contacts_composition(ip1,ip2,5)
            Traceback (most recent call last):
            ...
            ValueError: Invalid composition parameter
            sage: [TamariIntervalPosets.lower_contacts_composition(ip1,ip2,p) for p in xrange(ip2.number_of_contacts()+1)] == TamariIntervalPosets.lower_contacts_composition(ip1,ip2)
            True

        composition with an empty interval-poset::

            sage: ip0 = TamariIntervalPoset(0,[])
            sage: TamariIntervalPosets.lower_contacts_composition(ip1,ip0)
            [The tamari interval of size 5 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (3, 1), (2, 1)]]
            sage: TamariIntervalPosets.lower_contacts_composition(ip1,ip0,0)
            The tamari interval of size 5 induced by relations [(1, 5), (2, 4), (3, 4), (4, 5), (3, 1), (2, 1)]
            sage: TamariIntervalPosets.lower_contacts_composition(ip0,ip2)
            [The tamari interval of size 6 induced by relations [(2, 3), (6, 5), (4, 3)],
             The tamari interval of size 6 induced by relations [(2, 3), (6, 5), (4, 3), (2, 1)],
             The tamari interval of size 6 induced by relations [(2, 3), (6, 5), (4, 3), (3, 1), (2, 1)],
             The tamari interval of size 6 induced by relations [(2, 3), (6, 5), (5, 1), (4, 3), (3, 1), (2, 1)]]
            sage: TamariIntervalPosets.lower_contacts_composition(ip0,ip2,2)
            The tamari interval of size 6 induced by relations [(2, 3), (6, 5), (4, 3), (3, 1), (2, 1)]
            sage: TamariIntervalPosets.lower_contacts_composition(ip0,ip0)
            [The tamari interval of size 1 induced by relations []]

        TESTS::

            sage: set(TamariIntervalPosets(4)) == set([ip for comp in [TamariIntervalPosets.lower_contacts_composition(ip1,ip2) for i in xrange(4) for ip1 in TamariIntervalPosets(i) for ip2 in TamariIntervalPosets(3-i)] for ip in comp])
            True

        """
        if p is not None and (p<0 or p>ip2.number_of_contacts()):
            raise ValueError, "Invalid composition parameter"%()

        k = ip1.size()+1
        size = k + ip2.size()

        ### left part ###

        # we add all relations of ip1
        relations = list(ip1._cover_relations)
        # we add all relations i-->k for i in ip1
        relations.extend([(i,k) for i in xrange(1,k)])
        
        ### right part ###

        # we add all relations of ip2
        relations.extend([(i+k,j+k) for (i,j) in ip2._cover_relations])
        # we add the decreasing relations depending on p
        result = []
        if p is None:
            result.append(TamariIntervalPoset(size,relations))
        elif p==0:
            return TamariIntervalPoset(size,relations)
        c = 0
        for j in ip2.decreasing_roots():
            relations.append((j+k,k))
            c+=1
            if p is None:
                result.append(TamariIntervalPoset(size,relations))
            elif p==c:
                return TamariIntervalPoset(size,relations)
        return result

    @staticmethod
    def initial_rise_composition(ip1,ip2,p=None):
        r"""
        Return the initial rise composition of ``ip1`` and ``ip2`` as 
        described in [CCP]_. If ``p`` is not ``None``, then the partial 
        composition is performed.

        If ``ip1.size()`` is `n` and ``ip2.size()`` is `m`, the complete initial 
        rise composition of ``ip1`` and ``ip2`` consists of a set of 
        interval-posets of size `n+m+1`. If ``ip1.initial_rise()`` is `ir_1` 
        and ``ip2.initial_rise()`` is `ir_2`, then the number of interval-posets 
        of the composition is `ir_2+ir_1` and the initial rises of the 
        resulting interval-posets are `ir_1 + ir_2 + 1, ir_1 + ir_2, ir_1 + ir_2 - 1, \dots, ir_1 + 1`.
        
        If ``p`` is not ``None``, then only a partial composition is performed 
        and a unique interval-poset is returned with its initial rise
        equal to `ir_1 + 1 + p`.

        For a detailed description of the composition, please refer to [CCP]_.

        REFERENCES:

        .. [CCP] Two bijections on Tamari intervals, G. Chatel, V. Pons

        INPUT:

        - ``ip1`` -- an interval-poset of size `n`
        - ``ip2`` -- an interval-poset of size `m`
        - ``p`` -- (default: ``None``) a Integer between 0 and ``ip2.initial_rise()``

        OUTPUT:

        A list of ``ip2.initial_rise()`` interval-posets of size `n+m+1`
        if ``p`` is ``None`` or a unique interval-poset of size `n+m+1` if 
        ``p`` is a number.

        EXAMPLES::
        
            sage: ip1 = TamariIntervalPoset(3,[(1,2),(3,2)])
            sage: ip2 = TamariIntervalPoset(5,[(2,3),(3,5),(4,5),(2,1),(3,1),(4,3)])
            sage: TamariIntervalPosets.initial_rise_composition(ip1,ip2)
            [The tamari interval of size 9 induced by relations [(1, 8), (2, 8), (3, 5), (4, 5), (5, 7), (6, 7), (7, 8), (9, 8), (7, 1), (6, 5), (5, 1), (4, 2), (3, 2), (2, 1)],
             The tamari interval of size 9 induced by relations [(1, 8), (2, 8), (3, 4), (4, 5), (5, 7), (6, 7), (7, 8), (9, 8), (7, 1), (6, 5), (5, 1), (4, 2), (3, 2), (2, 1)],
             The tamari interval of size 9 induced by relations [(1, 8), (2, 3), (3, 8), (4, 5), (5, 7), (6, 7), (7, 8), (9, 8), (7, 1), (6, 5), (5, 3), (4, 3), (3, 1), (2, 1)]]
            sage: ip1.initial_rise()
            1
            sage: ip2.initial_rise()
            2
            sage: [ip.initial_rise() for ip in TamariIntervalPosets.initial_rise_composition(ip1,ip2)]
            [4, 3, 2]

        partial composition::

            sage: TamariIntervalPosets.initial_rise_composition(ip1,ip2,2)
            The tamari interval of size 9 induced by relations [(1, 8), (2, 3), (3, 8), (4, 5), (5, 7), (6, 7), (7, 8), (9, 8), (7, 1), (6, 5), (5, 3), (4, 3), (3, 1), (2, 1)]
            sage: TamariIntervalPosets.initial_rise_composition(ip1,ip2,0)
            The tamari interval of size 9 induced by relations [(1, 8), (2, 8), (3, 5), (4, 5), (5, 7), (6, 7), (7, 8), (9, 8), (7, 1), (6, 5), (5, 1), (4, 2), (3, 2), (2, 1)]
            sage: TamariIntervalPosets.initial_rise_composition(ip1,ip2,4)
            Traceback (most recent call last):
            ...
            ValueError: Invalid composition parameter
            sage: [TamariIntervalPosets.initial_rise_composition(ip1,ip2,p) for p in xrange(ip2.initial_rise()+1)] == TamariIntervalPosets.initial_rise_composition(ip1,ip2)
            True

        composition with an empty interval-poset::

            sage: ip0 = TamariIntervalPoset(0,[])
            sage: TamariIntervalPosets.initial_rise_composition(ip1,ip0)
            [The tamari interval of size 4 induced by relations [(1, 3), (2, 3), (4, 3), (2, 1)]]
            sage: TamariIntervalPosets.initial_rise_composition(ip1,ip0,0)
            The tamari interval of size 4 induced by relations [(1, 3), (2, 3), (4, 3), (2, 1)]
            sage: TamariIntervalPosets.initial_rise_composition(ip0,ip2)
            [The tamari interval of size 6 induced by relations [(2, 4), (3, 4), (4, 6), (5, 6), (5, 4), (3, 1), (2, 1)],
             The tamari interval of size 6 induced by relations [(2, 3), (3, 4), (4, 6), (5, 6), (5, 4), (3, 1), (2, 1)],
             The tamari interval of size 6 induced by relations [(1, 2), (3, 4), (4, 6), (5, 6), (5, 4), (4, 2), (3, 2)]]
            sage: TamariIntervalPosets.initial_rise_composition(ip0,ip2,2)
            The tamari interval of size 6 induced by relations [(1, 2), (3, 4), (4, 6), (5, 6), (5, 4), (4, 2), (3, 2)]
            sage: TamariIntervalPosets.initial_rise_composition(ip0,ip0)
            [The tamari interval of size 1 induced by relations []]

        TESTS::

            sage: set(TamariIntervalPosets(4)) == set([ip for comp in [TamariIntervalPosets.initial_rise_composition(ip1,ip2) for i in xrange(4) for ip1 in TamariIntervalPosets(i) for ip2 in TamariIntervalPosets(3-i)] for ip in comp])
            True

        """
        ir1 = ip1.initial_rise()
        ir2 = ip2.initial_rise()
        if p is None:
            return [TamariIntervalPosets.initial_rise_composition(ip1,ip2,i) for i in xrange(ir2+1)]
        if p<0 or p>ir2:
            raise ValueError, "Invalid composition parameter"%()
        relations = []
        n2 = ip2.size()
        n1 = ip1.size()

        #### relations from ip2 ####

        root = ir2 + 1 -p

        ## increasing relations ##
        for (i,j) in ip2.increasing_cover_relations():
            if i>=root:
                i+=1
            if j>=root:
                j+=1
            relations.append((i+ir1,j+ir1))

        if root!= ip2.size() +1:
            relations.append((root+ir1,root+ir1+1))

        ## decreasing relations ##

        children_decreasing = []
        for i in xrange(n2,root-1,-1):
            while len(children_decreasing)>0:
                j = children_decreasing.pop()
                if ip2.le(j,i):
                    relations.append((j+1+ir1,i+1+ir1))
                else:
                    children_decreasing.append(j)
                    break 
            children_decreasing.append(i)
        new_children = [i+1 for i in children_decreasing]
        if(root!= n2+1): new_children.append(root)
        for i in xrange(root-1,0,-1):
            while len(children_decreasing)>0:
                j = children_decreasing.pop()
                nj = new_children.pop()
                if ip2.le(j,i):
                    relations.append((nj+ir1,i+ir1))
                else:
                    children_decreasing.append(j)
                    new_children.append(nj)
                    break
            children_decreasing.append(i)
            new_children.append(i)
            
        #### relations from ip1 ####

        for (i,j) in ip1._cover_relations:
            if i>ir1:
                i+=n2+1
            if j>ir1:
                j+=n2+1
            relations.append((i,j))

        # extra increasing relations
        par = ip1.increasing_parent(ir1)
        if par!=None:
            for i in xrange(ir1+1,ir1+n2+2):
                relations.append((i,par+n2+1))
                
        # extra decreasing relations
        if(ir1!=0):
            for i in xrange(ir1+1,ir1+n2+2):
                relations.append((i,ir1))
        
        return TamariIntervalPoset(n1+n2+1,relations)

    def __call__(self, *args, **keywords):
        r"""
        Allows for a poset to be directly transformed into an interval-poset.

        It is some kind of coercion but cannot be made through the coercion
        system because posets do not have parents.

        EXAMPLES::

            sage: TIP = TamariIntervalPosets()
            sage: p = Poset( ([1,2,3], [(1,2)]))
            sage: TIP(p)
            The tamari interval of size 3 induced by relations [(1, 2)]
            sage: TIP(TIP(p))
            The tamari interval of size 3 induced by relations [(1, 2)]
            sage: TIP(3,[(1,2)])
            The tamari interval of size 3 induced by relations [(1, 2)]
            sage: p = Poset(([1,2,3],[(1,3)]))
            sage: TIP(p)
            Traceback (most recent call last):
            ...
            ValueError: This does not satisfy the Tamari interval-poset condition.
        """
        if isinstance(args[0], TamariIntervalPoset):
            return args[0]
        if len(args) == 1 and isinstance(args[0], FinitePoset):
            return self.element_class(args[0].cardinality(), args[0].cover_relations())

        return super(TamariIntervalPosets, self).__call__(*args, **keywords)

    def le(self, el1, el2):
        r"""
        Poset stucture on the set of interval-posets through interval
        containment.

        Return whether the interval represented by ``el1`` is contained in
        the interval represented by ``el2``.

        INPUT:

        - ``el1`` -- an interval-poset
        - ``el2`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: TamariIntervalPosets().le(ip1,ip2)
            True
            sage: TamariIntervalPosets().le(ip2,ip1)
            False
        """
        return el2.contains_interval(el1)

    global_options = TamariIntervalPosetOptions


#################################################################
# Enumerated set of all Tamari Interval-posets
#################################################################
class TamariIntervalPosets_all(DisjointUnionEnumeratedSets, TamariIntervalPosets):

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.combinat.interval_posets import TamariIntervalPosets_all
            sage: S = TamariIntervalPosets_all()
            sage: S.cardinality()
            +Infinity

            sage: it = iter(S)
            sage: [it.next() for i in xrange(5)]
            [The tamari interval of size 0 induced by relations [],
             The tamari interval of size 1 induced by relations [],
             The tamari interval of size 2 induced by relations [],
             The tamari interval of size 2 induced by relations [(2, 1)],
             The tamari interval of size 2 induced by relations [(1, 2)]]
            sage: it.next().parent()
            Interval-posets
            sage: S(0,[])
            The tamari interval of size 0 induced by relations []

            sage: S is TamariIntervalPosets_all()
            True
            sage: TestSuite(S).run()
            """
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), TamariIntervalPosets_size),
            facade=True, keepkey=False, category=(Posets(), EnumeratedSets()))

    def _repr_(self):
        r"""
        TEST::

            sage: TamariIntervalPosets()  # indirect doctest
            Interval-posets
        """
        return "Interval-posets"

    def _element_constructor_(self, size, relations):
        r"""
        EXAMPLES::

            sage: TIP = TamariIntervalPosets()
            sage: TIP(3,[(1,2)])  # indirect doctest
            The tamari interval of size 3 induced by relations [(1, 2)]
        """
        return self.element_class(size, relations)

    def __contains__(self, x):
        r"""
        TESTS::

            sage: S = TamariIntervalPosets()
            sage: 1 in S
            False
            sage: S(0,[]) in S
            True
        """
        return isinstance(x, self.element_class)

    Element = TamariIntervalPoset


#################################################################
# Enumerated set of Tamari interval-posets of a given size
#################################################################
class TamariIntervalPosets_size(TamariIntervalPosets):
    r"""
    The enumerated set of interval-posets of a given size.

    TESTS::

        sage: from sage.combinat.interval_posets import TamariIntervalPosets_size
        sage: for i in xrange(6): TestSuite(TamariIntervalPosets_size(i)).run()
    """
    def __init__(self, size):
        r"""
        TESTS::

            sage: S = TamariIntervalPosets(3)
            sage: S == loads(dumps(S))
            True

            sage: S is TamariIntervalPosets(3)
            True
        """
        # there is a natural order on interval-posets through inclusions
        # that is why we use the FinitePosets category
        super(TamariIntervalPosets_size, self).__init__(category=(FinitePosets(), FiniteEnumeratedSets()))

        self._size = size

    def _repr_(self):
        r"""
        TESTS::

            sage: TamariIntervalPosets(3)
            Interval-posets of size 3
        """
        return "Interval-posets of size %s" % self._size

    def __contains__(self, x):
        r"""
        TESTS::

            sage: S = TamariIntervalPosets(3)
            sage: 1 in S
            False
            sage: S([]) in S
            True
        """
        return isinstance(x, self.element_class) and x.size() == self._size

    def cardinality(self):
        r"""
        The cardinality of ``self``. That is, the number of
        interval-posets of size `n`.

        The formula was given in [ChapTamari08]: `\frac{2(4n+1)!}{(n+1)!(3n+2)!}`.

        EXAMPLES::

            sage: [TamariIntervalPosets(i).cardinality() for i in range(6)]
            [1, 1, 3, 13, 68, 399]
        """
        from sage.rings.arith import binomial
        n = self._size
        if n == 0:
            return Integer(1)
        return (2 * binomial(4 * n + 1, n - 1)) // (n * (n + 1))
        # return Integer(2 * factorial(4*n+1)/(factorial(n+1)*factorial(3*n+2)))

    def __iter__(self):
        r"""
        Recursive generation: we iterate through all interval-posets of
        size ``size - 1`` and add all possible relations to the last
        vertex.

        TESTS::

            sage: TIP1 = TamariIntervalPosets(1)
            sage: list(TIP1)
            [The tamari interval of size 1 induced by relations []]
            sage: TIP2 = TamariIntervalPosets(2)
            sage: list(TIP2)
            [The tamari interval of size 2 induced by relations [],
             The tamari interval of size 2 induced by relations [(2, 1)],
             The tamari interval of size 2 induced by relations [(1, 2)]]
            sage: TIP3 = TamariIntervalPosets(3)
            sage: list(TIP3)
            [The tamari interval of size 3 induced by relations [],
             The tamari interval of size 3 induced by relations [(3, 2)],
             The tamari interval of size 3 induced by relations [(2, 3)],
             The tamari interval of size 3 induced by relations [(1, 3), (2, 3)],
             The tamari interval of size 3 induced by relations [(2, 1)],
             The tamari interval of size 3 induced by relations [(3, 2), (2, 1)],
             The tamari interval of size 3 induced by relations [(3, 1), (2, 1)],
             The tamari interval of size 3 induced by relations [(2, 3), (2, 1)],
             The tamari interval of size 3 induced by relations [(2, 3), (3, 1), (2, 1)],
             The tamari interval of size 3 induced by relations [(1, 3), (2, 3), (2, 1)],
             The tamari interval of size 3 induced by relations [(1, 2)],
             The tamari interval of size 3 induced by relations [(1, 2), (3, 2)],
             The tamari interval of size 3 induced by relations [(1, 2), (2, 3)]]
            sage: all([len(list(TamariIntervalPosets(i)))==TamariIntervalPosets(i).cardinality() for i in xrange(6)])
            True
        """
        n = self._size
        if n == 0:
            yield TamariIntervalPoset(0, [])
        elif n == 1:
            yield TamariIntervalPoset(1, [])
        else:
            for tip in TamariIntervalPosets(n - 1):
                new_tip = TamariIntervalPoset(n, tip._cover_relations)
                yield new_tip  # we have added an extra vertex but no relations

                # adding a decreasing relation n>>m2 with m2<n
                for m2 in xrange(n - 1, 0, -1):
                    if new_tip.le(n - 1, m2):
                        yield TamariIntervalPoset(n, new_tip._cover_relations + ((n, m2),))

                for m in xrange(n - 1, 0, -1):

                    # adding an increasing relation m>>n
                    if not new_tip.le(m, n):
                        new_tip = TamariIntervalPoset(n, new_tip._cover_relations + ((m, n),))
                        yield new_tip
                    else:
                        continue

                    # adding a decreasing relation n>>m2 with m2<m
                    for m2 in xrange(m - 1, 0, -1):
                        if new_tip.le(n - 1, m2):
                            yield TamariIntervalPoset(n, new_tip._cover_relations + ((n, m2),))

    @lazy_attribute
    def _parent_for(self):
        r"""
        The parent of the element generated by ``self``

        TESTS::

            sage: TIP3 = TamariIntervalPosets(3)
            sage: TIP3._parent_for
            Interval-posets
        """
        return TamariIntervalPosets_all()

    @lazy_attribute
    def element_class(self):
        r"""
        TESTS::

            sage: S = TamariIntervalPosets(3)
            sage: S.element_class
            <class 'sage.combinat.interval_posets.TamariIntervalPosets_all_with_category.element_class'>
            sage: S.first().__class__ == TamariIntervalPosets().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, relations):
        r"""
        EXAMPLES::

            sage: TIP3 = TamariIntervalPosets(3)
            sage: TIP3([(1,2)])  # indirect doctest
            The tamari interval of size 3 induced by relations [(1, 2)]
            sage: TIP3([(3,4)])
            Traceback (most recent call last):
            ...
            ValueError: The relations do not correspond to the size of the poset.
        """
        return self.element_class(self._size, relations)
